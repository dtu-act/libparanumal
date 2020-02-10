/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "acoustics.h"
#include <stdio.h>

acoustics_t *acousticsSetup(mesh_t *mesh, setupAide &newOptions, char* boundaryHeaderFileName){
  acoustics_t *acoustics = (acoustics_t*) calloc(1, sizeof(acoustics_t));

  newOptions.getArgs("MESH DIMENSION", acoustics->dim);
  newOptions.getArgs("ELEMENT TYPE", acoustics->elementType);
  
  mesh->Nfields = (acoustics->dim==3) ? 4:3;
  acoustics->Nfields = mesh->Nfields;
  
  acoustics->mesh = mesh;

  dlong Ntotal = mesh->Nelements*mesh->Np*mesh->Nfields;
  acoustics->Nblock = (Ntotal+blockSize-1)/blockSize;
  
  hlong localElements = (hlong) mesh->Nelements;
  
  MPI_Allreduce(&localElements, &(acoustics->totalElements), 1, MPI_HLONG, MPI_SUM, mesh->comm);

  // viscosity
  int check;

  // compute samples of q at interpolation nodes
  acoustics->q = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  acoustics->rhsq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
  
  if (newOptions.compareArgs("TIME INTEGRATOR","LSERK4")){
    acoustics->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Nfields,
		  		sizeof(dfloat));
  }

  if (newOptions.compareArgs("TIME INTEGRATOR","DOPRI5")){
    int NrkStages = 7;
    acoustics->rkq  = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
          sizeof(dfloat));
    acoustics->rkrhsq = (dfloat*) calloc(NrkStages*mesh->Nelements*mesh->Np*mesh->Nfields,
          sizeof(dfloat));
    acoustics->rkerr  = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
          sizeof(dfloat));

    acoustics->errtmp = (dfloat*) calloc(acoustics->Nblock, sizeof(dfloat));

    // Dormand Prince -order (4) 5 with PID timestep control
    int Nrk = 7;
    dfloat rkC[7] = {0.0, 0.2, 0.3, 0.8, 8.0/9.0, 1.0, 1.0};
    dfloat rkA[7*7] ={             0.0,             0.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                                   0.2,             0.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                              3.0/40.0,        9.0/40.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                             44.0/45.0,      -56.0/15.0,       32.0/9.0,          0.0,             0.0,       0.0, 0.0,
                        19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0,             0.0,       0.0, 0.0,
                         9017.0/3168.0,     -355.0/33.0, 46732.0/5247.0,   49.0/176.0, -5103.0/18656.0,       0.0, 0.0, 
                            35.0/384.0,             0.0,   500.0/1113.0,  125.0/192.0,  -2187.0/6784.0, 11.0/84.0, 0.0 };
    dfloat rkE[7] = {71.0/57600.0,  0.0, -71.0/16695.0, 71.0/1920.0, -17253.0/339200.0, 22.0/525.0, -1.0/40.0 }; 

    acoustics->Nrk = Nrk;
    acoustics->rkC = (dfloat*) calloc(acoustics->Nrk, sizeof(dfloat));
    acoustics->rkE = (dfloat*) calloc(acoustics->Nrk, sizeof(dfloat));
    acoustics->rkA = (dfloat*) calloc(acoustics->Nrk*acoustics->Nrk, sizeof(dfloat));

    memcpy(acoustics->rkC, rkC, acoustics->Nrk*sizeof(dfloat));
    memcpy(acoustics->rkE, rkE, acoustics->Nrk*sizeof(dfloat));
    memcpy(acoustics->rkA, rkA, acoustics->Nrk*acoustics->Nrk*sizeof(dfloat));
    
    acoustics->dtMIN = 1E-9; //minumum allowed timestep
    acoustics->ATOL = 1E-6;  //absolute error tolerance
    acoustics->RTOL = 1E-6;  //relative error tolerance
    acoustics->safe = 0.8;   //safety factor

    //error control parameters
    acoustics->beta = 0.05;
    acoustics->factor1 = 0.2;
    acoustics->factor2 = 10.0;


    acoustics->exp1 = 0.2 - 0.75*acoustics->beta;
    acoustics->invfactor1 = 1.0/acoustics->factor1;
    acoustics->invfactor2 = 1.0/acoustics->factor2;
    acoustics->facold = 1E-4;
    
  }
  dfloat sloc[3];
  dfloat sxyz;
  newOptions.getArgs("SX", sloc[0]);
  newOptions.getArgs("SY", sloc[1]);
  newOptions.getArgs("SZ", sloc[2]);
  newOptions.getArgs("SXYZ", sxyz);
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat t = 0;
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];
      dfloat z = mesh->z[n + mesh->Np*e];

      dlong qbase = e*mesh->Np*mesh->Nfields + n;

      dfloat u = 0, v = 0, w = 0, r = 0;
      
      acousticsGaussianPulse(x, y, z, t, &r, &u, &v, &w, sloc, sxyz);
      acoustics->q[qbase+0*mesh->Np] = r;
      acoustics->q[qbase+1*mesh->Np] = u;
      acoustics->q[qbase+2*mesh->Np] = v;
      if(acoustics->dim==3)
	      acoustics->q[qbase+3*mesh->Np] = w;
    }
  }

  // set penalty parameter
  mesh->Lambda2 = 0.5;
  
  // set time step
  dfloat hmin = 1e9;
  for(dlong e=0;e<mesh->Nelements;++e){  

    for(int f=0;f<mesh->Nfaces;++f){
      dlong sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];

      if(invJ<0) printf("invJ = %g\n", invJ);
      
      // sJ = L/2, J = A/2,   sJ/J = L/A = L/(0.5*h*L) = 2/h
      // h = 0.5/(sJ/J)
      
      dfloat hest = .5/(sJ*invJ);

      hmin = mymin(hmin, hest);
    }
  }
  // need to change cfl and defn of dt
  //dfloat cfl = 0.5; // depends on the stability region size
  dfloat cfl;
  newOptions.getArgs("CFL", cfl);

  dfloat dtAdv  = hmin/(343.0*(mesh->N+1.)*(mesh->N+1.));
  //dfloat dtAdv  = hmin/(343.0);
  dfloat dt = cfl*dtAdv;
  
  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(mesh->dt), 1, MPI_DFLOAT, MPI_MIN, mesh->comm);
  
  //
  newOptions.getArgs("FINAL TIME", mesh->finalTime);

  mesh->NtimeSteps = mesh->finalTime/mesh->dt;
  if (newOptions.compareArgs("TIME INTEGRATOR","LSERK4")){
    mesh->dt = mesh->finalTime/mesh->NtimeSteps;
  }
    if (newOptions.compareArgs("TIME INTEGRATOR","EIRK4")){
    mesh->dt = mesh->finalTime/mesh->NtimeSteps;
  }
  
  if (mesh->rank ==0) printf("dtAdv = %lg (before cfl), dt = %lg\n",
   dtAdv, dt);

  //---------RECEIVER---------
  acoustics->qRecvCounter = 0; // Counter needed for later
  acoustics->qRecvCopyCounter = 0;


  // Read from receiver locations file
  string RecvDATAFileName;
  newOptions.getArgs("RECEIVER", RecvDATAFileName);

  FILE * RecvDATAFILE = fopen((char*)RecvDATAFileName.c_str(),"r");
  if (RecvDATAFILE == NULL) {
    printf("Could not find Reciver Locations file: %s\n", (char*)RecvDATAFileName.c_str());
    exit(-1);
  }

  fscanf(RecvDATAFILE,"%d",&acoustics->NReceivers);
  acoustics->recvXYZ = (dfloat*) calloc(acoustics->NReceivers*3, sizeof(dfloat));

  for(int iRead = 0; iRead < acoustics->NReceivers*3; iRead+=3){
    fscanf(RecvDATAFILE,"%lf %lf %lf",&acoustics->recvXYZ[iRead],
                                      &acoustics->recvXYZ[iRead+1],
                                      &acoustics->recvXYZ[iRead+2]);
  }
  fclose(RecvDATAFILE);
  acoustics->recvElements = (dlong*) calloc(acoustics->NReceivers, sizeof(dlong));
  // Assume that no receivers on this core
  for(int i = 0; i < acoustics->NReceivers;i++){
    acoustics->recvElements[i] = -1;
  } 
  acoustics->NReceiversLocal = 0; 
  
  acoustics->recvElementsIdx = (dlong*) calloc(acoustics->NReceivers, sizeof(dlong));

  //---------RECEIVER---------

  
  acoustics->frame = 0;
  // errorStep
  mesh->errorStep = 1000;

  if (mesh->rank ==0) printf("dt = %g\n", mesh->dt);

  // OCCA build stuff
  
  occa::properties kernelInfo;
 kernelInfo["defines"].asObject();
 kernelInfo["includes"].asArray();
 kernelInfo["header"].asArray();
 kernelInfo["flags"].asObject();

  if(acoustics->dim==3)
    meshOccaSetup3D(mesh, newOptions, kernelInfo);
  else
    meshOccaSetup2D(mesh, newOptions, kernelInfo);

  //add boundary data to kernel info
  kernelInfo["includes"] += boundaryHeaderFileName;
 
  acoustics->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), acoustics->q);

  acoustics->o_saveq =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), acoustics->q);
  
  acoustics->o_rhsq =
    mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), acoustics->rhsq);
  
  // [EA] Snapshot of solution q
  newOptions.getArgs("SNAPSHOT", acoustics->snapshot);
  if(acoustics->snapshot){
    newOptions.getArgs("SNAPSHOTMAX", acoustics->snapshotMax);
    acoustics->writeSnapshotEvery = 10; // Hardcoded.
    acoustics->Snapshott = (dfloat*) calloc(acoustics->writeSnapshotEvery, sizeof(dfloat));
    acoustics->qSnapshot = (dfloat*) calloc(acoustics->writeSnapshotEvery*mesh->Np*mesh->Nelements*mesh->Nfields, sizeof(dfloat));
    //acoustics->o_qSnapshot = 
    //    mesh->device.malloc(acoustics->writeSnapshotEvery*mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), acoustics->qSnapshot);
  }
  // [EA] Read and allocate space for LR/ER accumulators
  dlong Nangles = 91; // [EA] Number of angles for ER
  if(mesh->NERFaces){
    string ERDATAFileName;
    newOptions.getArgs("ERVECTFIT", ERDATAFileName);

    FILE * ERDATAFILE = fopen((char*)ERDATAFileName.c_str(),"r");
    if (ERDATAFILE == NULL) {
      printf("Could not find ERDATA file: %s\n", (char*)ERDATAFileName.c_str());
      exit(-1);
    }

    

    fscanf(ERDATAFILE,"%d %d %d",&acoustics->ERNpoles, &acoustics->ERNRealPoles,&acoustics->ERNImagPoles);

    acoustics->ER = (dfloat*) calloc(Nangles + acoustics->ERNRealPoles+acoustics->ERNRealPoles*Nangles 
    + 2*acoustics->ERNImagPoles*Nangles+2*acoustics->ERNImagPoles, sizeof(dfloat));
    acoustics->ERInfo = (dlong*) calloc(3,sizeof(dlong));

    acoustics->ERInfo[0] = acoustics->ERNpoles;
    acoustics->ERInfo[1] = acoustics->ERNRealPoles;
    acoustics->ERInfo[2] = acoustics->ERNImagPoles;


    dlong iReadOffset = 0;
    for(int iRead = 0; iRead < acoustics->ERNRealPoles*Nangles; iRead++){
      fscanf(ERDATAFILE,"%lf",&acoustics->ER[iReadOffset+iRead]);
    }
    iReadOffset += acoustics->ERNRealPoles*Nangles;
    for(int iRead = 0; iRead < acoustics->ERNImagPoles*Nangles; iRead++){
      fscanf(ERDATAFILE,"%lf",&acoustics->ER[iReadOffset+iRead]);
    }
    iReadOffset += acoustics->ERNImagPoles*Nangles;
    for(int iRead = 0; iRead < acoustics->ERNImagPoles*Nangles; iRead++){
      fscanf(ERDATAFILE,"%lf",&acoustics->ER[iReadOffset+iRead]);
    }
    iReadOffset += acoustics->ERNImagPoles*Nangles;
    for(int iRead = 0; iRead < acoustics->ERNRealPoles; iRead++){
      fscanf(ERDATAFILE,"%lf",&acoustics->ER[iReadOffset+iRead]);
    }
    iReadOffset += acoustics->ERNRealPoles;
    for(int iRead = 0; iRead < acoustics->ERNImagPoles; iRead++){
      fscanf(ERDATAFILE,"%lf",&acoustics->ER[iReadOffset+iRead]);
    }
    iReadOffset += acoustics->ERNImagPoles;
    for(int iRead = 0; iRead < acoustics->ERNImagPoles; iRead++){
      fscanf(ERDATAFILE,"%lf",&acoustics->ER[iReadOffset+iRead]);
      //printf("beta = %f\n",acoustics->ER[iReadOffset+iRead]);
    }
    iReadOffset += acoustics->ERNImagPoles;
    for(int iRead = 0; iRead < Nangles; iRead++){
      fscanf(ERDATAFILE,"%lf",&acoustics->ER[iReadOffset+iRead]);
    }
    fclose(ERDATAFILE);


    acoustics->o_ER = mesh->device.malloc((Nangles + acoustics->ERNRealPoles+acoustics->ERNRealPoles*Nangles 
    + 2*acoustics->ERNImagPoles*Nangles+2*acoustics->ERNImagPoles)*sizeof(dfloat), acoustics->ER);
    acoustics->o_ERInfo = mesh->device.malloc(3*sizeof(dlong), acoustics->ERInfo);

  } else {
      // [EA] To avoid empty pointers as this causes crashes on gpu
      acoustics->ERNpoles = 1;
      acoustics->ERNRealPoles = 0;
      acoustics->ERNImagPoles = 0;

      acoustics->ERInfo = (dlong*) calloc(3,sizeof(dlong));
      acoustics->ERInfo[0] = 0;
      acoustics->ERInfo[1] = 0;
      acoustics->ERInfo[2] = 0;
      acoustics->o_ER = mesh->device.malloc(1*sizeof(dfloat));
      acoustics->o_ERInfo = mesh->device.malloc(3*sizeof(dlong), acoustics->ERInfo);
  }
  newOptions.getArgs("BCCHANGETIME", acoustics->BCChangeTime);
  if(mesh->NLRFaces || (acoustics->BCChangeTime > 0.0 && mesh->NERFaces)){
    // [EA] LR BCs
    
    string LRDATAFileName;
    newOptions.getArgs("LRVECTFIT", LRDATAFileName);

    FILE * LRDATAFILE = fopen((char*)LRDATAFileName.c_str(),"r");
    if (LRDATAFILE == NULL) {
      printf("Could not find LRDATA file: %s\n", (char*)LRDATAFileName.c_str());
      exit(-1);
    }

    fscanf(LRDATAFILE,"%d %d %d",&acoustics->LRNpoles, &acoustics->LRNRealPoles,&acoustics->LRNImagPoles);
    
    acoustics->LR = (dfloat*) calloc(1 + 2*acoustics->LRNRealPoles + 4*acoustics->LRNImagPoles,sizeof(dfloat));
    acoustics->LRInfo = (dlong*) calloc(3,sizeof(dlong));

    acoustics->LRInfo[0] = acoustics->LRNpoles;
    acoustics->LRInfo[1] = acoustics->LRNRealPoles;
    acoustics->LRInfo[2] = acoustics->LRNImagPoles;

    
    dlong iReadOffset = 0;
    for(int iRead = 0; iRead < acoustics->LRNRealPoles; iRead++){
      fscanf(LRDATAFILE,"%lf",&acoustics->LR[iReadOffset+iRead]);
    }
    iReadOffset += acoustics->LRNRealPoles;
    for(int iRead = 0; iRead < acoustics->LRNImagPoles; iRead++){
      fscanf(LRDATAFILE,"%lf",&acoustics->LR[iReadOffset+iRead]);
    }
    iReadOffset += acoustics->LRNImagPoles;
    for(int iRead = 0; iRead < acoustics->LRNImagPoles; iRead++){
      fscanf(LRDATAFILE,"%lf",&acoustics->LR[iReadOffset+iRead]);
    }
    iReadOffset += acoustics->LRNImagPoles;
    for(int iRead = 0; iRead < acoustics->LRNRealPoles; iRead++){
      fscanf(LRDATAFILE,"%lf",&acoustics->LR[iReadOffset+iRead]);
    }

    iReadOffset += acoustics->LRNRealPoles;
    for(int iRead = 0; iRead < acoustics->LRNImagPoles; iRead++){
      fscanf(LRDATAFILE,"%lf",&acoustics->LR[iReadOffset+iRead]);
    }
    iReadOffset += acoustics->LRNImagPoles;
    for(int iRead = 0; iRead < acoustics->LRNImagPoles; iRead++){
      fscanf(LRDATAFILE,"%lf",&acoustics->LR[iReadOffset+iRead]);
    }
    iReadOffset += acoustics->LRNImagPoles;
    fscanf(LRDATAFILE,"%lf",&acoustics->LR[iReadOffset]);
    fclose(LRDATAFILE);
    acoustics->o_LR = mesh->device.malloc((1 + 2*acoustics->LRNRealPoles + 4*acoustics->LRNImagPoles)*sizeof(dfloat), acoustics->LR);
    acoustics->o_LRInfo = mesh->device.malloc(3*sizeof(dlong), acoustics->LRInfo);
    
  } else {
      // [EA] To avoid empty pointers as this causes crashes on gpu

      acoustics->LRNpoles = 1;
      acoustics->LRNRealPoles = 0;
      acoustics->LRNImagPoles = 0;


      acoustics->LRInfo = (dlong*) calloc(3,sizeof(dlong));
      acoustics->LRInfo[0] = 0;
      acoustics->LRInfo[1] = 0;
      acoustics->LRInfo[2] = 0;
      acoustics->o_LR = mesh->device.malloc(1*sizeof(dfloat));
      acoustics->o_LRInfo = mesh->device.malloc(3*sizeof(dlong), acoustics->LRInfo);

  }

  mesh->NboundaryPointsLocal = mesh->Nfp*mesh->NboundaryFacesLocal;
  mesh->NLRPoints = mesh->Nfp*mesh->NLRFaces;
  mesh->NERPoints = mesh->Nfp*mesh->NERFaces;


  dlong RPointsAlloc = mesh->NLRPoints+mesh->NERPoints ? mesh->NLRPoints+mesh->NERPoints:1; // [EA] Occa cannot have pointers to empty arrays
  dlong NERPointsAlloc = mesh->NERPoints ? mesh->NERPoints:1;
  mesh->mapAccToQ = (dlong*) calloc(RPointsAlloc, sizeof(dlong));
  mesh->mapAccToXYZ = (dlong*) calloc(NERPointsAlloc, sizeof(dlong));
  mesh->mapAccToN = (dlong*) calloc(NERPointsAlloc, sizeof(dlong));

  hlong counter = 0;
  hlong counter2 = 0;
  for(hlong e = 0; e < mesh->Nelements;e++){
    for(hlong n = 0; n < mesh->Nfp*mesh->Nfaces; n++){
      hlong id  = e*mesh->Nfp*mesh->Nfaces + n;
      if(mesh->mapAcc[id] == -2){ // The point id is a domain boundary point
        int face = n/mesh->Nfp;
        int bc = mesh->EToB[face+mesh->Nfaces*e];
        if(bc == 3){
          mesh->mapAcc[id] = counter;


          dlong idM = mesh->vmapM[id];
          int vidM = idM%mesh->Np;
          dlong qbaseM = e*mesh->Np*mesh->Nfields + vidM;
          mesh->mapAccToQ[counter] = qbaseM;
          counter++;
        }
        if(bc == 4){
          mesh->mapAcc[id] = counter2;

          dlong idM = mesh->vmapM[id];
          int vidM = idM%mesh->Np;
          dlong qbaseM = e*mesh->Np*mesh->Nfields + vidM;
          mesh->mapAccToQ[mesh->NLRPoints+counter2] = qbaseM;
          mesh->mapAccToXYZ[counter2] = e*mesh->Np + vidM;
          mesh->mapAccToN[counter2] = mesh->Nsgeo*(e*mesh->Nfaces+face);
          counter2++;
        }
      }
    }
  }
  #if 0
  printf("---MAPACC---\n");
  for(int i1 = 0; i1 < mesh->Nfp*mesh->Nfaces*mesh->Nelements; i1++){
    printf("%d\n",mesh->mapAcc[i1]);
  }
  printf("---MAPACC---\n");
  printf("---MAPACCQ---\n");
  for(int i1 = 0; i1 < RPointsAlloc; i1++){
    printf("%d\n",mesh->mapAccToQ[i1]);
  }
  printf("---MAPACCQ---\n");
  #endif

  // [EA] mapAcc to device
  mesh->o_mapAcc =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(dlong),
                        mesh->mapAcc);
  mesh->o_mapAccToQ =
    mesh->device.malloc(RPointsAlloc*sizeof(dlong),
                        mesh->mapAccToQ);
  mesh->o_mapAccToXYZ =
    mesh->device.malloc(NERPointsAlloc*sizeof(dlong),
                        mesh->mapAccToXYZ);
  mesh->o_mapAccToN =
    mesh->device.malloc(NERPointsAlloc*sizeof(dlong),
                        mesh->mapAccToN);

  // [EA] Moved these from below
  kernelInfo["defines/" "p_blockSize"]= blockSize;
  kernelInfo["defines/" "p_Nfields"]= mesh->Nfields; 

  
  // [EA] Build interpolation operators for ER wave-splitting points
  
  dfloat dx = 0.01; // [EA] Distance to move wave-splitting points
  kernelInfo["defines/" "p_ERdx"]= dx; 

  acoustics->NERPointsTotal = 0;
  MPI_Allreduce(&mesh->NERPoints, &acoustics->NERPointsTotal, 1, MPI_DLONG, MPI_SUM, mesh->comm);


  // [EA] BCChangeTime error checking
  if(acoustics->BCChangeTime > 0.0){
    if(acoustics->NERPointsTotal == 0){
      printf("BCCHANGETIME is greater than 0 with no Extended Reaction boundaries!\n");
      exit(-1);
    }
    if(acoustics->LRInfo[0] != acoustics->ERInfo[0]){
      printf("BCCHANGETIME is greater than 0, but number of poles in Local Reaction and Extended Reaction are not equal!\n");
      exit(-1);
    }
  }


  if(acoustics->NERPointsTotal){
    occa::kernel ERInterpolationOperators = 
              mesh->device.buildKernel(DACOUSTICS "/okl/acousticsERKernel.okl",
				      "ERInterpolationOperators",
				      kernelInfo);

    occa::memory o_invVB, o_EX, o_EY, o_EZ;

    o_EX = mesh->device.malloc(mesh->Nverts*mesh->Nelements*sizeof(dfloat),mesh->EX);
    o_EY = mesh->device.malloc(mesh->Nverts*mesh->Nelements*sizeof(dfloat),mesh->EY);
    o_EZ = mesh->device.malloc(mesh->Nverts*mesh->Nelements*sizeof(dfloat),mesh->EZ);
    o_invVB = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat),mesh->invVB);

    acoustics->o_ERintpol = mesh->device.malloc(mesh->Np*2*mesh->NERPoints*sizeof(dfloat));
    acoustics->o_ERintpolElements = mesh->device.malloc(2*mesh->NERPoints*sizeof(dlong));
    if(mesh->NERPoints){
      ERInterpolationOperators(mesh->NERPoints,
                                mesh->Nelements,
                                o_EX,
                                o_EY,
                                o_EZ,
                                acoustics->o_ERintpol,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_z,
                                mesh->o_sgeo,
                                mesh->o_mapAccToXYZ,
                                mesh->o_mapAccToN,
                                dx,
                                o_invVB,
                                acoustics->o_ERintpolElements);
    }
    
    #if 1
    // [EA] Assume no communication needed, check in next if statement.
    acoustics->NERComPoints = 0;
    acoustics->NComPointsToSendAllRanks = 0;
    if(mesh->size > 1){
      // [EA] Find wave-splitting points belonging to other cores

      // Copy ERintpolElements to host and find the number of wave-splitting points belonging to other cores
      dlong *ERintpolElements;
      ERintpolElements = (dlong*) calloc(2*mesh->NERPoints, sizeof(dlong));    
      if(mesh->NERPoints){
        acoustics->o_ERintpolElements.copyTo(ERintpolElements);
      }

      for(dlong itt = 0; itt < 2*mesh->NERPoints; itt++){
        if(ERintpolElements[itt] < 0){
          acoustics->NERComPoints++;
        }
      }
      acoustics->ERComPointsIdx = (dlong*) calloc(acoustics->NERComPoints, sizeof(dlong)); // Stores the idx of the ws points on other ranks
      acoustics->ERComPoints = (dfloat*) calloc(acoustics->NERComPoints*3, sizeof(dfloat)); // Stores xyz for ws points on other ranks

      // Save xyz of the missing wave-splitting points
      dlong counteritt = 0;
      for(dlong itt = 0; itt < 2*mesh->NERPoints; itt++){
        if(ERintpolElements[itt] < 0){
          acoustics->ERComPointsIdx[counteritt] = itt;
          
          dlong CPXYZidx = mesh->mapAccToXYZ[itt/2];
          dlong CPNidx = mesh->mapAccToN[itt/2];

          dfloat nxi = mesh->sgeo[CPNidx+NXID];
          dfloat nyi = mesh->sgeo[CPNidx+NYID];
          dfloat nzi = mesh->sgeo[CPNidx+NZID];

          dfloat xi = mesh->x[CPXYZidx];
          dfloat yi = mesh->y[CPXYZidx];
          dfloat zi = mesh->z[CPXYZidx];

          // ERintpolElements[itt] will either be -1 or -2 based on which WS of the two points are missing
          xi += ERintpolElements[itt]*dx*nxi; 
          yi += ERintpolElements[itt]*dx*nyi;
          zi += ERintpolElements[itt]*dx*nzi;

          acoustics->ERComPoints[counteritt*3+0] = xi;
          acoustics->ERComPoints[counteritt*3+1] = yi;
          acoustics->ERComPoints[counteritt*3+2] = zi;

          counteritt++;
        }
      }
      
      // Find total amount of communication points
      dlong NERComPointsAllRanks = 0; // Total number of communication points spread across all ranks.
      dlong *recvCounts, *recvCountsCum; // Array of NERComPoints*3 (xyz) for all ranks, and its cummulative, used for MPI comm
      recvCounts = (dlong*) calloc(mesh->size, sizeof(dlong));
      recvCountsCum = (dlong*) calloc(mesh->size, sizeof(dlong));
      MPI_Allgather(&acoustics->NERComPoints, 1, MPI_DLONG, recvCounts, 1, MPI_DLONG, mesh->comm);

      for(dlong itt = 0; itt < mesh->size; itt++){
        NERComPointsAllRanks += recvCounts[itt];
        recvCounts[itt] *= 3;
        if(itt > 0){
          recvCountsCum[itt] = recvCountsCum[itt-1] + recvCounts[itt-1];
        }
      }

      if(NERComPointsAllRanks){
        // Send unfound wave-splitting points to all ranks
        dfloat *ERComPointsAllRanks; // Contains missing xyz for missing ws points for all ranks. In order of ranks.
        ERComPointsAllRanks = (dfloat*) calloc(NERComPointsAllRanks*3, sizeof(dfloat));
        MPI_Allgatherv(acoustics->ERComPoints, acoustics->NERComPoints*3, MPI_DFLOAT,
        ERComPointsAllRanks, recvCounts, recvCountsCum, MPI_DFLOAT, mesh->comm);
        
        // Find element of previosly unfound WS points.
        acoustics->ERintpolElementsCom = (dlong*) calloc(NERComPointsAllRanks, sizeof(dlong));
        acoustics->o_ERintpolElementsCom = 
              mesh->device.malloc(NERComPointsAllRanks*sizeof(dlong),acoustics->ERintpolElementsCom);
        occa::memory o_recvCountsCum, o_recvCounts, o_ERComPointsAllRanks;
        o_recvCountsCum = 
              mesh->device.malloc(mesh->size*sizeof(dlong),recvCountsCum);
        o_recvCounts = 
              mesh->device.malloc(mesh->size*sizeof(dlong),recvCounts);    
        o_ERComPointsAllRanks = 
              mesh->device.malloc(NERComPointsAllRanks*3*sizeof(dfloat),ERComPointsAllRanks);  
        occa::kernel ERFindElementsCom = 
              mesh->device.buildKernel(DACOUSTICS "/okl/acousticsERKernel.okl",
              "ERFindElementsCom",
              kernelInfo);
        ERFindElementsCom(NERComPointsAllRanks,o_recvCountsCum,o_recvCounts,mesh->rank,o_ERComPointsAllRanks,
                          o_EX,o_EY,o_EZ,mesh->Nelements,acoustics->o_ERintpolElementsCom);


        // Count how many points current rank has to send to each other rank
        dlong NComPointsToSendAllRanks = 0;
        dlong *NComPointsToSendToRanks;
        NComPointsToSendToRanks = (dlong*) calloc(mesh->size, sizeof(dlong));
        acoustics->o_ERintpolElementsCom.copyTo(acoustics->ERintpolElementsCom);
        dlong currRank = 0;
        dlong currRankCum = recvCounts[currRank]/3;
        for(dlong itt = 0; itt < NERComPointsAllRanks; itt++){
          while(itt >= currRankCum){
              currRankCum += recvCounts[currRank+1]/3;
              currRank++;
              if(currRank > mesh->size){
                printf("Infinite while-loop\n");
                exit(-1);
              }
          }
          if(acoustics->ERintpolElementsCom[itt] != -1){
            NComPointsToSendToRanks[currRank]++;
            NComPointsToSendAllRanks++;
          }
        }
        dlong *NComPointsToSendToRanksCum;
        NComPointsToSendToRanksCum = (dlong*) calloc(mesh->size, sizeof(dlong));
        for(int itt = 1; itt < mesh->size; itt++){
          NComPointsToSendToRanksCum[itt] = NComPointsToSendToRanksCum[itt-1] + NComPointsToSendToRanks[itt-1];
        }

        // Find idx of WS points to send to each rank
        dlong *comPointsIdx, *comPointsToSend;
        comPointsIdx = (dlong*) calloc(NComPointsToSendAllRanks, sizeof(dlong));
        comPointsToSend = (dlong*) calloc(NComPointsToSendAllRanks, sizeof(dlong));

        currRank = 0;
        currRankCum = recvCounts[currRank]/3;;
        dlong currRankCounter = 0;
        dlong counterk = 0;
        for(dlong itt = 0; itt < NERComPointsAllRanks; itt++){
          while(itt >= currRankCum){
              currRankCum += recvCounts[currRank+1]/3;
              currRank++;
              currRankCounter = 0;
              if(currRank > mesh->size){
                printf("Infinite while-loop\n");
                exit(-1);
              }
          }
          if(acoustics->ERintpolElementsCom[itt] != -1){
            comPointsIdx[counterk] = currRankCounter;
            comPointsToSend[counterk] = itt;
            counterk++;
            
          }
          currRankCounter++;
        }

        acoustics->recvCountsArray = (dlong*) calloc(mesh->size*mesh->size, sizeof(dlong));
        MPI_Allgather(NComPointsToSendToRanks, mesh->size, MPI_DLONG, acoustics->recvCountsArray, mesh->size, MPI_DLONG, mesh->comm);

        dlong NComPointsIdx = 0;
        dlong *NComPointsIdxPerRank;
        dlong *NComPointsIdxPerRankCum;
        NComPointsIdxPerRank = (dlong*) calloc(mesh->size, sizeof(dlong));
        NComPointsIdxPerRankCum = (dlong*) calloc(mesh->size, sizeof(dlong));
        for(int itt = 0; itt < mesh->size; itt++){
          NComPointsIdx += acoustics->recvCountsArray[mesh->rank + mesh->size*itt]; // Total amount of WS points to recv from all other ranks
          NComPointsIdxPerRank[itt] = acoustics->recvCountsArray[mesh->rank + mesh->size*itt]; // Array of number of WS points to recveive from rank itt
          if(itt > 0){
            NComPointsIdxPerRankCum[itt] = NComPointsIdxPerRankCum[itt-1] + acoustics->recvCountsArray[mesh->rank + mesh->size*(itt-1)];
          }
        }

        if(acoustics->NERComPoints != NComPointsIdx){
          printf("Unable to find all Extended Reaction Wave-splitting points. Check if mesh is correct.\n");
          exit(-1);
        }

        dlong *comPointsIdxAll;
        comPointsIdxAll = (dlong*) calloc(acoustics->NERComPoints, sizeof(dlong));

        MPI_Alltoallv(comPointsIdx, NComPointsToSendToRanks, NComPointsToSendToRanksCum, MPI_DLONG, 
                      comPointsIdxAll, NComPointsIdxPerRank, NComPointsIdxPerRankCum, MPI_DLONG, mesh->comm);

        
        acoustics->o_ERComPointsIdx = 
              mesh->device.malloc(acoustics->NERComPoints*sizeof(dlong),acoustics->ERComPointsIdx);
        
        acoustics->o_ERintpolCom =
              mesh->device.malloc(mesh->Np*NComPointsToSendAllRanks*sizeof(dfloat));

        acoustics->o_recvCountsArray = 
              mesh->device.malloc(mesh->size*mesh->size*sizeof(dlong),acoustics->recvCountsArray);


        acoustics->o_comPointsIdxAll = 
              mesh->device.malloc(acoustics->NERComPoints*sizeof(dlong),comPointsIdxAll);

        acoustics->o_comPointsToSend = 
              mesh->device.malloc(NComPointsToSendAllRanks*sizeof(dlong),comPointsToSend);
        
        acoustics->NComPointsToSendAllRanks = NComPointsToSendAllRanks;

        occa::kernel ERmakeInterpolatorCom;
        ERmakeInterpolatorCom = 
            mesh->device.buildKernel(DACOUSTICS "/okl/acousticsERKernel.okl",
                  "ERmakeInterpolatorCom",
                  kernelInfo);
        if(NComPointsToSendAllRanks){
          ERmakeInterpolatorCom(NComPointsToSendAllRanks,
                                acoustics->o_comPointsToSend,
                                o_ERComPointsAllRanks,
                                acoustics->o_ERintpolElementsCom,
                                o_EX,
                                o_EY,
                                o_EZ,
                                acoustics->o_ERintpolCom,
                                o_invVB);
        }
        acoustics->vtSend = (dfloat*) calloc(NComPointsToSendAllRanks*3, sizeof(dfloat));
        acoustics->vtRecv = (dfloat*) calloc(acoustics->NERComPoints*3, sizeof(dfloat));
        acoustics->o_vtSend = 
              mesh->device.malloc(NComPointsToSendAllRanks*3*sizeof(dfloat),acoustics->vtSend);
        acoustics->o_vtRecv = 
              mesh->device.malloc(acoustics->NERComPoints*3*sizeof(dfloat),acoustics->vtRecv);
    
    free(ERintpolElements);
    free(recvCounts);
    free(recvCountsCum);
    free(ERComPointsAllRanks);
    free(NComPointsToSendToRanks);
    free(NComPointsToSendToRanksCum);
    free(comPointsIdx);
    free(comPointsToSend);
    free(NComPointsIdxPerRank);
    free(NComPointsIdxPerRankCum);
    free(comPointsIdxAll);
    // [DEBUG] PRINT EVERYTHING Delete!
    #if 0
    for(int itt2 = 0; itt2 < mesh->size; itt2++){
      
      if(itt2 == mesh->rank){
        printf("-------------- r = %d --------------\n",mesh->rank);
    printf("ERCOMPOINTSIDX:\n");
    for(int itt = 0; itt < acoustics->NERComPoints;itt++){
      printf("%d\n",acoustics->ERComPointsIdx[itt]);
    }
    printf("---------------------------------------------\n");

    printf("ERintpolCom:\n");
    dfloat *t123;
    t123 = (dfloat*) calloc(mesh->Np*NComPointsToSendAllRanks, sizeof(dfloat));
    acoustics->o_ERintpolCom.copyTo(t123);
    for(int itt = 0; itt < mesh->Np*NComPointsToSendAllRanks;itt++){
      printf("%.15f\n",t123[itt]);
    }
    printf("---------------------------------------------\n");

    printf("recvCountsArray:\n");
    for(int itt = 0; itt < mesh->size*mesh->size;itt++){
      printf("%d\n",acoustics->recvCountsArray[itt]);
    }
    printf("---------------------------------------------\n"); 

   printf("comPointsIdxAll:\n");
    for(int itt = 0; itt < acoustics->NERComPoints;itt++){
      printf("%d\n",comPointsIdxAll[itt]);
    }
    printf("---------------------------------------------\n"); 

    printf("comPointsToSend:\n");
    for(int itt = 0; itt < NComPointsToSendAllRanks;itt++){
      printf("%d\n",comPointsToSend[itt]);
    }
    printf("---------------------------------------------\n"); 
  printf("-------------- r = %d --------------\n",mesh->rank);
      }
      MPI_Barrier(mesh->comm);
    }

    #endif

      }
      
    }
   #endif

    dfloat *tempvt, *tempvi; 
    dlong *tempanglei;
    tempvt = (dfloat*) calloc(mesh->NERPoints*4*3*3, sizeof(dfloat));
    tempvi = (dfloat*) calloc(mesh->NERPoints*3, sizeof(dfloat));
    tempanglei = (dlong*) calloc(mesh->NERPoints, sizeof(dlong));
    
    acoustics->o_vt = mesh->device.malloc(mesh->NERPoints*4*3*3*sizeof(dfloat),tempvt);
    acoustics->o_vi = mesh->device.malloc(mesh->NERPoints*3*sizeof(dfloat),tempvi);
    acoustics->o_anglei = mesh->device.malloc(mesh->NERPoints*sizeof(dlong),tempanglei);

    free(tempvt);
    free(tempvi);
    free(tempanglei);

    // [EA] DEBUG: Print intpol
    #if 0
    dfloat * tempIntPol;
    tempIntPol =(dfloat*) calloc(mesh->Np*2*mesh->NERPoints, sizeof(dfloat));
    acoustics->o_ERintpol.copyTo(tempIntPol);
    for(int ti = 0; ti < 2*mesh->NERPoints;ti++){
      for(int tj = 0; tj < mesh->Np;tj++){
        printf("%f ",tempIntPol[ti*mesh->Np+tj]);
      }
      printf("\n");
    }
    #endif

    o_invVB.free();
    o_EX.free();
    o_EY.free();
    o_EZ.free();
    mesh->o_mapAccToXYZ.free();
    free(mesh->mapAccToXYZ);
  } else {
    // [EA] To avoid empty pointers, if no ER points.
    acoustics->o_ERintpol = mesh->device.malloc(1*sizeof(dfloat));
    acoustics->o_vt = mesh->device.malloc(1*sizeof(dfloat));
    acoustics->o_vi = mesh->device.malloc(1*sizeof(dfloat));
    acoustics->o_anglei = mesh->device.malloc(1*sizeof(dlong));
  }
  
  // [EA] erk and esdirk to device
  mesh->o_erka = 
    mesh->device.malloc(36*sizeof(dfloat), mesh->erka);
  mesh->o_erkb = 
    mesh->device.malloc(6*sizeof(dfloat), mesh->erkb);
  mesh->o_erke = 
    mesh->device.malloc(6*sizeof(dfloat), mesh->erke);
  mesh->o_esdirka = 
    mesh->device.malloc(36*sizeof(dfloat), mesh->esdirka);
  mesh->o_esdirkb = 
    mesh->device.malloc(6*sizeof(dfloat), mesh->esdirkb);
  mesh->o_esdirke = 
    mesh->device.malloc(6*sizeof(dfloat), mesh->esdirke);


  
  acoustics->acc = 
    (dfloat*) calloc(acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1, sizeof(dfloat));
  acoustics->rhsacc = 
    (dfloat*) calloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1), sizeof(dfloat));
  acoustics->resacc = 
    (dfloat*) calloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1), sizeof(dfloat));
  acoustics->o_acc =
    mesh->device.malloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1)*sizeof(dfloat), acoustics->acc);
  acoustics->o_rhsacc =
    mesh->device.malloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1)*sizeof(dfloat), acoustics->rhsacc);
  acoustics->o_resacc =
    mesh->device.malloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1)*sizeof(dfloat), acoustics->resacc);

  if(!mesh->rank){
    cout << "TIME INTEGRATOR (" << newOptions.getArgs("TIME INTEGRATOR") << ")" << endl;
  }
  if (newOptions.compareArgs("TIME INTEGRATOR","LSERK4")){
    acoustics->o_resq =
      mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), acoustics->resq);
  }

  if (newOptions.compareArgs("TIME INTEGRATOR","DOPRI5")){
    printf("setting up DOPRI5\n");
    int NrkStages = 7;
    acoustics->o_rkq =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), acoustics->rkq);
    acoustics->o_rkrhsq =
      mesh->device.malloc(NrkStages*mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), acoustics->rkrhsq);
    acoustics->o_rkerr =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), acoustics->rkerr);
  
    acoustics->o_errtmp = mesh->device.malloc(acoustics->Nblock*sizeof(dfloat), acoustics->errtmp);

    acoustics->o_rkA = mesh->device.malloc(acoustics->Nrk*acoustics->Nrk*sizeof(dfloat), acoustics->rkA);
    acoustics->o_rkE = mesh->device.malloc(  acoustics->Nrk*sizeof(dfloat), acoustics->rkE);
  }

  if (newOptions.compareArgs("TIME INTEGRATOR","EIRK4") || newOptions.compareArgs("TIME INTEGRATOR","EIRK4Adap")){
    acoustics->k1acc = (dfloat*) calloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1), sizeof(dfloat));
    acoustics->k2acc = (dfloat*) calloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1), sizeof(dfloat));
    acoustics->k3acc = (dfloat*) calloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1), sizeof(dfloat));
    acoustics->k4acc = (dfloat*) calloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1), sizeof(dfloat));
    acoustics->k5acc = (dfloat*) calloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1), sizeof(dfloat));
    acoustics->k6acc = (dfloat*) calloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1), sizeof(dfloat));

    acoustics->Xacc = (dfloat*) calloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1), sizeof(dfloat));

    acoustics->k1rhsq = (dfloat*) calloc(mesh->Np*mesh->Nelements*mesh->Nfields, sizeof(dfloat));
    acoustics->k2rhsq = (dfloat*) calloc(mesh->Np*mesh->Nelements*mesh->Nfields, sizeof(dfloat));
    acoustics->k3rhsq = (dfloat*) calloc(mesh->Np*mesh->Nelements*mesh->Nfields, sizeof(dfloat));
    acoustics->k4rhsq = (dfloat*) calloc(mesh->Np*mesh->Nelements*mesh->Nfields, sizeof(dfloat));
    acoustics->k5rhsq = (dfloat*) calloc(mesh->Np*mesh->Nelements*mesh->Nfields, sizeof(dfloat));
    acoustics->k6rhsq = (dfloat*) calloc(mesh->Np*mesh->Nelements*mesh->Nfields, sizeof(dfloat));

    acoustics->resq = (dfloat*) calloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields, sizeof(dfloat));

    acoustics->o_k1acc = 
          mesh->device.malloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1)*sizeof(dfloat), acoustics->k1acc);
    acoustics->o_k2acc = 
          mesh->device.malloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1)*sizeof(dfloat), acoustics->k2acc);
    acoustics->o_k3acc = 
          mesh->device.malloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1)*sizeof(dfloat), acoustics->k3acc);
    acoustics->o_k4acc = 
          mesh->device.malloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1)*sizeof(dfloat), acoustics->k4acc);
    acoustics->o_k5acc = 
          mesh->device.malloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1)*sizeof(dfloat), acoustics->k5acc);
    acoustics->o_k6acc = 
          mesh->device.malloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1)*sizeof(dfloat), acoustics->k6acc);

    acoustics->o_Xacc = 
          mesh->device.malloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1)*sizeof(dfloat), acoustics->resacc);

    acoustics->o_k1rhsq = 
          mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), acoustics->k1rhsq);
    acoustics->o_k2rhsq = 
          mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), acoustics->k2rhsq);
    acoustics->o_k3rhsq = 
          mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), acoustics->k3rhsq);
    acoustics->o_k4rhsq = 
          mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), acoustics->k4rhsq);
    acoustics->o_k5rhsq = 
          mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), acoustics->k5rhsq);
    acoustics->o_k6rhsq = 
          mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), acoustics->k6rhsq);

    acoustics->o_resq = 
          mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), acoustics->resq);
          
  }
  if(newOptions.compareArgs("TIME INTEGRATOR","EIRK4Adap")){
    acoustics->rkq = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
				sizeof(dfloat));
    acoustics->rkerr = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
				sizeof(dfloat));

    acoustics->rkAcc = (dfloat*) calloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1),
				sizeof(dfloat));

    acoustics->rkerrAcc = (dfloat*) calloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1),
				sizeof(dfloat));

    acoustics->o_rkq =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), acoustics->rkq);
    acoustics->o_rkerr =
      mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), acoustics->rkerr);
    acoustics->o_rkAcc = 
          mesh->device.malloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1)*sizeof(dfloat), acoustics->rkAcc);
    acoustics->o_rkerrAcc = 
          mesh->device.malloc((acoustics->LRInfo[0]*mesh->NLRPoints + acoustics->ERInfo[0]*mesh->NERPoints + 1)*sizeof(dfloat), acoustics->rkerrAcc);
  }
  
  if(mesh->totalHaloPairs>0){
    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    mesh->o_haloBuffer =
      mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat));

    // MPI send buffer
    acoustics->haloBytes = mesh->totalHaloPairs*mesh->Np*acoustics->Nfields*sizeof(dfloat);

    acoustics->o_haloBuffer = mesh->device.malloc(acoustics->haloBytes);

#if 0
    occa::memory o_sendBuffer = mesh->device.mappedAlloc(acoustics->haloBytes, NULL);
    occa::memory o_recvBuffer = mesh->device.mappedAlloc(acoustics->haloBytes, NULL);

    acoustics->sendBuffer = (dfloat*) o_sendBuffer.getMappedPointer();
    acoustics->recvBuffer = (dfloat*) o_recvBuffer.getMappedPointer();
#endif

    acoustics->sendBuffer = (dfloat*) occaHostMallocPinned(mesh->device, acoustics->haloBytes, NULL, acoustics->o_sendBuffer);
    acoustics->recvBuffer = (dfloat*) occaHostMallocPinned(mesh->device, acoustics->haloBytes, NULL, acoustics->o_recvBuffer);    
  }

  //  p_RT, p_rbar, p_ubar, p_vbar
  // p_half, p_two, p_third, p_Nstresses
  
  //kernelInfo["defines/" "p_Nfields"]= mesh->Nfields; // [EA] Moved further up
  const dfloat p_one = 1.0, p_two = 2.0, p_half = 1./2., p_third = 1./3., p_zero = 0;

  kernelInfo["defines/" "p_two"]= p_two;
  kernelInfo["defines/" "p_one"]= p_one;
  kernelInfo["defines/" "p_half"]= p_half;
  kernelInfo["defines/" "p_third"]= p_third;
  kernelInfo["defines/" "p_zero"]= p_zero;
  
  
  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int NblockV = 1024/mesh->Np; // works for CUDA
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  int cubMaxNodes = mymax(mesh->Np, (mesh->intNfp*mesh->Nfaces));
  kernelInfo["defines/" "p_cubMaxNodes"]= cubMaxNodes;
  int cubMaxNodes1 = mymax(mesh->Np, (mesh->intNfp));
  kernelInfo["defines/" "p_cubMaxNodes1"]= cubMaxNodes1;

  kernelInfo["defines/" "p_Lambda2"]= 0.5f;

  kernelInfo["defines/" "p_abstol"]= 1.0e-5;
  kernelInfo["defines/" "p_reltol"]= 1.0e-5;

  //kernelInfo["defines/" "p_blockSize"]= blockSize; // [EA] Moved further up

  // [EA] LR offset defines
  kernelInfo["defines/" "p_LRA"]= 0;
  kernelInfo["defines/" "p_LRB"]= acoustics->LRInfo[1];
  kernelInfo["defines/" "p_LRC"]= acoustics->LRInfo[1] + acoustics->LRInfo[2];
  kernelInfo["defines/" "p_LRLambda"]= acoustics->LRInfo[1] + 2*acoustics->LRInfo[2];
  kernelInfo["defines/" "p_LRAlpha"]= 2*acoustics->LRInfo[1] + 2*acoustics->LRInfo[2];
  kernelInfo["defines/" "p_LRBeta"]= 2*acoustics->LRInfo[1] + 3*acoustics->LRInfo[2];
  kernelInfo["defines/" "p_LRYinf"]= 2*acoustics->LRInfo[1] + 4*acoustics->LRInfo[2];
  
  // [EA] ER offset defines
  kernelInfo["defines/" "p_ERA"]= 0;
  kernelInfo["defines/" "p_ERB"]= acoustics->ERInfo[1]*Nangles;
  kernelInfo["defines/" "p_ERC"]= acoustics->ERInfo[1]*Nangles + acoustics->ERInfo[2]*Nangles;
  kernelInfo["defines/" "p_ERLambda"]= acoustics->ERInfo[1]*Nangles + 2*acoustics->ERInfo[2]*Nangles;
  kernelInfo["defines/" "p_ERAlpha"]= acoustics->ERInfo[1]+acoustics->ERInfo[1]*Nangles + 2*acoustics->ERInfo[2]*Nangles;
  kernelInfo["defines/" "p_ERBeta"]= acoustics->ERInfo[1]+acoustics->ERInfo[1]*Nangles + acoustics->ERInfo[2]+2*acoustics->ERInfo[2]*Nangles;
  kernelInfo["defines/" "p_ERYinf"]= acoustics->ERInfo[1]+acoustics->ERInfo[1]*Nangles + 2*acoustics->ERInfo[2]+2*acoustics->ERInfo[2]*Nangles;





  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  // set kernel name suffix
  char *suffix;
  
  if(acoustics->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(acoustics->elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(acoustics->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(acoustics->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  // kernels from volume file
  sprintf(fileName, DACOUSTICS "/okl/acousticsVolume%s.okl", suffix);
  sprintf(kernelName, "acousticsVolume%s", suffix);

  //printf("fileName=[ %s ] \n", fileName);
  //printf("kernelName=[ %s ]\n", kernelName);
  
  acoustics->volumeKernel =  mesh->device.buildKernel(fileName, kernelName, kernelInfo);

  // kernels from surface file
  sprintf(fileName, DACOUSTICS "/okl/acousticsSurface%s.okl", suffix);
  sprintf(kernelName, "acousticsSurface%s", suffix);
  
  acoustics->surfaceKernel = mesh->device.buildKernel(fileName, kernelName, kernelInfo);

  // kernels from update file
  acoustics->updateKernel =
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsUpdate.okl",
				       "acousticsUpdate",
				       kernelInfo);

  // [EA] LR update kernel
  acoustics->updateKernelLR = 
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsUpdate.okl",
				       "acousticsUpdateLRAcc",
				       kernelInfo);
  acoustics->updateKernelER = 
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsUpdate.okl",
				       "acousticsUpdateERAcc",
				       kernelInfo);
  acoustics->acousticsUpdateEIRK4 = 
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsUpdate.okl",
				       "acousticsUpdateEIRK4",
				       kernelInfo);
  acoustics->acousticsUpdateEIRK4AccLR = 
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsUpdate.okl",
				       "acousticsUpdateEIRK4AccLR",
				       kernelInfo);
  acoustics->acousticsUpdateEIRK4AccER = 
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsUpdate.okl",
				       "acousticsUpdateEIRK4AccER",
				       kernelInfo);

  acoustics->acousticsWSComInterpolation = 
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsERKernel.okl",
				       "acousticsWSComInterpolation",
				       kernelInfo);

  acoustics->acousticsReceiverInterpolation = 
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsReceiverKernel.okl",
				       "acousticsReceiverInterpolation",
				       kernelInfo);
  
  acoustics->ERInsertComVT = 
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsERKernel.okl",
				       "ERInsertComVT",
				       kernelInfo);
  
  acoustics->rkUpdateKernel =
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsUpdate.okl",
				       "acousticsRkUpdate",
				       kernelInfo);
  acoustics->rkStageKernel =
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsUpdate.okl",
				       "acousticsRkStage",
				       kernelInfo);

  acoustics->rkErrorEstimateKernel =
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsUpdate.okl",
				       "acousticsErrorEstimate",
				       kernelInfo);

  // [EA] Copy from q to qRecv (Currently not in use, changed to using copyTo instead)
  acoustics->receiverKernel =
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsReceiverKernel.okl",
              "acousticsReceiverKernel",
              kernelInfo);

  acoustics->ERMoveVT = 
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsERKernel.okl",
              "ERMoveVT",
              kernelInfo);


  // [EA] Detect angle of incoming wave for Extended Reaction boundary condition
  acoustics->ERangleDetection =
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsERKernel.okl",
                "ERangleDetection",
                kernelInfo);
  
  acoustics->acousticsErrorEIRK4 = 
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsUpdate.okl",
                "acousticsErrorEIRK4",
                kernelInfo);

  acoustics->acousticsErrorEIRK4r = 
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsUpdate.okl",
                "acousticsErrorEIRK4r",
                kernelInfo);

  acoustics->acousticsErrorEIRK4Acc = 
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsUpdate.okl",
                "acousticsErrorEIRK4Acc",
                kernelInfo);

  acoustics->acousticsErrorEIRK4Accr = 
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsUpdate.okl",
                "acousticsErrorEIRK4Accr",
                kernelInfo);

  // fix this later
  mesh->haloExtractKernel =
    mesh->device.buildKernel(DHOLMES "/okl/meshHaloExtract3D.okl",
				       "meshHaloExtract3D",
				       kernelInfo);


  #if 0
  // MAX TEST!
  printf("hejsa1\n");
  occa::kernel maxArray = 
          mesh->device.buildKernel(DACOUSTICS "/okl/acousticsERKernel.okl",
		      "maxarray",
		      kernelInfo);

printf("hejsa2\n");
  dlong goddagLen = 100000000;
  dlong *goddag;
  goddag = (dlong*) calloc(goddagLen, sizeof(dlong));
  printf("hejsa3\n");
  for(int itt = 0; itt < goddagLen; itt++){
    goddag[itt] = itt;
  }
  printf("hejsa4\n");
  goddag[70] = 123123123;
  
  occa::memory o_goddag;
  o_goddag = mesh->device.malloc(goddagLen*sizeof(dlong),goddag);

  double start, end;
  start = MPI_Wtime();
  dlong *out;
  dlong outLen = (goddagLen+blockSize-1)/blockSize;
  out = (dlong*) calloc(outLen, sizeof(dlong));
  occa::memory o_out = mesh->device.malloc(outLen*sizeof(dlong),out);

  maxArray(goddagLen, o_goddag, o_out);

  o_out.copyTo(out);

  dlong outmax = -5000;
  for(int itt = 0; itt < outLen; itt++){
    outmax = (outmax < out[itt]) ? out[itt] : outmax;
  }
  printf("time taken for max gpu : %f, outmax = %d\n", MPI_Wtime() - start,outmax);

  start = MPI_Wtime();
  o_goddag.copyTo(goddag);
  outmax = -5000;
  for(int itt = 0; itt < goddagLen; itt++){
    outmax = (outmax < goddag[itt]) ? goddag[itt] : outmax;
  }
  printf("time taken for max cpu : %f, outmax = %d\n", MPI_Wtime() - start,outmax);

  printf("outLen = %d\n",outLen);
  printf("MAX SKAL VÆRE 123123123:  DEN ER %d\n",outmax);

  #endif


  

  // Print setup parameters to out file
  if(!mesh->rank){
    string mshFileName;
    newOptions.getArgs("MESH FILE", mshFileName);
    string timeInt;
    newOptions.getArgs("TIME INTEGRATOR", timeInt);
    string recvFile;
    newOptions.getArgs("RECEIVER", recvFile);
    string LRFile;
    newOptions.getArgs("LRVECTFIT", LRFile);
    string ERFile;
    newOptions.getArgs("ERVECTFIT", ERFile);
    printf("-----SETUP PARAMETERS START-----\n");
    printf("N = %d\n",mesh->N);
    printf("Mesh file: %s\n",(char*)mshFileName.c_str());
    printf("Final Time: %g\n",mesh->finalTime);
    printf("CFL = %g\n",cfl);
    printf("dt = %g\n",mesh->dt);
    printf("Time integrator: %s\n",(char*)timeInt.c_str());
    printf("Receiver file: %s\n",(char*)recvFile.c_str());
    printf("Local Reaction file: %s\n",(char*)LRFile.c_str());
    printf("Extended Reaction file: %s\n",(char*)ERFile.c_str());
    printf("BCCHANGETIME = %g\n",acoustics->BCChangeTime);
    printf("Boundary conditions on surface indices from .msh file:\n");
    printf("Rigid = [ ");
    for(int jj = 1; jj < 1000; jj+=2){ // Hardcoded for 500 surfaces, see meshParallelReaderTet3D.c
      if(1 == mesh->mshPrint[jj]){
        printf("%d ",mesh->mshPrint[jj-1]);
      }
    }
    printf("]\nFrequency Independent = [ ");
    for(int jj = 1; jj < 1000; jj+=2){ // Hardcoded for 500 surfaces, see meshParallelReaderTet3D.c
      if(2 == mesh->mshPrint[jj]){
        printf("%d ",mesh->mshPrint[jj-1]);
      }
    }
    printf("]\nLocal Reaction = [ ");
    for(int jj = 1; jj < 1000; jj+=2){ // Hardcoded for 500 surfaces, see meshParallelReaderTet3D.c
      if(3 == mesh->mshPrint[jj]){
        printf("%d ",mesh->mshPrint[jj-1]);
      }
    }
    printf("]\nExtended Reaction = [ ");
    for(int jj = 1; jj < 1000; jj+=2){ // Hardcoded for 500 surfaces, see meshParallelReaderTet3D.c
      if(4 == mesh->mshPrint[jj]){
        printf("%d ",mesh->mshPrint[jj-1]);
      }
    }
  printf("]\n");
  printf("-----SETUP PARAMETERS END-----\n");
  }
  free(mesh->mshPrint);  
    
  

  
  return acoustics;
}
