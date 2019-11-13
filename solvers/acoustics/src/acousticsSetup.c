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

  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat t = 0;
      dfloat x = mesh->x[n + mesh->Np*e];
      dfloat y = mesh->y[n + mesh->Np*e];
      dfloat z = mesh->z[n + mesh->Np*e];

      dlong qbase = e*mesh->Np*mesh->Nfields + n;

      dfloat u = 0, v = 0, w = 0, r = 0;
      
      acousticsGaussianPulse(x, y, z, t, &r, &u, &v, &w);
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

  // [EA] Read and allocate space for LR accumulators
  string LRDATAFileName;
  newOptions.getArgs("VECTFIT", LRDATAFileName);

  FILE * LRDATAFILE = fopen((char*)LRDATAFileName.c_str(),"r");
  if (LRDATAFILE == NULL) {
    printf("Could not find LRDATA file: %s\n", (char*)LRDATAFileName.c_str());
    exit(-1);
  }

  fscanf(LRDATAFILE,"%d %d %d",&acoustics->Npoles, &acoustics->NRealPoles,&acoustics->NImagPoles);
  acoustics->LRA = (dfloat*) calloc(acoustics->NRealPoles, sizeof(dfloat));
  acoustics->LRLambda = (dfloat*) calloc(acoustics->NRealPoles, sizeof(dfloat));

  acoustics->LRB = (dfloat*) calloc(acoustics->NImagPoles, sizeof(dfloat));
  acoustics->LRC = (dfloat*) calloc(acoustics->NImagPoles, sizeof(dfloat));
  acoustics->LRAlpha = (dfloat*) calloc(acoustics->NImagPoles, sizeof(dfloat));
  acoustics->LRBeta = (dfloat*) calloc(acoustics->NImagPoles, sizeof(dfloat));

  for(int iRead = 0; iRead < acoustics->NRealPoles; iRead++){
    fscanf(LRDATAFILE,"%lf",&acoustics->LRA[iRead]);
  }
  for(int iRead = 0; iRead < acoustics->NImagPoles; iRead++){
    fscanf(LRDATAFILE,"%lf",&acoustics->LRB[iRead]);
  }
  for(int iRead = 0; iRead < acoustics->NImagPoles; iRead++){
    fscanf(LRDATAFILE,"%lf",&acoustics->LRC[iRead]);
  }
  for(int iRead = 0; iRead < acoustics->NRealPoles; iRead++){
    fscanf(LRDATAFILE,"%lf",&acoustics->LRLambda[iRead]);
  }
  for(int iRead = 0; iRead < acoustics->NImagPoles; iRead++){
    fscanf(LRDATAFILE,"%lf",&acoustics->LRAlpha[iRead]);
  }
  for(int iRead = 0; iRead < acoustics->NImagPoles; iRead++){
    fscanf(LRDATAFILE,"%lf",&acoustics->LRBeta[iRead]);
  }
  fscanf(LRDATAFILE,"%lf",&acoustics->LRYinf);
  fclose(LRDATAFILE);

  acoustics->o_LRA = mesh->device.malloc(acoustics->NRealPoles*sizeof(dfloat), acoustics->LRA);
  acoustics->o_LRB = mesh->device.malloc(acoustics->NImagPoles*sizeof(dfloat), acoustics->LRB);
  acoustics->o_LRC = mesh->device.malloc(acoustics->NImagPoles*sizeof(dfloat), acoustics->LRC);
  acoustics->o_LRLambda = mesh->device.malloc(acoustics->NRealPoles*sizeof(dfloat), acoustics->LRLambda);
  acoustics->o_LRAlpha = mesh->device.malloc(acoustics->NImagPoles*sizeof(dfloat), acoustics->LRAlpha);
  acoustics->o_LRBeta = mesh->device.malloc(acoustics->NImagPoles*sizeof(dfloat), acoustics->LRBeta);

  mesh->NboundaryPointsLocal = mesh->Nfp*mesh->NboundaryFacesLocal;
  mesh->NLRPoints = mesh->Nfp*mesh->NLRFaces;
  
  dlong NLRPointsAlloc = mesh->NLRPoints ? mesh->NLRPoints:1; // [EA] Occa cannot have pointers to empty arrays

  mesh->mapAccToQ = (dlong*) calloc(NLRPointsAlloc, sizeof(dlong));

  hlong counter = 0;
  for(hlong e = 0; e < mesh->Nelements;e++){
    for(hlong n = 0; n < mesh->Nfp*mesh->Nfaces; n++){
      hlong id  = e*mesh->Nfp*mesh->Nfaces + n;
      if(mesh->mapAcc[id] >= 0){
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
      }
    }
  }
  printf("counter:%d, NLRPoints:%d\n",counter,NLRPointsAlloc);
  // [EA] mapAcc to device
  mesh->o_mapAcc =
    mesh->device.malloc(mesh->Nelements*mesh->Nfp*mesh->Nfaces*sizeof(dlong),
                        mesh->mapAcc);

  
  mesh->o_mapAccToQ =
    mesh->device.malloc(NLRPointsAlloc*sizeof(dlong),
                        mesh->mapAccToQ);


  // [EA] erk and esdirk to device
  mesh->o_erka = 
    mesh->device.malloc(36*sizeof(dfloat), mesh->erka);
  mesh->o_erkb = 
    mesh->device.malloc(6*sizeof(dfloat), mesh->erkb);
  mesh->o_esdirka = 
    mesh->device.malloc(36*sizeof(dfloat), mesh->esdirka);
  mesh->o_esdirkb = 
    mesh->device.malloc(6*sizeof(dfloat), mesh->esdirkb);





  acoustics->acc = 
    (dfloat*) calloc(NLRPointsAlloc*acoustics->Npoles, sizeof(dfloat));
  acoustics->rhsacc = 
    (dfloat*) calloc(NLRPointsAlloc*acoustics->Npoles, sizeof(dfloat));
  acoustics->resacc = 
    (dfloat*) calloc(NLRPointsAlloc*acoustics->Npoles, sizeof(dfloat));
  acoustics->o_acc =
    mesh->device.malloc(NLRPointsAlloc*acoustics->Npoles*sizeof(dfloat), acoustics->acc);
  acoustics->o_rhsacc =
    mesh->device.malloc(NLRPointsAlloc*acoustics->Npoles*sizeof(dfloat), acoustics->rhsacc);
  acoustics->o_resacc =
    mesh->device.malloc(NLRPointsAlloc*acoustics->Npoles*sizeof(dfloat), acoustics->resacc);

  cout << "TIME INTEGRATOR (" << newOptions.getArgs("TIME INTEGRATOR") << ")" << endl;
  
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

  if (newOptions.compareArgs("TIME INTEGRATOR","EIRK4")){
    acoustics->k1acc = (dfloat*) calloc(NLRPointsAlloc*acoustics->Npoles, sizeof(dfloat));
    acoustics->k2acc = (dfloat*) calloc(NLRPointsAlloc*acoustics->Npoles, sizeof(dfloat));
    acoustics->k3acc = (dfloat*) calloc(NLRPointsAlloc*acoustics->Npoles, sizeof(dfloat));
    acoustics->k4acc = (dfloat*) calloc(NLRPointsAlloc*acoustics->Npoles, sizeof(dfloat));
    acoustics->k5acc = (dfloat*) calloc(NLRPointsAlloc*acoustics->Npoles, sizeof(dfloat));
    acoustics->k6acc = (dfloat*) calloc(NLRPointsAlloc*acoustics->Npoles, sizeof(dfloat));

    acoustics->Xacc = (dfloat*) calloc(NLRPointsAlloc*acoustics->Npoles, sizeof(dfloat));

    acoustics->k1rhsq = (dfloat*) calloc(mesh->Np*mesh->Nelements*mesh->Nfields, sizeof(dfloat));
    acoustics->k2rhsq = (dfloat*) calloc(mesh->Np*mesh->Nelements*mesh->Nfields, sizeof(dfloat));
    acoustics->k3rhsq = (dfloat*) calloc(mesh->Np*mesh->Nelements*mesh->Nfields, sizeof(dfloat));
    acoustics->k4rhsq = (dfloat*) calloc(mesh->Np*mesh->Nelements*mesh->Nfields, sizeof(dfloat));
    acoustics->k5rhsq = (dfloat*) calloc(mesh->Np*mesh->Nelements*mesh->Nfields, sizeof(dfloat));
    acoustics->k6rhsq = (dfloat*) calloc(mesh->Np*mesh->Nelements*mesh->Nfields, sizeof(dfloat));

    acoustics->resq = (dfloat*) calloc(mesh->Np*mesh->Nelements*mesh->Nfields, sizeof(dfloat));

    acoustics->o_k1acc = 
          mesh->device.malloc(NLRPointsAlloc*acoustics->Npoles*sizeof(dfloat), acoustics->k1acc);
    acoustics->o_k2acc = 
          mesh->device.malloc(NLRPointsAlloc*acoustics->Npoles*sizeof(dfloat), acoustics->k2acc);
    acoustics->o_k3acc = 
          mesh->device.malloc(NLRPointsAlloc*acoustics->Npoles*sizeof(dfloat), acoustics->k3acc);
    acoustics->o_k4acc = 
          mesh->device.malloc(NLRPointsAlloc*acoustics->Npoles*sizeof(dfloat), acoustics->k4acc);
    acoustics->o_k5acc = 
          mesh->device.malloc(NLRPointsAlloc*acoustics->Npoles*sizeof(dfloat), acoustics->k5acc);
    acoustics->o_k6acc = 
          mesh->device.malloc(NLRPointsAlloc*acoustics->Npoles*sizeof(dfloat), acoustics->k6acc);

    acoustics->o_Xacc = 
          mesh->device.malloc(NLRPointsAlloc*acoustics->Npoles*sizeof(dfloat), acoustics->resacc);

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
          mesh->device.malloc(mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat), acoustics->resq);
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
  
  kernelInfo["defines/" "p_Nfields"]= mesh->Nfields;
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

  int NblockS = 1024/maxNodes; // works for CUDA
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  int cubMaxNodes = mymax(mesh->Np, (mesh->intNfp*mesh->Nfaces));
  kernelInfo["defines/" "p_cubMaxNodes"]= cubMaxNodes;
  int cubMaxNodes1 = mymax(mesh->Np, (mesh->intNfp));
  kernelInfo["defines/" "p_cubMaxNodes1"]= cubMaxNodes1;

  kernelInfo["defines/" "p_Lambda2"]= 0.5f;

  kernelInfo["defines/" "p_blockSize"]= blockSize;


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

  printf("fileName=[ %s ] \n", fileName);
  printf("kernelName=[ %s ] \n", kernelName);
  
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
  acoustics->acousticsUpdateEIRK4 = 
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsUpdate.okl",
				       "acousticsUpdateEIRK4",
				       kernelInfo);
  acoustics->acousticsUpdateEIRK4Acc = 
    mesh->device.buildKernel(DACOUSTICS "/okl/acousticsUpdate.okl",
				       "acousticsUpdateEIRK4Acc",
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
              
  // fix this later
  mesh->haloExtractKernel =
    mesh->device.buildKernel(DHOLMES "/okl/meshHaloExtract3D.okl",
				       "meshHaloExtract3D",
				       kernelInfo);
  return acoustics;
}
