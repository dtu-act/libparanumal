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

void acousticsDopriStep(acoustics_t *acoustics, setupAide &newOptions, const dfloat time){

  mesh_t *mesh = acoustics->mesh;
  
  //RK step
  for(int rk=0;rk<acoustics->Nrk;++rk){
    
    // t_rk = t + C_rk*dt
    dfloat currentTime = time + acoustics->rkC[rk]*mesh->dt;
    
    //compute RK stage 
    // rkq = q + dt sum_{i=0}^{rk-1} a_{rk,i}*rhsq_i
    acoustics->rkStageKernel(mesh->Nelements,
		       rk,
		       mesh->dt,
		       acoustics->o_rkA,
		       acoustics->o_q,
		       acoustics->o_rkrhsq,
		       acoustics->o_rkq);
    
    //compute RHS
    // rhsq = F(currentTIme, rkq)

    // extract q halo on DEVICE
    if(mesh->totalHaloPairs>0){
      int Nentries = mesh->Np*acoustics->Nfields;          
      mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, acoustics->o_rkq, acoustics->o_haloBuffer);
      
      // copy extracted halo to HOST 
      acoustics->o_haloBuffer.copyTo(acoustics->sendBuffer);      
      
      // start halo exchange
      meshHaloExchangeStart(mesh, mesh->Np*acoustics->Nfields*sizeof(dfloat), acoustics->sendBuffer, acoustics->recvBuffer);
    }

    acoustics->volumeKernel(mesh->Nelements, 
		      mesh->o_vgeo, 
		      mesh->o_Dmatrices,
		      acoustics->o_rkq, 
		      acoustics->o_rhsq);

    // wait for q halo data to arrive
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);
          
      // copy halo data to DEVICE
      size_t offset = mesh->Np*acoustics->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      acoustics->o_rkq.copyFrom(acoustics->recvBuffer, acoustics->haloBytes, offset);
    }

    acoustics->surfaceKernel(mesh->Nelements, 
			     mesh->o_sgeo, 
			     mesh->o_LIFTT, 
			     mesh->o_vmapM, 
			     mesh->o_vmapP, 
			     mesh->o_EToB,
			     currentTime, 
			     mesh->o_x, 
			     mesh->o_y,
			     mesh->o_z, 
			     acoustics->o_rkq, 
			     acoustics->o_rhsq);
    
    // update solution using Runge-Kutta
    // rkrhsq_rk = rhsq
    // if rk==6 
    //   q = q + dt*sum_{i=0}^{rk} rkA_{rk,i}*rkrhs_i
    //   rkerr = dt*sum_{i=0}^{rk} rkE_{rk,i}*rkrhs_i
    acoustics->rkUpdateKernel(mesh->Nelements, 
			rk,
			mesh->dt, 
			acoustics->o_rkA, 
			acoustics->o_rkE, 
			acoustics->o_q,
			acoustics->o_rhsq, 
			acoustics->o_rkrhsq, 
			acoustics->o_rkq,
			acoustics->o_rkerr);
  }
}


void acousticsLserkStep(acoustics_t *acoustics, setupAide &newOptions, const dfloat time){

  mesh_t *mesh = acoustics->mesh;
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(int rk=0;rk<mesh->Nrk;++rk){
    dfloat currentTime = time + mesh->rkc[rk]*mesh->dt;
      
    // extract q halo on DEVICE
    if(mesh->totalHaloPairs>0){
      int Nentries = mesh->Np*acoustics->Nfields;
        
      mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, acoustics->o_q, acoustics->o_haloBuffer);
        
      // copy extracted halo to HOST 
      acoustics->o_haloBuffer.copyTo(acoustics->sendBuffer);      
        
      // start halo exchange
      meshHaloExchangeStart(mesh, mesh->Np*acoustics->Nfields*sizeof(dfloat), acoustics->sendBuffer, acoustics->recvBuffer);
    }
    acoustics->volumeKernel(mesh->Nelements, 
		      mesh->o_vgeo, 
		      mesh->o_Dmatrices,
		      acoustics->o_q, 
		      acoustics->o_rhsq);
    // wait for q halo data to arrive
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);
      
      // copy halo data to DEVICE
      size_t offset = mesh->Np*acoustics->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      acoustics->o_q.copyFrom(acoustics->recvBuffer, acoustics->haloBytes, offset);
    }
    acoustics->surfaceKernel(mesh->Nelements, 
		       mesh->o_sgeo, 
		       mesh->o_LIFTT, 
		       mesh->o_vmapM, 
		       mesh->o_vmapP, 
		       mesh->o_EToB,
		       currentTime, 
		       mesh->o_x, 
		       mesh->o_y,
		       mesh->o_z, 
		       acoustics->o_q, 
		       acoustics->o_rhsq,
           acoustics->o_acc,
           acoustics->o_rhsacc,
           mesh->o_mapAcc,
           acoustics->o_LRA,
           acoustics->o_LRB,
           acoustics->o_LRC,
           acoustics->o_LRLambda,
           acoustics->o_LRAlpha,
           acoustics->o_LRBeta,
           acoustics->LRYinf,
           acoustics->Npoles,
           acoustics->NRealPoles,
           acoustics->NImagPoles);
    // update solution using Runge-Kutta
    acoustics->updateKernel(mesh->Nelements, 
		      mesh->dt, 
		      mesh->rka[rk], 
		      mesh->rkb[rk], 
		      acoustics->o_rhsq, 
		      acoustics->o_resq, 
		      acoustics->o_q);
    if(mesh->NLRPoints){
      acoustics->updateKernelLR(mesh->NLRPoints,
            acoustics->Npoles,
            mesh->dt,  
            mesh->rka[rk],
            mesh->rkb[rk],
            acoustics->o_rhsacc,
            acoustics->o_resacc,
            acoustics->o_acc);
    }
  }

  //---------RECEIVER---------


  for(int iRecv = 0; iRecv < acoustics->NReceiversLocal; iRecv++){
    dlong qRecvOffset = acoustics->qRecvCounter*mesh->Np + iRecv*mesh->Np*mesh->NtimeSteps; // Offset writing location by # time steps taken*Np
    dlong qOffset = acoustics->recvElements[acoustics->recvElementsIdx[iRecv]]*mesh->Np*acoustics->Nfields; // Offset in q to get the correct element pres
    acoustics->o_q.copyTo(acoustics->o_qRecv,
          mesh->Np*sizeof(dfloat), 
          qRecvOffset*sizeof(dfloat),
          qOffset*sizeof(dfloat));
  }
  acoustics->qRecvCounter++;
  #if 0 
  // Old code that copied from device to host every time step.
  // copy data back to host 
  acoustics->o_q.copyTo(acoustics->q);
  // Offset in qRecv by Np*counter
  size_t qRecvOffset = acoustics->qRecvCounter*mesh->Np; // Offset writing location by # time steps taken*Np
  size_t qOffset = acoustics->recvElement*mesh->Np*acoustics->Nfields; // Offset in q to get the correct element pres

  for(int n=0;n<mesh->Np;++n){
    acoustics->qRecv[qRecvOffset+n] = acoustics->q[qOffset+n]; // See in acousticsReport where we print how to find a specific element. Currently just saves the 0th element.
  }
  acoustics->qRecvCounter++;
  #endif
  //---------RECEIVER---------
}

void acousticsEirkStep(acoustics_t *acoustics, setupAide &newOptions, const dfloat time){

  mesh_t *mesh = acoustics->mesh;
  
  
  //----------------------------------------- STAGE 1 -----------------------------------------
  dfloat currentTime = time + mesh->erkc[0]*mesh->dt;
  // extract q halo on DEVICE
  if(mesh->totalHaloPairs>0){
    int Nentries = mesh->Np*acoustics->Nfields;
      
    mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, acoustics->o_q, acoustics->o_haloBuffer);
      
    // copy extracted halo to HOST 
    acoustics->o_haloBuffer.copyTo(acoustics->sendBuffer);      
      
    // start halo exchange
    meshHaloExchangeStart(mesh, mesh->Np*acoustics->Nfields*sizeof(dfloat), acoustics->sendBuffer, acoustics->recvBuffer);
  }
  acoustics->volumeKernel(mesh->Nelements, 
	     mesh->o_vgeo, 
	     mesh->o_Dmatrices,
	     acoustics->o_q, 
	     acoustics->o_k1rhsq);

  if(mesh->totalHaloPairs>0){
    meshHaloExchangeFinish(mesh);
      
    // copy halo data to DEVICE
    size_t offset = mesh->Np*acoustics->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
    acoustics->o_q.copyFrom(acoustics->recvBuffer, acoustics->haloBytes, offset);
  }
  
  acoustics->surfaceKernel(mesh->Nelements, 
		       mesh->o_sgeo, 
		       mesh->o_LIFTT, 
		       mesh->o_vmapM, 
		       mesh->o_vmapP, 
		       mesh->o_EToB,
		       currentTime, 
		       mesh->o_x, 
		       mesh->o_y,
		       mesh->o_z, 
		       acoustics->o_q, 
		       acoustics->o_k1rhsq,
           acoustics->o_acc,
           acoustics->o_k1acc,
           mesh->o_mapAcc,
           acoustics->o_LRA,
           acoustics->o_LRB,
           acoustics->o_LRC,
           acoustics->o_LRLambda,
           acoustics->o_LRAlpha,
           acoustics->o_LRBeta,
           acoustics->LRYinf,
           acoustics->Npoles,
           acoustics->NRealPoles,
           acoustics->NImagPoles);
  
  acoustics->acousticsUpdateEIRK4(mesh->Nelements,
		      mesh->dt,  
		      mesh->o_erka,
		      mesh->o_erkb,
		      acoustics->o_k1rhsq,
          acoustics->o_k2rhsq,
          acoustics->o_k3rhsq,
          acoustics->o_k4rhsq,
          acoustics->o_k5rhsq,
          acoustics->o_k6rhsq,
		      acoustics->o_resq,
		      acoustics->o_q,
          1);
  if(mesh->NLRPoints){
  acoustics->acousticsUpdateEIRK4Acc(mesh->NLRPoints,
		      mesh->dt,  
		      mesh->o_esdirka,
		      mesh->o_esdirkb,
          acoustics->Npoles,
          acoustics->NRealPoles,
          mesh->o_mapAccToQ,
          acoustics->o_LRAlpha,
          acoustics->o_LRBeta,
          acoustics->o_LRLambda,
		      acoustics->o_k1acc,
          acoustics->o_k2acc,
          acoustics->o_k3acc,
          acoustics->o_k4acc,
          acoustics->o_k5acc,
          acoustics->o_k6acc,
		      acoustics->o_resq,
          acoustics->o_acc,
          acoustics->o_Xacc,
          1);
  }
          
  //----------------------------------------- STAGE 2 -----------------------------------------
  
  currentTime = time + mesh->erkc[1]*mesh->dt;

  // extract q halo on DEVICE
  if(mesh->totalHaloPairs>0){
    int Nentries = mesh->Np*acoustics->Nfields;
      
    mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, acoustics->o_resq, acoustics->o_haloBuffer);
      
    // copy extracted halo to HOST 
    acoustics->o_haloBuffer.copyTo(acoustics->sendBuffer);      
      
    // start halo exchange
    meshHaloExchangeStart(mesh, mesh->Np*acoustics->Nfields*sizeof(dfloat), acoustics->sendBuffer, acoustics->recvBuffer);
  }
  acoustics->volumeKernel(mesh->Nelements, 
	     mesh->o_vgeo, 
	     mesh->o_Dmatrices,
	     acoustics->o_resq, 
	     acoustics->o_k2rhsq);

  if(mesh->totalHaloPairs>0){
    meshHaloExchangeFinish(mesh);
      
    // copy halo data to DEVICE
    size_t offset = mesh->Np*acoustics->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
    acoustics->o_resq.copyFrom(acoustics->recvBuffer, acoustics->haloBytes, offset);
  }
  acoustics->surfaceKernel(mesh->Nelements, 
		       mesh->o_sgeo, 
		       mesh->o_LIFTT, 
		       mesh->o_vmapM, 
		       mesh->o_vmapP, 
		       mesh->o_EToB,
		       currentTime, 
		       mesh->o_x, 
		       mesh->o_y,
		       mesh->o_z, 
		       acoustics->o_resq, 
		       acoustics->o_k2rhsq,
           acoustics->o_Xacc,
           acoustics->o_k2acc,
           mesh->o_mapAcc,
           acoustics->o_LRA,
           acoustics->o_LRB,
           acoustics->o_LRC,
           acoustics->o_LRLambda,
           acoustics->o_LRAlpha,
           acoustics->o_LRBeta,
           acoustics->LRYinf,
           acoustics->Npoles,
           acoustics->NRealPoles,
           acoustics->NImagPoles);

  acoustics->acousticsUpdateEIRK4(mesh->Nelements,
		      mesh->dt,  
		      mesh->o_erka,
		      mesh->o_erkb,
		      acoustics->o_k1rhsq,
          acoustics->o_k2rhsq,
          acoustics->o_k3rhsq,
          acoustics->o_k4rhsq,
          acoustics->o_k5rhsq,
          acoustics->o_k6rhsq,
		      acoustics->o_resq,
		      acoustics->o_q,
          2);
  if(mesh->NLRPoints){
  acoustics->acousticsUpdateEIRK4Acc(mesh->NLRPoints,
		      mesh->dt,  
		      mesh->o_esdirka,
		      mesh->o_esdirkb,
          acoustics->Npoles,
          acoustics->NRealPoles,
          mesh->o_mapAccToQ,
          acoustics->o_LRAlpha,
          acoustics->o_LRBeta,
          acoustics->o_LRLambda,
		      acoustics->o_k1acc,
          acoustics->o_k2acc,
          acoustics->o_k3acc,
          acoustics->o_k4acc,
          acoustics->o_k5acc,
          acoustics->o_k6acc,
		      acoustics->o_resq,
          acoustics->o_acc,
          acoustics->o_Xacc,
          2);
  }
  //----------------------------------------- STAGE 3 -----------------------------------------
  
  currentTime = time + mesh->erkc[2]*mesh->dt;

  // extract q halo on DEVICE
  if(mesh->totalHaloPairs>0){
    int Nentries = mesh->Np*acoustics->Nfields;
      
    mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, acoustics->o_resq, acoustics->o_haloBuffer);
      
    // copy extracted halo to HOST 
    acoustics->o_haloBuffer.copyTo(acoustics->sendBuffer);      
      
    // start halo exchange
    meshHaloExchangeStart(mesh, mesh->Np*acoustics->Nfields*sizeof(dfloat), acoustics->sendBuffer, acoustics->recvBuffer);
  }
  acoustics->volumeKernel(mesh->Nelements, 
	     mesh->o_vgeo, 
	     mesh->o_Dmatrices,
	     acoustics->o_resq, 
	     acoustics->o_k3rhsq);

  if(mesh->totalHaloPairs>0){
    meshHaloExchangeFinish(mesh);
      
    // copy halo data to DEVICE
    size_t offset = mesh->Np*acoustics->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
    acoustics->o_resq.copyFrom(acoustics->recvBuffer, acoustics->haloBytes, offset);
  }
  acoustics->surfaceKernel(mesh->Nelements, 
		       mesh->o_sgeo, 
		       mesh->o_LIFTT, 
		       mesh->o_vmapM, 
		       mesh->o_vmapP, 
		       mesh->o_EToB,
		       currentTime, 
		       mesh->o_x, 
		       mesh->o_y,
		       mesh->o_z, 
		       acoustics->o_resq, 
		       acoustics->o_k3rhsq,
           acoustics->o_Xacc,
           acoustics->o_k3acc,
           mesh->o_mapAcc,
           acoustics->o_LRA,
           acoustics->o_LRB,
           acoustics->o_LRC,
           acoustics->o_LRLambda,
           acoustics->o_LRAlpha,
           acoustics->o_LRBeta,
           acoustics->LRYinf,
           acoustics->Npoles,
           acoustics->NRealPoles,
           acoustics->NImagPoles);

  acoustics->acousticsUpdateEIRK4(mesh->Nelements,
		      mesh->dt,  
		      mesh->o_erka,
		      mesh->o_erkb,
		      acoustics->o_k1rhsq,
          acoustics->o_k2rhsq,
          acoustics->o_k3rhsq,
          acoustics->o_k4rhsq,
          acoustics->o_k5rhsq,
          acoustics->o_k6rhsq,
		      acoustics->o_resq,
		      acoustics->o_q,
          3);
  if(mesh->NLRPoints){
  acoustics->acousticsUpdateEIRK4Acc(mesh->NLRPoints,
		      mesh->dt,  
		      mesh->o_esdirka,
		      mesh->o_esdirkb,
          acoustics->Npoles,
          acoustics->NRealPoles,
          mesh->o_mapAccToQ,
          acoustics->o_LRAlpha,
          acoustics->o_LRBeta,
          acoustics->o_LRLambda,
		      acoustics->o_k1acc,
          acoustics->o_k2acc,
          acoustics->o_k3acc,
          acoustics->o_k4acc,
          acoustics->o_k5acc,
          acoustics->o_k6acc,
		      acoustics->o_resq,
          acoustics->o_acc,
          acoustics->o_Xacc,
          3);
  }
  //----------------------------------------- STAGE 4 -----------------------------------------
  
  currentTime = time + mesh->erkc[3]*mesh->dt;

  // extract q halo on DEVICE
  if(mesh->totalHaloPairs>0){
    int Nentries = mesh->Np*acoustics->Nfields;
      
    mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, acoustics->o_resq, acoustics->o_haloBuffer);
      
    // copy extracted halo to HOST 
    acoustics->o_haloBuffer.copyTo(acoustics->sendBuffer);      
      
    // start halo exchange
    meshHaloExchangeStart(mesh, mesh->Np*acoustics->Nfields*sizeof(dfloat), acoustics->sendBuffer, acoustics->recvBuffer);
  }
  acoustics->volumeKernel(mesh->Nelements, 
	     mesh->o_vgeo, 
	     mesh->o_Dmatrices,
	     acoustics->o_resq, 
	     acoustics->o_k4rhsq);

  if(mesh->totalHaloPairs>0){
    meshHaloExchangeFinish(mesh);
      
    // copy halo data to DEVICE
    size_t offset = mesh->Np*acoustics->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
    acoustics->o_resq.copyFrom(acoustics->recvBuffer, acoustics->haloBytes, offset);
  }
  acoustics->surfaceKernel(mesh->Nelements, 
		       mesh->o_sgeo, 
		       mesh->o_LIFTT, 
		       mesh->o_vmapM, 
		       mesh->o_vmapP, 
		       mesh->o_EToB,
		       currentTime, 
		       mesh->o_x, 
		       mesh->o_y,
		       mesh->o_z, 
		       acoustics->o_resq, 
		       acoustics->o_k4rhsq,
           acoustics->o_Xacc,
           acoustics->o_k4acc,
           mesh->o_mapAcc,
           acoustics->o_LRA,
           acoustics->o_LRB,
           acoustics->o_LRC,
           acoustics->o_LRLambda,
           acoustics->o_LRAlpha,
           acoustics->o_LRBeta,
           acoustics->LRYinf,
           acoustics->Npoles,
           acoustics->NRealPoles,
           acoustics->NImagPoles);

  acoustics->acousticsUpdateEIRK4(mesh->Nelements,
		      mesh->dt,  
		      mesh->o_erka,
		      mesh->o_erkb,
		      acoustics->o_k1rhsq,
          acoustics->o_k2rhsq,
          acoustics->o_k3rhsq,
          acoustics->o_k4rhsq,
          acoustics->o_k5rhsq,
          acoustics->o_k6rhsq,
		      acoustics->o_resq,
		      acoustics->o_q,
          4);
  if(mesh->NLRPoints){
  acoustics->acousticsUpdateEIRK4Acc(mesh->NLRPoints,
		      mesh->dt,  
		      mesh->o_esdirka,
		      mesh->o_esdirkb,
          acoustics->Npoles,
          acoustics->NRealPoles,
          mesh->o_mapAccToQ,
          acoustics->o_LRAlpha,
          acoustics->o_LRBeta,
          acoustics->o_LRLambda,
		      acoustics->o_k1acc,
          acoustics->o_k2acc,
          acoustics->o_k3acc,
          acoustics->o_k4acc,
          acoustics->o_k5acc,
          acoustics->o_k6acc,
		      acoustics->o_resq,
          acoustics->o_acc,
          acoustics->o_Xacc,
          4);
  }
  //----------------------------------------- STAGE 5 -----------------------------------------
  
  currentTime = time + mesh->erkc[4]*mesh->dt;

  // extract q halo on DEVICE
  if(mesh->totalHaloPairs>0){
    int Nentries = mesh->Np*acoustics->Nfields;
      
    mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, acoustics->o_resq, acoustics->o_haloBuffer);
      
    // copy extracted halo to HOST 
    acoustics->o_haloBuffer.copyTo(acoustics->sendBuffer);      
      
    // start halo exchange
    meshHaloExchangeStart(mesh, mesh->Np*acoustics->Nfields*sizeof(dfloat), acoustics->sendBuffer, acoustics->recvBuffer);
  }
  acoustics->volumeKernel(mesh->Nelements, 
	     mesh->o_vgeo, 
	     mesh->o_Dmatrices,
	     acoustics->o_resq, 
	     acoustics->o_k5rhsq);

  if(mesh->totalHaloPairs>0){
    meshHaloExchangeFinish(mesh);
      
    // copy halo data to DEVICE
    size_t offset = mesh->Np*acoustics->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
    acoustics->o_resq.copyFrom(acoustics->recvBuffer, acoustics->haloBytes, offset);
  }
  acoustics->surfaceKernel(mesh->Nelements, 
		       mesh->o_sgeo, 
		       mesh->o_LIFTT, 
		       mesh->o_vmapM, 
		       mesh->o_vmapP, 
		       mesh->o_EToB,
		       currentTime, 
		       mesh->o_x, 
		       mesh->o_y,
		       mesh->o_z, 
		       acoustics->o_resq, 
		       acoustics->o_k5rhsq,
           acoustics->o_Xacc,
           acoustics->o_k5acc,
           mesh->o_mapAcc,
           acoustics->o_LRA,
           acoustics->o_LRB,
           acoustics->o_LRC,
           acoustics->o_LRLambda,
           acoustics->o_LRAlpha,
           acoustics->o_LRBeta,
           acoustics->LRYinf,
           acoustics->Npoles,
           acoustics->NRealPoles,
           acoustics->NImagPoles);

  acoustics->acousticsUpdateEIRK4(mesh->Nelements,
		      mesh->dt,  
		      mesh->o_erka,
		      mesh->o_erkb,
		      acoustics->o_k1rhsq,
          acoustics->o_k2rhsq,
          acoustics->o_k3rhsq,
          acoustics->o_k4rhsq,
          acoustics->o_k5rhsq,
          acoustics->o_k6rhsq,
		      acoustics->o_resq,
		      acoustics->o_q,
          5);
  if(mesh->NLRPoints){
  acoustics->acousticsUpdateEIRK4Acc(mesh->NLRPoints,
		      mesh->dt,  
		      mesh->o_esdirka,
		      mesh->o_esdirkb,
          acoustics->Npoles,
          acoustics->NRealPoles,
          mesh->o_mapAccToQ,
          acoustics->o_LRAlpha,
          acoustics->o_LRBeta,
          acoustics->o_LRLambda,
		      acoustics->o_k1acc,
          acoustics->o_k2acc,
          acoustics->o_k3acc,
          acoustics->o_k4acc,
          acoustics->o_k5acc,
          acoustics->o_k6acc,
		      acoustics->o_resq,
          acoustics->o_acc,
          acoustics->o_Xacc,
          5);
  }
  //----------------------------------------- STAGE 6 -----------------------------------------
  
  currentTime = time + mesh->erkc[5]*mesh->dt;

  // extract q halo on DEVICE
  if(mesh->totalHaloPairs>0){
    int Nentries = mesh->Np*acoustics->Nfields;
      
    mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, acoustics->o_resq, acoustics->o_haloBuffer);
      
    // copy extracted halo to HOST 
    acoustics->o_haloBuffer.copyTo(acoustics->sendBuffer);      
      
    // start halo exchange
    meshHaloExchangeStart(mesh, mesh->Np*acoustics->Nfields*sizeof(dfloat), acoustics->sendBuffer, acoustics->recvBuffer);
  }
  acoustics->volumeKernel(mesh->Nelements, 
	     mesh->o_vgeo, 
	     mesh->o_Dmatrices,
	     acoustics->o_resq, 
	     acoustics->o_k6rhsq);

  if(mesh->totalHaloPairs>0){
    meshHaloExchangeFinish(mesh);
      
    // copy halo data to DEVICE
    size_t offset = mesh->Np*acoustics->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
    acoustics->o_resq.copyFrom(acoustics->recvBuffer, acoustics->haloBytes, offset);
  }
  acoustics->surfaceKernel(mesh->Nelements, 
		       mesh->o_sgeo, 
		       mesh->o_LIFTT, 
		       mesh->o_vmapM, 
		       mesh->o_vmapP, 
		       mesh->o_EToB,
		       currentTime, 
		       mesh->o_x, 
		       mesh->o_y,
		       mesh->o_z, 
		       acoustics->o_resq, 
		       acoustics->o_k6rhsq,
           acoustics->o_Xacc,
           acoustics->o_k6acc,
           mesh->o_mapAcc,
           acoustics->o_LRA,
           acoustics->o_LRB,
           acoustics->o_LRC,
           acoustics->o_LRLambda,
           acoustics->o_LRAlpha,
           acoustics->o_LRBeta,
           acoustics->LRYinf,
           acoustics->Npoles,
           acoustics->NRealPoles,
           acoustics->NImagPoles);

  acoustics->acousticsUpdateEIRK4(mesh->Nelements,
		      mesh->dt,  
		      mesh->o_erka,
		      mesh->o_erkb,
		      acoustics->o_k1rhsq,
          acoustics->o_k2rhsq,
          acoustics->o_k3rhsq,
          acoustics->o_k4rhsq,
          acoustics->o_k5rhsq,
          acoustics->o_k6rhsq,
		      acoustics->o_resq,
		      acoustics->o_q,
          6);
  if(mesh->NLRPoints){
  acoustics->acousticsUpdateEIRK4Acc(mesh->NLRPoints,
		      mesh->dt,  
		      mesh->o_esdirka,
		      mesh->o_esdirkb,
          acoustics->Npoles,
          acoustics->NRealPoles,
          mesh->o_mapAccToQ,
          acoustics->o_LRAlpha,
          acoustics->o_LRBeta,
          acoustics->o_LRLambda,
		      acoustics->o_k1acc,
          acoustics->o_k2acc,
          acoustics->o_k3acc,
          acoustics->o_k4acc,
          acoustics->o_k5acc,
          acoustics->o_k6acc,
		      acoustics->o_resq,
          acoustics->o_acc,
          acoustics->o_Xacc,
          6);
    }

  //---------RECEIVER---------


  for(int iRecv = 0; iRecv < acoustics->NReceiversLocal; iRecv++){
    dlong qRecvOffset = acoustics->qRecvCounter*mesh->Np + iRecv*mesh->Np*mesh->NtimeSteps; // Offset writing location by # time steps taken*Np
    dlong qOffset = acoustics->recvElements[acoustics->recvElementsIdx[iRecv]]*mesh->Np*acoustics->Nfields; // Offset in q to get the correct element pres
    acoustics->o_q.copyTo(acoustics->o_qRecv,
          mesh->Np*sizeof(dfloat), 
          qRecvOffset*sizeof(dfloat),
          qOffset*sizeof(dfloat));
  }
  acoustics->qRecvCounter++;

  //---------RECEIVER---------
  
}




