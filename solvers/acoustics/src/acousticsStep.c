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

void acousticsVolumeKernel(acoustics_t *acoustics, occa::memory qPtr, occa::memory rhsqPtr){
  mesh_t *mesh = acoustics->mesh;
  if(!mesh->Ncurv){
    acoustics->volumeKernel(mesh->Nelements, 
          mesh->o_vgeo, 
          mesh->o_Dmatrices,
          qPtr, 
          rhsqPtr);
  } else {
    acoustics->volumeKernelCurv(mesh->Nelements,
        mesh->o_vgeo,
        mesh->o_vgeoCurv,
        mesh->o_mapCurv,
        mesh->o_Dmatrices,
        qPtr,
        rhsqPtr);
  }
}

void acousticsSurfaceKernel(acoustics_t *acoustics, occa::memory qPtr, occa::memory rhsqPtr,
                      occa::memory accPtr, occa::memory rhsaccPtr, const dfloat currentTime){
  mesh_t *mesh = acoustics->mesh;
  if(!mesh->Ncurv){
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
		       qPtr, 
		       rhsqPtr,
           accPtr,
           rhsaccPtr,
           mesh->o_mapAcc,
           acoustics->o_LR,
           acoustics->o_LRInfo,
           mesh->NLRPoints,
           acoustics->o_anglei,
           acoustics->o_ER,
           acoustics->o_ERInfo);
    } else {
      acoustics->surfaceKernelCurv(mesh->Nelements, 
		       mesh->o_sgeo, 
           mesh->o_sgeoCurv,
           mesh->o_mapCurv,
		       mesh->o_LIFTT, 
		       mesh->o_vmapM, 
		       mesh->o_vmapP, 
		       mesh->o_EToB,
		       currentTime, 
		       mesh->o_x, 
		       mesh->o_y,
		       mesh->o_z, 
		       qPtr, 
		       rhsqPtr,
           accPtr,
           rhsaccPtr,
           mesh->o_mapAcc,
           acoustics->o_LR,
           acoustics->o_LRInfo,
           mesh->NLRPoints,
           acoustics->o_anglei,
           acoustics->o_ER,
           acoustics->o_ERInfo);
    }
}

void acousticsDopriStep(acoustics_t *acoustics, const dfloat time){

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

void acousticsLserkStep(acoustics_t *acoustics, const dfloat time) {
  mesh_t *mesh = acoustics->mesh;
  
  // [EA] Angle detection using wave-splitting
  if(acoustics->NERPointsTotal){
    // Do interpolation and send/receive wave-splitting points from other ranks 
    if(acoustics->NERComPoints || acoustics->NComPointsToSendAllRanks){
      
      if(acoustics->NComPointsToSendAllRanks){
        acoustics->acousticsWSComInterpolation(acoustics->NComPointsToSendAllRanks,
																		acoustics->o_comPointsToSend,
																		acoustics->o_ERintpolElementsCom,
																		acoustics->o_ERintpolCom,
																		acoustics->o_q,
																		acoustics->o_vtSend);
      }
      
      if(acoustics->NERComPoints){
        acoustics->ERInsertComVT(acoustics->NERComPoints,
													acoustics->o_comPointsIdxAll,
													acoustics->o_ERComPointsIdx,
													acoustics->o_vtRecv,
													acoustics->o_vt,
													mesh->rank);
      }
    }
    acoustics->ERangleDetection(mesh->NERPoints,
														 mesh->NLRPoints,
														 acoustics->o_vt,
														 acoustics->o_vi,
														 acoustics->o_ERintpolElements,
														 acoustics->o_q,
														 mesh->o_sgeo,
														 acoustics->o_ERintpol,
														 acoustics->o_anglei,
														 mesh->o_mapAccToQ,
														 mesh->o_mapAccToN,
														 mesh->dt,
                             mesh->rank);


    // Move vt time steps
    acoustics->ERMoveVT(mesh->NERPoints, acoustics->o_vt);
  }

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

    acousticsVolumeKernel(acoustics, acoustics->o_q, acoustics->o_rhsq);

    // wait for q halo data to arrive
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);
      
      // copy halo data to DEVICE
      size_t offset = mesh->Np*acoustics->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      acoustics->o_q.copyFrom(acoustics->recvBuffer, acoustics->haloBytes, offset);
    }

    acousticsSurfaceKernel(acoustics, acoustics->o_q, acoustics->o_rhsq,
                      acoustics->o_acc, acoustics->o_rhsacc, currentTime);
    
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
            acoustics->LRInfo[0],
            mesh->dt,  
            mesh->rka[rk],
            mesh->rkb[rk],
            acoustics->o_rhsacc,
            acoustics->o_resacc,
            acoustics->o_acc);
    }
  }
}

void acousticsEirkStep(acoustics_t *acoustics, const dfloat time) {
  mesh_t *mesh = acoustics->mesh;

  // [EA] Angle detection using wave-splitting
  if(acoustics->NERPointsTotal){
    // Do interpolation and send/receive wave-splitting points from other ranks 
    if(acoustics->NERComPoints || acoustics->NComPointsToSendAllRanks){
      
      if(acoustics->NComPointsToSendAllRanks){
      acoustics->acousticsWSComInterpolation(acoustics->NComPointsToSendAllRanks,
																		acoustics->o_comPointsToSend,
																		acoustics->o_ERintpolElementsCom,
																		acoustics->o_ERintpolCom,
																		acoustics->o_q,
																		acoustics->o_vtSend);
      }
      
      if(acoustics->NERComPoints){
      acoustics->ERInsertComVT(acoustics->NERComPoints,
													acoustics->o_comPointsIdxAll,
													acoustics->o_ERComPointsIdx,
													acoustics->o_vtRecv,
													acoustics->o_vt,
													mesh->rank);
      }
    }

    if(mesh->NERPoints){
      acoustics->ERangleDetection(mesh->NERPoints,
                                mesh->NLRPoints,
                                acoustics->o_vt,
                                acoustics->o_vi,
                                acoustics->o_ERintpolElements,
                                acoustics->o_q,
                                mesh->o_sgeo,
                                acoustics->o_ERintpol,
                                acoustics->o_anglei,
                                mesh->o_mapAccToQ,
                                mesh->o_mapAccToN,
                                mesh->dt,
                                mesh->rank);


      // Move vt time steps
      acoustics->ERMoveVT(mesh->NERPoints, acoustics->o_vt);    
    }
  }
  
  for(int s = 0; s < 6; s++){
    dfloat currentTime = time + mesh->erkc[s]*mesh->dt;
    
    // [EA] Make pointers for each stage
    occa::memory qPtr;
    occa::memory rhsqPtr;
    occa::memory accPtr;
    occa::memory rhsaccPtr;
    switch(s){
      case 0:
        qPtr = acoustics->o_q;
        rhsqPtr = acoustics->o_k1rhsq;
        accPtr = acoustics->o_acc;
        rhsaccPtr = acoustics->o_k1acc;
        break;
      case 1:
        qPtr = acoustics->o_resq;
        rhsqPtr = acoustics->o_k2rhsq;
        accPtr = acoustics->o_Xacc;
        rhsaccPtr = acoustics->o_k2acc;
        break;
      case 2:
        qPtr = acoustics->o_resq;
        rhsqPtr = acoustics->o_k3rhsq;
        accPtr = acoustics->o_Xacc;
        rhsaccPtr = acoustics->o_k3acc;
        break;
      case 3:
        qPtr = acoustics->o_resq;
        rhsqPtr = acoustics->o_k4rhsq;
        accPtr = acoustics->o_Xacc;
        rhsaccPtr = acoustics->o_k4acc;
        break;
      case 4:
        qPtr = acoustics->o_resq;
        rhsqPtr = acoustics->o_k5rhsq;
        accPtr = acoustics->o_Xacc;
        rhsaccPtr = acoustics->o_k5acc;
        break;
      case 5:
        qPtr = acoustics->o_resq;
        rhsqPtr = acoustics->o_k6rhsq;
        accPtr = acoustics->o_Xacc;
        rhsaccPtr = acoustics->o_k6acc;
        break;
    }
    
    // extract q halo on DEVICE
    if(mesh->totalHaloPairs>0){
      int Nentries = mesh->Np*acoustics->Nfields;
      mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, qPtr, acoustics->o_haloBuffer);
        
      // copy extracted halo to HOST 
      acoustics->o_haloBuffer.copyTo(acoustics->sendBuffer);      
        
      // start halo exchange
      meshHaloExchangeStart(mesh, mesh->Np*acoustics->Nfields*sizeof(dfloat), acoustics->sendBuffer, acoustics->recvBuffer);
    }

    acousticsVolumeKernel(acoustics, qPtr, rhsqPtr);

    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);
        
      // copy halo data to DEVICE
      size_t offset = mesh->Np*acoustics->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      qPtr.copyFrom(acoustics->recvBuffer, acoustics->haloBytes, offset);
    }
    
    acousticsSurfaceKernel(acoustics, qPtr, rhsqPtr, accPtr, rhsaccPtr, currentTime);

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
            s+1);

    if(mesh->NLRPoints){
    acoustics->acousticsUpdateEIRK4AccLR(mesh->NLRPoints,
            mesh->dt,  
            mesh->o_esdirka,
            mesh->o_esdirkb,
            mesh->o_mapAccToQ,
            acoustics->o_LR,
            acoustics->o_LRInfo,
            acoustics->o_k1acc,
            acoustics->o_k2acc,
            acoustics->o_k3acc,
            acoustics->o_k4acc,
            acoustics->o_k5acc,
            acoustics->o_k6acc,
            acoustics->o_resq,
            acoustics->o_acc,
            acoustics->o_Xacc,
            s+1);

    }
  }
}
