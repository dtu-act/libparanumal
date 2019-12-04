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
  
  // [EA] Angle detection using wave-splitting
  if(mesh->NERPoints){
    // Do interpolation and send/receive wave-splitting points from other ranks 
    if(acoustics->NERComPoints || acoustics->NComPointsToSendAllRanks){
      
      acoustics->acousticsWSComInterpolation(acoustics->NComPointsToSendAllRanks,
																		acoustics->o_comPointsToSend,
																		acoustics->o_ERintpolElementsCom,
																		acoustics->o_ERintpolCom,
																		acoustics->o_q,
																		acoustics->o_vtSend);

      
      // Communicate wave-splitting points between ranks
      acousticsWSExchange(acoustics); 
      acoustics->ERInsertComVT(acoustics->NERComPoints,
													acoustics->o_comPointsIdxAll,
													acoustics->o_ERComPointsIdx,
													acoustics->o_vtRecv,
													acoustics->o_vt,
                          mesh->rank);
    }

    //printf("r = %d, nerpoints = %d\n",mesh->rank,mesh->NERPoints);
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

    #if 0
    dfloat ttt[91];
    for(int itt = 0; itt < 91; itt++){
      ttt[itt] = -1;
      printf("ttt = %g\n",ttt[itt]);
    }    
    acoustics->o_ERYinf.copyTo(ttt);
    for(int itt = 0; itt < 91; itt++){
      printf("yinf = %g\n",ttt[itt]);
    }
    #endif


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
           acoustics->LRNpoles,
           acoustics->LRNRealPoles,
           acoustics->LRNImagPoles,
           mesh->NLRPoints,
           acoustics->o_anglei,
           acoustics->o_ERYinf,
           acoustics->ERNRealPoles,
           acoustics->ERNImagPoles,
           acoustics->ERNpoles,
           acoustics->o_ERLambda,
           acoustics->o_ERA,
           acoustics->o_ERAlpha,
           acoustics->o_ERBeta,
           acoustics->o_ERB,
           acoustics->o_ERC);

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
            acoustics->LRNpoles,
            mesh->dt,  
            mesh->rka[rk],
            mesh->rkb[rk],
            acoustics->o_rhsacc,
            acoustics->o_resacc,
            acoustics->o_acc);
    }
    if(mesh->NERPoints){
      acoustics->updateKernelER(mesh->NERPoints,
            acoustics->ERNpoles,
            mesh->dt,  
            mesh->rka[rk],
            mesh->rkb[rk],
            acoustics->o_rhsacc,
            acoustics->o_resacc,
            acoustics->o_acc,
            acoustics->LRNpoles*mesh->NLRPoints);
    }
  }

  //---------RECEIVER---------
  acoustics->acousticsReceiverInterpolation(acoustics->NReceiversLocal,
																		acoustics->o_qRecv,
																		acoustics->o_recvElements,
                                    acoustics->o_recvElementsIdx,
																		acoustics->o_recvintpol,
																		acoustics->o_q,
																		acoustics->qRecvCounter);
  acoustics->qRecvCounter++;
  if(acoustics->qRecvCounter == recvCopyRate){
    for(int iRecv = 0; iRecv < acoustics->NReceiversLocal; iRecv++){
      dlong offset = recvCopyRate*acoustics->qRecvCopyCounter + mesh->NtimeSteps*iRecv;

      acoustics->o_qRecv.copyTo(acoustics->qRecv+offset,
            acoustics->qRecvCounter*sizeof(dfloat), 
            recvCopyRate*iRecv*sizeof(dfloat));  
    }
    acoustics->qRecvCounter = 0;
    acoustics->qRecvCopyCounter++;  
  } 
  //---------RECEIVER---------
}

void acousticsEirkStep(acoustics_t *acoustics, setupAide &newOptions, const dfloat time){

  mesh_t *mesh = acoustics->mesh;
  
  // [EA] Angle detection using wave-splitting
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
														 mesh->dt);


    // Move vt time steps
    acoustics->ERMoveVT(mesh->NERPoints, acoustics->o_vt);    
  }


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
           acoustics->LRNpoles,
           acoustics->LRNRealPoles,
           acoustics->LRNImagPoles,
           mesh->NLRPoints,
           acoustics->o_anglei,
           acoustics->o_ERYinf,
           acoustics->ERNRealPoles,
           acoustics->ERNImagPoles,
           acoustics->ERNpoles,
           acoustics->o_ERLambda,
           acoustics->o_ERA,
           acoustics->o_ERAlpha,
           acoustics->o_ERBeta,
           acoustics->o_ERB,
           acoustics->o_ERC);
  
 
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
  acoustics->acousticsUpdateEIRK4AccLR(mesh->NLRPoints,
		      mesh->dt,  
		      mesh->o_esdirka,
		      mesh->o_esdirkb,
          acoustics->LRNpoles,
          acoustics->LRNRealPoles,
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
  if(mesh->NERPoints){
    acoustics->acousticsUpdateEIRK4AccER(mesh->NERPoints,
					mesh->NLRPoints,
					acoustics->LRNpoles,
		      mesh->dt,  
		      mesh->o_esdirka,
		      mesh->o_esdirkb,
          acoustics->ERNpoles,
          acoustics->ERNRealPoles,
          mesh->o_mapAccToQ,
          acoustics->o_ERAlpha,
          acoustics->o_ERBeta,
          acoustics->o_ERLambda,
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
           acoustics->LRNpoles,
           acoustics->LRNRealPoles,
           acoustics->LRNImagPoles,
           mesh->NLRPoints,
           acoustics->o_anglei,
           acoustics->o_ERYinf,
           acoustics->ERNRealPoles,
           acoustics->ERNImagPoles,
           acoustics->ERNpoles,
           acoustics->o_ERLambda,
           acoustics->o_ERA,
           acoustics->o_ERAlpha,
           acoustics->o_ERBeta,
           acoustics->o_ERB,
           acoustics->o_ERC);
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
  acoustics->acousticsUpdateEIRK4AccLR(mesh->NLRPoints,
		      mesh->dt,  
		      mesh->o_esdirka,
		      mesh->o_esdirkb,
          acoustics->LRNpoles,
          acoustics->LRNRealPoles,
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
    if(mesh->NERPoints){
    acoustics->acousticsUpdateEIRK4AccER(mesh->NERPoints,
					mesh->NLRPoints,
					acoustics->LRNpoles,
		      mesh->dt,  
		      mesh->o_esdirka,
		      mesh->o_esdirkb,
          acoustics->ERNpoles,
          acoustics->ERNRealPoles,
          mesh->o_mapAccToQ,
          acoustics->o_ERAlpha,
          acoustics->o_ERBeta,
          acoustics->o_ERLambda,
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
           acoustics->LRNpoles,
           acoustics->LRNRealPoles,
           acoustics->LRNImagPoles,
           mesh->NLRPoints,
           acoustics->o_anglei,
           acoustics->o_ERYinf,
           acoustics->ERNRealPoles,
           acoustics->ERNImagPoles,
           acoustics->ERNpoles,
           acoustics->o_ERLambda,
           acoustics->o_ERA,
           acoustics->o_ERAlpha,
           acoustics->o_ERBeta,
           acoustics->o_ERB,
           acoustics->o_ERC);

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
  acoustics->acousticsUpdateEIRK4AccLR(mesh->NLRPoints,
		      mesh->dt,  
		      mesh->o_esdirka,
		      mesh->o_esdirkb,
          acoustics->LRNpoles,
          acoustics->LRNRealPoles,
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
  if(mesh->NERPoints){
  acoustics->acousticsUpdateEIRK4AccER(mesh->NERPoints,
        mesh->NLRPoints,
        acoustics->LRNpoles,
        mesh->dt,  
        mesh->o_esdirka,
        mesh->o_esdirkb,
        acoustics->ERNpoles,
        acoustics->ERNRealPoles,
        mesh->o_mapAccToQ,
        acoustics->o_ERAlpha,
        acoustics->o_ERBeta,
        acoustics->o_ERLambda,
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
           acoustics->LRNpoles,
           acoustics->LRNRealPoles,
           acoustics->LRNImagPoles,
           mesh->NLRPoints,
           acoustics->o_anglei,
           acoustics->o_ERYinf,
           acoustics->ERNRealPoles,
           acoustics->ERNImagPoles,
           acoustics->ERNpoles,
           acoustics->o_ERLambda,
           acoustics->o_ERA,
           acoustics->o_ERAlpha,
           acoustics->o_ERBeta,
           acoustics->o_ERB,
           acoustics->o_ERC);

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
  acoustics->acousticsUpdateEIRK4AccLR(mesh->NLRPoints,
		      mesh->dt,  
		      mesh->o_esdirka,
		      mesh->o_esdirkb,
          acoustics->LRNpoles,
          acoustics->LRNRealPoles,
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
  if(mesh->NERPoints){
  acoustics->acousticsUpdateEIRK4AccER(mesh->NERPoints,
        mesh->NLRPoints,
        acoustics->LRNpoles,
        mesh->dt,  
        mesh->o_esdirka,
        mesh->o_esdirkb,
        acoustics->ERNpoles,
        acoustics->ERNRealPoles,
        mesh->o_mapAccToQ,
        acoustics->o_ERAlpha,
        acoustics->o_ERBeta,
        acoustics->o_ERLambda,
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
           acoustics->LRNpoles,
           acoustics->LRNRealPoles,
           acoustics->LRNImagPoles,
           mesh->NLRPoints,
           acoustics->o_anglei,
           acoustics->o_ERYinf,
           acoustics->ERNRealPoles,
           acoustics->ERNImagPoles,
           acoustics->ERNpoles,
           acoustics->o_ERLambda,
           acoustics->o_ERA,
           acoustics->o_ERAlpha,
           acoustics->o_ERBeta,
           acoustics->o_ERB,
           acoustics->o_ERC);

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
  acoustics->acousticsUpdateEIRK4AccLR(mesh->NLRPoints,
		      mesh->dt,  
		      mesh->o_esdirka,
		      mesh->o_esdirkb,
          acoustics->LRNpoles,
          acoustics->LRNRealPoles,
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
  if(mesh->NERPoints){
  acoustics->acousticsUpdateEIRK4AccER(mesh->NERPoints,
        mesh->NLRPoints,
        acoustics->LRNpoles,
        mesh->dt,  
        mesh->o_esdirka,
        mesh->o_esdirkb,
        acoustics->ERNpoles,
        acoustics->ERNRealPoles,
        mesh->o_mapAccToQ,
        acoustics->o_ERAlpha,
        acoustics->o_ERBeta,
        acoustics->o_ERLambda,
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
           acoustics->LRNpoles,
           acoustics->LRNRealPoles,
           acoustics->LRNImagPoles,
           mesh->NLRPoints,
           acoustics->o_anglei,
           acoustics->o_ERYinf,
           acoustics->ERNRealPoles,
           acoustics->ERNImagPoles,
           acoustics->ERNpoles,
           acoustics->o_ERLambda,
           acoustics->o_ERA,
           acoustics->o_ERAlpha,
           acoustics->o_ERBeta,
           acoustics->o_ERB,
           acoustics->o_ERC);

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
    acoustics->acousticsUpdateEIRK4AccLR(mesh->NLRPoints,
		      mesh->dt,  
		      mesh->o_esdirka,
		      mesh->o_esdirkb,
          acoustics->LRNpoles,
          acoustics->LRNRealPoles,
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
    if(mesh->NERPoints){
      acoustics->acousticsUpdateEIRK4AccER(mesh->NERPoints,
					mesh->NLRPoints,
					acoustics->LRNpoles,
		      mesh->dt,  
		      mesh->o_esdirka,
		      mesh->o_esdirkb,
          acoustics->ERNpoles,
          acoustics->ERNRealPoles,
          mesh->o_mapAccToQ,
          acoustics->o_ERAlpha,
          acoustics->o_ERBeta,
          acoustics->o_ERLambda,
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
acoustics->acousticsReceiverInterpolation(acoustics->NReceiversLocal,
																		acoustics->o_qRecv,
																		acoustics->o_recvElements,
                                    acoustics->o_recvElementsIdx,
																		acoustics->o_recvintpol,
																		acoustics->o_q,
																		acoustics->qRecvCounter);
  acoustics->qRecvCounter++;
  if(acoustics->qRecvCounter == recvCopyRate){
    for(int iRecv = 0; iRecv < acoustics->NReceiversLocal; iRecv++){
      dlong offset = recvCopyRate*acoustics->qRecvCopyCounter + mesh->NtimeSteps*iRecv;

      acoustics->o_qRecv.copyTo(acoustics->qRecv+offset,
            acoustics->qRecvCounter*sizeof(dfloat), 
            recvCopyRate*iRecv*sizeof(dfloat));  
    }
    acoustics->qRecvCounter = 0;
    acoustics->qRecvCopyCounter++;  
  } 
  //---------RECEIVER---------
  
}




