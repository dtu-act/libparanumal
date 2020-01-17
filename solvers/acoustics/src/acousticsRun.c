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


void acousticsBCChange(acoustics_t *acoustics, dfloat time){

  mesh_t *mesh = acoustics->mesh;
  if(acoustics->BCChangeTime > 0.0 && acoustics->BCChangeTime < time){
  
    // Change ER BC to LR BC ( 4 to 3 )
    for(int ei = 0; ei < mesh->Nelements*mesh->Nfaces; ei++){
      if(mesh->EToB[ei] == 4){
        mesh->EToB[ei] = 3;
      }
    }
    mesh->o_EToB.copyFrom(mesh->EToB);

    printf("Swapping from ER to LR at time = %f\n",time);

    // Reset acc to 0
    acoustics->o_acc.copyFrom(acoustics->acc);

    // Swap ER points to LR
    mesh->NLRPoints = mesh->NERPoints;
    mesh->NERPoints = 0;

    // Make sure this swap only happens once
    acoustics->BCChangeTime = 0.0;
  }

}

void acousticsRun(acoustics_t *acoustics, setupAide &newOptions){

  mesh_t *mesh = acoustics->mesh;

  //acousticsReport(acoustics, 0, newOptions);

  //occa::timer timer; // OCCA TIMER CAUSES CRASHES
  
  //timer.initTimer(mesh->device);
  //timer.tic("Run");
  if (newOptions.compareArgs("TIME INTEGRATOR","DOPRI5")) {
    
    // hard code this for the moment
    dfloat outputInterval;
    newOptions.getArgs("OUTPUT INTERVAL", outputInterval);
    
    dfloat nextOutputTime = outputInterval;
    dfloat outputNumber = 0;
    
    //initial time
    dfloat time = 0.0;
    int tstep=0, allStep = 0;

    int done =0;
    while (!done) {

      acoustics->advSwitch = 1;
      
      if (mesh->dt<acoustics->dtMIN){
        printf("ERROR: Time step became too small at time step=%d\n", tstep);
        exit (-1);
      }
      if (isnan(mesh->dt)) {
        printf("ERROR: Solution became unstable at time step=%d\n", tstep);
        exit (-1);
      }

      //check for final timestep
      if (time+mesh->dt > mesh->finalTime){
	mesh->dt = mesh->finalTime-time;
	done = 1;
      }

      // try a step with the current time step
      acousticsDopriStep(acoustics, newOptions, time);

      // compute Dopri estimator
      dfloat err = acousticsDopriEstimate(acoustics);
					 
      // build controller
      dfloat fac1 = pow(err,acoustics->exp1);
      dfloat fac = fac1/pow(acoustics->facold,acoustics->beta);

      fac = mymax(acoustics->invfactor2, mymin(acoustics->invfactor1,fac/acoustics->safe));
      dfloat dtnew = mesh->dt/fac;

      if (err<1.0) { //dt is accepted

	// check for output during this step and do a mini-step
	if(time<nextOutputTime && time+mesh->dt>nextOutputTime){
	  dfloat savedt = mesh->dt;
	  
	  // save rkq
	  acoustics->o_saveq.copyFrom(acoustics->o_rkq);

	  // change dt to match output
	  mesh->dt = nextOutputTime-time;

	  // print
	  printf("Taking output mini step: %g\n", mesh->dt);
	  
	  // time step to output
	  acousticsDopriStep(acoustics, newOptions, time);	  

	  // shift for output
	  acoustics->o_rkq.copyTo(acoustics->o_q);
	  
	  // output  (print from rkq)
	  acousticsReport(acoustics, nextOutputTime, newOptions);

	  // restore time step
	  mesh->dt = savedt;

	  // increment next output time
	  nextOutputTime += outputInterval;

	  // accept saved rkq
	  acoustics->o_q.copyFrom(acoustics->o_saveq);
	}
	else{
	  // accept rkq
	  acoustics->o_q.copyFrom(acoustics->o_rkq);
	}

        time += mesh->dt;

        acoustics->facold = mymax(err,1E-4); // hard coded factor ?

	printf("\r time = %g (%d), dt = %g accepted                      ", time, allStep,  mesh->dt);
        tstep++;
      } else {
        dtnew = mesh->dt/(mymax(acoustics->invfactor1,fac1/acoustics->safe));
	printf("\r time = %g (%d), dt = %g rejected, trying %g", time, allStep, mesh->dt, dtnew);

	done = 0;
      }
      mesh->dt = dtnew;
      allStep++;

    }
    
    mesh->device.finish();
    
    double elapsed  = 1.0;//timer.toc("Run");
    printf("run took %lg seconds for %d accepted steps and %d total steps\n", elapsed, tstep, allStep);
    
  } else if (newOptions.compareArgs("TIME INTEGRATOR","LSERK4")) {
    dfloat time = 0.0;
    dlong snapshotCounter = 0;
    dlong snapshotFlag = 1;
    dlong snapshotTotal = 0;
    for(int tstep=0;tstep<mesh->NtimeSteps;++tstep){

      // [EA] Snapshot solution
      if(acoustics->snapshot){
        acoustics->Snapshott[snapshotCounter] = time;
        if(tstep % acoustics->snapshot == 0 && snapshotTotal < acoustics->snapshotMax){
          dlong offset = mesh->Np*mesh->Nelements*mesh->Nfields*snapshotCounter; //where to write to qSnapshot
          acoustics->o_q.copyTo(acoustics->qSnapshot+offset,
              mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat),
              0);
          snapshotCounter++;
          snapshotTotal++;
        }
        if(snapshotCounter == acoustics->writeSnapshotEvery || (tstep == mesh->NtimeSteps-1 && snapshotCounter > 0)){
          acousticsSnapshot(acoustics, time, newOptions, snapshotFlag, snapshotCounter);
          snapshotCounter = 0;
          snapshotFlag = 0;
        }
      }

      time = tstep*mesh->dt;
      
      // [EA] Change BC from ER to LR at time acoustics->BCChangeTime
      acousticsBCChange(acoustics, time);

      acousticsLserkStep(acoustics, newOptions, time);


      if(tstep % 500 == 0 && !mesh->rank){
        printf("LSERK4 - Step: %d, out of: %d\n",tstep, mesh->NtimeSteps);
      }
    }
  } else if (newOptions.compareArgs("TIME INTEGRATOR","EIRK4")){
    dfloat time = 0.0;
    dlong snapshotCounter = 0;
    dlong snapshotFlag = 1;
    dlong snapshotTotal = 0;
    for(int tstep=0;tstep<mesh->NtimeSteps;++tstep){

      // [EA] Snapshot solution
      if(acoustics->snapshot){
        acoustics->Snapshott[snapshotCounter] = time;
        if(tstep % acoustics->snapshot == 0 && snapshotTotal < acoustics->snapshotMax){
          dlong offset = mesh->Np*mesh->Nelements*mesh->Nfields*snapshotCounter; //where to write to qSnapshot
          acoustics->o_q.copyTo(acoustics->qSnapshot+offset,
              mesh->Np*mesh->Nelements*mesh->Nfields*sizeof(dfloat),
              0);
          snapshotCounter++;
          snapshotTotal++;
        }
        if(snapshotCounter == acoustics->writeSnapshotEvery || (tstep == mesh->NtimeSteps-1 && snapshotCounter > 0)){
          acousticsSnapshot(acoustics, time, newOptions, snapshotFlag, snapshotCounter);
          snapshotCounter = 0;
          snapshotFlag = 0;
        }
      }


      time = tstep*mesh->dt;

      // [EA] Change BC from ER to LR at time acoustics->BCChangeTime
      acousticsBCChange(acoustics, time);

      acousticsEirkStep(acoustics, newOptions, time);

      if(tstep % 500 == 0 && !mesh->rank){
        printf("EIRK4 - Step: %d, out of: %d\n",tstep, mesh->NtimeSteps);
      }
    }
  } else if(newOptions.compareArgs("TIME INTEGRATOR","EIRK4ADAP")){
    dfloat time = 0.0;
    dfloat dt = mesh->dt; // Initial dt
    while(time < mesh->finalTime){
      
      if(time+dt > mesh->finalTime){
        dt = mesh->finalTime - time;
      }

      // Solver
        // New solver write to rkq instead of q and rkacc instead of acc or do some pointer magic stuff

      // Time step controller
        // Calculate error from k1k2k3k4k5k6
      acoustics->acousticsErrorEIRK4(mesh->Nelements,
		            dt,  
		            mesh->o_erke,
		            acoustics->o_k1rhsq,
                acoustics->o_k2rhsq,
                acoustics->o_k3rhsq,
                acoustics->o_k4rhsq,
                acoustics->o_k5rhsq,
                acoustics->o_k6rhsq,
                acoustics->o_rkerr,
                acoustics->o_q);

      dlong accLength = mesh->NLRPoints*acoustics->LRNpoles + mesh->NERPoints*acoustics->ERNpoles;
      acoustics->acousticsErrorEIRK4Acc(accLength,
		            dt,  
		            mesh->o_esdirke,
		            acoustics->o_k1acc,
                acoustics->o_k2acc,
                acoustics->o_k3acc,
                acoustics->o_k4acc,
                acoustics->o_k5acc,
                acoustics->o_k6acc,
                acoustics->o_rkerrAcc,
                acoustics->o_acc);


      acoustics->acousticsErrorEIRK4r(mesh->Nelements,
                acoustics->o_rkerr,
                acoustics->o_rkq);
      
      acoustics->acousticsErrorEIRK4Accr(accLength,
                acoustics->o_rkerrAcc,
                acoustics->o_rkAcc);

    }
  }
  // [EA] Copy remaining o_qRecv from device to host
  if(acoustics->NReceiversLocal > 0){
    for(int iRecv = 0; iRecv < acoustics->NReceiversLocal; iRecv++){
      dlong offset = recvCopyRate*acoustics->qRecvCopyCounter + mesh->NtimeSteps*iRecv;
      
      acoustics->o_qRecv.copyTo(acoustics->qRecv+offset,
            acoustics->qRecvCounter*sizeof(dfloat), 
            recvCopyRate*iRecv*sizeof(dfloat));  
    }

  }
  if(acoustics->snapshot){
    acousticsSnapshotXYZ(acoustics, newOptions);
  }
}
