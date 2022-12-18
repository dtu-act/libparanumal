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

#include "acousticsWriters.h"
#include "acoustics.h"
#include <limits>
#include <algorithm>
#include <stdexcept>

void dopri5(acoustics_t *acoustics, dfloat outputInterval);
void lseRK4(acoustics_t *acoustics, std::vector<dfloat> timeStepsWrite, std::shared_ptr<IAcousticWriter> writer);
void eiRK4(acoustics_t *acoustics, std::vector<dfloat> timeStepsWrite, std::shared_ptr<IAcousticWriter> writer);
void eiRK4adap(acoustics_t *acoustics);

// helper
void copyReceiverToHost(acoustics_t *acoustics) {
  for(int iRecv = 0; iRecv < acoustics->NReceiversLocal; iRecv++){
        dlong offset = recvCopyRate*acoustics->qRecvCopyCounter + acoustics->qRecvCounter*iRecv;

        acoustics->o_qRecv.copyTo(acoustics->qRecv + offset,
              acoustics->qRecvCounter * sizeof(dfloat), 
              recvCopyRate * iRecv * sizeof(dfloat));  
      }
}

void acousticsBCChange(acoustics_t *acoustics, dfloat time)
{
  mesh_t *mesh = acoustics->mesh;
  if (acoustics->BCChangeTime > 0.0 && acoustics->BCChangeTime < time)
  {

    // Change ER BC to LR BC ( 4 to 3 )
    for (int ei = 0; ei < mesh->Nelements * mesh->Nfaces; ei++)
    {
      if (mesh->EToB[ei] == 4)
      {
        mesh->EToB[ei] = 3;
      }
    }
    mesh->o_EToB.copyFrom(mesh->EToB);

    printf("Swapping from ER to LR at time = %f\n", time);

    // Reset acc to 0
    acoustics->o_acc.copyFrom(acoustics->acc);

    // Swap ER points to LR
    mesh->NLRPoints = mesh->NERPoints;
    mesh->NERPoints = 0;

    // Make sure this swap only happens once
    acoustics->BCChangeTime = 0.0;
  }
}

void acousticsRun(acoustics_t *acoustics, setupAide &newOptions)
{
  WriteWaveFieldType waveFieldWriteType;
  int ppwWaveField;

  if (newOptions.compareArgs("WRITE_WAVE_FIELD", "NONE")) {
    waveFieldWriteType = None;
  } 
  else if (newOptions.compareArgs("WRITE_WAVE_FIELD", "XDMF")) {
    waveFieldWriteType = Xdmf;
  }
  else if (newOptions.compareArgs("WRITE_WAVE_FIELD", "H5")) {
    waveFieldWriteType = H5;
  }
  else if (newOptions.compareArgs("WRITE_WAVE_FIELD", "H5Compact")) {
    waveFieldWriteType = H5Compact;
  }
  else {
    printf("[WRITE_WAVE_FIELD] not valid [NONE | XDMF | H5 | H5Compact]: defaulting to None\n");
    waveFieldWriteType = None;
  }

  if (waveFieldWriteType != None) {
    if (newOptions.getArgs("TEMPORAL_PPW_OUTPUT", ppwWaveField) == 0) {
      throw std::invalid_argument("[TEMPORAL_PPW_OUTPUT] tag missing");
    }
  }

  std::vector<dfloat> timeStepsWrite;
  if (waveFieldWriteType != None) {
    dfloat dt_write = (1.0/acoustics->fmax)/ppwWaveField;
    int tstepsWrite = std::max((int)floor(dt_write/acoustics->mesh->dt), 1);
    
    for (int tstep = 0; tstep < acoustics->mesh->NtimeSteps; tstep++) {
      if (tstep % tstepsWrite == 0) {
        timeStepsWrite.push_back(tstep * acoustics->mesh->dt);
      }
    }
  }

  std::shared_ptr<IAcousticWriter> writer;

  if (waveFieldWriteType == Xdmf) {
    writer.reset(new AcousticXdmfWriter(acoustics, timeStepsWrite));
  }
  else if (waveFieldWriteType == H5Compact || waveFieldWriteType == H5) {
    auto conn = std::vector<std::vector<uint>>();
    auto x1d = std::vector<dfloat>();
    auto y1d = std::vector<dfloat>();
    auto z1d = std::vector<dfloat>();
    auto p1d = std::vector<dfloat>();

    auto writeConnTable = waveFieldWriteType == H5;
    
    extractUniquePoints(acoustics->mesh, acoustics, conn, x1d, y1d, z1d, p1d);
    writer.reset(new AcousticH5CompactWriter(acoustics, p1d.size(), timeStepsWrite, writeConnTable));
  }
  else {
    throw std::exception();
  }

  if (newOptions.compareArgs("TIME INTEGRATOR", "DOPRI5"))
  {
    dfloat outputInterval;
    newOptions.getArgs("OUTPUT INTERVAL", outputInterval);
    dopri5(acoustics, outputInterval); // NOT WORKING
  }
  else if (newOptions.compareArgs("TIME INTEGRATOR", "LSERK4"))
  {
    lseRK4(acoustics, timeStepsWrite, writer);
  }
  else if (newOptions.compareArgs("TIME INTEGRATOR", "EIRK4"))
  {
    eiRK4(acoustics, timeStepsWrite, writer);
  }
  else if (newOptions.compareArgs("TIME INTEGRATOR", "EIRK4ADAP"))
  {
    eiRK4adap(acoustics);
  }
  else {
    throw std::invalid_argument("TIME INTEGRATOR type not supported. Supported types: DOPRI5, LSERK4, EIRK4, EIRK4ADAP.");
  }
  //timer.toc("Run");

  // [EA] Copy remaining o_qRecv from device to host
  if (acoustics->NReceiversLocal > 0) {
    copyReceiverToHost(acoustics);
  }
}

void updateReceivers(acoustics_t *acoustics) {
  if (acoustics->NReceiversLocal > 0)  {
    acoustics->acousticsReceiverInterpolation(acoustics->NReceiversLocal,
                                      acoustics->o_qRecv,
                                      acoustics->o_recvElements,
                                      acoustics->o_recvElementsIdx,
                                      acoustics->o_recvintpol,
                                      acoustics->o_q,
                                      acoustics->qRecvCounter);
    acoustics->qRecvCounter++;

    if(acoustics->qRecvCounter == recvCopyRate) {
      copyReceiverToHost(acoustics);
      acoustics->qRecvCounter = 0;
      acoustics->qRecvCopyCounter++;  
    }
  }
}

void lseRK4(acoustics_t *acoustics, std::vector<dfloat> timeStepsWrite, std::shared_ptr<IAcousticWriter> writer) {
  mesh_t *mesh = acoustics->mesh;
  
  updateReceivers(acoustics); // update receivers at t=0 (IC)

  int tstepWrite = 0;
  for (int tstep = 0; tstep < mesh->NtimeSteps; ++tstep)
  {
    dfloat time = tstep*mesh->dt;

    if (timeStepsWrite.size() > 0 && abs(timeStepsWrite[tstepWrite] - time) < 1e-10) {
      acoustics->o_q.copyTo(acoustics->q); // copy from GPU to CPU
      writer->write(acoustics, tstepWrite);
      tstepWrite++;
    }
        
    acousticsBCChange(acoustics, time);
    acousticsLserkStep(acoustics, time);
    updateReceivers(acoustics);

    if (tstep % 500 == 0 && !mesh->rank)
    {
      printf("LSERK4 - Step: %d, out of: %d\n", tstep, mesh->NtimeSteps);
    }
  }
}

void eiRK4(acoustics_t *acoustics, std::vector<dfloat> timeStepsWrite, std::shared_ptr<IAcousticWriter> writer) {
  mesh_t *mesh = acoustics->mesh;

  updateReceivers(acoustics); // update receivers at t=0 (IC)

  int tstepWrite = 0;
  for (int tstep = 0; tstep < mesh->NtimeSteps; ++tstep)
  {
    dfloat time = tstep*mesh->dt;

    if (timeStepsWrite.size() > 0 && abs(timeStepsWrite[tstepWrite] - time) < 1e-10) {
      acoustics->o_q.copyTo(acoustics->q); // copy from GPU to CPU
      writer->write(acoustics, tstepWrite);
      tstepWrite++;
    }

    acousticsBCChange(acoustics, time);
    acousticsEirkStep(acoustics, time);
    updateReceivers(acoustics);    

    if (tstep % 500 == 0 && !mesh->rank)
    {
      printf("EIRK4 - Step: %d, out of: %d\n", tstep, mesh->NtimeSteps);
    }
  }    
}

void eiRK4adap(acoustics_t *acoustics) {
  dfloat time = 0.0;
  mesh_t *mesh = acoustics->mesh;

  updateReceivers(acoustics); // update receivers at t=0 (IC)

  dfloat dt = mesh->dt; // Initial dt
  while (time < mesh->finalTime)
  {
    if (time + dt > mesh->finalTime)
    {
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

    dlong accLength = mesh->NLRPoints * acoustics->LRNpoles + mesh->NERPoints * acoustics->ERNpoles;
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

// NOT WORKING (according to comments...)
void dopri5(acoustics_t *acoustics, dfloat outputInterval) {
  mesh_t *mesh = acoustics->mesh;

  // hard code outputInterval for the moment
  dfloat nextOutputTime = outputInterval;
  dfloat outputNumber = 0;

  //initial time
  dfloat time = 0.0;
  int tstep = 0, allStep = 0;

  int done = 0;
  while (!done)
  {

    acoustics->advSwitch = 1;

    if (mesh->dt < acoustics->dtMIN)
    {
      printf("ERROR: Time step became too small at time step=%d\n", tstep);
      exit(-1);
    }
    if (isnan(mesh->dt))
    {
      printf("ERROR: Solution became unstable at time step=%d\n", tstep);
      exit(-1);
    }

    //check for final timestep
    if (time + mesh->dt > mesh->finalTime)
    {
      mesh->dt = mesh->finalTime - time;
      done = 1;
    }

    // try a step with the current time step
    acousticsDopriStep(acoustics, time);

    // compute Dopri estimator
    dfloat err = acousticsDopriEstimate(acoustics);

    // build controller
    dfloat fac1 = pow(err, acoustics->exp1);
    dfloat fac = fac1 / pow(acoustics->facold, acoustics->beta);

    fac = mymax(acoustics->invfactor2, mymin(acoustics->invfactor1, fac / acoustics->safe));
    dfloat dtnew = mesh->dt / fac;

    if (err < 1.0)
    { //dt is accepted

      // check for output during this step and do a mini-step
      if (time < nextOutputTime && time + mesh->dt > nextOutputTime)
      {
        dfloat savedt = mesh->dt;

        // save rkq
        acoustics->o_saveq.copyFrom(acoustics->o_rkq);

        // change dt to match output
        mesh->dt = nextOutputTime - time;

        // print
        printf("Taking output mini step: %g\n", mesh->dt);

        // time step to output
        acousticsDopriStep(acoustics, time);

        // shift for output
        acoustics->o_rkq.copyTo(acoustics->o_q);

        // restore time step
        mesh->dt = savedt;

        // increment next output time
        nextOutputTime += outputInterval;

        // accept saved rkq
        acoustics->o_q.copyFrom(acoustics->o_saveq);
      }
      else
      {
        // accept rkq
        acoustics->o_q.copyFrom(acoustics->o_rkq);
      }

      time += mesh->dt;

      acoustics->facold = mymax(err, 1E-4); // hard coded factor ?

      printf("\r time = %g (%d), dt = %g accepted                      ", time, allStep, mesh->dt);
      tstep++;
    }
    else
    {
      dtnew = mesh->dt / (mymax(acoustics->invfactor1, fac1 / acoustics->safe));
      printf("\r time = %g (%d), dt = %g rejected, trying %g", time, allStep, mesh->dt, dtnew);

      done = 0;
    }
    mesh->dt = dtnew;
    allStep++;
  }

  mesh->device.finish();
}