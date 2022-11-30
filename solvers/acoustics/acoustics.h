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

#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"

#include "mesh.h"
#include "mesh2D.h"
#include "mesh3D.h"

// block size for reduction (hard coded)
#define blockSize 256
// include Gaussian random fields
#define INCLUDE_GRF 0
#define AMPLITUDE 2
#define OFFSET_BC 2.7 // only used when no IC mesh is loaded

#define TRIANGLES 3
#define QUADRILATERALS 4
#define TETRAHEDRA 6
#define HEXAHEDRA 12

enum WriteWaveFieldType { None, Vtu, Xdmf, Txt };
enum SourceType { GaussianFunction, GRF };

typedef struct {  
  string outDir;
  string simulationID;
  SourceType sourceType;
  
  dfloat fmax;   // max frequency
  dfloat sigma0; // source width
  dfloat sourcePosition[3];
  
  int dim;
  int elementType; // number of edges (3=tri, 4=quad, 6=tet, 12=hex)
  
  int Nfields;

  hlong totalElements;
  dlong Nblock;

  dfloat *q, *rhsq, *resq;

  dfloat xminmax[2], yminmax[2], zminmax[2];
  
  std::vector<dfloat> ic_rectilinear;
  std::vector<dfloat> x1d_rectilinear;
  std::vector<dfloat> y1d_rectilinear;
  std::vector<dfloat> z1d_rectilinear;
  dfloat ic_rectilinear_shape[3];
  
  //---------RECEIVER---------
  dfloat *qRecv; // Saves pres in receiver element in each timestep
  dlong qRecvCounter; // To keep track of which timestep we are on
  dlong qRecvCopyCounter;
  dlong *recvElements; // Index to elements where the receivers are located
  dlong *recvElementsIdx; // Index into recvElements
  dfloat *recvXYZ; // XYZ coordinates of receiver
  dlong NReceivers; // Total number of receivers
  dlong NReceiversLocal; // Total number of receivers

  occa::memory o_recvElements;
  occa::memory o_recvElementsIdx;
  //---------RECEIVER---------

  dfloat *acc;
  dfloat *rhsacc;
  dfloat *resacc;
  dlong LRNpoles;
  dlong LRNRealPoles;
  dlong LRNImagPoles; // Number of complex conjugate pairs
  dlong NBoundaryPoints;
  dlong ERNpoles;
  dlong ERNRealPoles;
  dlong ERNImagPoles;

  dfloat BCChangeTime;
  dfloat *LR;
  dlong *LRInfo;
  occa::memory o_LR;
  occa::memory o_LRInfo;
  dfloat *ER;
  dlong *ERInfo;
  occa::memory o_ER;
  occa::memory o_ERInfo;


  dlong NERComPoints;
  dfloat *ERComPoints;
  dlong NERPointsTotal;
  dlong *ERComPointsIdx;
  dlong *ERintpolElementsCom;
  dlong comPointsCounter; // Unused?
  dlong NComPointsToSendAllRanks;
  dfloat *vtSend;
  dfloat *vtRecv;
  dlong *recvCountsArray;
  

  occa::memory o_acc;
  occa::memory o_rhsacc;
  occa::memory o_resacc;
  occa::memory o_ERintpol;
  occa::memory o_recvintpol;
  occa::memory o_ERintpolElements;
  occa::memory o_ERintpolElementsCom;
  occa::memory o_ERintpolCom;
  occa::memory o_vtSend;
  occa::memory o_vtRecv;
  occa::memory o_ERComPointsIdx;
  occa::memory o_recvCountsArray;
  occa::memory o_comPointsIdxAll;
  occa::memory o_comPointsToSend;
  occa::memory o_vt;
  occa::memory o_vi;
  occa::memory o_anglei;


  // EIRK4 storage
  dfloat *k1acc;
  dfloat *k2acc;
  dfloat *k3acc;
  dfloat *k4acc;
  dfloat *k5acc;
  dfloat *k6acc;

  dfloat *Xacc;

  dfloat *k1rhsq;
  dfloat *k2rhsq;
  dfloat *k3rhsq;
  dfloat *k4rhsq;
  dfloat *k5rhsq;
  dfloat *k6rhsq;

  occa::memory o_k1acc;
  occa::memory o_k2acc;
  occa::memory o_k3acc;
  occa::memory o_k4acc;
  occa::memory o_k5acc;
  occa::memory o_k6acc;

  occa::memory o_Xacc;

  occa::memory o_k1rhsq;
  occa::memory o_k2rhsq;
  occa::memory o_k3rhsq;
  occa::memory o_k4rhsq;
  occa::memory o_k5rhsq;
  occa::memory o_k6rhsq;

  dfloat *Vort;

  dfloat *rkq, *rkrhsq, *rkerr;
  dfloat *errtmp;
  int frame;

  mesh_t *mesh;

  occa::kernel volumeKernel;
  occa::kernel surfaceKernel;
  occa::kernel updateKernel;
  occa::kernel rkStageKernel;
  occa::kernel rkUpdateKernel;
  occa::kernel rkErrorEstimateKernel;
  occa::kernel receiverKernel;
  occa::kernel updateKernelLR;
  occa::kernel updateKernelER;
  occa::kernel acousticsUpdateEIRK4;
  occa::kernel acousticsUpdateEIRK4AccLR;
  occa::kernel acousticsUpdateEIRK4AccER;
  occa::kernel ERangleDetection;
  occa::kernel ERMoveVT;
  occa::kernel ERInsertComVT;
  occa::kernel acousticsWSComInterpolation;
  occa::kernel acousticsReceiverInterpolation;
  occa::kernel volumeKernelCurv;
  occa::kernel surfaceKernelCurv;

  occa::memory o_q;
  occa::memory o_rhsq;
  occa::memory o_resq;
  occa::memory o_saveq;
  
  //[EA] 
  occa::memory o_qRecv;
  occa::kernel acousticsErrorEIRK4;
  occa::kernel acousticsErrorEIRK4r;
  occa::kernel acousticsErrorEIRK4Acc;
  occa::kernel acousticsErrorEIRK4Accr;
  dfloat *rkAcc, *rkerrAcc;
  occa::memory o_rkAcc, o_rkerrAcc;
  occa::memory o_rkq, o_rkrhsq, o_rkerr;
  occa::memory o_errtmp;
  
  //halo data
  dlong haloBytes;
  dfloat *sendBuffer;
  dfloat *recvBuffer;
  occa::memory o_sendBuffer;
  occa::memory o_recvBuffer;
  occa::memory o_haloBuffer;

  // DOPRI5 RK data
  int advSwitch;
  int Nrk;
  dfloat ATOL, RTOL;
  dfloat factor1, invfactor1;
  dfloat factor2, invfactor2;
  dfloat exp1, facold,  dtMIN, safe, beta;
  dfloat *rkA, *rkC, *rkE;
  occa::memory o_rkA, o_rkC, o_rkE;
} acoustics_t;

int acousticsSetupMain(setupAide &newOptions);

void acousticsRun(acoustics_t *acoustics, setupAide &newOptions);

acoustics_t *acousticsSetup(mesh_t *mesh, setupAide &newOptions, char* boundaryHeaderFileName);
void sourceSetup(mesh_t *mesh, acoustics_t *acoustics, setupAide &newOptions);

void acousticsError(acoustics_t *acoustics, dfloat time);

void acousticsCavitySolution(dfloat x, dfloat y, dfloat z, dfloat t,
		       dfloat *u, dfloat *v, dfloat *w, dfloat *p);

dfloat gaussianSource(dfloat x, dfloat y, dfloat z, dfloat t, dfloat *sloc, dfloat sxyz, dfloat amplitude = AMPLITUDE);
void gaussianSource(vector<dfloat> x1d, vector<dfloat> y1d, vector<dfloat> z1d, dfloat *sloc, dfloat sxyz, 
    vector<dfloat> &pressures, dfloat ampl = AMPLITUDE);
#if INCLUDE_GRF
void grfWindowed(vector<dfloat> x1d, vector<dfloat> y1d, vector<dfloat> z1d, 
    dfloat xminmax[2], dfloat yminmax[2], dfloat zminmax[2], 
    dfloat sigma_0, dfloat l_0, dfloat sigma0_window, vector<dfloat> &samples_out, dfloat amplitude = AMPLITUDE);
#endif

void acousticsWritePressureField(acoustics_t *acoustics, WriteWaveFieldType waveFieldWriteType, std::vector<dfloat> timeSteps, int iter);
void acousticsWriteXdmf(acoustics_t *acoustics, std::vector<dfloat> timeSteps, int iter);
void acousticsWriteVTU(acoustics_t *acoustics, bool writeVelocity = false);
void acousticsWritePressureFieldTxt(acoustics_t *acoustics, dfloat time);
int acousticsWriteIRs(acoustics_t *acoustics, setupAide &newOptions);
void acousticsWriteSimulationSettings(acoustics_t *acoustics, string filename);

void acousticsDopriStep(acoustics_t *acoustics, const dfloat time);
void acousticsLserkStep(acoustics_t *acoustics, const dfloat time);
void acousticsEirkStep(acoustics_t *acoustics, const dfloat time);
dfloat acousticsDopriEstimate(acoustics_t *acoustics);

void acousticsReceiverInterpolation(acoustics_t *acoustics);
void acousticsFindReceiverElement(acoustics_t *acoustics);
void acousticsRecvIntpolOperators(acoustics_t *acoustics);

int createDir(string path, bool deleteIfExists);
