#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh3D.h"

typedef struct {

  occa::memory o_vmapPP;
  occa::memory o_faceNodesP;

  occa::memory o_oasForward;
  occa::memory o_oasBack;
  occa::memory o_oasDiagInvOp;
  occa::memory o_invDegreeP;

  occa::memory o_oasForwardDg;
  occa::memory o_oasBackDg;
  occa::memory o_oasDiagInvOpDg;
  occa::memory o_invDegreeDGP;
  
  occa::kernel restrictKernel;
  occa::kernel preconKernel;

  occa::kernel coarsenKernel;
  occa::kernel prolongateKernel;  

  ogs_t *ogsP, *ogsDg;

  occa::memory o_diagA;

  // coarse grid basis for preconditioning
  occa::memory o_V1, o_Vr1, o_Vs1, o_Vt1;
  occa::memory o_r1, o_z1;
  dfloat *r1, *z1;
  void *xxt, *amg, *almond;

  occa::memory o_coarseInvDegree;
  occa::memory o_ztmp;

  iint coarseNp;
  iint coarseTotal;
  iint *coarseOffsets;
  dfloat *B, *tmp2;
  occa::memory *o_B, o_tmp2;
  void *xxt2;
  void *parAlmond;

  
} precon_t;

void ellipticRunHex3D(mesh3D *mesh);

void ellipticOccaRunHex3D(mesh3D *mesh);

void ellipticSetupHex3D(mesh3D *mesh, occa::kernelInfo &kernelInfo);

void ellipticVolumeHex3D(mesh3D *mesh);

void ellipticSurfaceHex3D(mesh3D *mesh, dfloat time);

void ellipticUpdateHex3D(mesh3D *mesh, dfloat rka, dfloat rkb);

void ellipticErrorHex3D(mesh3D *mesh, dfloat time);

void ellipticParallelGatherScatter(mesh3D *mesh, ogs_t *ogs, occa::memory &o_v, occa::memory &o_gsv,
				    const char *type, const char *op);

precon_t *ellipticPreconditionerSetupHex3D(mesh3D *mesh, ogs_t *ogs, dfloat lambda, const char *options);

void ellipticCoarsePreconditionerHex3D(mesh_t *mesh, precon_t *precon, dfloat *x, dfloat *b);

void ellipticCoarsePreconditionerSetupHex3D(mesh_t *mesh, precon_t *precon, ogs_t *ogs, dfloat lambda, const char *options);

typedef struct {

  mesh_t *mesh;

  precon_t *precon;

  ogs_t *ogs;

  ogs_t *ogsDg;

  // C0 halo gather-scatter info
  ogs_t *halo;

  // C0 nonhalo gather-scatter info
  ogs_t *nonHalo;
  
  
  iint Nblock;
  
  occa::memory o_p; // search direction
  occa::memory o_z; // preconditioner solution
  occa::memory o_zP; // extended OAS preconditioner patch solution
  occa::memory o_Ax; // A*initial guess
  occa::memory o_Ap; // A*search direction
  occa::memory o_tmp; // temporary
  occa::memory o_grad; // temporary gradient storage (part of A*)
  occa::memory o_rtmp;
  occa::memory o_invDegree;
  occa::memory o_pAp;

  dfloat *sendBuffer, *recvBuffer;

  // HOST shadow copies
  dfloat *Ax, *p, *r, *z, *zP, *Ap, *tmp, *grad;

  // integration storage for BP3
  int gNq;
  occa::memory o_gggeo; // Jacobian matrix at integration nodes
  occa::memory o_gI;    // interpolate from GLL to integration nodes
  occa::memory o_gD;    // differentiate and interpolate from GLL to integration nodes

  // list of elements that are needed for global gather-scatter
  iint NglobalGatherElements;
  iint *globalGatherElementList;
  occa::memory o_globalGatherElementList;

  // list of elements that are not needed for global gather-scatter
  iint NlocalGatherElements;
  iint *localGatherElementList;
  occa::memory o_localGatherElementList;
  
  occa::kernel AxKernel;
  occa::kernel partialAxKernel;
  
  occa::kernel gradientKernel;
  occa::kernel partialGradientKernel;

  occa::kernel ipdgKernel;
  occa::kernel partialIpdgKernel;
  
  occa::stream defaultStream;
  occa::stream dataStream;

  
  
}solver_t;

// block size for reduction (hard coded)
#define blockSize 256

void ellipticMatrixFreeAx(void **args, occa::memory o_q, occa::memory o_Aq, const char* options);

int ellipticSolveHex3D(solver_t *solver, dfloat lambda, occa::memory &o_r, occa::memory &o_x, iint maxIterations, const char *options);

solver_t *ellipticSolveSetupHex3D(mesh_t *mesh, dfloat lambda, occa::kernelInfo &kernelInfo, const char *options);


void ellipticStartHaloExchange3D(solver_t *solver, occa::memory &o_q, dfloat *sendBuffer, dfloat *recvBuffer);

void ellipticInterimHaloExchange3D(solver_t *solver, occa::memory &o_q, dfloat *sendBuffer, dfloat *recvBuffer);

void ellipticEndHaloExchange3D(solver_t *solver, occa::memory &o_q, dfloat *recvBuffer);

void ellipticParallelGatherScatterHex3D(mesh3D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type, const char *op);

void ellipticHaloGatherScatter(solver_t *solver, 
			       ogs_t *halo, 
			       occa::memory &o_v,
			       const char *type,
			       const char *op);

void ellipticNonHaloGatherScatter(solver_t *solver, 
				  ogs_t *nonHalo, 
				  occa::memory &o_v,
				  const char *type,
				  const char *op);


void ellipticParallelGatherScatterSetup(mesh_t *mesh,    // provides DEVICE
					iint Nlocal,     // number of local nodes
					iint Nbytes,     // number of bytes per node
					iint *gatherLocalIds,  // local index of nodes
					iint *gatherBaseIds,   // global index of their base nodes
					iint *gatherHaloFlags,
					ogs_t **halo,
					ogs_t **nonHalo);   // 1 for halo node, 0 for not
