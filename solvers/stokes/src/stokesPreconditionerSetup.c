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

#include "stokes.h"

static void stokesJacobiPreconditionerSetup(stokes_t *stokes);
static void stokesSCBlockPreconditionerSetup(stokes_t *stokes, occa::properties &kernelInfo);

static void stokesBuildLocalContinuousDiagQuad2D(stokes_t* stokes, dlong e, dfloat *diagA);

void stokesPreconditionerSetup(stokes_t *stokes, occa::properties &kernelInfo)
{
  if (stokes->options.compareArgs("PRECONDITIONER", "NONE")) {
    stokes->precon = NULL;
    return;
  } else if (stokes->options.compareArgs("PRECONDITIONER", "JACOBI")) {
    stokesJacobiPreconditionerSetup(stokes);
  } else if (stokes->options.compareArgs("PRECONDITIONER", "SCBLOCK")) {
    stokesSCBlockPreconditionerSetup(stokes, kernelInfo);
  } else {
    printf("ERROR:  Invalid value %s for [PRECONDITIONER] option.\n",
           stokes->options.getArgs("PRECONDITIONER").c_str());
    exit(-1);
  }

  return;
}

/*****************************************************************************/

static void stokesJacobiPreconditionerSetup(stokes_t *stokes)
{
  int elementType;

  stokes->precon = new stokesPrecon_t();

  stokes->options.getArgs("JACOBI BOOST PARAMETER", stokes->precon->boost);
  if (stokes->precon->boost == 0.0) {
    printf("ERROR:  Jacobi boost must be nonzero.\n");
    exit(-1);
  }

  stokesVecAllocate(stokes, &stokes->precon->invDiagA);

  switch (stokes->elementType) {
  case QUADRILATERALS:
    if (stokes->mesh->dim == 2) {
      for (dlong e = 0; e < stokes->mesh->Nelements; e++) {
        /* TODO:  This is wasteful---should just build once and then copy. */
        stokesBuildLocalContinuousDiagQuad2D(stokes, e, stokes->precon->invDiagA.x + e*stokes->mesh->Np);
        stokesBuildLocalContinuousDiagQuad2D(stokes, e, stokes->precon->invDiagA.y + e*stokes->mesh->Np);
      }

      for (int i = 0; i < stokes->Ntotal; i++){
        stokes->precon->invDiagA.p[i] = stokes->precon->boost;
      }
    } else {
      printf("ERROR:  Not implemented.\n");
      exit(-1);
    }
    break;
  default:
    printf("ERROR:  Not implemented.\n");
    exit(-1);
  }

  // BADNESS
  for (int i = 0; i < stokes->Ndof; i++)
    stokes->precon->invDiagA.v[i] = 1.0/stokes->precon->invDiagA.v[i];
  
  stokesVecCopyHostToDevice(stokes->precon->invDiagA);

  if (stokes->options.compareArgs("VELOCITY DISCRETIZATION", "CONTINUOUS")) {
    stokesVecGatherScatter(stokes, stokes->precon->invDiagA);
  }

  return;
}

static void stokesSCBlockPreconditionerSetup(stokes_t *stokes, occa::properties &kernelInfo)
{
  mesh_t     *mesh;
  elliptic_t *elliptic;
  setupAide  ellipticOptions;

  mesh = stokes->mesh;
  elliptic = new elliptic_t();

  /* Set up the elliptic sub-solver. */
  elliptic->mesh = mesh;
  elliptic->dim = elliptic->mesh->dim;
  elliptic->elementType = stokes->elementType;

  /* TODO:  Map these to the Stokes setup file. */
  ellipticOptions.setArgs("BASIS", "NODAL");
  ellipticOptions.setArgs("DISCRETIZATION", "CONTINUOUS");
  ellipticOptions.setArgs("DEBUG ENABLE OGS", "1");
  ellipticOptions.setArgs("DEBUG ENABLE REDUCTIONS", "1");
  ellipticOptions.setArgs("KRYLOV SOLVER", "PCG");
  ellipticOptions.setArgs("PRECONDITIONER", "MULTIGRID");
  ellipticOptions.setArgs("MULTIGRID COARSENING", "HALFDEGREES");
  ellipticOptions.setArgs("MULTIGRID SMOOTHER", "DAMPEDJACOBI+CHEBYSHEV");
  ellipticOptions.setArgs("MULTIGRID CHEBYSHEV DEGREE", "2");
  ellipticOptions.setArgs("PARALMOND CYCLE", "VCYCLE");
  ellipticOptions.setArgs("PARALMOND SMOOTHER", "CHEBYSHEV+DAMPEDJACOBI");
  ellipticOptions.setArgs("PARALMOND CHEBYSHEV DEGREE", "2");

  ellipticOptions.setArgs("VERBOSE", "FALSE");
  elliptic->options = ellipticOptions;

  elliptic->BCType = stokes->BCType;

  elliptic->r = (dfloat*)calloc(stokes->Ntotal, sizeof(dfloat));
  elliptic->o_r = mesh->device.malloc(stokes->Ntotal*sizeof(dfloat), elliptic->r);
  elliptic->x = (dfloat*)calloc(stokes->Ntotal, sizeof(dfloat));
  elliptic->o_x = mesh->device.malloc(stokes->Ntotal*sizeof(dfloat), elliptic->x);

  /* TODO:  This allocates a whole lot of extra stuff---we may be able to share
   * some scratch arrays with the stokes_t.
   */
  ellipticSolveSetup(elliptic, 0.0, kernelInfo);

  stokes->precon = new stokesPrecon_t();
  stokes->precon->elliptic = elliptic;

  stokesVecAllocate(stokes, &stokes->precon->invMM);

  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      stokes->precon->invMM.p[e*mesh->Np+n] = mesh->ggeo[e*mesh->Np*mesh->Nggeo+GWJID*mesh->Np+n];
    }
  }

#if 1
  // assemble mass matrix (suitable for C0 pressure)
  stokesVecCopyHostToDevice(stokes->precon->invMM);
  
  stokesVecGatherScatter(stokes, stokes->precon->invMM);
  
  stokesVecCopyDeviceToHost(stokes->precon->invMM);
#endif
  
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat val = stokes->precon->invMM.p[e*mesh->Np+n];
      if(val)
	stokes->precon->invMM.p[e*mesh->Np+n] = 1./val;
    }
  }
  
  stokesVecCopyHostToDevice(stokes->precon->invMM);
   
  return;
}

/*****************************************************************************/

/* TODO:  This was basically copied from the elliptic solver.
 *
 * TODO:  This may need modification for boundary conditions.
 */
static void stokesBuildLocalContinuousDiagQuad2D(stokes_t* stokes, dlong e, dfloat *diagA)
{
  mesh_t *mesh = stokes->mesh;

  for (int ny=0;ny<mesh->Nq;ny++) {
    for (int nx=0;nx<mesh->Nq;nx++) {
      int iid = nx+ny*mesh->Nq;
      diagA[iid] = 0;

      for (int k=0;k<mesh->Nq;k++) {
        int id = k+ny*mesh->Nq;
        dfloat eta = stokes->eta[e*mesh->Np + id];
        dfloat Grr = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G00ID*mesh->Np];
        diagA[iid] += Grr*mesh->D[nx+k*mesh->Nq]*mesh->D[nx+k*mesh->Nq]*eta;
      }

      for (int k=0;k<mesh->Nq;k++) {
        int id = nx+k*mesh->Nq;
        dfloat eta = stokes->eta[e*mesh->Np + id];
        dfloat Gss = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G11ID*mesh->Np];
        diagA[iid] += Gss*mesh->D[ny+k*mesh->Nq]*mesh->D[ny+k*mesh->Nq]*eta;
      }

      int id = nx+ny*mesh->Nq;
      dfloat eta = stokes->eta[e*mesh->Np + id];
      dfloat Grs = mesh->ggeo[e*mesh->Np*mesh->Nggeo + id + G01ID*mesh->Np];
      diagA[iid] += 2*Grs*mesh->D[nx+nx*mesh->Nq]*mesh->D[ny+ny*mesh->Nq]*eta;
    }
  }

  return;
}