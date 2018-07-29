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

#include "ellipticQuad2D.h"

void BuildLocalIpdgPatchAx(mesh2D *mesh, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *B, dfloat *Br, dfloat* Bs, dlong eM, dfloat *A);

void BuildLocalContinuousPatchAx(solver_t* solver, dfloat lambda,
                                  dlong eM, dfloat *B, dfloat *Br, dfloat* Bs, dfloat *A);

void ellipticBuildJacobiQuad2D(solver_t* solver, dfloat tau, dfloat lambda,
                                   int *BCType, dfloat **invDiagA,
                                   const char *options){

  mesh2D *mesh = solver->mesh;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  // build some monolithic basis arrays
  dfloat *B  = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Br = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bs = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

  int mode = 0;
  for(int nj=0;nj<mesh->N+1;++nj){
    for(int ni=0;ni<mesh->N+1;++ni){

      int node = 0;

      for(int j=0;j<mesh->N+1;++j){
        for(int i=0;i<mesh->N+1;++i){

          if(nj==j && ni==i)
            B[mode*mesh->Np+node] = 1;
          if(nj==j)
            Br[mode*mesh->Np+node] = mesh->D[ni+mesh->Nq*i]; 
          if(ni==i)
            Bs[mode*mesh->Np+node] = mesh->D[nj+mesh->Nq*j]; 
          
          ++node;
        }
      }
      ++mode;
    }
  }

  dlong diagNnum = mesh->Np*mesh->Nelements;

  dfloat *diagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));

  if(rank==0) printf("Building diagonal...");fflush(stdout);

  // loop over all elements
  #pragma omp parallel
  {
    //temp patch storage
    dfloat *patchA = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

    #pragma omp for 
    for(dlong eM=0;eM<mesh->Nelements;++eM){
      //build the patch A matrix for this element
      if (strstr(options,"IPDG")) {
        BuildLocalIpdgPatchAx(mesh, tau, lambda, BCType, B, Br, Bs, eM, patchA);
      } else if (strstr(options,"CONTINUOUS")) {
        BuildLocalContinuousPatchAx(solver, lambda, eM, B, Br, Bs, patchA);
      }

      // compute the diagonal entries
      for(int n=0;n<mesh->Np;++n){
        diagA[eM*mesh->Np + n] = patchA[n*mesh->Np+n]; //store the inverse diagonal entry
      }
    }

    free(patchA);
  }

  if (strstr(options,"CONTINUOUS")) 
    gsParallelGatherScatter(mesh->hostGsh, diagA, dfloatString, "add"); 
    
  *invDiagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));
  for (dlong n=0;n<mesh->Nelements*mesh->Np;n++) {
    (*invDiagA)[n] = 1/diagA[n];
  }

  if(rank==0) printf("done.\n");

  free(B); free(Br); free(Bs);
  free(diagA);
}
