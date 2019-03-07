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

static void stokesTestSolutionConstantViscosityQuad2D(dfloat x, dfloat y, dfloat *ux, dfloat *uy);
static void stokesTestSolutionVariableViscosityQuad2D(dfloat x, dfloat y, dfloat *ux, dfloat *uy);
static void stokesTestSolutionDirichletQuad2D(dfloat x, dfloat y, dfloat *ux, dfloat *uy);
static void stokesTestSolutionConstantViscosityHex3D(dfloat x, dfloat y, dfloat z, dfloat *ux, dfloat *uy, dfloat *uz);

int main(int argc, char **argv)
{
  stokes_t         *stokes;
  occa::properties kernelInfo;

  // Start up MPI.
  MPI_Init(&argc, &argv);

  if (argc != 2) {
    printf("usage: ./stokesMain setupfile\n");
    MPI_Finalize();
    exit(-1);
  }

  setupAide options(argv[1]);

  stokes = stokesSetup(kernelInfo, options);
  stokesSolve(stokes);

  stokesVecCopyDeviceToHost(stokes->u);

#if 1
  /* Compute error (if applicable.) */
  dfloat errxInf = 0.0, erryInf = 0.0, errzInf = 0.0;
  dfloat errxDL2 = 0.0, erryDL2 = 0.0, errzDL2 = 0.0;
  for (int e = 0; e < stokes->mesh->Nelements; e++) {
    for (int i = 0; i < stokes->mesh->Np; i++) {
      int    ind;
      dfloat x, y, z;
      dfloat errx, erry, errz;
      dfloat ux_exact, uy_exact, uz_exact;

      ind = e*stokes->mesh->Np + i;
      x = stokes->mesh->x[ind];
      y = stokes->mesh->y[ind];
      z = stokes->mesh->z[ind];

      if (stokes->mesh->dim == 2) {
        //stokesTestSolutionConstantViscosityQuad2D(x, y, &ux_exact, &uy_exact);
        //stokesTestSolutionVariableViscosityQuad2D(x, y, &ux_exact, &uy_exact);
        if (stokes->mapB[e*stokes->mesh->Np + i] == 1) {
          stokes->u.x[ind] = cos(y);  // Manually insert the boundary data.
          stokes->u.y[ind] = sin(x);

          ux_exact = cos(y);
          uy_exact = sin(x);
        } else {
          stokesTestSolutionDirichletQuad2D(x, y, &ux_exact, &uy_exact);
        }
      } else if (stokes->mesh->dim == 3) {
        stokesTestSolutionConstantViscosityHex3D(x, y, z, &ux_exact, &uy_exact, &uz_exact);
      }

      errx = stokes->u.x[ind] - ux_exact;
      erry = stokes->u.y[ind] - uy_exact;
      if (stokes->mesh->dim == 3)
        errz = stokes->u.z[ind] - uz_exact;

      if (fabs(errx) > errxInf)
        errxInf = fabs(errx);
      if (fabs(erry) > erryInf)
        erryInf = fabs(erry);
      errxDL2 += errx*errx;
      erryDL2 += erry*erry;

      if (stokes->mesh->dim == 3) {
        if (fabs(errz) > errzInf)
          errzInf = fabs(errz);
        errzDL2 += errz*errz;
      }
    }
  }

  errxDL2 = sqrt(errxDL2);
  erryDL2 = sqrt(erryDL2);
  if (stokes->mesh->dim == 3)
    errzDL2 = sqrt(errzDL2);

  printf("-----\n");

  printf("errxInf = % .15e\n", errxInf);
  printf("erryInf = % .15e\n", erryInf);
  if (stokes->mesh->dim == 3)
    printf("errzInf = % .15e\n", errzInf);
  printf("errxDL2 = % .15e\n", errxDL2);
  printf("erryDL2 = % .15e\n", erryDL2);
  if (stokes->mesh->dim == 3)
    printf("errzDL2 = % .15e\n", errzDL2);
#endif

#if 0
  printf("-----\n");

  printf("u = [");
  stokesVecPrint(stokes, stokes->u);
  printf("];\n");

  printf("x = [");
  for (int i = 0; i < stokes->Ntotal; i++) {
    printf("% .15e\n", stokes->mesh->x[i]);
  }
  printf("];\n");

  printf("y = [");
  for (int i = 0; i < stokes->NtotalV; i++) {
    printf("% .15e\n", stokes->mesh->y[i]);
  }
  printf("];\n");

  printf("Ntotal = %d\n", stokes->Ntotal);
#endif

  /* Export solution. */

  /* Report runtime statistics. */

  // Shut down MPI.
  MPI_Finalize();

  return 0;
}

/*****************************************************************************/

static void stokesTestSolutionConstantViscosityQuad2D(dfloat x, dfloat y, dfloat *ux, dfloat *uy)
{
  *ux = 6.0*pow(1.0 - x*x, 3.0)*pow(1.0 - y*y, 2.0)*y;
  *uy = -6.0*pow(1.0 - y*y, 3.0)*pow(1.0 - x*x, 2.0)*x;
  return;
}

static void stokesTestSolutionVariableViscosityQuad2D(dfloat x, dfloat y, dfloat *ux, dfloat *uy)
{
  *ux = 6.0*pow(1.0 - x*x, 3.0)*pow(1.0 - y*y, 2.0)*y;
  *uy = -6.0*pow(1.0 - y*y, 3.0)*pow(1.0 - x*x, 2.0)*x;
  return;
}

static void stokesTestSolutionDirichletQuad2D(dfloat x, dfloat y, dfloat *ux, dfloat *uy)
{
  *ux = cos(y);
  *uy = sin(x);
  return;
}

static void stokesTestSolutionConstantViscosityHex3D(dfloat x, dfloat y, dfloat z, dfloat *ux, dfloat *uy, dfloat *uz)
{
  *ux = -6.0*z*pow(1.0 - z*z, 2.0);
  *uy = -6.0*x*pow(1.0 - x*x, 2.0);
  *uz = -6.0*y*pow(1.0 - y*y, 2.0);
  return;
}
