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

void acousticsReport(acoustics_t *acoustics, dfloat time, setupAide &newOptions){

  mesh3D *mesh = acoustics->mesh;

  // copy data back to host
  acoustics->o_q.copyTo(acoustics->q);

// Print output to txt files
#if 1
  int nprocs, procid;
  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  hlong globalNelements;
  MPI_Reduce(&mesh->Nelements, &globalNelements, 1, MPI_HLONG, MPI_SUM, 0, mesh->comm);
  
  FILE * xFP, * yFP, * zFP, * pFP;
  char xFN[50] = "data/x.txt";
  char yFN[50] = "data/y.txt";
  char zFN[50] = "data/z.txt";
  char pFN[50] = "data/p.txt";  
  for(int i = 0; i < nprocs; i++){
    if(procid == i){
      if(procid == 0){
          xFP = fopen(xFN, "w");
          yFP = fopen(yFN, "w");
          zFP = fopen(zFN, "w");
          pFP = fopen(pFN, "w");

          // Write Np and Nelements as first line in x output
          fprintf(xFP, "%d %d\n",mesh->Np, globalNelements);
      } else {
          xFP = fopen(xFN, "a");
          yFP = fopen(yFN, "a");
          zFP = fopen(zFN, "a");
          pFP = fopen(pFN, "a");
      }
      for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat xEA = mesh->x[n + mesh->Np*e];
      dfloat yEA = mesh->y[n + mesh->Np*e];
      dfloat zEA = mesh->z[n + mesh->Np*e];

      dlong qbaseEA = e*mesh->Np*mesh->Nfields + n;

      dfloat rEA = acoustics->q[qbaseEA+0*mesh->Np];

      fprintf(xFP, "%.15lf ",xEA);
      fprintf(yFP, "%.15lf ",yEA);
      fprintf(zFP, "%.15lf ",zEA);
      fprintf(pFP, "%.15lf ",rEA);
        }
        fprintf(xFP, "\n");
        fprintf(yFP, "\n");
        fprintf(zFP, "\n");
        fprintf(pFP, "\n");
      }

      fclose(xFP);
      fclose(yFP);
      fclose(zFP);
      fclose(pFP);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif


  // do error stuff on host
  acousticsError(acoustics, time);

  // output field files
  char fname[BUFSIZ];

  sprintf(fname, "foo_%04d_%04d.vtu", mesh->rank, acoustics->frame++);

  acousticsPlotVTU(acoustics, fname);
  
}
