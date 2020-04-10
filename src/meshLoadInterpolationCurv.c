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

#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"

void interpolateToGL(mesh3D *mesh){
	dfloat *intpol = mesh->GLIntpolCurv;
	int cnt = 0;
	for(int e = 0; e < mesh->Nelements; e++){
		for(int n = 0; n < mesh->Np; n++){
			int id = e*mesh->Np + n; 
			for(int i = 0; i < mesh->Np; i++){
				mesh->x[id] += mesh->EXCurv[e*mesh->Np + i]*intpol[n*mesh->Np + i];
				mesh->y[id] += mesh->EYCurv[e*mesh->Np + i]*intpol[n*mesh->Np + i];
				mesh->z[id] += mesh->EZCurv[e*mesh->Np + i]*intpol[n*mesh->Np + i];
				
			}
		}
	}
}

void meshLoadInterpolationCurv(mesh3D *mesh, int N){

  char fname[BUFSIZ];
  sprintf(fname, DHOLMES "/nodes/gmshToGLInterpolationMatrix.dat");

  FILE *fp = fopen(fname, "r");

  if (!fp) {
    printf("ERROR: Cannot open file: '%s'\n", fname);
    exit(-1);
  }

  mesh->GLIntpolCurv = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));


	int Nrows, Ncols;
	char mname[BUFSIZ];
	sprintf(mname,"Interpolate gmsh to GL N%d",N);	
  readDfloatArray(fp, mname, &(mesh->GLIntpolCurv),&Nrows,&Ncols);
	fclose(fp);
  mesh->x = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->y = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->z = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
	
	interpolateToGL(mesh);

	free(mesh->GLIntpolCurv);
	free(mesh->EXCurv);
	free(mesh->EYCurv);
	free(mesh->EZCurv);

}


