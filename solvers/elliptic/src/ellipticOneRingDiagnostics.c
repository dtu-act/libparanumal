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

#include "elliptic.h"

void ellipticOneRingDiagnostics(elliptic_t *elliptic, elliptic_t *elliptic1){

  mesh_t *mesh  = elliptic->mesh;
  mesh_t *mesh1 = elliptic1->mesh;

  char fname[BUFSIZ];
  sprintf(fname, "diagnostics%04d.dat", mesh->rank);
  FILE *fp = fopen(fname, "w");
  fprintf(fp, "EToV=[\n");
  for(int e=0;e<mesh1->Nelements;++e){
    for(int v=0;v<mesh1->Nverts;++v)
      fprintf(fp, "%d ", mesh1->EToV[e*mesh1->Nverts+v]);
    if(e<mesh->Nelements) fprintf(fp, "%% original \n");
    else      fprintf(fp, "%% overlap \n");
  }
  fprintf(fp, "];\n");

  fprintf(fp, "EToE=[\n");
  for(int e=0;e<mesh1->Nelements;++e){
    for(int f=0;f<mesh1->Nfaces;++f)
      fprintf(fp, "%d ", mesh1->EToE[e*mesh1->Nfaces+f]);
    if(e<mesh->Nelements) fprintf(fp, "%% original \n");
    else      fprintf(fp, "%% overlap \n");
  }
  fprintf(fp, "];\n");

  fprintf(fp, "EToB=[\n");
  for(int e=0;e<mesh1->Nelements;++e){
    for(int f=0;f<mesh1->Nfaces;++f)
      fprintf(fp, "%d ", mesh1->EToB[e*mesh1->Nfaces+f]);
    if(e<mesh->Nelements) fprintf(fp, "%% original \n");
    else      fprintf(fp, "%% overlap \n");
  }
  fprintf(fp, "];\n");
  
  fclose(fp);

  string outName;
  elliptic1->options.getArgs("OUTPUT FILE NAME", outName);
  sprintf(fname, "%s_oneRing_%04d",(char*)outName.c_str(), mesh->rank);
  ellipticPlotVTUHex3D(mesh1, fname, 0);
  
}
