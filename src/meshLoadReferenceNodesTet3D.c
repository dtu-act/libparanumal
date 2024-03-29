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

void meshLoadReferenceNodesTet3D(mesh3D *mesh, int N){

  char fname[BUFSIZ];
  sprintf(fname, DHOLMES "/nodes/tetN%02d.dat", N);

  FILE *fp = fopen(fname, "r");

  if (!fp) {
    printf("ERROR: Cannot open file: '%s'\n", fname);
    exit(-1);
  }

  mesh->N = N;
  mesh->Np = ((N+1)*(N+2)*(N+3))/6;
  mesh->Nfp = ((N+1)*(N+2))/2;

  int Nrows, Ncols;

  /* Nodal Data */
  readDfloatArray(fp, "Nodal r-coordinates", &(mesh->r),&Nrows,&Ncols);
  readDfloatArray(fp, "Nodal s-coordinates", &(mesh->s),&Nrows,&Ncols);
  readDfloatArray(fp, "Nodal t-coordinates", &(mesh->t),&Nrows,&Ncols);
  readDfloatArray(fp, "Nodal Dr differentiation matrix", &(mesh->Dr), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal Ds differentiation matrix", &(mesh->Ds), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal Dt differentiation matrix", &(mesh->Dt), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal Mass Matrix", &(mesh->MM), &Nrows, &Ncols);
  readIntArray   (fp, "Nodal Face nodes", &(mesh->faceNodes), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal Lift Matrix", &(mesh->LIFT), &Nrows, &Ncols);
  //readIntArray   (fp, "Nodal rotation permutations", &(mesh->rmapP), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal degree raise matrix", &(mesh->interpRaise), &Nrows, &Ncols);
  readDfloatArray(fp, "Nodal degree lower matrix", &(mesh->interpLower), &Nrows, &Ncols);

  /* Plotting data */ 
  readDfloatArray(fp, "Plotting r-coordinates", &(mesh->plotR),&Nrows,&Ncols);
  readDfloatArray(fp, "Plotting s-coordinates", &(mesh->plotS),&Nrows,&Ncols);
  readDfloatArray(fp, "Plotting t-coordinates", &(mesh->plotT),&Nrows,&Ncols);
  mesh->plotNp = Nrows;

  readDfloatArray(fp, "Plotting Interpolation Matrix", &(mesh->plotInterp),&Nrows,&Ncols);
  readIntArray   (fp, "Plotting triangulation", &(mesh->plotEToV), &Nrows, &Ncols);
  mesh->plotNelements = Nrows;
  mesh->plotNverts = Ncols;

  readIntArray(fp,"Contour plot EToV", &(mesh->contourEToV), &Nrows, &Ncols);
  readDfloatArray(fp,"Contour plot VX", &(mesh->contourVX), &Nrows, &Ncols);
  readDfloatArray(fp,"Contour plot VY", &(mesh->contourVY), &Nrows, &Ncols);
  readDfloatArray(fp,"Contour plot VZ", &(mesh->contourVZ), &Nrows, &Ncols);

  readDfloatArray(fp, "Contour plot Interpolation",&(mesh->contourInterp), &Nrows, &Ncols);
  readDfloatArray(fp, "Contour plot Linear Interpolation",&(mesh->contourInterp1), &Nrows, &Ncols);
  readDfloatArray(fp, "Contour plot Filter",&(mesh->contourFilter), &Nrows, &Ncols);

  /* Cubature data */ 
  if (N<7) {
    readDfloatArray(fp, "Cubature r-coordinates", &(mesh->cubr),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature s-coordinates", &(mesh->cubs),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature t-coordinates", &(mesh->cubt),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature weights", &(mesh->cubw),&Nrows,&Ncols);
    mesh->cubNp = Nrows;

    readDfloatArray(fp, "Cubature Interpolation Matrix", &(mesh->cubInterp),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature Weak Dr Differentiation Matrix", &(mesh->cubDrW),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature Weak Ds Differentiation Matrix", &(mesh->cubDsW),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature Weak Dt Differentiation Matrix", &(mesh->cubDtW),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature Projection Matrix", &(mesh->cubProject),&Nrows,&Ncols);
    readDfloatArray(fp, "Cubature Surface Interpolation Matrix", &(mesh->intInterp),&Nrows,&Ncols);
    mesh->intNfp = Nrows/mesh->Nfaces; //number of interpolation points per face

    readDfloatArray(fp, "Cubature Surface Lift Matrix", &(mesh->intLIFT),&Nrows,&Ncols);
  }


  /* Bernstein-Bezier data */ 
  readDfloatArray(fp, "Bernstein-Bezier Vandermonde Matrix", &(mesh->VB),&Nrows,&Ncols);
  readDfloatArray(fp, "Bernstein-Bezier Inverse Vandermonde Matrix", &(mesh->invVB),&Nrows,&Ncols);
  readIntArray   (fp, "Bernstein-Bezier sparse D0 differentiation ids", &(mesh->D0ids), &Nrows, &Ncols);  //Ncols should be 4
  readIntArray   (fp, "Bernstein-Bezier sparse D1 differentiation ids", &(mesh->D1ids), &Nrows, &Ncols);  //Ncols should be 4
  readIntArray   (fp, "Bernstein-Bezier sparse D2 differentiation ids", &(mesh->D2ids), &Nrows, &Ncols);  //Ncols should be 4
  readIntArray   (fp, "Bernstein-Bezier sparse D3 differentiation ids", &(mesh->D3ids), &Nrows, &Ncols);  //Ncols should be 4
  readDfloatArray(fp, "Bernstein-Bezier sparse D differentiation values", &(mesh->Dvals), &Nrows, &Ncols);//Ncols should be 4

  readIntArray   (fp, "Bernstein-Bezier sparse D0T transpose differentiation ids", &(mesh->D0Tids), &Nrows, &Ncols);  //Ncols should be 4
  readIntArray   (fp, "Bernstein-Bezier sparse D1T transpose differentiation ids", &(mesh->D1Tids), &Nrows, &Ncols);  //Ncols should be 4
  readIntArray   (fp, "Bernstein-Bezier sparse D2T transpose differentiation ids", &(mesh->D2Tids), &Nrows, &Ncols);  //Ncols should be 4
  readIntArray   (fp, "Bernstein-Bezier sparse D3T transpose differentiation ids", &(mesh->D3Tids), &Nrows, &Ncols);  //Ncols should be 4
  readDfloatArray(fp, "Bernstein-Bezier sparse DT transpose differentiation values", &(mesh->DTvals), &Nrows, &Ncols);//Ncols should be 4

  readIntArray   (fp, "Bernstein-Bezier L0 Matrix ids", &(mesh->L0ids), &Nrows, &Ncols);  
  readDfloatArray(fp, "Bernstein-Bezier L0 Matrix values", &(mesh->L0vals), &Nrows, &Ncols); //Ncols should be 7
  readIntArray   (fp, "Bernstein-Bezier EL lift ids", &(mesh->ELids), &Nrows, &Ncols);  
  readDfloatArray(fp, "Bernstein-Bezier EL lift values", &(mesh->ELvals), &Nrows, &Ncols); 
  mesh->max_EL_nnz = Ncols;

  readIntArray   (fp, "Bernstein-Bezier sparse 2D degree raise ids", &(mesh->BBRaiseids), &Nrows, &Ncols);     //Ncols should be 3
  readDfloatArray(fp, "Bernstein-Bezier sparse 2D degree raise values", &(mesh->BBRaiseVals), &Nrows, &Ncols); //Ncols should be 3 
  readDfloatArray(fp, "Bernstein-Bezier sparse 2D degree lower matrix", &(mesh->BBLower), &Nrows, &Ncols); 

  /* IPDG patch data */ 
  readDfloatArray(fp, "IPDG overlapping patch forward matrix", &(mesh->oasForwardDg), &Nrows, &Ncols);   
  readDfloatArray(fp, "IPDG overlapping patch diagonal scaling", &(mesh->oasDiagOpDg), &Nrows, &Ncols);   
  readDfloatArray(fp, "IPDG overlapping patch backward matrix", &(mesh->oasBackDg), &Nrows, &Ncols);   
  mesh->NpP = Nrows; //overlapping patch size


  /* SEMFEM data */ 
  readDfloatArray(fp, "SEMFEM r-coordinates", &(mesh->rFEM),&Nrows,&Ncols);
  readDfloatArray(fp, "SEMFEM s-coordinates", &(mesh->sFEM),&Nrows,&Ncols);
  readDfloatArray(fp, "SEMFEM t-coordinates", &(mesh->tFEM),&Nrows,&Ncols);
  mesh->NpFEM = Nrows;

  readIntArray   (fp, "SEMFEM reference mesh", &(mesh->FEMEToV), &Nrows, &Ncols);
  mesh->NelFEM = Nrows;

  readDfloatArray(fp, "SEMFEM interpolation matrix", &(mesh->SEMFEMInterp),&Nrows,&Ncols);


  fclose(fp);

  // find node indices of vertex nodes
  dfloat NODETOL = 1e-12;
  mesh->vertexNodes = (int*) calloc(mesh->Nverts, sizeof(int));
  for(int n=0;n<mesh->Np;++n){
    if( (mesh->r[n]+1)*(mesh->r[n]+1)+(mesh->s[n]+1)*(mesh->s[n]+1)+(mesh->t[n]+1)*(mesh->t[n]+1)<NODETOL)
      mesh->vertexNodes[0] = n;
    if( (mesh->r[n]-1)*(mesh->r[n]-1)+(mesh->s[n]+1)*(mesh->s[n]+1)+(mesh->t[n]+1)*(mesh->t[n]+1)<NODETOL)
      mesh->vertexNodes[1] = n;
    if( (mesh->r[n]+1)*(mesh->r[n]+1)+(mesh->s[n]-1)*(mesh->s[n]-1)+(mesh->t[n]+1)*(mesh->t[n]+1)<NODETOL)
      mesh->vertexNodes[2] = n;
    if( (mesh->r[n]+1)*(mesh->r[n]+1)+(mesh->s[n]+1)*(mesh->s[n]+1)+(mesh->t[n]-1)*(mesh->t[n]-1)<NODETOL)
      mesh->vertexNodes[3] = n;
  }
}


