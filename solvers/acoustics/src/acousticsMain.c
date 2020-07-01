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


int main(int argc, char **argv){

  if(argc!=2){
    printf("usage2: ./acousticsMain setupfile\n");
    exit(-1);
  }

  // start up MPI
  MPI_Init(&argc, &argv);

  // if argv > 2 then should load input data from argv
  setupAide newOptions(argv[1]);
  
  // set up mesh stuff
  string fileName;
  int N, dim, elementType, curv;
  
  newOptions.getArgs("MESH FILE", fileName);
  newOptions.getArgs("POLYNOMIAL DEGREE", N);
  newOptions.getArgs("ELEMENT TYPE", elementType);
  newOptions.getArgs("MESH DIMENSION", dim);
  newOptions.getArgs("CURVILINEAR MESH", curv);
  // set up mesh
  mesh_t *mesh;
  switch(elementType){
  case TRIANGLES:
    mesh = meshSetupTri2D((char*)fileName.c_str(), N); break;
  case QUADRILATERALS:
    mesh = meshSetupQuad2D((char*)fileName.c_str(), N); break;
  case TETRAHEDRA:
    if(!curv){
      mesh = meshSetupTet3D((char*)fileName.c_str(), N); break;
    } else {
      mesh = meshSetupTet3DCurv((char*)fileName.c_str(), N); break;
    }
  case HEXAHEDRA:
    mesh = meshSetupHex3D((char*)fileName.c_str(), N); break;
  }
  char *boundaryHeaderFileName; // could sprintf
  if(dim==2)
    boundaryHeaderFileName = strdup(DACOUSTICS "/acousticsUniform2D.h"); // default
  if(dim==3)
    boundaryHeaderFileName = strdup(DACOUSTICS "/acousticsUniform3D.h"); // default

  // set up acoustics stuff
  acoustics_t *acoustics = acousticsSetup(mesh, newOptions, boundaryHeaderFileName);
  acousticsFindReceiverElement(acoustics);
  // If receiver is on this core, allocate array for storage
  if(acoustics->NReceiversLocal > 0){
    acoustics->o_recvElements = 
      mesh->device.malloc(acoustics->NReceivers*sizeof(dlong), acoustics->recvElements);
    acoustics->o_recvElementsIdx = 
      mesh->device.malloc(acoustics->NReceivers*sizeof(dlong), acoustics->recvElementsIdx);
    acoustics->qRecv = (dfloat*) calloc(acoustics->NReceiversLocal*mesh->NtimeSteps, sizeof(dfloat));
    acoustics->o_qRecv =
    mesh->device.malloc(acoustics->NReceiversLocal*recvCopyRate*sizeof(dfloat));
    acousticsRecvIntpolOperators(acoustics);
  }

  // run
  double startTime, endTime;
  startTime = MPI_Wtime();
  acousticsRun(acoustics, newOptions);
  endTime = MPI_Wtime();
  if(!mesh->rank){printf("Execution time: %lf\n",endTime-startTime);}
  acousticsReport(acoustics, mesh->finalTime, newOptions);


  //---------RECEIVER---------
  acousticsPrintReceiversToFile(acoustics, newOptions);
  
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
