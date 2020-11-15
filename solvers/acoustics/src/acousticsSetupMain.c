#include "acoustics.h"

int acousticsSetupMain(setupAide &newOptions) {
  
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
  switch (elementType) {
    case TRIANGLES:
      mesh = meshSetupTri2D((char*)fileName.c_str(), N); break;
    case QUADRILATERALS:
      mesh = meshSetupQuad2D((char*)fileName.c_str(), N); break;
    case TETRAHEDRA:
      if (!curv) {
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

  return 0;
}