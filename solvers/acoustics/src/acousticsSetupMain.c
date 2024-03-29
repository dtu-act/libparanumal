#include "acoustics.h"
#include <limits.h>

int acousticsSetupMain(setupAide &newOptions)
{  
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
  switch (elementType)
  {
  case TRIANGLES:
    mesh = meshSetupTri2D((char *)fileName.c_str(), N);
    break;
  case QUADRILATERALS:
    mesh = meshSetupQuad2D((char *)fileName.c_str(), N);
    break;
  case TETRAHEDRA:
    if (!curv)
    {
      mesh = meshSetupTet3D((char *)fileName.c_str(), N);
      break;
    }
    else
    {
      mesh = meshSetupTet3DCurv((char *)fileName.c_str(), N);
      break;
    }
  case HEXAHEDRA:
    mesh = meshSetupHex3D((char *)fileName.c_str(), N);
    break;
  }

  char *boundaryHeaderFileName; // could sprintf
  if (dim == 2)
    boundaryHeaderFileName = strdup(DACOUSTICS "/acousticsUniform2D.h"); // default
  if (dim == 3)
    boundaryHeaderFileName = strdup(DACOUSTICS "/acousticsUniform3D.h"); // default

  // set up acoustics stuff
  acoustics_t *acoustics = acousticsSetup(mesh, newOptions, boundaryHeaderFileName);
  acousticsFindReceiverElement(acoustics);
  // If receiver is on this core, allocate array for storage
  if (acoustics->NReceiversLocal > 0)
  {
    acoustics->o_recvElements =
        mesh->device.malloc(acoustics->NReceivers * sizeof(dlong), acoustics->recvElements);
    acoustics->o_recvElementsIdx =
        mesh->device.malloc(acoustics->NReceivers * sizeof(dlong), acoustics->recvElementsIdx);
    acoustics->qRecv = (dfloat *)calloc(acoustics->NReceiversLocal * acoustics->timeStepsOut.size(), sizeof(dfloat));
    acoustics->o_qRecv =
        mesh->device.malloc(acoustics->NReceiversLocal * RECV_COPY_RATE * sizeof(dfloat));
    acousticsRecvIntpolOperators(acoustics);
  }

  double startTime, endTime;
  startTime = MPI_Wtime();

  printf("Acoustics run...\n");
  acousticsRun(acoustics, newOptions);

  endTime = MPI_Wtime();
  if (!mesh->rank) {
    printf("Execution time: %lf\n", endTime - startTime);
  }

  std::string filepathJson = acoustics->outDir + "/simulation_parameters.json";
  acousticsWriteSimulationSettings(acoustics, filepathJson);

  return EXIT_SUCCESS;
}