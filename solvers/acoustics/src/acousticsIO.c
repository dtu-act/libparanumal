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
#include <limits.h>
#include<filesystem>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace fs = std::filesystem;

void acousticsWriteSimulationSettings(acoustics_t *acoustics, string filename) {
  // https://github.com/nlohmann/json
  std::ofstream file(filename.c_str());
  json j;

  j["SimulationParameters"]["fmax"] = acoustics->fmax;
  j["SimulationParameters"]["c"] = acoustics->mesh->c;
  j["SimulationParameters"]["rho"] = acoustics->mesh->rho;
  j["SimulationParameters"]["sigma"] = acoustics->sigma0;
  j["SimulationParameters"]["dt"] = acoustics->mesh->dt;

  if (acoustics->sourceType == GaussianFunction) {
    j["SimulationParameters"]["SourcePosition"] = acoustics->sourcePosition;
  }

  file << j.dump(4) << std::endl;
}

// Print interpolated receiver to file
int acousticsWriteIRs(acoustics_t *acoustics, setupAide &newOptions) {      
  if (acoustics->NReceiversLocal > 0) {
    mesh_t *mesh = acoustics->mesh;  
    char fname[BUFSIZ];
    
    for (dlong iRecv = 0; iRecv < acoustics->NReceiversLocal; iRecv++) {
      sprintf(fname, "%s/%s_%02d.txt", (char*)acoustics->outDir.c_str(), (char*)acoustics->simulationID.c_str(), acoustics->recvElementsIdx[iRecv]);

      FILE *iFP = fopen(fname,"w");
      if (iFP == NULL) {
        printf("ERROR: receiver output file could not be opened %s)\n", fname);      
        return 1;
      }

      dfloat time = 0;

      for(int i = 0; i < mesh->NtimeSteps + 1; i++) { // +1 for including IC
        fprintf(iFP, "%.15lf %.15le\n", time, acoustics->qRecv[i+iRecv*mesh->NtimeSteps]);
        time += mesh->dt;
      }
      
      fclose(iFP);
    }
  }

  return 0;
}

int createDir(string path, bool deleteIfExists) {
  if (deleteIfExists) {
    fs::remove_all(path);
  }  

  std::error_code ec;
  fs::create_directories(path, ec);
  if (ec) {
    printf("Creating directory failed!\n"); 
    return -1;
  }

  return 0;
}