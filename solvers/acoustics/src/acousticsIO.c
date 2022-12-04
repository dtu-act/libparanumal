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

void acousticsWritePressureField(acoustics_t *acoustics, WriteWaveFieldType waveFieldWriteType, std::vector<dfloat> timeSteps, int iter) {    
  if (waveFieldWriteType == Vtu) {    
    acousticsWriteVTU(acoustics);
  }
  else if  (waveFieldWriteType == Xdmf) {
    acousticsWriteXdmf(acoustics, timeSteps, iter);
  }
  else if  (waveFieldWriteType == Txt) {
    dfloat time = timeSteps[iter];
    acousticsWritePressureFieldTxt(acoustics, time);
  }
  else {
    throw std::exception();
  }
}

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

// write acoustic pressure as simple text file
// NOTE: not working properly, overwrites new file for each time step
void acousticsWritePressureFieldTxt(acoustics_t *acoustics, dfloat time)
{
  mesh3D *mesh = acoustics->mesh;

  std::string xFN = acoustics->outDir + "/x.txt";
  std::string yFN = acoustics->outDir + "/y.txt";
  std::string zFN = acoustics->outDir + "/z.txt";
  std::string pFN = acoustics->outDir + "/p.txt";

  // Print x,y,z, and pressure to txt files
  int nprocs, procid;
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  hlong globalNelements;
  MPI_Reduce(&mesh->Nelements, &globalNelements, 1, MPI_HLONG, MPI_SUM, 0, mesh->comm);
  string PREFIX = acoustics->simulationID;

  FILE *xFP, *yFP, *zFP, *pFP;

  for (int i = 0; i < nprocs; i++)
  {
    if (procid == i)
    {
      xFP = fopen(xFN.c_str(), "a");
      yFP = fopen(yFN.c_str(), "a");
      zFP = fopen(zFN.c_str(), "a");
      pFP = fopen(pFN.c_str(), "a");

      if (xFP == NULL || yFP == NULL || zFP == NULL || pFP == NULL) {
        throw std::invalid_argument("Could not open files for writing wave field [TXT].");
      }

      if (time < 1e-10) {
        // Write Np and Nelements as first line
        fprintf(xFP, "%d %d\n", mesh->Np, globalNelements);
        fprintf(yFP, "%d %d\n", mesh->Np, globalNelements);
        fprintf(zFP, "%d %d\n", mesh->Np, globalNelements);
        fprintf(pFP, "%d %d\n", mesh->Np, globalNelements);
      }

      for (dlong e = 0; e < mesh->Nelements; ++e)
      {
        for (int n = 0; n < mesh->Np; ++n)
        {
          dfloat xEA = mesh->x[n + mesh->Np * e];
          dfloat yEA = mesh->y[n + mesh->Np * e];
          dfloat zEA = mesh->z[n + mesh->Np * e];

          dlong qbaseEA = e * mesh->Np * mesh->Nfields + n;

          dfloat rEA = acoustics->q[qbaseEA + 0 * mesh->Np];

          fprintf(xFP, "%.15lf ", xEA);
          fprintf(yFP, "%.15lf ", yEA);
          fprintf(zFP, "%.15lf ", zEA);
          fprintf(pFP, "%.15lf ", rEA);
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

  // caluclate errors on the host
  // acousticsError(acoustics, time);
}

// void acousticsSnapshotXYZ(acoustics_t *acoustics, setupAide &newOptions)
// {
//   mesh3D *mesh = acoustics->mesh;
//   int nprocs, procid;
//   MPI_Comm_rank(MPI_COMM_WORLD, &procid);
//   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

//   hlong globalNelements;
//   MPI_Reduce(&mesh->Nelements, &globalNelements, 1, MPI_HLONG, MPI_SUM, 0, mesh->comm);

//   string PREFIX;
//   newOptions.getArgs("SNAPSHOTPREFIX", PREFIX);

//   FILE *xFP, *yFP, *zFP;
//   char xFN[BUFSIZ];
//   char yFN[BUFSIZ];
//   char zFN[BUFSIZ];

//   sprintf(xFN, "data/snapshot/%s_Snapshot_x.txt", (char *)PREFIX.c_str());
//   sprintf(yFN, "data/snapshot/%s_Snapshot_y.txt", (char *)PREFIX.c_str());
//   sprintf(zFN, "data/snapshot/%s_Snapshot_z.txt", (char *)PREFIX.c_str());

//   for (int i = 0; i < mesh->size; i++)
//   {
//     if (mesh->rank == i)
//     {
//       if (mesh->rank == 0)
//       {
//         xFP = fopen(xFN, "w");
//         yFP = fopen(yFN, "w");
//         zFP = fopen(zFN, "w");
//       }
//       else
//       {
//         xFP = fopen(xFN, "a");
//         yFP = fopen(yFN, "a");
//         zFP = fopen(zFN, "a");
//       }
//       for (dlong e = 0; e < mesh->Nelements; ++e)
//       {
//         for (int n = 0; n < mesh->Np; ++n)
//         {
//           dfloat xEA = mesh->x[n + mesh->Np * e];
//           dfloat yEA = mesh->y[n + mesh->Np * e];
//           dfloat zEA = mesh->z[n + mesh->Np * e];

//           fprintf(xFP, "%.15lf ", xEA);
//           fprintf(yFP, "%.15lf ", yEA);
//           fprintf(zFP, "%.15lf ", zEA);
//         }
//       }

//       fclose(xFP);
//       fclose(yFP);
//       fclose(zFP);
//     }
//     MPI_Barrier(MPI_COMM_WORLD);
//   }
// }

// void acousticsSnapshot(acoustics_t *acoustics, dfloat time, setupAide &newOptions, dlong openNewFile, dlong snapshotCounter)
// {

//   mesh3D *mesh = acoustics->mesh;
//   FILE *snapshotFilex;
//   FILE *snapshotFiley;
//   FILE *snapshotFilez;
//   FILE *snapshotFilep;

//   string PREFIX;
//   newOptions.getArgs("SNAPSHOTPREFIX", PREFIX);
//   char fnamex[BUFSIZ];
//   char fnamey[BUFSIZ];
//   char fnamez[BUFSIZ];
//   char fnamep[BUFSIZ];
//   char fnamet[BUFSIZ];
//   sprintf(fnamex, "data/snapshot/%s_Snapshot_vx.txt", (char *)PREFIX.c_str());
//   sprintf(fnamey, "data/snapshot/%s_Snapshot_vy.txt", (char *)PREFIX.c_str());
//   sprintf(fnamez, "data/snapshot/%s_Snapshot_vz.txt", (char *)PREFIX.c_str());
//   sprintf(fnamep, "data/snapshot/%s_Snapshot_p.txt", (char *)PREFIX.c_str());
//   sprintf(fnamet, "data/snapshot/%s_Snapshot_t.txt", (char *)PREFIX.c_str());

//   if (!mesh->rank)
//   {
//     FILE *snapshotFilet;
//     if (openNewFile == 1)
//     {
//       snapshotFilet = fopen(fnamet, "w");
//     }
//     else
//     {
//       snapshotFilet = fopen(fnamet, "a");
//     }

//     for (int itt = 0; itt < snapshotCounter; itt++)
//     {
//       fprintf(snapshotFilet, "%.15lf ", acoustics->Snapshott[itt]);
//     }
//     fclose(snapshotFilet);
//   }

//   for (int jtt = 0; jtt < snapshotCounter; jtt++)
//   {
//     for (int i = 0; i < mesh->size; i++)
//     {
//       if (mesh->rank == i)
//       {
//         if (mesh->rank == 0 && openNewFile == 1)
//         {
//           snapshotFilex = fopen(fnamex, "w");
//           snapshotFiley = fopen(fnamey, "w");
//           snapshotFilez = fopen(fnamez, "w");
//           snapshotFilep = fopen(fnamep, "w");
//           openNewFile = 0;
//         }
//         else
//         {
//           snapshotFilex = fopen(fnamex, "a");
//           snapshotFiley = fopen(fnamey, "a");
//           snapshotFilez = fopen(fnamez, "a");
//           snapshotFilep = fopen(fnamep, "a");
//         }
//         for (int itt = 0; itt < mesh->Nfields; itt++)
//         {
//           for (dlong e = 0; e < mesh->Nelements; ++e)
//           {
//             for (int n = 0; n < mesh->Np; ++n)
//             {

//               dlong qbaseEA = e * mesh->Np * mesh->Nfields + n;
//               dlong offsetSnapshot = mesh->Np * mesh->Nelements * mesh->Nfields * jtt;
//               dfloat out = acoustics->qSnapshot[offsetSnapshot + qbaseEA + itt * mesh->Np];

//               if (itt == 0)
//               {
//                 fprintf(snapshotFilep, "%.15lf ", out);
//               }
//               else if (itt == 1)
//               {
//                 fprintf(snapshotFilex, "%.15lf ", out);
//               }
//               else if (itt == 2)
//               {
//                 fprintf(snapshotFiley, "%.15lf ", out);
//               }
//               else if (itt == 3)
//               {
//                 fprintf(snapshotFilez, "%.15lf ", out);
//               }
//             }
//           }
//         }
//         if (mesh->rank + 1 == mesh->size)
//         {
//           fprintf(snapshotFilex, "\n");
//           fprintf(snapshotFiley, "\n");
//           fprintf(snapshotFilez, "\n");
//           fprintf(snapshotFilep, "\n");
//         }
//         fclose(snapshotFilex);
//         fclose(snapshotFiley);
//         fclose(snapshotFilez);
//         fclose(snapshotFilep);
//       }
//       MPI_Barrier(MPI_COMM_WORLD);
//     }
//   }
// }