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

void acousticsWritePressureField(acoustics_t *acoustics) {  
  char fname[BUFSIZ];
  mesh3D *mesh = acoustics->mesh;

  sprintf(fname, "%s/pressure_field_%04d_%04d.vtu", (char*)acoustics->outDir.c_str(), mesh->rank, acoustics->frame++);
  acousticsPlotVTU(acoustics, fname);
}

// Print interpolated receiver to file
int acousticsWriteIRs(acoustics_t *acoustics, setupAide &newOptions) {  
  mesh_t *mesh = acoustics->mesh;
  string PREFIX;  
  char fname[BUFSIZ];
  
  if (acoustics->NReceiversLocal == 0) {
    printf("ERROR: No receivers defined\n");
    return 0;
  }   
  if (newOptions.getArgs("SIMULATION_ID", PREFIX) == 0) {
    printf("ERROR: [SIMULATION_ID] tag missing");
    return 1;
  }  

  for (dlong iRecv = 0; iRecv < acoustics->NReceiversLocal; iRecv++){    
    sprintf(fname, "%s/%s_%02d.txt", (char*)acoustics->outDir.c_str(), (char*)PREFIX.c_str(), acoustics->recvElementsIdx[iRecv]);

    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) == NULL) {
      perror("getcwd() error");      
      return 1;
    }

    FILE *iFP = fopen(fname,"w");
    if (iFP == NULL) {
      printf("ERROR: receiver output file could not be opened %s/%s)\n", cwd, fname);      
      return 1;
    }

    dfloat time = 0;

    for(int i = 0; i < mesh->NtimeSteps + 1; i++) { // +1 for including IC
      fprintf(iFP, "%.15lf %.15le\n", time, acoustics->qRecv[i+iRecv*mesh->NtimeSteps]);
      time += mesh->dt;
    }

    printf("Receiver impulse response was written to disk: %s/%s\n", cwd, fname);
    fclose(iFP);
  }

  return 0;
}

int createDir(string path) {
  char *outDir = (char*)path.c_str();

  // create output folder if not existing
  struct stat st = {0};
  if (stat(outDir, &st) == -1) {
    return mkdir(outDir, 0700);
  }

  return 0;
}

// write acoustic pressure as simple text file instead of VTU (for Matlab?)
void acousticsWritePressureFieldTxt(acoustics_t *acoustics, dfloat time, setupAide &newOptions)
{
  mesh3D *mesh = acoustics->mesh;

  // Print x,y,z, and pressure to txt files
  int nprocs, procid;
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  hlong globalNelements;
  MPI_Reduce(&mesh->Nelements, &globalNelements, 1, MPI_HLONG, MPI_SUM, 0, mesh->comm);
  string PREFIX;
  newOptions.getArgs("RECEIVERPREFIX", PREFIX);

  FILE *xFP, *yFP, *zFP, *pFP;
  char xFN[BUFSIZ];
  char yFN[BUFSIZ];
  char zFN[BUFSIZ];
  char pFN[BUFSIZ];
  sprintf(xFN, "data/x_%s.txt", (char *)PREFIX.c_str());
  sprintf(yFN, "data/y_%s.txt", (char *)PREFIX.c_str());
  sprintf(zFN, "data/z_%s.txt", (char *)PREFIX.c_str());
  sprintf(pFN, "data/p_%s.txt", (char *)PREFIX.c_str());
  for (int i = 0; i < nprocs; i++)
  {
    if (procid == i)
    {
      if (procid == 0)
      {
        xFP = fopen(xFN, "w");
        yFP = fopen(yFN, "w");
        zFP = fopen(zFN, "w");
        pFP = fopen(pFN, "w");

        // Write Np and Nelements as first line in x output
        fprintf(xFP, "%d %d\n", mesh->Np, globalNelements);
      }
      else
      {
        xFP = fopen(xFN, "a");
        yFP = fopen(yFN, "a");
        zFP = fopen(zFN, "a");
        pFP = fopen(pFN, "a");
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

  // do error stuff on host
  acousticsError(acoustics, time);
  // output field files
  char fname[BUFSIZ];

  sprintf(fname, "foo_%04d_%04d.vtu", mesh->rank, acoustics->frame++);
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