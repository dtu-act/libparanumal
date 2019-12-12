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
#if 0
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

  //acousticsPlotVTU(acoustics, fname);
}




void acousticsSnapshotXYZ(acoustics_t *acoustics, setupAide &newOptions){
  mesh3D *mesh = acoustics->mesh;
  int nprocs, procid;
  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  hlong globalNelements;
  MPI_Reduce(&mesh->Nelements, &globalNelements, 1, MPI_HLONG, MPI_SUM, 0, mesh->comm);
  

  string PREFIX;
  newOptions.getArgs("SNAPSHOTPREFIX", PREFIX);

  FILE * xFP, * yFP, * zFP;
  char xFN[BUFSIZ];
  char yFN[BUFSIZ];
  char zFN[BUFSIZ];

  sprintf(xFN, "data/snapshot/%s_Snapshot_x.txt", (char*)PREFIX.c_str());
  sprintf(yFN, "data/snapshot/%s_Snapshot_y.txt", (char*)PREFIX.c_str());
  sprintf(zFN, "data/snapshot/%s_Snapshot_z.txt", (char*)PREFIX.c_str());
    
  for(int i = 0; i < mesh->size; i++){
    if(mesh->rank == i){
      if(mesh->rank == 0){
          xFP = fopen(xFN, "w");
          yFP = fopen(yFN, "w");
          zFP = fopen(zFN, "w");
      } else {
          xFP = fopen(xFN, "a");
          yFP = fopen(yFN, "a");
          zFP = fopen(zFN, "a");
          
      }
      for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat xEA = mesh->x[n + mesh->Np*e];
      dfloat yEA = mesh->y[n + mesh->Np*e];
      dfloat zEA = mesh->z[n + mesh->Np*e];


      fprintf(xFP, "%.15lf ",xEA);
      fprintf(yFP, "%.15lf ",yEA);
      fprintf(zFP, "%.15lf ",zEA);
      
        }
        //fprintf(xFP, "\n");
        //fprintf(yFP, "\n");
        //fprintf(zFP, "\n");
        
      }

      fclose(xFP);
      fclose(yFP);
      fclose(zFP);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

}



void acousticsSnapshot(acoustics_t *acoustics, dfloat time, setupAide &newOptions, dlong openNewFile, dlong snapshotCounter){
  
  mesh3D *mesh = acoustics->mesh;
  //acoustics->o_qSnapshot.copyTo(acoustics->qSnapshot);
  FILE *snapshotFilex;
  FILE *snapshotFiley;
  FILE *snapshotFilez;
  FILE *snapshotFilep;

  string PREFIX;
  newOptions.getArgs("SNAPSHOTPREFIX", PREFIX);
  char fnamex[BUFSIZ];
  char fnamey[BUFSIZ];
  char fnamez[BUFSIZ];
  char fnamep[BUFSIZ];
  char fnamet[BUFSIZ];
  sprintf(fnamex, "data/snapshot/%s_Snapshot_vx.txt", (char*)PREFIX.c_str());
  sprintf(fnamey, "data/snapshot/%s_Snapshot_vy.txt", (char*)PREFIX.c_str());
  sprintf(fnamez, "data/snapshot/%s_Snapshot_vz.txt", (char*)PREFIX.c_str());
  sprintf(fnamep, "data/snapshot/%s_Snapshot_p.txt", (char*)PREFIX.c_str());
  sprintf(fnamet, "data/snapshot/%s_Snapshot_t.txt", (char*)PREFIX.c_str());

  if(!mesh->rank){
    FILE * snapshotFilet;
    if(openNewFile == 1){
      snapshotFilet = fopen(fnamet,"w");
    } else {
       snapshotFilet = fopen(fnamet,"a");
    }

    for(int itt = 0; itt < snapshotCounter; itt++){
      fprintf(snapshotFilet, "%.15lf ", acoustics->Snapshott[itt]);
    }
    fclose(snapshotFilet);
  }

  for(int i = 0; i < mesh->size; i++){
    if(mesh->rank == i){
      if(mesh->rank == 0 && openNewFile == 1){
        snapshotFilex = fopen(fnamex,"w");
        snapshotFiley = fopen(fnamey,"w");
        snapshotFilez = fopen(fnamez,"w");
        snapshotFilep = fopen(fnamep,"w");
      } else {
        snapshotFilex = fopen(fnamex,"a");
        snapshotFiley = fopen(fnamey,"a");
        snapshotFilez = fopen(fnamez,"a");
        snapshotFilep = fopen(fnamep,"a");
      }
      for(int jtt = 0; jtt < snapshotCounter; jtt++){
        for(int itt = 0; itt < mesh->Nfields; itt++){
          for(dlong e=0;e<mesh->Nelements;++e){
            for(int n=0;n<mesh->Np;++n){


              dlong qbaseEA = e*mesh->Np*mesh->Nfields + n;
              dlong offsetSnapshot = mesh->Np*mesh->Nelements*mesh->Nfields*jtt;
              dfloat out = acoustics->qSnapshot[offsetSnapshot+qbaseEA+itt*mesh->Np];

              //if(offsetSnapshot+qbaseEA+itt*mesh->Np >= mesh->Np*mesh->Nelements*mesh->Nfields*10 - 2){
              //  printf("Hejsa her fra %d\n",offsetSnapshot+qbaseEA+itt*mesh->Np);
              //}

              if(itt == 0){
                fprintf(snapshotFilep, "%.15lf ",out);
              } else if(itt == 1){
                fprintf(snapshotFilex, "%.15lf ",out);
              } else if(itt == 2){
                fprintf(snapshotFiley, "%.15lf ",out);
              } else if(itt == 3){
                fprintf(snapshotFilez, "%.15lf ",out);
              }
            }
          }
        }
        if(mesh->rank+1 == mesh->size){
          fprintf(snapshotFilex, "\n");
          fprintf(snapshotFiley, "\n");
          fprintf(snapshotFilez, "\n");
          fprintf(snapshotFilep, "\n");
        }
      }

      fclose(snapshotFilex);
      fclose(snapshotFiley);
      fclose(snapshotFilez);
      fclose(snapshotFilep);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  #if 0
  FILE * temp;
  if(openNewFile){
    temp = fopen("data/qsnap.txt","w");
  }else{
    temp = fopen("data/qsnap.txt","a");
  }
  printf("h = %d\n",snapshotCounter);
  for(int i = 0; i < mesh->Np*mesh->Nelements*mesh->Nfields*snapshotCounter; i++){
    fprintf(temp,"%.15f\n",acoustics->qSnapshot[i]);
  }
  fclose(temp);
  #endif
}
