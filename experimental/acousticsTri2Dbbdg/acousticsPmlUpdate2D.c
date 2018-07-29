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

#include "acoustics2D.h"


void acousticsMRABpmlUpdate2D(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, int lev, dfloat dt){

  for(int et=0;et<mesh->MRABpmlNelements[lev];et++){
    int e     = mesh->MRABpmlElementIds[lev][et];
    int pmlId = mesh->MRABpmlIds[lev][et];

    for(int n=0;n<mesh->Np;++n){
      int id = mesh->Nfields*(e*mesh->Np + n);
      int pid = mesh->pmlNfields*(pmlId*mesh->Np + n);

      int rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      int rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      int rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      int pmlrhsId1 = 3*pid + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->pmlNfields;
      int pmlrhsId2 = 3*pid + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->pmlNfields;
      int pmlrhsId3 = 3*pid + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->pmlNfields;
      
      mesh->pmlq[pid+0] += dt*(a1*mesh->pmlrhsq[pmlrhsId1+0] + a2*mesh->pmlrhsq[pmlrhsId2+0] + a3*mesh->pmlrhsq[pmlrhsId3+0]);
      mesh->pmlq[pid+1] += dt*(a1*mesh->pmlrhsq[pmlrhsId1+1] + a2*mesh->pmlrhsq[pmlrhsId2+1] + a3*mesh->pmlrhsq[pmlrhsId3+1]);
      mesh->q[id+0] += dt*(a1*mesh->rhsq[rhsId1+0] + a2*mesh->rhsq[rhsId2+0] + a3*mesh->rhsq[rhsId3+0]);
      mesh->q[id+1] += dt*(a1*mesh->rhsq[rhsId1+1] + a2*mesh->rhsq[rhsId2+1] + a3*mesh->rhsq[rhsId3+1]);
      mesh->q[id+2] = mesh->pmlq[pid+0]+mesh->pmlq[pid+1];
    }

    //write new traces to fQ
    for (int f =0;f<mesh->Nfaces;f++) {
      for (int n=0;n<mesh->Nfp;n++) {
        int id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
        int qidM = mesh->Nfields*mesh->vmapM[id];

        int qid = mesh->Nfields*id;
        // save trace node values of q
        for (int fld=0; fld<mesh->Nfields;fld++) {
          mesh->fQP[qid+fld] = mesh->q[qidM+fld];
          mesh->fQM[qid+fld] = mesh->q[qidM+fld];
        }
      }
    }
  }
}

void acousticsMRABpmlUpdateTrace2D(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, int lev, dfloat dt){

  for(int et=0;et<mesh->MRABpmlNhaloElements[lev];et++){
    int e = mesh->MRABpmlHaloElementIds[lev][et];
    int pmlId = mesh->MRABpmlHaloIds[lev][et];

    dfloat s_q[mesh->Np*mesh->Nfields];

    for(int n=0;n<mesh->Np;++n){
      int id = mesh->Nfields*(e*mesh->Np + n);
      int pid = mesh->pmlNfields*(pmlId*mesh->Np + n);

      int rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      int rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      int rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      int pmlrhsId1 = 3*pid + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->pmlNfields;
      int pmlrhsId2 = 3*pid + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->pmlNfields;
      int pmlrhsId3 = 3*pid + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->pmlNfields;
      
      // Increment solutions
      dfloat px = mesh->pmlq[pid+0] + dt*(a1*mesh->pmlrhsq[pmlrhsId1+0] + a2*mesh->pmlrhsq[pmlrhsId2+0] + a3*mesh->pmlrhsq[pmlrhsId3+0]);
      dfloat py = mesh->pmlq[pid+1] + dt*(a1*mesh->pmlrhsq[pmlrhsId1+1] + a2*mesh->pmlrhsq[pmlrhsId2+1] + a3*mesh->pmlrhsq[pmlrhsId3+1]);
      s_q[n*mesh->Nfields+0] = mesh->q[id+0] + dt*(a1*mesh->rhsq[rhsId1+0] + a2*mesh->rhsq[rhsId2+0] + a3*mesh->rhsq[rhsId3+0]);
      s_q[n*mesh->Nfields+1] = mesh->q[id+1] + dt*(a1*mesh->rhsq[rhsId1+1] + a2*mesh->rhsq[rhsId2+1] + a3*mesh->rhsq[rhsId3+1]);
      s_q[n*mesh->Nfields+2] = px+py;
    }

    //write new traces to fQ
    for (int f =0;f<mesh->Nfaces;f++) { 
      for (int n=0;n<mesh->Nfp;n++) {
        int id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
        int qidM = mesh->Nfields*(mesh->vmapM[id]-e*mesh->Np);

        int qid = n*mesh->Nfields + f*mesh->Nfp*mesh->Nfields + e*mesh->Nfaces*mesh->Nfp*mesh->Nfields;
        // save trace node values of q
        for (int fld=0; fld<mesh->Nfields;fld++) {
          mesh->fQM[qid+fld] = s_q[qidM+fld];
          mesh->fQP[qid+fld] = s_q[qidM+fld];
        }
      }
    }
  }
}


void acousticsMRABpmlUpdate2D_wadg(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, int lev, dfloat dt){

  for(int et=0;et<mesh->MRABpmlNelements[lev];et++){
    int e = mesh->MRABpmlElementIds[lev][et];
    int pmlId = mesh->MRABpmlIds[lev][et];

    dfloat p[mesh->Np];
    dfloat cubp[mesh->cubNp];

    for(int n=0;n<mesh->Np;++n){
      int id = mesh->Nfields*(e*mesh->Np + n);
      int pid = mesh->pmlNfields*(pmlId*mesh->Np + n);

      int rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      int rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      int rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      int pmlrhsId1 = 3*pid + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->pmlNfields;
      int pmlrhsId2 = 3*pid + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->pmlNfields;
      int pmlrhsId3 = 3*pid + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->pmlNfields;
      
      mesh->pmlq[pid+0] += dt*(a1*mesh->pmlrhsq[pmlrhsId1+0] + a2*mesh->pmlrhsq[pmlrhsId2+0] + a3*mesh->pmlrhsq[pmlrhsId3+0]);
      mesh->pmlq[pid+1] += dt*(a1*mesh->pmlrhsq[pmlrhsId1+1] + a2*mesh->pmlrhsq[pmlrhsId2+1] + a3*mesh->pmlrhsq[pmlrhsId3+1]);
      mesh->q[id+0] += dt*(a1*mesh->rhsq[rhsId1+0] + a2*mesh->rhsq[rhsId2+0] + a3*mesh->rhsq[rhsId3+0]);
      mesh->q[id+1] += dt*(a1*mesh->rhsq[rhsId1+1] + a2*mesh->rhsq[rhsId2+1] + a3*mesh->rhsq[rhsId3+1]);
      p[n] = mesh->pmlq[pid+0]+mesh->pmlq[pid+1];
    }

    // Interpolate rhs to cubature nodes
    for(int n=0;n<mesh->cubNp;++n){
      cubp[n] = 0.f;
      for (int i=0;i<mesh->Np;++i){
        cubp[n] += mesh->cubInterp[n*mesh->Np + i] * p[i];
      }
      // Multiply result by wavespeed c2 at cubature node
      cubp[n] *= mesh->c2[n + e*mesh->cubNp];
    }

    // Increment solution, project result back down
    for(int n=0;n<mesh->Np;++n){
      // Project scaled rhs down
      dfloat c2p = 0.f;
      for (int i=0;i<mesh->cubNp;++i){
        c2p += mesh->cubProject[n*mesh->cubNp + i] * cubp[i];
      }
      int id = mesh->Nfields*(e*mesh->Np + n);
      mesh->q[id+2] = c2p;
    }

    //write new traces to fQ
    for (int f =0;f<mesh->Nfaces;f++) {
      for (int n=0;n<mesh->Nfp;n++) {
        int id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
        int qidM = mesh->Nfields*mesh->vmapM[id];

        int qid = mesh->Nfields*id;
        // save trace node values of q
        for (int fld=0; fld<mesh->Nfields;fld++) {
          mesh->fQP[qid+fld] = mesh->q[qidM+fld];
          mesh->fQM[qid+fld] = mesh->q[qidM+fld];
        }
      }
    }
  }
}

void acousticsMRABpmlUpdateTrace2D_wadg(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, int lev, dfloat dt){

  for(int et=0;et<mesh->MRABpmlNhaloElements[lev];et++){
    int e = mesh->MRABpmlHaloElementIds[lev][et];
    int pmlId = mesh->MRABpmlHaloIds[lev][et];

    dfloat cubp[mesh->cubNp];
    dfloat s_q[mesh->Np*mesh->Nfields];

    for(int n=0;n<mesh->Np;++n){
      int id = mesh->Nfields*(e*mesh->Np + n);
      int pid = mesh->pmlNfields*(pmlId*mesh->Np + n);

      int rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->Nfields;
      int rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->Nfields;
      int rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->Nfields;
      int pmlrhsId1 = 3*pid + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->pmlNfields;
      int pmlrhsId2 = 3*pid + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->pmlNfields;
      int pmlrhsId3 = 3*pid + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->pmlNfields;
      
      // Increment solutions
      dfloat px = mesh->pmlq[pid+0] + dt*(a1*mesh->pmlrhsq[pmlrhsId1+0] + a2*mesh->pmlrhsq[pmlrhsId2+0] + a3*mesh->pmlrhsq[pmlrhsId3+0]);
      dfloat py = mesh->pmlq[pid+1] + dt*(a1*mesh->pmlrhsq[pmlrhsId1+1] + a2*mesh->pmlrhsq[pmlrhsId2+1] + a3*mesh->pmlrhsq[pmlrhsId3+1]);
      s_q[n*mesh->Nfields+0] = mesh->q[id+0] + dt*(a1*mesh->rhsq[rhsId1+0] + a2*mesh->rhsq[rhsId2+0] + a3*mesh->rhsq[rhsId3+0]);
      s_q[n*mesh->Nfields+1] = mesh->q[id+1] + dt*(a1*mesh->rhsq[rhsId1+1] + a2*mesh->rhsq[rhsId2+1] + a3*mesh->rhsq[rhsId3+1]);
      s_q[n*mesh->Nfields+2] = px+py;
    }

    // Interpolate rhs to cubature nodes
    for(int n=0;n<mesh->cubNp;++n){
      cubp[n] = 0.f;
      for (int i=0;i<mesh->Np;++i){
        cubp[n] += mesh->cubInterp[n*mesh->Np + i] * s_q[i*mesh->Nfields+2];
      }
      // Multiply result by wavespeed c2 at cubature node
      cubp[n] *= mesh->c2[n + e*mesh->cubNp];
    }

    // Increment solution, project result back down
    for(int n=0;n<mesh->Np;++n){
      // Project scaled rhs down
      s_q[n*mesh->Nfields+2] = 0.f;
      for (int i=0;i<mesh->cubNp;++i){
        s_q[n*mesh->Nfields+2] += mesh->cubProject[n*mesh->cubNp + i] * cubp[i];
      }
    }

    //write new traces to fQ
    for (int f =0;f<mesh->Nfaces;f++) {
      for (int n=0;n<mesh->Nfp;n++) {
        int id  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
        int qidM = mesh->Nfields*(mesh->vmapM[id]-e*mesh->Np);

        int qid = n*mesh->Nfields + f*mesh->Nfp*mesh->Nfields + e*mesh->Nfaces*mesh->Nfp*mesh->Nfields;
        // save trace node values of q
        for (int fld=0; fld<mesh->Nfields;fld++) {
          mesh->fQP[qid+fld] = s_q[qidM+fld];
          mesh->fQM[qid+fld] = s_q[qidM+fld];
        }
      }
    }
  }
}


void acousticsMRABpmlUpdate2D2(mesh2D *mesh,  
                           dfloat a1,
                           dfloat a2,
                           dfloat a3, int lev, dfloat dt){

  for(int m=0;m<mesh->MRABpmlNelements[lev];m++){
    int pmlId = mesh->MRABpmlIds[lev][m];
    for(int n=0;n<mesh->Np;++n){
      int id = mesh->pmlNfields*(pmlId*mesh->Np + n);

      int rhsId1 = 3*id + ((mesh->MRABshiftIndex[lev]+0)%3)*mesh->pmlNfields;
      int rhsId2 = 3*id + ((mesh->MRABshiftIndex[lev]+1)%3)*mesh->pmlNfields;
      int rhsId3 = 3*id + ((mesh->MRABshiftIndex[lev]+2)%3)*mesh->pmlNfields;
      for (int fld=0;fld<mesh->pmlNfields;fld++)
        mesh->pmlq[id+fld] += dt*(a1*mesh->pmlrhsq[rhsId1+fld] + a2*mesh->pmlrhsq[rhsId2+fld] + a3*mesh->pmlrhsq[rhsId3+fld]);
    }
  }
}



