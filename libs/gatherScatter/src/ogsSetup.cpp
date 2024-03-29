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

#include "ogs.hpp"
#include "ogsKernels.hpp"
#include "ogsInterface.h"

typedef struct{

  dlong localId;    // local node id
  hlong baseId;     // original global index

  dlong newId;         // new global id
  int owned;

}parallelNode_t;

// compare on baseId then by localId
int compareBaseId(const void *a, const void *b){

  parallelNode_t *fa = (parallelNode_t*) a;
  parallelNode_t *fb = (parallelNode_t*) b;

  if(abs(fa->baseId) < abs(fb->baseId)) return -1; //group by abs(baseId)
  if(abs(fa->baseId) > abs(fb->baseId)) return +1;

  if(fa->localId < fb->localId) return -1; //sort by local id
  if(fa->localId > fb->localId) return +1;

  return 0;
}

// compare on haloOwned then localId
int compareLocalId(const void *a, const void *b){

  parallelNode_t *fa = (parallelNode_t*) a;
  parallelNode_t *fb = (parallelNode_t*) b;

  if(fa->localId < fb->localId) return -1;
  if(fa->localId > fb->localId) return +1;

  return 0;
}

ogs_t *ogsSetup(dlong N, hlong *ids, MPI_Comm &comm,
                int ogsUnique, int verbose, occa::device device){

  ogs_t *ogs = (ogs_t*) calloc(1, sizeof(ogs_t));

  //Keep track of how many gs handles we've created, and
  // build kernels if this is the first
  if (!ogs::Nrefs) ogs::initKernels(comm, device);
  ogs::Nrefs++;

  ogs->N = N;
  ogs->comm = comm;

  int rank, size;
  MPI_Comm_rank(ogs->comm, &rank);
  MPI_Comm_size(ogs->comm, &size);

  //use the host gs to find what nodes are local to this rank
  int *minRank = (int *) calloc(N,sizeof(int));
  int *maxRank = (int *) calloc(N,sizeof(int));
  hlong *flagIds   = (hlong *) calloc(N,sizeof(hlong));
  for (dlong i=0;i<N;i++) {
    minRank[i] = rank;
    maxRank[i] = rank;
    flagIds[i] = abs(ids[i]);
  }

  //make a symmetric host gs handle (calls gslib)
  ogs->hostGsh = ogsHostSetup(comm, N, flagIds, 0, 0);

  ogsHostGatherScatter(minRank, ogsInt, ogsMin, ogs->hostGsh); //minRank[n] contains the smallest rank taking part in the gather of node n
  ogsHostGatherScatter(maxRank, ogsInt, ogsMax, ogs->hostGsh); //maxRank[n] contains the largest rank taking part in the gather of node n

  if (ogsUnique) {
    //remake the host gs handle respecting the nonsymmetric behavior
    ogsHostFree(ogs->hostGsh);
    ogs->hostGsh = ogsHostSetup(comm, N, ids, 0, 0);
  }

  //count local and halo nodes
  ogs->Nlocal=0; ogs->Nhalo=0;
  for (dlong i=0;i<N;i++) {
    if (ids[i]==0) continue;

    if ((minRank[i]!=rank)||(maxRank[i]!=rank))
      ogs->Nhalo++;
    else
      ogs->Nlocal++;
  }

  //-----------Local GS setup -------------

  parallelNode_t *localNodes = (parallelNode_t*) calloc(ogs->Nlocal,
                                                  sizeof(parallelNode_t));

  ogs->NlocalGather = 0;
  dlong cnt=0;
  for (dlong i=0;i<N;i++) {
    if (ids[i]==0) continue;

    if ((minRank[i]==rank)&&(maxRank[i]==rank)) {
      localNodes[cnt].localId = i;
      localNodes[cnt].baseId  = ids[i];
      localNodes[cnt].owned   = 0;
      cnt++;
    }
  }

  // sort based on base ids then local id
  qsort(localNodes, ogs->Nlocal, sizeof(parallelNode_t), compareBaseId);

  if (ogs->Nlocal) {
    localNodes[0].newId = 0;
    localNodes[0].owned = 1;
  }
  for (dlong i=1;i<ogs->Nlocal;i++) {
    int s = 0;
    if (abs(localNodes[i].baseId)!=abs(localNodes[i-1].baseId)) {
      ogs->NlocalGather++;
      s = 1;
    }
    localNodes[i].newId = ogs->NlocalGather;
    localNodes[i].owned = s;
  }
  if (ogs->Nlocal) ogs->NlocalGather++;

  // sort based on local ids
  qsort(localNodes, ogs->Nlocal, sizeof(parallelNode_t), compareLocalId);

  //tally up how many nodes are being gathered to each gatherNode and
  //  map to a local ordering
  dlong *localGatherCounts = (dlong*) calloc(ogs->NlocalGather,sizeof(dlong));
  dlong *localGatherMap    = (dlong*) calloc(ogs->NlocalGather,sizeof(dlong));
  cnt = 0;
  for (dlong i=0;i<ogs->Nlocal;i++) {
    dlong newId = localNodes[i].newId; //get the ordered id

    if (localNodes[i].owned)
      localGatherMap[newId] = cnt++; //record a new index if this is a new gatherNode

    localNodes[i].newId = localGatherMap[newId]; //reorder
    localGatherCounts[localGatherMap[newId]]++;  //tally
  }
  free(localGatherMap);

  ogs->localGatherOffsets = (dlong*) calloc(ogs->NlocalGather+1,sizeof(dlong));
  for (dlong i=0;i<ogs->NlocalGather;i++) {
    ogs->localGatherOffsets[i+1] = ogs->localGatherOffsets[i] + localGatherCounts[i];
    localGatherCounts[i] = 0;
  }

  ogs->localGatherIds = (dlong*) calloc(ogs->Nlocal,sizeof(dlong));
  for (dlong i=0;i<ogs->Nlocal;i++) {
    dlong gatherId = localNodes[i].newId;
    dlong offset = ogs->localGatherOffsets[gatherId];
    int index  = localGatherCounts[gatherId];

    ogs->localGatherIds[offset+index] = localNodes[i].localId;
    localGatherCounts[gatherId]++;
  }
  free(localGatherCounts);

  ogs->o_localGatherOffsets = device.malloc((ogs->NlocalGather+1)*sizeof(dlong), ogs->localGatherOffsets);
  ogs->o_localGatherIds     = device.malloc((ogs->Nlocal)*sizeof(dlong), ogs->localGatherIds);

  free(localNodes);

  //-----------Halo GS setup -------------

  //set up the halo gatherScatter
  parallelNode_t *haloNodes = (parallelNode_t*) calloc(ogs->Nhalo,sizeof(parallelNode_t));

  cnt=0;
  for (dlong i=0;i<N;i++) {
    if (ids[i]==0) continue; //skip masked ids

    //find and record shared nodes
    if ((minRank[i]!=rank)||(maxRank[i]!=rank)) {
      haloNodes[cnt].localId = i;          //original id
      haloNodes[cnt].baseId  = ids[i]; //global id
      haloNodes[cnt].owned   = 0;
      cnt++;
    }
  }

  // sort based on base ids then local id
  qsort(haloNodes, ogs->Nhalo, sizeof(parallelNode_t), compareBaseId);

  ogs->NhaloGather=0;
  if (ogs->Nhalo) { //set the first node id
    haloNodes[0].newId = 0;
    haloNodes[0].owned = 1;
  }
  for (dlong i=1;i<ogs->Nhalo;i++) {
    int s = 0;
    if (abs(haloNodes[i].baseId)!=abs(haloNodes[i-1].baseId)) { //new gather node
      ogs->NhaloGather++;
      s = 1;
    }

    haloNodes[i].owned = s;
    haloNodes[i].newId = ogs->NhaloGather;
  }
  if (ogs->Nhalo) ogs->NhaloGather++;

  // sort based on local ids
  qsort(haloNodes, ogs->Nhalo, sizeof(parallelNode_t), compareLocalId);

  //make an array of just the gathered halo globalIds
  hlong *haloFlagIds = (hlong *) calloc(ogs->NhaloGather,sizeof(hlong));
  cnt = 0;
  for (dlong i=0;i<ogs->Nhalo;i++) {
    if (haloNodes[i].owned)
      haloFlagIds[cnt++] = haloNodes[i].baseId;
  }

  if (!ogsUnique) {
    //use gslib to uniquely flag a single id one node in each group
    //i.e. one unique node in each group is 'flagged' (kept positive),
    // while others are turned negative.
    ogsGsUnique(haloFlagIds, ogs->NhaloGather, comm);
  }

  //count how many node this rank will actually 'own'
  ogs->NownedHalo=0;
  for (dlong i=0;i<ogs->NhaloGather;i++) {
    if(haloFlagIds[i]>0) ogs->NownedHalo++;
  }

  //tally up how many nodes are being gathered to each gatherNode and
  //  map to a local ordering
  dlong *haloGatherCounts = (dlong*) calloc(ogs->NhaloGather,sizeof(dlong));
  dlong *haloGatherMap    = (dlong*) calloc(ogs->NhaloGather,sizeof(dlong));
  hlong *symIds      = (hlong *) calloc(ogs->NhaloGather,sizeof(hlong));
  hlong *nonSymIds   = (hlong *) calloc(ogs->NhaloGather,sizeof(hlong));

  cnt = 0;
  dlong cnt_owned = 0;
  dlong cnt_local = ogs->NownedHalo;
  for (dlong i=0;i<ogs->Nhalo;i++) {
    dlong newId = haloNodes[i].newId; //get the ordered id

    if (haloNodes[i].owned) { //gathered node
      dlong c;
      if (haloFlagIds[cnt]>0)
        c = cnt_owned++;
      else
        c = cnt_local++;

      symIds[c]    = abs(haloFlagIds[cnt]); //record the base id
      nonSymIds[c] =     haloFlagIds[cnt];  //record the base id
      haloGatherMap[newId] = c; //record a new index if this is a new gatherNode
      cnt++;
    }

    haloNodes[i].newId = haloGatherMap[newId];  //reorder
    haloGatherCounts[haloGatherMap[newId]]++;  //tally
  }
  free(haloGatherMap);

  ogs->haloGatherOffsets = (dlong*) calloc(ogs->NhaloGather+1,sizeof(dlong));
  for (dlong i=0;i<ogs->NhaloGather;i++) {
    ogs->haloGatherOffsets[i+1] = ogs->haloGatherOffsets[i] + haloGatherCounts[i];
    haloGatherCounts[i] = 0;
  }

  ogs->haloGatherIds = (dlong*) calloc(ogs->Nhalo,sizeof(dlong));
  for (dlong i=0;i<ogs->Nhalo;i++) {
    dlong gatherId = haloNodes[i].newId;
    dlong offset = ogs->haloGatherOffsets[gatherId];
    int index  = haloGatherCounts[gatherId];

    ogs->haloGatherIds[offset+index] = haloNodes[i].localId;
    haloGatherCounts[gatherId]++;
  }
  free(haloGatherCounts);

  ogs->o_haloGatherOffsets = device.malloc((ogs->NhaloGather+1)*sizeof(dlong), ogs->haloGatherOffsets);
  ogs->o_haloGatherIds     = device.malloc((ogs->Nhalo)*sizeof(dlong), ogs->haloGatherIds);

  //make a host gs handle
  ogs->haloGshSym    = ogsHostSetup(comm, ogs->NhaloGather, symIds,    0,0);
  ogs->haloGshNonSym = ogsHostSetup(comm, ogs->NhaloGather, nonSymIds, 0,0);

  free(symIds); free(nonSymIds);
  free(haloNodes);
  free(minRank); free(maxRank); free(flagIds);

  //total number of owned gathered nodes
  ogs->Ngather = ogs->NlocalGather+ogs->NownedHalo;

  ogs->device = device;

  // build degree vectors
  ogs->invDegree = (dfloat*) calloc(N, sizeof(dfloat));
  ogs->gatherInvDegree = (dfloat*) calloc(ogs->Ngather, sizeof(dfloat));
  for(dlong n=0;n<N;++n) ogs->invDegree[n] = 1;

  ogs->o_invDegree = device.malloc(N*sizeof(dfloat), ogs->invDegree);
  ogs->o_gatherInvDegree = device.malloc(ogs->Ngather*sizeof(dfloat), ogs->gatherInvDegree);

  ogsGather(ogs->o_gatherInvDegree, ogs->o_invDegree, ogsDfloat, ogsAdd, ogs);

  if(ogs->Ngather)
    ogs->o_gatherInvDegree.copyTo(ogs->gatherInvDegree);

  ogsScatter(ogs->o_invDegree, ogs->o_gatherInvDegree, ogsDfloat, ogsAdd, ogs);

  if (N) ogs->o_invDegree.copyTo(ogs->invDegree);

  for(dlong n=0;n<ogs->N;++n)
    ogs->invDegree[n] = 1./ogs->invDegree[n];

  for(dlong n=0;n<ogs->Ngather;++n)
    ogs->gatherInvDegree[n] = 1./ogs->gatherInvDegree[n];

  if(ogs->Ngather)
    ogs->o_gatherInvDegree.copyFrom(ogs->gatherInvDegree);

  if(ogs->N)
    ogs->o_invDegree.copyFrom(ogs->invDegree);

  return ogs;
}


void ogsFree(ogs_t *ogs) {

  if (ogs->Nlocal) {
    free(ogs->localGatherOffsets);
    free(ogs->localGatherIds);
    ogs->o_localGatherOffsets.free();
    ogs->o_localGatherIds.free();
  }

  if (ogs->Nhalo) {
    free(ogs->haloGatherOffsets);
    free(ogs->haloGatherIds);
    ogs->o_haloGatherOffsets.free();
    ogs->o_haloGatherIds.free();
    ogsHostFree(ogs->haloGshSym);
    ogsHostFree(ogs->haloGshNonSym);
  }

  if (ogs->N) {
    free(ogs->invDegree);
    ogs->o_invDegree.free();
    ogsHostFree(ogs->hostGsh);
  }

  if (ogs->Ngather) {
    free(ogs->gatherInvDegree);
    ogs->o_gatherInvDegree.free();
  }

  free(ogs);

  ogs::Nrefs--;
  if (!ogs::Nrefs) ogs::freeKernels();
}