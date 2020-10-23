#include "acoustics.h"
#include <limits.h>

dlong factorial(dlong x){
  dlong res = 1;
  for(int i = 2; i <= x; i++){
    res *= i;
  }
  return res;
}

// OLD - From when interpolation was done on host.
void acousticsReceiverInterpolation(acoustics_t *acoustics){
	
  mesh_t *mesh = acoustics->mesh;
  for(int iRecv = 0; iRecv < acoustics->NReceiversLocal; iRecv++){

    // Receriver element
    dlong rIdx = acoustics->recvElementsIdx[iRecv];
    dlong recvElement = acoustics->recvElements[rIdx];
    dfloat xRecvElement[4];
    dfloat yRecvElement[4];
    dfloat zRecvElement[4];

    xRecvElement[0] = mesh->EX[recvElement*mesh->Nverts+0];
    xRecvElement[1] = mesh->EX[recvElement*mesh->Nverts+1];
    xRecvElement[2] = mesh->EX[recvElement*mesh->Nverts+2];
    xRecvElement[3] = mesh->EX[recvElement*mesh->Nverts+3];
    
    yRecvElement[0] = mesh->EY[recvElement*mesh->Nverts+0];
    yRecvElement[1] = mesh->EY[recvElement*mesh->Nverts+1];
    yRecvElement[2] = mesh->EY[recvElement*mesh->Nverts+2];
    yRecvElement[3] = mesh->EY[recvElement*mesh->Nverts+3];

    zRecvElement[0] = mesh->EZ[recvElement*mesh->Nverts+0];
    zRecvElement[1] = mesh->EZ[recvElement*mesh->Nverts+1];
    zRecvElement[2] = mesh->EZ[recvElement*mesh->Nverts+2];
    zRecvElement[3] = mesh->EZ[recvElement*mesh->Nverts+3];

    dfloat L1_rec = -(xRecvElement[1]*yRecvElement[2]*zRecvElement[3] - xRecvElement[1]*yRecvElement[2]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[1]*yRecvElement[3]*zRecvElement[2] + xRecvElement[1]*yRecvElement[3]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[1]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[2] - xRecvElement[1]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[3] - xRecvElement[2]*yRecvElement[1]*zRecvElement[3] + xRecvElement[2]*yRecvElement[1]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[2]*yRecvElement[3]*zRecvElement[1] - xRecvElement[2]*yRecvElement[3]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[2]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[1] + xRecvElement[2]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[3] + xRecvElement[3]*yRecvElement[1]*zRecvElement[2] - xRecvElement[3]*yRecvElement[1]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[3]*yRecvElement[2]*zRecvElement[1] + xRecvElement[3]*yRecvElement[2]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[3]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[1] - xRecvElement[3]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[2] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[1]*zRecvElement[2] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[1]*zRecvElement[3] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[2]*zRecvElement[1] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[2]*zRecvElement[3] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[3]*zRecvElement[1] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[3]*zRecvElement[2])/(xRecvElement[0]*yRecvElement[1]*zRecvElement[2] - xRecvElement[0]*yRecvElement[1]*zRecvElement[3] - xRecvElement[0]*yRecvElement[2]*zRecvElement[1] + xRecvElement[0]*yRecvElement[2]*zRecvElement[3] + xRecvElement[0]*yRecvElement[3]*zRecvElement[1] - xRecvElement[0]*yRecvElement[3]*zRecvElement[2] - xRecvElement[1]*yRecvElement[0]*zRecvElement[2] + xRecvElement[1]*yRecvElement[0]*zRecvElement[3] + xRecvElement[1]*yRecvElement[2]*zRecvElement[0] - xRecvElement[1]*yRecvElement[2]*zRecvElement[3] - xRecvElement[1]*yRecvElement[3]*zRecvElement[0] + xRecvElement[1]*yRecvElement[3]*zRecvElement[2] + xRecvElement[2]*yRecvElement[0]*zRecvElement[1] - xRecvElement[2]*yRecvElement[0]*zRecvElement[3] - xRecvElement[2]*yRecvElement[1]*zRecvElement[0] + xRecvElement[2]*yRecvElement[1]*zRecvElement[3] + xRecvElement[2]*yRecvElement[3]*zRecvElement[0] - xRecvElement[2]*yRecvElement[3]*zRecvElement[1] - xRecvElement[3]*yRecvElement[0]*zRecvElement[1] + xRecvElement[3]*yRecvElement[0]*zRecvElement[2] + xRecvElement[3]*yRecvElement[1]*zRecvElement[0] - xRecvElement[3]*yRecvElement[1]*zRecvElement[2] - xRecvElement[3]*yRecvElement[2]*zRecvElement[0] + xRecvElement[3]*yRecvElement[2]*zRecvElement[1]);
    dfloat L2_rec = (xRecvElement[0]*yRecvElement[2]*zRecvElement[3] - xRecvElement[0]*yRecvElement[2]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[0]*yRecvElement[3]*zRecvElement[2] + xRecvElement[0]*yRecvElement[3]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[0]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[2] - xRecvElement[0]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[3] - xRecvElement[2]*yRecvElement[0]*zRecvElement[3] + xRecvElement[2]*yRecvElement[0]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[2]*yRecvElement[3]*zRecvElement[0] - xRecvElement[2]*yRecvElement[3]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[2]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[0] + xRecvElement[2]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[3] + xRecvElement[3]*yRecvElement[0]*zRecvElement[2] - xRecvElement[3]*yRecvElement[0]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[3]*yRecvElement[2]*zRecvElement[0] + xRecvElement[3]*yRecvElement[2]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[3]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[0] - xRecvElement[3]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[2] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[0]*zRecvElement[2] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[0]*zRecvElement[3] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[2]*zRecvElement[0] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[2]*zRecvElement[3] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[3]*zRecvElement[0] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[3]*zRecvElement[2])/(xRecvElement[0]*yRecvElement[1]*zRecvElement[2] - xRecvElement[0]*yRecvElement[1]*zRecvElement[3] - xRecvElement[0]*yRecvElement[2]*zRecvElement[1] + xRecvElement[0]*yRecvElement[2]*zRecvElement[3] + xRecvElement[0]*yRecvElement[3]*zRecvElement[1] - xRecvElement[0]*yRecvElement[3]*zRecvElement[2] - xRecvElement[1]*yRecvElement[0]*zRecvElement[2] + xRecvElement[1]*yRecvElement[0]*zRecvElement[3] + xRecvElement[1]*yRecvElement[2]*zRecvElement[0] - xRecvElement[1]*yRecvElement[2]*zRecvElement[3] - xRecvElement[1]*yRecvElement[3]*zRecvElement[0] + xRecvElement[1]*yRecvElement[3]*zRecvElement[2] + xRecvElement[2]*yRecvElement[0]*zRecvElement[1] - xRecvElement[2]*yRecvElement[0]*zRecvElement[3] - xRecvElement[2]*yRecvElement[1]*zRecvElement[0] + xRecvElement[2]*yRecvElement[1]*zRecvElement[3] + xRecvElement[2]*yRecvElement[3]*zRecvElement[0] - xRecvElement[2]*yRecvElement[3]*zRecvElement[1] - xRecvElement[3]*yRecvElement[0]*zRecvElement[1] + xRecvElement[3]*yRecvElement[0]*zRecvElement[2] + xRecvElement[3]*yRecvElement[1]*zRecvElement[0] - xRecvElement[3]*yRecvElement[1]*zRecvElement[2] - xRecvElement[3]*yRecvElement[2]*zRecvElement[0] + xRecvElement[3]*yRecvElement[2]*zRecvElement[1]);
    dfloat L3_rec = -(xRecvElement[0]*yRecvElement[1]*zRecvElement[3] - xRecvElement[0]*yRecvElement[1]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[0]*yRecvElement[3]*zRecvElement[1] + xRecvElement[0]*yRecvElement[3]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[0]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[1] - xRecvElement[0]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[3] - xRecvElement[1]*yRecvElement[0]*zRecvElement[3] + xRecvElement[1]*yRecvElement[0]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[1]*yRecvElement[3]*zRecvElement[0] - xRecvElement[1]*yRecvElement[3]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[1]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[0] + xRecvElement[1]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[3] + xRecvElement[3]*yRecvElement[0]*zRecvElement[1] - xRecvElement[3]*yRecvElement[0]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[3]*yRecvElement[1]*zRecvElement[0] + xRecvElement[3]*yRecvElement[1]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[3]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[0] - xRecvElement[3]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[1] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[0]*zRecvElement[1] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[0]*zRecvElement[3] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[1]*zRecvElement[0] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[1]*zRecvElement[3] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[3]*zRecvElement[0] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[3]*zRecvElement[1])/(xRecvElement[0]*yRecvElement[1]*zRecvElement[2] - xRecvElement[0]*yRecvElement[1]*zRecvElement[3] - xRecvElement[0]*yRecvElement[2]*zRecvElement[1] + xRecvElement[0]*yRecvElement[2]*zRecvElement[3] + xRecvElement[0]*yRecvElement[3]*zRecvElement[1] - xRecvElement[0]*yRecvElement[3]*zRecvElement[2] - xRecvElement[1]*yRecvElement[0]*zRecvElement[2] + xRecvElement[1]*yRecvElement[0]*zRecvElement[3] + xRecvElement[1]*yRecvElement[2]*zRecvElement[0] - xRecvElement[1]*yRecvElement[2]*zRecvElement[3] - xRecvElement[1]*yRecvElement[3]*zRecvElement[0] + xRecvElement[1]*yRecvElement[3]*zRecvElement[2] + xRecvElement[2]*yRecvElement[0]*zRecvElement[1] - xRecvElement[2]*yRecvElement[0]*zRecvElement[3] - xRecvElement[2]*yRecvElement[1]*zRecvElement[0] + xRecvElement[2]*yRecvElement[1]*zRecvElement[3] + xRecvElement[2]*yRecvElement[3]*zRecvElement[0] - xRecvElement[2]*yRecvElement[3]*zRecvElement[1] - xRecvElement[3]*yRecvElement[0]*zRecvElement[1] + xRecvElement[3]*yRecvElement[0]*zRecvElement[2] + xRecvElement[3]*yRecvElement[1]*zRecvElement[0] - xRecvElement[3]*yRecvElement[1]*zRecvElement[2] - xRecvElement[3]*yRecvElement[2]*zRecvElement[0] + xRecvElement[3]*yRecvElement[2]*zRecvElement[1]);
    dfloat L4_rec = (xRecvElement[0]*yRecvElement[1]*zRecvElement[2] - xRecvElement[0]*yRecvElement[1]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[0]*yRecvElement[2]*zRecvElement[1] + xRecvElement[0]*yRecvElement[2]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[0]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[1] - xRecvElement[0]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[2] - xRecvElement[1]*yRecvElement[0]*zRecvElement[2] + xRecvElement[1]*yRecvElement[0]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[1]*yRecvElement[2]*zRecvElement[0] - xRecvElement[1]*yRecvElement[2]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[1]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[0] + xRecvElement[1]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[2] + xRecvElement[2]*yRecvElement[0]*zRecvElement[1] - xRecvElement[2]*yRecvElement[0]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[2]*yRecvElement[1]*zRecvElement[0] + xRecvElement[2]*yRecvElement[1]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[2]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[0] - xRecvElement[2]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[1] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[0]*zRecvElement[1] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[0]*zRecvElement[2] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[1]*zRecvElement[0] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[1]*zRecvElement[2] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[2]*zRecvElement[0] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[2]*zRecvElement[1])/(xRecvElement[0]*yRecvElement[1]*zRecvElement[2] - xRecvElement[0]*yRecvElement[1]*zRecvElement[3] - xRecvElement[0]*yRecvElement[2]*zRecvElement[1] + xRecvElement[0]*yRecvElement[2]*zRecvElement[3] + xRecvElement[0]*yRecvElement[3]*zRecvElement[1] - xRecvElement[0]*yRecvElement[3]*zRecvElement[2] - xRecvElement[1]*yRecvElement[0]*zRecvElement[2] + xRecvElement[1]*yRecvElement[0]*zRecvElement[3] + xRecvElement[1]*yRecvElement[2]*zRecvElement[0] - xRecvElement[1]*yRecvElement[2]*zRecvElement[3] - xRecvElement[1]*yRecvElement[3]*zRecvElement[0] + xRecvElement[1]*yRecvElement[3]*zRecvElement[2] + xRecvElement[2]*yRecvElement[0]*zRecvElement[1] - xRecvElement[2]*yRecvElement[0]*zRecvElement[3] - xRecvElement[2]*yRecvElement[1]*zRecvElement[0] + xRecvElement[2]*yRecvElement[1]*zRecvElement[3] + xRecvElement[2]*yRecvElement[3]*zRecvElement[0] - xRecvElement[2]*yRecvElement[3]*zRecvElement[1] - xRecvElement[3]*yRecvElement[0]*zRecvElement[1] + xRecvElement[3]*yRecvElement[0]*zRecvElement[2] + xRecvElement[3]*yRecvElement[1]*zRecvElement[0] - xRecvElement[3]*yRecvElement[1]*zRecvElement[2] - xRecvElement[3]*yRecvElement[2]*zRecvElement[0] + xRecvElement[3]*yRecvElement[2]*zRecvElement[1]);

    //Vandermonde Berstein matrix in receiver point
    dfloat *VB_rec;
    VB_rec = (dfloat*) calloc(mesh->Np, sizeof(dfloat));

    int sk = 0;
    for(int l = 0; l <= mesh->N; l++){
      for(int k = 0; k <= mesh->N - l; k++){
        for(int j = 0; j <= mesh->N - k - l; j++){
          int i = mesh->N - j - k - l;
          dfloat temp = factorial(mesh->N)/(factorial(i)*factorial(j)*factorial(k)*factorial(l));
          VB_rec[sk] = temp*pow(L1_rec,i)*pow(L2_rec,j)*pow(L3_rec,k)*pow(L4_rec,l);
          sk++;
        }
      }
    }
    // interpolation
    dfloat *intpol;
    intpol = (dfloat*) calloc(mesh->Np, sizeof(dfloat));

    for(int j = 0; j < mesh->Np; j++){
      for(int i = 0; i < mesh->Np; i++){
        intpol[i] += VB_rec[j]*mesh->invVB[i+j*mesh->Np];
      }
    }
    #if 0
    printf("Original, irecv = %d, ele = %d\n",iRecv,recvElement);
    for(int i = 0; i < mesh->Np; i++){
      printf("%.15lf \n",intpol[i]);
    }
    printf("\n");
    #endif


    #if 0
    // Interpolated receiver
    dfloat *intRecv;
    intRecv = (dfloat*) calloc(mesh->NtimeSteps, sizeof(dfloat)); // DOES NOT HAVE INITIAL CONDITION! Add to qRecv also!

    for(int i = 0; i < mesh->NtimeSteps; i++){
      dfloat interpolated = 0;
      dlong qRecvOffset = mesh->Np*i + iRecv*mesh->Np*mesh->NtimeSteps;
      for(int j = 0; j < mesh->Np;j++){
        interpolated += intpol[j]*acoustics->qRecv[qRecvOffset+j];
      }
      intRecv[i] = interpolated;
    }

    // Print interpolated receiver to file
    FILE *iFP;
    char fname[BUFSIZ];

    sprintf(fname, "data/interpolatedRecvPoint_%02d.txt", acoustics->recvElementsIdx[iRecv]);
    iFP = fopen(fname,"w");
    for(int i = 0; i < mesh->NtimeSteps; i++){
      fprintf(iFP, "%.15lf\n", intRecv[i]);
    }
    fclose(iFP);
    #endif
    free(VB_rec);
    free(intpol);
    //free(intRecv);
  }
}

// Create interpolation operators
void acousticsRecvIntpolOperators(acoustics_t *acoustics){
	
  mesh_t *mesh = acoustics->mesh;
  dfloat *intpol;
  intpol = (dfloat*) calloc(acoustics->NReceiversLocal*mesh->Np, sizeof(dfloat));
  
  for(int iRecv = 0; iRecv < acoustics->NReceiversLocal; iRecv++){

    // Receriver element
    dlong rIdx = acoustics->recvElementsIdx[iRecv];
    dlong recvElement = acoustics->recvElements[rIdx];
    dfloat xRecvElement[4];
    dfloat yRecvElement[4];
    dfloat zRecvElement[4];

    xRecvElement[0] = mesh->EX[recvElement*mesh->Nverts+0];
    xRecvElement[1] = mesh->EX[recvElement*mesh->Nverts+1];
    xRecvElement[2] = mesh->EX[recvElement*mesh->Nverts+2];
    xRecvElement[3] = mesh->EX[recvElement*mesh->Nverts+3];
    
    yRecvElement[0] = mesh->EY[recvElement*mesh->Nverts+0];
    yRecvElement[1] = mesh->EY[recvElement*mesh->Nverts+1];
    yRecvElement[2] = mesh->EY[recvElement*mesh->Nverts+2];
    yRecvElement[3] = mesh->EY[recvElement*mesh->Nverts+3];

    zRecvElement[0] = mesh->EZ[recvElement*mesh->Nverts+0];
    zRecvElement[1] = mesh->EZ[recvElement*mesh->Nverts+1];
    zRecvElement[2] = mesh->EZ[recvElement*mesh->Nverts+2];
    zRecvElement[3] = mesh->EZ[recvElement*mesh->Nverts+3];

    dfloat L1_rec = -(xRecvElement[1]*yRecvElement[2]*zRecvElement[3] - xRecvElement[1]*yRecvElement[2]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[1]*yRecvElement[3]*zRecvElement[2] + xRecvElement[1]*yRecvElement[3]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[1]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[2] - xRecvElement[1]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[3] - xRecvElement[2]*yRecvElement[1]*zRecvElement[3] + xRecvElement[2]*yRecvElement[1]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[2]*yRecvElement[3]*zRecvElement[1] - xRecvElement[2]*yRecvElement[3]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[2]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[1] + xRecvElement[2]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[3] + xRecvElement[3]*yRecvElement[1]*zRecvElement[2] - xRecvElement[3]*yRecvElement[1]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[3]*yRecvElement[2]*zRecvElement[1] + xRecvElement[3]*yRecvElement[2]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[3]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[1] - xRecvElement[3]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[2] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[1]*zRecvElement[2] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[1]*zRecvElement[3] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[2]*zRecvElement[1] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[2]*zRecvElement[3] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[3]*zRecvElement[1] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[3]*zRecvElement[2])/(xRecvElement[0]*yRecvElement[1]*zRecvElement[2] - xRecvElement[0]*yRecvElement[1]*zRecvElement[3] - xRecvElement[0]*yRecvElement[2]*zRecvElement[1] + xRecvElement[0]*yRecvElement[2]*zRecvElement[3] + xRecvElement[0]*yRecvElement[3]*zRecvElement[1] - xRecvElement[0]*yRecvElement[3]*zRecvElement[2] - xRecvElement[1]*yRecvElement[0]*zRecvElement[2] + xRecvElement[1]*yRecvElement[0]*zRecvElement[3] + xRecvElement[1]*yRecvElement[2]*zRecvElement[0] - xRecvElement[1]*yRecvElement[2]*zRecvElement[3] - xRecvElement[1]*yRecvElement[3]*zRecvElement[0] + xRecvElement[1]*yRecvElement[3]*zRecvElement[2] + xRecvElement[2]*yRecvElement[0]*zRecvElement[1] - xRecvElement[2]*yRecvElement[0]*zRecvElement[3] - xRecvElement[2]*yRecvElement[1]*zRecvElement[0] + xRecvElement[2]*yRecvElement[1]*zRecvElement[3] + xRecvElement[2]*yRecvElement[3]*zRecvElement[0] - xRecvElement[2]*yRecvElement[3]*zRecvElement[1] - xRecvElement[3]*yRecvElement[0]*zRecvElement[1] + xRecvElement[3]*yRecvElement[0]*zRecvElement[2] + xRecvElement[3]*yRecvElement[1]*zRecvElement[0] - xRecvElement[3]*yRecvElement[1]*zRecvElement[2] - xRecvElement[3]*yRecvElement[2]*zRecvElement[0] + xRecvElement[3]*yRecvElement[2]*zRecvElement[1]);
    dfloat L2_rec = (xRecvElement[0]*yRecvElement[2]*zRecvElement[3] - xRecvElement[0]*yRecvElement[2]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[0]*yRecvElement[3]*zRecvElement[2] + xRecvElement[0]*yRecvElement[3]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[0]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[2] - xRecvElement[0]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[3] - xRecvElement[2]*yRecvElement[0]*zRecvElement[3] + xRecvElement[2]*yRecvElement[0]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[2]*yRecvElement[3]*zRecvElement[0] - xRecvElement[2]*yRecvElement[3]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[2]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[0] + xRecvElement[2]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[3] + xRecvElement[3]*yRecvElement[0]*zRecvElement[2] - xRecvElement[3]*yRecvElement[0]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[3]*yRecvElement[2]*zRecvElement[0] + xRecvElement[3]*yRecvElement[2]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[3]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[0] - xRecvElement[3]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[2] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[0]*zRecvElement[2] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[0]*zRecvElement[3] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[2]*zRecvElement[0] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[2]*zRecvElement[3] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[3]*zRecvElement[0] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[3]*zRecvElement[2])/(xRecvElement[0]*yRecvElement[1]*zRecvElement[2] - xRecvElement[0]*yRecvElement[1]*zRecvElement[3] - xRecvElement[0]*yRecvElement[2]*zRecvElement[1] + xRecvElement[0]*yRecvElement[2]*zRecvElement[3] + xRecvElement[0]*yRecvElement[3]*zRecvElement[1] - xRecvElement[0]*yRecvElement[3]*zRecvElement[2] - xRecvElement[1]*yRecvElement[0]*zRecvElement[2] + xRecvElement[1]*yRecvElement[0]*zRecvElement[3] + xRecvElement[1]*yRecvElement[2]*zRecvElement[0] - xRecvElement[1]*yRecvElement[2]*zRecvElement[3] - xRecvElement[1]*yRecvElement[3]*zRecvElement[0] + xRecvElement[1]*yRecvElement[3]*zRecvElement[2] + xRecvElement[2]*yRecvElement[0]*zRecvElement[1] - xRecvElement[2]*yRecvElement[0]*zRecvElement[3] - xRecvElement[2]*yRecvElement[1]*zRecvElement[0] + xRecvElement[2]*yRecvElement[1]*zRecvElement[3] + xRecvElement[2]*yRecvElement[3]*zRecvElement[0] - xRecvElement[2]*yRecvElement[3]*zRecvElement[1] - xRecvElement[3]*yRecvElement[0]*zRecvElement[1] + xRecvElement[3]*yRecvElement[0]*zRecvElement[2] + xRecvElement[3]*yRecvElement[1]*zRecvElement[0] - xRecvElement[3]*yRecvElement[1]*zRecvElement[2] - xRecvElement[3]*yRecvElement[2]*zRecvElement[0] + xRecvElement[3]*yRecvElement[2]*zRecvElement[1]);
    dfloat L3_rec = -(xRecvElement[0]*yRecvElement[1]*zRecvElement[3] - xRecvElement[0]*yRecvElement[1]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[0]*yRecvElement[3]*zRecvElement[1] + xRecvElement[0]*yRecvElement[3]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[0]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[1] - xRecvElement[0]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[3] - xRecvElement[1]*yRecvElement[0]*zRecvElement[3] + xRecvElement[1]*yRecvElement[0]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[1]*yRecvElement[3]*zRecvElement[0] - xRecvElement[1]*yRecvElement[3]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[1]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[0] + xRecvElement[1]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[3] + xRecvElement[3]*yRecvElement[0]*zRecvElement[1] - xRecvElement[3]*yRecvElement[0]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[3]*yRecvElement[1]*zRecvElement[0] + xRecvElement[3]*yRecvElement[1]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[3]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[0] - xRecvElement[3]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[1] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[0]*zRecvElement[1] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[0]*zRecvElement[3] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[1]*zRecvElement[0] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[1]*zRecvElement[3] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[3]*zRecvElement[0] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[3]*zRecvElement[1])/(xRecvElement[0]*yRecvElement[1]*zRecvElement[2] - xRecvElement[0]*yRecvElement[1]*zRecvElement[3] - xRecvElement[0]*yRecvElement[2]*zRecvElement[1] + xRecvElement[0]*yRecvElement[2]*zRecvElement[3] + xRecvElement[0]*yRecvElement[3]*zRecvElement[1] - xRecvElement[0]*yRecvElement[3]*zRecvElement[2] - xRecvElement[1]*yRecvElement[0]*zRecvElement[2] + xRecvElement[1]*yRecvElement[0]*zRecvElement[3] + xRecvElement[1]*yRecvElement[2]*zRecvElement[0] - xRecvElement[1]*yRecvElement[2]*zRecvElement[3] - xRecvElement[1]*yRecvElement[3]*zRecvElement[0] + xRecvElement[1]*yRecvElement[3]*zRecvElement[2] + xRecvElement[2]*yRecvElement[0]*zRecvElement[1] - xRecvElement[2]*yRecvElement[0]*zRecvElement[3] - xRecvElement[2]*yRecvElement[1]*zRecvElement[0] + xRecvElement[2]*yRecvElement[1]*zRecvElement[3] + xRecvElement[2]*yRecvElement[3]*zRecvElement[0] - xRecvElement[2]*yRecvElement[3]*zRecvElement[1] - xRecvElement[3]*yRecvElement[0]*zRecvElement[1] + xRecvElement[3]*yRecvElement[0]*zRecvElement[2] + xRecvElement[3]*yRecvElement[1]*zRecvElement[0] - xRecvElement[3]*yRecvElement[1]*zRecvElement[2] - xRecvElement[3]*yRecvElement[2]*zRecvElement[0] + xRecvElement[3]*yRecvElement[2]*zRecvElement[1]);
    dfloat L4_rec = (xRecvElement[0]*yRecvElement[1]*zRecvElement[2] - xRecvElement[0]*yRecvElement[1]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[0]*yRecvElement[2]*zRecvElement[1] + xRecvElement[0]*yRecvElement[2]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[0]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[1] - xRecvElement[0]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[2] - xRecvElement[1]*yRecvElement[0]*zRecvElement[2] + xRecvElement[1]*yRecvElement[0]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[1]*yRecvElement[2]*zRecvElement[0] - xRecvElement[1]*yRecvElement[2]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[1]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[0] + xRecvElement[1]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[2] + xRecvElement[2]*yRecvElement[0]*zRecvElement[1] - xRecvElement[2]*yRecvElement[0]*acoustics->recvXYZ[2+rIdx*3] - xRecvElement[2]*yRecvElement[1]*zRecvElement[0] + xRecvElement[2]*yRecvElement[1]*acoustics->recvXYZ[2+rIdx*3] + xRecvElement[2]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[0] - xRecvElement[2]*acoustics->recvXYZ[1+rIdx*3]*zRecvElement[1] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[0]*zRecvElement[1] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[0]*zRecvElement[2] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[1]*zRecvElement[0] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[1]*zRecvElement[2] - acoustics->recvXYZ[0+rIdx*3]*yRecvElement[2]*zRecvElement[0] + acoustics->recvXYZ[0+rIdx*3]*yRecvElement[2]*zRecvElement[1])/(xRecvElement[0]*yRecvElement[1]*zRecvElement[2] - xRecvElement[0]*yRecvElement[1]*zRecvElement[3] - xRecvElement[0]*yRecvElement[2]*zRecvElement[1] + xRecvElement[0]*yRecvElement[2]*zRecvElement[3] + xRecvElement[0]*yRecvElement[3]*zRecvElement[1] - xRecvElement[0]*yRecvElement[3]*zRecvElement[2] - xRecvElement[1]*yRecvElement[0]*zRecvElement[2] + xRecvElement[1]*yRecvElement[0]*zRecvElement[3] + xRecvElement[1]*yRecvElement[2]*zRecvElement[0] - xRecvElement[1]*yRecvElement[2]*zRecvElement[3] - xRecvElement[1]*yRecvElement[3]*zRecvElement[0] + xRecvElement[1]*yRecvElement[3]*zRecvElement[2] + xRecvElement[2]*yRecvElement[0]*zRecvElement[1] - xRecvElement[2]*yRecvElement[0]*zRecvElement[3] - xRecvElement[2]*yRecvElement[1]*zRecvElement[0] + xRecvElement[2]*yRecvElement[1]*zRecvElement[3] + xRecvElement[2]*yRecvElement[3]*zRecvElement[0] - xRecvElement[2]*yRecvElement[3]*zRecvElement[1] - xRecvElement[3]*yRecvElement[0]*zRecvElement[1] + xRecvElement[3]*yRecvElement[0]*zRecvElement[2] + xRecvElement[3]*yRecvElement[1]*zRecvElement[0] - xRecvElement[3]*yRecvElement[1]*zRecvElement[2] - xRecvElement[3]*yRecvElement[2]*zRecvElement[0] + xRecvElement[3]*yRecvElement[2]*zRecvElement[1]);

    //Vandermonde Berstein matrix in receiver point
    dfloat *VB_rec;
    VB_rec = (dfloat*) calloc(mesh->Np, sizeof(dfloat));

    int sk = 0;
    for(int l = 0; l <= mesh->N; l++){
      for(int k = 0; k <= mesh->N - l; k++){
        for(int j = 0; j <= mesh->N - k - l; j++){
          int i = mesh->N - j - k - l;
          dfloat temp = factorial(mesh->N)/(factorial(i)*factorial(j)*factorial(k)*factorial(l));
          VB_rec[sk] = temp*pow(L1_rec,i)*pow(L2_rec,j)*pow(L3_rec,k)*pow(L4_rec,l);
          sk++;
        }
      }
    }
    // interpolation   
    for(int j = 0; j < mesh->Np; j++){
      for(int i = 0; i < mesh->Np; i++){
        intpol[iRecv*mesh->Np+i] += VB_rec[j]*mesh->invVB[i+j*mesh->Np];
      }
    }
    #if 0
    printf("new, irecv = %d, ele = %d\n",iRecv,recvElement);
    for(int i = 0; i < mesh->Np; i++){
      printf("%.15lf \n",intpol[iRecv*mesh->Np+i]);
    }
    printf("\n");
    #endif
    free(VB_rec);
  }

  acoustics->o_recvintpol = 
        mesh->device.malloc(acoustics->NReceiversLocal*mesh->Np*sizeof(dfloat),intpol);

  free(intpol);
}



void acousticsFindReceiverElement(acoustics_t *acoustics){
  // Only works for tets!
  // [TODO] If receiver point is on the boundary of two cores both will find it! Do some mpi stuff to fix!
  mesh_t *mesh = acoustics->mesh;

  for(dlong k = 0; k < acoustics->NReceivers; k++){

    dfloat recvLoc[3];
    recvLoc[0] = acoustics->recvXYZ[k*3+0];
    recvLoc[1] = acoustics->recvXYZ[k*3+1];
    recvLoc[2] = acoustics->recvXYZ[k*3+2];

    dlong faceVertices[4][4] = {{0,1,2,3},{0,1,3,2},{1,2,3,0},{2,0,3,1}};

    for(dlong i = 0; i < mesh->Nelements; i++){
      // Assume receiver is in element
      dlong isInside = 1;
      for(int j = 0; j < mesh->Nfaces; j++){
        dlong fv1 = faceVertices[j][0];
        dlong fv2 = faceVertices[j][1];
        dlong fv3 = faceVertices[j][2];
        dlong fv4 = faceVertices[j][3];

        // b r s defines the plane
        dfloat b[3] = {mesh->EX[i*mesh->Nverts + fv1],
                      mesh->EY[i*mesh->Nverts + fv1],
                      mesh->EZ[i*mesh->Nverts + fv1]};
        dfloat r[3] = {mesh->EX[i*mesh->Nverts + fv2],
                      mesh->EY[i*mesh->Nverts + fv2],
                      mesh->EZ[i*mesh->Nverts + fv2]};
        dfloat s[3] = {mesh->EX[i*mesh->Nverts + fv3],
                      mesh->EY[i*mesh->Nverts + fv3],
                      mesh->EZ[i*mesh->Nverts + fv3]};
        
        // d is the control point (the last point in the tet)
        dfloat d[3] = {mesh->EX[i*mesh->Nverts + fv4],
                      mesh->EY[i*mesh->Nverts + fv4],
                      mesh->EZ[i*mesh->Nverts + fv4]};


        // Cross product to get normal vector of plane brs
        dfloat n[3] = {(r[1]-b[1])*(s[2]-b[2]) - (r[2]-b[2])*(s[1]-b[1]),
        (r[2]-b[2])*(s[0]-b[0]) - (r[0]-b[0])*(s[2]-b[2]),
        (r[0]-b[0])*(s[1]-b[1]) - (r[1]-b[1])*(s[0]-b[0])};


        // Calculate the plane equation for both receiver and leftover tet point.
        dfloat planeEqRecv = n[0]*(recvLoc[0]-b[0]) + n[1]*(recvLoc[1]-b[1]) + n[2]*(recvLoc[2]-b[2]);
        dfloat planeEqOther = n[0]*(d[0]-b[0]) + n[1]*(d[1]-b[1]) + n[2]*(d[2]-b[2]);

        // Check if the two points are on the same side of the plane
        dlong recvSide = planeEqRecv > 0 ? 1 : -1;
        dlong otherSide = planeEqOther > 0 ? 1 : -1;
        planeEqRecv = planeEqRecv >= 0 ? planeEqRecv:-1.0*planeEqRecv;
        if(recvSide != otherSide && planeEqRecv > 1.0e-15){
          // Recv is not inside element i
          isInside = 0;
          break;
        }
      }
      //Check if found recv point inside element i
      if(isInside == 1 ){
        acoustics->recvElements[k] = i;
        acoustics->recvElementsIdx[acoustics->NReceiversLocal] = k;
        acoustics->NReceiversLocal++;
        break;
      } else if(i == mesh->Nelements-1){
          // [EA] Currently prints if the receiver is not on this core! Even if found by another core!
          //printf("RECEIVER LOCATION NOT FOUND!!,(x,y,z) = (%f,%f,%f)\n", recvLoc[0],recvLoc[1],recvLoc[2]);
      }
    }
  }
}

void acousticsPrintReceiversToFile(acoustics_t *acoustics, setupAide &newOptions){
  mesh_t *mesh = acoustics->mesh;
  dfloat sloc[3];
  dfloat sxyz;
  string outDir_cppstring;

  newOptions.getArgs("OUTPUT DIRECTORY", outDir_cppstring);
  newOptions.getArgs("SX", sloc[0]);
  newOptions.getArgs("SY", sloc[1]);
  newOptions.getArgs("SZ", sloc[2]);
  newOptions.getArgs("SXYZ", sxyz);

  char *outDir = (char*)outDir_cppstring.c_str();

  string PREFIX;
  newOptions.getArgs("RECEIVERPREFIX", PREFIX);
  for(dlong iRecv = 0; iRecv < acoustics->NReceiversLocal; iRecv++){
    // Print interpolated receiver to file
    
    char fname[BUFSIZ];

    // create output folder if not existing
    struct stat st = {0};
    if (stat(outDir, &st) == -1) {
      mkdir(outDir, 0700);
    }

    sprintf(fname, "%s/%s_RecvPoint_%02d.txt", outDir, (char*)PREFIX.c_str(), acoustics->recvElementsIdx[iRecv]);

    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) == NULL) {
      perror("getcwd() error");
      return;
    }

    FILE *iFP = fopen(fname,"w");
    if (iFP == NULL) {
      printf("ERROR: receiver output file could not be opened %s/%s)\n", cwd, fname);      
      return;
    }

    dfloat time = 0;
    dfloat u = 0, v = 0, w = 0, r = 0;

    dlong rIdx = acoustics->recvElementsIdx[iRecv];
    dfloat x = acoustics->recvXYZ[rIdx*3+0];
    dfloat y = acoustics->recvXYZ[rIdx*3+1];
    dfloat z = acoustics->recvXYZ[rIdx*3+2];
    
    acousticsGaussianPulse(x, y, z, 0, &r, &u, &v, &w, sloc, sxyz);
    
    fprintf(iFP, "%.15lf %.15le\n", time, r);
    
    for(int i = 0; i < mesh->NtimeSteps; i++){
      time += mesh->dt;
      fprintf(iFP, "%.15lf %.15le\n", time, acoustics->qRecv[i+iRecv*mesh->NtimeSteps]);
    }

    printf("Receiver impulse response was written to disk: %s/%s\n", cwd, fname);

    fclose(iFP);
  }
}


