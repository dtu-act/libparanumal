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

THE SOFTWARE IS PROVIDCED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

  
#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"


void meshSurfaceGeometricFactorsTet3DCurv(mesh3D *mesh){
	
	mesh->NsgeoCurv = 5;
	mesh->sgeoCurv = (dfloat*) calloc(mesh->Ncurv*mesh->NsgeoCurv*mesh->Nfaces*mesh->Nfp, sizeof(dfloat));
	
	for(int e = 0; e < mesh->Ncurv; e++){
		for(int n = 0; n < mesh->Nfp; n++){
			
			// Face 1
			int f = 0;
			int fn = mesh->faceNodes[f*mesh->Nfp+n];
			int id = mesh->NvgeoCurv*e*mesh->Np + fn*mesh->NvgeoCurv;
			dfloat J = mesh->vgeoCurv[id+JIDC];
			dfloat rx =  mesh->vgeoCurv[id+RXIDC], ry = mesh->vgeoCurv[id+RYIDC], rz = mesh->vgeoCurv[id+RZIDC];
			dfloat sx =  mesh->vgeoCurv[id+SXIDC], sy = mesh->vgeoCurv[id+SYIDC], sz = mesh->vgeoCurv[id+SZIDC];
			dfloat tx =  mesh->vgeoCurv[id+TXIDC], ty = mesh->vgeoCurv[id+TYIDC], tz = mesh->vgeoCurv[id+TZIDC];

			dfloat nx = -tx;
			dfloat ny = -ty;
			dfloat nz = -tz;
			dfloat sJ = norm3(nx,ny,nz);
			
			dlong base = e*mesh->NsgeoCurv*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp*mesh->NsgeoCurv + n*mesh->NsgeoCurv;
			mesh->sgeoCurv[base+NXIDC] = nx/sJ;
			mesh->sgeoCurv[base+NYIDC] = ny/sJ;
			mesh->sgeoCurv[base+NZIDC] = nz/sJ;
			mesh->sgeoCurv[base+SJIDC] = sJ*J;
			mesh->sgeoCurv[base+IJIDC] = 1./J;

			
			// Face 2
			f = 1;
			fn = mesh->faceNodes[f*mesh->Nfp+n];
			id = mesh->NvgeoCurv*e*mesh->Np + fn*mesh->NvgeoCurv;
			J = mesh->vgeoCurv[id+JIDC];
			rx =  mesh->vgeoCurv[id+RXIDC], ry = mesh->vgeoCurv[id+RYIDC], rz = mesh->vgeoCurv[id+RZIDC];
			sx =  mesh->vgeoCurv[id+SXIDC], sy = mesh->vgeoCurv[id+SYIDC], sz = mesh->vgeoCurv[id+SZIDC];
			tx =  mesh->vgeoCurv[id+TXIDC], ty = mesh->vgeoCurv[id+TYIDC], tz = mesh->vgeoCurv[id+TZIDC];

			nx = -sx;
			ny = -sy;
			nz = -sz;
			sJ = norm3(nx,ny,nz);
			
			base = e*mesh->NsgeoCurv*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp*mesh->NsgeoCurv + n*mesh->NsgeoCurv;
			mesh->sgeoCurv[base+NXIDC] = nx/sJ;
			mesh->sgeoCurv[base+NYIDC] = ny/sJ;
			mesh->sgeoCurv[base+NZIDC] = nz/sJ;
			mesh->sgeoCurv[base+SJIDC] = sJ*J;
			mesh->sgeoCurv[base+IJIDC] = 1./J;

			// Face 3
			f = 2;
			fn = mesh->faceNodes[f*mesh->Nfp+n];
			id = mesh->NvgeoCurv*e*mesh->Np + fn*mesh->NvgeoCurv;
			J = mesh->vgeoCurv[id+JIDC];
			rx =  mesh->vgeoCurv[id+RXIDC], ry = mesh->vgeoCurv[id+RYIDC], rz = mesh->vgeoCurv[id+RZIDC];
			sx =  mesh->vgeoCurv[id+SXIDC], sy = mesh->vgeoCurv[id+SYIDC], sz = mesh->vgeoCurv[id+SZIDC];
			tx =  mesh->vgeoCurv[id+TXIDC], ty = mesh->vgeoCurv[id+TYIDC], tz = mesh->vgeoCurv[id+TZIDC];
			
			nx = rx+sx+tx;
			ny = ry+sy+ty;
			nz = rz+sz+tz;
			sJ = norm3(nx,ny,nz);
			
			base = e*mesh->NsgeoCurv*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp*mesh->NsgeoCurv + n*mesh->NsgeoCurv;
			mesh->sgeoCurv[base+NXIDC] = nx/sJ;
			mesh->sgeoCurv[base+NYIDC] = ny/sJ;
			mesh->sgeoCurv[base+NZIDC] = nz/sJ;
			mesh->sgeoCurv[base+SJIDC] = sJ*J;
			mesh->sgeoCurv[base+IJIDC] = 1./J;

			// Face 4
			f = 3;
			fn = mesh->faceNodes[f*mesh->Nfp+n];
			id = mesh->NvgeoCurv*e*mesh->Np + fn*mesh->NvgeoCurv;
			J = mesh->vgeoCurv[id+JIDC];
			rx =  mesh->vgeoCurv[id+RXIDC], ry = mesh->vgeoCurv[id+RYIDC], rz = mesh->vgeoCurv[id+RZIDC];
			sx =  mesh->vgeoCurv[id+SXIDC], sy = mesh->vgeoCurv[id+SYIDC], sz = mesh->vgeoCurv[id+SZIDC];
			tx =  mesh->vgeoCurv[id+TXIDC], ty = mesh->vgeoCurv[id+TYIDC], tz = mesh->vgeoCurv[id+TZIDC];

			nx = -rx;
			ny = -ry;
			nz = -rz;
			sJ = norm3(nx,ny,nz);
			
			base = e*mesh->NsgeoCurv*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp*mesh->NsgeoCurv + n*mesh->NsgeoCurv;
			mesh->sgeoCurv[base+NXIDC] = nx/sJ;
			mesh->sgeoCurv[base+NYIDC] = ny/sJ;
			mesh->sgeoCurv[base+NZIDC] = nz/sJ;
			mesh->sgeoCurv[base+SJIDC] = sJ*J;
			mesh->sgeoCurv[base+IJIDC] = 1./J;

		}
	}
}
