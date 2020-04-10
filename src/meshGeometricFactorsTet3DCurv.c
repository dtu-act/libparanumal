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

#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"


void calcVGeoFacs(mesh3D *mesh, int e, dfloat *geo){
	dfloat *Dr = mesh->Dr;
	dfloat *Ds = mesh->Ds;
	dfloat *Dt = mesh->Dt;

	int id = e*mesh->Np;
	// Dr
	int offset = 0;
	for(int i = 0; i < mesh->Np; i++){
		int idx = offset + i;
		int idy = offset + i + 3*mesh->Np;
		int idz = offset + i + 6*mesh->Np;
		
		geo[idx] = 0;
		geo[idy] = 0;
		geo[idz] = 0;
		for(int j = 0; j < mesh->Np; j++){
			geo[idx] += Dr[j+i*mesh->Np] * mesh->x[id+j];
			geo[idy] += Dr[j+i*mesh->Np] * mesh->y[id+j];
			geo[idz] += Dr[j+i*mesh->Np] * mesh->z[id+j];
		}
	}

	// Ds
	offset = mesh->Np;
	for(int i = 0; i < mesh->Np; i++){
		int idx = offset + i;
		int idy = offset + i + 3*mesh->Np;
		int idz = offset + i + 6*mesh->Np;
		
		geo[idx] = 0;
		geo[idy] = 0;
		geo[idz] = 0;
		for(int j = 0; j < mesh->Np; j++){
			geo[idx] += Ds[j+i*mesh->Np] * mesh->x[id+j];
			geo[idy] += Ds[j+i*mesh->Np] * mesh->y[id+j];
			geo[idz] += Ds[j+i*mesh->Np] * mesh->z[id+j];
		}
	}

	// Dt
	offset = 2*mesh->Np;
	for(int i = 0; i < mesh->Np; i++){
		int idx = offset + i;
		int idy = offset + i + 3*mesh->Np;
		int idz = offset + i + 6*mesh->Np;
		
		geo[idx] = 0;
		geo[idy] = 0;
		geo[idz] = 0;
		for(int j = 0; j < mesh->Np; j++){
			geo[idx] += Dt[j+i*mesh->Np] * mesh->x[id+j];
			geo[idy] += Dt[j+i*mesh->Np] * mesh->y[id+j];
			geo[idz] += Dt[j+i*mesh->Np] * mesh->z[id+j];
		}
	}

}


void meshGeometricFactorsTet3DCurv(mesh3D *mesh){

	// Initialized to max possible size, reallocated to correct size at the end of the function
	mesh->NvgeoCurv = 10;
	mesh->vgeoCurv = (dfloat*) calloc(mesh->Nelements*mesh->NvgeoCurv*mesh->Np, sizeof(dfloat));
	
	mesh->mapCurv = (dlong*) calloc(mesh->Nelements, sizeof(dlong));
	mesh->Ncurv = 0;	

	// Temp storage
	dfloat *geo = (dfloat*) calloc(9*mesh->Np, sizeof(dfloat));
	dfloat *Jtmp = (dfloat*) calloc(mesh->Np, sizeof(dfloat));

	for(dlong e = 0; e < mesh->Nelements; e++){
		calcVGeoFacs(mesh, e, geo);
		for(int n = 0; n < mesh->Np;n++){
				dfloat xr = geo[n], xs = geo[mesh->Np + n], xt = geo[2*mesh->Np + n];
				dfloat yr = geo[3*mesh->Np + n], ys = geo[4*mesh->Np + n], yt = geo[5*mesh->Np + n];
				dfloat zr = geo[6*mesh->Np + n], zs = geo[7*mesh->Np + n], zt = geo[8*mesh->Np + n];

				Jtmp[n] = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);
		}
		dfloat maxDif = 0;
		for(dlong n = 0; n < mesh->Np; n++){
			for(dlong nn = 0; nn < mesh->Np; nn++){
				dfloat tmp = fabs(Jtmp[n] - Jtmp[nn]);
				if(tmp > maxDif){
					maxDif = tmp;
				}
			}
		}
		mesh->mapCurv[e] = -1;
		if(maxDif > 1.0e-13){ // Check if element is curvilinear
			mesh->mapCurv[e] = mesh->Ncurv;
			mesh->Ncurv++;
		}
		dlong eCurv = mesh->mapCurv[e];
		if(eCurv > -1){
			for(dlong n = 0; n < mesh->Np; n++){
				dfloat xr = geo[n], xs = geo[mesh->Np + n], xt = geo[2*mesh->Np + n];
				dfloat yr = geo[3*mesh->Np + n], ys = geo[4*mesh->Np + n], yt = geo[5*mesh->Np + n];
				dfloat zr = geo[6*mesh->Np + n], zs = geo[7*mesh->Np + n], zt = geo[8*mesh->Np + n];

				dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);
				dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
				dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
				dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;

				const dlong vid = mesh->NvgeoCurv*eCurv*mesh->Np + n*mesh->NvgeoCurv;
				mesh->vgeoCurv[vid + RXIDC] = rx;
				mesh->vgeoCurv[vid + RYIDC] = ry;
				mesh->vgeoCurv[vid + RZIDC] = rz;
				mesh->vgeoCurv[vid + SXIDC] = sx;
				mesh->vgeoCurv[vid + SYIDC] = sy;
				mesh->vgeoCurv[vid + SZIDC] = sz;
				mesh->vgeoCurv[vid + TXIDC] = tx;
				mesh->vgeoCurv[vid + TYIDC] = ty;
				mesh->vgeoCurv[vid + TZIDC] = tz;
				mesh->vgeoCurv[vid +  JIDC] = J;
				
				if(J<0) printf("bugger: got negative curvilinear geofac\n");
				
			}
		}
	}

	// Resize array to the correct size
	mesh->vgeoCurv = (dfloat*) realloc(mesh->vgeoCurv, mesh->Ncurv*mesh->NvgeoCurv*mesh->Np*sizeof(dfloat));
	free(geo);
	free(Jtmp);
}
