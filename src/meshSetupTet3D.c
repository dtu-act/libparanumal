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

#include "mesh3D.h"

mesh3D *meshSetupTet3D(char *filename, int N){

  // read chunk of elements
  mesh3D *mesh = meshParallelReaderTet3D(filename);

  // partition elements using Morton ordering & parallel sort
  meshGeometricPartition3D(mesh);
  
  // connect elements using parallel sort
  meshParallelConnect(mesh);

  // print out connectivity statistics
  meshPartitionStatistics(mesh);

  // connect elements to boundary faces
  meshConnectBoundary(mesh);

  // load reference (r,s,t) element nodes
  meshLoadReferenceNodesTet3D(mesh, N);

  // compute physical (x,y,z) locations of the element nodes
  meshPhysicalNodesTet3D(mesh);

  // compute geometric factors
  meshGeometricFactorsTet3D(mesh);

  // set up halo exchange info for MPI (do before connect face nodes)
  // [AE] Halo nodes are u^+ on other procs
  meshHaloSetup(mesh);
  
  // connect face nodes (find trace indices)
  meshConnectFaceNodes3D(mesh);

  // compute surface geofacs
  meshSurfaceGeometricFactorsTet3D(mesh);

  // global nodes
  meshParallelConnectNodes(mesh);

  // initialize LSERK4 time stepping coefficients
  int Nrk = 5;

  dfloat rka[5] = {0.0,
		   -567301805773.0/1357537059087.0 ,
		   -2404267990393.0/2016746695238.0 ,
		   -3550918686646.0/2091501179385.0  ,
		   -1275806237668.0/842570457699.0};
  dfloat rkb[5] = { 1432997174477.0/9575080441755.0 ,
		    5161836677717.0/13612068292357.0 ,
		    1720146321549.0/2090206949498.0  ,
		    3134564353537.0/4481467310338.0  ,
		    2277821191437.0/14882151754819.0};
  dfloat rkc[5] = {0.0  ,
		   1432997174477.0/9575080441755.0 ,
		   2526269341429.0/6820363962896.0 ,
		   2006345519317.0/3224310063776.0 ,
		   2802321613138.0/2924317926251.0}; 

  mesh->Nrk = Nrk;
  memcpy(mesh->rka, rka, Nrk*sizeof(dfloat));
  memcpy(mesh->rkb, rkb, Nrk*sizeof(dfloat));
  memcpy(mesh->rkc, rkc, Nrk*sizeof(dfloat));

  // [EA] initialize explicit-implicit time stepping coefficients for Local Reaction, EIRK4 (Explicit Implicit Runge-Kutta 4th order)
  int INrk = 6;
  dfloat erka[36] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        1.0/2.0,0.0,0.0,0.0,0.0,0.0,
        13861.0/62500.0,6889.0/62500.0,0.0,0.0,0.0,0.0,
        -116923316275.0/2393684061468.0,-2731218467317.0/15368042101831.0,9408046702089.0/11113171139209.0,0.0,0.0,0.0,
        -451086348788.0/2902428689909.0,-2682348792572.0/7519795681897.0,12662868775082.0/11960479115383.0,3355817975965.0/11060851509271.0,0.0,0.0,
        647845179188.0/3216320057751.0,73281519250.0/8382639484533.0,552539513391.0/3454668386233.0,3354512671639.0/8306763924573.0,4040.0/17871.0,0.0
        };
  dfloat erkb[6] = {82889.0/524892.0,0.0,15625.0/83664.0,69875.0/102672.0,-2260.0/8211.0,1.0/4.0};
  dfloat erkc[6] = {0.0,1.0/2.0,83.0/250.0,31.0/50.0,17.0/20.0,1.0};
  dfloat erke[6] = {31666707/9881966720, 0, -256875/105007616, -2768025/128864768, 169839/3864644, -5247/225920};
  dfloat esdirka[36] = {
        0.0,0.0,0.0,0.0,0.0,0.0,
        1.0/4.0,1.0/4.0,0.0,0.0,0.0,0,
        8611.0/62500.0,-1743.0/31250.0,1.0/4.0,0.0,0.0,0,
        5012029.0/34652500.0,-654441.0/2922500.0,174375.0/388108.0,1.0/4.0,0.0,0,
        15267082809.0/155376265600.0,-71443401.0/120774400.0,730878875.0/902184768.0,2285395.0/8070912.0,1.0/4.0,0,
        82889.0/524892.0,0.0,15625.0/83664.0,69875.0/102672.0,-2260.0/8211.0,1.0/4
        };
  dfloat esdirkb[6] = {82889.0/524892.0,0.0,15625.0/83664.0,69875.0/102672.0,-2260.0/8211.0,1.0/4.0};
  dfloat esdirkc[6] = {0.0,1.0/2.0,83.0/250.0,31.0/50.0,17.0/20.0,1.0};
  dfloat esdirke[6] = {31666707/9881966720, 0, -256875/105007616, -2768025/128864768, 169839/3864644, -5247/225920};

  mesh->INrk = INrk;
  memcpy(mesh->erka, erka, INrk*INrk*sizeof(dfloat));
  memcpy(mesh->erkb, erkb, INrk*sizeof(dfloat));
  memcpy(mesh->erkc, erkc, INrk*sizeof(dfloat));
  memcpy(mesh->erke, erke, INrk*sizeof(dfloat));
  memcpy(mesh->esdirka, esdirka, INrk*INrk*sizeof(dfloat));
  memcpy(mesh->esdirkb, esdirkb, INrk*sizeof(dfloat));
  memcpy(mesh->esdirkc, esdirkc, INrk*sizeof(dfloat));
  memcpy(mesh->esdirke, esdirke, INrk*sizeof(dfloat));


#if 0
  for(int eM=0;eM<mesh->Nelements;++eM){
    for(int fM=0;fM<mesh->Nfaces;++fM){

      int eP = mesh->EToE[eM*mesh->Nfaces+fM];
      int fP = mesh->EToF[eM*mesh->Nfaces+fM];

      if(eP<mesh->Nelements){
	if(eP<0){ eP = eM; fP = fM; }
      
	dfloat sc = 1;
	if(eM==eP) sc = -1;
      
	int baseM = (eM*mesh->Nfaces + fM)*mesh->Nsgeo;
	int baseP = (eP*mesh->Nfaces + fP)*mesh->Nsgeo;

	dfloat avenx = mesh->sgeo[baseM+NXID] + sc*mesh->sgeo[baseP+NXID];
	dfloat aveny = mesh->sgeo[baseM+NYID] + sc*mesh->sgeo[baseP+NYID];
	dfloat avenz = mesh->sgeo[baseM+NZID] + sc*mesh->sgeo[baseP+NZID];
	dfloat dsJ   = mesh->sgeo[baseM+SJID] - mesh->sgeo[baseP+SJID];
	dfloat aven  = norm3(avenx,aveny,avenz);

	if(aven>1e-4)
	  printf("(%d,%d=>%d,%d) aven = %g,%g,%g", eM,fM, eP,fP, avenx,aveny,avenz);

	if(fabs(dsJ)>1e-4)
	  printf("sJ mismatch %g\n", dsJ);
	
	for(int n=0;n<mesh->Nfp;++n){
	  int id = n + fM*mesh->Nfp + eM*mesh->Nfp*mesh->Nfaces;
	  int idM = mesh->vmapM[id];
	  int idP = mesh->vmapP[id];

	  dfloat d = norm3(mesh->x[idM]-mesh->x[idP],mesh->y[idM]-mesh->y[idP],mesh->z[idM]-mesh->z[idP]);

	  if(d>1e-4)
	    printf(", %g\n", d);
	}
      }
    }
  }
#endif  
  return mesh;
}
