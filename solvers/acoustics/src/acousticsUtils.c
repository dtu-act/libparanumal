#include <random>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <unordered_map>
#include <map>
#include <cmath>

#include "acoustics.h"
#include "acousticsUtils.h"

std::string generateUUID(int length)
{
    std::ostringstream out;

    for (int i=0; i<length; i++) {
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> dist6(0,9);

        out << dist6(rng);
    }
    return out.str();
}

/*
   purpose: read gmsh nodes (only)
*/
void pointReader3D(char *fileName, std::vector<dfloat> &VX, std::vector<dfloat> &VY, std::vector<dfloat> &VZ)
{
  int Nnodes;

  FILE *fp = fopen(fileName, "r");

  char *status;

  if (fp == NULL)
  {
    printf("meshReaderTet3D: could not load file %s\n", fileName);
    exit(0);
  }

  char buf[BUFSIZ];
  do
  {
    status = fgets(buf, BUFSIZ, fp);
  } while (!strstr(buf, "$Nodes"));

  /* read number of nodes in mesh */
  status = fgets(buf, BUFSIZ, fp);
  sscanf(buf, hlongFormat, &Nnodes);

  VX.resize(Nnodes);
  VY.resize(Nnodes);
  VZ.resize(Nnodes);

  dfloat VX_;
  dfloat VY_;
  dfloat VZ_;

  /* load nodes */
  for (hlong n = 0; n < Nnodes; ++n)
  {
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d" dfloatFormat dfloatFormat dfloatFormat,
           &VX_, &VY_, &VZ_);
    
    VX[n] = VX_;
    VY[n] = VY_;
    VZ[n] = VZ_;
  }

  fclose(fp);
}

void extractUniquePoints(mesh_t *mesh, acoustics_t *acoustics, 
  std::vector<uint> &conn, std::vector<dfloat> &x1d, std::vector<dfloat> &y1d, std::vector<dfloat> &z1d, std::vector<dfloat> &p1d)
{  
  struct Coord3D { dfloat x, y, z; };

  // https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key
  struct Coord3DHasher
  {
    std::size_t operator()(const Coord3D &coord) const
    {
      return ((std::hash<std::string>()(toString(coord.x)) ^ 
        (std::hash<std::string>()(toString(coord.y)) << 1)) >> 1) ^ 
        (std::hash<std::string>()(toString(coord.z)) << 1);
    }
  };

  dfloat t = 0;

  auto comp = [](const Coord3D &c1, const Coord3D &c2)
  {
    const dfloat eps = 1e-5; // std::numeric_limits<float>::epsilon();
    return (std::abs(c1.x - c2.x) < eps) && (std::abs(c1.y - c2.y) < eps) && (std::abs(c1.z - c2.z) < eps);
  };

  std::unordered_map<Coord3D, int, Coord3DHasher, decltype(comp)> coordsToIndx;
  conn.resize(mesh->Nelements * mesh->Np);

  for (dlong e = 0; e < mesh->Nelements; ++e)
  {
    for (int n = 0; n < mesh->Np; ++n)
    {
      dfloat x = mesh->x[n + mesh->Np * e];
      dfloat y = mesh->y[n + mesh->Np * e];
      dfloat z = mesh->z[n + mesh->Np * e];

      auto coord_ = Coord3D{x, y, z};
      auto ret = coordsToIndx.insert({coord_, x1d.size()});

      if (ret.second)
      {
        // a new coordinate set was found
        x1d.push_back(x);
        y1d.push_back(y);
        z1d.push_back(z);

        dlong qbase = e * mesh->Np * mesh->Nfields + n;
        p1d.push_back(acoustics->q[qbase + 0 * mesh->Np]);
      }

      conn[n + mesh->Np * e] = ret.first->second;
    }
  }
}

void meshMinMax(mesh_t *mesh, dfloat xminmax[2], dfloat yminmax[2], dfloat zminmax[2]) {
  xminmax[0] = std::numeric_limits<float>::max();
  xminmax[1] = std::numeric_limits<float>::min();
  yminmax[0] = std::numeric_limits<float>::max();
  yminmax[1] = std::numeric_limits<float>::min();
  zminmax[0] = std::numeric_limits<float>::max();
  zminmax[1] = std::numeric_limits<float>::min();

  for (dlong e = 0; e < mesh->Nelements; ++e)
  {
    for (int n = 0; n < mesh->Np; ++n)
    {
      dfloat x = mesh->x[n + mesh->Np * e];
      dfloat y = mesh->y[n + mesh->Np * e];
      dfloat z = mesh->z[n + mesh->Np * e];

      xminmax[0] = (x < xminmax[0]) ? x : xminmax[0];
      xminmax[1] = (x > xminmax[1]) ? x : xminmax[1];
      yminmax[0] = (y < xminmax[0]) ? y : yminmax[0];
      yminmax[1] = (y > yminmax[1]) ? y : yminmax[1];
      zminmax[0] = (z < zminmax[0]) ? z : zminmax[0];
      zminmax[1] = (z > zminmax[1]) ? z : zminmax[1];
    }
  }
}