#ifndef ACOUSTICS_UTILS_H
#define ACOUSTICS_UTILS_H

#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sstream>

std::string generateUUID(int length);
void extractUniquePoints(mesh_t *mesh, acoustics_t *acoustics, 
  std::vector<std::vector<uint>> &conn, std::vector<float> &x1d, 
  std::vector<float> &y1d, std::vector<float> &z1d, std::vector<float> &p1d);
void meshMinMax(mesh_t *mesh, dfloat xminmax[2], dfloat yminmax[2], dfloat zminmax[2]);
void pointReader3D(char *fileName, std::vector<dfloat> &VX, std::vector<dfloat> &VY, std::vector<dfloat> &VZ);

template <typename T>
std::string toString(const T a_value, const int n = 5)
{
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

#endif /* ACOUSTICS_UTILS_H */