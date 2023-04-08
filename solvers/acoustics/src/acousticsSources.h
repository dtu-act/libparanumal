#pragma once

#include "acoustics.h"
#include <vector>
// #include "acousticsUtils.h"

dfloat gaussianSource(dfloat x, dfloat y, dfloat z, dfloat t, dfloat *sloc, dfloat sxyz, dfloat amplitude = AMPLITUDE);
// public template functions (not part of a class) need to be in the header file
template <typename T, typename U>
void gaussianSource(const std::vector<T>& x1d, const std::vector<T>& y1d, const std::vector<T>& z1d, 
                    dfloat *sloc, dfloat sxyz, std::vector<U> &pressures, dfloat amplitude = AMPLITUDE)
{
  pressures.resize(x1d.size());
  for (uint i=0; i < x1d.size(); ++i) { 
    pressures[i] = (U)gaussianSource((dfloat)x1d[i],(dfloat)y1d[i],(dfloat)z1d[i],0.0,sloc,sxyz,amplitude);
  }
}
#if INCLUDE_GRF
// Only double vector types are supported. Change to template types if needed.
void grfWindowed(std::vector<dfloat> x1d, std::vector<dfloat> y1d, std::vector<dfloat> z1d, 
    dfloat xminmax[2], dfloat yminmax[2], dfloat zminmax[2], 
    dfloat sigma_0, dfloat l_0, dfloat sigma0_window, std::vector<dfloat> &samples_out, dfloat amplitude = AMPLITUDE);
#endif