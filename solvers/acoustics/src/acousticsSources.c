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
#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

#if INCLUDE_GRF
#include <armadillo>
using namespace arma;
#endif

dfloat gaussianSource(dfloat x, dfloat y, dfloat z, dfloat t, dfloat *sloc, dfloat sxyz, dfloat amplitude) {
  return amplitude*std::exp(-((x - sloc[0])*(x - sloc[0]) + (y - sloc[1])*(y - sloc[1]) + (z - sloc[2])*(z - sloc[2]))/(sxyz*sxyz));
}

void gaussianSource(const vector<dfloat> x1d, const vector<dfloat> y1d, const vector<dfloat> z1d, dfloat *sloc, dfloat sxyz, 
    vector<dfloat> &pressures, dfloat amplitude)
{
  pressures.resize(x1d.size());
  for (uint i=0; i < x1d.size(); ++i) { 
    pressures[i] = gaussianSource(x1d[i],y1d[i],z1d[i],0,sloc,sxyz,amplitude);
  }
}

#if INCLUDE_GRF
void applyWindowFunction(vector<dfloat> x1d, vector<dfloat> y1d, vector<dfloat> z1d, 
    dfloat xminmax[], dfloat yminmax[], dfloat zminmax[], dfloat sigma0_window, vector<dfloat> &samples_out, dfloat amplitude);
void grf(vector<dfloat> x1d_, vector<dfloat> y1d_, vector<dfloat> z1d_, dfloat sigma_0, dfloat l_0, vector<dfloat> &samples_out);

void grfWindowed(vector<dfloat> x1d, vector<dfloat> y1d, vector<dfloat> z1d, dfloat xminmax[2], dfloat yminmax[2], dfloat zminmax[2], 
    dfloat sigma_0, dfloat l_0, dfloat sigma0_window, vector<dfloat> &samples_out, dfloat amplitude) {

    grf(x1d, y1d, z1d, sigma_0, l_0, samples_out);
    applyWindowFunction(x1d, y1d, z1d, xminmax, yminmax, zminmax, sigma0_window, samples_out, amplitude);
}

void grf(vector<dfloat> x1d_, vector<dfloat> y1d_, vector<dfloat> z1d_, dfloat sigma_0, dfloat l_0, vector<dfloat> &samples_out) {
    arma_rng::set_seed_random();

    int num_samples = 1;
    vec x1d = vec(x1d_);
    vec y1d = vec(y1d_);
    vec z1d = vec(z1d_);

    int num_nodes = x1d.n_elem;

    mat distances_squared = mat(num_nodes, num_nodes);

    for (int i = 0; i < num_nodes; i++) {
        for (int j = 0; j < num_nodes; j++) {
            distances_squared(i,j)  = pow(x1d(i) - x1d(j), 2) + pow(y1d(i) - y1d(j), 2) + pow(z1d(i) - z1d(j), 2);
        }
    }

    mat covariance_matrix = (1.0/(sigma_0*sqrt(2*datum::pi)) * exp(- 0.5 / pow(l_0,2) * distances_squared));    
    vec mu = zeros<vec>(num_nodes);
    
    mat samples = mvnrnd(mu, covariance_matrix, num_samples);    

    samples_out.resize(samples.n_elem);

    // copy into output matrix
    for (int i=0; i<samples.n_elem; i++) {
        samples_out[i] = samples(i);
    }
}

// HELPER FUNCTIONS
void clamp(vector<dfloat> &v) {
    dfloat max = *max_element(v.begin(), v.end());
    dfloat min = *min_element(v.begin(), v.end());

    if (min < -1) {
        for (int i=0; i<v.size(); i++) {
            v[i] /= std::abs(min);
        }
    }
        
    if (max > 1) {
        for (int i=0;i<v.size();i++) {
            v[i] /= std::abs(max);
        }
    }
}

void clampPerSample(mat &samples) {
    samples.each_col(
        [](vec& a) {
            if (a.min() < -1)
                a /= abs(a).min();
            if (a.max() > 1)
                a /= abs(a).max();
        }
    );
}

void applyWindowFunction(vector<dfloat> x1d, vector<dfloat> y1d, vector<dfloat> z1d, 
    dfloat xminmax[], dfloat yminmax[], dfloat zminmax[], dfloat sigma0_window, vector<dfloat> &samples_out, dfloat amplitude) {

    dfloat offset = sigma0_window*3;
    dfloat x0 = xminmax[0] + offset;
    dfloat y0 = yminmax[0] + offset;
    dfloat z0 = zminmax[0] + offset;
    dfloat xhat0 = xminmax[1] - offset;
    dfloat yhat0 = yminmax[1] - offset;
    dfloat zhat0 = zminmax[1] - offset;

    auto gauss_fn = [sigma0_window](dfloat x, dfloat y, dfloat z, dfloat x0, dfloat y0, dfloat z0) {
        return std::exp(-((std::pow(x-x0,2) + std::pow(y-y0,2) + std::pow(z-z0,2))/std::pow(sigma0_window,2)));
    };

    auto mask_fun = [x0,y0,z0,xhat0,yhat0,zhat0,gauss_fn](dfloat x, dfloat y, dfloat z) {
        // 6 squares
        // x
        if (x <= x0 and (y >= y0 and y <= yhat0) and (z >= z0 and z <= zhat0))
            return gauss_fn(x,0,0,x0,0,0);           
        else if (x >= xhat0 and (y >= y0 and y <= yhat0) and (z >= z0 and z <= zhat0))
            return gauss_fn(x,0,0,xhat0,0,0);
        // y
        else if (y <= y0 and (x >= x0 and x <= xhat0) and (z >= z0 and z <= zhat0))
            return gauss_fn(0,y,0,0,y0,0);
        else if (y >= yhat0 and (x >= x0 and x <= xhat0) and (z >= z0 and z <= zhat0))
            return gauss_fn(0,y,0,0,yhat0,0);
        // z
        else if (z <= z0 and (x >= x0 and x <= xhat0) and (y >= y0 and y <= yhat0))
            return gauss_fn(0,0,z,0,0,z0);
        else if (z >= zhat0 and (x >= x0 and x <= xhat0) and (y >= y0 and y <= yhat0))
            return gauss_fn(0,0,z,0,0,zhat0);
        
        // 12 edges
        // along y-axis lower z-plane
        else if (x <= x0 and (y >= y0 and y <= yhat0) and z <= z0)
            return gauss_fn(x,0,z,x0,0,z0);
        else if (x >= xhat0 and (y >= y0 and y <= yhat0) and z <= z0)
            return gauss_fn(x,0,z,xhat0,0,z0);
        // along y-axis upper z-plane
        else if (x <= x0 and (y >= y0 and y <= yhat0) and z >= zhat0)
            return gauss_fn(x,0,z,x0,0,zhat0);
        else if (x >= xhat0 and (y >= y0 and y <= yhat0) and z >= zhat0)
            return gauss_fn(x,0,z,xhat0,0,zhat0);
        // along x-axis lower z-plane
        else if (y <= y0 and (x >= x0 and x <= xhat0) and z <= z0)
            return gauss_fn(0,y,z,0,y0,z0);
        else if (y >= yhat0 and (x >= x0 and x <= xhat0) and z <= z0)
            return gauss_fn(0,y,z,0,yhat0,z0);
        // along x-axis upper z-plane
        else if (y <= y0 and (x >= x0 and x <= xhat0) and z >= zhat0)
            return gauss_fn(0,y,z,0,y0,zhat0);
        else if (y >= yhat0 and (x >= x0 and x <= xhat0) and z >= zhat0)
            return gauss_fn(0,y,z,0,yhat0,zhat0);
        // along z-axis nearest y-plane
        else if (x <= x0 and (z >= z0 and z <= zhat0) and y <= y0)
            return gauss_fn(x,y,0,x0,y0,0);
        else if (x >= xhat0 and (z >= z0 and z <= xhat0) and y <= y0)
            return gauss_fn(x,y,0,xhat0,y0,0);
        // along z-axis furthest y-plane
        else if (x <= x0 and (z >= z0 and z <= zhat0) and y >= yhat0)
            return gauss_fn(x,y,0,x0,yhat0,0);
        else if (x >= xhat0 and (z >= z0 and z <= xhat0) and y >= yhat0)
            return gauss_fn(x,y,0,xhat0,yhat0,0);

        // 8 corners
        // lower z-plane
        else if (x <= x0 and y <= y0 and z <= z0)
            return gauss_fn(x,y,z,x0,y0,z0);
        else if (x >= xhat0 and y <= y0 and z <= z0)
            return gauss_fn(x,y,z,xhat0,y0,z0);
        else if (x <= x0 and y >= yhat0 and z <= z0)
            return gauss_fn(x,y,z,x0,yhat0,z0);
        else if (x >= xhat0 and y >= yhat0 and z <= z0)
            return gauss_fn(x,y,z,xhat0,yhat0,z0);
        // upper z-plane
        else if (x <= x0 and y <= y0 and z >= z0)
            return gauss_fn(x,y,z,x0,y0,zhat0);
        else if (x >= xhat0 and y <= y0 and z >= z0)
            return gauss_fn(x,y,z,xhat0,y0,zhat0);
        else if (x <= x0 and y >= yhat0 and z >= z0)
            return gauss_fn(x,y,z,x0,yhat0,zhat0);
        else if (x >= xhat0 and y >= yhat0 and z >= z0)
            return gauss_fn(x,y,z,xhat0,yhat0,zhat0);
        else
            return 1.0; // inside
    };
    
    int N = samples_out.size();
    std::vector<dfloat> mask = vector<dfloat>(N);

    for (int i=0; i<N; i++) {
        mask[i] = mask_fun(x1d[i],y1d[i], z1d[i]);
        samples_out[i] *= mask[i];
    }

    // clip all samples between -1 and 1
    clamp(samples_out);

    for (int i=0; i<N; i++) {
        samples_out[i] *= amplitude;
    }
}
#endif