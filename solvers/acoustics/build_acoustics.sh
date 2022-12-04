#!/bin/bash

# DTU HPC
module load mpi/4.1.2-gcc-10.3.0-binutils-2.36.1
module load cuda/11.6
module load hdf5/1.12.1-gcc-10.3.0
module load armadillo/11.2.4

# OCCA environment variables
export OCCA_DIR=~/libparanumal/occa
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib

# Build OCCA - OCCA must be build on a system with GPU access for GPU support (on HPC switch to 'voltash')
cd ../../occa
make clean
make -j 4
cd ../solvers/acoustics/

# Build project
make clean
make -j 4