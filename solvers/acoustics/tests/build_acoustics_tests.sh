#!/bin/bash

# Load modules
module purge
module load mpi/4.0.3-gcc-8.4.0
module load cuda/10.0

# OCCA environment variables
export OCCA_DIR=~/libparanumal/occa
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib

# Build OCCA - OCCA must be build on a system with GPU access for GPU support (on HPC switch to 'voltash')
#cd ../../occa
#make clean
#make -j 4
#cd ../solvers/acoustics/

# Build project
#make clean
make tests -j 4