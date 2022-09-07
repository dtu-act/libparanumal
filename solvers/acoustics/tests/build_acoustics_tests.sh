#!/bin/bash

# DTU HPC
module load mpi/3.1.4-gcc-9.2.0
module load cuda/11.3

# OSCAR
# module load gcc/10.2
# module load mpi/openmpi_3.1.6_gcc
# module load cuda/11.3.1
# module load armadillo/9.200.4

# OCCA environment variables
export OCCA_DIR=~/libparanumal/occa
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib

# Build OCCA - OCCA must be build on a system with GPU access for GPU support (on HPC switch to 'voltash')
# cd ../../occa
# make clean
# make -j 4
# cd ../solvers/acoustics/

# Build project
# make clean
make tests -j 4