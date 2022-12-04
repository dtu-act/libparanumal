#!/bin/bash
#BSUB -q gpuv100
#BSUB -n 1
#BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -R "span[hosts=1]"
#BSUB -W 01:00
#BSUB -R "rusage[mem=15GB]"
#BSUB -J "logs/acoustics_example"

### -- Notify me by email when execution begins --
#BSUB -B
### -- Notify me by email when execution ends   --
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
#BSUB -o "logs/acoustics_example-%J.out"
#BSUB -e "logs/acoustics_example-%J.err"

module load mpi/4.1.2-gcc-10.3.0-binutils-2.36.1
module load cuda/11.6
module load hdf5/1.12.1-gcc-10.3.0
module load armadillo/11.2.4

export OCCA_DIR=~/libparanumal/occa
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib

# Run solver
mpirun -n 1 ./acousticsMain examples/setups/setup_cube_1000Hz_p4_perf_refl_test

