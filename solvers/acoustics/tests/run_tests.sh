#!/bin/bash
#BSUB -q gpuv100
#BSUB -n 2
#BSUB -gpu "num=2:mode=exclusive_process"
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=15GB]"
#BSUB -W 01:00
#BSUB -J "logs/testing"

### -- Notify me by email when execution begins --
#BSUB -B
### -- Notify me by email when execution ends   --
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
#BSUB -o "logs/testing-%J.out"
#BSUB -e "logs/testing-%J.err"

# DTU HPC
module load mpi/3.1.4-gcc-9.2.0
module load cuda/11.3

# OSCAR
# module load gcc/10.2
# module load mpi/openmpi_3.1.6_gcc
# module load cuda/11.3.1
# module load armadillo/9.200.4

export OCCA_DIR=~/libparanumal/occa
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib

# Run solver
mpirun -n 1 ./acousticsTestMain # [cylinder][freq_indep]