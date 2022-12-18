#!/bin/bash

#BSUB -W 00:10
#BSUB -q gpuv100
#BSUB -n 1
#BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=1GB]"
#BSUB -J Array_Gaussians[2310-3000]%5

### -- Notify me by email when execution begins --
##BSUB -B
### -- Notify me by email when execution ends   --
##BSUB -N
### -- Specify the output and error file. %J is the job-id --
#BSUB -o "/work3/nibor/data/logs/libparanumal_%J_%I.out"
#BSUB -e "/work3/nibor/data/logs/libparanumal_%J_%I.err"

module load mpi/4.1.2-gcc-10.3.0-binutils-2.36.1
module load cuda/11.6
module load hdf5/1.12.1-gcc-10.3.0
module load armadillo/11.2.4

export OCCA_DIR=~/libparanumal/occa
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib

src_index=`expr $LSB_JOBINDEX - 1`
echo BASH: Simulation using source index $src_index from source position mesh.

# Run solver
mpirun -np 1 ./acousticsMain simulationSetups/deeponet/setup_1000Hz_p4_rand_src $src_index
