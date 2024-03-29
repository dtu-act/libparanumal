#!/bin/bash

#BSUB -W 00:05
#BSUB -q gpuv100
#BSUB -n 4
#BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=1GB]"
#BSUB -J furn_train[1-1000]%5

### -- Notify me by email when execution begins --
##BSUB -B
### -- Notify me by email when execution ends   --
##BSUB -N
### -- Specify the output and error file. %J is the job-id --
#BSUB -o "/work3/nibor/data/logs/libparanumal_furn_train_%J_%I.out"
#BSUB -e "/work3/nibor/data/logs/libparanumal_furn_train_%J_%I.err"

module load mpi/4.1.2-gcc-10.3.0-binutils-2.36.1
module load cuda/11.6
module load hdf5/1.12.1-gcc-10.3.0
module load armadillo/11.2.4

export OCCA_DIR=~/libparanumal/occa
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib

src_index=`expr $LSB_JOBINDEX - 1`
echo "BASH: Simulation using source index $src_index from source position mesh."

# Run solver
mpirun -np 1 ./acousticsMain simulationSetups/deeponet/setup_train_val_srcpos_mesh/setup_furnished_train $src_index
