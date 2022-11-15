#!/bin/bash

#SBATCH -t 00:05:00
#SBATCH --mem=8gb
#SBATCH -p a6000-gcondo
#SBATCH --gres=gpu:1
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o /users/nborrelj/data/nborrelj/logs/libparanumal%j.out
#SBATCH -e /users/nborrelj/data/nborrelj/logs/libparanumal%j.err
#SBATCH --job-name=libp_gaussians
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nikolas_borrel-jensen@brown.edu
#SBATCH -a 0-4

# OSCAR
module load gcc/10.2
module load mpi/openmpi_4.0.1_gcc
module load cuda/11.3.1
module load armadillo/11.2.4
module load hdf5/1.10.5

export OCCA_DIR=~/libparanumal/occa
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib

SAMPLE_LIST=($(<simulationSetups/deeponet/joblist.list))
SAMPLE=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}

echo This is task $SLURM_ARRAY_TASK_ID, which will run $SAMPLE

# Run solver
mpirun -np 1 ./acousticsMain simulationSetups/deeponet/$SAMPLE