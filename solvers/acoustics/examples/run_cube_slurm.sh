#!/bin/bash

#SBATCH -t 05:00:00
#SBATCH --mem=32gb
#SBATCH -p a6000-gcondo
#SBATCH --gres=gpu:1
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o /users/nborrelj/data/nborrelj/logs/lipparanumal%j.out
#SBATCH -e /users/nborrelj/data/nborrelj/logs/lipparanumal%j.err
#SBATCH --job-name=lipparanumal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nikolas_borrel-jensen@brown.edu
#SBATCH -a 1-20%10

# OSCAR
module load gcc/10.2
module load mpi/openmpi_4.0.1_gcc
module load cuda/11.3.1
module load armadillo/11.2.4
module load hdf5/1.10.5

export OCCA_DIR=~/libparanumal/occa
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib

# Run solver
mpirun -np 1 ./acousticsMain examples/setups/setup_cube_1000Hz_p4_perf_refl
