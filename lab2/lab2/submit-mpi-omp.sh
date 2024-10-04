#!/bin/bash

#SBATCH --job-name=submit-mpi-omp.sh
#SBATCH -D .
#SBATCH --output=submit-mpi-omp.sh.o%j
#SBATCH --error=submit-mpi-omp.sh.e%j
#SBATCH --nodes=3
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=8

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

PROGRAM=$1
size=$2

srun --mpi=pmi2 --cpu-bind=socket ./$PROGRAM $size

