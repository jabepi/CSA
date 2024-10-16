#!/bin/bash

#SBATCH --job-name=submit-mpi.sh
#SBATCH -D .
#SBATCH --output=submit-mpi.sh.o%j
#SBATCH --error=submit-mpi.sh.e%j
#SBATCH --nodes=1              # Use only 1 node
#SBATCH --cpus-per-task=1      # 1 CPU per task

make cholesky_blocked_omp