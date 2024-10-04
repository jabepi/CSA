#!/bin/bash

#SBATCH --job-name=submit-omp.sh
#SBATCH -D .
#SBATCH --output=timing-omp.sh.o%j
#SBATCH --error=timing-omp.sh.e%j

export OMP_NUM_THREADS=$2
/usr/bin/time -o time-$1-$2 ../BackSubs_omp.o $1
