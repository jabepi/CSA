#!/bin/bash

#SBATCH --job-name=submit-omp.sh
#SBATCH -D .

export OMP_NUM_THREADS=$2
/usr/bin/time -o time-$1-$2 ../BackSubs_omp.o $1
