#!/bin/bash

#SBATCH --job-name=submit-omp.sh
#SBATCH -D .
#SBATCH --output=submit-omp.sh.o%j
#SBATCH --error=submit-omp.sh.e%j

/usr/bin/time -o time-$1 ../BackSubsOmpTime.o $1
