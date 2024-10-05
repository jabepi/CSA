#!/bin/bash

#SBATCH --job-name=submit-omp.sh
#SBATCH -D .

/usr/bin/time ../BackSubsOmpTime.o $1
