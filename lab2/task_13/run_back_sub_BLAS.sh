#!/bin/bash

#SBATCH --job-name=submit-omp.sh
#SBATCH -D .

/usr/bin/time ../BackSubsBLAS.o $1
