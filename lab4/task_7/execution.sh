#!/bin/bash

#SBATCH --job-name=submit-omp.sh
#SBATCH -D .
#SBATCH --output=submit-omp.sh.o%j
#SBATCH --error=submit-omp.sh.e%j

cd .. 
make clean
make
cd task_6

for n in 100 200 400; do
    ../cgp3d.x -n $n
done
