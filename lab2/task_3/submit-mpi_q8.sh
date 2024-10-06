#!/bin/bash

#SBATCH --job-name=submit-mpi.sh
#SBATCH -D .
#SBATCH --output=.ubmit-mpi.sh.o%j
#SBATCH --error=submit-mpi.sh.e%j
#SBATCH --nodes=1              # Use only 1 node
#SBATCH --cpus-per-task=1      # 1 CPU per task

cd ..
PROGRAM=$1
size=$2
make $1

# Loop to test with 1 to 20 cores
for ntasks in $(seq 1 20); do
    echo "Running with $ntasks tasks"
    sbatch -W --nodes=1 --ntasks=$ntasks --cpus-per-task=1 --wrap="srun --mpi=pmi2 --ntasks=$ntasks --cpu-bind=cores ./$PROGRAM $size"
    wait
done
