#!/bin/bash

#SBATCH --job-name=submit-mpi.sh
#SBATCH -D .
#SBATCH --output=submit-mpi.sh.o%j
#SBATCH --error=submit-mpi.sh.e%j
#SBATCH --nodes=1-3             # Use 1 to 3 nodes
#SBATCH --cpus-per-task=1       # 1 CPU per task

cd ../lab2/
PROGRAM=$1
size=$2
make $1

# Loop to test with 1 to 3 nodes, 1 MPI process per node

for nodes in $(seq 1 3); do
    echo "Running with $nodes nodes, 1 MPI process per node"
    sbatch -W --nodes=$nodes --ntasks-per-node=1 --cpus-per-task=1 --wrap="srun --mpi=pmi2 --ntasks=$nodes --cpu-bind=cores ./$PROGRAM $size"
    wait
done
