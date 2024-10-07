#!/bin/bash

#SBATCH --job-name=submit-mpi.sh
#SBATCH -D .
#SBATCH --output=submit-mpi.sh.o%j
#SBATCH --error=submit-mpi.sh.e%j
#SBATCH --nodes=1-3            # Use 1 to 3 nodes
#SBATCH --cpus-per-task=20     # 20 CPUs per task

cd ..
PROGRAM=$1
size=$2
make $1

# Loop to test with 1 to 3 nodes, with 1 process per node and 20 threads per process
for nodes in $(seq 1 3); do
    echo "Running with $nodes nodes, 1 process per node, and 20 threads per process"
    sbatch -W --nodes=$nodes --ntasks=$nodes --cpus-per-task=20 --wrap="srun --mpi=pmi2 --ntasks=$nodes --cpus-per-task=20 --cpu-bind=cores ./$PROGRAM $size"
    wait
done
