#!/bin/bash

#SBATCH --job-name=submit-mpi.sh
#SBATCH -D .
#SBATCH --output=submit-mpi.sh.o%j
#SBATCH --error=submit-mpi.sh.e%j
#SBATCH --nodes=1-3             # Use 1 to 3 nodes
#SBATCH --cpus-per-task=1       # 1 CPU per task

cd ..
PROGRAM=$1
size=$2
make $1

# Loop to test with 1 to 3 nodes, 2 and 4 MPI processes per node

for nodes in $(seq 1 3); do
    for procs_per_node in 2 4; do
        total_procs=$((nodes * procs_per_node))
        
        if [ $total_procs -le $((3 * 4)) ]; then  # Check if total processes exceed 12 (3 nodes * 4 procs)
            echo "Running with $nodes nodes, $procs_per_node MPI processes per node (Total: $total_procs)"
            sbatch -W --nodes=$nodes --ntasks-per-node=$procs_per_node --cpus-per-task=1 --wrap="srun --mpi=pmi2 --ntasks=$total_procs --cpu-bind=cores ./$PROGRAM $size"
            wait
        fi
    done
done
