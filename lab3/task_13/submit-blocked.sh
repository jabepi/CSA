#!/bin/bash

#SBATCH --job-name=submit-seq.sbatch
#SBATCH -D .

export size=4096

make -C "/scratch/nas/1/sca1011/CSA/lab3" clean
make -C "/scratch/nas/1/sca1011/CSA/lab3" cholesky_blocked_deps_mkl
make -C "/scratch/nas/1/sca1011/CSA/lab3" cholesky_blocked_deps
make -C "/scratch/nas/1/sca1011/CSA/lab3" cholesky_blocked_deps_optimized

for threads in 4 8 16 20 32
do
    export OMP_NUM_THREADS=$threads
    export MKL_NUM_THREADS=$threads
    printf "Using %d threads\n" $threads
    ../cholesky_blocked_deps_mkl.o $size 0
done

printf "\n\n"

for threads in 4 8 16 20 32
do
    export OMP_NUM_THREADS=$threads
    printf "Using %d threads\n" $threads
    ../cholesky_blocked_deps.o $size 0
done

printf "\n\n"

for threads in 4 8 16 20 32
do
    export OMP_NUM_THREADS=$threads
    printf "Using %d threads\n" $threads
    ../cholesky_blocked_deps_optimized.o $size 0
done