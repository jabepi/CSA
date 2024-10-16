#!/bin/bash

#SBATCH --job-name=submit-seq.sbatch
#SBATCH -D .

export size=4096

make -C "/scratch/nas/1/sca1011/CSA/lab3" clean
make -C "/scratch/nas/1/sca1011/CSA/lab3" cholesky_blocked_deps_mkl
make -C "/scratch/nas/1/sca1011/CSA/lab3" cholesky_blocked_deps

unset MKL_NUM_THREADS
../cholesky_blocked_deps_mkl.o $size 0

printf "\n\n\n\n\n\n\n"

unset OMP_NUM_THREADS
../cholesky_blocked_deps.o $size 0
