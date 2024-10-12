#!/bin/bash

#SBATCH --job-name=submit-seq.sbatch
#SBATCH -D .

export size=3200

printf "cholesky_blocked omp NON optimized\n"
/usr/bin/time ../cholesky_blockseq_omp $size 1


printf "cholesky_blocked omp optimized\n"
/usr/bin/time ../cholesky_blockseq_omp_optimized $size 1


printf "cholesky_blocked omp NON optimized\n"
/usr/bin/time ../cholesky_blockseq_omp $size 0


printf "cholesky_blocked omp optimized\n"
/usr/bin/time ../cholesky_blockseq_omp_optimized $size 0