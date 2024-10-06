#!/bin/bash

export OMP_NUM_THREADS=4
sbatch --output=timing_${OMP_NUM_THREADS}_4000 run_back_sub.sh 4000
sleep 5
sbatch --output=timing_BLAS_4000 run_back_sub_BLAS.sh 4000
