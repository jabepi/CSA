#!/bin/bash

export OMP_NUM_THREADS=1
sbatch --output=timing_${OMP_NUM_THREADS}_2000 run_back_sub.sh 2000
sleep 5
export OMP_NUM_THREADS=4
sbatch --output=timing_${OMP_NUM_THREADS}_4000 run_back_sub.sh 4000
sleep 5
export OMP_NUM_THREADS=16
sbatch --output=timing_${OMP_NUM_THREADS}_8000 run_back_sub.sh 8000
sleep 10
export OMP_NUM_THREADS=64
sbatch --output=timing_${OMP_NUM_THREADS}_16000 run_back_sub.sh 16000
sleep 20
export OMP_NUM_THREADS=256
sbatch --output=timing_${OMP_NUM_THREADS}_32000 run_back_sub.sh 32000
sleep 20