#!/bin/bash

for i in {1..20}
do
  export OMP_NUM_THREADS=$i
  sbatch --output=timing-omp.sh.o%j_${OMP_NUM_THREADS} timing_omp.sh 12000 $i
  sleep 5
done