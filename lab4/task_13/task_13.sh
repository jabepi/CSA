#!/bin/bash

#SBATCH --job-name=submit-seq.sbatch
#SBATCH -D .
#SBATCH --output=execution_result.out

make -C "/scratch/nas/1/sca1011/CSA/lab4" clean -s
make -C "/scratch/nas/1/sca1011/CSA/lab4" CFLAGS="-std=gnu99 -O3 -march=native -Wall -m64 -g -DCLOCK=CLOCK_REALTIME" -s

for size in 50 100 150 200; do
  printf -- "---- Running with size: %d ----\n" "$size"
  for o in $(seq 0 10 200); do
    echo -n "$o "
    ../cgp3d.x -o "$o" -M 500 -p as -n "$size" |
      grep "residual reduction" |
      sed -n 's/\([0-9]*\) steps.*time \([0-9.]*\) seconds.*/Iterations: \1, Time: \2 seconds/p'
  done
done