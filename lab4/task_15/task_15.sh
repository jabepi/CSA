#!/bin/bash

#SBATCH --job-name=submit-seq.sbatch
#SBATCH -D .
#SBATCH --output=execution_result.out

make -C "/scratch/nas/1/sca1011/CSA/lab4" clean -s
make -C "/scratch/nas/1/sca1011/CSA/lab4" CFLAGS="-std=gnu99 -O3 -march=native -Wall -m64 -g -DCLOCK=CLOCK_REALTIME -DUSE_RHS0" -s

printf -- "---- Running with RHS0 ----\n"

for size in 50 100 150 200; do
  printf -- "---- Running with size: %d ----\n" "$size"
    echo -n "ID: "
    ../cgp3d.x -M 500 -p id -n "$size" |
      grep "residual reduction" |
      sed -n 's/\([0-9]*\) steps.*time \([0-9.]*\) seconds.*/Iterations: \1, Time: \2 seconds/p'

    echo -n "SSOR: "
    ../cgp3d.x -w 1.90 -M 500 -p ssor -n "$size" |
      grep "residual reduction" |
      sed -n 's/\([0-9]*\) steps.*time \([0-9.]*\) seconds.*/Iterations: \1, Time: \2 seconds/p'

    echo -n "SCHWARZ : "
    ../cgp3d.x -o 2 -M 500 -p as -n "$size" |
      grep "residual reduction" |
      sed -n 's/\([0-9]*\) steps.*time \([0-9.]*\) seconds.*/Iterations: \1, Time: \2 seconds/p'
done

make -C "/scratch/nas/1/sca1011/CSA/lab4" clean -s
make -C "/scratch/nas/1/sca1011/CSA/lab4" CFLAGS="-std=gnu99 -O3 -march=native -Wall -m64 -g -DCLOCK=CLOCK_REALTIME" -s

printf -- "---- Running with RHS1 ----\n"

for size in 50 100 150 200; do
  printf -- "---- Running with size: %d ----\n" "$size"
    echo -n "ID: "
    ../cgp3d.x -M 500 -p id -n "$size" |
      grep "residual reduction" |
      sed -n 's/\([0-9]*\) steps.*time \([0-9.]*\) seconds.*/Iterations: \1, Time: \2 seconds/p'

    echo -n "SSOR: "
    ../cgp3d.x -w 1.90 -M 500 -p ssor -n "$size" |
      grep "residual reduction" |
      sed -n 's/\([0-9]*\) steps.*time \([0-9.]*\) seconds.*/Iterations: \1, Time: \2 seconds/p'

    echo -n "SCHWARZ : "
    ../cgp3d.x -o 2 -M 500 -p as -n "$size" |
      grep "residual reduction" |
      sed -n 's/\([0-9]*\) steps.*time \([0-9.]*\) seconds.*/Iterations: \1, Time: \2 seconds/p'
done
