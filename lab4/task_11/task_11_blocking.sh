#!/bin/bash

#SBATCH --job-name=submit-seq.sbatch
#SBATCH -D .
#SBATCH --output=blocking_result.out

Bk_values=(4 8 16 32)
Bj_values=(4 8 16 32)
Bi_values=(4 8 16 32)

make -C "/scratch/nas/1/sca1011/CSA/lab4" clean -s
make -C "/scratch/nas/1/sca1011/CSA/lab4" CFLAGS="-std=gnu99 -O3 -march=native -Wall -m64 -g -DCLOCK=CLOCK_REALTIME" -s

echo "non blocked"
../cgp3d.x -n 100

for Bk in "${Bk_values[@]}"; do
  for Bj in "${Bj_values[@]}"; do
    for Bi in "${Bi_values[@]}"; do
      echo "Testing Bk=$Bk, Bj=$Bj, Bi=$Bi"

      # Clean previous build
      make -C "/scratch/nas/1/sca1011/CSA/lab4" clean -s

      # Compile with current block sizes
      make -C "/scratch/nas/1/sca1011/CSA/lab4" \
        CFLAGS="-DUSE_BLOCKING -DBk=$Bk -DBj=$Bj -DBi=$Bi -std=gnu99 -O3 -march=native -Wall -m64 -g -DCLOCK=CLOCK_REALTIME" -s

      ../cgp3d.x
    done
  done
done