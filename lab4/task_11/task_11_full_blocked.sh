#!/bin/bash

#SBATCH --job-name=submit-seq.sbatch
#SBATCH -D .
#SBATCH --output=full_block_result.out

make -C "/scratch/nas/1/sca1011/CSA/lab4" clean
make -C "/scratch/nas/1/sca1011/CSA/lab4" CFLAGS="-std=gnu99 -O3 -march=native -Wall -m64 -g -DCLOCK=CLOCK_REALTIME"

../cgp3d.x -n 200 -p ssor -M 2000

make -C "/scratch/nas/1/sca1011/CSA/lab4" clean
make -C "/scratch/nas/1/sca1011/CSA/lab4" CFLAGS="-DUSE_BLOCKED_SSOR -DUSE_BLOCKING -DBk=4 -DBj=4 -DBi=32 -std=gnu99 -O3 -march=native -Wall -m64 -g -DCLOCK=CLOCK_REALTIME"

../cgp3d.x -n 200 -p ssor -M 2000
