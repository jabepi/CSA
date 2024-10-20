#!/bin/bash

#SBATCH --job-name=submit-seq.sbatch
#SBATCH -D .
#SBATCH --output=execution_result.out

make -C "/scratch/nas/1/sca1011/CSA/lab4" clean
make -C "/scratch/nas/1/sca1011/CSA/lab4" CFLAGS="-DUSE_RHS0"

../cgp3d.x

../cgp3d.x -p as

../cgp3d.x -p ssor

make -C "/scratch/nas/1/sca1011/CSA/lab4" clean
make -C "/scratch/nas/1/sca1011/CSA/lab4"

../cgp3d.x

../cgp3d.x -p as

../cgp3d.x -p ssor