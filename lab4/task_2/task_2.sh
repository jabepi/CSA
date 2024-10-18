#!/bin/bash

#SBATCH --job-name=submit-seq.sbatch
#SBATCH -D .
#SBATCH --output=execution_result.out

make -C "/scratch/nas/1/sca1011/CSA/lab4" clean
make -C "/scratch/nas/1/sca1011/CSA/lab4"

../cgp3d.x -n 200

../cgp3d.x -n 200 -M 500