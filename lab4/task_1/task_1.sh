#!/bin/bash

#SBATCH --job-name=submit-seq.sbatch
#SBATCH -D .
#SBATCH --output=execution_result.out

make -C "/scratch/nas/1/sca1011/CSA/lab4" clean
make -C "/scratch/nas/1/sca1011/CSA/lab4"

../cgp3d.x -h

../cgp3d.x

../cgp3d.x -M 500

../cgp3d.x -M 199