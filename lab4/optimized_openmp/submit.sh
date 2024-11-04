#!/bin/bash

#SBATCH --job-name=submit-seq.sbatch
#SBATCH -D .
#SBATCH --output=execution_result.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --export=OMP_DYNAMIC=false
#SBATCH --export=OMP_NUM_THREADS=40

make -C "/scratch/nas/1/sca1011/CSA/lab4/optimized_openmp" clean -s
make -C "/scratch/nas/1/sca1011/CSA/lab4/optimized_openmp" CFLAGS="-std=gnu99 -O3 -march=native -Wall -m64 -g -DCLOCK=CLOCK_REALTIME -DUSE_BLOCKING" -s

printf -- "---- 100 ----\n"
./cgp3d.x -n 100 -M 2000

printf -- "---- 200 ----\n"
./cgp3d.x -n 200 -M 2000

printf -- "---- 300 ----\n"
./cgp3d.x -n 300 -M 2000

printf -- "------------- Sequential -------------\n"

printf -- "---- 100 ----\n"
../cgp3d.x -n 100 -M 2000

printf -- "---- 200 ----\n"
../cgp3d.x -n 200 -M 2000

printf -- "---- 300 ----\n"
../cgp3d.x -n 300 -M 2000