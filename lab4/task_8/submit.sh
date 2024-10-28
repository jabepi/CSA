#!/bin/bash

#SBATCH --job-name=task7.sbatch
#SBATCH -D .
#SBATCH --output=task7-%j.out

cd ..
make -s clean
make -s
cd task_8

for n in 100 200 400; do
    echo "Executing CGP3D with n=$n" 
    sbatch -W --nodes=1 --ntasks=1 --cpus-per-task=40 --export=OMP_DYNAMIC=false --export=OMP_NUM_THREADS=40 --wrap="../cgp3d.x -n $n"
    wait
done
