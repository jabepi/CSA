#!/bin/bash

#SBATCH --job-name=task9.sbatch
#SBATCH -D .
#SBATCH --output=task9-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --export=OMP_DYNAMIC=false
#SBATCH --export=OMP_NUM_THREADS=40


cd ..
make -s clean
make -s
cd task_9

for n in 100 200 400; do
    echo "Executing CGP3D with n=$n" 
    ../cgp3d.x -n $n -p as
done
