#!/bin/bash

#SBATCH --job-name=task10.sbatch
#SBATCH -D .
#SBATCH --output=task10-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --export=OMP_DYNAMIC=false
#SBATCH --export=OMP_NUM_THREADS=20


cd ..
make -s clean
make -s
cd task_10

for n in 100 200 400; do
    echo "Executing CGP3D with n=$n" 
    ../cgp3d.x -n $n -p ssor
done
