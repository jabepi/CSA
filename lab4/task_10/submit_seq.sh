#!/bin/bash

#SBATCH --job-name=task10_seq
#SBATCH -D .
#SBATCH --output=task10_seq_%j.out

cd ..
make -s clean
make -s
cd task_10

for n in 100 200 400; do
    echo "Executing CGP3D with n=$n" 
    ../cgp3d.x -n $n -p ssor
done