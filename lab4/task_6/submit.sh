#!/bin/bash

#SBATCH --job-name=task6
#SBATCH -D .
#SBATCH --output=task6_%j.out

cd ..
make -s clean
make -s
cd task_6

for n in 100 200 400; do
    ../cgp3d.x -n $n
done