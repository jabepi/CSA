#!/bin/bash

#SBATCH --job-name=task7.sbatch
#SBATCH -D .
#SBATCH --output=task7-%j.out

cd ..
make -s clean
make -s
cd task_7

for n in 100 200 400; do
    for p in id ssor as; do

        echo "Executing CGP3D with n=$n and p=$p"
        ../cgp3d.x -n $n -p $p
    done
done