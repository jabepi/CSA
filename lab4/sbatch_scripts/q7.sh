#!/bin/bash

#SBATCH --job-name=q7.sh
#SBATCH -D .
#SBATCH --output=q7.sh.o%j
#SBATCH --error=q7.sh.e%j

cd .. 
#Recompile the code
make

PROGRAM=cgp3d.x
precondition=id


./$PROGRAM -p $precondition -n 100
