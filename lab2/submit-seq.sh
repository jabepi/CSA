#!/bin/bash

#SBATCH --job-name=submit-seq.sbatch
#SBATCH -D .
#SBATCH --output=submit-seq.sbatch.o%j
#SBATCH --error=submit-seq.sbatch.e%j

export PROG=BackSubs
make $PROG

export size=10000


/usr/bin/time -o ${PROG}_time.txt ./$PROG $size
