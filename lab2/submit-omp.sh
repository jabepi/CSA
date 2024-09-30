#!/bin/bash

#SBATCH --job-name=submit-omp.sh
#SBATCH -D .
#SBATCH --output=submit-omp.sh.o%j
#SBATCH --error=submit-omp.sh.e%j


USAGE="\n USAGE: submit-omp.sh PROG size n_th\n
        PROG   -> omp program name\n
        size   -> size of the problem\n
	n_th   -> number of threads\n"

if (test $# -lt 3 || test $# -gt 3)
then
	echo -e $USAGE
        exit 0
fi

make $1

HOST=$(echo $HOSTNAME | cut -f 1 -d'.')

export OMP_NUM_THREADS=$3

/usr/bin/time -o time-$1-$3-${HOST}.txt ./$1 $2
