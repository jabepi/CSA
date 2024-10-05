#!/bin/bash

#SBATCH --job-name=submit-omp-i.sh
#SBATCH -D .
#SBATCH --output=submit-omp-i.sh.o%j
#SBATCH --error=submit-omp-i.sh.e%j

export OMP_NUM_THREADS=4
export size=10

HOST=$(echo $HOSTNAME | cut -f 1 -d'.')

export LD_PRELOAD=${EXTRAE_HOME}/lib/libomptrace.so
../BackSubs_omp_trace.o $size
unset LD_PRELOAD

mpi2prv -f TRACE.mpits -o BackSubs_omp_trace-10-${HOST}.prv -e BackSubs_omp_trace.o -paraver
rm -rf  TRACE.mpits set-0 >& /dev/null
