#!/bin/bash

#SBATCH --job-name=submit-seq.sbatch
#SBATCH -D .

export size=3200

printf "cholesky_blocked\n"
/usr/bin/time ../cholesky_blocked $size false # MKL one

printf "cholesky_blockseq\n"
/usr/bin/time ../cholesky_blockseq $size false # Manual one
