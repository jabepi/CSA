#!/bin/bash

#SBATCH --job-name=submit-seq.sbatch
#SBATCH -D .
#SBATCH --output=partial_predication.out

make -C "/scratch/nas/1/sca1011/CSA/lab4" clean -s
make -C "/scratch/nas/1/sca1011/CSA/lab4" CFLAGS="-std=gnu99 -O3 -march=native -Wall -m64 -g -DCLOCK=CLOCK_REALTIME"

../cgp3d.x

make -C "/scratch/nas/1/sca1011/CSA/lab4" clean -s
make -C "/scratch/nas/1/sca1011/CSA/lab4" CFLAGS="-DUSE_PARTIAL_PREDICATION -std=gnu99 -O3 -march=native -Wall -m64 -g -DCLOCK=CLOCK_REALTIME" -s

../cgp3d.x

make -C "/scratch/nas/1/sca1011/CSA/lab4" clean -s
make -C "/scratch/nas/1/sca1011/CSA/lab4" CFLAGS="-std=gnu99 -O3 -march=native -Wall -m64 -g -DCLOCK=CLOCK_REALTIME"

../cgp3d.x -n 200 -M 2000

make -C "/scratch/nas/1/sca1011/CSA/lab4" clean -s
make -C "/scratch/nas/1/sca1011/CSA/lab4" CFLAGS="-DUSE_PARTIAL_PREDICATION -std=gnu99 -O3 -march=native -Wall -m64 -g -DCLOCK=CLOCK_REALTIME" -s

../cgp3d.x -n 200 -M 2000

make -C "/scratch/nas/1/sca1011/CSA/lab4" clean -s
make -C "/scratch/nas/1/sca1011/CSA/lab4" CFLAGS="-std=gnu99 -O3 -march=native -Wall -m64 -g -DCLOCK=CLOCK_REALTIME"

../cgp3d.x -n 300 -M 2000

make -C "/scratch/nas/1/sca1011/CSA/lab4" clean -s
make -C "/scratch/nas/1/sca1011/CSA/lab4" CFLAGS="-DUSE_PARTIAL_PREDICATION -std=gnu99 -O3 -march=native -Wall -m64 -g -DCLOCK=CLOCK_REALTIME" -s

../cgp3d.x -n 300 -M 2000
