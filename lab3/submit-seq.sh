
export PROG=cholesky
make $PROG

for size in 1024 2048 4096; do
    echo "Running with size $size"
    sbatch -W --ntasks=1 --wrap=" srun ./$PROG $size"
    wait
done

