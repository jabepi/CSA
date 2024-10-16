export PROG=cholesky_blocked_omp
make $PROG

# Define thread counts to test
THREADS_LIST="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 20 24 28 32"

# File to store all outputs
OUTPUT_FILE="results.txt"
echo "" > $OUTPUT_FILE  # Clear the output file

for threads in $THREADS_LIST; do
    export OMP_NUM_THREADS=$threads  # Set OpenMP threads
    echo "Running using $threads threads" | tee -a $OUTPUT_FILE
    JOB_ID=$(sbatch --ntasks=1 --nodes=1 --cpus-per-task=$threads -W --wrap="srun ./$PROG 4096 0")
    wait
    
    JOB_ID=$(echo $JOB_ID | awk '{print $4}')
    OUT_FILE="slurm-$JOB_ID.out"

    grep "Execution time (secs):" $OUT_FILE | awk '{print $4}' | tee -a $OUTPUT_FILE
done

