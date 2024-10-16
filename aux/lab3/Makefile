CC       = gcc
SEQ      = -DUSE_MANUAL
OPENMP   = -fopenmp
#MKLROOT  = /scratch/nas/1/Soft/intel_parallel_studio_xe.2018/mkl
#MKLROOT = /scratch/nas/1/Soft/intel_parallel_studio_xe/compilers_and_libraries_2020.2.254/linux/mkl
MKLROOT = /scratch/nas/1/Soft/intel_parallel_studio_xe/mkl
CINCL   = -I. -I$(MKLROOT)/include
CLIBS   = -L. -L$(MKLROOT)/lib/intel64
RPATH   = -Wl,-rpath=$(MKLROOT)/lib/intel64
CFLAGS   = -Wall -Wimplicit-function-declaration -O3 -std=c99 -march=native -DUSE_MKL -DUSE_DOUBLE -lmkl_core -lmkl_sequential -lmkl_rt -lm $(CINCL) $(CLIBS) $(RPATH)

TARGETS = cholesky cholesky_seq cholesky_omp cholesky_blocked cholesky_blockseq cholesky_blockseq_omp cholesky_blocked_omp cholesky_blockseq_deps cholesky_blocked_deps

all: $(TARGETS)

cholesky: cholesky.c
	$(CC) $< $(CFLAGS) -o $@

cholesky_seq: cholesky.c
	$(CC) $< $(CFLAGS) $(SEQ) -o $@

cholesky_omp: cholesky_omp.c
	$(CC) $(CFLAGS) $(SEQ) $(OPENMP) $< -o $@

cholesky_blocked: cholesky_blocked.c
	$(CC) $(CFLAGS) $(OPENMP) $< -o $@

cholesky_blockseq: cholesky_blocked.c
	$(CC) $(CFLAGS) $(SEQ) $(OPENMP) $< -o $@

cholesky_blockseq_omp: cholesky_blocked_omp.c
	$(CC) $(CFLAGS) $(SEQ) $(OPENMP) $< -o $@

cholesky_blocked_omp: cholesky_blocked_omp.c
	$(CC) $(CFLAGS) $(OPENMP) $< -o $@

cholesky_blockseq_deps: cholesky_blocked_deps.c
	$(CC) $(CFLAGS) $(SEQ) $(OPENMP) $< -o $@

cholesky_blocked_deps: cholesky_blocked_deps.c
	$(CC) $(CFLAGS) $(OPENMP) $< -o $@

clean:
	rm -rf $(TARGETS) *.o *.pe* *.po* *.o* *.e*

ultraclean:
	rm -rf $(TARGETS) *.o *.pe* *.po* *.o* *.e* *.txt *.ps *.prv *.pcf *.row

