CC       = gcc
MPICC    = mpicc.mpich
# Needed if user events are introduced in source file
MKLROOT = /scratch/nas/1/Soft/intel_parallel_studio_xe/mkl
IINCL    = -I${EXTRAE_HOME}/include
ILIBS    = -L$(EXTRAE_HOME)/lib -lomptrace
CINCL   = -I. -I$(MKLROOT)/include
CLIBS   = -L. -L$(MKLROOT)/lib/intel64
RPATH   = -Wl,-rpath=$(MKLROOT)/lib/intel64
MKLBLASFLAGS   = -Wall -Wimplicit-function-declaration -O3 -std=c99 -march=native -DUSE_MKL -DUSE_DOUBLE -lmkl_core -lmkl_sequential -lmkl_rt -lm $(CINCL) $(CLIBS) $(RPATH)
CFLAGS   = -Wall -Wimplicit-function-declaration -O3 -std=c99 -march=native -fopenmp
BLASFLAGS   = -Wall -O3 -std=c99 -march=native -lblas -llapack -lm

TARGETS = BackSubs BackSubs_i BackSubs_omp BackSubs_omp_i BackSubs_mpi BackSubs_blas MKLBackSubs_blas

all: $(TARGETS)

BackSubs: BackSubs.c
	$(CC) $< $(CFLAGS) -o $@.o

BackSubsOmp: BackSubs_omp.c
	$(CC) $< $(CFLAGS) -o $@.o

BackSubsOmpTime: BackSubs_omp_time.c
	$(CC) $< $(CFLAGS) -o $@.o

BackSubsBLAS: BackSubs_blas.c
	$(CC) $< $(BLASFLAGS) -o $@.o

BackSubs_i: BackSubs.c
	$(CC) -D_EXTRAE_ $(CFLAGS) $(OPENMP) $(IINCL) $< -o $@ $(ILIBS)

BackSubs_omp: BackSubs_omp.c
	$(CC) $(CFLAGS) $(OPENMP) $< -o $@.o
	
BackSubs_omp_trace: BackSubs_omp.c
	$(CC) -D_EXTRAE_ $(CFLAGS) $(OPENMP) $(IINCL) $< -o $@.o $(ILIBS)

BackSubs_omp_i: BackSubs_omp.c
	$(CC) -D_EXTRAE_ $(CFLAGS) $(OPENMP) $(IINCL) $< -o $@ $(ILIBS)

BackSubs_mpi: BackSubs_mpi.c
	$(MPICC) $< $(CFLAGS) $(CINCL) -o $@ $(CLIBS)

BackSubs_mpi_q8: BackSubs_mpi_q8.c
	$(MPICC) $< $(CFLAGS) $(CINCL) -o $@ $(CLIBS)

BackSubs_blas: BackSubs_blas.c
	$(CC) $< $(BLASFLAGS) -o $@

MKLBackSubs_blas: BackSubs_blas.c
	$(CC) $< $(MKLBLASFLAGS) -o $@

clean:
	rm -rf $(TARGETS) *.o *.pe* *.po* *.o* *.e*

ultraclean:
	rm -rf $(TARGETS) *.o *.pe* *.po* *.o* *.e* *.txt *.ps *.prv *.pcf *.row *.sym

