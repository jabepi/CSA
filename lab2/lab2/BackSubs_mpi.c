/* C Example of Backward Substitution */
#define _XOPEN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
   
/* header files for getting hostname and process id */ 
#include <unistd.h>
#include <sys/types.h>
   
/* Get time and resources */     
#include <sys/time.h>   
#include <sys/resource.h> 
   
//#define _DEBUG_  
//#include "constants.h"   
   
#if _EXTRAE_   
#include "extrae_user_events.h"    
// Extrae Constants    
#define  PROGRAM    1000   
#define  END    0  
#define  SERIAL 1  
#define  PARALLEL   2  
#else  
double getusec_() {    
    struct timeval time;   
    gettimeofday(&time, NULL);     
    return ((double)time.tv_sec * (double)1e6 + (double)time.tv_usec); 
}    
   
#define START_COUNT_TIME stamp = getusec_();   
#define STOP_COUNT_TIME(_m) stamp = getusec_() - stamp;\
    stamp = stamp/1e6;\
    printf ("Time: %s%0.6f\n",(_m), stamp);  
#endif

#ifdef SINGLE_PRECISION

typedef float fp_t;

#else

typedef double fp_t;

#endif 
  
int m; /* order of symmetric input matrix */
int n; /* number of right-hand sides */
int b; /* block size */
int reps; /* number of repetitions */
int check; /* check result */
fp_t alpha;
int lda;
int ldb;
char uplo; /* lower/upper triangular */
char side;
char trans;
char diag;

fp_t *A;
fp_t *B;
fp_t *Bchk;   
   
void GENMAT_IP(fp_t *A, int m, int n, int scale, int seed) 
{
	srand48(time(NULL)+seed);

	int j;
	for (j = 0; j < n; ++j ) {
		int i;
		for( i = 0; i < m; ++i ) {
			A[j*m+i] = (fp_t) drand48();
			if ( j == i )
				A[j*m+i] += scale;
		}
  	}
}

int trsm_setup(int check, int m, int n, int b, int lda, int ldb, fp_t **A, fp_t **B, fp_t **Bchk) 
{
	fp_t *lA = *A = malloc(lda*m*sizeof(fp_t));
	if ( lA == NULL )
		return 1;

	GENMAT_IP(lA, lda, m, 1,0);
    int j, i;
	for (j = 0; j < m; ++j ) 
		for( i = 0; i < j; ++i ) 
			lA[j*m+i] = (fp_t) 0.0;

	fp_t *lB = *B = malloc(ldb * n * sizeof(fp_t));
	if (lB == NULL)
		return 2;
	GENMAT_IP(lB, ldb, n, 1,1);

	fp_t *lBchk = *Bchk = malloc(ldb * n * sizeof(fp_t));
	if (lBchk == NULL)
		return 3;

	for ( i = 0; i < ldb*n; ++i )
		lBchk[i] = lB[i];

	return 0;
}

void trsm_shutdown(fp_t *A, fp_t *B, fp_t *X) 
{
	free(A);
	free(B);
	free(X);
}    
    

int main(int argc, char *argv[]) {
    int  myid, numprocs; 
#ifdef _DEBUG_ 
    char hostname[128];
#endif
    MPI_Status status;

    const char Usage[] = "Usage: BackSubs <size> (try 10000)\n";
    if (argc < 2) {
    fprintf(stderr, Usage);
    exit(1);
    }
    unsigned long long int size = atoll(argv[1]);
    
    side = 'L';
    uplo = 'U';

    m = size;
    n = size;
    b = 1;
    lda = size;
    ldb = 1;

    alpha = 1.0;
    trans = 'N';
    diag = 'N';
    reps = 1;
    check = 1;
    
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs); 
    MPI_Comm_rank(MPI_COMM_WORLD, &myid); 
    
    //Print the number of processes
    // if (myid == 0) {
    //     printf("Number of processes: %d\n", numprocs);
    // }
    // //Exit 
    // MPI_Finalize();
    // return 0;

    if (myid==0) { // Only Id 0 creates the matrix to share it
       if ( trsm_setup(check, m, n, b, lda, ldb, &A, &B, &Bchk) ) {
		fprintf(stderr, "err: allocating matrix\n");
		return 2;
       }
    } else {
       A = malloc(lda*m*sizeof(fp_t));
       B = malloc(ldb * n * sizeof(fp_t));
       if (A==NULL || B==NULL ) 
          return 1;
    }

#if _EXTRAE_
    Extrae_event (PROGRAM, SERIAL);
#endif

#if _EXTRAE_
    Extrae_event (PROGRAM, END);
#endif

    /* do computation -- using all available threads */
#if _EXTRAE_
    Extrae_event (PROGRAM, PARALLEL);
#endif

    MPI_Bcast(A, lda*m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(B, ldb * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#if _EXTRAE_
#else
    double stamp;
    if (myid == 0) START_COUNT_TIME;
#endif

    // Calculate the number of lines each process will handle
    int linesperprocess = size / numprocs;
    // Determine the starting line for this process
    int mylinesinit = myid * linesperprocess;
    // Determine the ending line for this process
    int mylinesend = (myid + 1) * linesperprocess;
    // Adjust the ending line if it exceeds the matrix size
    if ((mylinesend + linesperprocess) > size)
       mylinesend = size;

    int i, j;

    // Loop over each row in reverse order
    for (i = size - 1; i >= 0; i--) {
        // Determine which process is responsible for the current row
        int sender = i / linesperprocess;
        if (sender >= numprocs) sender = numprocs - 1;

        // If the current row is within the range of this process
        if (i >= mylinesinit && i < mylinesend) {
            // Perform the division step of backward substitution
            B[i] = B[i] / A[i * m + i];
            // Send the updated value to all previous processes
            for (int receiver = sender - 1; receiver >= 0; receiver--)
                MPI_Send(&B[i], 1, MPI_DOUBLE, receiver, 0, MPI_COMM_WORLD);
        }

        // If the current row is outside the range of this process
        if (i >= mylinesend)
            // Receive the updated value from the responsible process
            MPI_Recv(&B[i], 1, MPI_DOUBLE, sender, 0, MPI_COMM_WORLD, &status);

        // Update the remaining rows within the range of this process
        for (j = mylinesinit; j < i && j < mylinesend; j++) {
            B[j] = B[j] - B[i] * A[j * m + i];
        }
    }

 
#if _EXTRAE_
    Extrae_event (PROGRAM, END);
    Extrae_event (PROGRAM, SERIAL);
#else
    if (myid == 0) { STOP_COUNT_TIME("");}
#endif

#if _EXTRAE_
    Extrae_event (PROGRAM, END);
#endif

    if (myid == 0) {
        if (size<=1000) {
           if (size<=10) { // Simple printing code for debug purposes    
              for ( i = 0; i < lda; ++i ) { 
                for ( j = 0; j < m; ++j ) 
                  fprintf(stdout,"%lf ",A[i*m+j]);
              fprintf(stdout,"\n");
	      }
              fprintf(stdout,"\n");
              for ( i = 0; i < ldb*n; ++i ) 
                  fprintf(stdout,"%lf ",B[i]);
              fprintf(stdout,"\n");
 	      for ( i = 0; i < ldb*n; ++i ) 
                  fprintf(stdout,"%lf ",Bchk[i]);
              fprintf(stdout,"\n");
           }
        fp_t error=0;
        fp_t suma;
        for (i=size-1; i>=0; i--) {
            suma=0.0;
            for (j=i; j<=size-1; j++)
                suma+=B[j]*A[i*m+j];
            if (error<fabs(Bchk[i]-suma)) {
                error=fabs(Bchk[i]-suma);
		fprintf(stdout,"%d:%lf vs %lf -> %lf\n",i,suma,Bchk[i],error);
	    }
        }
        fprintf(stdout,"Error: %lf\n",error);
        }
    }

    if (myid == 0) {
        trsm_shutdown(A, B, Bchk);
    }else {
        free(A);
        free(B);
    }
    
#ifdef _DEBUG_ 
    gethostname(hostname, 126);
    printf( "Hello world from process %d of %d at hostname %s Processes lines %d to %d\n",myid, numprocs, hostname, mylinesinit, mylinesend );
#endif

    MPI_Finalize(); 

    return EXIT_SUCCESS;
}
