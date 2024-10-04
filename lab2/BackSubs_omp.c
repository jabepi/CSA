/* C Example of Backward Substitution */
#define _XOPEN_SOURCE
#include <stdio.h>
#include <stdlib.h>                                                                          
#include <time.h>
#include <math.h>
#include <omp.h>

/* header files for getting hostname and process id */                                                                                                                             
#include <unistd.h>
#include <sys/types.h>

/* Get time and resources */

#include <sys/time.h>
#include <sys/resource.h>

#if _EXTRAE_                                                                                                                                                                       
#include "extrae_user_events.h"                                                                                                                                                    
// Extrae Constants                                                                                                                                                                
#define  PROGRAM    1000                                                                                                                                                           
#define  END        0                                                                                                                                                              
#define  SERIAL     1                                                                                                                                                              
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
                        printf ("%s%0.6f\n",(_m), stamp);                                                                                                                          
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


    #pragma omp parallel
    {
        #pragma omp master
        {
            printf("Running with %d threads\n", omp_get_num_threads());
        }
    }

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

#if _EXTRAE_
    Extrae_event (PROGRAM, SERIAL);
#else
    double stamp;
    START_COUNT_TIME;
#endif

    if ( trsm_setup(check, m, n, b, lda, ldb, &A, &B, &Bchk) ) {
		fprintf(stderr, "err: allocating matrix\n");
		return 2;
    }
    
#if _EXTRAE_
    Extrae_event (PROGRAM, END);
#endif

    /* do computation -- using all available threads */
#if _EXTRAE_
    Extrae_event (PROGRAM, PARALLEL);
#endif


    int i,j;
    for (i=size-1; i>=0; i--) {
        B[i] = B[i] / A[i * m + i];
        #pragma omp parallel for private(j)
        for (j = 0; j < i; j++) {
           B[j] = B[j] - B[i] * A[j * m + i];
        }
    }
    
#if _EXTRAE_
    Extrae_event (PROGRAM, END);
    Extrae_event (PROGRAM, SERIAL);
#endif

   if (size<=1000) {
        if (size<=10) { // Simple printing code for debug purposes    
            for ( i = 0; i < lda*m; ++i ) 
                fprintf(stdout,"%lf ",A[i]);
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
            if (error<fabs(Bchk[i]-suma))
                error=fabs(Bchk[i]-suma);
        }
    fprintf(stdout,"Error: %lf\n",error);
    }

#if _EXTRAE_
    Extrae_event (PROGRAM, END);
#else
    STOP_COUNT_TIME("");
#endif       
 
    trsm_shutdown(A, B, Bchk);
    
    return EXIT_SUCCESS;
}
