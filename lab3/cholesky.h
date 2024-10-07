/*
* Copyright (c) 2017, BSC (Barcelona Supercomputing Center)
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of the <organization> nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY BSC ''AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL <copyright holder> BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <sys/time.h>
#include <sys/times.h>
#include <math.h>
#if USE_MKL
# include <mkl_lapacke.h>
# include <mkl.h>
#else
// USE_OPENBLAS
# include <clapack.h>
# include <cblas.h>
#endif

const int ts = 32; // tile size

#if defined(USE_DOUBLE)
#  define type_t     double
#  define ELEM_T_STR "double"
#  define gemm       cblas_dgemm
#  define trsm       cblas_dtrsm
#  define trmm       cblas_dtrmm
#  define syrk       cblas_dsyrk
#  if USE_MKL
#    define potrf    dpotrf
#    define lacpy    dlacpy
#    define lange    dlange
#    define larnv    dlarnv
#  else
#    define potrf    LAPACK_dpotrf
#    define lacpy    LAPACK_dlacpy
#    define lange    LAPACK_dlange
#    define larnv    LAPACK_dlarnv
#  endif
#else
#  define type_t     float
#  define ELEM_T_STR "float"
#  define gemm       cblas_sgemm
#  define trsm       cblas_strsm
#  define trmm       cblas_strmm
#  define syrk       cblas_ssyrk
#  if USE_MKL
#    define potrf    spotrf
#    define lacpy    slacpy
#    define lange    slange
#    define larnv    slarnv
#  else
#    define potrf    LAPACK_spotrf
#    define lacpy    LAPACK_slacpy
#    define lange    LAPACK_slange
#    define larnv    LAPACK_slarnv
#  endif
#endif
#define CBLAS_MAT_ORDER   CblasColMajor
#define CBLAS_T           CblasTrans
#define CBLAS_NT          CblasNoTrans
#define CBLAS_LO          CblasLower
#define CBLAS_UP          CblasUpper
#define CBLAS_LF          CblasLeft
#define CBLAS_RI          CblasRight
#define CBLAS_NU          CblasNonUnit

float get_time()
{
   static double gtod_ref_time_sec = 0.0;

   struct timeval tv;
   gettimeofday(&tv, NULL);

   // If this is the first invocation of through dclock(), then initialize the
   // "reference time" global variable to the seconds field of the tv struct.
   if (gtod_ref_time_sec == 0.0)
      gtod_ref_time_sec = (double) tv.tv_sec;

   // Normalize the seconds field of the tv struct so that it is relative to the
   // "reference time" that was recorded during the first invocation of dclock().
   const double norm_sec = (double) tv.tv_sec - gtod_ref_time_sec;

   // Compute the number of seconds since the reference time.
   const double t = norm_sec + tv.tv_usec * 1.0e-6;

   return (float) t;
}

static type_t pow_di(type_t x, int n)
{
   type_t rv = 1.0;

   if (n < 0) {
      n = -n;
      x = 1.0 / x;
   }

   for (; n; n >>= 1, x *= x) {
      if (n & 1) rv *= x;
   }

   return rv;
}

static void gather_block(const int N, type_t *Alin, type_t *A)
{
   for (int i = 0; i < ts; i++) {
      for (int j = 0; j < ts; j++) {
         A[i*ts + j] = Alin[i*N + j];
      }
   }
}

static void scatter_block(const int N, type_t *A, type_t *Alin)
{
   for (int i = 0; i < ts; i++) {
      for (int j = 0; j < ts; j++) {
         Alin[i*N + j] = A[i*ts + j];
      }
   }
}

static void convert_to_blocks(const int DIM, const int N, type_t (*Alin)[N], type_t *A[DIM][DIM])
{
   for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < DIM; j++) {
         gather_block(N, &Alin[i*ts][j*ts], A[i][j]);
      }
   }
}

static void convert_to_linear(const int DIM, const int N, type_t *A[DIM][DIM], type_t (*Alin)[N])
{
   for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < DIM; j++) {
         scatter_block(N, A[i][j], (type_t *) &Alin[i*ts][j*ts]);
      }
   }
}
