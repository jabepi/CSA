#ifndef PCG_H
#define PCG_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "timing.h"

typedef void (*mul_fun_t)(int n, void* data, double* Ax, double* x);

extern float time_in_poisson;
extern float time_in_SSOR;

double dot(int n, const double* x, const double* y)
{
    double result = 0;
    for (int i = 0; i < n; ++i)
        result += x[i]*y[i];
    return result;
}

double pcg(int n,
           mul_fun_t Mfun, void* Mdata,
           mul_fun_t Afun, void* Adata,
           double* x,
           const double* b,
           int maxit,
           double rtol)
{
    double* r = malloc(n*sizeof(double));
    double* z = malloc(n*sizeof(double));
    double* q = malloc(n*sizeof(double));
    double* p = malloc(n*sizeof(double));

    double rho0     = 0;
    double rho      = 0;
    double rho_prev = 0;
    double rtol2 = rtol*rtol;
    int is_converged = 0;
    int step;

    tic(0);

    /* Form residual */
    Afun(n, Adata, r, x);

    for (int i = 0; i < n; ++i) r[i] = b[i]-r[i];

    for (step = 0; step < maxit && !is_converged; ++step) {
        Mfun(n, Mdata, z, r);
        rho_prev = rho;

        rho = dot(n, r, z);

        if (step == 0) {
            rho0 = rho;
            memcpy(p, z, n*sizeof(double));
        } else {
            double beta = rho/rho_prev;
            for (int i = 0; i < n; ++i) p[i] = z[i] + beta*p[i];
        }
        Afun(n, Adata, q, p);

        double alpha = rho/dot(n, p, q);
        for (int i = 0; i < n; ++i) x[i] += alpha*p[i];
        for (int i = 0; i < n; ++i) r[i] -= alpha*q[i];
        is_converged = (rho/rho0 < rtol2);
    }

    printf("%d steps, residual reduction %g (%s tol %g); time %g seconds.\n " \
           "Time in Poisson: %g seconds\nTime in SSOR: %g seconds\n", step,
            sqrt(rho/rho0), is_converged ? "<=" : ">", rtol, toc(0),
            time_in_poisson, time_in_SSOR);
    free(p);
    free(q);
    free(z);
    free(r);

    return rho/rho0;
}

#endif /* PCG_H */
