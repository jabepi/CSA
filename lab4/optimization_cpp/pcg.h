#ifndef PCG_H
#define PCG_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "timing.h"
#include <vector>

typedef void (*mul_fun_t)(int n, void* data, double* Ax, double* x);
typedef void (*mul_fun_vec_t)(int n, void* data, std::vector<double>& Ax, const std::vector<double>& x);

extern float time_in_poisson;
extern float time_in_SSOR;

double dot(int n, const std::vector<double>& x, const std::vector<double>& y)
{
    double result = 0;
    for (int i = 0; i < n; ++i)
        result += x[i]*y[i];
    return result;
}

double pcg(int n,
           mul_fun_vec_t Mfun, void* Mdata,
           mul_fun_t Afun, void* Adata,
           std::vector<double>& x,
           const std::vector<double>& b,
           int maxit,
           double rtol)
{
    std::vector<double> r(n, 0.0);
    std::vector<double> z(n, 0.0);
    std::vector<double> q(n, 0.0);
    std::vector<double> p(n, 0.0);

    double rho0     = 0;
    double rho      = 0;
    double rho_prev = 0;
    double rtol2 = rtol*rtol;
    int is_converged = 0;
    int step;

    tic(0);

    /* Form residual */
    Afun(n, Adata, r.data(), x.data());

    for (int i = 0; i < n; ++i) r[i] = b[i]-r[i];

    for (step = 0; step < maxit && !is_converged; ++step) {
        Mfun(n, Mdata, z, r);
        rho_prev = rho;

        rho = dot(n, r, z);

        if (step == 0) {
            rho0 = rho;
            p = z;
        } else {
            double beta = rho/rho_prev;
            for (int i = 0; i < n; ++i) p[i] = z[i] + beta*p[i];
        }
        Afun(n, Adata, q.data(), p.data());

        double alpha = rho/dot(n, p, q);
        for (int i = 0; i < n; ++i) x[i] += alpha*p[i];
        for (int i = 0; i < n; ++i) r[i] -= alpha*q[i];
        is_converged = (rho/rho0 < rtol2);
    }

    printf("%d steps, residual reduction %g (%s tol %g); time %g seconds.\n " \
           "Time in Poisson: %g seconds\nTime in SSOR: %g seconds\n", step,
            sqrt(rho/rho0), is_converged ? "<=" : ">", rtol, toc(0),
            time_in_poisson, time_in_SSOR);

    return rho/rho0;
}

#endif /* PCG_H */
