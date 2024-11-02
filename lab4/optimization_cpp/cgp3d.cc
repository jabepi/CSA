#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <vector>

#include "timing.h"
#include "pcg.h"
#include "params.h"

float time_in_poisson = 0;
float time_in_SSOR = 0;

#ifdef USE_NO_BRANCH_SSOR
void mul_poisson3d(int N, void* data, double* Ax, double* x)
{
    tic(1);
    #define X(i,j,k) (x[((k)*n+(j))*n+(i)])
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])

    int n = *(int*) data;
    int inv_h2 = (n - 1) * (n - 1);

    // Handle interior points where all neighbors are within bounds
    for (int k = 1; k < n - 1; ++k) {
        for (int j = 1; j < n - 1; ++j) {
            for (int i = 1; i < n - 1; ++i) {
                AX(i, j, k) = (6 * X(i, j, k)
                                - X(i - 1, j, k) - X(i + 1, j, k)
                                - X(i, j - 1, k) - X(i, j + 1, k)
                                - X(i, j, k - 1) - X(i, j, k + 1)) * inv_h2;
            }
        }
    }

    // Handle boundary faces k = 0 and k = n - 1
    for (int k = 0; k < n; k += n - 1) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                double xn = (i > 0)     ? X(i - 1, j, k) : 0;
                double xs = (i < n - 1) ? X(i + 1, j, k) : 0;
                double xe = (j > 0)     ? X(i, j - 1, k) : 0;
                double xw = (j < n - 1) ? X(i, j + 1, k) : 0;
                double xu = (k > 0)     ? X(i, j, k - 1) : 0;
                double xd = (k < n - 1) ? X(i, j, k + 1) : 0;
                AX(i, j, k) = (6 * X(i, j, k) - xn - xs - xe - xw - xu - xd) * inv_h2;
            }
        }
    }

    // Handle boundary faces j = 0 and j = n - 1
    for (int k = 1; k < n - 1; ++k) {
        for (int j = 0; j < n; j += n - 1) {
            for (int i = 0; i < n; ++i) {
                double xn = (i > 0)     ? X(i - 1, j, k) : 0;
                double xs = (i < n - 1) ? X(i + 1, j, k) : 0;
                double xu = X(i, j, k - 1);
                double xd = X(i, j, k + 1);
                double xe = (j > 0)     ? X(i, j - 1, k) : 0;
                double xw = (j < n - 1) ? X(i, j + 1, k) : 0;
                AX(i, j, k) = (6 * X(i, j, k) - xn - xs - xe - xw - xu - xd) * inv_h2;
            }
        }
    }

    // Handle boundary faces i = 0 and i = n - 1
    for (int k = 1; k < n - 1; ++k) {
        for (int j = 1; j < n - 1; ++j) {
            for (int i = 0; i < n; i += n - 1) {
                double xn = (i > 0)     ? X(i - 1, j, k) : 0;
                double xs = (i < n - 1) ? X(i + 1, j, k) : 0;
                double xe = X(i, j - 1, k);
                double xw = X(i, j + 1, k);
                double xu = X(i, j, k - 1);
                double xd = X(i, j, k + 1);
                AX(i, j, k) = (6 * X(i, j, k) - xn - xs - xe - xw - xu - xd) * inv_h2;
            }
        }
    }

    // Handle edges where two indices are at the boundary
    // Edges along k-direction
    for (int k = 1; k < n - 1; ++k) {
        for (int ij = 0; ij < 4; ++ij) {
            int i = (ij & 1) * (n - 1); // i = 0 or n - 1
            int j = ((ij >> 1) & 1) * (n - 1); // j = 0 or n - 1
            double xn = (i > 0)     ? X(i - 1, j, k) : 0;
            double xs = (i < n - 1) ? X(i + 1, j, k) : 0;
            double xe = (j > 0)     ? X(i, j - 1, k) : 0;
            double xw = (j < n - 1) ? X(i, j + 1, k) : 0;
            double xu = X(i, j, k - 1);
            double xd = X(i, j, k + 1);
            AX(i, j, k) = (6 * X(i, j, k) - xn - xs - xe - xw - xu - xd) * inv_h2;
        }
    }

    // Edges along j-direction
    for (int k = 0; k < n; k += n - 1) {
        for (int i = 1; i < n - 1; ++i) {
            for (int j = 0; j < n; j += n - 1) {
                double xn = X(i - 1, j, k);
                double xs = X(i + 1, j, k);
                double xe = (j > 0)     ? X(i, j - 1, k) : 0;
                double xw = (j < n - 1) ? X(i, j + 1, k) : 0;
                double xu = (k > 0)     ? X(i, j, k - 1) : 0;
                double xd = (k < n - 1) ? X(i, j, k + 1) : 0;
                AX(i, j, k) = (6 * X(i, j, k) - xn - xs - xe - xw - xu - xd) * inv_h2;
            }
        }
    }

    // Edges along i-direction
    for (int k = 0; k < n; k += n - 1) {
        for (int j = 1; j < n - 1; ++j) {
            for (int i = 0; i < n; i += n - 1) {
                double xn = (i > 0)     ? X(i - 1, j, k) : 0;
                double xs = (i < n - 1) ? X(i + 1, j, k) : 0;
                double xe = X(i, j - 1, k);
                double xw = X(i, j + 1, k);
                double xu = (k > 0)     ? X(i, j, k - 1) : 0;
                double xd = (k < n - 1) ? X(i, j, k + 1) : 0;
                AX(i, j, k) = (6 * X(i, j, k) - xn - xs - xe - xw - xu - xd) * inv_h2;
            }
        }
    }

    // Handle corner points where all three indices are at the boundary
    for (int k = 0; k < n; k += n - 1) {
        for (int j = 0; j < n; j += n - 1) {
            for (int i = 0; i < n; i += n - 1) {
                double xn = (i > 0)     ? X(i - 1, j, k) : 0;
                double xs = (i < n - 1) ? X(i + 1, j, k) : 0;
                double xe = (j > 0)     ? X(i, j - 1, k) : 0;
                double xw = (j < n - 1) ? X(i, j + 1, k) : 0;
                double xu = (k > 0)     ? X(i, j, k - 1) : 0;
                double xd = (k < n - 1) ? X(i, j, k + 1) : 0;
                AX(i, j, k) = (6 * X(i, j, k) - xn - xs - xe - xw - xu - xd) * inv_h2;
            }
        }
    }

    #undef AX
    #undef X
    time_in_poisson += toc(1);
}


#elif defined(USE_LOOP_ORDER)
void mul_poisson3d(int N, void* data, double* Ax, double* x)
{
    tic(1);
    #define X(i,j,k) (x[((k)*n+(j))*n+(i)])
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])

    int n = *(int*) data;
    int inv_h2 = (n - 1) * (n - 1);
#if defined(USE_JKI)
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
#elif defined(USE_KIJ)
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
#elif defined(USE_IKJ)
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            for (int j = 0; j < n; ++j) {
#elif defined(USE_IJK)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
#elif defined(USE_JIK)
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < n; ++k) {
#else
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
#endif
                double xx = X(i,j,k);
                double xn = (i > 0)   ? X(i-1,j,k) : 0;
                double xs = (i < n-1) ? X(i+1,j,k) : 0;
                double xe = (j > 0)   ? X(i,j-1,k) : 0;
                double xw = (j < n-1) ? X(i,j+1,k) : 0;
                double xu = (k > 0)   ? X(i,j,k-1) : 0;
                double xd = (k < n-1) ? X(i,j,k+1) : 0;
                AX(i,j,k) = (6*xx - xn - xs - xe - xw - xu - xd)*inv_h2;
            }
        }
    }

    #undef AX
    #undef X
    time_in_poisson += toc(1);
}
#elif defined(USE_PARTIAL_PREDICATION)
void mul_poisson3d(int N, void* data, double* Ax, double* x)
{
    tic(1);
    #define X(i,j,k) (x[((k)*n+(j))*n+(i)])
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])
    int n = *(int*) data;
    int inv_h2 = (n-1) * (n-1);
    for (int k = 0; k < n; ++k) {
        int k_gt_0 = k > 0;
        int k_lt_n1 = k < n - 1;
        int ku = k - k_gt_0;
        int kd = k + k_lt_n1;
        for (int j = 0; j < n; ++j) {
            int j_gt_0 = j > 0;
            int j_lt_n1 = j < n - 1;
            int je = j - j_gt_0;
            int jw = j + j_lt_n1;
            for (int i = 0; i < n; ++i) {
                int i_gt_0 = i > 0;
                int i_lt_n1 = i < n - 1;
                int in = i - i_gt_0;
                int is = i + i_lt_n1;
                double xx = X(i, j, k);
                double xn = X(in, j, k) * i_gt_0;
                double xs = X(is, j, k) * i_lt_n1;
                double xe = X(i, je, k) * j_gt_0;
                double xw = X(i, jw, k) * j_lt_n1;
                double xu = X(i, j, ku) * k_gt_0;
                double xd = X(i, j, kd) * k_lt_n1;
                AX(i, j, k) = (6 * xx - xn - xs - xe - xw - xu - xd) * inv_h2;
            }
        }
    }
    #undef AX
    #undef X
    time_in_poisson += toc(1);
}
#elif defined(USE_BLOCKING)
void mul_poisson3d(int N, void* data, double* Ax, double* x)
{
    tic(1);
    #define X(i,j,k) (x[((k)*n+(j))*n+(i)])
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])

    int n = *(int*) data;
    int inv_h2 = (n - 1) * (n - 1);

    // Set default block sizes if not defined
    #ifndef Bk
    #define Bk 8
    #endif

    #ifndef Bj
    #define Bj 8
    #endif

    #ifndef Bi
    #define Bi 32
    #endif

    for (int kk = 0; kk < n; kk += Bk) {
        for (int jj = 0; jj < n; jj += Bj) {
            for (int ii = 0; ii < n; ii += Bi) {
                int k_end = (kk + Bk < n) ? kk + Bk : n;
                int j_end = (jj + Bj < n) ? jj + Bj : n;
                int i_end = (ii + Bi < n) ? ii + Bi : n;

                for (int k = kk; k < k_end; ++k) {
                    for (int j = jj; j < j_end; ++j) {
                        for (int i = ii; i < i_end; ++i) {
                            double xx = X(i, j, k);
                            double xn = (i > 0)     ? X(i - 1, j, k) : 0;
                            double xs = (i < n - 1) ? X(i + 1, j, k) : 0;
                            double xe = (j > 0)     ? X(i, j - 1, k) : 0;
                            double xw = (j < n - 1) ? X(i, j + 1, k) : 0;
                            double xu = (k > 0)     ? X(i, j, k - 1) : 0;
                            double xd = (k < n - 1) ? X(i, j, k + 1) : 0;
                            AX(i, j, k) = (6 * xx - xn - xs - xe - xw - xu - xd) * inv_h2;
                        }
                    }
                }
            }
        }
    }

    #undef AX
    #undef X

    // Optionally undefine block sizes to avoid conflicts elsewhere
    #undef Bk
    #undef Bj
    #undef Bi

    time_in_poisson += toc(1);
}
#else
void mul_poisson3d(int N, void* data, double* Ax, double* x)
{
    tic(1);
    #define X(i,j,k) (x[((k)*n+(j))*n+(i)])
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])

    int n = *(int*) data;
    int inv_h2 = (n-1) * (n-1);
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                double xx = X(i,j,k);
                double xn = (i > 0)   ? X(i-1,j,k) : 0;
                double xs = (i < n-1) ? X(i+1,j,k) : 0;
                double xe = (j > 0)   ? X(i,j-1,k) : 0;
                double xw = (j < n-1) ? X(i,j+1,k) : 0;
                double xu = (k > 0)   ? X(i,j,k-1) : 0;
                double xd = (k < n-1) ? X(i,j,k+1) : 0;
                AX(i,j,k) = (6*xx - xn - xs - xe - xw - xu - xd)*inv_h2;
            }
        }
    }

    #undef AX
    #undef X
    time_in_poisson += toc(1);
}
#endif

void pc_identity(int n, void* data, std::vector<double>& Ax, const std::vector<double>& x)
{
    Ax.assign(x.begin(), x.end());
}

typedef struct pc_ssor_p3d_t {
    int n;          /* Number of points in one direction on mesh */
    double omega;   /* SSOR relaxation parameter */
} pc_ssor_p3d_t;

#ifdef USE_BLOCKED_SSOR
#define BLOCK_SIZE 8
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

void ssor_forward_sweep(int n, int i1, int i2, int j1, int j2, int k1, int k2,
                                double* Ax, double w)
{
    tic(2);
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])
    for (int kk = k1; kk < k2; kk += BLOCK_SIZE) {
        int k_max = MIN(kk + BLOCK_SIZE, k2);
        for (int jj = j1; jj < j2; jj += BLOCK_SIZE) {
            int j_max = MIN(jj + BLOCK_SIZE, j2);
            for (int ii = i1; ii < i2; ii += BLOCK_SIZE) {
                int i_max = MIN(ii + BLOCK_SIZE, i2);
                for (int k = kk; k < k_max; ++k) {
                    for (int j = jj; j < j_max; ++j) {
                        for (int i = ii; i < i_max; ++i) {
                            double xx = AX(i,j,k);
                            double xn = (i > 0)   ? AX(i-1,j,k) : 0.0;
                            double xe = (j > 0)   ? AX(i,j-1,k) : 0.0;
                            double xu = (k > 0)   ? AX(i,j,k-1) : 0.0;
                            AX(i,j,k) = (xx + xn + xe + xu) / 6.0 * w;
                        }
                    }
                }
            }
        }
    }
    #undef AX
    time_in_SSOR += toc(2);
}

void ssor_backward_sweep(int n, int i1, int i2, int j1, int j2, int k1, int k2,
                                 double* Ax, double w)
{
    tic(2);
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])
    for (int kk = k2 - 1; kk >= k1; kk -= BLOCK_SIZE) {
        int k_min = MAX(kk - BLOCK_SIZE + 1, k1);
        for (int jj = j2 - 1; jj >= j1; jj -= BLOCK_SIZE) {
            int j_min = MAX(jj - BLOCK_SIZE + 1, j1);
            for (int ii = i2 - 1; ii >= i1; ii -= BLOCK_SIZE) {
                int i_min = MAX(ii - BLOCK_SIZE + 1, i1);
                for (int k = kk; k >= k_min; --k) {
                    for (int j = jj; j >= j_min; --j) {
                        for (int i = ii; i >= i_min; --i) {
                            double xx = AX(i,j,k);
                            double xs = (i < n-1) ? AX(i+1,j,k) : 0.0;
                            double xw = (j < n-1) ? AX(i,j+1,k) : 0.0;
                            double xd = (k < n-1) ? AX(i,j,k+1) : 0.0;
                            AX(i,j,k) = (xx + xs + xw + xd) / 6.0 * w;
                        }
                    }
                }
            }
        }
    }
    #undef AX
    time_in_SSOR += toc(2);
}
#else
void ssor_forward_sweep(int n, int i1, int i2, int j1, int j2, int k1, int k2,
                        double* Ax, double w)
{
    tic(2);
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])
    for (int k = k1; k < k2; ++k) {
        for (int j = j1; j < j2; ++j) {
            for (int i = i1; i < i2; ++i) {
                double xx = AX(i,j,k);
                double xn = (i > 0)   ? AX(i-1,j,k) : 0;
                double xe = (j > 0)   ? AX(i,j-1,k) : 0;
                double xu = (k > 0)   ? AX(i,j,k-1) : 0;
                AX(i,j,k) = (xx+xn+xe+xu)/6*w;
            }
        }
    }
    #undef AX
    time_in_SSOR += toc(2);
}

void ssor_backward_sweep(int n, int i1, int i2, int j1, int j2, int k1, int k2,
                         double* Ax, double w)
{
    tic(2);
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])
    for (int k = k2-1; k >= k1; --k) {
        for (int j = j2-1; j >= j1; --j) {
            for (int i = i2-1; i >= i1; --i) {
                double xx = AX(i,j,k);
                double xs = (i < n-1) ? AX(i+1,j,k) : 0;
                double xw = (j < n-1) ? AX(i,j+1,k) : 0;
                double xd = (k < n-1) ? AX(i,j,k+1) : 0;
                AX(i,j,k) = (xx+xs+xw+xd)/6*w;
            }
        }
    }
    #undef AX
    time_in_SSOR += toc(2);
}
#endif

void ssor_diag_sweep(int n, int i1, int i2, int j1, int j2, int k1, int k2,
                     double* Ax, double w)
{
    tic(2);
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])
    for (int k = k1; k < k2; ++k)
        for (int j = j1; j < j2; ++j)
            for (int i = i1; i < i2; ++i)
                AX(i,j,k) *= (6*(2-w)/w);
    #undef AX
    time_in_SSOR += toc(2);
}

void pc_ssor_poisson3d(int N, void* data, std::vector<double>& Ax, const std::vector<double>& x)
{
    pc_ssor_p3d_t* ssor_data = (pc_ssor_p3d_t*) data;
    int n = ssor_data->n;
    double w = ssor_data->omega;

    std::copy(x.begin(), x.end(), Ax.begin());
    ssor_forward_sweep(n, 0, n, 0, n, 0, n, Ax.data(), w);
    ssor_diag_sweep(n, 0, n, 0, n, 0, n, Ax.data(), w);
    ssor_backward_sweep(n, 0, n, 0, n, 0, n, Ax.data(), w);
}

typedef struct pc_schwarz_p3d_t {
    int n;           /* Number of mesh points on a side */
    int overlap;     /* Number of points through the overlap region */
    double omega;    /* SSOR relaxation parameter */
    double* scratch; /* Scratch space used by the preconditioner */
} pc_schwarz_p3d_t;

void schwarz_get(int n, int i1, int i2, int j1, int j2, int k1, int k2,
                 double* x_local,
                 double* x)
{
    #define X(i,j,k) (x[((k)*n+(j))*n+(i)])
    #define XL(i,j,k) (x_local[((k)*n+(j))*n+(i)])

    for (int k = k1; k < k2; ++k)
        for (int j = j1; j < j2; ++j)
            for (int i = i1; i < i2; ++i)
                XL(i,j,k) = X(i,j,k);

    if (k1 > 0)
        for (int j = j1; j < j2; ++j)
            for (int i = i1; i < i2; ++i)
                XL(i,j,k1-1) = 0;
    if (j1 > 0)
        for (int k = k1; k < k2; ++k)
            for (int i = i1; i < i2; ++i)
                XL(i,j1-1,k) = 0;
    if (i1 > 0)
        for (int k = k1; k < k2; ++k)
            for (int j = j1; j < j2; ++j)
                XL(i1-1,j,k) = 0;

    if (k2 < n)
        for (int j = j1; j < j2; ++j)
            for (int i = i1; i < i2; ++i)
                XL(i,j,k2) = 0;
    if (j2 < n)
        for (int k = k1; k < k2; ++k)
            for (int i = i1; i < i2; ++i)
                XL(i,j2,k) = 0;
    if (i2 < n)
        for (int k = k1; k < k2; ++k)
            for (int j = j1; j < j2; ++j)
                XL(i2,j,k) = 0;

    #undef XL
    #undef X
}

void schwarz_add(int n, int i1, int i2, int j1, int j2, int k1, int k2,
                 double* Ax_local,
                 double* Ax)
{
    #define AX(i,j,k) (Ax[((k)*n+(j))*n+(i)])
    #define AXL(i,j,k) (Ax_local[((k)*n+(j))*n+(i)])
    for (int k = k1; k < k2; ++k)
        for (int j = j1; j < j2; ++j)
            for (int i = i1; i < i2; ++i)
                AX(i,j,k) += AXL(i,j,k);
    #undef AXL
    #undef AX
}

void pc_schwarz_poisson3d(int N, void* data, std::vector<double>& Ax, const std::vector<double>& x)
{
    pc_schwarz_p3d_t* ssor_data = (pc_schwarz_p3d_t*) data;
    double* scratch = ssor_data->scratch;
    int n = ssor_data->n;
    int o = ssor_data->overlap / 2;
    double w = ssor_data->omega;
    std::fill(Ax.begin(), Ax.end(), 0.0);

    int n1 = n / 2 + o;
    int n2 = n / 2 - o;

    schwarz_get(n, 0, n, 0, n, 0, n1, scratch, x.data());
    ssor_forward_sweep(n, 0, n, 0, n, 0, n1, scratch, w);
    ssor_diag_sweep(n, 0, n, 0, n, 0, n1, scratch, w);
    ssor_backward_sweep(n, 0, n, 0, n, 0, n1, scratch, w);
    schwarz_add(n, 0, n, 0, n, 0, n1, scratch, Ax.data());

    schwarz_get(n, 0, n, 0, n, n2, n, scratch, x.data());
    ssor_forward_sweep(n, 0, n, 0, n, n2, n, scratch, w);
    ssor_diag_sweep(n, 0, n, 0, n, n2, n, scratch, w);
    ssor_backward_sweep(n, 0, n, 0, n, n2, n, scratch, w);
    schwarz_add(n, 0, n, 0, n, n2, n, scratch, Ax.data());
}

void setup_rhs0(int n, double* b)
{
    int N = n*n*n;
    memset(b, 0, N*sizeof(double));
    b[0] = 1;
}

void setup_rhs1(int n, double* b)
{
    int N = n*n*n;
    memset(b, 0, N*sizeof(double));
    for (int i = 0; i < n; ++i) {
        double x = 1.0*(i+1)/(n+1);
        for (int j = 0; j < n; ++j) {
            double y = 1.0*(i+1)/(n+1);
            for (int k = 0; k < n; ++k) {
                double z = 1.0*(i+1)/(n+1);
                b[(k*n+j)*n+i] = x*(1-x) * y*(1-y) * z*(1-z);
            }
        }
    }
}

int main(int argc, char** argv)
{
    solve_param_t params;
    if (get_params(argc, argv, &params))
        return -1;

    int n = params.n;
    int N = n*n*n;

    std::vector<double> b(N, 0.0);
    std::vector<double> x(N, 0.0);
    std::vector<double> r(N, 0.0);

    /* Set up right hand side */
#ifdef USE_RHS0
    setup_rhs0(n, b.data());
#else
    setup_rhs1(n, b.data());
#endif

    /* Solve via PCG */
    int maxit   = params.maxit;
    double rtol = params.rtol;

    if (params.ptype == PC_SCHWARZ) {
        std::vector<double> scratch(N, 0.0);
        pc_schwarz_p3d_t pcdata = {n, params.overlap, params.omega, scratch.data()};
        pcg(N, pc_schwarz_poisson3d, &pcdata, mul_poisson3d, &n, x, b, maxit, rtol);
    } else if (params.ptype == PC_SSOR) {
        pc_ssor_p3d_t ssor_data = {n, params.omega};
        pcg(N, pc_ssor_poisson3d, &ssor_data, mul_poisson3d, &n, x, b, maxit, rtol);
    } else {
        pcg(N, pc_identity, NULL, mul_poisson3d, &n, x, b, maxit, rtol);
    }

    /* Check answer */
    mul_poisson3d(N, &n, r.data(), x.data());
    double rnorm2 = 0;
    for (int i = 0; i < n; ++i) r[i] = b[i]-r[i];
    for (int i = 0; i < n; ++i) rnorm2 += r[i]*r[i];
    printf("rnorm = %g\n", sqrt(rnorm2));
}
