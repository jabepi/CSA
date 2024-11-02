#include <vector>
#include <cstring>
#include <cmath>
#include <cstdio>

int main(int argc, char** argv)
{
    solve_param_t params;
    if (get_params(argc, argv, &params))
        return -1;

    int n = params.n;
    int N = n * n * n;

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
        pcg(N, pc_schwarz_poisson3d, &pcdata, mul_poisson3d, &n, x.data(), b.data(),
            maxit, rtol);
    } else if (params.ptype == PC_SSOR) {
        pc_ssor_p3d_t ssor_data = {n, params.omega};
        pcg(N, pc_ssor_poisson3d, &ssor_data, mul_poisson3d, &n, x.data(), b.data(),
            maxit, rtol);
    } else {
        pcg(N, pc_identity, nullptr, mul_poisson3d, &n, x.data(), b.data(), maxit, rtol);
    }

    /* Check answer */
    mul_poisson3d(N, &n, r.data(), x.data());
    double rnorm2 = 0.0;
    for (int i = 0; i < n; ++i) {
        r[i] = b[i] - r[i];
    }
    for (int i = 0; i < n; ++i) {
        rnorm2 += r[i] * r[i];
    }
    printf("rnorm = %g\n", std::sqrt(rnorm2));

    return 0;
}
