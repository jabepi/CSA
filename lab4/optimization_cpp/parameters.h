#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <string>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>

enum class PreconditionerType {
    IdentityPreconditioner,
    SSORPreconditioner,
    SchwarzPreconditioner
};

struct SolverParameters {
    int meshSize = 100;
    int maximumIterations = 200;
    double relativeTolerance = 1e-6;
    PreconditionerType preconditionerType = PreconditionerType::IdentityPreconditioner;
    double ssorRelaxationParameter = 1.9;
    int overlapSize = 10;
};

inline void print_usage()
{
    fprintf(stderr,
            "cgp3d\n"
            "\t-h: print this message\n"
            "\t-n: mesh size (100)\n"
            "\t-M: maximum iteration count (200)\n"
            "\t-r: relative residual tolerance (1e-6)\n"
            "\t-p: preconditioner: as, ssor, or id (id)\n"
            "\t-w: SSOR relaxation parameter (1.9)\n"
            "\t-o: Schwarz method overlap (10)\n");
}

inline void set_default_parameters(SolverParameters& params)
{
    params.meshSize = 100;
    params.maximumIterations = 200;
    params.relativeTolerance = 1e-6;
    params.preconditionerType = PreconditionerType::IdentityPreconditioner;
    params.ssorRelaxationParameter = 1.9;
    params.overlapSize = 10;
}

inline int get_params(int argc, char** argv, SolverParameters& params)
{
    const char* optstring = "hn:M:r:p:w:o:";
    int option;

    set_default_parameters(params);

    while ((option = getopt(argc, argv, optstring)) != -1) {
        switch (option) {
            case 'h':
                print_usage();
                return -1;
            case 'p':
                if (strcmp(optarg, "id") == 0) {
                    params.preconditionerType = PreconditionerType::IdentityPreconditioner;
                } else if (strcmp(optarg, "ssor") == 0) {
                    params.preconditionerType = PreconditionerType::SSORPreconditioner;
                } else if (strcmp(optarg, "as") == 0) {
                    params.preconditionerType = PreconditionerType::SchwarzPreconditioner;
                } else {
                    fprintf(stderr, "Unknown preconditioner type: %s\n", optarg);
                    return -1;
                }
                break;
            case 'n':
                params.meshSize = std::atoi(optarg);
                break;
            case 'M':
                params.maximumIterations = std::atoi(optarg);
                break;
            case 'r':
                params.relativeTolerance = std::atof(optarg);
                break;
            case 'w':
                params.ssorRelaxationParameter = std::atof(optarg);
                break;
            case 'o':
                params.overlapSize = std::atoi(optarg);
                break;
            default:
                fprintf(stderr, "Unknown option\n");
                return -1;
        }
    }
    return 0;
}

#endif // PARAMS_HPP
