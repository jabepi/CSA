#ifndef TIMING_H
#define TIMING_H

#include <time.h>

#define NWATCHES 8

#ifdef CLOCK
static struct timespec watches[NWATCHES];
#else
static clock_t watches[NWATCHES];
#endif

static void tic(int watch)
{
#ifdef CLOCK
    clock_gettime(CLOCK, watches + watch);
#else
    watches[watch] = clock();
#endif
}

static double toc(int watch)
{
    double elapsed;
#ifdef CLOCK
    struct timespec now;
    clock_gettime(CLOCK, &now);
    elapsed = now.tv_nsec - (double)watches[watch].tv_nsec;
    elapsed *= 1.0E-9L;
    elapsed += now.tv_sec - (double)watches[watch].tv_sec;
#else
    clock_t now = clock();
    elapsed = (double)(now - watches[watch]) / CLOCKS_PER_SEC;
#endif
    return elapsed;
}

#endif /* TIMING_H */
