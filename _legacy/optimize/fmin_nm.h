#ifndef FMIN_NM_H_INCLUDED
#define FMIN_NM_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define FMIN_NM_CONV_FUNCTION 0
#define FMIN_NM_CONV_SIMPLEX_SIZE 1

typedef struct {
    double reflection;
    double expansion;
    double contraction;
    double reduction;

    int max_iterations;
    FILE *log;
    double simplex_min;
    double simplex_max;
    double *init;

    int convergence_criterion;
} fmin_nm_args_t;

typedef struct {
    int function_calls;
    int reflections;
    int expansions;
    int contractions;
    int reductions;

    double tol;
} fmin_nm_stats_t;

double *fmin_nm(double *x,
               int n,
               double (*val)(double *, int, void *),
               double tol,
               void *instance,

               fmin_nm_args_t *args,
               fmin_nm_stats_t *stats,
               double *work_array);


#endif // FMIN_NM_H_INCLUDED
