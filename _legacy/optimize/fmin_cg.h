#ifndef FMIN_CG_H_INCLUDED
#define FMIN_CG_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "fmin_linesearch.h"

#define FMIN_CG_HS  0
#define FMIN_CG_FR  1
#define FMIN_CG_PRP 2
#define FMIN_CG_CD  3
#define FMIN_CG_LS  4
#define FMIN_CG_DY  5
#define FMIN_CG_N   6

#define FMIN_CG_HEURISTICS_NONE 0
#define FMIN_CG_HEURISTICS_PRP  1
#define FMIN_CG_HEURISTICS_P    2
#define FMIN_CG_HEURISTICS_N    3

typedef struct {
    int max_iterations;
    FILE *log;

    int cg_type;
    int heuristic_type;
    fmin_linesearch_t linesearch;
} fmin_cg_args_t;

typedef struct {
    int function_calls;
    int gradient_calls;

    double tol;
} fmin_cg_stats_t;

double *fmin_cg(double *x,
                int n,

                double (*val) (double *, int, void *),
                void (*grad) (double *, double *, int, void *),
                double (*valgrad) (double *, double *, int, void *),

                double tol,
                void *instance,

                fmin_cg_args_t *args,
                fmin_cg_stats_t *stats,
                double *work_array);

#endif // FMIN_CG_H_INCLUDED
