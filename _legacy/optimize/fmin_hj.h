#ifndef AURO_OPTIMIZE_FMIN_HOOKE_JEEVES_H_INCLUDED
#define AURO_OPTIMIZE_FMIN_HOOKE_JEEVES_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

typedef struct {
    long int max_iterations;
    double max_epsilon;

    double lambda;
    double alpha;
    double *delta;
    int ndelta;
    double default_delta;
} fmin_hj_args_t;

fmin_hj_args_t *new_fmin_hj_args(int ndelta, double default_delta);
void delete_fmin_hj_args(fmin_hj_args_t *args);

double *fmin_hj(double (*func)(double *, int, void *),
                double *result,
                int n,
                double *x0,
                void *instance,
                void *method_instance,
                void (*callback)(double *, int, void *, int, double, void*),
                void *callback_instance);

#endif // AURO_OPTIMIZE_FMIN_HOOKE_JEEVES_H_INCLUDED
