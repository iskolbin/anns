#ifndef FMIN_CALLBACK_H_INCLUDED
#define FMIN_CALLBACK_H_INCLUDED

#include <stdio.h>

typedef struct {
    double *x;
    int n;
    void *instance;

    int iteration;
    double f;
    double gnorm;
    double gnorm2;

    FILE *fout;
    double (*mse) (double *, int, void *);
} fmin_callback_args_t;

typedef void (* fmin_callback_t) (fmin_callback_args_t *args);

void fmin_callback_fprintf(fmin_callback_args_t *args);
void fmin_callback_fprintf_mse(fmin_callback_args_t *args);

#endif // FMIN_CALLBACK_H_INCLUDED
