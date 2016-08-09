#ifndef ANNS_UTILS_H_INCLUDED
#define ANNS_UTILS_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#ifndef uniform
#define uniform(a,b) (((b)-(a))*(rand()/(RAND_MAX + 1.0))+(a))
#endif

double **d2alloc(int rows, int cols);
double **d2allocs(int rows, int *cols);
double ***d3alloc(int rows, int cols, int depth);
void d2free(double **d2, int rows);
void d3free(double ***d3, int rows, int cols);

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

#define ANNS_CND(name) static double name(double *x, double *u, double **u_x, double **u_xx, double ***u_xy, double *v, void *instance)
#define ANNS_CND_G(name) static double name(double *x, double *u, double **u_x, double **u_xx, double ***u_xy, double u_a, double *u_ax, double *u_axx, double **u_axy, double *v, void *instance)

#define ANNS_FUN(name) static double name(double *x, int n, void *instance)

#define ANNS_FUNPTR(name) double (* name )(double *, int, void *)

#endif // ANNS_UTILS_H_INCLUDED
