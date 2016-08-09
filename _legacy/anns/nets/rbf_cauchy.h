#ifndef ANNS_NETS_CAUCHY_H_INCLUDED
#define ANNS_NETS_CAUCHY_H_INCLUDED

// пая-й (RBF-C)
// CLASSIC | RADIAL | CAUCHY_FUNCTION

#include <math.h>
#include "../anns_types.h"

void rbf_cauchy_val(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
void rbf_cauchy_val_d1(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
void rbf_cauchy_val_d2(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);

void rbf_cauchy_valgrad(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
void rbf_cauchy_valgrad_d1(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
void rbf_cauchy_valgrad_d2(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);

void rbf_cauchy_loadfuncs(int t, void *anns_instance);

#endif // ANNS_NETS_CAUCHY_H_INCLUDED
