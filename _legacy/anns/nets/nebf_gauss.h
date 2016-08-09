#ifndef ANNS_NETS_NEBF_GAUSS_H_INCLUDED
#define ANNS_NETS_NEBF_GAUSS_H_INCLUDED

// НЭБС-Г (NEBF-G)
// NORMALIZED | ELLIPTIC | GAUSSIAN

#include <math.h>
#include "../anns_types.h"

// Про дополнительную нормировку см. nrbf_gauss.h

void nebf_gauss_val(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
void nebf_gauss_val_d1(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
void nebf_gauss_val_d2(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);

void nebf_gauss_valgrad(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
void nebf_gauss_valgrad_d1(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
void nebf_gauss_valgrad_d2(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);

void nebf_gauss_loadfuncs(int t, void *anns_instance);

#endif // ANNS_NETS_NEBF_GAUSS_H_INCLUDED
