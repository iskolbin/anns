#ifndef ANNS_NETS_EBF_GAUSS_H_INCLUDED
#define ANNS_NETS_EBF_GAUSS_H_INCLUDED

// щая-ц (EBF-G)
// CLASSIC | ELLIPTIC | GAUSSIAN

#include <math.h>
#include "../anns_types.h"

void ebf_gauss_val(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
void ebf_gauss_val_d1(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
void ebf_gauss_val_d2(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);

void ebf_gauss_valgrad(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
void ebf_gauss_valgrad_d1(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
void ebf_gauss_valgrad_d2(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);

void ebf_gauss_loadfuncs(int t, void *anns_instance);

#endif // ANNS_NETS_EBF_GAUSS_H_INCLUDED
