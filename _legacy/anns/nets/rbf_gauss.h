#ifndef ANNS_NETS_RBF_GAUSS_H_INCLUDED
#define ANNS_NETS_RBF_GAUSS_H_INCLUDED

// тву-з (RBF-G)
// CLASSIC | RADIAL | GAUSSIAN

#include <math.h>
#include "../anns_types.h"

void rbf_gauss_val(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void rbf_gauss_val_d1(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void rbf_gauss_val_d2(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);

void rbf_gauss_valgrad(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void rbf_gauss_valgrad_d1(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void rbf_gauss_valgrad_d2(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);

void rbf_gauss_val_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void rbf_gauss_val_d1_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void rbf_gauss_val_d2_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);

void rbf_gauss_valgrad_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void rbf_gauss_valgrad_d1_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void rbf_gauss_valgrad_d2_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);


//void rbf_gauss_loadfuncs(int t, void *anns_instance);
void rbf_gauss_loadfuncs(anns_net_t *net);
#endif // ANNS_NETS_RBF_GAUSS_H_INCLUDED
