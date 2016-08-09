#ifndef ANNS_H_INCLUDED
#define ANNS_H_INCLUDED

#include <math.h>

#include "anns_instance.h"
#include "anns_const.h"
#include "anns_net.h"
#include "anns_cond.h"
#include "anns_loadfuncs.h"

double anns_val_f(double *x, int A, void *anns_instance);
void anns_grad_f(double *g, double *x, int A, void *anns_instance);
double anns_valgrad_f(double *g, double *x, int A, void *anns_instance);

double anns_val_f_check(double *x, int A, void *anns_instance);
double anns_evalsingle(int t, double *point, double *x, void *anns_instance);
void anns_csum(double *x, int A, void *anns_instance);
double anns_eval(double *point, double *x, void *anns_instance);

double anns_val_f_wonly(double *x, int A, void *instance);
void anns_grad_f_wonly(double *g, double *x, int A, void *anns_instance);
double anns_valgrad_f_wonly(double *g, double *x, int A, void *anns_instance);

double *anns_refine(anns_instance_t *ai);

#endif // ANNS_H_INCLUDED
