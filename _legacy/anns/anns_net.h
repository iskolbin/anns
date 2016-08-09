#ifndef ANNS_NET_H_INCLUDED
#define ANNS_NET_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include "anns_types.h"

anns_net_t *an_new(double *z, int nc, int nm, int nsize, int dim, int *map);
void an_delete(anns_net_t *an);

double an_eval(anns_net_t *an, double *point);
double *an_ninsert(anns_net_t *an, int from, int count);
double *an_nremove(anns_net_t *an, int from, int count);
void an_fprintf(anns_net_t *an, FILE *fout, char *valf, char *valsep, char *nsep, int aux);
double an_distance(anns_net_t *an, int i, int j);

double *an_cut(anns_net_t *an, double *buffer, int i, size_t n);
void an_paste(anns_net_t *an, double *buffer, int i, size_t n);
void an_overwrite(anns_net_t *an, double *buffer, int i, size_t n);

double an_min_distance(anns_net_t *an, int i, int *nearest);
double an_max_distance(anns_net_t *an, int i, int *farthest);
void an_swap(anns_net_t *an, int i, int j);
void an_setwidth_nearest(anns_net_t *an, double beta);

#define an_neuron(an,i) ((an)->z + (i)*(an)->nsize)
#define an_aux(an,T_,I_,zeroinsert_) {(an)->T = (T_); (an)->I = (I_); (an)->zeroinsert = (zeroinsert_);}
#define an_nappend(an,count) (an_ninsert((an), (an)->nc, (count)))
#define an_ntrim(an,count) (an_nremove((an), (an)->nc - (count), (count)) )
#define an_printf(an) (an_fprintf((an), stdout, "%g", " ", "\n", 1) )

void an_test(void);

#endif // ANNS_NET_H_INCLUDED
