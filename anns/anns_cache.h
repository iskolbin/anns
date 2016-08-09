#ifndef ANNS_CACHE_H_INCLUDED
#define ANNS_CACHE_H_INCLUDED

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "anns_const.h"
#include "anns_utils.h"
#include "anns_types.h"

anns_cache_t *ac_new(int A, int dim, int nconds, int nnets, int nmem, anns_cond_t **conds, anns_net_t **nets);
void ac_push(anns_cache_t *cache, double *x);
void ac_put(anns_cache_t *cache, double *x, int i);
double *ac_get(anns_cache_t *ac, int i);
void ac_delete(anns_cache_t *ac);
void ac_genpoints(anns_cache_t *ac, double bounds[]);
void ac_fprintfpoints(FILE *fout, anns_cache_t *ac, int mode);
#define ac_printfpoints(ac,mode) ac_fprintfpoints(stdout,(ac),mode)
void ac_eval_v(anns_cache_t *ac, int c, int iv, ANNS_FUNPTR(f), void *instance);
void ac_load_v(anns_cache_t *ac, int c, int iv, double *darr);
void ac_set_v(anns_cache_t *ac, int c, int iv, double val);

void ac_v_fout(FILE *fout, anns_cache_t *ac);
void ac_p_fout(FILE *fout, anns_cache_t *ac);

void ac_switchpoints(anns_cache_t *ac);

#endif // ANNS_CACHE_H_INCLUDED
