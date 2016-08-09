#ifndef ANNS_COND_H_INCLUDED
#define ANNS_COND_H_INCLUDED

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#include "anns_types.h"

anns_cond_t *ad_new(fbvp_t A, gbvp_t *A_g, int nvars, int m, int *eval, int nt);
void ad_delete(anns_cond_t *ad);

void ad_fprintf(anns_cond_t *ad, FILE *fout);
#define ad_printf(ad) (ad_fprintf(ad,stdout))

#endif // ANNS_COND_H_INCLUDED
