#include "anns_cond.h"

anns_cond_t *ad_new(fbvp_t A, gbvp_t *A_g, int nvars, int m, int *eval, int nt) {
    anns_cond_t *ad = NULL;

    ad = malloc(sizeof *ad);

    assert(ad);

    ad->A = A;
    ad->A_g = A_g;
    ad->nvars = nvars;
    ad->m = m;
    ad->eval = eval;
    ad->nt = nt;
    ad->delta = 1;

    return ad;
}


void ad_delete(anns_cond_t *ad) {
    if (ad) {
        free(ad);
    }
}


void ad_fprintf(anns_cond_t *ad, FILE *fout) {
    int t;

    assert(ad);
    assert(fout);

    fprintf(fout, "nvars=%d m=%d nt=%d delta=%g eval={", ad->nvars, ad->m, ad->nt, ad->delta);
    fprintf(fout, "%d", ad->eval[0]);
    for (t = 1; t < ad->nt; t++) {
        fprintf(fout, ",%d", ad->eval[t]);
    }
    fprintf(fout, "}\n");

    fprintf(fout, "A=<function %p>\n", ad->A);
    fprintf(fout, "A_g={<function %p>", ad->A_g[0]);
    for (t = 1; t < ad->nt; t++) {
        fprintf(fout, ",<function %p>", ad->A_g[t]);
    }
    fprintf(fout, "}\n");
}
