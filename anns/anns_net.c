#include "anns_net.h"
#include <math.h>

anns_net_t *an_new(double *z, int nc, int nm, int nsize, int dim, int *map) {
    anns_net_t *an = NULL;
    int j;

    assert(nc > 0);
    assert(nm >= nc);
    assert(nsize > 0);

    an = malloc(sizeof *an);

    assert(an);

    an->z = z;
    an->nc = nc;
    an->nm = nm;
    an->nsize = nsize;
    an->lenc = nc * nsize;
    an->lenm = nm * nsize;
    an->T = 0;
    an->I = 0;
    an->dim = dim;
    an->zeroinsert = 0;

    an->map = map;

    return an;
}

void an_delete(anns_net_t *an) {
    if (an) {
        free(an);
    }
}

double an_eval(anns_net_t *an, double *point) {
    return 0;
}

double *an_ninsert(anns_net_t *an, int from, int count) {
    int j;

    assert(an);
    assert(count > 0);
    assert(an->nc + count <= an->nm);
    assert(from >= 0);

    if (from < an->nc) {
        for (j = an->lenc-1; j >= from * an->nsize; j--) {
            an->z[j+an->nsize * count] = an->z[j];
        }
    }

    if (an->zeroinsert) {
        for (j = from * an->nsize; j < (from + count)*an->nsize; j++) {
            an->z[j] = 0;
        }
    }

    an->nc += count;
    an->lenc += count*an->nsize;

    return an->z + from;
}

double an_distance(anns_net_t *an, int i, int j) {
    double d = 0, tmp, *c1, *c2;
    int k;

    c1 = an_neuron(an, i) + 2;
    c2 = an_neuron(an, j) + 2;

    for (k = 0; k < an->dim; k++) {
        tmp = c1[k] - c2[k];
        d += tmp * tmp;
    }

    return sqrt(d);
}

void an_swap(anns_net_t *an, int i, int j) {
    double d = 0, tmp, *c1, *c2;
    int k;

    c1 = an_neuron(an, i);
    c2 = an_neuron(an, j);

    for (k = 0; k < an->nsize; k++) {
        tmp = c1[k];
        c1[k] = c2[k];
        c2[k] = tmp;
    }
}

double an_min_distance(anns_net_t *an, int i, int *nearest) {
    double min_d = DBL_MAX, d;
    int j, k;

    for (j = 0; j < an->nc; j++) {
        if (i == j) continue;
        d = an_distance(an, i, j);
        if (d < min_d) {
            min_d = d;
            k = j;
        }
    }

    if (nearest) *nearest = k;

    return min_d;
}

double an_max_distance(anns_net_t *an, int i, int *farthest) {
    double max_d = 0, d;
    int j, k;

    for (j = 0; j < an->nc; j++) {
        if (i == j) continue;
        d = an_distance(an, i, j);
        if (d > max_d) {
            max_d = d;
            k = j;
        }
    }

    if (farthest) *farthest = k;

    return max_d;
}

void an_setwidth_nearest(anns_net_t *an, double beta) {
    int i;
    for (i = 0; i < an->nc; i++) {
        an_neuron(an,i)[1] = beta * an_min_distance(an, i, NULL);
    }
}

double *an_nremove(anns_net_t *an, int from, int count) {
    int j;

    assert(an);
    assert(count > 0);
    assert(from >= 0);
    assert(from + count <= an->nc);

    if (from-1 <= an->nc) {
        for (j = from * an->nsize; j <= (from + count)*an->nsize; j++) {
            an->z[j] = an->z[j + an->nsize * count];
        }
    }

    an->nc -= count;
    an->lenc -= count*an->nsize;

    return an->z + from;
}

void an_fprintf(anns_net_t *an, FILE *fout, char *valf, char *valsep, char *nsep, int aux) {
    int i, j;

    assert(an);
    assert(fout);

    if (fout) {
        if (aux) {
            fprintf(fout, "nc=%d nm=%d nsize=%d lenc=%d lenm=%d dim=%d T=%d I=%d zeroinsert=%d map= ",
                    an->nc, an->nm, an->nsize, an->lenc, an->lenm, an->dim, an->T, an->I, an->zeroinsert);

            if (an->map) {
                fprintf(fout, "{%d", an->map[0]);
                for (j = 1; j < an->dim; j++) {
                    fprintf(fout, ",%d", an->map[j]);
                }
                fprintf(fout, "}");
            } else {
                fprintf(fout, "NULL");
            }
            fprintf(fout, "\n");

            fprintf(fout, "eval=<functions %p>\n", an->ann_val);
            fprintf(fout, "eval_der=<functions %p>\n", an->ann_valgrad);
        }

        for (i = 0; i < an->nc; i++) {
            fprintf(fout, valf, an->z[i*an->nsize]);
            for (j = 1; j < an->nsize; j++) {
                fprintf(fout, valsep);
                fprintf(fout, valf, an->z[i*an->nsize + j]);

            }
            fprintf(fout, nsep);
        }
    }
}

double *an_cut(anns_net_t *an, double *buffer, int i, size_t n) {
    int j;

    assert(an);
    assert(buffer);
    assert(i >= 0);
    assert(n > 0);
    assert(i + n <= an->nc);

//    memcpy(buffer, an->z + i*an->nsize, n*sizeof(*an->z)*an->nsize);
//    memmove(an->z + i*an->nsize, an->z + i*an->nsize + n*an->nsize, (an->nc-i-n)*sizeof(*an->z)*an->nsize);

    for (j = 0; j < n*an->nsize ; j++) {
        buffer[j] = an->z[i*an->nsize + j];
    }

    for (j = i*an->nsize ; j < (an->nc-n)*an->nsize ; j++) {
        an->z[j] = an->z[j+n*an->nsize];
    }

    an->nc = an->nc - n;
    an->lenc = an->nc * an->nsize;

    return buffer;
}

void an_paste(anns_net_t *an, double *buffer, int i, size_t n) {
    int j;

    assert(an);
    assert(buffer);
    assert(i >= 0);
    assert(n > 0);
    assert(an->nc + n <= an->nm);

//    if (i < an->nc) {
//        memmove(an->z + i*an->nsize + n*an->nsize, an->z + i*an->nsize, n*sizeof(*an->z)*an->nsize);
//    }
//    memmove(an->z + i*an->nsize, buffer,  n*sizeof(*an->z)*an->nsize);

    if (i < an->nc) {
        for (j = (an->nc-1)*an->nsize; j < i*an->nsize ; j++) {
            an->z[j + n*an->nsize] = an->z[j];
        }
    }

    for (j = i*an->nsize ; j < (i+n)*an->nsize ; j++) {
        an->z[j] = buffer[j - i*an->nsize];
    }

    an->nc = an->nc + n;
    an->lenc = an->nc * an->nsize;
}

void an_overwrite(anns_net_t *an, double *buffer, int i, size_t n) {
    int j;

    assert(an);
    assert(buffer);
    assert(i >= 0);
    assert(n > 0);
    assert(i + n < an->nm);

//    memcpy(an->z + i*an->nsize, buffer, n*sizeof(*an->z)*an->nsize);
    for (j = i*an->nsize ; j < (i+n)*an->nsize ; j++) {
        an->z[j] = buffer[j - i*an->nsize];
    }

    if (i + n > an->nc) {
        an->nc = i + n;
        an->lenc = an->nc * an->nsize;
    }
}

void an_test(void) {
    int nsize = 4, nneurons = 8, i, j;
    double z[nsize * nneurons];
    double buffer[nsize * nneurons];

    int cut_from = 2, cut_n = 2, paste_at = 6;

    anns_net_t *an = an_new(z, nneurons, nneurons, nsize, 2, NULL);

    for (i = 0; i < an->lenc; i++) z[i] = i / nsize;

    printf("[Testing an_net]\n[Net after init]\n");
    an_printf(an);

    printf("[Cut %d neurons from %d position and copy to buffer]\n[Net after cut]\n", cut_n, cut_from);
    an_cut(an, buffer, cut_from, cut_n);
    an_printf(an);
    printf("[Buffer]\n");
    for (i = 0; i < nneurons; i++) {
        for (j = 0; j < nsize; j++) {
            printf("%g ", buffer[i*nsize + j]);
        }
        printf("\n");
    }

    printf("[Paste %d neurons at %d position from buffer]\n[Net after pasting]\n", cut_n, paste_at);
    an_paste(an, buffer, paste_at, cut_n);
    an_printf(an);
    printf("[Buffer]\n");
    for (i = 0; i < nneurons; i++) {
        for (j = 0; j < nsize; j++) {
            printf("%g ", buffer[i*nsize + j]);
        }
        printf("\n");
    }

    printf("[Overwrite %d neurons at %d position from buffer]\n[Net after overwrite]\n", cut_n, cut_from);
    an_overwrite(an, buffer, cut_from, cut_n);
    an_printf(an);
    printf("[Buffer]\n");
    for (i = 0; i < nneurons; i++) {
        for (j = 0; j < nsize; j++) {
            printf("%g ", buffer[i*nsize + j]);
        }
        printf("\n");
    }



    an_delete(an);
}
