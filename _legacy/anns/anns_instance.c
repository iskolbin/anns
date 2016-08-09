#include "anns_instance.h"


static char *generate_file_name(int problemtag, int typetag, int index, int subindex);
static char typetags[][8] = {"dir", "isrc", "ibnd", "iicd", "icft"};
static char problemtags[][8] = {"ell", "par", "hyp", "app"};

static char *generate_file_name(int typetag, int problemtag, int index, int subindex) {
    char *tmp = calloc(1024, sizeof(char));
    char *name = NULL;
    int len;

    len = sprintf(tmp, "%s/%s_%s_%d_%d.txt", OUTPUT_DIR, typetags[typetag%100], problemtags[problemtag%200], index, subindex);
    name = calloc(len + 1, sizeof(char));
    strcpy(name, tmp);

    return name;
}


anns_instance_t *ai_new(int data[], void *instance) {
    anns_instance_t *ai = malloc(sizeof(anns_instance_t));
    anns_net_t *net;
    anns_cond_t *cond;
    int c, t, d, cursor = 0, I = 0;;

    assert(data);
    assert(ai);

    // Считываем описательные параметры и генерируем на их основе название файла
    ai->typetag = data[cursor++];
    ai->problemtag = data[cursor++];
    ai->index = data[cursor++];
    ai->subindex = data[cursor++];
    ai->filename = generate_file_name(ai->typetag, ai->problemtag,  ai->index, ai->subindex);

    // Считываем глобальные переменные: размерность пространства, число уравнений, число сетей,
    // число запоминаемых наборов сетей
    ai->F = data[cursor++];
    ai->tdelta = data[cursor++];

    ai->dim = data[cursor++];
    ai->nconds = data[cursor++];
    ai->nnets = data[cursor++];
    ai->nmem = data[cursor++];

    ai->nets = calloc(ai->nnets, sizeof *ai->nets);
    ai->conds = calloc(ai->nconds, sizeof *ai->conds);

    ai->M = 0;
    for (c = 0; c < ai->nconds; c++) {
        ai->conds[c] = malloc(sizeof **ai->conds);
        cond = ai->conds[c];
        cond->m = data[cursor++];
        cond->mtmp = cond->m;
        cond->m_ = data[cursor++];
        ai->M += cond->m;
        cond->nvars = data[cursor++];
        cond->nt = ai->nnets;
        cond->eval = calloc(cond->nt, sizeof *cond->eval);
        cond->A_g = calloc(cond->nt, sizeof *cond->A_g);
    }

    for (c = 0; c < ai->nconds; c++) {
        switch (ai->tdelta) {
            case DELTA_MANUAL: ai->conds[c]->delta = 1; break;
            case DELTA_INV_COND: ai->conds[c]->delta = 1. / ai->conds[c]->m; break;
            case DELTA_INV_ALL: ai->conds[c]->delta = 1. / ai->M; break;
        }
    }

    for (t = 0; t < ai->nnets; t++) {
        ai->nets[t] = malloc(sizeof **ai->nets);
        net = ai->nets[t];
        net->T = data[cursor++];
        net->nm = data[cursor++];
        net->nc = net->nm;
        net->dim = data[cursor++];
        if (data[cursor++]) {
            net->map = calloc(net->dim, sizeof *net->map);
            for (d = 0; d < net->dim; d++) {
                net->map[d] = data[cursor++];
            }
        } else {
            net->map = NULL;
        }

        for (c = 0; c < ai->nconds; c++) {
            ai->conds[c]->eval[t] = data[cursor++];
        }
    }

    // Вычисляем вспомогательные переменные - общее число переменных сетей, индексы сетей, всеобщее число переменных
    ai->A = 0;
    ai->Aw = 0;


    for (t = 0; t < ai->nnets; t++) {
        net = ai->nets[t];
        net->I = I;

        if (net->T & RADIAL) {
            net->nsize = net->dim + 2;
        } else if (net->T & ELLIPTIC) {
            net->nsize = 2*net->dim + 1;
        } else {
            // Очевидно ошибка, пока не введены другие сети
        }
        net->lenc = net->nsize * net->nc;
        net->lenm = net->nsize * net->nm;
        I += net->lenm;
        ai->A += net->lenm;
        ai->Aw += net->nm;
    }

    ai->z = calloc(ai->A, sizeof *ai->z);
    ai->zbest = calloc(ai->A, sizeof *ai->zbest);
    ai->zw = calloc(ai->Aw, sizeof *ai->zw);

    // Даем смещения
    ai->nets[0]->z = ai->z;
    ai->nets[0]->zw = ai->zw;
    for (t = 1; t < ai->nnets; t++) {
        net = ai->nets[t];
        net->z = ai->z + ai->nets[t-1]->lenm;
        net->zw = ai->zw + ai->nets[t-1]->nm;
    }

    // Вычисляем число точек
    ai->M = 0;
    for (c = 0; c < ai->nconds; c++) {
        ai->M += ai->conds[c]->m;
    }

    // Создаём контейнер для вычислений
    ai->cache = ac_new(ai->A, ai->dim, ai->nconds, ai->nnets, ai->nmem, ai->conds, ai->nets);

    // Создаём точку для свёртки
    ai->dpoint = calloc(ai->dim, sizeof *ai->dpoint);

    ai->instance = instance;

    ai->t = 0;
    ai->nsteps = ai->nmem;
    ai->start = 0;
    ai->step = 1./ (ai->nsteps - 1);

    ai->preeval = NULL;

    return ai;
}

void ai_delete(anns_instance_t *ai) {
    int i, t, c;

    if (ai) {
        ac_delete(ai->cache);

        free(ai->dpoint);
        free(ai->z);
        free(ai->zbest);
        free(ai->filename);

        for (t = 0; t < ai->nnets; t++) {
            if (ai->nets[t]->map) free(ai->nets[t]->map);
            free(ai->nets[t]);
        }
        free(ai->nets);

        for (c = 0; c < ai->nconds; c++) {
            free(ai->conds[c]->eval);
            free(ai->conds[c]->A_g);
            free(ai->conds[c]);
        }
        free(ai->conds);

        free(ai);
    }
}

void ai_switchz(anns_instance_t *ai) {
    int t;

    if (ai_isbestz(ai)) {
        ai->z = ai->ztmp;
        for (t = 0; t < ai->nnets; t++) {
            ai->nets[t]->nc = ai->nets[t]->nctmp;
        }
    } else {
        ai->z = ai->zbest;
        for (t = 0; t < ai->nnets; t++) {
            ai->nets[t]->nc = ai->nets[t]->ncbest;
        }
    }
}

void ai_savez(anns_instance_t *ai) {
    int i, t;

    for (i = 0; i < ai->A; i++) {
        ai->zbest[i] = ai->z[i];
    }

    for (t = 0; t < ai->nnets; t++) {
        ai->nets[t]->ncbest = ai->nets[t]->nc;
        ai->nets[t]->zbest = ai->zbest + ai->nets[t]->I;
    }
}

void ai_conditions(anns_instance_t *ai, fbvp_t *fbvp, gbvp_t *gbvp, double *delta) {
    int c, t, cursor1 = 0, cursor2 = 0;

    for (c = 0; c < ai->nconds; c++) {
        ai->conds[c]->A = fbvp[cursor1++];

        for (t = 0; t < ai->nnets; t++) {
            ai->conds[c]->A_g[t] = gbvp[cursor2++];
        }

        if (delta) {
            ai->conds[c]->delta *= delta[c];
        }
    }
}

void ai_updatenets(anns_instance_t *ai, double *newz, int *nc) {
    int l = 0, k, t;

    assert(ai);
    assert(newz);

    for (t = 0; t < ai->nnets; t++) {
        for (k = 0; k < ai->nets[t]->lenm; k++) {
            ai->z[l] = newz[l];
            l++;
        }

        if (nc) {
            ai->nets[t]->lenc = nc[t] * ai->nets[t]->nsize;
            ai->nets[t]->nc = nc[t];
        }
    }
}

void ai_fprintf(anns_instance_t *ai, FILE *fout, char *valf, char *valsep, char *nsep, int aux) {
    int c, t;

    assert(ai);
    assert(fout);

    if (fout) {
        fprintf(fout, "[Instance]\n");
        fprintf(fout, "tags=%d.%d index=%d.%d filename=\"%s\"\n", ai->problemtag, ai->typetag, ai->index, ai->subindex, ai->filename);
        fprintf(fout, "dim=%d nconds=%d nnets=%d nmem=%d A=%d M=%d F=%d\n", ai->dim, ai->nconds, ai->nnets, ai->nmem, ai->A, ai->M, ai->F);
        fprintf(fout, "[Nets (%d)]\n", ai->nnets);
        for (t = 0; t < ai->nnets; t++) {
            fprintf(fout, "[Net %d]\n", t);
            an_fprintf(ai->nets[t], fout, valf, valsep, nsep, aux);
        }
        fprintf(fout, "\n[Conditions (%d)]\n", ai->nconds);
        for (c = 0; c < ai->nconds; c++) {
            fprintf(fout, "[Condition %d]\n", c);
            ad_fprintf(ai->conds[c], fout);
        }
    }
}

double *ai_collectweights(anns_instance_t *ai) {
    int t, i, i_ = 0;

    assert(ai);

    for (t = 0; t < ai->nnets; t++) {
        for (i = 0; i < ai->nets[t]->nc; i++) {
            ai->zw[i_++] = ai->nets[t]->z[i*ai->nets[t]->nsize];
        }
    }

    return ai->zw;
}

double *ai_expandweights(anns_instance_t *ai) {
    int t, i, i_ = 0;

    assert(ai);

    for (t = 0; t < ai->nnets; t++) {
        for (i = 0; i < ai->nets[t]->nc; i++) {
            ai->nets[t]->z[i*ai->nets[t]->nsize] = ai->zw[i_++];
        }
    }

    return ai->z;
}

