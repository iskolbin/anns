#include "anns_cache.h"

anns_cache_t *ac_new(int A, int dim, int nconds, int nnets, int nmem, anns_cond_t **conds, anns_net_t **nets) {
    int t, i, c, p, j;
    int nmax; // максимальное число нейронов
    anns_cache_t *ac = malloc(sizeof(anns_cache_t));

    assert(dim > 0);
    assert(nconds > 0);
    assert(nnets > 0);
    assert(nmem >= 0);
    assert(ac);
    nmax = nets[0]->nc;
    for (t = 1; t < nnets; t++) {
        if (nets[t]->nc > nmax) {
            nmax = nets[t]->nc;
        }
    }

    ac->A = A;
    ac->dim = dim;
    ac->nconds = nconds;
    ac->nnets = nnets;
    ac->nmem = nmem;
    ac->nmax = nmax;
    ac->conds = conds;
    ac->nets = nets;
    ac->pointsselected = 0;

    ac->xm = (nmem > 0) ? d2alloc(nmem, A) : NULL;

    ac->w = calloc(nmax, sizeof(*ac->w));
    ac->a = d2alloc(nmax, dim);
    ac->b = d2alloc(nmax, dim);
    ac->xc = d2alloc(nmax, dim);
    ac->r = calloc(nmax, sizeof(*ac->r));

    ac->phi = calloc(nconds, sizeof *ac->phi);
    ac->phi_x = calloc(nconds, sizeof *ac->phi_x);
    ac->phi_xx = calloc(nconds, sizeof *ac->phi_xx);
    ac->q = calloc(nconds, sizeof *ac->q);
    ac->q_x = calloc(nconds, sizeof *ac->q_x);
    ac->q_xx = calloc(nconds, sizeof *ac->q_xx);

    for (c = 0; c < nconds; c++) {
        ac->phi[c] = calloc(conds[c]->m, sizeof **ac->phi);
        ac->phi_x[c] = calloc(conds[c]->m, sizeof **ac->phi_x);
        ac->phi_xx[c] = calloc(conds[c]->m, sizeof **ac->phi_xx);
        ac->q[c] = calloc(conds[c]->m, sizeof **ac->q);
        ac->q_x[c] = calloc(conds[c]->m, sizeof **ac->q_x);
        ac->q_xx[c] = calloc(conds[c]->m, sizeof **ac->q_xx);

        for (p = 0; p < conds[c]->m; p++) {
            ac->phi[c][p] = calloc(nnets, sizeof ***ac->phi);
            ac->phi_x[c][p] = calloc(nnets, sizeof ***ac->phi_x);
            ac->phi_xx[c][p] = calloc(nnets, sizeof ***ac->phi_xx);
            ac->q[c][p] = calloc(nnets, sizeof ***ac->q);
            ac->q_x[c][p] = calloc(nnets, sizeof ***ac->q_x);
            ac->q_xx[c][p] = calloc(nnets, sizeof ***ac->q_xx);

            for (t = 0; t < nnets; t++) {
                ac->phi[c][p][t] = calloc(nets[t]->nm, sizeof ****ac->phi);
                ac->phi_x[c][p][t] = calloc(nets[t]->nm, sizeof ****ac->phi_x);
                ac->phi_xx[c][p][t] = calloc(nets[t]->nm, sizeof ****ac->phi_xx);
                ac->q_x[c][p][t] = calloc(nets[t]->nm, sizeof ****ac->q_x);
                ac->q_xx[c][p][t] = calloc(nets[t]->nm, sizeof ****ac->q_xx);

                for (i = 0; i < nets[t]->nm; i++) {
                    ac->phi_x[c][p][t][i] = calloc(nets[t]->dim, sizeof *****ac->phi_x);
                    ac->phi_xx[c][p][t][i] = calloc(nets[t]->dim, sizeof *****ac->phi_xx);
                    ac->q_xx[c][p][t][i] = calloc(nets[t]->dim, sizeof *****ac->q_xx);

                    for (j = 0; j < nets[t]->dim; j++) {
                        ac->phi_xx[c][p][t][i][j] = calloc(nets[t]->dim, sizeof ******ac->phi_xx);
                    }
                }
            }
        }
    }

    ac->phi_b = d2alloc(nmax, dim);
    ac->phi_bx = calloc(nmax, sizeof *ac->phi_bx);
    for (i = 0; i < nmax; i++) {
        ac->phi_bx[i] = d2alloc(dim, dim);
    }

    ac->u = calloc(nnets, sizeof(*ac->u));
    ac->u_a = calloc(nnets, sizeof *ac->u_a);
    for (t = 0; t < nnets; t++) {
        ac->u_a[t] = calloc(nets[t]->lenm, sizeof **ac->u_a);
    }

    ac->u_x = calloc(nnets, sizeof *ac->u_x);
    ac->u_xx = calloc(nnets, sizeof *ac->u_xx);
    for (t = 0; t < nnets; t++) {
        ac->u_x[t] = calloc(nets[t]->dim, sizeof **ac->u_x);
        ac->u_xx[t] = calloc(nets[t]->dim, sizeof **ac->u_xx);
    }

    ac->u_xy = calloc(nnets, sizeof *ac->u_xy);
    ac->u_ax = calloc(nnets, sizeof *ac->u_ax);
    ac->u_axx = calloc(nnets, sizeof *ac->u_axx);
    ac->u_axy = calloc(nnets, sizeof *ac->u_axy);

    for (t = 0; t < nnets; t++) {
        ac->u_ax[t] = d2alloc(nets[t]->lenm, nets[t]->dim);
        ac->u_axx[t] = d2alloc(nets[t]->lenm, nets[t]->dim);
        ac->u_xy[t] = d2alloc(nets[t]->dim, nets[t]->dim);
        ac->u_axy[t] = d3alloc(nets[t]->lenm, nets[t]->dim, nets[t]->dim);
    }

    // Создаём контейнеры для предвычисленных членов
    ac->v = calloc(ac->nconds, sizeof *ac->v);
    ac->v_ = calloc(ac->nconds, sizeof *ac->v_);
    ac->vtmp = ac->v;
    for (c = 0; c < ac->nconds; c++) {
        if (conds[c]->nvars > 0) {
            ac->v[c] = calloc(conds[c]->m, sizeof **ac->v);
            for (p = 0; p < conds[c]->m; p++) {
                ac->v[c][p] = calloc(conds[c]->nvars, sizeof ***ac->v);
            }

            ac->v_[c] = calloc(conds[c]->m_, sizeof **ac->v_);
            for (p = 0; p < conds[c]->m_; p++) {
                ac->v_[c][p] = calloc(conds[c]->nvars, sizeof ***ac->v_);
            }
        } else {
            ac->v[c] = NULL;
            ac->v_[c] = NULL;
        }
    }

    // Выделяем память для точек
    ac->points = calloc(ac->nconds, sizeof *ac->points);
    ac->points_ = calloc(ac->nconds, sizeof *ac->points_);
    ac->pointstmp = ac->points;
    for (c = 0; c < ac->nconds; c++) {
        ac->points[c] = d2alloc(conds[c]->m, ac->dim);
        ac->points_[c] = d2alloc(conds[c]->m_, ac->dim);
    }

    ac->csum = calloc(nconds, sizeof *ac->csum);

    ac->g = calloc(A, sizeof *ac->g);

    return ac;
}

void ac_delete(anns_cache_t *ac) {
    int t, c;
    int A, dim, nconds, nnets, nvars, nmem, nmax;
    int *nn, *nd, *na, *np;
    anns_net_t *net;

    if (ac) {
        A = ac->A;
        dim = ac->dim;
        nconds = ac->nconds;
        nnets = ac->nnets;
        nvars = ac->nvars;
        nmem = ac->nmem;
        nmax = ac->nmax;

        free(ac->csum);

        for (t = 0; t < nnets; t++) {
            net = ac->nets[t];
            d3free(ac->u_axy[t], net->lenm, net->dim);
            d2free(ac->u_xy[t], net->dim);
            d2free(ac->u_axx[t], net->lenm);
            d2free(ac->u_ax[t], net->lenm);
        }
        free(ac->u_axy);
        free(ac->u_axx);
        free(ac->u_ax);
        free(ac->u_xy);

        d2free(ac->u_xx, nnets);
        d2free(ac->u_x, nnets);

        d2free(ac->u_a, nnets);
        free(ac->u);

        d3free(ac->phi_bx, nmax, dim);
        d2free(ac->phi_b, nmax);

        d2free(ac->b, nmax);
        d2free(ac->a, nmax);

        free(ac->w);

        if (nmem > 0) d2free(ac->xm, nmem);

        if (ac->pointsselected) {
            ac_switchpoints(ac);
        }

        for (c = 0; c < ac->nconds; c++) {
            d2free(ac->points[c], ac->conds[c]->m);
            d2free(ac->points_[c], ac->conds[c]->m_);

            if (ac->conds[c]->nvars > 0) {

            }
        }
        free(ac->v);
        free(ac->v_);
        free(ac->points);
        free(ac->points_);

        free(ac->g);

        free(ac);
    }
}

void ac_switchpoints(anns_cache_t *ac) {
    int c;

    assert(ac);

    if (ac->pointsselected) {
        ac->points = ac->pointstmp;
        ac->v = ac->vtmp;
        for (c = 0; c < ac->nconds; c++) {
            ac->conds[c]->m = ac->conds[c]->mtmp;
        }
    } else {
        ac->points = ac->points_;
        ac->v = ac->v_;
        for (c = 0; c < ac->nconds; c++) {
            ac->conds[c]->m = ac->conds[c]->m_;
        }
    }
    ac->pointsselected = !ac->pointsselected;
}

void ac_push(anns_cache_t *ac, double *x) {
    double **xm = NULL, *tmp;
    int i;

    if (ac && ac->xm) {
        xm = ac->xm;

        // Сохраняем в последний вектор
        for (i = 0; i < ac->A; i++) {
            xm[ac->nmem-1][i] = x[i];
        }

        // Обновляем ссылки - последний вектор теперь первый
        tmp = xm[ac->nmem-1];
        for (i = 1; i < ac->nmem; i++){
            xm[i] = xm[i-1];
        }
        xm[0] = tmp;
    }
}

void ac_put(anns_cache_t *ac, double *x, int i) {
    int j;

    assert(ac);
    assert(x);
    assert(i >= 0 && i < ac->nmem);

    for (j = 0;j < ac->A; j++) {
        ac->xm[i][j] = x[j];
    }
}

double *ac_get(anns_cache_t *ac, int i) {
    assert(ac);
    assert(i >= 0 && i < ac->nmem);
    return ac->xm[i];
}

void ac_eval_v(anns_cache_t *ac, int c, int iv, ANNS_FUNPTR(f), void *instance) {
    int p;

    assert(f);
    assert(ac);
    assert(ac->nconds);
//    assert(ac->nv[c]);
//    assert(iv >= 0 && iv < ac->nv[c]);

//    for (p = 0; p < ac->np[c]; p++) {
    for (p = 0; p < ac->conds[c]->m; p++) {
        ac->v[c][p][iv] = f(ac->points[c][p], ac->dim, instance);
    }

    for (p = 0; p < ac->conds[c]->m_; p++) {
        ac->v_[c][p][iv] = f(ac->points_[c][p], ac->dim, instance);
    }
}

void ac_load_v(anns_cache_t *ac, int c, int iv, double *darr) {
    int p;

    assert(ac);
    assert(c >= 0 && c < ac->nconds);
    assert(iv >= 0 && iv < ac->conds[c]->nvars);
    assert(darr);

    for (p = 0; p < ac->conds[c]->m; p++) {
        ac->v[c][p][iv] = darr[p];
    }


}

void ac_set_v(anns_cache_t *ac, int c, int iv, double val) {
    int p;

    assert(ac);
    assert(c >= 0 && c < ac->nconds);
    assert(iv >= 0 && iv < ac->conds[c]->nvars);

    for (p = 0; p < ac->conds[c]->m; p++) {
        ac->v[c][p][iv] = val;
    }
}

void ac_v_fout(FILE *fout, anns_cache_t *ac) {
    int c, p, iv;

    for (c = 0; c < ac->nconds; c++) {
        fprintf(fout, "[Condition %d]\n", c);
//        for (p = 0; p < ac->np[c]; p++) {
        for (p = 0; p < ac->conds[c]->m; p++) {
            fprintf(fout, "[");
//            for (iv = 0; iv < ac->nv[c]; iv++) {
            for (iv = 0; iv < ac->conds[c]->nvars; iv++) {
                fprintf(fout, "%g ", ac->v[c][p][iv]);
            }
            fprintf(fout, "]");
        }
        fprintf(fout,"\n");
    }
}

void ac_p_fout(FILE *fout, anns_cache_t *ac) {
    int c, p, j;

    for (c = 0; c < ac->nconds; c++) {
        fprintf(fout, "[Condition %d]\n", c);
        for (p = 0; p < ac->conds[c]->m; p++) {
            fprintf(fout, "(");
            fprintf(fout, "%g", ac->points[c][p][0]);
            for (j = 1; j < ac->dim; j++) {
                fprintf(fout, ",%g", ac->points[c][p][j]);
            }
            fprintf(fout, ") ");
        }
        fprintf(fout,"\n");
    }
}
