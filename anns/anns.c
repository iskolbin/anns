#include "anns.h"

double anns_val_f(double *x, int A, void *instance) {
    anns_instance_t *ai = instance;
    anns_cache_t *ac = ai->cache;
    anns_cond_t *cond;
    anns_net_t *net;
    double *point = NULL, Fp, functional = 0., functional_c, *dpoint = ai->dpoint;
    double *u = ac->u, **u_x = ac->u_x, **u_xx = ac->u_xx, ***u_xy = ac->u_xy;
    int c, p, t, eval, j, *k = NULL;

    // Проходим по уравнениям
    ai->Jmax = -1e10;
        if (ai->preeval) {
        ai->preeval(x, ai);
    }


    for (c = 0; c < ai->nconds; c++) {
        cond = ac->conds[c];

        // Проходим по точкам, ассоцированным с уравнением
        functional_c = 0;
        for (p = 0; p < cond->m; p++) {
            point = ac->points[c][p];

            // Вычисляем нейросети в точке (если они участвуют в вычислениях)
            for (t = 0; t < ac->nnets; t++) {
                net = ac->nets[t];
                eval = cond->eval[t];

                // Проверяем условие - используется ли сеть в уравнении
                if (eval >= 0) {
                    net->ann_val[eval % MIXED_DERIVATIVES](x, c, p, t, eval & MIXED_DERIVATIVES, ai);
//                    if (net->map) {
//                        for (j = 0; j < net->dim; j++) {
//                            dpoint[j] = point[net->map[j]];
//                        }
//                        net->ann_val[eval % MIXED_DERIVATIVES](dpoint, x, t, eval & MIXED_DERIVATIVES, ai);
//                    } else {
//                        net->ann_val[eval % MIXED_DERIVATIVES](point, x, t, eval & MIXED_DERIVATIVES, ai);
//                    }
                }
            }

            // Вычисляем уравенение с нейросетями в точке
            Fp = cond->A(point, u, u_x, u_xx, u_xy, ac->v[c] ? ac->v[c][p] : NULL, ai);

            // Обновляем значение функционала
            switch (ai->F) {
                case FUNCTIONAL_SQUARE: case FUNCTIONAL_SQUAREROOT: functional_c += Fp*Fp; break;
                case FUNCTIONAL_ABSOLUTE: functional_c += fabs(Fp); break;
                case FUNCTIONAL_LINEAR: functional_c += Fp; break;
            }

            if (functional_c > ai->Jmax) {
                ai->Jmax = functional_c;
                ai->cJmax = c;
                ai->pJmax = p;
            }
        }

        functional += cond->delta*functional_c;
    }

    if (ai->F == FUNCTIONAL_SQUAREROOT) {
        functional = sqrt(functional);
    }

    return functional;
}

double anns_val_f_check(double *x, int A, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    double f_;

    ai_switchpoints(ai);
    f_ = anns_val_f(x, A, anns_instance);
    ai_switchpoints(ai);

    return f_;
//    ai->cache->points
}

void anns_grad_f(double *g, double *x, int A, void *instance) {
    // Вычисление градиента невозможно без вычисления значения
    anns_valgrad_f(g, x, A, instance);
}

double anns_valgrad_f(double *g, double *x, int A, void *instance) {
    anns_instance_t *ai = instance;
    anns_cache_t *ac = ai->cache;
    anns_cond_t *cond;
    anns_net_t *net;
    double *point = NULL, Fp, Fp_a, functional = 0., functional_c, *dpoint = ai->dpoint;
    double *u = ac->u, **u_x = ac->u_x, **u_xx = ac->u_xx, ***u_xy = ac->u_xy, *g_c = ac->g,
        **u_a = ac->u_a, ***u_ax = ac->u_ax, ***u_axx = ac->u_axx, ****u_axy = ac->u_axy;
    int i, j, c, p, t, eval, I, *k = NULL;

    for (i = 0; i < ai->A; i++) {
        g[i] = 0;
    }

    // Проходим по уравнениям
    ai->Jmax = -1e10;

    if (ai->preeval) {
        ai->preeval(x, ai);
    }

    for (c = 0; c < ai->nconds; c++) {
        cond = ac->conds[c];
        functional_c = 0;

        for (i = 0; i < ai->A; i++) {
            g_c[i] = 0;
        }

        // Проходим по точкам, ассоцированным с уравнением
        for (p = 0; p < cond->m; p++) {
            point = ac->points[c][p];

            // Вычисляем нейросети в точке (если они участвуют в вычислениях)
            for (t = 0; t < ai->nnets; t++) {
                net = ac->nets[t];
                eval = cond->eval[t];

                // Проверяем условие - используется ли сеть в уравнении
                if (eval >= 0) {

                    net->ann_valgrad[eval % MIXED_DERIVATIVES](x, c, p, t, eval & MIXED_DERIVATIVES, ai);
//                    if (net->map) {
//                        for (j = 0; j < net->dim; j++) {
//                            dpoint[j] = point[net->map[j]];
//                        }
//                        net->ann_valgrad[eval % MIXED_DERIVATIVES](dpoint, x, t, eval & MIXED_DERIVATIVES, ai);
//                    } else {
//                        net->ann_valgrad[eval % MIXED_DERIVATIVES](point, x, t, eval & MIXED_DERIVATIVES, ai);
//                    }
                }
            }

            // Вычисляем уравенение с нейросетями в точке
            Fp = cond->A(point, u, u_x, u_xx, u_xy, ac->v[c] ? ac->v[c][p] : NULL, ai);

            // Обновляем значение функционала
            switch (ai->F) {
                case FUNCTIONAL_SQUARE: case FUNCTIONAL_SQUAREROOT: functional_c += Fp*Fp; break;
                case FUNCTIONAL_ABSOLUTE: functional_c += fabs(Fp); break;
                case FUNCTIONAL_LINEAR: functional_c += Fp; break;
            }

            if (functional_c > ai->Jmax) {
                ai->Jmax = functional_c;
                ai->cJmax = c;
                ai->pJmax = p;
            }

            // Обновляем градиент
            for (t = 0; t < ai->nnets; t++) {
                eval = cond->eval[t];
                if (eval >= 0) {
                    I = ai->nets[t]->I;
                    for (i = I; i < I+ai->nets[t]->lenc; i++) {
//                        printf("[[%p]]", cond->A_g[t]);
                        Fp_a = cond->A_g[t](point, u, u_x, u_xx, u_xy, u_a[t][i-I], u_ax[t][i-I], u_axx[t][i-I], u_axy[t][i-I], ac->v[c] ? ac->v[c][p] : NULL, ai);

                        switch (ai->F) {
                            case FUNCTIONAL_SQUARE: case FUNCTIONAL_SQUAREROOT: g_c[i] += 2 * Fp_a * Fp; break;
                            case FUNCTIONAL_ABSOLUTE: g_c[i] += (Fp > 0 ? Fp_a: Fp < 0 ? -Fp_a : 0); break;
                            case FUNCTIONAL_LINEAR: g_c[i] += Fp_a; break;
                        }
                    }
                }
            }
        }

        functional += cond->delta*functional_c;

        for (i = 0; i < ai->A; i++) {
            g[i] += cond->delta * g_c[i];
        }
    }

    if (ai->F == FUNCTIONAL_SQUAREROOT) {
        functional = sqrt(functional);
        for (i = 0; i < ai->A; i++) {
            g[i] *= (0.5/functional);
        }
    }

    return functional;
}

double anns_eval(double *point, double *x, void *anns_instance) {
    int t;
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    double s = 0.0;

    for (t = 0; t < ai->nnets; t++) {
//        ai->nets[t]->ann_val[0](c, p, t, 0, anns_instance);
//        s += ac->u[t];
        s += ai->nets[t]->ann_eval(x, point, t, anns_instance);
//        ai->nets[t]->ann_eval(x, point, t, anns_instance);
    }

    return s;
}


void anns_csum(double *x, int A, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net;
    double *point = NULL, Fp, functional_c, *dpoint = ai->dpoint;
    double *u = ac->u, **u_x = ac->u_x, **u_xx = ac->u_xx, ***u_xy = ac->u_xy, *g_c = ac->g,
        **u_a = ac->u_a, ***u_ax = ac->u_ax, ***u_axx = ac->u_axx, ****u_axy = ac->u_axy;
    int i, c, p, t, eval, j, *k = NULL;


    for (c = 0; c < ai->nconds; c++) {
        ac->csum[c] = 0;
        functional_c = 0;

        // Проходим по точкам, ассоцированным с уравнением
        for (p = 0; p < ai->conds[c]->m; p++) {
            point = ac->points[c][p];

            // Вычисляем нейросети в точке (если они участвуют в вычислениях)
            for (t = 0; t < ai->nnets; t++) {
                net = ai->nets[t];
                eval = ai->conds[c]->eval[t];

                // Проверяем условие - используется ли сеть в уравнении
                if (eval >= 0) {
                    net->ann_val[eval % MIXED_DERIVATIVES](x, c, p, t, eval & MIXED_DERIVATIVES, ai);
//                    if (net->map) {
//                        for (j = 0; j < net->dim; j++) {
//                            dpoint[j] = point[net->map[j]];
//                        }
//                        net->ann_val[eval % MIXED_DERIVATIVES](dpoint, x, t, eval & MIXED_DERIVATIVES, ai);
//                    } else {
//                        net->ann_val[eval % MIXED_DERIVATIVES](point, x, t, eval & MIXED_DERIVATIVES, ai);
//                    }
                }
            }

            // Вычисляем уравенение с нейросетями в точке
            Fp = ai->conds[c]->A(point, u, u_x, u_xx, u_xy, ac->v[c] ? ac->v[c][p] : NULL, ai);

            // Обновляем значение функционала
            switch (ai->F) {
                case FUNCTIONAL_SQUARE: case FUNCTIONAL_SQUAREROOT: functional_c += Fp*Fp; break;
                case FUNCTIONAL_ABSOLUTE: functional_c += fabs(Fp); break;
                case FUNCTIONAL_LINEAR: functional_c += Fp; break;
            }


        }

        ac->csum[c] = ai->conds[c]->delta*functional_c;
    }
}

double anns_evalsingle(int t, double *point, double *x, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
//    ai->nets[t]->ann_val[0](point, x, t, 0, ai);
    return ai->cache->u[t];
}


void anns_grad_f_wonly(double *g, double *x, int A, void *anns_instance) {
    anns_valgrad_f_wonly(g, x, A, anns_instance);
}


double anns_val_f_wonly(double *x, int A, void *instance) {
    anns_instance_t *ai = instance;
    anns_cache_t *ac = ai->cache;
    anns_cond_t *cond;
    anns_net_t *net;
    double *point = NULL, Fp, functional = 0., functional_c, *dpoint = ai->dpoint;
    double *u = ac->u, **u_x = ac->u_x, **u_xx = ac->u_xx, ***u_xy = ac->u_xy;
    int c, p, t, eval, j, *k = NULL;

    // Проходим по уравнениям
    ai->Jmax = -1e10;
    for (c = 0; c < ai->nconds; c++) {
        cond = ac->conds[c];

        // Проходим по точкам, ассоцированным с уравнением
        functional_c = 0;
        for (p = 0; p < cond->m; p++) {
            point = ac->points[c][p];

            // Вычисляем нейросети в точке (если они участвуют в вычислениях)
            for (t = 0; t < ac->nnets; t++) {
                net = ac->nets[t];
                eval = cond->eval[t];

                // Проверяем условие - используется ли сеть в уравнении
                if (eval >= 0) {
                    net->ann_val_wonly[eval % MIXED_DERIVATIVES](x, c, p, t, eval & MIXED_DERIVATIVES, ai);
//                    if (net->map) {
//                        for (j = 0; j < net->dim; j++) {
//                            dpoint[j] = point[net->map[j]];
//                        }
//                        net->ann_val[eval % MIXED_DERIVATIVES](dpoint, x, t, eval & MIXED_DERIVATIVES, ai);
//                    } else {
//                        net->ann_val[eval % MIXED_DERIVATIVES](point, x, t, eval & MIXED_DERIVATIVES, ai);
//                    }
                }
            }

            // Вычисляем уравенение с нейросетями в точке
            Fp = cond->A(point, u, u_x, u_xx, u_xy, ac->v[c] ? ac->v[c][p] : NULL, ai);

            // Обновляем значение функционала
            switch (ai->F) {
                case FUNCTIONAL_SQUARE: functional_c += Fp*Fp; break;
                case FUNCTIONAL_ABSOLUTE: functional_c += fabs(Fp); break;
                case FUNCTIONAL_LINEAR: functional_c += Fp; break;
            }

            if (functional_c > ai->Jmax) {
                ai->Jmax = functional_c;
                ai->cJmax = c;
                ai->pJmax = p;
            }
        }

        functional += cond->delta*functional_c;
    }

    return functional;
}

double anns_valgrad_f_wonly(double *g, double *x, int A, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_cond_t *cond;
    anns_net_t *net;
    double *point = NULL, Fp, Fp_a, functional = 0., functional_c, *dpoint = ai->dpoint;
    double *u = ac->u, **u_x = ac->u_x, **u_xx = ac->u_xx, ***u_xy = ac->u_xy, *g_c = ac->g,
        **u_a = ac->u_a, ***u_ax = ac->u_ax, ***u_axx = ac->u_axx, ****u_axy = ac->u_axy;
    int i, j, c, p, t, eval, I, *k = NULL;

    for (i = 0; i < ai->A; i++) {
        g[i] = 0;
    }

    // Проходим по уравнениям
    ai->Jmax = -1e10;
    for (c = 0; c < ai->nconds; c++) {
        cond = ac->conds[c];
        functional_c = 0;

        for (i = 0; i < ai->A; i++) {
            g_c[i] = 0;
        }

        // Проходим по точкам, ассоцированным с уравнением
        for (p = 0; p < cond->m; p++) {
            point = ac->points[c][p];

            // Вычисляем нейросети в точке (если они участвуют в вычислениях)
            for (t = 0; t < ai->nnets; t++) {
                net = ac->nets[t];
                eval = cond->eval[t];

                // Проверяем условие - используется ли сеть в уравнении
                if (eval >= 0) {
                    net->ann_valgrad_wonly[eval % MIXED_DERIVATIVES](x, c, p, t, eval & MIXED_DERIVATIVES, ai);
//                    if (net->map) {
//                        for (j = 0; j < net->dim; j++) {
//                            dpoint[j] = point[net->map[j]];
//                        }
//                        net->ann_valgrad[eval % MIXED_DERIVATIVES](c, p, t, eval & MIXED_DERIVATIVES, ai);
//                    } else {
//                        net->ann_valgrad[eval % MIXED_DERIVATIVES](c, p, t,  eval & MIXED_DERIVATIVES, ai);
//                    }
                }
            }

            // Вычисляем уравенение с нейросетями в точке
            Fp = cond->A(point, u, u_x, u_xx, u_xy, ac->v[c] ? ac->v[c][p] : NULL, ai);

            // Обновляем значение функционала
            switch (ai->F) {
                case FUNCTIONAL_SQUARE: functional_c += Fp*Fp; break;
                case FUNCTIONAL_ABSOLUTE: functional_c += fabs(Fp); break;
                case FUNCTIONAL_LINEAR: functional_c += Fp; break;
            }

            if (functional_c > ai->Jmax) {
                ai->Jmax = functional_c;
                ai->cJmax = c;
                ai->pJmax = p;
            }

            // Обновляем градиент
            for (t = 0; t < ai->nnets; t++) {
                eval = cond->eval[t];
                if (eval >= 0) {
                    I = ai->nets[t]->I;
                    for (i = I; i < I+ai->nets[t]->lenc; i++) {
                        Fp_a = cond->A_g[t](point, u, u_x, u_xx, u_xy, u_a[t][i-I], u_ax[t][i-I], u_axx[t][i-I], u_axy[t][i-I], ac->v[c] ? ac->v[c][p] : NULL, ai);

                        switch (ai->F) {
                            case FUNCTIONAL_SQUARE: g_c[i] += 2 * Fp_a * Fp; break;
                            case FUNCTIONAL_ABSOLUTE: g_c[i] += (Fp > 0 ? Fp_a: Fp < 0 ? -Fp_a : 0); break;
                            case FUNCTIONAL_LINEAR: g_c[i] += Fp_a; break;
                        }
                    }
                }
            }
        }

        functional += cond->delta*functional_c;

        for (i = 0; i < ai->A; i++) {
            g[i] += cond->delta * g_c[i];
        }
    }

    return functional;
}


double *anns_refine(anns_instance_t *ai) {
    int i, t;
    double *buffer = NULL;
    double J, J1;
    double jtol = 0;

    assert(ai);

    for (t = 0; t < ai->nnets; t++) {
        buffer = calloc(ai->nets[t]->nsize, sizeof(*buffer) );

        J = anns_val_f(ai->z, ai->A, ai);

        for (i = 0; i < ai->nets[t]->nc; i++) {
            an_cut(ai->nets[t], buffer, i, 1);
            J1 = anns_val_f(ai->z, ai->A, ai);

            if (J1 - J > jtol) {
                an_paste(ai->nets[t], buffer, i, 1);
            }
        }
        free(buffer);
    }
}
