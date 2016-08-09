#include "rbf_gauss.h"
#include "../anns_tablecalc.h"

static anns_tablecalc_t *atexp = NULL;
#define MIN_EPOW (-40.0)


double rbf_gauss_eval(double *z, double *y, int t, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_net_t *net = ai->nets[t];
    int i, j, I, nc = net->nc, from = net->I, nsize = net->nsize;
    double p = 0, phi, epow, r, w, b, a, xc;

    for (i = 0, I = from; i < nc; i++, I += nsize) {
        r = 0;
        w = z[I];
        b = z[I + 1];
        a = 1 / (b*b);

        for (j = 0; j < net->dim; j++) {
            xc = y[j] - z[I + 2 + j];
            r += xc * xc;
        }

        epow = -0.5 * (r * a);
        if (epow > MIN_EPOW) {
            phi = at_fasteval(atexp, epow);
            p += w * phi;
        }
    }

    return p;
}

void rbf_gauss_val(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *an = ac->nets[t];
    int i, j, I, dim = an->dim, nc = an->nc, nsize = an->nsize;
    double p_ = 0., *point = ac->points[c][p];
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi[c][p][t];

    for (i = 0, I = an->I; i < nc; i++, I += nsize) {
        r[i] = 0.;
        w[i] = x[I];
        b[i][0] = x[I + 1];
        a[i][0] = 1 / (b[i][0]*b[i][0]);

        for (j = 0; j < dim; j++) {
            xc[i][j] = point[j] - x[I + 2 + j];
            r[i] += xc[i][j] * xc[i][j];
            xc[i][j] *= a[i][0];
        }
    }

    for (i = 0, I = an->I; i < nc; i++, I += nsize) {
        phi[i] = exp(-0.5 * (r[i] * a[i][0]));

        if (phi[i] != 0) {
            p_ += w[i] * phi[i];
        }
    }

    ac->u[t] = p_;
}

void rbf_gauss_val_d1(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *an = ac->nets[t];
    int i, j, dim = an->dim, nc = an->nc;
    double u, p_x;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi[c][p][t],
        **phi_x = ac->phi_x[c][p][t], **phi_b = ac->phi_b, *u_x = ac->u_x[t];

    rbf_gauss_val(x, c, p, t, mixed_derivatives, anns_instance);

    u = ac->u[t];

    for (j = 0; j < dim; j++) {
        p_x = 0.;
        for (i = 0; i < nc; i++) {
            if (phi[i] != 0) {
                phi_x[i][j] = -xc[i][j] * phi[i];
                p_x += w[i] * phi_x[i][j];
            }
        }
        u_x[j] = p_x;
    }
}

void rbf_gauss_val_d2(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *an = ac->nets[t];
    int i, j, k, nc = an->nc, dim = an->dim;
    double q = 0., u, p_xx;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b;
    double **xc = ac->xc;


    double *phi = ac->phi[c][p][t];


     double   **phi_x = ac->phi_x[c][p][t], **phi_b = ac->phi_b, *u_x = ac->u_x[t];
     double ***phi_xx = ac->phi_xx[c][p][t],
        *u_xx = ac->u_xx[t], **u_xy = ac->u_xy[t];

    rbf_gauss_val_d1(x, c, p, t, mixed_derivatives, anns_instance);

    u = ac->u[t];

    for (j = 0; j < dim; j++) {
        p_xx = 0;
        for (i = 0; i < nc; i++) {
            if (phi[i] != 0) {
                phi_xx[i][j][j] = -phi_x[i][j]*xc[i][j] - a[i][0]*phi[i];
                p_xx += w[i] * phi_xx[i][j][j];
            }
        }
        u_xx[j] = p_xx;

        if (mixed_derivatives) {
            for (k = 0; k < j; k++) {
                p_xx = 0;
                for (i = 0; i < nc; i++) {
                    if (phi[i] != 0) {
                        phi_xx[i][j][k] = xc[i][j] * xc[i][k] * phi[i];
                        p_xx += w[i] * phi_xx[i][j][k];
                    }
                }
                u_xy[j][k] = p_xx;
            }
        }
    }
}


void rbf_gauss_valgrad(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *an = ac->nets[t];
    int i, j, m, s, I, nc = an->nc, dim = an->dim, nsize = an->nsize;
    double p_ = 0., *point = ac->points[c][p];;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi[c][p][t],
        *u_a = ac->u_a[t], **phi_x = ac->phi_x[c][p][t], **phi_b = ac->phi_b;

    for (i = 0, I = an->I; i < nc; i++, I += nsize) {
        r[i] = 0.;
        w[i] = x[I];
        b[i][0] = x[I + 1];
        a[i][0] = 1 / (b[i][0]*b[i][0]);

        for (j = 0; j < dim; j++) {
            xc[i][j] = point[j] - x[I + 2 + j];
            r[i] += xc[i][j] * xc[i][j];
            xc[i][j] *= a[i][0];
        }
    }

    for (i = 0, I = an->I; i < an->nc; i++, I += nsize) {

        phi[i] = exp(-0.5 * (r[i] * a[i][0]));

        if (phi[i] != 0) {
            p_ += w[i] * phi[i];

            for (j = 0; j < dim; j++) {
                phi_x[i][j] = -xc[i][j] * phi[i];
            }
        }
    }

    ac->u[t] = p_;

    // Вычисление градиента
    for (i = 0, s = 0; i < nc; i++, s+= nsize) {
        if (phi[i] != 0) {
            u_a[s] = phi[i];

            phi_b[i][0] = phi[i] * r[i] * a[i][0] / b[i][0];
            u_a[s+1] = phi_b[i][0] * w[i];

            for (m = 0; m < dim; m++) {
                u_a[s+2+m] = -phi_x[i][m] * w[i];
            }
        } else {
            u_a[s] = 0;
            u_a[s+1] = 0;
            for (m = 0; m < dim; m++) {
                u_a[s+2+m] = 0;
            }
        }
    }
}

void rbf_gauss_valgrad_d1(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *an = ac->nets[t];
    int i, j, m, s, k, dim = an->dim, nc = an->nc, nsize = an->nsize;
    double u, p_x;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi[c][p][t],
        *u_a = ac->u_a[t], **u_ax = ac->u_ax[t], **phi_x = ac->phi_x[c][p][t], **phi_b = ac->phi_b,
        *u_x = ac->u_x[t], ***phi_xx = ac->phi_xx[c][p][t], *u_xx = ac->u_xx[t], ***phi_bx = ac->phi_bx;

    rbf_gauss_valgrad(x, c, p, t, mixed_derivatives, anns_instance);

    for (i = 0; i < nc; i++) {
        if (phi[i] != 0) {
            for (j = 0; j < dim; j++) {
                phi_xx[i][j][j] = -phi_x[i][j]*xc[i][j]-a[i][0]*phi[i];
                for (k = 0; k < j; k++) {
                    phi_xx[i][j][k] = -xc[i][k] * phi_x[i][j];
                    phi_xx[i][k][j] = phi_xx[i][j][k];
                }
            }
        }
    }

    u = ac->u[t];

    for (j = 0; j < dim; j++) {
        p_x = 0.;
        for (i = 0; i < nc; i++) {
            if (phi[i] != 0) {
                p_x += w[i] * phi_x[i][j];
            }
        }
        u_x[j] = p_x;
    }

    // Вычисление градиента
    for (i = 0, s = 0; i < nc; i++, s+= nsize) {
        if (phi[i] != 0) {
            for (j = 0; j < dim; j++) {
                u_ax[s][j] = phi_x[i][j];

                phi_bx[i][j][0] = phi_x[i][j] * (r[i]*a[i][0] - 2) / b[i][0];
                u_ax[s+1][j] = phi_bx[i][j][0]*w[i];

                for (m = 0; m < dim; m++) {
                    u_ax[s+2+m][j] = -phi_xx[i][j][m]*w[i];
                }
            }
        } else {
            for (j = 0; j < dim; j++) {
                u_ax[s][j] = 0;
                u_ax[s+1][j] = 0;
                for (m = 0; m < dim; m++) {
                    u_ax[s+2+m][j] = 0;
                }
            }
        }
    }
}

void rbf_gauss_valgrad_d2(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *an = ac->nets[t];
    int i, j, m, s, k, dim = an->dim, nc = an->nc, nsize = an->nsize;
    double u, p_xx, phi_bxx, phi_xxx;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi[c][p][t],
        *u_a = ac->u_a[t], **u_ax = ac->u_ax[t], **phi_x = ac->phi_x[c][p][t], **phi_b = ac->phi_b,
        *u_x = ac->u_x[t], ***phi_xx = ac->phi_xx[c][p][t], *u_xx = ac->u_xx[t],
        **u_xy = ac->u_xy[t], ***phi_bx = ac->phi_bx, **u_axx = ac->u_axx[t], ***u_axy = ac->u_axy[t];

    rbf_gauss_valgrad_d1(x, c, p, t, mixed_derivatives, anns_instance);

    u = ac->u[t];
    for (j = 0; j < dim; j++) {
        p_xx = 0;
        for (i = 0; i < nc; i++) {
            if (phi[i] != 0) {
                p_xx += w[i] * phi_xx[i][j][j];
            }
        }
        u_xx[j] = p_xx;

        if (mixed_derivatives) {
            for (k = 0; k < j; k++) {
                p_xx = 0;
                for (i = 0; i < nc; i++) {
                    if (phi[i] != 0) {
                        p_xx += w[i] * phi_xx[i][j][k];
                    }
                }
                u_xy[j][k] = p_xx;
            }
        }
    }

    // Вычисление градиента
    for (i = 0, s = 0; i < nc; i++, s += nsize) {
        if (phi[i] != 0) {
            for (j = 0; j < dim; j++) {
                u_axx[s][j] = phi_xx[i][j][j];

                phi_bxx = phi_bx[i][j][0]*phi_x[i][j]/phi[i] - phi_b[i][0]*a[i][0] - 2*phi_xx[i][j][j]/b[i][0];

                u_axx[s+1][j] = phi_bxx*w[i];

                for (m = 0; m < dim; m++) {
                    phi_xxx = (-2*phi_xx[i][j][m]*phi_x[i][j] + phi_x[i][m]*phi_xx[i][j][j]) / phi[i] + 2*phi_x[i][m]*a[i][0];

                    u_axx[s+2+m][j] = phi_xxx*w[i];
                }

                if (mixed_derivatives) {
                    for (k = 0; k < j; k++) {
                        u_axy[s][j][k] = phi_xx[i][j][k];

                        phi_bxx = (phi_bx[i][j][0]*phi_x[i][k] + phi_x[i][j]*phi_bx[i][k][0] - phi_b[i][0]*phi_xx[i][j][k]) / phi[i];

                        u_axy[s+1][j][k] = phi_bxx*w[i];

                        for (m = 0; m < dim; m++) {
                            phi_xxx = (-phi_xx[i][m][j]*phi_x[i][k] - phi_xx[i][m][k]*phi_x[i][j] + phi_x[i][m]*phi_xx[i][j][k] ) / phi[i];

                            u_axy[s+2+m][j][k] = phi_xxx*w[i];
                        }
                    }
                }
            }
        } else {
            for (j = 0; j < dim; j++) {
                u_axx[s][j] = 0;
                u_axx[s+1][j] = 0;
                for (m = 0; m < dim; m++) {
                    u_axx[s+2+m][j] = 0;
                }

                if (mixed_derivatives) {
                    for (k = 0; k < j; k++) {
                        u_axy[s][j][k] = 0;
                        u_axy[s+1][j][k] = 0;
                        for (m = 0; m < ai->dim; m++) {
                            u_axy[s+2+m][j][k] = 0;
                        }
                    }
                }
            }
        }
    }
    return;
}

void rbf_gauss_val_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, I, nc = net->nc;
    double *w = ac->w, *phi = ac->phi[c][p][t], p_ = 0;

    for (i = 0; i < nc; i++) {
        p_ += x[i] * phi[i];
    }

    ac->u[t] = p_;
}


void rbf_gauss_val_d1_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, j, dim = net->dim, nc = net->nc;
    double u, p_x;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi[c][p][t],
        **phi_x = ac->phi_x[c][p][t], **phi_b = ac->phi_b, *u_x = ac->u_x[t];

    rbf_gauss_val_wonly(x, c, p, t, mixed_derivatives, anns_instance);

    u = ac->u[t];

    for (j = 0; j < dim; j++) {
        p_x = 0.;
        for (i = 0; i < nc; i++) {
            if (phi[i] != 0) {
                p_x += x[i] * phi_x[i][j];
            }
        }
        u_x[j] = p_x;
    }
}

void rbf_gauss_val_d2_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, j, k, dim = net->dim, nc = net->nc;
    double u, p_x, *r = ac->r, **a = ac->a, *w = ac->w, **xc = ac->xc, *phi = ac->phi[c][p][t], p_xx,
        **phi_x = ac->phi_x[c][p][t], *u_x = ac->u_x[t], ***phi_xx = ac->phi_xx[c][p][t],
        *u_xx = ac->u_xx[t], **u_xy = ac->u_xy[t];

    rbf_gauss_val_d1_wonly(x, c, p, t, mixed_derivatives, anns_instance);

    u = ac->u[t];

    for (j = 0; j < dim; j++) {
        p_xx = 0;
        for (i = 0; i < nc; i++) {
            if (phi[i] != 0) {
                p_xx += x[i] * phi_xx[i][j][j];
            }
        }
        u_xx[j] = p_xx;

        if (mixed_derivatives) {
            for (k = 0; k < j; k++) {
                p_xx = 0;
                for (i = 0; i < nc; i++) {
                    if (phi[i] != 0) {
                        p_xx += x[i] * phi_xx[i][j][k];
                    }
                }
                u_xy[j][k] = p_xx;
            }
        }
    }
}

void rbf_gauss_valgrad_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, nc = ai->nets[t]->nc;
    double q, *phi = ac->phi[c][p][t], *u_a = ac->u_a[t];

    rbf_gauss_val_wonly(x, c, p, t, mixed_derivatives, anns_instance);

    for (i = 0; i < nc; i++) {
        if (phi[i] != 0) {
            u_a[i] = phi[i];
        } else {
            u_a[i] = 0;
        }
    }
}


void rbf_gauss_valgrad_d1_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, j, nc = net->nc, dim = net->dim;
    double *phi = ac->phi[c][p][t], **phi_x = ac->phi_x[c][p][t], *u_a = ac->u_a[t], **u_ax = ac->u_ax[t];

    rbf_gauss_val_d1_wonly(x, c, p, t, mixed_derivatives, anns_instance);

    for (i = 0; i < nc; i++) {
        if (phi[i] != 0) {
            u_a[i] = phi[i];
            for (j = 0; j < dim; j++) {
                u_ax[i][j] = phi_x[i][j];
            }
        } else {
            u_a[i] = 0;
            for (j = 0; j < dim; j++) {
                u_ax[i][j] = 0;
            }
        }
    }
}


void rbf_gauss_valgrad_d2_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, j, k, nc = net->nc, dim = net->dim;
    double q = ac->q[c][p][t], *phi = ac->phi[c][p][t], **phi_x = ac->phi_x[c][p][t], *u_a = ac->u_a[t], **u_ax = ac->u_ax[t],
        **u_axx = ac->u_axx[t], ***u_axy = ac->u_axy[t], ***phi_xx = ac->phi_xx[c][p][t];

    rbf_gauss_val_d2_wonly(x, c, p, t, mixed_derivatives, anns_instance);

    for (i = 0; i < nc; i++) {
        if (phi[i] != 0) {
            u_a[i] = phi[i];
            for (j = 0; j < dim; j++) {
                u_ax[i][j] = phi_x[i][j];
                u_axx[i][j] = phi_xx[i][j][j];
                if (mixed_derivatives) {
                    for (k = 0; k < j; k++) {
                        u_axy[i][j][k] = phi_xx[i][j][k];
                    }
                }
            }
        } else {
            u_a[i] = 0;
            for (j = 0; j < dim; j++) {
                u_ax[i][j] = 0;
                u_axx[i][j] = 0;
                if (mixed_derivatives) {
                    for (k = 0; k < j; k++) {
                        u_axy[i][j][k] = 0;
                    }
                }
            }
        }
    }
}

static anns_func_t val_functions[] = {rbf_gauss_val, rbf_gauss_val_d1, rbf_gauss_val_d2, };
static anns_func_t valgrad_functions[] = {rbf_gauss_valgrad, rbf_gauss_valgrad_d1, rbf_gauss_valgrad_d2, };
static anns_func_t val_functions_wonly[] = {rbf_gauss_val_wonly, rbf_gauss_val_d1_wonly, rbf_gauss_val_d2_wonly, };
static anns_func_t valgrad_functions_wonly[] = {rbf_gauss_valgrad_wonly, rbf_gauss_valgrad_d1_wonly, rbf_gauss_valgrad_d2_wonly, };



void rbf_gauss_loadfuncs(anns_net_t *an) {
    if (!atexp) {
        atexp = at_new(-40, 1, 1e-3, exp);
    }

    an->ann_eval = rbf_gauss_eval;
    an->ann_val = val_functions;
    an->ann_valgrad = valgrad_functions;
    an->ann_val_wonly = val_functions_wonly;
    an->ann_valgrad_wonly = valgrad_functions_wonly;
}

