#include "nrbf_gauss.h"
#include "../anns_tablecalc.h"

static anns_tablecalc_t *atexp = NULL;
#define MIN_EPOW (-40.0)

double nrbf_gauss_eval(double *z, double *y, int t, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_net_t *net = ai->nets[t];
    int i, j, I, nc = net->nc, from = net->I, nsize = net->nsize;
    double p = 0, q = 0, phi, epow, r, w, b, a, xc;

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
            q += phi;
            p += w * phi;
        }
    }

    return p / q;
}

void nrbf_gauss_val(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, j, I, from = net->I, nsize = net->nsize, nc = net->nc;
    double p_ = 0., q = 0, lambda, min_lambda = 1e100, epow, *point = ac->points[c][p];
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc;
    double *phi = ac->phi[c][p][t];

    for (i = 0, I = from; i < nc; i++, I += nsize) {
        r[i] = 0.;
        w[i] = x[I];
        b[i][0] = x[I + 1];
        a[i][0] = 1 / (b[i][0]*b[i][0]);

        for (j = 0; j < net->dim; j++) {
            xc[i][j] = point[j] - x[I + 2 + j];
            r[i] += xc[i][j] * xc[i][j];
            xc[i][j] *= a[i][0];
        }

        // Переменная для нормировки - минимальное значение для экспоненты
        lambda = r[i] * a[i][0];
        if (lambda < min_lambda) {
            min_lambda = lambda;
        }
    }

    for (i = 0, I = from; i < nc; i++, I += nsize) {
        epow = -0.5 * (r[i] * a[i][0] - min_lambda);
        if (epow > MIN_EPOW) {
            phi[i] = at_fasteval(atexp, epow);
//            phi[i] = exp(epow);

            q += phi[i];
            p_ += w[i] * phi[i];
        } else {
            phi[i] = 0;
        }
    }

    ac->u[t] = p_ / q;
    ac->q[c][p][t] = q;
}

void nrbf_gauss_val_d1(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, j, dim = net->dim, nc = net->nc;
    double q, u, p_x;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi[c][p][t],
        *q_x = ac->q_x[c][p][t], **phi_x = ac->phi_x[c][p][t], **phi_b = ac->phi_b, *u_x = ac->u_x[t];

    nrbf_gauss_val(x, c, p, t, mixed_derivatives, anns_instance);

    u = ac->u[t];
    q = ac->q[c][p][t];

    for (j = 0; j < dim; j++) {
        p_x = 0.;
        q_x[j] = 0.;
        for (i = 0; i < nc; i++) {
            if (phi[i] != 0) {
                phi_x[i][j] = -xc[i][j] * phi[i];
                q_x[j] += phi_x[i][j];
                p_x += w[i] * phi_x[i][j];
            }
        }
        u_x[j] = (p_x - u*q_x[j]) / q;
    }
}

void nrbf_gauss_val_d2(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, j, k, dim = net->dim, nc = net->nc;
    double q, u, p_xx;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi[c][p][t],
        *q_x = ac->q_x[c][p][t], **phi_x = ac->phi_x[c][p][t], **phi_b = ac->phi_b, *u_x = ac->u_x[t],
        ***phi_xx = ac->phi_xx[c][p][t], **q_xx = ac->q_xx[c][p][t], *u_xx = ac->u_xx[t], **u_xy = ac->u_xy[t];

    nrbf_gauss_val_d1(x, c, p, t, mixed_derivatives, anns_instance);

    u = ac->u[t];
    q = ac->q[c][p][t];

    for (j = 0; j < dim; j++) {
        p_xx = 0;
        q_xx[j][j] = 0;
        for (i = 0; i < nc; i++) {
            if (phi[i] != 0) {
                phi_xx[i][j][j] = -phi_x[i][j]*xc[i][j] - a[i][0]*phi[i];
                q_xx[j][j] += phi_xx[i][j][j];
                p_xx += w[i] * phi_xx[i][j][j];
            }
        }
        u_xx[j] = (p_xx - 2*q_x[j]*u_x[j] - u*q_xx[j][j]) / q;

        if (mixed_derivatives) {
            for (k = 0; k < j; k++) {
                p_xx = 0;
                q_xx[j][k] = 0;
                for (i = 0; i < nc; i++) {
                    if (phi[i] != 0) {
                        phi_xx[i][j][k] = xc[i][j] * xc[i][k] * phi[i];
                        q_xx[j][k] += phi_xx[i][j][k];
                        p_xx += w[i] * phi_xx[i][j][k];
                    }
                }
                u_xy[j][k] = (p_xx - q_x[j]*u_x[k] - q_x[k]*u_x[j] - u*q_xx[j][k]) / q;
            }
        }
    }
}


void nrbf_gauss_valgrad(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, j, m, I, s, from = net->I, nsize = net->nsize, nc = net->nc, dim = net->dim;
    double p_ = 0., q = 0., u, lambda, min_lambda = 1e100, epow, *point = ac->points[c][p];
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi[c][p][t],
        *u_a = ac->u_a[t], *q_x = ac->q_x[c][p][t], **phi_x = ac->phi_x[c][p][t], **phi_b = ac->phi_b;

    for (i = 0, I = from; i < nc; i++, I += nsize) {
        r[i] = 0.;
        w[i] = x[I];
        b[i][0] = x[I + 1];
        a[i][0] = 1 / (b[i][0]*b[i][0]);

        for (j = 0; j < dim; j++) {
            xc[i][j] = point[j] - x[I + 2 + j];
            r[i] += xc[i][j] * xc[i][j];
            xc[i][j] *= a[i][0];
        }

        // Переменная для нормировки - минимальное значение для экспоненты
        lambda = r[i] * a[i][0];
        if (lambda < min_lambda) {
            min_lambda = lambda;
        }
    }

    for (i = 0, I = from; i < nc; i++, I += nsize) {
        epow = -0.5 * (r[i] * a[i][0] - min_lambda);
        if (epow > MIN_EPOW) {
            phi[i] = at_fasteval(atexp, epow);
//            phi[i] = exp(epow);

            q += phi[i];
            p_ += w[i] * phi[i];

            for (j = 0; j < dim; j++) {
                phi_x[i][j] = -xc[i][j] * phi[i];
            }
        } else {
            phi[i] = 0;
        }
    }

    u =  p_ / q;
    ac->u[t] = u;
    ac->q[c][p][t] = q;

    // Вычисление градиента
    for (i = 0, s = 0; i < nc; i++, s += nsize) {
        if (phi[i] != 0) {
            u_a[s] = phi[i] / q;

            phi_b[i][0] = phi[i] * r[i] * a[i][0] / b[i][0];
            u_a[s+1] = phi_b[i][0] * (w[i] - u) / q;

            for (m = 0; m < dim; m++) {
                u_a[s+2+m] = -phi_x[i][m] * (w[i] - u) / q;
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

void nrbf_gauss_valgrad_d1(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, j, m, s, k, nsize = net->nsize, nc = net->nc, dim = net->dim;
    double q = 0., u, p_x;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi[c][p][t],
        *u_a = ac->u_a[t], **u_ax = ac->u_ax[t], *q_x = ac->q_x[c][p][t], **phi_x = ac->phi_x[c][p][t], **phi_b = ac->phi_b,
        *u_x = ac->u_x[t], ***phi_xx = ac->phi_xx[c][p][t], *u_xx = ac->u_xx[t], ***phi_bx = ac->phi_bx;

    nrbf_gauss_valgrad(x, c, p, t, mixed_derivatives, anns_instance);

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
    q = ac->q[c][p][t];

    for (j = 0; j < dim; j++) {
        p_x = 0.;
        q_x[j] = 0.;
        for (i = 0; i < nc; i++) {
            if (phi[i] != 0) {
                q_x[j] += phi_x[i][j];
                p_x += w[i] * phi_x[i][j];
            }
        }
        u_x[j] = (p_x - u*q_x[j]) / q;
    }

    // Вычисление градиента
    for (i = 0, s = 0; i < nc; i++, s += nsize) {
        if (phi[i] != 0) {
             for (j = 0; j < dim; j++) {
                u_ax[s][j] = (phi_x[i][j] - u_a[s] * q_x[j]) / q;

                phi_bx[i][j][0] = phi_x[i][j] * (r[i]*a[i][0] - 2) / b[i][0];
                u_ax[s+1][j] = (phi_bx[i][j][0]*(w[i]-u) - u_a[s+1]*q_x[j] - u_x[j]*phi_b[i][0]) / q;

                 for (m = 0; m < dim; m++) {
                    u_ax[s+2+m][j] = (-phi_xx[i][j][m]*(w[i] - u) + phi_x[i][m]*u_x[j] - u_a[s+2+m]*q_x[j]) / q;
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

void nrbf_gauss_valgrad_d2(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, j, m, s, k, nsize = net->nsize, nc = net->nc, dim = net->dim;
    double q = 0., u, p_xx, phi_bxx, phi_xxx;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi[c][p][t],
        *u_a = ac->u_a[t], **u_ax = ac->u_ax[t], *q_x = ac->q_x[c][p][t], **phi_x = ac->phi_x[c][p][t], **phi_b = ac->phi_b,
        *u_x = ac->u_x[t], ***phi_xx = ac->phi_xx[c][p][t], **q_xx = ac->q_xx[c][p][t], *u_xx = ac->u_xx[t],
        **u_xy = ac->u_xy[t], ***phi_bx = ac->phi_bx, **u_axx = ac->u_axx[t], ***u_axy = ac->u_axy[t];

    nrbf_gauss_valgrad_d1(x, c, p, t, mixed_derivatives, anns_instance);

    u = ac->u[t];
    q = ac->q[c][p][t];
    for (j = 0; j < dim; j++) {
        p_xx = 0;
        q_xx[j][j] = 0;
        for (i = 0; i < nc; i++) {
            if (phi[i] != 0) {
                q_xx[j][j] += phi_xx[i][j][j];
                p_xx += w[i] * phi_xx[i][j][j];
            }
        }
        u_xx[j] = (p_xx - 2*q_x[j]*u_x[j] - u*q_xx[j][j]) / q;

        if (mixed_derivatives) {
            for (k = 0; k < j; k++) {
                p_xx = 0;
                q_xx[j][k] = 0;
                for (i = 0; i < nc; i++) {
                    if (phi[i] != 0) {
                        q_xx[j][k] += phi_xx[i][j][k];
                        p_xx += w[i] * phi_xx[i][j][k];
                    }
                }
                u_xy[j][k] = (p_xx - q_x[j]*u_x[k] - q_x[k]*u_x[j] - u*q_xx[j][k]) / q;
            }
        }
    }

    // Вычисление градиента
    for (i = 0, s = 0; i < nc; i++, s += nsize) {
        if (phi[i] != 0) {
            for (j = 0; j < dim; j++) {
                u_axx[s][j] = (phi_xx[i][j][j] - 2*u_ax[s][j]*q_x[j] - u_a[s]*q_xx[j][j]) / q;

                phi_bxx = phi_bx[i][j][0]*phi_x[i][j]/phi[i] - phi_b[i][0]*a[i][0] - 2*phi_xx[i][j][j]/b[i][0];

                if (phi_bxx != phi_bxx) {phi_bxx = 0;}

                u_axx[s+1][j] = (phi_bxx*(w[i]-u) - 2*u_ax[s+1][j]*q_x[j] - 2*u_x[j]*phi_bx[i][j][0] - u_a[s+1]*q_xx[j][j] - u_xx[j]*phi_b[i][0])/q;

                for (m = 0; m < dim; m++) {
                    phi_xxx = (-2*phi_xx[i][j][m]*phi_x[i][j] + phi_x[i][m]*phi_xx[i][j][j]) / phi[i] + 2*phi_x[i][m]*a[i][0];

                    if (phi_xxx != phi_xxx) {phi_xxx = 0;}

                    u_axx[s+2+m][j] = (phi_xxx*(w[i]-u) - 2*u_ax[s+2+m][j]*q_x[j] - 2*u_x[j]*(-phi_xx[i][m][j]) - u_xx[j]*(-phi_x[i][m]) - u_a[s+2+m]*q_xx[j][j]) / q;
                }

                if (mixed_derivatives) {
                    for (k = 0; k < j; k++) {
                        u_axy[s][j][k] = (phi_xx[i][j][k] - u_ax[s][j]*q_x[k] - u_ax[s][k]*q_x[j] - u_a[s]*q_xx[j][k]) / q;

                        phi_bxx = (phi_bx[i][j][0]*phi_x[i][k] + phi_x[i][j]*phi_bx[i][k][0] - phi_b[i][0]*phi_xx[i][j][k]) / phi[i];

                        if (phi_bxx != phi_bxx) {phi_bxx = 0;}

                        u_axy[s+1][j][k] = (phi_bxx*(w[i]-u) - u_ax[s+1][j]*q_x[k] - u_ax[s+1][k]*q_x[j] - u_xy[j][k]*phi_b[i][0] - u_x[j]*phi_bx[i][k][0] - u_x[k]*phi_bx[i][j][0] - u_a[s+1]*q_xx[j][k])/q;

                        for (m = 0; m < dim; m++) {
                            phi_xxx = (-phi_xx[i][m][j]*phi_x[i][k] - phi_xx[i][m][k]*phi_x[i][j] + phi_x[i][m]*phi_xx[i][j][k] ) / phi[i];

                            if (phi_xxx != phi_xxx) {phi_xxx = 0;}

                            u_axy[s+2+m][j][k] = (phi_xxx*(w[i]-u) - u_ax[s+2+m][j]*q_x[k] - u_ax[s+2+m][k]*q_x[j] - u_xy[j][k]*(-phi_x[i][m]) - u_x[j]*(-phi_xx[i][m][k]) - u_x[k]*(-phi_xx[i][m][j]) - u_a[s+2+m]*q_xx[j][k] ) / q;
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
                        for (m = 0; m < dim; m++) {
                            u_axy[s+2+m][j][k] = 0;
                        }
                    }
                }
            }
        }
    }
    return;
}

void nrbf_gauss_val_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, I, nc = net->nc;
    double *w = ac->w, *phi = ac->phi[c][p][t], p_ = 0;

    for (i = 0; i < nc; i++) {
//        p_ += w[i] * phi[i];
        p_ += x[i] * phi[i];
    }

    ac->u[t] = p_ / ac->q[c][p][t];
}

void nrbf_gauss_val_d1_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, j, dim = net->dim, nc = net->nc;
    double q = 0., u, p_x;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi[c][p][t],
        *q_x = ac->q_x[c][p][t], **phi_x = ac->phi_x[c][p][t], **phi_b = ac->phi_b, *u_x = ac->u_x[t];

    nrbf_gauss_val_wonly(x, c, p, t, mixed_derivatives, anns_instance);

    u = ac->u[t];
    q = ac->q[c][p][t];

    for (j = 0; j < dim; j++) {
        p_x = 0.;
        for (i = 0; i < nc; i++) {
            if (phi[i] != 0) {
                p_x+= x[i] * phi_x[i][j];
//                p_x += w[i] * phi_x[i][j];
            }
        }
        u_x[j] = (p_x - u*q_x[j]) / q;
    }
}

void nrbf_gauss_val_d2_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, j, k, dim = net->dim, nc = net->nc;
    double u, q, p_x, *r = ac->r, **a = ac->a, *w = ac->w, **xc = ac->xc, *phi = ac->phi[c][p][t], p_xx,
        *q_x = ac->q_x[c][p][t], **phi_x = ac->phi_x[c][p][t], *u_x = ac->u_x[t], ***phi_xx = ac->phi_xx[c][p][t],
        *u_xx = ac->u_xx[t], **q_xx = ac->q_xx[c][p][t], **u_xy = ac->u_xy[t];

    nrbf_gauss_val_d1_wonly(x, c, p, t, mixed_derivatives, anns_instance);

    u = ac->u[t];
    q = ac->q[c][p][t];

    for (j = 0; j < dim; j++) {
        p_xx = 0;
        for (i = 0; i < nc; i++) {
            if (phi[i] != 0) {
//                p_xx += w[i] * phi_xx[i][j][j];
                p_xx += x[i] * phi_xx[i][j][j];
            }
        }
        u_xx[j] = (p_xx - 2*q_x[j]*u_x[j] - u*q_xx[j][j]) / q;

        if (mixed_derivatives) {
            for (k = 0; k < j; k++) {
                p_xx = 0;
                for (i = 0; i < nc; i++) {
                    if (phi[i] != 0) {
//                        p_xx += w[i] * phi_xx[i][j][k];
                        p_xx += x[i] * phi_xx[i][j][k];
                    }
                }
                u_xy[j][k] = (p_xx - q_x[j]*u_x[k] - q_x[k]*u_x[j] - u*q_xx[j][k]) / q;
            }
        }
    }
}

void nrbf_gauss_valgrad_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, nc = ai->nets[t]->nc;
    double q, *phi = ac->phi[c][p][t], *u_a = ac->u_a[t];

    nrbf_gauss_val_wonly(x, c, p, t, mixed_derivatives, anns_instance);

    q = ac->q[c][p][t];
    for (i = 0; i < nc; i++) {
        if (phi[i] != 0) {
            u_a[i] = phi[i] / q;
        } else {
            u_a[i] = 0;
        }
    }

}

void nrbf_gauss_valgrad_d1_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, j, nc = net->nc, dim = net->dim;
    double q = ac->q[c][p][t], *phi = ac->phi[c][p][t], **phi_x = ac->phi_x[c][p][t], *u_a = ac->u_a[t], **u_ax = ac->u_ax[t],
        *q_x = ac->q_x[c][p][t];

    nrbf_gauss_val_d1_wonly(x, c, p, t, mixed_derivatives, anns_instance);

    for (i = 0; i < nc; i++) {
        if (phi[i] != 0) {
            u_a[i] = phi[i] / q;
            for (j = 0; j < dim; j++) {
                u_ax[i][j] = (phi_x[i][j] - u_a[i] * q_x[j]) / q;
            }
        } else {
            u_a[i] = 0;
            for (j = 0; j < dim; j++) {
                u_ax[i][j] = 0;
            }
        }
    }
}

void nrbf_gauss_valgrad_d2_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    anns_net_t *net = ai->nets[t];
    int i, j, k, nc = net->nc, dim = net->dim;
    double q = ac->q[c][p][t], *phi = ac->phi[c][p][t], **phi_x = ac->phi_x[c][p][t], *u_a = ac->u_a[t], **u_ax = ac->u_ax[t],
        *q_x = ac->q_x[c][p][t], **u_axx = ac->u_axx[t], ***u_axy = ac->u_axy[t], ***phi_xx = ac->phi_xx[c][p][t],
        **q_xx = ac->q_xx[c][p][t];

    nrbf_gauss_val_d2_wonly(x, c, p, t, mixed_derivatives, anns_instance);

    for (i = 0; i < nc; i++) {
        if (phi[i] != 0) {
            u_a[i] = phi[i] / q;
            for (j = 0; j < dim; j++) {
                u_ax[i][j] = (phi_x[i][j] - u_a[i] * q_x[j]) / q;
                u_axx[i][j] = (phi_xx[i][j][j] - 2*u_ax[i][j]*q_x[j] - u_a[i]*q_xx[j][j]) / q;
                if (mixed_derivatives) {
                    for (k = 0; k < j; k++) {
                        u_axy[i][j][k] = (phi_xx[i][j][k] - u_ax[i][j]*q_x[k] - u_ax[i][k]*q_x[j] - u_a[i]*q_xx[j][k]) / q;
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


static anns_func_t val_functions[] = {nrbf_gauss_val, nrbf_gauss_val_d1, nrbf_gauss_val_d2, };
static anns_func_t valgrad_functions[] = {nrbf_gauss_valgrad, nrbf_gauss_valgrad_d1, nrbf_gauss_valgrad_d2, };
static anns_func_t val_wonly_functions[] = {nrbf_gauss_val_wonly, nrbf_gauss_val_d1_wonly, nrbf_gauss_val_d2_wonly, };
static anns_func_t valgrad_wonly_functions[] = {nrbf_gauss_valgrad_wonly, nrbf_gauss_valgrad_d1_wonly, nrbf_gauss_valgrad_d2_wonly, };


void nrbf_gauss_loadfuncs(anns_net_t *net) {
    if (!atexp) {
        atexp = at_new(-40, 1, 1e-3, exp);
    }

    net->ann_eval = nrbf_gauss_eval;
    net->ann_val = val_functions;
    net->ann_valgrad = valgrad_functions;
    net->ann_val_wonly = val_wonly_functions;
    net->ann_valgrad_wonly = valgrad_wonly_functions;
}
