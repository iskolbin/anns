#include "rbf_cauchy.h"

void rbf_cauchy_val(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, j, I;
    double p = 0.;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi;

    for (i = 0, I = ai->I[t]; i < ai->nn[t]; i++, I += ai->S[t]) {
        r[i] = 0.;
        w[i] = x[I];
        b[i][0] = x[I + 1];
        a[i][0] = 1 / (b[i][0]*b[i][0]);

        for (j = 0; j < ai->nd[t]; j++) {
            xc[i][j] = point[j] - x[I + 2 + j];
            r[i] += xc[i][j] * xc[i][j];
            xc[i][j] *= a[i][0];
        }
    }

    for (i = 0, I = ai->I[t]; i < ai->nn[t]; i++, I += ai->S[t]) {
        phi[i] = 1/(1.0+0.5*r[i]*a[i][0]);
        if (phi[i] != 0) {
            p += w[i] * phi[i];
        }
    }

    ac->u[t] = p;
}

void rbf_cauchy_val_d1(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, j;
    double p = 0.,  u, p_x;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi,
        **phi_x = ac->phi_x, **phi_b = ac->phi_b, *u_x = ac->u_x[t];

    rbf_cauchy_val(point, x, t, mixed_derivatives, anns_instance);

    u = ac->u[t];

    for (j = 0; j < ai->nd[t]; j++) {
        p_x = 0.;
        for (i = 0; i < ai->nn[t]; i++) {
            if (phi[i] != 0) {
                phi_x[i][j] = -xc[i][j]*phi[i]*phi[i];
                p_x += w[i] * phi_x[i][j];
            }
        }
        u_x[j] = p_x;
    }
}

void rbf_cauchy_val_d2(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, j, k;
    double p = 0., u, p_xx;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi,
        **phi_x = ac->phi_x, **phi_b = ac->phi_b, *u_x = ac->u_x[t], ***phi_xx = ac->phi_xx,
        *u_xx = ac->u_xx[t], **u_xy = ac->u_xy[t];

    rbf_cauchy_val_d1(point, x, t, mixed_derivatives, anns_instance);

    u = ac->u[t];

    for (j = 0; j < ai->nd[t]; j++) {
        p_xx = 0;
        for (i = 0; i < ai->nn[t]; i++) {
            if (phi[i] != 0) {
                phi_xx[i][j][j] = -2*xc[i][j]*phi[i]*phi_x[i][j] - a[i][0]*phi[i]*phi[i];
                p_xx += w[i] * phi_xx[i][j][j];
            }
        }
        u_xx[j] = p_xx;

        if (mixed_derivatives) {
            for (k = 0; k < j; k++) {
                p_xx = 0;
                for (i = 0; i < ai->nn[t]; i++) {
                    if (phi[i] != 0) {
                        phi_xx[i][j][k] = -2*xc[i][j]*phi[i]*phi_x[i][k];
                        p_xx += w[i] * phi_xx[i][j][k];
                    }
                }
                u_xy[j][k] = p_xx;
            }
        }
    }
}


void rbf_cauchy_valgrad(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, j, m, I, s;
    double p = 0., u;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi,
        *u_a = ac->u_a[t], *q_x = ac->q_x, **phi_x = ac->phi_x, **phi_b = ac->phi_b;


    for (i = 0, I = ai->I[t]; i < ai->nn[t]; i++, I += ai->S[t]) {
        r[i] = 0.;
        w[i] = x[I];
        b[i][0] = x[I + 1];
        a[i][0] = 1 / (b[i][0]*b[i][0]);

        for (j = 0; j < ai->nd[t]; j++) {
            xc[i][j] = point[j] - x[I + 2 + j];
            r[i] += xc[i][j] * xc[i][j];
            xc[i][j] *= a[i][0];
        }
    }

    for (i = 0, I = ai->I[t]; i < ai->nn[t]; i++, I += ai->S[t]) {
        phi[i] = 1/(1.0+0.5*r[i]*a[i][0]);

        if (phi[i] != 0) {
            p += w[i] * phi[i];

            for (j = 0; j < ai->nd[t]; j++) {
                phi_x[i][j] = -xc[i][j]*phi[i]*phi[i];
            }
        }
    }

    u = p;
    ac->u[t] = u;

    // Вычисление градиента
    for (i = 0, s = 0; i < ai->nn[t]; i++, s+= ai->S[t]) {
        if (phi[i] != 0) {
            u_a[s] = phi[i];

            phi_b[i][0] = phi[i] * phi[i] * r[i] * a[i][0] / b[i][0];
            u_a[s+1] = phi_b[i][0] * w[i];

            for (m = 0; m < ai->nd[t]; m++) {
                u_a[s+2+m] = -phi_x[i][m] * w[i];
            }
        } else {
            u_a[s] = 0;
            u_a[s+1] = 0;
            for (m = 0; m < ai->nd[t]; m++) {
                u_a[s+2+m] = 0;
            }
        }
    }
}

void rbf_cauchy_valgrad_d1(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, j, m, s, k;
    double p = 0., u, p_x;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi,
        *u_a = ac->u_a[t], **u_ax = ac->u_ax[t], **phi_x = ac->phi_x, **phi_b = ac->phi_b,
        *u_x = ac->u_x[t], ***phi_xx = ac->phi_xx, *u_xx = ac->u_xx[t], ***phi_bx = ac->phi_bx;

    rbf_cauchy_valgrad(point, x, t, mixed_derivatives, anns_instance);

    for (i = 0; i < ai->nn[t]; i++) {
        if (phi[i] != 0) {
            for (j = 0; j < ai->nd[t]; j++) {
                phi_xx[i][j][j] = -2*xc[i][j]*phi[i]*phi_x[i][j] - a[i][0]*phi[i]*phi[i];
                for (k = 0; k < j; k++) {
                    phi_xx[i][j][k] = -2*xc[i][j]*phi[i]*phi_x[i][k];
                    phi_xx[i][k][j] = phi_xx[i][j][k];
                }
            }
        }
    }

    u = ac->u[t];

    for (j = 0; j < ai->nd[t]; j++) {
        p_x = 0.;
        for (i = 0; i < ai->nn[t]; i++) {
            if (phi[i] != 0) {
                p_x += w[i] * phi_x[i][j];
            }
        }
        u_x[j] = p_x;
    }

    // Вычисление градиента
    for (i = 0, s = 0; i < ai->nn[t]; i++, s+= ai->S[t]) {
        if (phi[i] != 0) {
            for (j = 0; j < ai->nd[t]; j++) {
                u_ax[s][j] = phi_x[i][j];

                phi_bx[i][j][0] = 2*phi_x[i][j]*(phi_b[i][0]/phi[i] - 1/b[i][0]);
                u_ax[s+1][j] = phi_bx[i][j][0]*w[i];

                for (m = 0; m < ai->nd[t]; m++) {
                    u_ax[s+2+m][j] = -phi_xx[i][j][m]*w[i];
                }
            }
        } else {
            for (j = 0; j < ai->nd[t]; j++) {
                u_ax[s][j] = 0;
                u_ax[s+1][j] = 0;
                for (m = 0; m < ai->nd[t]; m++) {
                    u_ax[s+2+m][j] = 0;
                }
            }
        }
    }
}

void rbf_cauchy_valgrad_d2(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, j, m, s, k;
    double p = 0., u, p_xx, phi_bxx, phi_xxx;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi,
        *u_a = ac->u_a[t], **u_ax = ac->u_ax[t], **phi_x = ac->phi_x, **phi_b = ac->phi_b,
        *u_x = ac->u_x[t], ***phi_xx = ac->phi_xx, *u_xx = ac->u_xx[t],
        **u_xy = ac->u_xy[t], ***phi_bx = ac->phi_bx, **u_axx = ac->u_axx[t], ***u_axy = ac->u_axy[t];

    rbf_cauchy_valgrad_d1(point, x, t, mixed_derivatives, anns_instance);

    u = ac->u[t];

    for (j = 0; j < ai->nd[t]; j++) {
        p_xx = 0;
        for (i = 0; i < ai->nn[t]; i++) {
            if (phi[i] != 0) {
                p_xx += w[i] * phi_xx[i][j][j];
            }
        }
        u_xx[j] = p_xx;

        if (mixed_derivatives) {
            for (k = 0; k < j; k++) {
                p_xx = 0;
                for (i = 0; i < ai->nn[t]; i++) {
                    if (phi[i] != 0) {
                        p_xx += w[i] * phi_xx[i][j][k];
                    }
                }
                u_xy[j][k] = p_xx;
            }
        }
    }

    // Вычисление градиента
    for (i = 0, s = 0; i < ai->nn[t]; i++, s+= ai->S[t]) {
        if (phi[i] != 0) {
            for (j = 0; j < ai->nd[t]; j++) {
                u_axx[s][j] = phi_xx[i][j][j];

                phi_bxx = phi_xx[i][j][j]*phi_bx[i][j][0] / phi_x[i][j] + 2*phi_x[i][j]*(phi_bx[i][j][0]/phi[i] + xc[i][j]*phi_b[i][0]);
                u_axx[s+1][j] = phi_bxx*w[i];

                for (m = 0; m < ai->nd[t]; m++) {
                    phi_xxx = -2 * (xc[i][j]*(phi_xx[i][m][j]*phi[i] + phi_x[i][m]*phi_x[i][j]) + phi[i]*phi_x[i][m]*a[i][0] + (m==j ? phi[i]*phi_x[i][j]*a[i][0] : 0) );
                    u_axx[s+2+m][j] = -phi_xxx*w[i];
                }

                if (mixed_derivatives) {
                    for (k = 0; k < j; k++) {
                        u_axy[s][j][k] = phi_xx[i][j][k];

                        phi_bxx = phi_xx[i][j][k]*phi_bx[i][j][0] / phi_x[i][j] + 2*phi_x[i][j]*(phi_bx[i][k][0]/phi[i] + xc[i][k]*phi_b[i][0]);
                        u_axy[s+1][j][k] = phi_bxx*w[i];

                        for (m = 0; m < ai->nd[t]; m++) {
                            phi_xxx = -2 * (xc[i][j]*(phi_xx[i][m][k]*phi[i] + phi_x[i][m]*phi_x[i][k]) + (m==j ? phi[i]*phi_x[i][k]*a[i][0] : 0) );
                            u_axy[s+2+m][j][k] = -phi_xxx*w[i];
                        }
                    }
                }
            }
        } else {
            for (j = 0; j < ai->nd[t]; j++) {
                u_axx[s][j] = 0;
                u_axx[s+1][j] = 0;
                for (m = 0; m < ai->nd[t]; m++) {
                    u_axx[s+2+m][j] = 0;
                }

                if (mixed_derivatives) {
                    for (k = 0; k < j; k++) {
                        u_axy[s][j][k] = 0;
                        u_axy[s+1][j][k] = 0;
                        for (m = 0; m < ai->nd[t]; m++) {
                            u_axy[s+2+m][j][k] = 0;
                        }
                    }
                }
            }
        }
    }
    return;
}

static anns_func_t val_functions[] = {rbf_cauchy_val, rbf_cauchy_val_d1, rbf_cauchy_val_d2, };
static anns_func_t valgrad_functions[] = {rbf_cauchy_valgrad, rbf_cauchy_valgrad_d1, rbf_cauchy_valgrad_d2, };

void rbf_cauchy_loadfuncs(int t, void *anns_instance) {
    anns_instance_t *ai = anns_instance;

    ai->ann_val[t] = val_functions;
    ai->ann_valgrad[t] = valgrad_functions;
}
