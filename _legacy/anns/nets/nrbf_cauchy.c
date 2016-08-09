#include "nrbf_cauchy.h"

void nrbf_cauchy_val(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, j, I;
    double p = 0., q = 0.;
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
            q += phi[i];
            p += w[i] * phi[i];
        }
    }

    ac->u[t] = p / q;
    ac->q = q;
}

void nrbf_cauchy_val_d1(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, j;
    double p = 0., q = 0., u, p_x;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi,
        *q_x = ac->q_x, **phi_x = ac->phi_x, **phi_b = ac->phi_b, *u_x = ac->u_x[t];

    nrbf_cauchy_val(point, x, t, mixed_derivatives, anns_instance);

    u = ac->u[t];
    q = ac->q;

    for (j = 0; j < ai->nd[t]; j++) {
        p_x = 0.;
        q_x[j] = 0.;
        for (i = 0; i < ai->nn[t]; i++) {
            if (phi[i] != 0) {
                phi_x[i][j] = -xc[i][j]*phi[i]*phi[i];
                q_x[j] += phi_x[i][j];
                p_x += w[i] * phi_x[i][j];
            }
        }
        u_x[j] = (p_x - u*q_x[j]) / q;
    }
}

void nrbf_cauchy_val_d2(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, j, k;
    double p = 0., q = 0., u, p_xx;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi,
        *q_x = ac->q_x, **phi_x = ac->phi_x, **phi_b = ac->phi_b, *u_x = ac->u_x[t], ***phi_xx = ac->phi_xx,
        **q_xx = ac->q_xx, *u_xx = ac->u_xx[t], **u_xy = ac->u_xy[t];

    nrbf_cauchy_val_d1(point, x, t, mixed_derivatives, anns_instance);

    u = ac->u[t];
    q = ac->q;

    for (j = 0; j < ai->nd[t]; j++) {
        p_xx = 0;
        q_xx[j][j] = 0;
        for (i = 0; i < ai->nn[t]; i++) {
            if (phi[i] != 0) {
                phi_xx[i][j][j] = -2*xc[i][j]*phi[i]*phi_x[i][j] - a[i][0]*phi[i]*phi[i];
                q_xx[j][j] += phi_xx[i][j][j];
                p_xx += w[i] * phi_xx[i][j][j];
            }
        }
        u_xx[j] = (p_xx - 2*q_x[j]*u_x[j] - u*q_xx[j][j]) / q;

        if (mixed_derivatives) {
            for (k = 0; k < j; k++) {
                p_xx = 0;
                q_xx[j][k] = 0;
                for (i = 0; i < ai->nn[t]; i++) {
                    if (phi[i] != 0) {
                        phi_xx[i][j][k] = -2*xc[i][j]*phi[i]*phi_x[i][k];
                        q_xx[j][k] += phi_xx[i][j][k];
                        p_xx += w[i] * phi_xx[i][j][k];
                    }
                }
                u_xy[j][k] = (p_xx - q_x[j]*u_x[k] - q_x[k]*u_x[j] - u*q_xx[j][k]) / q;
            }
        }
    }
}


void nrbf_cauchy_valgrad(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, j, m, I, s;
    double p = 0., q = 0., u;
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
            q += phi[i];
            p += w[i] * phi[i];

            for (j = 0; j < ai->nd[t]; j++) {
                phi_x[i][j] = -xc[i][j]*phi[i]*phi[i];
            }
        }
    }

    u =  p / q;
    ac->u[t] = u;
    ac->q = q;

    // Вычисление градиента
    for (i = 0, s = 0; i < ai->nn[t]; i++, s+= ai->S[t]) {
        if (phi[i] != 0) {
            u_a[s] = phi[i] / q;

            phi_b[i][0] = phi[i] * phi[i] * r[i] * a[i][0] / b[i][0];
            u_a[s+1] = phi_b[i][0] * (w[i] - u) / q;

            for (m = 0; m < ai->nd[t]; m++) {
                u_a[s+2+m] = -phi_x[i][m] * (w[i] - u) / q;
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

void nrbf_cauchy_valgrad_d1(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, j, m, s, k;
    double p = 0., q = 0., u, p_x;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi,
        *u_a = ac->u_a[t], **u_ax = ac->u_ax[t], *q_x = ac->q_x, **phi_x = ac->phi_x, **phi_b = ac->phi_b,
        *u_x = ac->u_x[t], ***phi_xx = ac->phi_xx, *u_xx = ac->u_xx[t], ***phi_bx = ac->phi_bx;

    nrbf_cauchy_valgrad(point, x, t, mixed_derivatives, anns_instance);

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
    q = ac->q;

    for (j = 0; j < ai->nd[t]; j++) {
        p_x = 0.;
        q_x[j] = 0.;
        for (i = 0; i < ai->nn[t]; i++) {
            if (phi[i] != 0) {
                q_x[j] += phi_x[i][j];
                p_x += w[i] * phi_x[i][j];
            }
        }
        u_x[j] = (p_x - u*q_x[j]) / q;
    }

    // Вычисление градиента
    for (i = 0, s = 0; i < ai->nn[t]; i++, s+= ai->S[t]) {
        if (phi[i] != 0) {
            for (j = 0; j < ai->nd[t]; j++) {
                u_ax[s][j] = (phi_x[i][j] - u_a[s] * q_x[j]) / q;

                phi_bx[i][j][0] = 2*phi_x[i][j]*(phi_b[i][0]/phi[i] - 1/b[i][0]);
                u_ax[s+1][j] = (phi_bx[i][j][0]*(w[i]-u) - u_a[s+1]*q_x[j] - u_x[j]*phi_b[i][0]) / q;

                for (m = 0; m < ai->nd[t]; m++) {
                    u_ax[s+2+m][j] = (-phi_xx[i][j][m]*(w[i] - u) + phi_x[i][m]*u_x[j] - u_a[s+2+m]*q_x[j]) / q;
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

void nrbf_cauchy_valgrad_d2(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, j, m, s, k;
    double p = 0., q = 0., u, p_xx, phi_bxx, phi_xxx;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi,
        *u_a = ac->u_a[t], **u_ax = ac->u_ax[t], *q_x = ac->q_x, **phi_x = ac->phi_x, **phi_b = ac->phi_b,
        *u_x = ac->u_x[t], ***phi_xx = ac->phi_xx, **q_xx = ac->q_xx, *u_xx = ac->u_xx[t],
        **u_xy = ac->u_xy[t], ***phi_bx = ac->phi_bx, **u_axx = ac->u_axx[t], ***u_axy = ac->u_axy[t];

    nrbf_cauchy_valgrad_d1(point, x, t, mixed_derivatives, anns_instance);

    u = ac->u[t];
    q = ac->q;
    for (j = 0; j < ai->nd[t]; j++) {
        p_xx = 0;
        q_xx[j][j] = 0;
        for (i = 0; i < ai->nn[t]; i++) {
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
                for (i = 0; i < ai->nn[t]; i++) {
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
    for (i = 0, s = 0; i < ai->nn[t]; i++, s+= ai->S[t]) {
        if (phi[i] != 0) {
            for (j = 0; j < ai->nd[t]; j++) {
                u_axx[s][j] = (phi_xx[i][j][j] - 2*u_ax[s][j]*q_x[j] - u_a[s]*q_xx[j][j]) / q;

                phi_bxx = phi_xx[i][j][j]*phi_bx[i][j][0] / phi_x[i][j] + 2*phi_x[i][j]*(phi_bx[i][j][0]/phi[i] + xc[i][j]*phi_b[i][0]);
                u_axx[s+1][j] = (phi_bxx*(w[i]-u) - 2*u_ax[s+1][j]*q_x[j] - 2*u_x[j]*phi_bx[i][j][0] - u_a[s+1]*q_xx[j][j] - u_xx[j]*phi_b[i][0])/q;

                for (m = 0; m < ai->nd[t]; m++) {
                    phi_xxx = -2 * (xc[i][j]*(phi_xx[i][m][j]*phi[i] + phi_x[i][m]*phi_x[i][j]) + phi[i]*phi_x[i][m]*a[i][0] + (m==j ? phi[i]*phi_x[i][j]*a[i][0] : 0) );
                    u_axx[s+2+m][j] = (-phi_xxx*(w[i]-u) - 2*u_ax[s+2+m][j]*q_x[j] - 2*u_x[j]*(-phi_xx[i][m][j]) - u_xx[j]*(-phi_x[i][m]) - u_a[s+2+m]*q_xx[j][j]) / q;
                }

                if (mixed_derivatives) {
                    for (k = 0; k < j; k++) {
                        u_axy[s][j][k] = (phi_xx[i][j][k] - u_ax[s][j]*q_x[k] - u_ax[s][k]*q_x[j] - u_a[s]*q_xx[j][k]) / q;

                        phi_bxx = phi_xx[i][j][k]*phi_bx[i][j][0] / phi_x[i][j] + 2*phi_x[i][j]*(phi_bx[i][k][0]/phi[i] + xc[i][k]*phi_b[i][0]);
                        u_axy[s+1][j][k] = (phi_bxx*(w[i]-u) - u_ax[s+1][j]*q_x[k] - u_ax[s+1][k]*q_x[j] - u_xy[j][k]*phi_b[i][0] - u_x[j]*phi_bx[i][k][0] - u_x[k]*phi_bx[i][j][0] - u_a[s+1]*q_xx[j][k])/q;

                        for (m = 0; m < ai->nd[t]; m++) {
                            phi_xxx = -2 * (xc[i][j]*(phi_xx[i][m][k]*phi[i] + phi_x[i][m]*phi_x[i][k]) + (m==j ? phi[i]*phi_x[i][k]*a[i][0] : 0) );
                            u_axy[s+2+m][j][k] = (-phi_xxx*(w[i]-u) - u_ax[s+2+m][j]*q_x[k] - u_ax[s+2+m][k]*q_x[j] - u_xy[j][k]*(-phi_x[i][m]) - u_x[j]*(-phi_xx[i][m][k]) - u_x[k]*(-phi_xx[i][m][j]) - u_a[s+2+m]*q_xx[j][k] ) / q;
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

static anns_func_t val_functions[] = {nrbf_cauchy_val, nrbf_cauchy_val_d1, nrbf_cauchy_val_d2, };
static anns_func_t valgrad_functions[] = {nrbf_cauchy_valgrad, nrbf_cauchy_valgrad_d1, nrbf_cauchy_valgrad_d2, };

void nrbf_cauchy_loadfuncs(int t, void *anns_instance) {
    anns_instance_t *ai = anns_instance;

    ai->ann_val[t] = val_functions;
    ai->ann_valgrad[t] = valgrad_functions;
}
