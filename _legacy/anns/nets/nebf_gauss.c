#include "nebf_gauss.h"

void nebf_gauss_val(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, j, I, Nd = ac->nd[t];
    double p = 0., q = 0., min_lambda = 1e100, x0;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi;

    for (i = 0, I = ai->I[t]; i < ai->nn[t]; i++, I += ai->S[t]) {
        r[i] = 0.;
        w[i] = x[I];

        for (j = 0; j < Nd; j++) {
            b[i][j] = x[I + j + 1];
            a[i][j] = 1.0/(b[i][j]*b[i][j]);

            x0 = point[j] - x[I + j + Nd + 1];
            xc[i][j] = x0*a[i][j];
            r[i] += x0*xc[i][j];
        }
        r[i] *= 0.5;

        // Переменная для нормировки - минимальное значение для экспоненты
        if (r[i] < min_lambda) {
            min_lambda = r[i];
        }
    }

    for (i = 0, I = ai->I[t]; i < ai->nn[t]; i++, I += ai->S[t]) {

        phi[i] = exp(-(r[i]-min_lambda));

        // В случае, если после нормировки функция все равно обращается в нуль, то нет смысла её вычислять
        if (phi[i] != 0) {
            q += phi[i];
            p += w[i] * phi[i];
        }
    }

    // q никогда не обращается в нуль из-за нормировки
    ac->u[t] = p / q;
    ac->q = q;
}

void nebf_gauss_val_d1(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, j;
    double p = 0., q = 0., u, p_x, min_lambda = 1e100;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi,
        *q_x = ac->q_x, **phi_x = ac->phi_x, **phi_b = ac->phi_b, *u_x = ac->u_x[t];

    nebf_gauss_val(point, x, t, mixed_derivatives, anns_instance);

    u = ac->u[t];
    q = ac->q;

    for (j = 0; j < ai->nd[t]; j++) {
        p_x = 0.;
        q_x[j] = 0.;
        for (i = 0; i < ai->nn[t]; i++) {
            if (phi[i] != 0) {
                phi_x[i][j] = -xc[i][j] * phi[i];
                q_x[j] += phi_x[i][j];
                p_x += w[i] * phi_x[i][j];
            }
        }
        u_x[j] = (p_x - u*q_x[j]) / q;
    }
}

void nebf_gauss_val_d2(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, j, k;
    double p = 0., q = 0., u, p_xx, min_lambda = 1e100;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi,
        *q_x = ac->q_x, **phi_x = ac->phi_x, **phi_b = ac->phi_b, *u_x = ac->u_x[t], ***phi_xx = ac->phi_xx,
        **q_xx = ac->q_xx, *u_xx = ac->u_xx[t], **u_xy = ac->u_xy[t];

    nebf_gauss_val_d1(point, x, t, mixed_derivatives, anns_instance);

    u = ac->u[t];
    q = ac->q;

    for (j = 0; j < ai->nd[t]; j++) {
        p_xx = 0;
        q_xx[j][j] = 0;
        for (i = 0; i < ai->nn[t]; i++) {
            if (phi[i] != 0) {
                phi_xx[i][j][j] = -phi_x[i][j]*xc[i][j] - a[i][j]*phi[i];
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
//                        phi_xx[i][j][k] = xc[i][j] * xc[i][k] * phi[i];
                        phi_xx[i][j][k] = -phi_x[i][j] * xc[i][k];
                        q_xx[j][k] += phi_xx[i][j][k];
                        p_xx += w[i] * phi_xx[i][j][k];
                    }
                }
                u_xy[j][k] = (p_xx - q_x[j]*u_x[k] - q_x[k]*u_x[j] - u*q_xx[j][k]) / q;
            }
        }
    }
}

void nebf_gauss_valgrad(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, m, s, Nd = ac->nd[t];
    double p = 0., q = 0., u, min_lambda = 1e100;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi,
        *u_a = ac->u_a[t], *q_x = ac->q_x, **phi_x = ac->phi_x, **phi_b = ac->phi_b;

    nebf_gauss_val(point, x, t, mixed_derivatives, anns_instance);

    u = ac->u[t];
    q = ac->q;
    // Вычисление градиента
    for (i = 0, s = 0; i < ai->nn[t]; i++, s+= ai->S[t]) {
        if (phi[i] != 0) {
            u_a[s] = phi[i] / q;

            for (m = 0; m < Nd; m++) {
                phi_x[i][m] = -xc[i][m] * phi[i];
                phi_b[i][m] = -b[i][m]*xc[i][m]*phi_x[i][m];

                u_a[s+1+m] = phi_b[i][m] * (w[i] - u) / q;
                u_a[s+1+m+Nd] = -phi_x[i][m] * (w[i] - u) / q;
            }
        } else {
            u_a[s] = 0;
            for (m = 0; m < Nd; m++) {
                u_a[s+1+m] = 0;
                u_a[s+1+Nd+m] = 0;
            }
        }
    }
}

void nebf_gauss_valgrad_d1(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, j, m, s, k, Nd = ac->nd[t];
    double p = 0., q = 0., u, p_x;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi,
        *u_a = ac->u_a[t], **u_ax = ac->u_ax[t], *q_x = ac->q_x, **phi_x = ac->phi_x, **phi_b = ac->phi_b,
        *u_x = ac->u_x[t], ***phi_xx = ac->phi_xx, *u_xx = ac->u_xx[t], ***phi_bx = ac->phi_bx;

    nebf_gauss_valgrad(point, x, t, mixed_derivatives, anns_instance);

    for (i = 0; i < ai->nn[t]; i++) {
        if (phi[i] != 0) {
            for (j = 0; j < Nd; j++) {
                phi_xx[i][j][j] = -phi_x[i][j]*xc[i][j]-a[i][j]*phi[i];
                for (k = 0; k < j; k++) {
                    phi_xx[i][j][k] = -xc[i][k] * phi_x[i][j];
                    phi_xx[i][k][j] = phi_xx[i][j][k];
                }
            }
        }
    }

    u = ac->u[t];
    q = ac->q;

    for (j = 0; j < Nd; j++) {
        p_x = 0.;
        q_x[j] = 0.;
        for (i = 0; i < ai->nn[t]; i++) {
            if (phi[i] != 0) {
                phi_x[i][j] = -xc[i][j] * phi[i];
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

                for (m = 0; m < ai->nd[t]; m++) {
                    phi_bx[i][j][m] = -b[i][m]*xc[i][m]*phi_xx[i][j][m] - (j == m ? phi_x[i][j]/b[i][j] : 0);
                    u_ax[s+1+m][j] = (phi_bx[i][j][m]*(w[i]-u) - u_a[s+1+m]*q_x[j] - u_x[j]*phi_b[i][m]) / q;

                    u_ax[s+1+Nd+m][j] = (-phi_xx[i][j][m]*(w[i] - u) + phi_x[i][m]*u_x[j] - u_a[s+1+Nd+m]*q_x[j]) / q;
                }
            }
        } else {
            for (j = 0; j < ai->nd[t]; j++) {
                u_ax[s][j] = 0;
                for (m = 0; m < ai->nd[t]; m++) {
                    u_ax[s+1+m][j] = 0;
                    u_ax[s+1+Nd+m][j] = 0;
                }
            }
        }
    }
}

void nebf_gauss_valgrad_d2(double *point, double *x, int t, int mixed_derivatives, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int i, j, m, s, k, Nd = ac->nd[t];
    double p = 0., q = 0., u, p_xx, phi_bxx, phi_xxx;
    double *r = ac->r, **a = ac->a, *w = ac->w, **b = ac->b, **xc = ac->xc, *phi = ac->phi,
        *u_a = ac->u_a[t], **u_ax = ac->u_ax[t], *q_x = ac->q_x, **phi_x = ac->phi_x, **phi_b = ac->phi_b,
        *u_x = ac->u_x[t], ***phi_xx = ac->phi_xx, **q_xx = ac->q_xx, *u_xx = ac->u_xx[t],
        **u_xy = ac->u_xy[t], ***phi_bx = ac->phi_bx, **u_axx = ac->u_axx[t], ***u_axy = ac->u_axy[t];

    nebf_gauss_valgrad_d1(point, x, t, mixed_derivatives, anns_instance);

    u = ac->u[t];
    q = ac->q;
    for (j = 0; j < Nd; j++) {
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
            for (j = 0; j < Nd; j++) {
                u_axx[s][j] = (phi_xx[i][j][j] - 2*u_ax[s][j]*q_x[j] - u_a[s]*q_xx[j][j]) / q;

                for (m = 0; m < Nd; m++) {
                    phi_bxx = b[i][m]*xc[i][m]*xc[i][m]*phi_xx[i][m][j] - (m == j ? 2*(2*phi_xx[i][j][j]+phi[i]*a[i][m])/b[i][m] : 0 );
                    u_axx[s+1+m][j] = (phi_bxx*(w[i]-u) - 2*u_ax[s+1+m][j]*q_x[j] - 2*u_x[j]*phi_bx[i][j][m] - u_a[s+1+m]*q_xx[j][j] - u_xx[j]*phi_b[i][m])/q;

                    phi_xxx = xc[i][m] * phi_xx[i][j][j] + (m == j ? 2*phi_x[i][j]*a[i][j] : 0);
                    u_axx[s+1+Nd+m][j] = (phi_xxx*(w[i]-u) - 2*u_ax[s+1+Nd+m][j]*q_x[j] - 2*u_x[j]*(-phi_xx[i][m][j]) - u_xx[j]*(-phi_x[i][m]) - u_a[s+Nd+1+m]*q_xx[j][j]) / q;
                }

                if (mixed_derivatives) {
                    for (k = 0; k < j; k++) {
                        u_axy[s][j][k] = (phi_xx[i][j][k] - u_ax[s][j]*q_x[k] - u_ax[s][k]*q_x[j] - u_a[s]*q_xx[j][k]) / q;

                        for (m = 0; m < Nd; m++) {
                            phi_bxx = b[i][m]*xc[i][m]*xc[i][m]*phi_xx[i][j][k] - ( (m == j ? 1 : 0) + (m == k ? 1 : 0) )*(2*phi_xx[i][j][k]/b[i][m]);
                            u_axy[s+1+m][j][k] = (phi_bxx*(w[i]-u) - u_ax[s+1+m][j]*q_x[k] - u_ax[s+1+m][k]*q_x[j] - u_xy[j][k]*phi_b[i][0] - u_x[j]*phi_bx[i][k][0] - u_x[k]*phi_bx[i][j][0] - u_a[s+1+m]*q_xx[j][k])/q;

//                            phi_xxx = (-phi_xx[i][m][j]*phi_x[i][k] - phi_xx[i][m][k]*phi_x[i][j] + phi_x[i][m]*phi_xx[i][j][k] ) / phi[i];

//                            if (phi_xxx != phi_xxx) {phi_xxx = 0;}
                            phi_xxx = xc[i][m]*phi_xx[i][j][k] + (m == j ? phi_x[i][k]*a[i][m] : 0) + (m == k ? phi_x[i][j]*a[i][m] : 0);

                            u_axy[s+Nd+1+m][j][k] = (phi_xxx*(w[i]-u) - u_ax[s+Nd+1+m][j]*q_x[k] - u_ax[s+Nd+1+m][k]*q_x[j] - u_xy[j][k]*(-phi_x[i][m]) - u_x[j]*(-phi_xx[i][m][k]) - u_x[k]*(-phi_xx[i][m][j]) - u_a[s+Nd+1+m]*q_xx[j][k] ) / q;
                        }
                    }
                }
            }
        } else {
            for (j = 0; j < Nd; j++) {
                u_axx[s][j] = 0;
                for (m = 0; m < Nd; m++) {
                    u_axx[s+1+m][j] = 0;
                    u_axx[s+1+Nd+m][j] = 0;
                }
                if (mixed_derivatives) {
                    for (k = 0; k < j; k++) {
                        u_axy[s][j][k] = 0;
                        for (m = 0; m < Nd; m++) {
                            u_axy[s+1+m][j][k] = 0;
                            u_axy[s+1+Nd+m][j][k] = 0;
                        }
                    }
                }
            }
        }
    }
    return;
}

static anns_func_t val_functions[] = {nebf_gauss_val, nebf_gauss_val_d1, nebf_gauss_val_d2, };
static anns_func_t valgrad_functions[] = {nebf_gauss_valgrad, nebf_gauss_valgrad_d1, nebf_gauss_valgrad_d2, };

void nebf_gauss_loadfuncs(int t, void *anns_instance) {
    anns_instance_t *ai = anns_instance;

    ai->ann_val[t] = val_functions;
    ai->ann_valgrad[t] = valgrad_functions;
}
