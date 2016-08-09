#include "fmin_cg.h"

fmin_cg_args_t fmin_cg_default_args = {
    -1, NULL,
    FMIN_CG_N, FMIN_CG_HEURISTICS_P,
    fmin_linesearch_wolfe,
};

fmin_linesearch_args_t fmin_linesearch_default_args_ = {
    1e-4, 0.1, 2.0,
    1e-3, 15
};

#define FMIN_CG_INF_NORM

//
// void linesearch_armijo(fmin_linesearch_vars_t* vars, fmin_linesearch_args_t *args) {
//    double alpha0 = 1, alpha1, alpha2;
//    double derphi0 = 0, derphi1;
//    int i, n = vars->n;
//    double *xtemp = vars->xtemp, *x = vars->x, *g = vars->g, *gtemp = vars->gtemp, *d = vars->d;
//    double phi0, phi_a0, phi_a1, phi_a2, factor, a, b;
//    void *instance = vars->instance;
//
//    args = args ? args : &fmin_linesearch_default_args_;
//
//    phi0 = vars->val(x, n, instance);
//    vars->nval++;
//
//    for (i = 0; i < vars->n; i++) {
//        xtemp[i] = x[i] + alpha0 * d[i];
//        derphi0 += g[i] * d[i];
//    }
//
//    phi_a0 = vars->val(xtemp, n, instance);
//    vars->nval++;
//
//    if (phi_a0 <= phi0 + args->c1*alpha0*derphi0) {
//        return;
//    }
//
//    alpha1 = -0.5 * derphi0 * alpha0 * alpha0 / (phi_a0 - phi0 - derphi0 * alpha0);
//    for (i = 0; i < vars->n; i++) {
//        xtemp[i] = x[i] + alpha1 * d[i];
//    }
//    phi_a1 = vars->val(xtemp, n, instance);
//    vars->nval++;
//
//    if (phi_a1 <= phi0 + args->c1*alpha1*derphi0) {
//        return;
//    }
//
//    // amin
//    while (alpha1 > 0) {
//        factor = alpha0*alpha0 * (alpha1-alpha0);
//
//        a = alpha0*alpha0 * (phi_a1 - phi0 - derphi0*alpha1) - alpha1*alpha1 * (phi_a0 - phi0 - derphi0*alpha0);
//        a /= factor;
//
//        b = -alpha0*alpha0*alpha0 * (phi_a1 - phi0 - derphi0*alpha1) + alpha1*alpha1*alpha1 * (phi_a0 - phi0 - derphi0*alpha0);
//        b /= factor;
//
//        alpha2 = (-b + sqrt(fabs(b*b - 3*a*derphi0))) / (3*a);
//        for (i = 0; i < vars->n; i++) {
//            xtemp[i] = x[i] + alpha2 * d[i];
//        }
//        phi_a2 = vars->val(x, n, instance);
//        vars->nval++;
//
//        if (phi_a2 <= phi0 + args->c1*alpha2*derphi0) {
//            return;
//        }
//
//        if ( ((alpha1 - alpha2) > 0.5*alpha1) || (1-alpha2 / alpha1) < 0.96 ) {
//            alpha2 = 0.5*alpha1;
//        }
//
//        alpha0 = alpha1;
//        alpha1 = alpha2;
//        phi_a0 = phi_a1;
//        phi_a1 = phi_a2;
//    }
//}

double *fmin_cg(double *x,
                int n,

                double (*val) (double *, int, void *),
                void (*grad) (double *, double *, int, void *),
                double (*valgrad) (double *, double *, int, void *),

                double tol,
                void *instance,

                fmin_cg_args_t *args,
                fmin_cg_stats_t *stats,
                double *work_array) {

    void *preallocated = work_array;
    double *d, *xtemp, *g, *gtemp;
    double beta;
    double yk, tmp, dnorm, gnorm;
    int max_iterations, i;
    fmin_vars_t vars = {val, grad, valgrad, x, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, n, instance, 0, 0, 0};

    work_array = preallocated ? work_array : calloc(4*n, sizeof *work_array);
    args = args ? args : &fmin_cg_default_args;

    max_iterations = args->max_iterations == -1 ? 200*n : args->max_iterations;

    d = work_array;
    xtemp = work_array + n;
    g = work_array + 2*n;
    gtemp = work_array + 3*n;

    grad(g, x, n, instance);
    vars.ngrad++;
    gnorm = 0;
    for (i = 0; i < n; i++) {
        d[i] = -g[i];
        gnorm += g[i]*g[i];
    }

    vars.x = x;
    vars.xtemp = xtemp;
    vars.d = d;
    vars.g = g;
    vars.gtemp = gtemp;

    for (vars.niterations = 0; gnorm > tol && vars.niterations < max_iterations; vars.niterations++) {
        vars.f0 = val(x, n, instance);
        vars.nval++;
        vars.old_fa = vars.f0 + 5000;
        vars.g0 = 0;
        for (i = 0; i < n; i++) {
            vars.g0 += g[i]*d[i];
        }

        args->linesearch(&vars, NULL);
        for (i = 0; i < vars.n; i++) {
            x[i] = xtemp[i];
        }

        if (vars.ga == FMIN_LINESEARCH_FAIL) {
            for (i = 0; i < vars.n; i++) {
                gtemp[i] = g[i];
            }
        }

        beta = 0;
        tmp = 0;
        switch (args->cg_type) {
            case FMIN_CG_HS: {
                for (i = 0; i < n; i++) {
                    yk = gtemp[i] - g[i];
                    beta += gtemp[i] * yk;
                    tmp += d[i]*yk;
                }
                beta /= tmp;
                break;
            }

            case FMIN_CG_FR: {
                for (i = 0; i < n; i++) {
                    beta += gtemp[i]*gtemp[i];
                    tmp += g[i]*g[i];
                }
                beta /= tmp;
                break;
            }

            case FMIN_CG_PRP: {
                for (i = 0; i < n; i++) {
                    beta += gtemp[i]*(gtemp[i] - g[i]);
                    tmp += g[i]*g[i];
                }
                beta /= tmp;
                break;
            }

            case FMIN_CG_CD: {
                for (i = 0; i < n; i++) {
                    beta += gtemp[i]*gtemp[i];
                    tmp += -d[i]*g[i];
                }
                beta /= tmp;
                break;
            }

            case FMIN_CG_LS: {
                for (i = 0; i < n; i++) {
                    beta += gtemp[i]*(gtemp[i]-g[i]);
                    tmp += -d[i]*g[i];
                }
                beta /= tmp;
                break;
            }

            case FMIN_CG_DY: {
                for (i = 0; i < n; i++) {
                    beta += gtemp[i]*gtemp[i];
                    tmp += d[i]*(gtemp[i]-g[i]);
                }
                beta /= tmp;
                break;
            }

            case FMIN_CG_N: {
                double dkyk = 0;
                tmp = 0;
                for (i = 0; i < n; i++) {
                    yk = (gtemp[i]-g[i]);
                    dkyk += d[i]*yk;
                    tmp += yk*yk;
                }
                tmp /= dkyk;

                for (i = 0; i < n; i++) {
                    beta += ((gtemp[i]-g[i]) - 2*d[i]*tmp) * gtemp[i];
                }
                beta /= dkyk;

                break;
            }
        }

        if (args->heuristic_type == FMIN_CG_HEURISTICS_P) {
            beta = beta >= 0 ? beta : 0;
        } else if (args->heuristic_type == FMIN_CG_HEURISTICS_N) {
            dnorm = 0;
            gnorm = 0;
            for (i = 0; i < n; i++) {
                dnorm += d[i]*d[i];
                gnorm += g[i]*g[i];
            }
            tmp = -1/(dnorm * (0.01 < gnorm ? 0.01 : gnorm));
            beta = beta > tmp ? beta : tmp;
        }

        gnorm = 0;
        for (i = 0; i < n; i++) {
            d[i] = -gtemp[i] + beta*d[i];
            g[i] = gtemp[i];
#ifdef FMIN_CG_INF_NORM
            tmp = fabs(g[i]);
            gnorm = tmp > gnorm ? tmp : gnorm;
#else
            gnorm += g[i]*g[i];
#endif
        }
    }

    if (preallocated) {
        free(work_array);
    }

    printf("i=%d f=%d g=%d", vars.niterations, vars.nval, vars.ngrad);

    return x;
}

