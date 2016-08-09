#include "fmin_linesearch.h"

fmin_linesearch_args_t fmin_linesearch_default_args = {
    1e-4, 0.1, 2.0,
    0.2, 0.1,
    1e-3, 10
};

double fmin_linesearch_goldenratio(fmin_vars_t *vars, fmin_linesearch_args_t *args) {
    double a = 0, b = 5, fa, fb, x1, x2, fx1, fx2, alpha;
    int i, step, n;
    double *x, *xtemp, *d, *gtemp;
    void *instance;
    double phi = (1+sqrt(5)) / 2;

    args = args ? args : &fmin_linesearch_default_args;

    assert(vars);
    assert(args);

    x = vars->x; xtemp = vars->xtemp;
    n = vars->n; gtemp = vars->gtemp;
    d = vars->d; instance = vars->instance;

    fa = vars->f0;

    for (i = 0; i < n; i++) {
        xtemp[i] = x[i] + b * d[i];
    }
    fb = vars->val(xtemp, n, instance);
    vars->nval++;

    for (step = 1; step < args->max_steps; step++) {
        x1 = b - (b-a) / phi;
        for (i = 0; i < n; i++) {
            xtemp[i] = x[i] + x1 * d[i];
        }
        fx1 = vars->val(xtemp, n, instance);
        vars->nval++;

        x2 = a + (b-a) / phi;
        for (i = 0; i < n; i++) {
            xtemp[i] = x[i] + x2 * d[i];
        }
        fx2 = vars->val(xtemp, n, instance);
        vars->nval++;

        if (fx1 > fx2) {
            a = x1;
        } else {
            b = x2;
        }

        if (fabs(b-a) < args->tol) {
            break;
        }
    }

    alpha = 0.5*(b+a);
    for (i = 0; i < n; i++) {
        xtemp[i] = x[i] + alpha * d[i];
    }
    vars->fa = vars->valgrad(gtemp, xtemp, n, instance);
    vars->ngrad++;
    vars->nval++;

    return alpha;
}

double fmin_linesearch_cubic(double a, double fa, double fpa, double b, double fb, double c, double fc) {
    // f(x) - A*(x-a)^3 + B*(x-a)^2 + C*(x-a) + D
    double A, B, C = fpa, D = fa;
    double db = b-a, dc = c-a, denom, radical;
    double tb = fb - D - C*db, tc = fc - D - C*dc;
    double db2 = db*db, dc2 = dc*dc;

    if (!db || !dc || b==c) return FMIN_LINESEARCH_FAIL;

    denom = (db2*dc2)*(b-c);

    A = dc2*tb - db2*tc;
    B = -dc*dc2*tb + db*db2*tc;

    A /= denom;
    B /= denom;

    radical = B*B - 3*A*C;

    if (radical < 0 || !A) return FMIN_LINESEARCH_FAIL;

    return a + (-B + sqrt(radical)) / (3*A);
}

double fmin_linesearch_quad(double a, double fa, double fpa, double b, double fb) {
    // f(x) = B*(x-a)^2 + C*(x-a) + D
    double B, D = fa, C= fpa;
    double db = b-a;

    if (!db) return FMIN_LINESEARCH_FAIL;

    B = (fb - D - C*db) / (db*db);
    if (B <= 0) return FMIN_LINESEARCH_FAIL;

    return a - C / (2*B);
}

double fmin_linesearch_zoom(double a_lo, double a_hi, double phi_lo, double phi_hi, double derphi_lo, fmin_vars_t *vars, fmin_linesearch_args_t *args) {
    int maxiter;
    double delta1, delta2;
    double phi_rec = vars->f0, a_rec = 0, phi0 = vars->f0, derphi0 = vars->g0;
    double *x = vars->x, *xtemp = vars->xtemp, *d = vars->d, *gtemp = vars->gtemp;
    void *instance = vars->instance;
    int n = vars->n, i = 0, j;
    double c1, c2, dalpha, a, b, cchk, qchk, a_j, phi_aj, derphi_aj;

    args = args ? args : &fmin_linesearch_default_args;

    c1 = args->c1; c2 = args->c2;
    delta1 = args->delta1; delta2 = args->delta2;
    maxiter = args->max_steps;

    for (;;) {
        dalpha = a_hi - a_lo;
        if (dalpha < 0) {
            a = a_hi; b = a_lo;
        } else {
            a = a_lo; b = a_hi;
        }

        if (i > 0) {
            cchk = delta1*dalpha;
            a_j = fmin_linesearch_cubic(a_lo, phi_lo, derphi_lo, a_hi, phi_hi, a_rec, phi_rec);
        }

        if (!i || (a_j == FMIN_LINESEARCH_FAIL) || (a_j > (b-cchk)) || (a_j < (a+cchk)) ) {
            qchk = delta2*dalpha;
            a_j = fmin_linesearch_quad(a_lo, phi_lo, derphi_lo, a_hi, phi_hi);

            if (a_j == FMIN_LINESEARCH_FAIL || (a_j > b - qchk) || (a_j < a + qchk)) {
                a_j = a_lo + 0.5*dalpha;
            }
        }

        for (j = 0; j < n; j++) xtemp[j] = x[j] + a_j*d[j];
        phi_aj = vars->val(xtemp, n, instance);
        vars->nval++;

        if ((phi_aj > phi0 + c1*a_j*derphi0) || (phi_aj >= phi_lo)) {
            phi_rec = phi_hi;
            a_rec = a_hi;
            a_hi = a_j;
            phi_hi = phi_aj;
        } else {
            vars->grad(gtemp, xtemp, n, instance);
            vars->ngrad++;

            derphi_aj = 0;
            for (j = 0; j < n; j++) {
                derphi_aj += d[j]*gtemp[j];
            }

            if (fabs(derphi_aj) <= -c2*derphi0) {
                vars->a = a_j;
                vars->fa = phi_aj;
                vars->ga = derphi_aj;
                break;
            }

            if (derphi_aj*(a_hi - a_lo) >= 0) {
                phi_rec = phi_hi;
                a_rec = a_hi;
                a_hi = a_lo;
                phi_hi = phi_lo;
            } else {
                phi_rec = phi_lo;
                a_rec = a_lo;
            }

            a_lo = a_j;
            phi_lo = phi_aj;
            derphi_lo = derphi_aj;
        }

        i++;

        if (i > maxiter) {
            vars->a = a_j;
            vars->fa = phi_aj;
            vars->ga = FMIN_LINESEARCH_FAIL;
            break;
        }
    }

    return vars->a;
}

double fmin_linesearch_wolfe(fmin_vars_t *vars, fmin_linesearch_args_t *args) {
    int maxiter, step = 0;
    double delta1, delta2;
    double phi0 = vars->f0, derphi0 = vars->g0;
    double *x = vars->x, *xtemp = vars->xtemp, *d = vars->d, *gtemp = vars->gtemp;
    void *instance = vars->instance;
    int n = vars->n, i = 0;
    double c1, c2, tmp;
    double alpha0, alpha1, alpha2, alpha_star, phi_star, derphi_star, phi_a1, phi_a0, derphi_a0, derphi_a1;

    args = args ? args : &fmin_linesearch_default_args;

    c1 = args->c1; c2 = args->c2;
    delta1 = args->delta1; delta2 = args->delta2;
    maxiter = args->max_steps;

    phi0 = vars->f0;
    derphi0 = vars->g0;

    alpha0 = 0;
    if (vars->old_fa != FMIN_LINESEARCH_FAIL) {
        tmp = 1.01*2*(phi0 - vars->old_fa) / derphi0;
        alpha1 = 1 < tmp ? 1 : tmp;
    } else {
        alpha1 = 1;
    }

    if (alpha1 < 0) {
        alpha1 = 1;
    }

    if (alpha1 == 0) {
        alpha_star = FMIN_LINESEARCH_FAIL;
        phi_star = phi0;
        phi0 = vars->old_fa;
        derphi_star = FMIN_LINESEARCH_FAIL;
    }

    for (i = 0; i < n; i++) {
        xtemp[i] = x[i] + alpha1*d[i];
    }
    phi_a1 = vars->val(xtemp, n, instance);
    vars->nval++;

    phi_a0 = phi0;
    derphi_a0 = derphi0;

    for (;;) {
        if (alpha1 == 0.) {
            break;
        } else {
            if ( (phi_a1 > phi0 + c1*alpha1*derphi0) || (phi_a1 >= phi_a0 && (i > 1)) ) {
                fmin_linesearch_zoom(alpha0, alpha1, phi_a0, phi_a1, derphi_a0, vars, args);
                break;
            }

//            if (vars->ga == FMIN_LINESEARCH_FAIL) {
            for (i = 0; i < n; i++) {
                xtemp[i] = x[i] + alpha1*d[i];
            }
            vars->grad(gtemp, xtemp, n, instance);
            vars->ngrad++;
//            }

            derphi_a1 = 0;
            for (i = 0; i < n; i++) {
                derphi_a1 += gtemp[i]*d[i];
            }

            if (fabs(derphi_a1) <= -c2*derphi0) {
                vars->a = alpha1;
                vars->fa = phi_a1;
                vars->ga = derphi_a1;
                break;
            }

            if (derphi_a1 >= 0) {
                fmin_linesearch_zoom(alpha1, alpha0, phi_a1, phi_a0, derphi_a1, vars, args);
                break;
            }

            alpha2 = 2*alpha1;

            step++;

            alpha0 = alpha1;
            alpha1 = alpha2;
            phi_a0 = phi_a1;

            for (i = 0; i < n; i++) {
                xtemp[i] = x[i] + alpha1*d[i];
            }
            phi_a1 = vars->val(xtemp, n, instance);
            vars->nval++;

            derphi_a0 = derphi_a1;

            if ( step > maxiter) {
                vars->a = alpha1;
                vars->fa = phi_a1;
                vars->ga = FMIN_LINESEARCH_FAIL;
                break;
            }
        }
    }
    return vars->a;
}
