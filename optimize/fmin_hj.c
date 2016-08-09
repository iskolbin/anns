#include "fmin_hj.h"

fmin_hj_args_t *new_fmin_hj_args(int ndelta, double default_delta) {
    fmin_hj_args_t *args = NULL;
    int i;

    args = malloc(sizeof(args));

    assert(args);

    args->max_iterations = 0xFFFF;
    args->max_epsilon = 1E-5;
    args->lambda = 1.0;
    args->alpha = 2.0;
    args->default_delta = default_delta;

    if (ndelta > 0) {
        args->ndelta = ndelta;
        args->delta = calloc(ndelta, sizeof *args->delta);
        for (i = 0; i < ndelta; i++) {
            args->delta[i] = default_delta;
        }
    } else {
        args->ndelta = 0;
        args->delta = NULL;
    }

    return args;
}

void delete_fmin_hj_args(fmin_hj_args_t *args) {
    if (args) {
        if (args->delta) {
            free(args->delta);
        }
        free(args);
    }
}

double *fmin_hj(double (*func)(double *, int, void *), double *result, int n, double *x0, void *instance, void *method_instance, void (*callback)(double *, int, void *, int, double, void*), void *callback_instance) {
    double lambda, alpha, f_yi, f_yi1, f_xk, max_epsilon;
    double *yi = NULL, *xk = NULL, xk1, *delta = NULL;
    int i, stop_key = 1;
    long int iteration = 0, max_iterations;
    fmin_hj_args_t *args = NULL;

    assert(result && n && x0);

    yi = result;
    xk = calloc(n, sizeof *xk);
    delta = calloc(n, sizeof *delta);

    if (method_instance) {
        args = (fmin_hj_args_t *) method_instance;
        max_iterations = args->max_iterations;
        max_epsilon = args->max_epsilon;
        lambda = args->lambda;
        alpha = args->alpha;
        if (args->delta && args->ndelta>0) {
            for (i = 0; i < n; i++) {
                delta[i] = args->delta[i%args->ndelta];
            }
        } else if (args) {
            for (i = 0; i < n; i++) {
                delta[i] = args->default_delta;
            }
        } else {
            for (i = 0; i < n; i++) {
                delta[i] = 0.25;
            }
        }
    } else {
        max_iterations = 0xFFFF;
        max_epsilon = 1E-8;
        lambda = 1.0;
        alpha = 2.0;
        for (i = 0; i < n; i++) {
            delta[i] = 1.0;
        }
    }

    for (i = 0; i < n; i++) {
        yi[i] = x0[i];
        xk[i] = x0[i];
    }

    f_xk = func(x0, n, instance);
    f_yi = f_xk;

    do {
        for (i = 0; i < n; i++) {
            yi[i] += delta[i]; /* evaluate y + delta*d */
            f_yi1 = func(yi, n, instance);
            if (f_yi <= f_yi1) {
                yi[i] -= 2*delta[i]; /* evaluate y - delta*d */
                f_yi1 = func(yi, n, instance);
                if (f_yi <= f_yi1) {
                    yi[i] += delta[i]; /* return y to default */
                } else {
                    f_yi = f_yi1;
                }
            } else {
                f_yi = f_yi1;
            }
        }

        if (f_yi < f_xk) {
            for (i = 0; i < n; i++) {
                xk1 = yi[i];
                yi[i] = (1 + lambda) * yi[i] - xk[i];
                xk[i] = xk1;
            }

            f_xk = f_yi;
            f_yi = func(yi, n, instance);
        } else {
            stop_key = 0;
            for (i = 0; i < n; i++) {
                if (delta[i] > max_epsilon) {
                    stop_key = 1;
                    delta[i] /= alpha;
                }
            }
            for (i = 0; i < n; i++) {
                yi[i] = xk[i];
            }

            f_yi = f_xk;
        }

        iteration++;

        if (callback) {
            callback(yi, n, instance, iteration, f_yi, callback_instance);
        }
    } while ((iteration < max_iterations) && stop_key);

    free(xk);
    free(delta);

    return result;
}
