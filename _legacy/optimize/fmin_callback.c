#include "fmin_callback.h"

void fmin_callback_fprintf(fmin_callback_args_t *args) {
    if (args && args->fout) {
        fprintf(args->fout, "i=%d f=%g gnorm=%g gnorm2=%g\n", args->iteration, args->f, args->gnorm, args->gnorm2);
    }
}

void fmin_callback_fprintf_mse(fmin_callback_args_t *args) {
    if (args && args->fout && args->mse) {
        fprintf(args->fout, "i=%d f=%g gnorm=%g gnorm2=%g mse=%g\n", args->iteration, args->f, args->gnorm, args->gnorm2, args->mse(args->x, args->n, args->instance));
    }
}
