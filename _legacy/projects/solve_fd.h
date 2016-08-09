#ifndef ANNS_SCHEME
#define ANNS_SCHEME

#include <time.h>
#include "../anns/anns.h"

#define SOLVE_OK 0
#define SOLVE_WRITE_ERROR 1

#define N_TRIALS 1

#define TOLERANCE 1e-5

#define DBL_FORMAT "%.16g"

//#define NFILEOUT

#ifndef NTRACE
#define tfprintf(...) fprintf(__VA_ARGS__)
#else
#define tfprintf(...)
#endif


#include "dir_par_6_fd.h"


int solve(int argc, char **argv) {
    int trial;
    anns_instance_t *ai = NULL;
    double *x = NULL, *x_best = NULL, *x_ = NULL, J_best = 1e100, J, J_after;
    double *Work = NULL;
    double s, mse, tmp,  *neuron = NULL;
    int i, j, k, p, t;
    int result = SOLVE_OK;
    long int seed = time(NULL);

    FILE *clog = stdout;
    FILE *fout = NULL;

    tfprintf(clog, "[Begin application]\n");
    tfprintf(clog, "[Finite difference time split]\n");
#ifdef NFILEOUT
    tfprintf(clog, "[File output is off]\n");
#endif
#ifndef NDEBUG
    tfprintf(clog, "[Debug is on]\n");
#else
    tfprintf(clog, "[Debug is off]\n");
#endif
    tfprintf(clog, "[Init random generator with seed %ld]\n", seed);
    srand(seed);

    tfprintf(clog, "[Allocate memory]\n");
    ai = ai_new(data, NULL);
    x = ai_createann(ai);
    Work = calloc(4*ai->A, sizeof(*Work));

    tfprintf(clog, "[Initialize ai]\n");
    ai_init(ai);
    ai_initann(ai, x);

    for (i = 1; i < FD_STEPS_EVAL; i++) {
        tfprintf(clog, "[Begin step %d of %d]\n", i, FD_STEPS_EVAL-1);
        ai_update(ai, x, i);
        cg_descent(x, ai->A, NULL, NULL, TOLERANCE, anns_val_f, anns_grad_f, anns_valgrad_f, Work, ai);
        tfprintf(clog, "[End step %d of %d]\n", i, FD_STEPS_EVAL-1);
    }
    ac_put(ai->cache, x, FD_STEPS_EVAL-1);

#ifndef NFILEOUT
        tfprintf(clog, "[Writing ouput to file \"%s\"]\n", ai->fileAme);

        fout = fopen(ai->fileAme, "w+");
        if (fout) {
            fprintf(fout, "tag=%d type=%d index=%d sub=%d nsteps=%d\n", ai->problemtag, ai->typetag, ai->index, ai->subindex, FD_STEPS);
            for (k = 0; k < FD_STEPS_EVAL; k++) {
                fprintf(fout, "k=%d\n", k);
                for (t = 0; t < ai->nnets; t++) {
                    fprintf(fout, "T=%d nn=%d nd=%d\n", ai->T[t], ai->nn[t], ai->nd[t]);
                    x_ = ai_get(ai, k);
                    for (i = 0; i < ai->nn[t]; i++) {
                        neuron = x_ + ai->I[t] + i*ai->S[t];
                        fprintf(fout, DBL_FORMAT " " DBL_FORMAT " " DBL_FORMAT, neuron[0], neuron[1], neuron[2]);
                        for (j = 1; j < ai->nd[t]; j++) {
                            fprintf(fout, " " DBL_FORMAT, neuron[2+j]);
                        }
                        fprintf(fout, "\n");
                    }
                }
            }
            tfprintf(clog, "[Done]\n");
        } else {
            tfprintf(clog, "[Error durning opening file]\n");
            result = SOLVE_WRITE_ERROR;
        }

        fclose(fout);
#endif

    tfprintf(clog, "[Release memory]\n");
    free(x);
    free(x_best);
    free(Work);
    ai_delete(ai);
    tfprintf(clog, "[End application]\n");

    return result;
}
#endif // ANNS_SCHEME
