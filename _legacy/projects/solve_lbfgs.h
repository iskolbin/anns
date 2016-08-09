#ifndef SOLVE_LBFGS_H_INCLUDED
#define SOLVE_LBFGS_H_INCLUDED


#include <time.h>

#include "../anns/anns.h"
#include "../optimize/fmin_lbfgs.h"


//#include "dir_ell_1.h"
#include "dir_ell_1_e.h"
//#include "dir_ell_1_2nets.h"
//#include "dir_ell_2.h"
//#include "dir_par_4.h"
//#include "ibnd_par_1.h"
//#include "ibnd_par_1.h"

#define N_TRIALS 1
#define NFILEOUT
#define RANDSEED (time(NULL))
//#define RANDSEED 123L

#define TOLERANCE 1e-4

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    ) {
        return anns_valgrad_f(g, x, n, instance);
//    int i;
//    double fx = 0.0;
//
//    for (i = 0;i < n;i += 2) {
//        lbfgsfloatval_t t1 = 1.0 - x[i];
//        lbfgsfloatval_t t2 = 10.0 * (x[i+1] - x[i] * x[i]);
//        g[i+1] = 20.0 * t2;
//        g[i] = -2.0 * (x[i] * g[i+1] + t1);
//        fx += t1 * t1 + t2 * t2;
//    }
//    return fx;
}

int solve(int argc, char **argv) {
    int trial;
    anns_instance_t *bi = NULL;
    double *x = NULL, *x_best = NULL, J_best = 1e100, J, J_after;
    double *Work = NULL;
    double s, mse, tmp,  *neuron = NULL;
    int i, j, p, t;
    long int seed = RANDSEED;
    FILE *fout = NULL;
    clock_t start,end;
    double dif;
    lbfgs_parameter_t param;

    param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE  ;
    param.max_linesearch = 50;
//    lbfgs_parameter_init(&param);

    printf("[Begin application]\n");
#ifdef NFILEOUT
    printf("[File output is off]\n");
#endif
    printf("[Init random generator with seed %ld]\n", seed);
    srand(seed);

    printf("[Allocate memory]\n");
    bi = ai_new(data, NULL);
//    x = lbfgs_malloc(bi->A);
//    x_best = lbfgs_malloc(bi->A);
    x = ai_createann(bi);
    x_best = ai_createann(bi);
//    x = calloc(1000, sizeof *x);
//    x_best = calloc(1000, sizeof *x);
//    Work = calloc(4*bi->A, sizeof(*Work));

    printf("[Initialize bi]\n");
    ai_init(bi);
    printf("[Reload bi]\n");
    ai_reload(bi);

    start = clock();
    for (trial = 0; trial < N_TRIALS; trial++) {
        printf("[Begin trial %d]\n", trial);
        printf("[Reload ANN]\n");
        ai_initann(bi, x);
        printf("[Begin optimization]\n");
        printf("[Result %s]\n",lbfgs_geterror(lbfgs(bi->A, x, x_best, evaluate, NULL, bi,NULL) ));
//        printf("[Result %s]\n",lbfgs_geterror(lbfgs(1000, x, x_best, evaluate, NULL, bi, NULL) ));
//        cg_descent(x, bi->A, NULL, NULL, TOLERANCE, anns_val_f, anns_grad_f_, anns_valgrad_f_, Work, bi);
//        cg_descent(x, bi->A, NULL, NULL, TOLERANCE, anns_val_f, anns_grad_f, anns_valgrad_f, Work, bi);
        printf("[End optimization]\n");
        J = anns_val_f(x, bi->A, bi);
        printf("[Reload bi]\n");
        ai_reload(bi);
        J_after = anns_val_f(x, bi->A, bi);
        printf("[J=%g J_after=%g J_best=%g]\n",J, J_after, J_best);
        if (J_after < J_best) {
            for (i = 0; i < bi->A; i++) {
                x_best[i] = x[i];
            }
            printf("[Update J (was %g, now %g)]\n", J_best, J_after);
            J_best = J_after;
        }
        printf("[End trial %d]\n", trial);
    }

    mse = ai_mse(bi, x_best);
    end = clock();
    dif = 1.0*(end - start)/CLOCKS_PER_SEC;

    printf("[Calculations took %g seconds]\n", dif);
    printf("[Best ANNs after %d trials]\n", N_TRIALS);
    printf("[mse=%g J=%g]\n", mse, anns_val_f(x_best, bi->A, bi));
    for (t = 0; t < bi->nnets; t++) {
        for (i = 0; i < bi->nn[t]; i++) {
            neuron = x_best + bi->I[t] + i*bi->S[t];
            printf("  %g %g %g", neuron[0], neuron[1], neuron[2]);
            for (j = 1; j < bi->nd[t]; j++) {
                printf(" %g", neuron[2+j]);
            }
            printf("\n");
        }
    }

#ifndef NFILEOUT
    printf("[Writing ouput to file \"%s\"]\n", bi->fileAme);

    fout = fopen(bi->fileAme, "w+");
    if (fout) {
        fprintf(fout, "tag=%d type=%d index=%d sub=%d time=%g\n", bi->problemtag, bi->typetag, bi->index, bi->subindex, dif);
        fprintf(fout, "mse=%g J=%g\n", mse, J_best);
        for (t = 0; t < bi->nnets; t++) {
            fprintf(fout, "T=%d nn=%d nd=%d\n", bi->T[t], bi->nn[t], bi->nd[t]);
            for (i = 0; i < bi->nn[t]; i++) {
                neuron = x_best + bi->I[t] + i*bi->S[t];
                fprintf(fout, "%g %g %g", neuron[0], neuron[1], neuron[2]);
                for (j = 1; j < bi->nd[t]; j++) {
                    fprintf(fout, " %g", neuron[2+j]);
                }
                fprintf(fout, "\n");
            }
        }

        printf("[Done]\n");
    } else {
        printf("[Error durning opening file]\n");
    }

    fclose(fout);
#endif

    printf("[Release memory]\n");
    free(x);
    free(x_best);
    free(Work);
    ai_delete(bi);
    printf("[End application]\n");

    return 0;
}

#endif // SOLVE_LBFGS_H_INCLUDED
