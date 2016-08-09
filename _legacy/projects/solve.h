#ifndef ANNS_SCHEME
#define ANNS_SCHEME

#include <time.h>
#include "../anns/anns.h"
#include "../optimize/fmin_cg_descent.h"
#include "../optimize/fmin_nm.h"
#include "../optimize/fmin_hj.h"
#include "../optimize/fmin_cg.h"

#include "dir_ell_1.h"
//#include "dir_ell_1_e.h"
//#include "dir_ell_1_2nets.h"
//#include "dir_par_1.h"
//#include "dir_ell_2.h"
//#include "dir_ell_3.h"
//#include "dir_ell_4.h"
//#include "dir_ell_4.h"
//#include "dir_ell_6.h"
//#include "dir_ell_7.h"
//#include "dir_ell_8.h"
//#include "dir_ell_9.h"
//#include "dir_ell_9_1.h"
//#include "dir_ell_10.h"
//#include "dir_ell_11.h"
//#include "dir_ell_12.h"
//#include "dir_par_4.h"
//#include "dir_par_5.h"
//#include "ibnd_par_1.h"
//#include "ibnd_par_2.h"
//#include "ibnd_par_3.h"
//#include "ibnd_par_4.h"
//#include "ibnd_par_4_err01.h"
//#include "ibnd_par_4_2.h"
//#include "appr_2.h"
//#include "isrc_ell_1.h"
//#include "dir_par_6.h"
//#include "isrc_ell_2.h"
//#include "isrc_ell_2_1.h"

#define N_TRIALS 1
#define NFILEOUT
//#define RANDSEED (time(NULL))
#define RANDSEED 123L

#define TOLERANCE 1e-5

int solve(int argc, char **argv) {
    int trial;
    anns_instance_t *ai = NULL;
    double *x = NULL, *x_best = NULL, J_best = 1e100, J, J_after;
    double *Work = NULL;
    double mse,  *neuron = NULL;
    int i, j, t;
    long int seed = RANDSEED;
    FILE *fout = NULL;
    clock_t start,end;
    double dif;

    printf("[Begin application]\n");
#ifdef NFILEOUT
    printf("[File output is off]\n");
#endif
    printf("[Init random generator with seed %ld]\n", seed);
    srand(seed);

    printf("[Allocate memory]\n");
    ai = ai_new(data, NULL);
    x = ai_createann(ai);
    x_best = ai_createann(ai);
    Work = calloc(4*ai->A, sizeof(*Work));

    printf("[Initialize ai]\n");
    ai_init(ai);
    printf("[Reload ai]\n");
    ai_reload(ai);

    start = clock();
    for (trial = 0; trial < N_TRIALS; trial++) {
        printf("[Begin trial %d]\n", trial);
        printf("[Reload ANN]\n");
        srand(seed);
        ai_initann(ai, x);
//        ac_v_fout(stdout, ai->cache);
        printf("[Begin optimization]\n");
//        fmin_cg(x, ai->A, anns_val_f, anns_grad_f, anns_valgrad_f, 1e-7, ai, NULL, NULL, NULL);
        cg_descent(x, ai->A, NULL, NULL, TOLERANCE, anns_val_f, anns_grad_f, anns_valgrad_f, Work, ai);
//        fmin_nm(x, ai->A, anns_val_f, TOLERANCE, ai, NULL, NULL, NULL);
        printf("[End optimization]\n");
        J = anns_val_f(x, ai->A, ai);
        printf("[Reload ai]\n");
        ai_reload(ai);
        J_after = anns_val_f(x, ai->A, ai);
        printf("[J=%g J_after=%g J_best=%g]\n",J, J_after, J_best);
        if (J_after < J_best) {
            for (i = 0; i < ai->A; i++) {
                x_best[i] = x[i];
            }
            printf("[Update J (was %g, now %g)]\n", J_best, J_after);
            J_best = J_after;
        }
        printf("[End trial %d]\n", trial);
    }

    mse = ai_mse(ai, x_best);
    end = clock();
    dif = 1.0*(end - start)/CLOCKS_PER_SEC;

    printf("[Calculations took %g seconds]\n", dif);
    printf("[Best ANNs after %d trials]\n", N_TRIALS);
    printf("[mse=%g J=%g]\n", mse, anns_val_f(x_best, ai->A, ai));
    for (t = 0; t < ai->nnets; t++) {
        for (i = 0; i < ai->nn[t]; i++) {
            neuron = x_best + ai->I[t] + i*ai->S[t];
            printf("%g", neuron[0]);
            for (j = 1; j < ai->S[t]; j++) {
                printf(" %g", neuron[j]);
            }
            printf("\n");
        }
    }
    printf("[Eval CSUM]\n");
    anns_csum(x_best, ai->A, ai);
    for (i = 0; i < ai->nconds; i++) {
        printf("[%d : %15.15f]\n", i, ai->cache->csum[i]);
    }
#ifndef NFILEOUT
    printf("[Writing ouput to file \"%s\"]\n", ai->fileAme);

    fout = fopen(ai->fileAme, "w+");
    if (fout) {
        fprintf(fout, "tag=%d type=%d index=%d sub=%d time=%g\n", ai->problemtag, ai->typetag, ai->index, ai->subindex, dif);
        fprintf(fout, "mse=%g J=%g\n", mse, J_best);
        for (t = 0; t < ai->nnets; t++) {
            fprintf(fout, "T=%d nn=%d nd=%d\n", ai->T[t], ai->nn[t], ai->nd[t]);
            for (i = 0; i < ai->nn[t]; i++) {
                neuron = x_best + ai->I[t] + i*ai->S[t];
                fprintf(fout, "%g", neuron[0]);
                for (j = 1; j < ai->S[t]; j++) {
                    fprintf(fout, " %g", neuron[j]);
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
    ai_delete(ai);
    printf("[End application]\n");

    return 0;
}

#endif // ANNS_SCHEME
