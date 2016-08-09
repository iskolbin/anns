#include "method_hybrid.h"


static void printf_callback(double *x, int n, void *ai, int iteration, double eps, void*inst) {
    printf("i=%d eps=%g\n", iteration, eps);
}

int method_hybrid(anns_methoddata_t *am) {
    int trial;
    anns_instance_t *ai = NULL;
    double *x = NULL, *x_best = NULL, J_pre, J_best = 1e100, J, J_after, *zw, tmp;
    double *Work = NULL;
    double mse,  *neuron = NULL;
    int i, j, t, k;
    FILE *fout = NULL;
    clock_t start, end, diff = 0;
    double dif, asum, J_mean = 0;
    int counter, times = 1;
    long int seed;

    printf("[Begin application]\n");

    if (am->nofileout) {
        printf("[File output is off]\n");
    }

    seed = am->randomseed ? time(NULL) : am->seed;
    printf("[Init random generator with seed %ld]\n", seed);
    srand(seed);

    printf("[Allocate memory]\n");
    ai = ai_new(am->data, &am->tau);
    Work = calloc(4*ai->A, sizeof(*Work));

    x = ai->z;
    x_best = ai->zbest;

    printf("[Reload ANN]\n");
    am->ai_initann(ai, x);

    printf("[Reload ai]\n");
    am->ai_newpoints(ai, am->bounds);

    if (am->cond_i && am->cond_i_a) {
        printf("[Initialize ai with approximation conditions]\n");
        am->ai_init(ai, am->cond_i, am->cond_i_a, am->delta);

        printf("[Begin approximation of initial condition]\n");
        ai_printf(ai);

        start = clock();
//            fmin_cg(x, ai->A, anns_val_f, anns_grad_f, anns_valgrad_f, am->gtol, ai, NULL, NULL, NULL);
        cg_descent(x, ai->A, NULL, NULL, am->gtol, anns_val_f, anns_grad_f, anns_valgrad_f, Work, ai);
        J_mean += anns_val_f(x, ai->A, ai);
        printf("[End approximation of initial condition]\n");
    } else {
        printf("[Skip approximation of initial condition]\n");
    }

    am->ai_init(ai, am->cond, am->cond_a, am->delta);
    ai->step = am->tau;

    printf("[Begin finite difference steps]\n");
    ai->preeval = am->preeval;

    printf("[Dump ai_instance]\n");
    ai_printf(ai);
    printf("[Points]\n");
    ac_printfpoints(ai->cache,0);

    for (i = 1; i < ai->nsteps; i++) {
        am->ai_update(ai, x, i, am->tau);

//        cg_descent(x, ai->A, NULL, NULL, am->gtol, anns_val_f, anns_grad_f, anns_valgrad_f, Work, ai);
        fmin_hj(anns_val_f, x, ai->A, x,ai,NULL, printf_callback, NULL );
        J_mean += anns_val_f(x, ai->A, ai);

        printf("[Eval CSUM]\n");
        anns_csum(x, ai->A, ai);
        asum = 0;
        for (j = 0; j < ai->nconds; j++) {
            asum += ai->cache->csum[j];
        }
        for (j = 0; j < ai->nconds; j++) {
            printf("[%d : %15.15g (%15.15g \%)]\n", i, ai->cache->csum[j], 100. * ai->cache->csum[j] / asum);
        }

    }
    J_mean /= ai->nsteps;
    ai_put(ai, x, ai->nsteps-1);
    printf("[End finite difference steps]\n");

    end = clock();
    mse = am->ai_mse(x_best, ai->A, ai );
    dif = (1.0*(end - start))/CLOCKS_PER_SEC;

    printf("[Calculations took %g seconds (%d/%d)]\n", dif,  start, end);
    printf("[Mean J = %15.15g]\n", J_mean);


    if (!am->nofileout) {
        printf("[Writing ouput to file \"%s\"]\n", ai->filename);

        fout = fopen(ai->filename, "w+");
        if (fout) {
            fprintf(fout, "tag=%d type=%d index=%d sub=%d time=%g nsteps=%d tau=%.15g Jmean=%.15g\n", ai->problemtag, ai->typetag, ai->index, ai->subindex, dif, am->nsteps, am->tau, J_mean);
            for (k= 0; k < ai->nmem; k++) {
                fprintf(fout, "k=%d\n", k);
                for (t = 0; t < ai->nnets; t++) {
                    fprintf(fout, "T=%d n=%d dim=%d\n", ai->nets[t]->T, ai->nets[t]->nc, ai->nets[t]->dim);
                    for (i = 0; i < ai->nets[t]->nc; i++) {
                        neuron = ai->cache->xm[k] + i*ai->nets[t]->nsize + ai->nets[t]->I ;
                        fprintf(fout, "%g", neuron[0]);
                        for (j = 1; j < ai->nets[t]->nsize; j++) {
                            fprintf(fout, " %g", neuron[j]);
                        }
                        fprintf(fout, "\n");
                    }
                }
            }

            printf("[Done]\n");
        } else {
            printf("[Error durning opening file]\n");
        }

        fclose(fout);
    }

    printf("[Release memory]\n");
    free(Work);
    ai_delete(ai);
    printf("[End application]\n");

    return 0;
}
