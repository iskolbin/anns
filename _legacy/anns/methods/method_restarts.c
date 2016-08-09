#include "method_restarts.h"

int method_restarts(anns_methoddata_t *am) {
    int trial;
    anns_instance_t *ai = NULL;
    double *x = NULL, *x_best = NULL, J_pre, J_best = 1e100, J, J_after, *zw, tmp;
    double *Work = NULL;
    double mse,  *neuron = NULL;
    int i, j, t;
    FILE *fout = NULL, *msefout = fopen("mse.txt", "w+");
    clock_t start, end, diff = 0;
    double dif, asum;
    int counter, times = 1;
    long int seed;
    fmin_callback_args_t args = {NULL, 0, NULL, 0, 0, 0, 0, stdout};

    printf("[Begin application]\n");

    if (am->nofileout) {
        printf("[File output is off]\n");
    }

    seed = am->randomseed ? time(NULL) : am->seed;
    printf("[Init random generator with seed %ld]\n", seed);
    srand(seed);

    printf("[Allocate memory]\n");
    ai = ai_new(am->data, NULL);
    Work = calloc(4*ai->A, sizeof(*Work));
    args.n = ai->A;
    args.instance = ai;
    args.mse = am->ai_mse;
    args.fout = msefout;

    x = ai->z;
    x_best = ai->zbest;
    args.x = x;

    printf("[Initialize ai]\n");
    am->ai_init(ai, am->cond, am->cond_a, am->delta);
    printf("[Reload ai]\n");
    am->ai_newpoints(ai, am->bounds);
    ai->preeval = am->preeval;
//    ac_printfpoints(ai->cache,1);
    ai_printf(ai);

    for (trial = 0; trial < am->maxTrials; trial++) {
        printf("[Begin trial %d]\n", trial);

        printf("[Reload ANN]\n");
        am->ai_initann(ai, x);

        printf("[Initial CSUM]\n");
        anns_csum(x, ai->A, ai);
        asum = 0;
        for (i = 0; i < ai->nconds; i++) {
            asum += ai->cache->csum[i];
        }
        for (i = 0; i < ai->nconds; i++) {
            printf("[%d : %15.15g (%15.15g \%)]\n", i, ai->cache->csum[i], 100. * ai->cache->csum[i] / asum);
        }

        printf("[Begin optimization]\n");
        J_pre = anns_val_f(x, ai->A, ai);

        start = clock();
//        fmin_cg(x, ai->A, anns_val_f, anns_grad_f, anns_valgrad_f, am->gtol, ai, NULL, NULL, NULL);


        cg_descent(x, ai->A, NULL, NULL, am->gtol, anns_val_f, anns_grad_f, anns_valgrad_f, Work, ai, fmin_callback_fprintf_mse, &args);
        end = clock();
        diff += (end-start);
        printf("[End optimization]\n");

        J = anns_val_f(x, ai->A, ai);

        printf("[Reload ai]\n");
        am->ai_newpoints(ai, am->bounds);

        J_after = anns_val_f(x, ai->A, ai);
        printf("[J_pre=%g J=%g J_after=%g J_best=%g]\n",J_pre, J, J_after, J_best);

        if (J_after < J_best) {
            for (i = 0; i < ai->A; i++) {
                x_best[i] = x[i];
            }
            printf("[Update J (was %g, now %g)]\n", J_best, J_after);
            J_best = J_after;
        }

        printf("[End trial %d]\n", trial);
    }

    mse = am->ai_mse(x_best, ai->A, ai);
    dif = 1.0*(end - start)/CLOCKS_PER_SEC;

    printf("[Best ANNs after %d trials]\n", am->maxTrials);
    printf("[mse=%g J=%g Jvg=%g]\n", mse, anns_val_f(x_best, ai->A, ai), anns_valgrad_f(Work, x_best, ai->A, ai));

    printf("[Calculations took %g seconds]\n", 1.0 * diff / (CLOCKS_PER_SEC*times) );
    for (t = 0; t < ai->nnets; t++) {
        for (i = 0; i < ai->nets[t]->nc; i++) {
            neuron = x_best + ai->nets[t]->I + i * ai->nets[t]->nsize;
            printf("%g", neuron[0]);
            for (j = 1; j < ai->nets[t]->nsize; j++) {
                printf(" %g", neuron[j]);
            }
            printf("\n");
        }
    }
    printf("[Eval CSUM]\n");
    anns_csum(x_best, ai->A, ai);
    asum = 0;
    for (i = 0; i < ai->nconds; i++) {
        asum += ai->cache->csum[i];
    }
    for (i = 0; i < ai->nconds; i++) {
        printf("[%d : %15.15g (%15.15g \%)]\n", i, ai->cache->csum[i], 100. * ai->cache->csum[i] / asum);
    }

//    printf("[Refining]\n");
//    for (i = 0; i < ai->A; i++) {
//        x_best[i] = x_best[i];
//    }
//    anns_refine(ai);
//    mse = am->ai_mse(x_best, ai->A, ai);
//    printf("[Best ANNs after %d trials]\n", am->maxTrials);
//    printf("[mse=%g J=%g Jvg=%g]\n", mse, anns_val_f(x_best, ai->A, ai), anns_valgrad_f(Work, x_best, ai->A, ai));

    if (!am->nofileout) {
        printf("[Writing ouput to file \"%s\"]\n", ai->filename);

        fout = fopen(ai->filename, "w+");
        if (fout) {
            fprintf(fout, "tag=%d type=%d index=%d sub=%d time=%g\n", ai->problemtag, ai->typetag, ai->index, ai->subindex, dif);
            fprintf(fout, "mse=%g J=%g\n", mse, J_best);
            for (t = 0; t < ai->nnets; t++) {
                for (i = 0; i < ai->nets[t]->nc; i++) {
                    neuron = x_best + ai->nets[t]->I + i*ai->nets[t]->nsize;

                    fprintf(fout, "%g", neuron[0]);
                    for (j = 1; j < ai->nets[t]->nsize; j++) {
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
    }

    printf("[Release memory]\n");
    free(Work);
    ai_delete(ai);
    printf("[End application]\n");
    fclose(msefout);

    return 0;
}
