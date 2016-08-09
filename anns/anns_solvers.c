#include "anns_solvers.h"
#include "anns_cache.h"

static void ai_default_update (anns_instance_t *ai, double *z, int step);
static double ai_default_mse  (anns_instance_t *ai, double *z);
static void ai_default_func   (anns_instance_t *ai);

void ai_default_update (anns_instance_t *ai, double *z, int step) {}
double ai_default_mse  (anns_instance_t *ai, double *z) {return -1;}
void ai_default_func   (anns_instance_t *ai) {}

anns_solverdata_t *as_new(void) {
    anns_solverdata_t *as = malloc(sizeof *as);

    assert(as);
    as->bounds = NULL;
    as->data = NULL;
    as->cond = NULL;
    as->cond_a = NULL;
    as->delta = NULL;

    as->lambda = 0.25;
    as->rho = 0.95;
    as->tol = 1e-4;
    as->gtol = 1e-4;
    as->randomseed = 1;
    as->seed = 0;
    as->maxTrials = 1;
//    as->N = 10000;
//    as->Nopt = 64;
//    as->Nremove = 4;
//    as->Add = 1;
//    as->Nrestart = 1;
//    as->Nnewpoints = 16;
    as->nofileout = 0;

    as->Nopt = 2;
    as->ticksGrow = 2;
    as->ticksReduce = 3;
    as->ticksRegeneration = 10;
    as->ticksTrial = 60;
    as->ticksMax = 100;

    as->n0 = 8;
    as->nmax = 32;

    as->ai_update = ai_default_update;
    as->ai_mse = ai_default_mse;

    return as;
}

anns_solverdata_t *as_new2(void (*ai_init)(anns_instance_t *, fbvp_t *, gbvp_t *, double *),
                            void (*ai_newpoints)(anns_instance_t *, double *),
                            void (*ai_initann)(anns_instance_t *, double *),
                            double (*ai_mse)(anns_instance_t *, double *),
                            int *data, double *bounds, fbvp_t *cond, gbvp_t *cond_a, double *delta) {

    anns_solverdata_t *as = as_new();

    as->ai_init = ai_init;
    as->ai_newpoints = ai_newpoints;
    as->ai_initann = ai_initann;
    if (ai_mse) {
        as->ai_mse = ai_mse;
    }

    as->data = data;
    as->bounds = bounds;
    as->cond = cond;
    as->cond_a = cond_a;
    as->delta = delta;

    return as;
}

//static void f

int anns_solve_restarts(anns_solverdata_t *as) {
    int trial;
    anns_instance_t *ai = NULL;
    double *x = NULL, *x_best = NULL, J_pre, J_best = 1e100, J, J_after, *zw, tmp;
		double *Work = NULL;
		double mse,  *neuron = NULL;
		int i, j, t;
		FILE *fout = NULL;
		clock_t start,end, diff = 0;
		double dif;
		int counter, times = 1;
		long int seed;
		cg_parameter uparam = {0};

		cg_default(&uparam);
		//    uparam.PrintFiAl = 0;
		//    uparam.maxit_fac = 0.125;

		printf("[Begin application]\n");
		//
		if (as->nofileout) {
			printf("[File output is off]\n");
		}

		//for (counter = 0; counter < times; counter++) {
		seed = as->randomseed ? time(NULL) : as->seed;
		printf("[Init random generator with seed %ld]\n", seed);
		srand(seed);

		printf("[Allocate memory]\n");
		ai = ai_new(as->data, NULL);
		x = ai->z;
		x_best = ai->zbest;

		Work = calloc(4*ai->A, sizeof(*Work));

		printf("[Initialize ai]\n");
		as->ai_init(ai, as->cond, as->cond_a, as->delta);
		printf("[Reload ai]\n");
		as->ai_newpoints(ai, as->bounds);

		//    for (trial = 0; trial < as->maxTrials; trial++) {
		printf("[Begin trial %d]\n", trial);
		printf("[Reload ANN]\n");
		//
		as->ai_initann(ai, x);
		//
		printf("[Begin optimization]\n");
		J_pre = anns_val_f(x, ai->A, ai);
		////        fmin_cg(x, ai->A, anns_val_f, anns_grad_f, anns_valgrad_f, 1e-7, ai, NULL, NULL, NULL);
		start = clock();
		cg_descent(x, ai->A, NULL, &uparam, as->gtol, anns_val_f, anns_grad_f, anns_valgrad_f, Work, ai, NULL, NULL);
		end = clock();
		diff += (end-start);
		////        fmin_nm(x, ai->A, anns_val_f, TOLERANCE, ai, NULL, NULL, NULL);
		printf("[End optimization]\n");
		J = anns_val_f(x, ai->A, ai);
		printf("[Reload ai]\n");
		//        as->ai_newpoints(ai, as->bounds);
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
		//    }
		//
		mse = as->ai_mse(ai, x_best);
		//
		dif = 1.0*(end - start)/CLOCKS_PER_SEC;


		printf("[Best ANNs after %d trials]\n", as->maxTrials);
		printf("[mse=%g J=%g Jvg=%g]\n", mse, anns_val_f(x_best, ai->A, ai), anns_valgrad_f(Work, x_best, ai->A, ai));

		//}
		//     anns_valgrad_f(ai->cache->g, x, ai->A, ai);
		//    zw = ai_collectweights(ai);
		//printf("[Weights]\n");
		//    for (t = 0; t < ai->nnets; t++) {
		//        for (i = 0; i < ai->nets[t]->nc; i++) {
		//            printf("%g ", zw[i]);
		//        }
		//    }
		//    printf("\n");
		//
		//    tmp = anns_val_f_wonly(zw, ai->Aw, ai);
		//    printf("[Weights f=%g]\n", tmp);
		//
		//    tmp = anns_valgrad_f_wonly(ai->cache->g, zw, ai->Aw, ai);
		//
		//    printf("[Weights grad (f=%g)]\n", tmp);
		//    for (i = 0; i < ai->Aw; i++) {
		//        printf("%g ", ai->cache->g[i]);
		//    }
		//    printf("\n");
		//
		//        tmp = anns_val_f(x, ai->A, ai);
		//    printf("[f=%g]\n", tmp);
		//
		//    tmp = anns_valgrad_f(ai->cache->g, x, ai->A, ai);
		//    printf("[Grad (f=%g)]\n", tmp);
		//    for (i = 0; i < ai->A; i++) {
		//        printf("%g ", ai->cache->g[i]);
		//        if (i && !((i+1) % ai->nets[0]->nsize)) {
		//            printf("\n");
		//        }
		//    }
		//    printf("\n");
		//as->ai_newpoints(ai, as->bounds);

		//    cg_descent(zw, ai->Aw, NULL, &uparam, 1e-6, anns_val_f_wonly, anns_grad_f_wonly, anns_valgrad_f_wonly, Work, ai);
		//
		//     printf("[Weights grad]\n");
		//    for (i = 0; i < ai->Aw; i++) {
		//        printf("%g ", ai->cache->g[i]);
		//    }
		//    printf("\n");
		//   x= ai_expandweights(ai);
		mse = as->ai_mse(ai, x);
		printf("[MSE at the end is %g]\n", mse);
		if ( !as->nofileout ) {
			printf("[Writing points to file %s]", ai->filename );
			fout = fopen(ai->filename, "w+");
			ac_fprintfpoints(fout, ai->cache, 1);
			fclose(fout);
		}
		// printf("[MSE after weights optimization is %g]\n", mse);

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
		for (i = 0; i < ai->nconds; i++) {
			printf("[%d : %15.15f]\n", i, ai->cache->csum[i]);
		}

		if (!as->nofileout) {
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

		return 0;
}

int anns_solve_growreduce(anns_solverdata_t *as) {
     int trial;
    anns_instance_t *ai = NULL, *aic = NULL;
    double *x = NULL, *x_best = NULL, J_pre, J_best = 1e100, J, J_, J_after, J_current, J_pred;
    double *Work = NULL;
    double mse,  *neuron = NULL;
    int i, j, t;
    FILE *fout = NULL;
    clock_t start,end,diff;
    double dif;
    long int seed;
    cg_parameter uparam = {0};
    int tick = 0;

    cg_default(&uparam);
    uparam.maxit_fac = as->Nopt;
//    uparam.maxit_fac = 0.125;

    printf("[Begin application]\n");

    if (as->nofileout) {
        printf("[File output is off]\n");
    }

    seed = as->randomseed ? time(NULL) : as->seed;
    printf("[Init random generator with seed %ld]\n", seed);
    srand(seed);

    printf("[Allocate memory]\n");
    ai = ai_new(as->data, NULL);
    ai->nets[0]->nc = as->n0;
    x = ai->z;
    x_best = ai->zbest;

    Work = calloc(4*ai->A, sizeof(*Work));

    printf("[Initialize ai]\n");
    as->ai_init(ai, as->cond, as->cond_a, as->delta);
    printf("[Reload ai]\n");
    as->ai_newpoints(ai, as->bounds);

    as->ai_initann(ai, x);

    start = clock();
    J_pred = J_best = anns_val_f_check(x, ai->A, ai);
    for (tick = 0; tick < as->ticksMax && J_best > as->tol; tick++) {


        cg_descent(x, ai->A, NULL, &uparam, as->gtol, anns_val_f, anns_grad_f, anns_valgrad_f, Work, ai, NULL, NULL);
        J_pred = J_current;
        J_current = anns_val_f_check(x, ai->A, ai);
        printf("\nJ_ = %g\n", J_current);
        if (J_current < J_best) {
            J_best = J_current;
            ai_savez(ai);
//            for (i = 0; i < ai->A; i++) {
//                x_best[i] = x[i];
//            }
        } else {
//                    if (ai->nets[0]->nc < ai->nets[0]->nm) {
//            ai->nets[0]->nc++;
//            ai->nets[0]->nctmp++;
//            }
//            as->ai_newpoints(ai, as->bounds);
        }

        if (J_current > as->rho * J_pred) {
            as->ai_newpoints(ai, as->bounds);
            printf("\n New points \n");
        }


    }
//
//    for (trial = 0; trial < as->Nrestart; trial++) {
//        printf("[Begin trial %d]\n", trial);
//        printf("[Reload ANN]\n");
//
//        as->ai_initann(ai, x);
//
//        printf("[Begin optimization]\n");
//        J_pre = anns_val_f(x, ai->A, ai);
////        fmin_cg(x, ai->A, anns_val_f, anns_grad_f, anns_valgrad_f, 1e-7, ai, NULL, NULL, NULL);
//        cg_descent(x, ai->A, NULL, &uparam, as->tol, anns_val_f, anns_grad_f, anns_valgrad_f, Work, ai);
////        fmin_nm(x, ai->A, anns_val_f, TOLERANCE, ai, NULL, NULL, NULL);
//        printf("[End optimization]\n");
//        J = anns_val_f(x, ai->A, ai);
//        printf("[Reload ai]\n");
//        as->ai_newpoints(ai, as->bounds);
//        J_after = anns_val_f(x, ai->A, ai);
//        printf("[J_pre=%g J=%g J_after=%g J_best=%g]\n",J_pre, J, J_after, J_best);
//        if (J_after < J_best) {
//            for (i = 0; i < ai->A; i++) {
//                x_best[i] = x[i];
//            }
//            printf("[Update J (was %g, now %g)]\n", J_best, J_after);
//            J_best = J_after;
//        }
//        printf("[End trial %d]\n", trial);
//    }

    mse = as->ai_mse(ai, x_best);
    end = clock();
    dif = 1.0*(end - start)/CLOCKS_PER_SEC;

    printf("[Calculations took %g seconds]\n", dif);
    printf("[Best ANNs after %d trials]\n", as->maxTrials);
    ai_switchz(ai);
    J = anns_val_f(x, ai->A, ai);
    J_ = anns_val_f_check(x, ai->A, ai);
    printf("[mse=%g J=%g J_=%g K=%g]\n", mse, J, J_, J_ / J);
    ai_switchz(ai);

//    ai_printf(ai);
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
    for (i = 0; i < ai->nconds; i++) {
        printf("[%d : %15.15f]\n", i, ai->cache->csum[i]);
    }

    if (!as->nofileout) {
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
    free(x);
    free(x_best);
    free(Work);
    ai_delete(ai);
    printf("[End application]\n");

    return 0;
}

int anns_solve_grow(anns_solverdata_t *as) {
    int trial;
    anns_instance_t *ai = NULL;
    double *x = NULL, *x_best = NULL, J_pre, J_best = 1e100, J, J_after, *zw, tmp;
    double *Work = NULL;
    double mse,  *neuron = NULL;
    int i, j, t;
    FILE *fout = NULL;
    clock_t start,end, diff = 0;
    double dif;
    int counter, times = 1;
    long int seed;
    cg_parameter uparam = {0};

    cg_default(&uparam);

    printf("[Begin application]\n");

    if (as->nofileout) {
        printf("[File output is off]\n");
    }

    seed = as->randomseed ? time(NULL) : as->seed;
    printf("[Init random generator with seed %ld]\n", seed);
    srand(seed);

    printf("[Allocate memory]\n");
    ai = ai_new(as->data, NULL);
    ai->nets[0]->nc = as->n0;
    x = ai->z;
    x_best = ai->zbest;

    Work = calloc(4*ai->A, sizeof(*Work));

    printf("[Initialize ai]\n");
    as->ai_init(ai, as->cond, as->cond_a, as->delta);
    as->ai_newpoints(ai, as->bounds);

    as->ai_initann(ai, x);

        printf("[Begin optimization]\n");
        J_pre = anns_val_f(x, ai->A, ai);
        start = clock();
        cg_descent(x, ai->A, NULL, &uparam, as->gtol, anns_val_f, anns_grad_f, anns_valgrad_f, Work, ai, NULL, NULL);
        end = clock();
        diff += (end-start);
        printf("[End optimization]\n");
        J = anns_val_f(x, ai->A, ai);
////        printf("[Reload ai]\n");
        as->ai_newpoints(ai, as->bounds);
        J_after = anns_val_f(x, ai->A, ai);
        printf("[J_pre=%g J=%g J_after=%g J_best=%g]\n",J_pre, J, J_after, J_best);
        if (J_after < J_best) {
            for (i = 0; i < ai->A; i++) {
                x_best[i] = x[i];
            }
            J_best = J_after;
        }

    mse = as->ai_mse(ai, x);
    printf("[MSE at the end is %g]\n", mse);
    fout = fopen(ai->filename, "w+");

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
}

#define DBL_FORMAT "%15.15f"

int anns_solve_hybrid(anns_solverdata_t *as) {
  int trial;
    anns_instance_t *ai = NULL;
    double *x = NULL, *x_best = NULL, *x_ = NULL, J_best = 1e100, J, J_after;
    double *Work = NULL;
    double s, mse, tmp,  *neuron = NULL;
    int i, j, k, p, t;
//    int result = SOLVE_OK;
    long int seed = time(NULL);

    FILE *clog = stdout;
    FILE *fout = NULL;

    fprintf(clog, "[Begin application]\n");
    fprintf(clog, "[Finite difference time split]\n");
#ifdef NFILEOUT
    fprintf(clog, "[File output is off]\n");
#endif
#ifndef NDEBUG
    fprintf(clog, "[Debug is on]\n");
#else
    fprintf(clog, "[Debug is off]\n");
#endif
    fprintf(clog, "[Init random generator with seed %ld]\n", seed);
    srand(seed);


    seed = as->randomseed ? time(NULL) : as->seed;
//    printf("[Init random generator with seed %ld]\n", seed);
    srand(seed);

//    printf("[Allocate memory]\n");
    ai = ai_new(as->data, NULL);



    fprintf(clog, "[Allocate memory]\n");
    ai = ai_new(as->data, NULL);
    x = ai->z;
    Work = calloc(4*ai->A, sizeof(*Work));

    fprintf(clog, "[Initialize bi]\n");
    as->ai_init(ai, as->cond, as->cond_a, as->delta);
    as->ai_initann(ai, x);

    for (i = 1; i < as->fdsteps; i++) {
        fprintf(clog, "[Begin step %d of %d]\n", i, as->fdsteps-1);
        as->ai_update(ai, x, i);
        cg_descent(x, ai->A, NULL, NULL, as->gtol, anns_val_f, anns_grad_f, anns_valgrad_f, Work, ai, NULL, NULL);
        fprintf(clog, "[End step %d of %d]\n", i, as->fdsteps-1);
    }
    ac_put(ai->cache, x, as->fdsteps-1);

#ifndef NFILEOUT
        fprintf(clog, "[Writing ouput to file \"%s\"]\n", ai->filename);

        fout = fopen(ai->filename, "w+");
        if (fout) {
            fprintf(fout, "tag=%d type=%d index=%d sub=%d nsteps=%d\n", ai->problemtag, ai->typetag, ai->index, ai->subindex, as->fdsteps);
            for (k = 0; k < as->fdsteps; k++) {
                fprintf(fout, "k=%d\n", k);
                for (t = 0; t < ai->nnets; t++) {
                    fprintf(fout, "T=%d nn=%d nd=%d\n", ai->nets[t]->T, ai->nets[t]->nc, ai->nets[t]->dim);
//                    x_ = ac_get(ai->cache, k);
                    x_ = ai->cache->xm[k];
                    for (i = 0; i < ai->nets[t]->nc; i++) {
                        neuron = x_ + ai->nets[t]->I + i*ai->nets[t]->nsize;
                        fprintf(fout, DBL_FORMAT " " DBL_FORMAT " " DBL_FORMAT, neuron[0], neuron[1], neuron[2]);
                        for (j = 1; j < ai->nets[t]->dim; j++) {
                            fprintf(fout, " " DBL_FORMAT, neuron[2+j]);
                        }
                        fprintf(fout, "\n");
                    }
                }
            }
            fprintf(clog, "[Done]\n");
        } else {
            fprintf(clog, "[Error durning opening file]\n");
//            result = SOLVE_WRITE_ERROR;
        }

        fclose(fout);
#endif

    fprintf(clog, "[Release memory]\n");
//    free(x);
//    free(x_best);
    free(Work);
    ai_delete(ai);
    fprintf(clog, "[End application]\n");

    return 0;
}
