#include "anns_methoddata.h"


static void ai_default_update (anns_instance_t *ai, double *z, int step, double tau);
static double ai_default_mse  (double *, int, void *);
static void ai_default_func   (anns_instance_t *ai);

void ai_default_update (anns_instance_t *ai, double *z, int step, double tau) {}
double ai_default_mse  (double *x, int n, void *instance) {return -1;}
void ai_default_func   (anns_instance_t *ai) {}

anns_methoddata_t *am_new(void) {
    anns_methoddata_t *am = malloc(sizeof *am);

    assert(am);
    am->bounds = NULL;
    am->data = NULL;
    am->cond = NULL;
    am->cond_a = NULL;
    am->cond_i = NULL;
    am->cond_i_a = NULL;
    am->delta = NULL;
    am->fdsteps = 0;

    am->lambda = 0.25;
    am->rho = 0.95;
    am->tol = 1e-4;
    am->gtol = 1e-4;
    am->randomseed = 1;
    am->seed = 0;
    am->maxTrials = 1;
    am->nofileout = 1;

    am->Nopt = 2;
    am->ticksGrow = 2;
    am->ticksReduce = 3;
    am->ticksRegeneration = 10;
    am->ticksTrial = 60;
    am->ticksMax = 100;

    am->n0 = 8;
    am->nmax = 32;

    am->ai_update = ai_default_update;
    am->ai_mse = ai_default_mse;

    am->preeval = NULL;

    return am;
}

anns_methoddata_t *am_new2(void (*ai_init)(anns_instance_t *, fbvp_t *, gbvp_t *, double *),
                            void (*ai_newpoints)(anns_instance_t *, double *),
                            void (*ai_initann)(anns_instance_t *, double *),
                            double (*ai_mse)(double *, int, void *),
                            void (*ai_update)(anns_instance_t *, double *, int , double),
                            int *data, double *bounds, fbvp_t *cond, gbvp_t *cond_a, double *delta) {

    anns_methoddata_t *am = am_new();

    am->ai_init = ai_init;
    am->ai_newpoints = ai_newpoints;
    am->ai_initann = ai_initann;
    am->ai_mse = ai_mse ? ai_mse : ai_default_mse;
    am->ai_update = ai_update ? ai_update : ai_default_update;

    am->data = data;
    am->bounds = bounds;
    am->cond = cond;
    am->cond_a = cond_a;
    am->delta = delta;

    return am;
}

anns_methoddata_t *am_new3(void (*ai_init)(anns_instance_t *, fbvp_t *, gbvp_t *, double *),
                            void (*ai_newpoints)(anns_instance_t *, double *),
                            void (*ai_initann)(anns_instance_t *, double *),
                            double (*ai_mse)(double *, int, void *),
                            void (*ai_update)(anns_instance_t *, double *, int, double ),
                            int *data, double *bounds, fbvp_t *cond, gbvp_t *cond_a, double *delta,
                            fbvp_t *cond_i, gbvp_t *cond_i_a, double T0, double Tm, int nsteps) {

    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints, ai_initann, ai_mse, ai_update, data, bounds, cond, cond_a, delta);
    am->cond_i = cond_i;
    am->cond_i_a = cond_i_a;
    am->T0 = T0;
    am->Tm = Tm;
    am->nsteps = nsteps;
    am->tau = (Tm - T0) / (nsteps-1);

    return am;
}

anns_methoddata_t *am_new_n(anns_methoddata_pointer_t* pointers, anns_methoddata_number_t* numbers) {
    int i = 0;
    anns_methoddata_t *am = am_new();
    double value = 0.0;
    void *pointer = NULL;

    if (pointers) {
        while (pointers->key != AM_NULL) {
            pointer = pointers->value;
            switch (pointers->key) {
                case AM_DATA: am->data = pointer; break;
                case AM_COND: am->cond = pointer; break;
                case AM_COND_A: am->cond_a = pointer; break;
                case AM_DELTA: am->delta = pointer; break;
                case AM_BOUNDS: am->bounds = pointer; break;
                case AM_COND_I: am->cond_i = pointer; break;
                case AM_COND_I_A: am->cond_i_a = pointer; break;
                case AM_UPDATE: am->ai_update = pointer; break;
                case AM_INIT: am->ai_init = pointer; break;
                case AM_INITANN: am->ai_initann = pointer; break;
                case AM_PREEVAL: am->preeval = pointer; break;
                case AM_MSE: am->ai_mse = pointer; break;
                case AM_NEWPOINTS: am->ai_newpoints = pointer; break;
                case AM_UPDATEPOINTS: am->updatepoints = pointer; break;
            }
            pointers++;
        }
    }

    if (numbers) {
        while (numbers->key != AM_NULL) {
            value = numbers->value;
            switch (numbers->key) {
                case AM_SEED:  am->randomseed = 0; am->seed = (int) value; break;
                case AM_TOL: am->tol = value; break;
                case AM_GTOL: am->gtol = value; break;
                case AM_NOFILEOUT: am->nofileout = (int) value; break;
                case AM_OPTSTEPS: am->Nopt = (int) value; break;
                case AM_GROW: am->ticksGrow = (int) value; break;
                case AM_REDUCE: am->ticksReduce = (int) value; break;
                case AM_REGENERATE: am->ticksRegeneration = (int) value; break;
                case AM_TRIAL: am->ticksTrial = (int) value; break;
                case AM_MAXTICKS: am->ticksMax = (int) value; break;
                case AM_INITIAL_N: am->n0 = (int) value; break;
                case AM_MAXIMAL_N: am->nmax = (int) value; break;
                case AM_TBEGIN: am->T0 = value; break;
                case AM_TEND: am->Tm = value; break;
                case AM_TSTEPS: am->fdsteps = (int) value; break;
                case AM_MAXTRIALS: am->maxTrials = (int) value; break;
                case AM_MAXITERS: am->maxIters = (int) value; break;

            }
            numbers++;
        }

        if (am->fdsteps) {
            am->tau = (am->Tm - am->T0) / am->fdsteps;
        }
    }

    return am;
}
