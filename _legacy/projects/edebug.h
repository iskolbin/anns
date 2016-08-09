#ifndef EDEBUG_H_INCLUDED
#define EDEBUG_H_INCLUDED

#include <time.h>
#include "../anns/anns.h"

#define RANDSEED 123L

int data[] = {
    PJ_ELLIPTIC_PROBLEM, PJ_DIRECT, 1, 0,

    3, 1, 1, 0,

    0, 0,

    EBF_GAUSS, 2, 3, 0,  2 | MIXED_DERIVATIVES,
};

double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID, 0.0, M_PI,  0.0, M_PI, 0.0, M_PI, // A
};


ANNS_CND(A) {return u_xx[0][0] + u_xx[0][1] + u_xx[0][2] + u_xy[0][0][1] + u_xy[0][1][0] + u_xy[0][2][0] + u_xy[0][2][1];}

ANNS_CND_G(A_a) {return u_axx[0] + u_axx[1] + u_axy[0][1] + u_axy[1][0] + u_axy[2][0] + u_axy[2][1];}


fbvp2_t cond[] = {A};
gbvp2_t cond_a[] = {A_a};

void ai_init(anns_instance_t *ai) {
    ai_conditions(ai, cond, cond_a, NULL);
    ai_loadfuncs(ai);
}

void ai_reload(anns_instance_t *ai) {
    ai_genpoints(ai, bounds);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (t = 0; t < ai->Nt; t++) {
        for (i = ai->I[t]; i < ai->I[t+1]; ) {
            x[i++] = uniform(-4., 4.);         // Веса
            x[i++] = uniform(0.5, 4.);        // Ширины
            x[i++] = uniform(0.5, 4.);        // Ширины
            x[i++] = uniform(0.5, 4.);        // Ширины
            x[i++] = uniform(-M_PI/2,M_PI+M_PI/2);       // Центры
            x[i++] = uniform(-M_PI/2,M_PI+M_PI/2);       // Центры
            x[i++] = uniform(-M_PI/2,M_PI+M_PI/2);       // Центры
        }
    }
}

int solve(int argc, char **argv) {
    int trial;
    anns_instance_t *ai = NULL;
    anns_cache_t *ac = NULL, *acn = NULL;
    double *x = NULL, *g = NULL;
    double s, mse, tmp,  *neuron = NULL;
    int i, j, k, t;
    long int seed = RANDSEED;
    double dif;
    double point[] = {M_PI_2, M_PI_2, M_PI_2};

    printf("[Begin debug application]\n");
    printf("[Init random generator with seed %ld]\n", seed);
    srand(seed);

    printf("[Allocate memory]\n");
    ai = ai_new(data, NULL);
    ac = ai->cache;
    acn = ac_new(ai->Na, ai->Nd, ai->Nc, ai->Nt, ai->Nm, ai->nn, ai->nd, ai->na, ai->np, ai->nv);
    x = ai_createann(ai);
    g = calloc(4*ai->Na, sizeof(*g));

    printf("[Initialize bi]\n");
    ai_init(ai);
    printf("[Reload bi]\n");
    ai_reload(ai);
    ai_initann(ai, x);

    ai->ann_valgrad[0][2](point, x, 0, 1, ai);
    ai->cache = acn;
    srand(seed);
    ai_init(ai);
    ai_reload(ai);
    ai_initann(ai, x);
    anns_numeric_loadfuncs(0, ai);
    ai->ann_valgrad[0][2](point, x, 0, 1, ai);
    printf("[Analytic-Numeric]\n");
    for (t = 0; t < ai->Nt; t++) {
        printf("[Net index %d]\n", t);
        printf("[u = (%.2g-%.2g)=%.2g]\n", ac->u[t], acn->u[t], ac->u[t]-acn->u[t]);
        printf("[u_x =");
        for (j = 0; j < ai->nd[t]; j++) {
            printf("\n  (%.2g-%.2g)=%.2g", ac->u_x[t][j], acn->u_x[t][j], ac->u_x[t][j]-acn->u_x[t][j]);
        }
        printf("\n]\n[u_xx =");
        for (j = 0; j < ai->nd[t]; j++) {
            printf("\n  (%.2g-%.2g)=%.2g", ac->u_xx[t][j], acn->u_xx[t][j], ac->u_xx[t][j]-acn->u_xx[t][j]);
        }

        printf("\n]\n[u_xy =");
        for (j = 0; j < ai->nd[t]; j++) {
            for (k = 0; k < j; k++) {
                printf("\n  (%.2g-%.2g)=%.2g", ac->u_xy[t][j][k], acn->u_xy[t][j][k], ac->u_xy[t][j][k]-acn->u_xy[t][j][k]);
            }
        }

        printf("\n]\n[Gradient]\n[u_a =");
        for (i = 0; i < ai->na[t]; i++) {
            printf("\n  (%.2g-%.2g)=%.2g", ac->u_a[t][i], acn->u_a[t][i], ac->u_a[t][i]-acn->u_a[t][i]);
        }
        printf("\n]\n[u_ax =");
        for (j = 0; j < ai->nd[t]; j++) {
            printf("\n  [%d]:[", j);
            for (i = 0; i < ai->na[t]; i++) {
                printf("\n    (%.2g-%.2g)=%.2g", ac->u_ax[t][i][j], acn->u_ax[t][i][j], ac->u_ax[t][i][j]-acn->u_ax[t][i][j]);
            }
            printf("\n  ]");
        }
        printf("\n]\n[u_axx =");
        for (j = 0; j < ai->nd[t]; j++) {
            printf("\n  [%d]:[", j);
            for (i = 0; i < ai->na[t]; i++) {
                printf("\n    (%.2g-%.2g)=%.2g", ac->u_axx[t][i][j], acn->u_axx[t][i][j], ac->u_axx[t][i][j]-acn->u_axx[t][i][j]);
            }
            printf("\n  ]");
        }
        printf("\n]");

        printf("\n[u_axy =");
        for (j = 0; j < ai->nd[t]; j++) {
            for (k = 0; k < j; k++) {
                printf("\n  [%d,%d]:[", j, k);
                for (i = 0; i < ai->na[t]; i++) {
                    printf("\n    (%.2g-%.2g)=%.2g", ac->u_axy[t][i][j][k], acn->u_axy[t][i][j][k], ac->u_axy[t][i][j][k]-acn->u_axy[t][i][j][k]);
                }
                printf("\n  ]");
            }
        }
        printf("\n]\n");
    }
    printf("[Release memory]\n");
    ac_delete(ac);
    ai_delete(ai);

    free(x);
    free(g);
    printf("[End application]\n");

    return 0;
}

#endif // EDEBUG_H_INCLUDED
