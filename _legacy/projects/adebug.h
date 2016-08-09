#ifndef ADEBUG_H_INCLUDED
#define ADEBUG_H_INCLUDED

#include <time.h>
#include "../anns/anns.h"

#define RANDSEED 11123L

int data[] = {
    PJ_ELLIPTIC_PROBLEM, PJ_DIRECT, 1, 0,

    FUNCTIOAL_SQUARE, 3, 1, 2, 0,

    1, 0,

    CLASSIC | RADIAL | CAUCHY_FUNCTION, 2, 3, 0,  2 | MIXED_DERIVATIVES,
    CLASSIC | RADIAL | CAUCHY_FUNCTION, 2, 3, 0,  2 | MIXED_DERIVATIVES,
};

double bounds[] = {
    GP_POINTS, M_PI_2, M_PI_2, M_PI_2, // A
};


ANNS_CND(A) {return u_xx[0][0] + u_xx[0][1] + u_xx[0][2] +
    u_xy[0][1][0] + u_xy[0][2][0] + u_xy[0][2][1] +
    u_xx[1][0] + u_xx[1][1] + u_xx[1][2] +
    u_xy[1][1][0] + u_xy[1][2][0] + u_xy[1][2][1];}

ANNS_CND_G(A_a) {return u_axx[0] + u_axx[1] + u_axx[2] + u_axy[1][0] + u_axy[2][0] + u_axy[2][1];}


fbvp_t cond[] = {A};
gbvp_t cond_a[] = {A_a, A_a};

void ai_init(anns_instance_t *ai) {
    ai_conditions(ai, cond, cond_a, NULL);
    ai_loadfuncs(ai);
}

void ai_reload(anns_instance_t *ai) {
    ai_genpoints(ai, bounds);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (t = 0; t < ai->nnets; t++) {
        for (i = ai->I[t]; i < ai->I[t+1]; ) {
            x[i++] = uniform(-4., 4.);         // Веса
            x[i++] = uniform(0.5, 4.);        // Ширины
            x[i++] = uniform(-M_PI/2,M_PI+M_PI/2);       // Центры
            x[i++] = uniform(-M_PI/2,M_PI+M_PI/2);       // Центры
            x[i++] = uniform(-M_PI/2,M_PI+M_PI/2);       // Центры
        }
    }
}


double ai_mse(anns_instance_t *ai, double *x) {
    double s, tmp, point[2];
    int i, p, M = 10000;

    s = 0;
    for (p = 0; p < M; p++) {
        point[0] = uniform(0, M_PI);
        point[1] = uniform(0, M_PI);

        tmp = 0;
        tmp += anns_eval(point, x, ai);
        tmp += 0.5*sin(point[0])*sin(point[1]);

        s += tmp * tmp;
    }

    return sqrt(s / (M-1) );
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
    acn = ac_new(ai->A, ai->dim, ai->nconds, ai->nnets, ai->nmem, ai->nn, ai->nd, ai->A, ai->np, ai->nv);
    x = ai_createann(ai);
    g = calloc(4*ai->A, sizeof(*g));

    t = 0;
    printf("[AAlytic-Numeric : AAlytic/Numeric-1]\n");
    for (t = 0; t < ai->nnets; t++) {
        printf("[Initialize bi]\n");
        srand(seed);
        ai_init(ai);
        printf("[Reload bi]\n");
        ai_reload(ai);
        ai_initann(ai, x);
        ai->cache = ac;
        ai->ann_valgrad[t][2](point, x, t, 1, ai);
        ai->cache = acn;
        srand(seed);
        ai_init(ai);
        ai_reload(ai);
        ai_initann(ai, x);
        anns_numeric_loadfuncs(t, ai);
        ai->ann_valgrad[t][2](point, x, t, 1, ai);
        anns_numeric_freefuncs(t, ai);
        printf("[Net index %d]\n", t);
        printf("[u = (%.2g-%.2g)=%.2g : %.2g]\n", ac->u[t], acn->u[t], ac->u[t]-acn->u[t], ac->u[t]/acn->u[t]-1);
        printf("[u_x =");
        for (j = 0; j < ai->nd[t]; j++) {
            printf("\n  (%.2g-%.2g)=%.2g : %.2g", ac->u_x[t][j], acn->u_x[t][j], ac->u_x[t][j]-acn->u_x[t][j], ac->u_x[t][j]/acn->u_x[t][j]-1);
        }
        printf("\n]\n[u_xx =");
        for (j = 0; j < ai->nd[t]; j++) {
            printf("\n  (%.2g-%.2g)=%.2g : %.2g", ac->u_xx[t][j], acn->u_xx[t][j], ac->u_xx[t][j]-acn->u_xx[t][j], ac->u_xx[t][j]/acn->u_xx[t][j]-1);
        }

        printf("\n]\n[u_xy =");
        for (j = 0; j < ai->nd[t]; j++) {
            for (k = 0; k < j; k++) {
                printf("\n  (%.2g-%.2g)=%.2g : %.2g", ac->u_xy[t][j][k], acn->u_xy[t][j][k], ac->u_xy[t][j][k]-acn->u_xy[t][j][k], ac->u_xy[t][j][k]/acn->u_xy[t][j][k]-1);
            }
        }

        printf("\n]\n[Gradient]\n[u_a =");
        for (i = 0; i < ai->A[t]; i++) {
            printf("\n  (%.2g-%.2g)=%.2g : %.2g", ac->u_a[t][i], acn->u_a[t][i], ac->u_a[t][i]-acn->u_a[t][i], ac->u_a[t][i]/acn->u_a[t][i]-1);
        }
        printf("\n]\n[u_ax =");
        for (j = 0; j < ai->nd[t]; j++) {
            printf("\n  [%d]:[", j);
            for (i = 0; i < ai->A[t]; i++) {
                printf("\n    (%.2g-%.2g)=%.2g : %.2g", ac->u_ax[t][i][j], acn->u_ax[t][i][j], ac->u_ax[t][i][j]-acn->u_ax[t][i][j], ac->u_ax[t][i][j]/acn->u_ax[t][i][j]-1);
            }
            printf("\n  ]");
        }
        printf("\n]\n[u_axx =");
        for (j = 0; j < ai->nd[t]; j++) {
            printf("\n  [%d]:[", j);
            for (i = 0; i < ai->A[t]; i++) {
                printf("\n    (%.2g-%.2g)=%.2g : %.2g", ac->u_axx[t][i][j], acn->u_axx[t][i][j], ac->u_axx[t][i][j]-acn->u_axx[t][i][j], ac->u_axx[t][i][j]/acn->u_axx[t][i][j]-1);
            }
            printf("\n  ]");
        }
        printf("\n]");

        printf("\n[u_axy =");
        for (j = 0; j < ai->nd[t]; j++) {
            for (k = 0; k < j; k++) {
                printf("\n  [%d,%d]:[", j, k);
                for (i = 0; i < ai->A[t]; i++) {
                    printf("\n    (%.2g-%.2g)=%.2g : %.2g", ac->u_axy[t][i][j][k], acn->u_axy[t][i][j][k], ac->u_axy[t][i][j][k]-acn->u_axy[t][i][j][k], ac->u_axy[t][i][j][k]/acn->u_axy[t][i][j][k]-1);
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

#endif // ADEBUG_H_INCLUDED
