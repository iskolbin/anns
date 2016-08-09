#include "stefan.h"

static void ai_init(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
static void ai_newpoints(anns_instance_t *, double *bounds);
static void ai_initann(anns_instance_t *, double *);
static double ai_mse(anns_instance_t *, double *);

#define C_C    1.0
#define C_K    1.0
#define C_lambda 1.0
#define C_To   0.0
#define C_Tm   C_To
#define C_Th   -1.0
#define C_xmax 5.0

ANNS_CND(A1d_i) {return u[0] - C_To;}
ANNS_CND_G(A1d_i_a) {return u_a;}

//ANNS_CND(A1d) {anns_instance_t *ai = instance; return x[0] < u[1] ? (0.5 * (ai->step) * C_K * u_xx[0][0] - u[0] + v[0]) : 0;}
ANNS_CND(A1d) {anns_instance_t *ai = instance; return x[0] < u[1] ? (0.5 * (ai->step) * (C_K/C_C) * u_xx[0][0] - u[0] + v[0]) : 0;}
ANNS_CND(B1d_left) {return u[0] - C_Th;}
ANNS_CND(B1d_front) {return u[0] - C_Tm;}
ANNS_CND(B1d_phase) {anns_instance_t *ai = instance; return (C_K*ai->step/C_lambda) * u_x[0][0] - u[1] + v[0];}

ANNS_CND_G(A1d_a_1) {anns_instance_t *ai = instance; return x[0] < u[1] ? (0.5 * (ai->step) * (C_K/C_C) * u_axx[0]- u_a) : 0;}
ANNS_CND_G(A1d_a_2) {return 0;}
ANNS_CND_G(B1d_a) {return u_a;}
ANNS_CND_G(B1d_phase_a_1) {anns_instance_t *ai = instance; return (C_K*ai->step/C_lambda) * u_ax[0]; }
ANNS_CND_G(B1d_phase_a_2) {anns_instance_t *ai = instance; return -u_a; }

ANNS_CND_G(F_zero) {return 0; }



static void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);

    ai->cache->v[3][0][0] = 0;
}

static void ai_newpoints_1d(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
}


static void ai_initann_1d(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;
    double b = 1/sqrt(ai->nets[0]->nc);

    // Сеть, приближающая решение
    for (i = 0; i < ai->nets[0]->lenc;) {
        x[i++] = uniform(C_To, C_To);   // Веса
        x[i++] = b;
        x[i++] = uniform(-0.05, C_xmax + 0.05);
    }

    // Сеть, приближающая точку фазового перехода
    x[i++] = 0;
    x[i++] = 1;
    x[i++] = 1;
}


static void ai_update_1d(anns_instance_t *ai, double *x, int i, double tau) {
    int p;
    anns_cache_t *ac = ai->cache;
    double point;

    // Сохраняем сеть i-1 шага в контейнер
    ai_put(ai, x, i-1);
//    s_ = s;

    // Обновляем предысторию (v[0]) для конечно-разностной схемы (явно-неявная, Кранк-Никлсон)
    for (p = 0; p < ai->conds[0]->m; p++) {
        ai->nets[0]->ann_val[2](x, 0, p, 0, 0, ai);
        point = ai->cache->points[0][p][0];
        ac->v[0][p][0] = 0.5*tau*(point < ac->u[1] ? C_K : 0) *ac->u_xx[0][0] + ac->u[0];
    }
    ac->v[3][0][0] = ac->u[1];
}

//#define C_near 1e-6

static void ai_preeval_1d(double *z, void *anns_instance) {
//    double u_x_left, u_x_right;

    anns_instance_t *ai = anns_instance;
//net->ann_val[eval % MIXED_DERIVATIVES](x, c, p, t, eval & MIXED_DERIVATIVES, ai);
//    ai->cache->points[2][0][0] -= C_near;
//    ai->cache->nets[0]->ann_val[1](z, 2, 0, 0, 0, ai);
//    u_x_left = ai->cache->u_x[0][0];
//
//    ai->cache->points[2][0][0] += 2*C_near;
//    ai->cache->nets[0]->ann_val[1](z, 2, 0, 0, 0, ai);
//    u_x_right = ai->cache->u_x[0][0];
//(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);

//    ai->cache->points[2][0][0];
//    printf("->s=%f s_=%f p=%f ux=%f\n ", s, s_, ai->cache->points[2][0][0], ai->cache->u_x[0][0]);
//    s = s_ + (ai->step / C_lambda) * (C_k1*u_x_left - C_K2*u_x_right);
//    printf("s=%g ", s);
//    printf("<-s=%f s_=%f p=%f ux=%f\n ", s, s_, ai->cache->points[2][0][0], ai->cache->u_x[0][0]);
    ai->cache->nets[1]->ann_val[0](z, 3, 0, 1, 0, ai);
    ai->cache->points[2][0][0] = ai->cache->u[1];
    ai->cache->points[3][0][0] = ai->cache->u[1];
}


//void stefan_1d(double Tm, int fdsteps) {
//
//}
void stefan1phase(double Tm, int fdsteps) {
    int na = 30;

    int data[] = {
        PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 5, 1,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 1, 4, 2, fdsteps,
        na, na,  1, /* Уравнение переноса */
        1, 1,   0, /* Левая граница */
        1, 1,   0, /* Граница фазового перехода */
        1, 1,   1, /* Граница фазового перехода  */
        NORMALIZED | RADIAL | GAUSSIAN, 4, 1, 0,  2, 1, 0, 0,
        NORMALIZED | RADIAL | GAUSSIAN, 1, 1, 0,  0, -1, -1, 0,};

    double bounds[] = {
        GP_GRID_1D, na, 0, C_xmax,
        GP_POINTS, 0,
        GP_POINTS, 0,
        GP_POINTS, 0,
    };

    fbvp_t cond[] = {A1d, B1d_left, B1d_front, B1d_phase };
    gbvp_t cond_a[] = {A1d_a_1, A1d_a_2,
        B1d_a, NULL,
        B1d_a, NULL,
        B1d_phase_a_1, B1d_phase_a_2};

    fbvp_t cond_i[] = {A1d_i, A1d_i, A1d_i, A1d_i};
    gbvp_t cond_i_a[] = {A1d_i_a, F_zero,
        A1d_i_a, F_zero,
        A1d_i_a, F_zero,
        F_zero, F_zero};

    anns_methoddata_t *am = am_new3(ai_init, ai_newpoints_1d, ai_initann_1d, NULL, ai_update_1d,
                                    data, bounds,
                                    cond, cond_a, NULL,
                                    cond_i, cond_i_a, 0, Tm, fdsteps);

    am->randomseed = 0;
    am->seed = 123L;
    am->Nopt = 5;
    am->gtol = 1e-5;
    am->tol = 1e-5;
    am->n0 = 4;
    am->nofileout = 0;
    am->preeval = ai_preeval_1d;

    method_hybrid(am);

    am_delete(am);
}
