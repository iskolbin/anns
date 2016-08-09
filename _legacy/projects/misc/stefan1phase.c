#include "stefan.h"

static void ai_init(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
static void ai_newpoints(anns_instance_t *, double *bounds);
static void ai_initann(anns_instance_t *, double *);
static double ai_mse(anns_instance_t *, double *);

static double s_pred = 0;
static double s_current = 0;

#define C_C      1.0
#define C_K      1.0
#define C_lambda 1.0
#define C_To     0.0
#define C_Tm     C_To
#define C_Th     -1.0
#define C_xmax   3.0

ANNS_CND(A1d_i) {return u[0] - C_To;}
ANNS_CND_G(A1d_i_a) {return u_a;}

//ANNS_CND(A1d) {anns_instance_t *ai = instance; return x[0] < u[1] ? (0.5 * (ai->step) * C_K * u_xx[0][0] - u[0] + v[0]) : 0;}
ANNS_CND(A1d) {
    anns_instance_t *ai = instance;
    double tau = ai->step;
return x[0] < s_current ? tau * (C_K/C_C) * u_xx[0][0] - u[0] + v[0] : 0;
//    if (x[0] < s_current) {
//        if (s_pred < x[0]) {
//            tau *= (s_current - x[0]) / (s_current - s_pred);
//            return tau * (C_K/C_C) * u_xx[0][0] - u[0] + C_Tm;
//        } else {
//            return tau * (C_K/C_C) * u_xx[0][0] - u[0] + v[0];
//        }
//    } else {
//        return u[0] - C_To;
//    }
}

ANNS_CND(B1d_left) {return u[0] - C_Th;}
ANNS_CND(B1d_front) {return u[0] - C_Tm;}

ANNS_CND_G(A1d_a) {
//    anns_instance_t *ai = instance; return x[0] < s_current ? ai->step) * (C_K/C_C) * u_axx[0] - u_a : u_a;
    anns_instance_t *ai = instance;
    double tau = ai->step;
return x[0] < s_current ? tau * (C_K/C_C) * u_axx[0] - u_a : 0;
//    if (x[0] < s_current) {
//        if (s_pred < x[0]) {
//            tau *= (s_current - x[0]) / (s_current - s_pred);
//            return tau * (C_K/C_C) * u_axx[0] - u_a;
//        } else {
//            return tau * (C_K/C_C) * u_axx[0] - u_a;
//        }
//    } else {
//        return u_a;
//    }
}

ANNS_CND_G(B1d_a) {return u_a;}
//ANNS_CND_G(B1d_phase_a_1) {anns_instance_t *ai = instance; return (C_K*ai->step/C_lambda) * u_ax[0]; }
//ANNS_CND_G(B1d_phase_a_2) {anns_instance_t *ai = instance; return -u_a; }

//ANNS_CND_G(F_zero) {return 0; }



static void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
    ai_set_v(ai, 0, 0, 0);
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
}


static void ai_update_1d(anns_instance_t *ai, double *x, int i, double tau) {
    int p;
    anns_cache_t *ac = ai->cache;
    double point;
printf("\npred=%g current=%g\n", s_pred, s_current);
    // Сохраняем сеть i-1 шага в контейнер
    ai_put(ai, x, i-1);
s_pred = s_current;
    // Обновляем предысторию (v[0]) для конечно-разностной схемы (явно-неявная, Кранк-Никлсон)
    for (p = 0; p < ai->conds[0]->m; p++) {
//        point = ai->cache->points[0][p][0];
//        if (point < s_current) {
            ai->nets[0]->ann_val[0](x, 0, p, 0, 0, ai);
            ac->v[0][p][0] = ac->u[0];
//        } else {
//            ac->v[0][p][0] = 0;
//        }
    }



}

static void ai_preeval_1d(double *z, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    ai->cache->nets[0]->ann_val[1](z, 2, 0, 0, 0, ai);
    s_current = s_pred - (ai->step * C_K / C_lambda) * ai->cache->u_x[0][0];
    ai->cache->points[2][0][0] = s_current;
}

void stefan1phase(double Tm, int fdsteps) {
    int na = 30;

    int data[] = {
        PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 5, 1,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 1, 3, 1, fdsteps,
        na, na,  1, /* Уравнение переноса */
        1, 1,   0, /* Левая граница */
        1, 1,   0, /* Граница фазового перехода */
        NORMALIZED | RADIAL | GAUSSIAN, 2, 1, 0,  2, 0, 0,};

    double bounds[] = {
        GP_GRID_1D, na, C_xmax/na, C_xmax,
        GP_POINTS, 0,
        GP_POINTS, 0,
    };

    fbvp_t cond[] = {A1d, B1d_left, B1d_front};
    gbvp_t cond_a[] = {A1d_a, B1d_a, B1d_a};

    fbvp_t cond_i[] = {A1d_i, A1d_i, A1d_i};
    gbvp_t cond_i_a[] = {A1d_i_a, A1d_i_a, A1d_i_a};

//    double delta[] = {1e5, 1, 1};

    anns_methoddata_t *am = am_new3(ai_init, ai_newpoints_1d, ai_initann_1d, NULL, ai_update_1d,
                                    data, bounds,
                                    cond, cond_a,
                                    NULL,
//                                    delta,
                                    cond_i, cond_i_a, 0, Tm, fdsteps);

    am->randomseed = 0;
    am->seed = 123L;
    am->Nopt = 5;
    am->gtol = 1e-5;
    am->tol = 1e-5;
    am->n0 = 4;
    am->nofileout = 1;
    am->preeval = ai_preeval_1d;

    method_hybrid(am);

    am_delete(am);
}

#undef C_C
#undef C_K
#undef C_lambda
#undef C_To
#undef C_Tm
#undef C_Th
#undef C_xmax
