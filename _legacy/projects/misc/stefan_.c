#include "stefan.h"

static void ai_init(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
static void ai_newpoints_1d(anns_instance_t *, double *bounds);
static void ai_initann_1d(anns_instance_t *, double *);
static double ai_mse(anns_instance_t *, double *);

static double s = 0;
static double s_ = 0;

#define C_beta 1.0
#define C_a 1.0
#define C_Xm 3.0
#define C_f 1.0

ANNS_CND(A1d_i) {return u[0];}
ANNS_CND_G(A1d_i_a) {return u_a;}

ANNS_CND(A1d) {anns_instance_t *ai = instance; return x[0] < s ? .5*(ai->step)*C_a*u_xx[0][0] - u[0] + v[0] : 0;}
ANNS_CND(B1d_left) {return u_x[0][0] + C_f;}
ANNS_CND(B1d_right) {return u[0];}

ANNS_CND_G(A1d_a) {anns_instance_t *ai = instance; return x[0] < s ? .5*(ai->step)*C_a*u_axx[0] - u_a : 0;}
ANNS_CND_G(B1d_left_a) {return u_ax[0];}
ANNS_CND_G(B1d_a) {return u_a;}


static void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

static void ai_newpoints_1d(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
}


static void ai_initann_1d(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;
    double b = 1/sqrt(ai->nets[0]->nc);

    for (i = 0; i < ai->nets[0]->lenc; i+=3) {
        x[i] = uniform(0, 0);   // Веса
        x[i+1] = b;
        x[i+2] = uniform(-0.05, C_Xm + 0.05);
    }
}


static void ai_update_1d(anns_instance_t *ai, double *x, int i, double tau) {
    int p;
    anns_cache_t *ac = ai->cache;
    // Сохраняем сеть i-1 шага в контейнер
    ai_put(ai, x, i-1);
    s_ = s;

    // Обновляем предысторию (v[0]) для конечно-разностной схемы (явно-неявная, Кранк-Никлсон)
    for (p = 0; p < ai->conds[0]->m; p++) {
        ai->nets[0]->ann_val[2](x, 0, p, 0, 0, ai);
        ac->v[0][p][0] = 0.5*tau*C_a*ac->u_xx[0][0] + ac->u[0];
    }
}

static void ai_preeval_1d(double *z, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
//net->ann_val[eval % MIXED_DERIVATIVES](x, c, p, t, eval & MIXED_DERIVATIVES, ai);
    ai->cache->nets[0]->ann_val[1](z, 2, 0, 0, 0, ai);
//    printf("->s=%f s_=%f p=%f ux=%f\n ", s, s_, ai->cache->points[2][0][0], ai->cache->u_x[0][0]);
    s = s_ - ai->step * ai->cache->u_x[0][0] / C_beta;
//    printf("<-s=%f s_=%f p=%f ux=%f\n ", s, s_, ai->cache->points[2][0][0], ai->cache->u_x[0][0]);
    ai->cache->points[2][0][0] = s;
}


void stefan_1d(double Tm, int fdsteps) {
    int na = 30;

    int data[] = {
        PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 4, 1,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 1, 3, 1, fdsteps,
        na, na,  1, /* Уравнение переноса */
        1, 1,   0, /* Левая граница */
        1, 1,   0, /* Граница фазового перехода */
        NORMALIZED | RADIAL | GAUSSIAN, 2, 1, 0,  2, 1, 0};

    double bounds[] = {
        GP_GRID_1D, na, 0, C_Xm,
        GP_POINTS, 0,
        GP_POINTS, 0,
    };

    fbvp_t cond[] = {A1d, B1d_left, B1d_right };
    gbvp_t cond_a[] = {A1d_a, B1d_left_a, B1d_a};

    fbvp_t cond_i[] = {A1d_i, A1d_i, A1d_i};
    gbvp_t cond_i_a[] = {A1d_i_a, A1d_i_a, A1d_i_a};

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
