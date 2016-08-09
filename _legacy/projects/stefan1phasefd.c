#include "stefan.h"

static void ai_init(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
static void ai_newpoints(anns_instance_t *, double *bounds);
static void ai_initann(anns_instance_t *, double *);
static double ai_mse(anns_instance_t *, double *);

#define C_K 1.0
#define C_c 1.0
#define C_sigma 1.0
#define C_alpha -1.0
#define C_lambda 1.0
#define C_a (C_K/(C_c*C_sigma))

#define C_Tmax 10.0
#define C_Xmax 10.0

static double h_current = 0;
static double h_pred = 0;


ANNS_CND_G(D_a) {return u_a;}

// 1 - Уравнение переноса
ANNS_CND(A) {double tau = ((anns_instance_t*)instance)->step;
//printf("%g ", x[0]);
return (h_current > x[0]) ? u[0] - (tau*C_a ) * u_xx[0][0] - v[0] : u[0];}
ANNS_CND_G(A_a) {double tau = ((anns_instance_t*)instance)->step; return (h_current > x[0]) ? u_a - (tau*C_a) * u_axx[0] : u_a;}

// 2 - Граничное условие
ANNS_CND(B) {return u[0] - C_alpha;}

// 3 - Условие на фронте фазового перехода
ANNS_CND(P) {return u[0];}


static void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

static void ai_newpoints(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
}

static void ai_update(anns_instance_t *ai, double *x, int i, double tau) {
    int p;
    double t = i * tau;
    anns_cache_t *ac = ai->cache;

    ai_put(ai, x, i-1);

    // Обновляем предысторию (v[0]) для конечно-разностной схемы (явно-неявная, Кранк-Никлсон)
    for (p = 0; p < ai->conds[0]->m; p++) {
        ai->nets[0]->ann_val[2](x, 0, p, 0, 0, ai);
        ac->v[0][p][0] = .5*tau*C_a*ac->u_xx[0][0]+ ac->u[0];
    }

    h_pred = h_current;

}

static void ai_initann(anns_instance_t *ai, double *x) {
    int i = 0, j, t, n, s;
    double b;

    // Сеть для первой фазы
    for (; i < ai->nets[0]->lenc;) {
        x[i++] = uniform(-1, -1);   // Веса
        x[i++] = 1;
        x[i++] = uniform(-0.05, C_Xmax + 0.05);
    }

    an_setwidth_nearest(ai->nets[0], 1.5);
}

static void ai_preeval(double *z, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    double tau = ai->step;
    int p;

    // Обновляем координаты точек фронта фазового перехода
    h_current = h_pred + (C_lambda*C_sigma*C_K) * ac->u_x[0][0];
    ac->points[2][0][0] = h_current;
    printf("[h_current = %g  h_pred = %g]\n", h_current, h_pred);
}

void stefan1phasefd(double Tmax, int fdsteps) {
    int nx = 5;

    // Входной вектор (t,x)
    // Сети 1 - фаза u1(t,x)
    //   2 - фронт фазового перехода s(t)
    int data[] = {
        PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 5, 2,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 1, 3, 1, fdsteps,
        nx , nx,  1, /* 1 - Уравнение переноса */
        1, 1,   0, /* 2 - Условие на границе x = 0 */
        1, 1,   0, /* 3 - Условие на фронте фазового перехода h */

        NORMALIZED | RADIAL | GAUSSIAN, 4, 1, 0,  2, 0, 0,};

    double bounds[] = {
        GP_GRID_1D,   nx, C_Xmax/nx, C_Xmax, // 1
        GP_POINTS, 0.0,                   // 2
        GP_POINTS, 0.0,                 // 3
    };

    fbvp_t cond[] = {A, B, P};
    gbvp_t cond_a[] = {A_a, D_a, D_a};

    double *delta = NULL;

    anns_methoddata_number_t mnumbers[] = {
        AM_SEED,        1231,
        AM_GTOL,        1e-4,
        AM_TOL,         1e-4,
        AM_NOFILEOUT,   1,
        AM_MAXITERS,    20,

        AM_TBEGIN,      0,
        AM_TEND,        Tmax,
        AM_TSTEPS,      fdsteps,

        AM_NULL,        0,
    };

    anns_methoddata_pointer_t mpointers[] = {
        AM_INIT,        ai_init,
        AM_NEWPOINTS,   ai_newpoints,
        AM_INITANN,     ai_initann,
        AM_DATA,        data,
        AM_BOUNDS,      bounds,
        AM_COND,        cond,
        AM_COND_A,      cond_a,
        AM_DELTA,       delta,
        AM_PREEVAL,     ai_preeval,
        AM_UPDATE,      ai_update,

        AM_NULL,        NULL,
    };


    anns_methoddata_t *am = am_new_n(mpointers, mnumbers);

    method_hybrid(am);

    am_delete(am);
}

#undef C_K
#undef C_c
#undef C_sigma
#undef C_alpha
#undef C_lambda
#undef C_a
#undef C_Tmax
#undef C_Xmax

