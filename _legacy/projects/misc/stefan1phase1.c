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
#define C_Xmax 1.0


ANNS_CND(F0) {return 0;}
ANNS_CND_G(F0_a) {return 0;}
ANNS_CND_G(D_a) {return u_a;}


// 1 - Уравнение переноса
ANNS_CND(A) {return (u[1] > x[1]) ? (u_x[0][0] - C_a * u_xx[0][1]) : u[0];}
ANNS_CND_G(A_a_1) {return (u[1] > x[1]) ? (u_ax[0] - C_a * u_axx[1]) : u_a;}

// 2 - Начальное условие
ANNS_CND(PI) {return u[1];}

// 3 - Граничное условие
ANNS_CND(B) {return u[0] - C_alpha;}
ANNS_CND_G(B_a_1) {return u_a;}

// 4 - s(t)
ANNS_CND(P1) {return u[0];}

// 5 - s(t)
ANNS_CND(P) {return (C_K/(C_lambda*C_sigma)) * u_x[0][1] - u_x[1][0];}
ANNS_CND_G(P_a_1) {return (C_K/(C_lambda*C_sigma)) * u_ax[1];}
ANNS_CND_G(P_a_s) {return -u_ax[0];}


static void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

static void ai_newpoints(anns_instance_t *ai, double *bounds) {
    anns_cache_t *ac = ai->cache;
    int p, m = ac->conds[3]->m;

    ai_genpoints(ai, bounds);

    for (p = 0; p < m; p++) {
        ac->points[3][p][1] = ((p+1.)/m)*(C_Xmax/C_Tmax);
        ac->points[4][p][1] = ((p+1.)/m)*(C_Xmax/C_Tmax);
        printf("%g ", ac->points[4][p][1]);
    }
}

static void ai_initann(anns_instance_t *ai, double *x) {
    int i = 0, j, t, n, s;
    double b;

    // Сеть для первой фазы
    for (; i < ai->nets[0]->lenc;) {
        x[i++] = uniform(0, 2);   // Веса
        x[i++] = 1;
        x[i++] = uniform(-0.05, C_Tmax + 0.05);
        x[i++] = uniform(-0.05, C_Xmax + 0.05);
    }

    // Сеть для фронта фазового перехода
    for (; i < ai->nets[1]->lenc + ai->nets[0]->lenc;) {
        x[i++] = uniform(0, 0);   // Веса
        x[i++] = 1;
        x[i++] = uniform(-0.05, C_Tmax + 0.05);
    }

    an_setwidth_nearest(ai->nets[0], 1.5);
    an_setwidth_nearest(ai->nets[1], 1.5);
}

static void ai_updatepoints(double *z, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int p;

    // Обновляем координаты точек фронта фазового перехода
    for (p = 0; p < ac->conds[3]->m; p++) {
        ac->nets[1]->ann_val[0](z, 3, p, 1, 0, ai);
        ac->points[3][p][1] = ac->u[1];
        ac->points[4][p][1] = ac->u[1];
        printf("%g ", ac->u[1]);
    }
    printf("\n");
}

void stefan1phase(void) {
    int nx = 20, nt = 20;


    // Входной вектор (t,x)
    // Сети 1 - фаза u1(t,x)
    //   2 - фронт фазового перехода s(t)
    int data[] = {
        PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 5, 1,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 2, 5, 2, 0,
        nx*nt , nx*nt,  0, /* 1 - Уравнение переноса */
        1, 1,   0, /* 2 - Начальное условие s(t) */
        nt, nt, 0, /* 3 - Условие на границе x = 0 */
        nt, nt, 0, /* 4 - s(t) для u */
        nt, nt, 0, /* 5 - s(t) для s */

        NORMALIZED | RADIAL | GAUSSIAN, 16, 2, 0,  2, -1, 0, 0,  1,
        NORMALIZED | RADIAL | GAUSSIAN, 3, 1, 0,  0, 0, -1, -1, 1,};

    double bounds[] = {
        GP_GRID_2D,   nt, C_Tmax/nt, C_Tmax,  nx, C_Xmax/nx, C_Xmax, // 1
        GP_POINTS, 0.0, 0.0,                         // 2
        GP_GRID_2D,   nt, C_Tmax/nt, C_Tmax,  1, 0, 0,       // 3
        GP_GRID_2D,   nt, C_Tmax/nt, C_Tmax,  1, 0, 0,       // 4
        GP_GRID_2D,   nt, C_Tmax/nt, C_Tmax,  1, 0, 0,       // 5
    };

    fbvp_t cond[] = {A, PI, B, P1, P};
    gbvp_t cond_a[] = {
        A_a_1, F0_a,
        NULL, D_a,
        B_a_1, NULL,
        D_a, NULL,
        P_a_1, P_a_s,
    };

    double *delta = NULL;

    anns_methoddata_number_t mnumbers[] = {
        AM_SEED,        1231,
        AM_GTOL,        1e-4,
        AM_TOL,         1e-4,
        AM_NOFILEOUT,   0,
        AM_MAXITERS,    20,

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
        AM_UPDATEPOINTS, ai_updatepoints,

        AM_NULL,        NULL,
    };


    anns_methoddata_t *am = am_new_n(mpointers, mnumbers);

    method_updatepoints(am);

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
