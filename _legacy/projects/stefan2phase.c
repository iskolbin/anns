#include "stefan.h"

static void ai_init(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
static void ai_newpoints(anns_instance_t *, double *bounds);
static void ai_initann(anns_instance_t *, double *);
static double ai_mse(anns_instance_t *, double *);

#define C_a1 1.2
#define C_a2 1.0
#define C_up 0.0
#define C_u0 -1.0
#define C_k1 1.2
#define C_k2 1.0
#define C_q  1.0
#define C_T  3.0

ANNS_CND(F0) {return 0;}
ANNS_CND_G(F0_a) {return 0;}
ANNS_CND_G(D_a) {return u_a;}

// 1 - Уравнение переноса
ANNS_CND(A) {return (u[2] > x[0]) ?
    (u_x[0][0] - C_a1*C_a1 * u_xx[0][1]) :
    (u_x[1][0] - C_a2*C_a2 * u_xx[1][1]);}

ANNS_CND_G(A_a_1) {return (u[2] > x[0]) ? u_ax[0] - C_a1*C_a1 * u_axx[1] : 0;}
ANNS_CND_G(A_a_2) {return (u[2] < x[0]) ? u_ax[0] - C_a2*C_a2 * u_axx[1] : 0;}

// 2 - Граничное условие x = 0 (phi = t - 1)
ANNS_CND(B0) {return u[0] - x[1] + 1;}

// 3 - s(t) для первой фазы
ANNS_CND(P1) {return u[0] - C_up;}

// 4 - Начальное условие
ANNS_CND(I) {return u[1] - C_u0;}

// 5 - Граничное условие x = 1 (psi = -1)
ANNS_CND(B1) {return u[1] + 1;}

// 6 - s(t) для второй фазы
ANNS_CND(P2) {return u[1] - C_up;}

// 7 - s(t) для двух фаз
ANNS_CND(P12) {return u[0] - u[1];}
ANNS_CND_G(P12_a_1) {return u_a;}
ANNS_CND_G(P12_a_2) {return -u_a;}

// 8 - s(t)
ANNS_CND(P) {return C_k1 * u_x[0][1] - C_k2 * u_x[1][1] - C_q*u_x[2][0];}

ANNS_CND_G(P_a_1) {return C_k1 * u_ax[1];}
ANNS_CND_G(P_a_2) {return -C_k2 * u_ax[1];}
ANNS_CND_G(P_a_s) {return -C_q * u_ax[0];}

ANNS_CND(PI) {return u[2];}

static void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

static void ai_newpoints(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
}

static void ai_initann(anns_instance_t *ai, double *x) {
    int i = 0, j, t, n, s;
    double b;

    // Сеть для первой фазы
    for (; i < ai->nets[0]->lenc;) {
        x[i++] = uniform(0, 1);   // Веса
        x[i++] = 1;
        x[i++] = uniform(-0.05, C_T + 0.05);
        x[i++] = uniform(-0.05, 1 + 0.05);
    }

    // Сеть для второй фазы
    for (; i < ai->nets[1]->lenc + ai->nets[0]->lenc;) {
        x[i++] = uniform(0, 1);   // Веса
        x[i++] = 1;
        x[i++] = uniform(-0.05, C_T + 0.05);
        x[i++] = uniform(-0.05, 1 + 0.05);
    }

    // Сеть для фронта фазового перехода
    for (; i < ai->nets[2]->lenc + ai->nets[1]->lenc + ai->nets[0]->lenc;) {
        x[i++] = uniform(0, 0);   // Веса
        x[i++] = 1;
        x[i++] = uniform(-0.05, C_T + 0.05);
    }

    an_setwidth_nearest(ai->nets[0], 1.5);
    an_setwidth_nearest(ai->nets[1], 1.5);
    an_setwidth_nearest(ai->nets[2], 1.5);
}

static void ai_updatepoints(double *z, void *anns_instance) {
    anns_instance_t *ai = anns_instance;
    anns_cache_t *ac = ai->cache;
    int p;

    // Обновляем координаты точек фронта фазового перехода
    for (p = 0; p < ac->conds[1]->m; p++) {
        ac->nets[2]->ann_val[0](z, 1, p, 2, 0, ai);
        ac->points[2][p][1] = ac->u[2];
        ac->points[5][p][1] = ac->u[2];
        ac->points[6][p][1] = ac->u[2];
        ac->points[7][p][1] = ac->u[2];
    }
}

static double ai_mse(anns_instance_t *ai, double *x) {
    return -1;
}

void stefan2phase1(double Tm, double Xm) {
    int nx = 20, nt = 20;

    // Входной вектор (t,x)
    // Сети 1 - первая фаза u1(t,x)
    //   2 - вторая фаза u2(t,x)
    //   3 - фронт фазового перехода s(t)
    int data[] = {
        PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 6, 1,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 2, 9, 3, 0,

        nx*nt , nx*nt,  0, /* 1 - Уравнение переноса */
        nt, nt, 0, /* 2 - x = 0 */
        nt, nt, 0,
        nx, nx, 0, /* 3 - t = 0 */
        nt, nt, 0, /* 4 - s(t) для u1 */
        nt, nt, 0, /* 5 - s(t) для u2 */
        nt, nt, 0, /* 6 - s(t) для s */
        nt, nt, 0,
        1,  1,  0,

        NORMALIZED | RADIAL | GAUSSIAN, 8, 2, 0,  2, 0, 0, -1, -1, -1, 0, 1, -1,
        NORMALIZED | RADIAL | GAUSSIAN, 8, 2, 0,  2, -1, -1, 0, 0, 0, 0, 1, -1,
        NORMALIZED | RADIAL | GAUSSIAN, 4, 1, 0,  0, -1, -1, -1, -1, -1, -1, 1, 0 };

    double bounds[] = {
        GP_GRID_2D,   nt, 0, C_T,  nx, 0, 1, // 1
        GP_GRID_2D,   nt, 0, C_T,  1, 0, 0,  // 2
        GP_GRID_2D,   nt, 0, C_T,  1, 0, 0,  // 3
        GP_GRID_2D,   1,  0, 0,    nx, 0, 1, // 4
        GP_GRID_2D,   nt, 0, C_T,  1, 1, 1,  // 5
        GP_GRID_2D,   nt, 0, C_T,  1, 0, 0,  // 6
        GP_GRID_2D,   nt, 0, C_T,  1, 0, 0,  // 7
        GP_GRID_2D,   nt, 0, C_T,  1, 0, 0,  // 8
        GP_POINTS,    0, 0,
    };

    fbvp_t cond[] = {A, B0, P1, I, B1, P2, P12, P, PI};
    gbvp_t cond_a[] = {
        A_a_1, A_a_2, F0_a,
        D_a, NULL, NULL,
        D_a, NULL, NULL,
        NULL, D_a, NULL,
        NULL, D_a, NULL,
        NULL, D_a, NULL,
        P12_a_1, P12_a_2, NULL,
        P_a_1, P_a_2, P_a_s,
        NULL, NULL, D_a
    };

    anns_methoddata_number_t mnumbers[] = {
        AM_SEED,        12311,
        AM_GTOL,        1e-3,
        AM_TOL,         1e-3,
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
        AM_UPDATEPOINTS, ai_updatepoints,

        AM_NULL,        NULL,
    };


    anns_methoddata_t *am = am_new_n(mpointers, mnumbers);

//    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints, ai_initann, ai_mse, NULL, data, bounds, cond, cond_a, delta);

//    am->randomseed = 0;
//    am->seed = 12311L;
//    am->Nopt = 5;
//    am->gtol = 1e-5;
//    am->tol = 1e-5;
//    am->n0 = 4;
//    am->nofileout = 0;
//    am->preeval = ai_preeval;
//
//    method_restarts(am);
    method_updatepoints(am);

    am_delete(am);
}


#undef C_T0
#undef C_Tc
#undef C_Tm
#undef C_C1
#undef C_C2
#undef C_a1
#undef C_a2
#undef C_u0
#undef C_k1
#undef C_k2
#undef C_s0
