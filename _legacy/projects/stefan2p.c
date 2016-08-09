#include "stefan.h"

static void ai_init(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
static void ai_newpoints(anns_instance_t *, double *bounds);
static void ai_initann(anns_instance_t *, double *);
static double ai_mse(anns_instance_t *, double *);

#define C_C_S 1.762
#define C_C_L 4.22
#define C_k_S 2.22
#define C_k_L 0.556
#define C_rho 1000.
#define C_a_S (C_k_S/(C_C_S*C_rho))
#define C_a_L (C_k_L/(C_C_L*C_rho))
#define C_h_SL 338.
#define C_Ti -5.
#define C_T0 20.0
#define C_Tm 0.0

#define C_T 40.0

// L - жидкая фаза
// S - твердая фаза

// Первая сеть - жидкая фаза u0(t, x)
// Вторая сеть - твердая фаза u1(t, x)
// Третья сеть - фронт u2(t)

// 1 - уравнение переноса для жидкой фазы
ANNS_CND(A_L) {return u_x[0][0] - C_a_L * u_xx[0][1];}

// 2 - уравнение переноса для твердой фазы
ANNS_CND(A_S) {return u_x[1][0] - C_a_S * u_xx[1][1];}

// 3 - граничное условие x = 0
ANNS_CND(B) {return u[0] - C_T0;}

// 4 - начальное условие t = 0
ANNS_CND(I) {return u[1] - C_Ti;}

// 5 - фронт фазового перехода
ANNS_CND(P) {return C_k_S * u_x[1][1] - C_k_L * u_x[0][1] - C_rho * C_h_SL * u_x[2][0];}

// 6 - согласование на подвижной границе
ANNS_CND(P1) {return u[0] - C_Tm;}

// 7 - согласование на подвижной границе
ANNS_CND(P2) {return u[1] - C_Tm;}

// 8 - изначальное положение границы
ANNS_CND(PI) {return u[2];}



ANNS_CND_G(A_L_a) {return u_ax[0] - C_a_L * u_axx[1];}
ANNS_CND_G(A_S_a) {return u_ax[0] - C_a_S * u_axx[1];}
ANNS_CND_G(B_a) {return u_a;}
ANNS_CND_G(I_a) {return u_a;}

ANNS_CND_G(P_a_0) {return -C_k_L * u_ax[1];}
ANNS_CND_G(P_a_1) {return C_k_S * u_ax[1];}
ANNS_CND_G(P_a_2) {return C_k_S * u_xx[1][1] * u_a - C_k_L * u_xx[0][1] * u_a - C_rho * C_h_SL * u_ax[0];}

ANNS_CND_G(P1_a_0) {return u_a;}
ANNS_CND_G(P1_a_2) {return u_x[0][1] * u_a;}

ANNS_CND_G(P2_a_1) {return u_a;}
ANNS_CND_G(P2_a_2) {return u_x[1][1] * u_a;}

ANNS_CND_G(PI_a) {return u_a;}

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


void stefan2phase(void) {
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

    fbvp_t cond[] = {A_L, A_S, B, I, P, P1, P2, PI};
    gbvp_t cond_a[] = {
        A_L_a, NULL, NULL,
        NULL, A_S_a, NULL,
        B_a, NULL, NULL,
        NULL, I_a, NULL,
        P_a_0, P_a_1, P_a_2,
        P1_a_0, NULL, P1_a_2,
        NULL, P2_a_1, P2_a_2,
        NULL, NULL, PI_a,
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

#undef C_C_S
#undef C_C_L
#undef C_k_S
#undef C_k_L
#undef C_rho
#undef C_a_S
#undef C_a_L
#undef C_h_SL
#undef C_Ti
#undef C_T0
#undef C_Tm
