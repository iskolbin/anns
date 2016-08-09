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

#define C_Tmax 1.0


ANNS_CND(F0) {return 0;}
ANNS_CND_G(F0_a) {return 0;}
ANNS_CND_G(D_a) {return u_a;}


// 1 - Уравнение переноса
//ANNS_CND(A) {return u_x[0][0] - u_xx[0][1]/(u[1]*u[1]) - x[1]*u_x[1][0]*u_x[0][1]/u[1];}
//ANNS_CND_G(A_a_z) {return u_ax[0] - u_axx[1]/(u[1]*u[1]) - x[1]*u_x[1][0]*u_ax[1]/u[1];}
//ANNS_CND_G(A_a_s) {return x[1]*u_x[0][1]*u_x[1][0]*u_ax[0]/u[1] - 2*u_xx[0][1]*u_a/(u[1]*u[1]) - x[1]*u_x[0][1]*u_ax[0];}
//ANNS_CND(A) {return u[1]*u[1]*u_x[0][0] - u_xx[0][1] - x[1]*u[1]*u_x[1][0]*u_x[0][1];}
//ANNS_CND_G(A_a_z) {return u[1]*u[1]*u_ax[0] - u_axx[1] - x[1]*u[1]*u_x[1][0]*u_ax[1];}
//ANNS_CND_G(A_a_s) {return 2*u[1]*u_a*u_x[0][0] - x[1]*u_x[0][1]*(u_a*u_x[1][0] + u[1]*u_ax[0]);}

ANNS_CND(A) {return u_x[0][0] - u_xx[0][1]/(u[1]*u[1]) - (x[1]/u[1])*u_x[1][0]*u_x[0][1];}
ANNS_CND_G(A_a_z) {return u_ax[0] - u_axx[1]/(u[1]*u[1]) - (x[1]/u[1])*u_x[1][0]*u_ax[1];}
ANNS_CND_G(A_a_s) {return 2*u_a*u_xx[0][1]/(u[1]*u[1]*u[1]) + (x[1]/(u[1]*u[1]))*u_a*u_x[1][0]*u_x[0][1] - (x[1]/u[1])*u_ax[0]*u_x[0][1];}


// 2 - Начальное условие s(0) = s0
ANNS_CND(PI) {return u[1];}

// 3 - Начальное условие u(0, x) = x - 1
//ANNS_CND(AI) {return u[0] - x[1] + 1;}
ANNS_CND(AI) {return u[0];}
// 4 - Граничное условие u(t, 0) = -1
ANNS_CND(B) {return u[0] + 1;}

// 5 - Граничное условие u(t, 1) = 0 (на самом деле - фронт фазового перехода)
ANNS_CND(P1) {return u[0];}

// 6 - Граничное условие u(t, 1) (на самом деле - фронт фазового перехода)
ANNS_CND(P) {return u_x[1][0] - u_x[0][1]/u[1];}
ANNS_CND_G(P_a_z) {return -u_ax[1]/u[1];}
ANNS_CND_G(P_a_s) {return u_ax[0] + u_a*u_x[0][1]/(u[1]*u[1]);}
//ANNS_CND(P) {return u[1]*u_x[1][0] - u_x[0][1];}
//ANNS_CND_G(P_a_z) {return -u_ax[1];}
//ANNS_CND_G(P_a_s) {return u_a*u_x[1][0] + u[1]*u_ax[0];}

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
        x[i++] = uniform(-1, 1);   // Веса
        x[i++] = 1;
        x[i++] = uniform(-0.05, C_Tmax + 0.05);
        x[i++] = uniform(-0.05, 1 + 0.05);
    }

    // Сеть для фронта фазового перехода
    for (; i < ai->nets[1]->lenc + ai->nets[0]->lenc;) {
        x[i++] = uniform(1, 1);   // Веса
        x[i++] = 1;
        x[i++] = uniform(-0.05, C_Tmax + 0.05);
    }

    an_setwidth_nearest(ai->nets[0], 1.5);
    an_setwidth_nearest(ai->nets[1], 1.5);
}

void stefan1phase(void) {
    int nx = 15, nt = 15;


    // Входной вектор (t,x)
    // Сети 1 - фаза u1(t,x)
    //   2 - фронт фазового перехода s(t)
    int data[] = {
        PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 5, 1,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 2, 6, 2, 0,
        nx*nt , nx*nt,  0,  /* 1 - Уравнение переноса */
        1, 1,   0,          /* 2 - Начальное условие s(0) = 0 */
        nx, nx, 0,          /* 3 - Начальное условие u(0, x) = x - 1 */
        nt, nt, 0,          /* 4 - Граничное условие u(t, 0) = -1 */
        nt, nt, 0,          /* 5 - Граничное условие u(t, 1) = 0 (на самом деле - фронт фазового перехода) */
        nt, nt, 0,          /* 6 - Граничное условие u(t, 1) (на самом деле - фронт фазового перехода) */

        NORMALIZED | RADIAL | GAUSSIAN, 16, 2, 0,  2, -1, 0, 0, 0, 1,
        NORMALIZED | RADIAL | GAUSSIAN, 5, 1, 0,   1, 0, -1, -1, -1, 1};

    double bounds[] = {
//        GP_UNIFORM_CUBOID, 0.05, C_Tmax,  0.05, 1.0,
        GP_GRID_2D,   nt, C_Tmax/nt, C_Tmax,    nx, 1.0/nx, 1.0, // 1
        GP_GRID_2D,   1, 0, 0,                  1, 0, 0,         // 2
        GP_GRID_2D,   1, 0, 0,                  nx, 1.0/nx, 1.0, // 3
        GP_GRID_2D,   nt, C_Tmax/nt, C_Tmax,    1, 0, 0,         // 4
        GP_GRID_2D,   nt, C_Tmax/nt, C_Tmax,    1, 0, 0,         // 5
        GP_GRID_2D,   nt, C_Tmax/nt, C_Tmax,    1, 0, 0,         // 6
    };

    fbvp_t cond[] = {A, PI, AI, B, P1, P};
    gbvp_t cond_a[] = {
        A_a_z, A_a_s,
        NULL, D_a,
        D_a, NULL,
        D_a, NULL,
        D_a, NULL,
        P_a_z, P_a_s,

    };

//    double *delta = NULL;
    double delta[] = {1,100,1,1,1,1};

    anns_methoddata_number_t mnumbers[] = {
        AM_SEED,        1231111,
        AM_GTOL,        1e-4,
        AM_TOL,         1e-4,
        AM_NOFILEOUT,   0,

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

        AM_NULL,        NULL,
    };


    anns_methoddata_t *am = am_new_n(mpointers, mnumbers);

    method_restarts(am);

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
