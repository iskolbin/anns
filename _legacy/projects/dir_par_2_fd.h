#ifndef ANNS_TASK
#define ANNS_TASK

#include "../bvp/bvp_ann.h"

#define FD_STEPS_EVAL 161
#define FD_STEPS 161

int data[] = {
    PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 2, 1,

    FUNCTIONAL_SQUARE, 1, 3, 1, FD_STEPS_EVAL,

    20, 1,
    1, 1,
    1, 1,

    NORMALIZED | RADIAL | GAUSSIAN, 4, 1, 0,   2, 0, 1,
};

#define C_T 7.0
#define TAU (C_T/(FD_STEPS-1))

double bounds[] = {
    GP_GRID_1D, 20, 0.01, M_PI - 0.01,
    GP_UNIFORM_TILE, 0.0, 0.0,
    GP_UNIFORM_TILE, M_PI, M_PI,
};

ANNS_CND(A) {return 0.5*TAU*u_xx[0][0] - u[0] + v[0];}
ANNS_CND(B) {return u[0] - v[0];}
ANNS_CND(C) {return u_x[0][0] - v[0];}

ANNS_CND_G(A_a) {return 0.5*TAU*u_axx[0] - u_a;}
ANNS_CND_G(B_a) {return u_a;}
ANNS_CND_G(C_a) {return u_ax[0];}

fbvp2_t cond[] = {A, B, C};
gbvp2_t cond_a[] = {A_a, B_a, C_a};

double delta[] = {1.0, 1.0, 1.0};

void init_bvp_instance (bvp_instance_t *bi) {
    bi_load_conditions(bi, cond, cond_a, delta);
    bi_load_netfuncs(bi);
}

void reload_bvp_instance(bvp_instance_t *bi) {
    bi_generate_points(bi, bounds);
//    bi_eval_v(bi, 1, 0, F_u);
}

void fd_step(bvp_instance_t *bi, double *x, int i) {
    int p;
    double pointX[1], x_;

    bvp_cache_t *bc = bi->cache;

    // Сохраняем сеть i-1 шага в контейнер
    bc_put(bc, x, i-1);

    bi->v[1][0][0] = sin(TAU*i);
    bi->v[2][0][0] = -sin(TAU*i);

    // Обновляем v для конечно-разностной схемы (Кранк-Никлсон)
    for (p = 0; p < bi->np[0]; p++) {
        x_ = bi->points[0][p][0];
        bi->ann_val[0][2](bi->points[0][p], x, 0, 2, bi);
        bi->v[0][p][0] = 0.5*TAU*(cos(x_)*(cos(TAU*i) + sin(TAU*i) + cos(TAU*(i-1)) + sin(TAU*(i-1))) + bc->u_xx[0][0]) + bc->u[0];
//        bi->ann_val[0][0](bi->points[1][p], x, 0, 0, bi);
//        bi->v[1][p][1] = -bi->v[1][p][0]*bc->u[0];
    }
}

void init_ann(double *x, bvp_instance_t *bi) {
    int i, j;

    for (i = bi->I[0]; i < bi->I[1]; ) {
        x[i++] = uniform(0.0, 0.0);        // Веса (в начальный момент T = 0)
        x[i++] = uniform(0.05, 1.05);        // Ширины
        x[i++] = uniform(-0.05, 1.05);     // Центры
    }
}

double eval_mse(double *x, bvp_instance_t *bi) {
    return 0.0;
}

#endif // ANNS_TASK
