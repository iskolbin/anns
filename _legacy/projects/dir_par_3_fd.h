#ifndef ANNS_TASK
#define ANNS_TASK

#include "../bvp/bvp_ann.h"

#define FD_STEPS_EVAL 21
#define FD_STEPS 21

int data[] = {
    PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 3, 1,

    FUNCTIONAL_SQUARE, 1, 3, 1, FD_STEPS_EVAL,

    20, 1,  // Уравнение переноса (v[0] - предыстория)
    1, 0,   // Изолированная сторона
    1, 1,   // Нагреваемая сторона (v[0] - мощность нагрева)


    NORMALIZED | RADIAL | GAUSSIAN, 8, 1, 0,   2, 1, 0,
};

#define TAU (1.0/(FD_STEPS-1))

double bounds[] = {
    GP_GRID_1D, 20, 0.0, 1.0,   // Уравнение переноса
    GP_UNIFORM_TILE, 0.0, 0.0,     // Изолированная сторона
    GP_UNIFORM_TILE, 1.0, 1.0, // Нагреваемая сторона
};

static double ksi(double t) {return (t < 0.5) ? 2*t : 2 - 2*t;}

ANNS_CND(A) {return 0.5*TAU*u_xx[0][0] - u[0] + v[0];}
ANNS_CND(B) {return u_x[0][0];}
ANNS_CND(C) {return u[0] - v[0];}

//ANNS_CND_G(A_a) {return TAU*u_axx[0] - F_u(x[0])*u_a;}
ANNS_CND_G(A_a) {return 0.5*TAU*u_axx[0] - u_a;}
ANNS_CND_G(B_a) {return u_ax[0];}
ANNS_CND_G(C_a) {return u_a;}

fbvp2_t cond[] = {A, B, C};
gbvp2_t cond_a[] = {A_a, B_a, C_a};

double delta[] = {1.0, 1.0, 1.0};

void init(bvp_instance_t *bi, double *x) {
    int i, p;

    bi_load_conditions(bi, cond, cond_a, delta);
    bi_load_netfuncs(bi);
    bi_generate_points(bi, bounds);

    for (i = bi->I[0]; i < bi->I[1]; ) {
        x[i++] = uniform(0.0, 0.0);          // Веса (в начальный момент t = 0)
        x[i++] = uniform(0.05, 1.05);        // Ширины
        x[i++] = uniform(-0.25, 1.25);  // Центры
    }

//    for (p = 0; p < bi->np[0]; p++) {
//        bi->v[0][p][0] = F_u(bi->points[0][p][0]);
//    }
}

void update(bvp_instance_t *bi, double *x, int i) {
    int p;
    bvp_cache_t *bc = bi->cache;

    // Сохраняем сеть i-1 шага в контейнер
    bc_put(bc, x, i-1);

    // Обновляем предысторию (v[1]) для конечно-разностной схемы (явно-неявная, Кранк-Никлсон)
    for (p = 0; p < bi->np[0]; p++) {
        bi->ann_val[0][2](bi->points[0][p], x, 0, 0, bi);
        bi->v[0][p][0] = 0.5*TAU*bc->u_xx[0][0] + bc->u[0];
//        bi->ann_val[0][0](bi->points[0][p], x, 0, 0, bi);
//        bi->v[0][p][1] = F_u(bi->points[0][p][0])*bc->u[0];
    }

    // Обновляем поток
    bi->v[2][0][0] = ksi(i * TAU);
}

double eval_mse(double *x, bvp_instance_t *bi) {
    return 0.0;
}

#endif // ANNS_TASK
