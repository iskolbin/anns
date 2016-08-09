#include "curvelinear.h"

static void ai_init(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
static void ai_newpoints(anns_instance_t *, double *bounds);
static void ai_initann(anns_instance_t *, double *);
static double ai_mse(anns_instance_t *, double *);

#define C_r 0.5
#define C_R 1.0

void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_newpoints(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;
    double b = 1/sqrt(ai->nets[0]->nc), r, angle;

    for (i = 0; i < ai->nets[0]->lenc; i+=4) {

        x[i] = uniform(-2, 2.);   // Веса
//            x[i++] = uniform(0.5, 2.);   // Ширины
        x[i+1] = b;

        r = uniform(C_r, C_R);
        angle = uniform(0, 0.5*M_PI);

        x[i+2] = r * cos(angle);
        x[i+3] = r * sin(angle);
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    double s, tmp, point[2], r, angle, x2, y2;
    int i, p, M = 5000;

    s = 0;
    for (p = 0; p < M; p++) {
        r = uniform(0.5, 1);
        angle = uniform(0, 0.5*M_PI);

        point[0] = r * cos(angle);
        point[1] = r * sin(angle);
        x2 = point[0] * point[0];
        y2 = point[1] * point[1];

//        tmp = anns_eval(point, x, ai) - 3.5*log(x2 + y2) / log(2) - 4;
        tmp = anns_eval(point, x, ai) - 0.5*log(x2 + y2) / log(2) - 1;
        s += tmp * tmp;
    }

    return sqrt(s / M);
}
ANNS_CND(A) {return u_xx[0][0] + u_xx[0][1];}
//ANNS_CND(BT) {return u[0] - 4;}
//ANNS_CND(BB) {return u[0] + 3;}
ANNS_CND(BT) {return u[0] - 1;}
ANNS_CND(BL) {return u_x[0][0];}
ANNS_CND(BR) {return u_x[0][1];}
ANNS_CND(BB) {return u[0];}

ANNS_CND_G(A_a) {return u_axx[0] + u_axx[1];}
ANNS_CND_G(BT_a) {return u_a;}
ANNS_CND_G(BL_a) {return u_ax[0];}
ANNS_CND_G(BR_a) {return u_ax[1];}
ANNS_CND_G(BB_a) {return u_a;}

void curvelinear_coil_quad(void) {
   int data[] = {
        PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 3, 1,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 2, 5, 1, 0,
        30, 15,  0, /* Уравнение переноса */
        10, 5,   0, /* Верхняя граница */
        10, 5,   0, /* Левая граница */
        10, 5,   0, /* Правая граница */
        10, 5,   0, /* Нижняя граница */
        NORMALIZED | RADIAL | GAUSSIAN, 4, 2, 0,  2, 0, 1, 1, 0,};

    double bounds[] = {
        GP_UNIFORM_COIL_SECTOR, 0, 0, C_r, C_R, 0, 0.5*M_PI,
        GP_UNIFORM_ARC, 0, 0, C_R, 0, 0.5*M_PI,
        GP_UNIFORM_SOLID_CUBOID, 0, 0, C_r, C_R,
        GP_UNIFORM_SOLID_CUBOID, C_r, C_R, 0, 0,
        GP_UNIFORM_ARC, 0, 0, C_r, 0, 0.5*M_PI,
    };
//        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, 0, M_PI,
//        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, M_PI, M_PI,
//        GP_UNIFORM_SOLID_CUBOID, 0, 0, 0, M_PI,
//        GP_UNIFORM_SOLID_CUBOID, M_PI, M_PI, 0, M_PI,
//        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, 0, 0};

    fbvp_t cond[] = {A, BT, BL, BR, BB};
    gbvp_t cond_a[] = {A_a, BT_a, BL_a, BR_a, BB_a};
    double delta[] = {1./1, 1, 1, 1, 1};

//    anns_solverdata_t *as = as_new2(ai_init, ai_newpoints, ai_initann, ai_mse, data, bounds, cond, cond_a, NULL);
    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints, ai_initann, ai_mse, NULL, data, bounds, cond, cond_a, delta);

    am->randomseed = 0;
    am->seed = 123L;
    am->Nopt = 5;
    am->gtol = 1e-4;
    am->tol = 1e-5;
    am->n0 = 16;

//    anns_solve_growreduce(as);
//    anns_solve_restarts(as);
//
//    as_delete(as);
    method_restarts(am);

    am_delete(am);
}

ANNS_CND(Aappr) {return u[0] - 0.5*log(x[0]*x[0] + x[1]*x[1]) / log(2) - 1;}
ANNS_CND_G(Aappr_a) {return u_a;}

void curvelinear_coil_quad_appr(void) {
   int data[] = {
        PJ_APPROXIMATION, PJ_ELLIPTIC_PROBLEM, 3, 1,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 2, 1, 1, 0,
        200, 15,  0,
        NORMALIZED | RADIAL | GAUSSIAN, 32, 2, 0,  0};

    double bounds[] = {
        GP_UNIFORM_COIL_SECTOR, 0, 0, C_r, C_R, 0, M_PI_2,
    };
//        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, 0, M_PI,
//        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, M_PI, M_PI,
//        GP_UNIFORM_SOLID_CUBOID, 0, 0, 0, M_PI,
//        GP_UNIFORM_SOLID_CUBOID, M_PI, M_PI, 0, M_PI,
//        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, 0, 0};

    fbvp_t cond[] = {Aappr};
    gbvp_t cond_a[] = {Aappr_a};
    double delta[] = {1./30, 1./10, 1./10, 1./10, 1./10};

//    anns_solverdata_t *as = as_new2(ai_init, ai_newpoints, ai_initann, ai_mse, data, bounds, cond, cond_a, NULL);
    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints, ai_initann, ai_mse, NULL, data, bounds, cond, cond_a, NULL);

    am->randomseed = 0;
    am->seed = 123L;
    am->Nopt = 5;
    am->gtol = 1e-5;
    am->tol = 1e-5;
    am->n0 = 16;

//    anns_solve_growreduce(as);
//    anns_solve_restarts(as);
    method_restarts(am);

    am_delete(am);
}

#undef C_r
#undef C_R
