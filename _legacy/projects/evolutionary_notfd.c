#include "evolutionary_notfd.h"

static double static_Tm = 1.0;

#define C_a 1.0
//#define C_x (2*M_PI)
#define C_x 1.0

ANNS_CND(A1d) {return u_x[0][1] - C_a*u_xx[0][0];}
ANNS_CND_G(A1d_a) {return u_ax[1] - C_a*u_axx[0];}

ANNS_CND(B1d) {return u[0];}
ANNS_CND_G(B1d_a) {return u_a;}

ANNS_CND(I1d) {return u[0] - v[0];}
ANNS_CND_G(I1d_a) {return u_a;}
ANNS_FUN(I1d_f) {return sin(C_x * x[0]);}

void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_newpoints(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 3, 0, I1d_f, NULL);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i = 0, j, t, n, s;
    double b = M_PI/sqrt(ai->nets[0]->nc);

    for (t = 0; t < ai->nnets; t++) {
        for (i = ai->nets[t]->I; i < ai->nets[t]->lenc;) {
            x[i++] = uniform(-2., 2.);         // Веса
            x[i++] = b;
            x[i++] = uniform(-0.1, 2*M_PI + 0.1);       // Центры x
            x[i++] = uniform(-0.05*static_Tm, 1.05*static_Tm);       // Центры t
        }

        an_setwidth_nearest(ai->nets[t], 1.5);
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    int i, M = 10000;
    double point[2], tmp, s = 0.0, c = 2*M_PI;

    for (i = 0; i < M; i++) {
        point[0] = uniform(0, 2*M_PI / C_x);
        point[1] = uniform(0, static_Tm);
        tmp = anns_eval(point, x, ai) - exp(-C_a*C_x*C_x*point[1])*sin(C_x * point[0]);
        s += tmp * tmp;
    }
    s = sqrt(s / (M-1));
    return s;
}


void evolutionary_1d_notfd(double Tm) {
    int data[] = {
        PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 1, 2,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 2, 4, 1, 0,
        100, 10,  0, /* Уравнение переноса */
        10, 10,  0, /* Граница (0,t) */
        10, 10,  0, /* Граница (2PI, t) */
        40, 10,  1, /* Начальное условие */
        NORMALIZED | RADIAL | GAUSSIAN, 8, 2, 0,  2, 0, 0, 0};

    double bounds[] = {
//        GP_UNIFORM_SOLID_CUBOID, 0, 2*M_PI, 0, Tm,
//        GP_UNIFORM_SOLID_CUBOID, 0, 0,      0, Tm,
//        GP_UNIFORM_SOLID_CUBOID, 2*M_PI, 2*M_PI,  0, Tm,
//        GP_UNIFORM_SOLID_CUBOID, 0, 2*M_PI, 0, 0,
        GP_GRID_2D, 10, 0, 2*M_PI,      10, 0, Tm,
        GP_GRID_2D, 1,  0, 0,           10, 0, Tm,
        GP_GRID_2D, 1,  2*M_PI, 2*M_PI, 10, 0, Tm,
        GP_GRID_2D, 40, 0, 2*M_PI,      1, 0, 0,
    };

    fbvp_t cond[] = {A1d, B1d, B1d, I1d};
    gbvp_t cond_a[] = {A1d_a, B1d_a, B1d_a, I1d_a};

    static_Tm = Tm;

    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints, ai_initann, ai_mse, NULL,
                                    data, bounds,
                                    cond, cond_a, NULL);


    am->randomseed = 0;
    am->seed = 1231L;
    am->Nopt = 5;
    am->gtol = 1e-5;
    am->tol = 1e-5;
//    am->n0 = 4;
    am->nofileout = 0;

    method_restarts(am);

    am_delete(am);
}

ANNS_CND(A2d) {return u_x[0][2] - C_a*(u_xx[0][0] + u_xx[0][1]);}
ANNS_CND(B2d_left) {return u[0] - v[0];}
ANNS_CND(B2d_right) {return u_x[0][0] - v[0];}
ANNS_CND(B2d_bottom) {return u_x[0][1] - v[0];}
ANNS_CND(B2d_top) {return u[0] - v[0];}
ANNS_CND(I2d) {return u[0] - v[0];}

ANNS_CND_G(A2d_a) {return u_ax[2] - C_a*(u_axx[0] + u_axx[1]);}
ANNS_CND_G(B2d_left_a) {return u_a;}
ANNS_CND_G(B2d_right_a) {return u_ax[0];}
ANNS_CND_G(B2d_bottom_a) {return u_ax[1];}
ANNS_CND_G(B2d_top_a) {return u_a;}
ANNS_CND_G(I2d_a) {return u_a;}

ANNS_FUN(B2d_left_f) {return sinh(x[1]) * exp(-3*C_a*x[2]);}
ANNS_FUN(B2d_right_f) {return -2*sinh(x[1]) * exp(-3*C_a*x[2]);}
ANNS_FUN(B2d_bottom_f) {return cos(2*x[0])*exp(-3*C_a*x[2]);}
ANNS_FUN(B2d_top_f) {return 0.75*cos(2*x[0])*exp(-3*C_a*x[2]);}
ANNS_FUN(I2d_f) {return cos(2*x[0])*sinh(x[1]);}

void ai_newpoints2d(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 1, 0, B2d_left_f, NULL);
    ai_eval_v(ai, 2, 0, B2d_right_f, NULL);
    ai_eval_v(ai, 3, 0, B2d_bottom_f, NULL);
    ai_eval_v(ai, 4, 0, B2d_top_f, NULL);
    ai_eval_v(ai, 5, 0, I2d_f, NULL);
}

void ai_initann2d(anns_instance_t *ai, double *x) {
    int i = 0, j, t, n, s;
    double b = M_PI/sqrt(ai->nets[0]->nc);

    for (t = 0; t < ai->nnets; t++) {
        for (i = ai->nets[t]->I; i < ai->nets[t]->lenc;) {
            x[i++] = uniform(-2., 2.);         // Веса
            x[i++] = b;
            x[i++] = uniform(-0.1, 1 + 0.1);       // Центры x
            x[i++] = uniform(-0.05, log(2) + 0.05);       // Центры y
            x[i++] = uniform(-0.05*static_Tm, 1.05*static_Tm);       // Центры t
        }

        an_setwidth_nearest(ai->nets[t], 1.5);
    }
}

double ai_mse2d(anns_instance_t *ai, double *x) {
    int i, M = 10000;
    double point[3], tmp, s = 0.0, c = 2*M_PI;

    for (i = 0; i < M; i++) {
        point[0] = uniform(0, M_PI/4);
        point[1] = uniform(0, log(2));
        point[2] = uniform(0, static_Tm);
        tmp = anns_eval(point, x, ai) - cos(2*point[0])*sinh(point[1])*exp(-3*C_a*point[2]);
        s += tmp * tmp;
    }
    s = sqrt(s / (M-1));
    return s;
}


void evolutionary_2d_notfd(double Tm) {
    int na = 5, nb = 5, ni = 5;

    int data[] = {
        PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 2, 2,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 3, 6, 1, 0,
        na*na*na, na*na*na,  0, /* Уравнение переноса */
        nb*nb, nb*nb,  1, /* Граница (0,y,t) */
        nb*nb, nb*nb,  1, /* Граница (pi/4,y,t)*/
        nb*nb, nb*nb,  1, /* Граница (x,0,t)*/
        nb*nb, nb*nb,  1, /* Граница (x,ln2,t)*/
        ni*ni, ni*ni,  1, /* Начальное условие (x,y,0)*/
        NORMALIZED | RADIAL | GAUSSIAN, 32, 3, 0,  2, 0, 1, 1, 0, 0};

    double bounds[] = {
        GP_GRID_3D, na, 0, M_PI/4,      na, 0, log(2),   na, 0, Tm,

        GP_GRID_3D, 1, 0, 0,           nb, 0, log(2),   nb, 0, Tm,
        GP_GRID_3D, 1,M_PI/4,M_PI/4,   nb, 0, log(2),   nb, 0, Tm,

        GP_GRID_3D, nb, 0, M_PI/4,     1, 0, 0,             nb, 0, Tm,
        GP_GRID_3D, nb, 0, M_PI/4,     1, log(2), log(2),   nb, 0, Tm,

        GP_GRID_3D, ni, 0, M_PI/4,     ni, 0, log(2),      1, 0, 0,
    };

    fbvp_t cond[] = {A2d, B2d_left, B2d_right, B2d_bottom, B2d_top, I2d};
    gbvp_t cond_a[] = {A2d_a, B2d_left_a, B2d_right_a, B2d_bottom_a, B2d_top_a, I2d_a};

    static_Tm = Tm;

    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints2d, ai_initann2d, ai_mse2d, NULL,
                                    data, bounds,
                                    cond, cond_a, NULL);


    am->randomseed = 0;
    am->seed = 1231L;
    am->Nopt = 5;
    am->gtol = 1e-4;
    am->tol = 1e-5;
//    am->n0 = 4;
    am->nofileout = 0;

    method_restarts(am);

    am_delete(am);
}

#undef C_A
#undef C_x
