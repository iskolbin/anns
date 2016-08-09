#include "special.h"

static void ai_init(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
static void ai_initann(anns_instance_t *, double *);

#define C_k 27.79

static void ai_newpoints_bl(anns_instance_t *, double *bounds);
static double ai_mse_bl(anns_instance_t *, double *);

void ai_newpoints_bl(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
}

double ai_mse_bl(anns_instance_t *ai, double *x) {
    double s=0, tmp, point[1];
    int i, p, Np = 101;
    double dx = 1./(Np-1);
    double ek = exp(C_k), ek_ = exp(-C_k), ekx, ekx_;

    s = 0;
    printf("[s=%g ", s);
    for (p = 0; p < Np; p++) {
        point[0] =  dx * p;

        tmp = 0;
        tmp += anns_eval(point, x, ai);
        printf(" nrbf tmp=%g ", tmp);

        ekx = exp(C_k * point[0]);
        ekx_ = exp(-C_k * point[0]);

        tmp -= ((ek - 1) * ekx_ + (1 - ek_) * ekx ) / (ek - ek_);

        s += tmp * tmp;
        printf("ek=%g ek_=%g ekx=%g ekx_=%g tmp=%g s=%g tmp^2=%g]\n", ek, ek_, ekx, ekx_, tmp, s, tmp*tmp);
    }

    return sqrt(s / (Np-1) );;
}

ANNS_CND(Abl) {return u_xx[0][0] - C_k*C_k*u[0];}
ANNS_CND(Bbl) {return u[0]-1;}

ANNS_CND_G(Abl_a) {return u_axx[0] - C_k*C_k*u_a;}
ANNS_CND_G(Bbl_a) {return u_a;}

void special_bl_problem(void) {
    int data[] = {
        PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 4, 1,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 1, 2, 1, 0,
        20, 10, 0, /* Уравнение переноса */
        2, 2,   0, /* Граница граница */
        NORMALIZED | RADIAL | GAUSSIAN, 8, 1, 0,  2, 0};

    double bounds[] = {
        GP_UNIFORM_SOLID_CUBOID, 0, 1,
        GP_POINTS, 0, 1};

    fbvp_t cond[] = {Abl, Bbl};
    gbvp_t cond_a[] = {Abl_a, Bbl_a};
    double delta[] = {1e-7, 1};

    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints_bl, ai_initann, ai_mse_bl, NULL, data, bounds, cond, cond_a, delta);

    am->randomseed = 0;
//    am->seed = 12311L;
    am->seed = 123111111L;
    am->Nopt = 5;
    am->gtol = 1e-6;
    am->tol = 1e-5;
//    am->nofileout = 0;

    method_restarts(am);

    as_delete(am);
}


static void ai_newpoints_spe(anns_instance_t *, double *bounds);
static double ai_mse_spe(anns_instance_t *, double *);

#define C_LAMBDA 2.0
#define C_MU     sqrt(C_LAMBDA*C_LAMBDA + M_PI*M_PI)

ANNS_FUN(Aspe_s) {return cos(M_PI * x[0]) * sinh(C_MU * (1.0 - x[1])) / sinh(C_MU);}

void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_newpoints_spe(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);

    ai_eval_v(ai, 1, 0, Aspe_s, NULL);
    ai_eval_v(ai, 2, 0, Aspe_s, NULL);
    ai_eval_v(ai, 3, 0, Aspe_s, NULL);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i = 0, j, t, n, s;
    double b = M_PI/sqrt(ai->nets[0]->nc), d;

    srand(123L);
    for (t = 0; t < ai->nnets; t++) {
        for (i = ai->nets[t]->I; i < ai->nets[t]->lenc;) {
            x[i++] = uniform(-2., 2.);         // Веса
            x[i++] = b;
            for (j = 0; j < ai->nets[t]->dim; j++) {
                x[i++] = uniform(-0.05, 1.05);       // Центры
            }
        }

        an_setwidth_nearest(ai->nets[t], 1.5);
    }
}

ANNS_CND(Aspe) {return u_xx[0][0] + u_xx[0][1] - C_LAMBDA*C_LAMBDA*u[0];}
ANNS_CND(Bspe_l) {return u[0] - sinh(C_MU * (1-x[1])) / sinh(C_MU);}
ANNS_CND(Bspe_r) {return u[0] + sinh(C_MU * (1-x[1])) / sinh(C_MU);}
ANNS_CND(Bspe_b) {return u[0] - cos(M_PI * x[0]);}
ANNS_CND(B)      {return u[0];}

ANNS_CND_G(Aspe_a) {return u_axx[0] + u_axx[1] - C_LAMBDA*C_LAMBDA*u_a;}
ANNS_CND_G(B_a) {return u_a;}

double ai_mse_spe(anns_instance_t *ai, double *x) {
    double s = 0, tmp, point[2];
    int i, npoints = 10000;

    for (i = 0; i < npoints; i++) {
        point[0] = uniform(0.0, 1.0);
        point[1] = uniform(0.0, 1.0);

        tmp = anns_eval(point, x, ai) - Aspe_s(point, 2, NULL);
//        tmp = anns_eval(point, x, ai) - cos(M_PI * point[0]) * sinh(C_MU * (1.0 - point[1])) / sinh(C_MU);
        s += tmp * tmp;
    }

    return sqrt(s / (npoints-1) );
}

void special_spe_problem(void) {
    int data[] = {
        PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 5, 1,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 2, 5, 1, 0,
        80, 15,  0, /* Уравнение переноса */
        20, 5,   1, /* Левая граница */
        20, 5,   1, /* Правая граница */
        20, 5,   1, /* Нижняя граница */
        20, 5,   0, /* Верхняя граница */
        NORMALIZED | RADIAL | GAUSSIAN, 16, 2, 0,  2, 0, 0, 0, 0,};

    double bounds[] = {
        GP_UNIFORM_SOLID_CUBOID, 0, 1, 0, 1,
        GP_UNIFORM_SOLID_CUBOID, 0, 0, 0, 1, // Левая
        GP_UNIFORM_SOLID_CUBOID, 1, 1, 0, 1, // Правая
        GP_UNIFORM_SOLID_CUBOID, 0, 1, 0, 0, // Нижняя
        GP_UNIFORM_SOLID_CUBOID, 0, 1, 1, 1}; // Верхняя



    fbvp_t cond[] = {Aspe, Bspe_l, Bspe_r, Bspe_b, B};
    gbvp_t cond_a[] = {Aspe_a, B_a, B_a, B_a, B_a};
    double delta[] = {1./ 500, 1, 1, 1, 1};

    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints_spe, ai_initann, ai_mse_spe, NULL, data, bounds, cond, cond_a, delta);

    am->randomseed = 0;
//    am->seed = 12311L;
    am->seed = 1231111111L;
    am->Nopt = 5;
    am->gtol = 1e-5;
    am->tol = 1e-5;
//    am->nofileout = 0;

    method_restarts(am);

    as_delete(am);
}


#define C_A 1.0
#define C_a 9.0
#define C_b 4.0
#define C_c 2.0
#define C_x0 0.25
#define C_y0 0.75

ANNS_FUN(Ae_s) {
    double x_ = x[0] - C_x0, y_ = x[1] - C_y0;
    return C_A * exp( -(C_a * x_ * x_ + 2 * C_b * x_ * y_ + C_c * y_ * y_ ) );}

ANNS_FUN(Ae_src) {
    double x_ = x[0] - C_x0, y_ = x[1] - C_y0;
    double e = Ae_s(x, 2, NULL);
    double bxcy = (C_b * x_ + C_c * y_), axby = (C_a * x_ + C_b * y_);
    return 4 * e * (axby*axby + bxcy*bxcy - 0.5*(C_a + C_c));
}

ANNS_CND(Ae) {return u_xx[0][0] + u_xx[0][1] - v[0];}
ANNS_CND(Be_l) {return u[0] - v[0];}
ANNS_CND(Be_r) {return u[0] - v[0];}
ANNS_CND(Be_b) {return u[0] - v[0];}
ANNS_CND(Be_t) {return u[0] - v[0];}

ANNS_CND_G(Ae_a) {return u_axx[0] + u_axx[1];}

void ai_newpoints_e(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);

    ai_eval_v(ai, 0, 0, Ae_src, NULL);
    ai_eval_v(ai, 1, 0, Ae_s, NULL);
    ai_eval_v(ai, 2, 0, Ae_s, NULL);
    ai_eval_v(ai, 3, 0, Ae_s, NULL);
    ai_eval_v(ai, 4, 0, Ae_s, NULL);
}

double ai_mse_e(anns_instance_t *ai, double *x) {
    double s = 0, tmp, point[2];
    int i, npoints = 10000;

    for (i = 0; i < npoints; i++) {
        point[0] = uniform(0.0, 1.0);
        point[1] = uniform(0.0, 1.0);

        tmp = anns_eval(point, x, ai) - Ae_s(point, 2, NULL);
        s += tmp * tmp;
    }

    return sqrt(s / (npoints-1) );
}

void special_e_problem(void) {
    int data[] = {
        PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 6, 1,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 2, 5, 1, 0,
        150, 15,  1, /* Уравнение переноса */
        20, 5,   1, /* Левая граница */
        20, 5,   1, /* Правая граница */
        20, 5,   1, /* Нижняя граница */
        20, 5,   1, /* Верхняя граница */
        CLASSIC | RADIAL | GAUSSIAN, 8, 2, 0,  2, 0, 0, 0, 0,};

    double bounds[] = {
        GP_UNIFORM_SOLID_CUBOID, 0, 1, 0, 1,
        GP_UNIFORM_SOLID_CUBOID, 0, 0, 0, 1, // Левая
        GP_UNIFORM_SOLID_CUBOID, 1, 1, 0, 1, // Правая
        GP_UNIFORM_SOLID_CUBOID, 0, 1, 0, 0, // Нижняя
        GP_UNIFORM_SOLID_CUBOID, 0, 1, 1, 1}; // Верхняя

    fbvp_t cond[] = {Ae, Be_l, Be_r, Be_b, Be_t};
    gbvp_t cond_a[] = {Ae_a, B_a, B_a, B_a, B_a};
    double delta[] = {1./ 1000, 1, 1, 1, 1};

    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints_e, ai_initann, ai_mse_e, NULL, data, bounds, cond, cond_a, delta);

    am->randomseed = 0;
    am->seed = 123111L;
    am->Nopt = 5;
    am->gtol = 1e-5;
    am->tol = 1e-5;
    am->nofileout = 0;

    method_restarts(am);

    as_delete(am);
}

#undef C_LAMBDA
#undef C_MU
#undef C_k
#undef C_A
#undef C_a
#undef C_b
#undef C_c
#undef C_x0
#undef C_y0
