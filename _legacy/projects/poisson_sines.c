#include "poisson_sines.h"

static void ai_init(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
static void ai_newpoints(anns_instance_t *, double *bounds);
static void ai_initann(anns_instance_t *, double *);
static double ai_mse(double *x, int n, void *instance);

ANNS_CND(A1d) {return u_xx[0][0] - v[0];}
ANNS_CND(A2d) {return u_xx[0][0] + u_xx[0][1] - v[0];}
ANNS_CND(A3d) {return u_xx[0][0] + u_xx[0][1] + u_xx[0][2] - v[0];}
ANNS_CND(B) {return u[0];}

ANNS_CND_G(A1d_a) {return u_axx[0];}
ANNS_CND_G(A2d_a) {return u_axx[0] + u_axx[1];}
ANNS_CND_G(A3d_a) {return u_axx[0] + u_axx[1] + u_axx[2];}
ANNS_CND_G(B_a) {return u_a;}

ANNS_FUN(A1d_f) {return sin(x[0]);}
ANNS_FUN(A2d_f) {return sin(x[0]) * sin(x[1]);}
ANNS_FUN(A3d_f) {return sin(x[0]) * sin(x[1]) * sin(x[2]);}

void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_newpoints(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);

    switch (ai->dim) {
        case 1: ai_eval_v(ai, 0, 0, A1d_f, NULL); break;
        case 2: ai_eval_v(ai, 0, 0, A2d_f, NULL); break;
        case 3: ai_eval_v(ai, 0, 0, A3d_f, NULL); break;
    }
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;
    double b = M_PI/sqrt(ai->nets[0]->nc);

    i = 0;
//    x[i++] = -0.5; x[i++] = b; x[i++] = 0;      x[i++] = 0;      x[i++] = 0;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI/2; x[i++] = 0;      x[i++] = 0;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI;   x[i++] = 0;      x[i++] = 0;
//    x[i++] = -0.5; x[i++] = b; x[i++] = 0;      x[i++] = M_PI/2; x[i++] = 0;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI/2; x[i++] = M_PI/2; x[i++] = 0;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI;   x[i++] = M_PI/2; x[i++] = 0;
//    x[i++] = -0.5; x[i++] = b; x[i++] = 0;      x[i++] = M_PI;   x[i++] = 0;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI/2; x[i++] = M_PI;   x[i++] = 0;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI;   x[i++] = M_PI;   x[i++] = 0;
//
//
//    x[i++] = -0.5; x[i++] = b; x[i++] = 0;      x[i++] = 0;      x[i++] = M_PI/2;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI/2; x[i++] = 0;      x[i++] = M_PI/2;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI;   x[i++] = 0;      x[i++] = M_PI/2;
//    x[i++] = -0.5; x[i++] = b; x[i++] = 0;      x[i++] = M_PI/2; x[i++] = M_PI/2;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI/2; x[i++] = M_PI/2; x[i++] = M_PI/2;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI;   x[i++] = M_PI/2; x[i++] = M_PI/2;
//    x[i++] = -0.5; x[i++] = b; x[i++] = 0;      x[i++] = M_PI;   x[i++] = M_PI/2;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI/2; x[i++] = M_PI;   x[i++] = M_PI/2;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI;   x[i++] = M_PI;   x[i++] = M_PI/2;
//
//    x[i++] = -0.5; x[i++] = b; x[i++] = 0;      x[i++] = 0;      x[i++] = M_PI;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI/2; x[i++] = 0;      x[i++] = M_PI;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI;   x[i++] = 0;      x[i++] = M_PI;
//    x[i++] = -0.5; x[i++] = b; x[i++] = 0;      x[i++] = M_PI/2; x[i++] = M_PI;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI/2; x[i++] = M_PI/2; x[i++] = M_PI;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI;   x[i++] = M_PI/2; x[i++] = M_PI;
//    x[i++] = -0.5; x[i++] = b; x[i++] = 0;      x[i++] = M_PI;   x[i++] = M_PI;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI/2; x[i++] = M_PI;   x[i++] = M_PI;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI;   x[i++] = M_PI;   x[i++] = M_PI;

//    x[i++] = -0.5; x[i++] = b; x[i++] = 0;      x[i++] = 0;      x[i++] = 0;
//    x[i++] = -0.5; x[i++] = b; x[i++] = 0;      x[i++] = 0;      x[i++] = M_PI;
//    x[i++] = -0.5; x[i++] = b; x[i++] = 0;      x[i++] = M_PI;   x[i++] = 0;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI;   x[i++] = 0;      x[i++] = M_PI;
//    x[i++] = -0.5; x[i++] = b; x[i++] = 0;      x[i++] = M_PI;   x[i++] = M_PI;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI;   x[i++] = M_PI;   x[i++] = 0;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI;   x[i++] = 0;      x[i++] = M_PI;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI;   x[i++] = M_PI;   x[i++] = M_PI;
//    x[i++] = -0.5; x[i++] = b; x[i++] = M_PI/2; x[i++] = M_PI/2; x[i++] = M_PI/2;
    srand(123L);
    for (t = 0; t < ai->nnets; t++) {
        for (i = ai->nets[t]->I; i < ai->nets[t]->lenc;) {
            x[i++] = uniform(-2., 2.);         // Веса
//            x[i++] = M_PI/sqrt(ai->nets[t]->nc);
            x[i++] = b;
            for (j = 0; j < ai->nets[t]->dim; j++) {
                x[i++] = uniform(-M_PI/8,M_PI+M_PI/8);       // Центры
            }
        }
    }
}

double ai_mse(double *x, int n, void *instance) {
    anns_instance_t *ai = instance;
    double s, tmp, point[] = {0, 0, 0}, h;
    int i, p, M, j, k, grid = 1;

    s = 0;
    if (grid) {
        M = 100;

        h = M_PI / M;
        switch (ai->dim) {
            case 1: {
                for (i = 0; i < M + 1; i++) {
                    point[0] = i * h;
                    tmp = anns_eval(point, x, ai) + sin(point[0]);
                    s += tmp * tmp;
                }
                s = sqrt(s / M);
                break;
            }
            case 2: {
                for (i = 0; i < M + 1; i++) {
                    point[0] = i * h;
                    for (j = 0; j < M + 1; j++) {
                        point[1] = j * h;
                        tmp = anns_eval(point, x, ai) + 0.5 * sin(point[0]) * sin(point[1]);
                        s += tmp * tmp;
                    }
                }
                s = sqrt(s / (M * M));
                break;
            }
            case 3: {
                for (i = 0; i < M + 1; i++) {
                    point[0] = i * h;
                    for (j = 0; j < M + 1; j++) {
                        point[1] = j * h;
                        for (k = 0; k < M + 1; k++) {
                            point[2] = k * h;
                            tmp = anns_eval(point, x, ai) + (1./3.) * sin(point[0]) * sin(point[1]) * sin(point[2]);
                            s += tmp * tmp;
                        }
                    }
                }
                s = sqrt(s / (M * M * M));
                break;
            }
        }
    } else {
        M = 10000;
        for (p = 0; p < M; p++) {
            tmp = 1./ai->dim;
            for (j = 0; j < ai->dim; j++) {
                point[j] = uniform(0, M_PI);
                tmp *= sin(point[j]);
            }
            tmp += anns_eval(point, x, ai);

            s += tmp * tmp;
        }
        s = sqrt(s / M);
    }

    return s;
}

void poisson_1d_sines(void) {
    int data[] = {
        PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 1, 1,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 1, 3, 1, 0,
        50, 15, 1, /* Уравнение переноса */
        1, 0,   0, /* Верхняя граница */
        1, 0,   0, /* Левая граница */
        CLASSIC | RADIAL | GAUSSIAN, 6, 1, 0,  2, 0, 0};

    double bounds[] = {
        GP_UNIFORM_SOLID_CUBOID, 0, M_PI,
        GP_POINTS, 0,
        GP_POINTS, M_PI};

    fbvp_t cond[] = {A1d, B, B};
    gbvp_t cond_a[] = {A1d_a, B_a, B_a};
    double delta[] = {1./30, 1./10, 1./10, 1./10, 1./10};

    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints, ai_initann, ai_mse, NULL, data, bounds, cond, cond_a, NULL);

    am->randomseed = 0;
    am->seed = 123111111L;
//    am->Nopt = 5;
    am->gtol = 1e-5;
    am->tol = 1e-5;
    am->n0 = 8;

//    anns_solve_growreduce(as);
    method_restarts(am);
//    anns_solve_grow(as);

    as_delete(am);
}

void poisson_2d_sines(void) {
    int data[] = {
        PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 1, 2,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 2, 5, 1, 0,
        60, 15,  1, /* Уравнение переноса */
        15, 5,   0, /* Верхняя граница */
        15, 5,   0, /* Левая граница */
        15, 5,   0, /* Правая граница */
        15, 5,   0, /* Нижняя граница */
        NORMALIZED | RADIAL | GAUSSIAN, 8, 2, 0,  2, 0, 0, 0, 0,};

    double bounds[] = {
        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, 0, M_PI,
        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, M_PI, M_PI,
        GP_UNIFORM_SOLID_CUBOID, 0, 0, 0, M_PI,
        GP_UNIFORM_SOLID_CUBOID, M_PI, M_PI, 0, M_PI,
        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, 0, 0};

    fbvp_t cond[] = {A2d, B, B, B, B};
    gbvp_t cond_a[] = {A2d_a, B_a, B_a, B_a, B_a};
    double delta[] = {1./30, 1./10, 1./10, 1./10, 1./10};

//    anns_solverdata_t *as = as_new2(ai_init, ai_newpoints, ai_initann, ai_mse, data, bounds, cond, cond_a, NULL);
//
//    as->randomseed = 0;
//    as->seed = 123L;
//    as->Nopt = 5;
//    as->gtol = 1e-5;
//    as->tol = 1e-5;
//    as->n0 = 8;
    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints, ai_initann, ai_mse, NULL, data, bounds, cond, cond_a, NULL);

    am->randomseed = 0;
    am->seed = 1231111L;
//    am->Nopt = 5;
    am->gtol = 1e-5;
    am->tol = 1e-5;
    am->n0 = 8;
    am->nofileout = 0;

//    anns_solve_growreduce(as);
    method_restarts(am);
//    anns_solve_grow(as);

    as_delete(am);
}

void poisson_3d_sines(void) {
    int data[] = {
        PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 1, 3,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 3, 7, 1, 0,
        100, 15,  1, /* Уравнение переноса */
        10, 5,   0, /* Верхняя граница */
        10, 5,   0, /* Левая граница */
        10, 5,   0, /* Правая граница */
        10, 5,   0, /* Нижняя граница */
        10, 5,   0, /* Правая граница */
        10, 5,   0, /* Нижняя граница */
        CLASSIC | RADIAL | GAUSSIAN, 16, 3, 0,  2, 0, 0, 0, 0, 0, 0};

    double bounds[] = {
        GP_UNIFORM_SOLID_CUBOID, 0, M_PI,     0, M_PI,     0, M_PI,
        GP_UNIFORM_SOLID_CUBOID, 0, M_PI,     0, M_PI,     0, 0,
        GP_UNIFORM_SOLID_CUBOID, 0, M_PI,     0, M_PI,     M_PI, M_PI,
        GP_UNIFORM_SOLID_CUBOID, 0, M_PI,     0, 0,        0, M_PI,
        GP_UNIFORM_SOLID_CUBOID, 0, M_PI,     M_PI, M_PI,  0, M_PI,
        GP_UNIFORM_SOLID_CUBOID, 0, 0,        0, M_PI,     0, M_PI,
        GP_UNIFORM_SOLID_CUBOID, M_PI, M_PI,  0, M_PI,     0, M_PI};

    fbvp_t cond[] = {A3d, B, B, B, B, B, B};
    gbvp_t cond_a[] = {A3d_a, B_a, B_a, B_a, B_a, B_a, B_a};
    double w[] = {1./1, 1, 1, 1, 1, 1, 1};

    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints, ai_initann, ai_mse, NULL, data, bounds, cond, cond_a, w);

    am->randomseed = 0;
    am->seed = 123111111L;
//    am->Nopt = 5;
    am->gtol = 1e-5;
    am->tol = 1e-5;
    am->n0 = 8;

//    anns_solve_growreduce(as);
    method_restarts(am);
//    anns_solve_grow(as);

    as_delete(am);
}


