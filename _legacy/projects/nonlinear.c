#include "nonlinear.h"

static void ai_init(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
static void ai_newpoints(anns_instance_t *, double *bounds);
static void ai_initann(anns_instance_t *, double *);
static double ai_mse(anns_instance_t *, double *);

ANNS_CND(A) {return u_xx[0][0] + u_xx[0][1] + u_xx[0][2] - sin(u[0]*u[0]) - v[0];}
ANNS_CND(B) {return u[0]- v[0];}

ANNS_CND_G(A_a) {return u_axx[0] + u_axx[1] + u_axx[2] - 2*u[0]*u_a*cos(u[0]*u[0]);}
ANNS_CND_G(B_a) {return u_a;}

ANNS_FUN(A_f) {double v = sin(x[0] + x[1] + x[2]); return -3*v - sin(v*v);}
ANNS_FUN(B_f) {return sin(x[0] + x[1] + x[2]);}

void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_newpoints(anns_instance_t *ai, double *bounds) {
    int i;

    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, A_f, NULL);
    for (i = 1; i < 7; i++) {
        ai_eval_v(ai, i, 0, B_f, NULL);
    }
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;
    double b = M_PI/sqrt(ai->nets[0]->nc);

    i = 0;
//    srand(123L);
    for (t = 0; t < ai->nnets; t++) {
        for (i = ai->nets[t]->I; i < ai->nets[t]->lenc;) {
//            x[i++] = uniform(-2., 2.);         // Веса
////            x[i++] = M_PI/sqrt(ai->nets[t]->nc);
//            x[i++] = 1;
//            for (j = 0; j < ai->nets[t]->dim; j++) {
//                x[i++] = uniform(-M_PI/8,M_PI+M_PI/8);       // Центры
//            }
            x[i++] = uniform(-4., 4.);   // Веса
//            x[i++] = uniform(0.5, 2.);   // Ширины
            x[i++] =1;
            x[i++] = uniform(-0.5, 1.5); // Центры
            x[i++] = uniform(-0.5, 1.5);
            x[i++] = uniform(-0.5, 1.5);
        }

        an_setwidth_nearest(ai->nets[t], 1.5);
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    double s, tmp, point[3], psum;
    int i, p, M = 10000, j;

    s = 0;
    for (p = 0; p < M; p++) {
//        tmp = 1./ai->dim;
        psum = 0;
        for (j = 0; j < ai->dim; j++) {
            point[j] = uniform(0, 1);
            psum += point[j];
//            tmp *= sin(point[j]);
        }
        tmp = anns_eval(point, x, ai) - sin(psum);

        s += tmp * tmp;
    }

    return sqrt(s / M);
}
void simple_nonlinear_3d(void) {
    int data[] = {
        PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 2, 1,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 3, 7, 1, 0,
        70, 15,  1, /* Уравнение переноса */
        15, 5,   1, /* Верхняя граница */
        15, 5,   1, /* Левая граница */
        15, 5,   1, /* Правая граница */
        15, 5,   1, /* Нижняя граница */
        15, 5,   1, /* Правая граница */
        15, 5,   1, /* Нижняя граница */
        NORMALIZED | RADIAL | GAUSSIAN, 2, 3, 0,  2, 0, 0, 0, 0, 0, 0};

    double bounds[] = {
        GP_UNIFORM_SOLID_CUBOID, 0, 1,     0, 1,     0, 1,
        GP_UNIFORM_SOLID_CUBOID, 0, 1,     0, 1,     0, 0,
        GP_UNIFORM_SOLID_CUBOID, 0, 1,     0, 1,     1, 1,
        GP_UNIFORM_SOLID_CUBOID, 0, 1,     0, 0,     0, 1,
        GP_UNIFORM_SOLID_CUBOID, 0, 1,     1, 1,     0, 1,
        GP_UNIFORM_SOLID_CUBOID, 0, 0,     0, 1,     0, 1,
        GP_UNIFORM_SOLID_CUBOID, 1, 1,     0, 1,     0, 1};

    fbvp_t cond[] = {A, B, B, B, B, B, B};
    gbvp_t cond_a[] = {A_a, B_a, B_a, B_a, B_a, B_a, B_a};
    double delta[] = {1./120, 100./120, 100./120, 100./120, 100./120, 100./120, 100./120};

    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints, ai_initann, ai_mse, NULL, data, bounds, cond, cond_a, NULL);

//    am->nofileout = 0;
    am->randomseed = 0;
    am->seed = 123L;
    am->Nopt = 5;
    am->gtol = 1e-5;
    am->tol = 1e-5;
    am->n0 = 4;

//    anns_solve_growreduce(as);
    method_restarts(am);

    as_delete(am);
}
