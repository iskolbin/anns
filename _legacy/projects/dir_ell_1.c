#include "projects.h"

static void ai_init(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
static void ai_newpoints(anns_instance_t *, double *bounds);
static void ai_initann(anns_instance_t *, double *);
static double ai_mse(anns_instance_t *, double *);

//ANNS_CND(A) {return (u_xx[0][0] + u_xx[0][1] - v[0])*(u_xx[0][0] + u_xx[0][1] - v[0]);}
//ANNS_CND(B) {return u[0]*u[0];}
//
//ANNS_CND_G(A_a) {return 2*(u_xx[0][0] + u_xx[0][1] - v[0])*(u_axx[0] + u_axx[1]);}
//ANNS_CND_G(B_a) {return 2*u_a*u[0];}

ANNS_CND(A) {return u_xx[0][0] + u_xx[0][1] - v[0];}
ANNS_CND(B) {return u[0];}

ANNS_CND_G(A_a) {return u_axx[0] + u_axx[1];}
ANNS_CND_G(B_a) {return u_a;}


ANNS_FUN(A_f) {return sin(x[0])*sin(x[1]);}

void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_newpoints(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, A_f, NULL);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (t = 0; t < ai->nnets; t++) {
        for (i = ai->nets[t]->I; i < ai->nets[t]->lenc; ) {
            x[i++] = uniform(-3., 3.);         // Веса
            x[i++] = M_PI/sqrt(ai->nets[t]->nc);
            x[i++] = uniform(-M_PI/4,M_PI+M_PI/4);       // Центры
            x[i++] = uniform(-M_PI/4,M_PI+M_PI/4);       // Центры
        }
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    double s, tmp, point[2];
    int i, p, M = 10000;

    s = 0;
    for (p = 0; p < M; p++) {
        point[0] = uniform(0, M_PI);
        point[1] = uniform(0, M_PI);

        tmp = 0;
        tmp += anns_eval(point, x, ai);
        tmp += 0.5*sin(point[0])*sin(point[1]);

        s += tmp * tmp;
    }

    return sqrt(s / (M-1) );
}

void dir_ell_1(void) {
    int data[] = {
        PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 1, 0,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 2, 5, 1, 0,
        30, 15,  1, /* Уравнение переноса */
        10, 5,   0, /* Верхняя граница */
        10, 5,   0, /* Левая граница */
        10, 5,   0, /* Правая граница */
        10, 5,   0, /* Нижняя граница */
        CLASSIC | RADIAL | GAUSSIAN, 8, 2, 0,  2, 0, 0, 0, 0,};

    double bounds[] = {
        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, 0, M_PI,
        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, M_PI, M_PI,
        GP_UNIFORM_SOLID_CUBOID, 0, 0, 0, M_PI,
        GP_UNIFORM_SOLID_CUBOID, M_PI, M_PI, 0, M_PI,
        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, 0, 0};

    fbvp_t cond[] = {A, B, B, B, B};
    gbvp_t cond_a[] = {A_a, B_a, B_a, B_a, B_a};
    double delta[] = {1./30, 1./10, 1./10, 1./10, 1./10};

    anns_solverdata_t *as = as_new2(ai_init, ai_newpoints, ai_initann, ai_mse, data, bounds, cond, cond_a, NULL);

    as->randomseed = 0;
    as->seed = 123L;
    as->Nopt = 5;
    as->gtol = 1e-4;
    as->tol = 1e-4;
    as->n0 = 8;

//    anns_solve_growreduce(as);
    anns_solve_restarts(as);

    as_delete(as);
}
