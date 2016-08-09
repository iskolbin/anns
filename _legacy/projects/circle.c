#include "curvelinear.h"

static void ai_init(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
static void ai_newpoints(anns_instance_t *, double *bounds);
static void ai_initann(anns_instance_t *, double *);
static double ai_mse(anns_instance_t *, double *);

ANNS_CND(A) {return u_xx[0][0] + u_xx[0][1] + 4;}
ANNS_CND(G) {return u[0];}

ANNS_CND_G(A_a) {return u_axx[0] + u_axx[1];}
ANNS_CND_G(G_a) {return u_a;}


void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_newpoints(anns_instance_t *ai, double *bounds) {
    int i;

    ai_genpoints(ai, bounds);
}

#define C_R 1.0

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;
    double b = 1/sqrt(ai->nets[0]->nc), r, angle;

    for (i = 0; i < ai->nets[0]->lenc; i+=4) {

        x[i] = uniform(-2, 2.);   // Веса
//            x[i++] = uniform(0.5, 2.);   // Ширины
        x[i+1] = b;

        r = uniform(0, C_R);
        angle = uniform(0, 2*M_PI);

        x[i+2] = r * cos(angle);
        x[i+3] = r * sin(angle);
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    double s, tmp, point[2], r, angle, x2, y2;
    int i, p, M = 500;

    s = 0;
    for (p = 0; p < M; p++) {
        r = uniform(0, C_R);
        angle = uniform(0, 2*M_PI);

        point[0] = r * cos(angle);
        point[1] = r * sin(angle);
        x2 = point[0] * point[0];
        y2 = point[1] * point[1];

//        tmp = anns_eval(point, x, ai) - 3.5*log(x2 + y2) / log(2) - 4;
        tmp = anns_eval(point, x, ai) - C_R * C_R + x2 + y2;
        s += tmp * tmp;
    }

    return sqrt(s / M);
}


void curvelinear_circle(void) {
   int data[] = {
        PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 3, 3,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 2, 2, 1, 0,
        30, 15,  0, /* Уравнение переноса */
        30, 5,   0, /* Граница */
        NORMALIZED | RADIAL | GAUSSIAN, 2, 2, 0,  2, 0,};

    double bounds[] = {
        GP_UNIFORM_COIL_SECTOR, 0, 0, 0, C_R, 0, 2*M_PI,
        GP_UNIFORM_ARC, 0, 0, C_R, 0, 2*M_PI,
    };
//        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, 0, M_PI,
//        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, M_PI, M_PI,
//        GP_UNIFORM_SOLID_CUBOID, 0, 0, 0, M_PI,
//        GP_UNIFORM_SOLID_CUBOID, M_PI, M_PI, 0, M_PI,
//        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, 0, 0};

    fbvp_t cond[] = {A, G, };
    gbvp_t cond_a[] = {A_a, G_a, };
    double delta[] = {1./30, 1./10, 1./10, 1./10, 1./10};

    anns_methoddata_t *as = am_new2(ai_init, ai_newpoints, ai_initann, ai_mse, NULL, data, bounds, cond, cond_a, NULL);
//anns_methoddata_t *am_new2(void (*ai_init)(anns_instance_t *, fbvp_t *, gbvp_t *, double *),
//                            void (*ai_newpoints)(anns_instance_t *, double *),
//                            void (*ai_initann)(anns_instance_t *, double *),
//                            double (*ai_mse)(anns_instance_t *, double *),
//                            void (*ai_update)(anns_instance_t *, double *, int , double),
//                            int *data, double *bounds, fbvp_t *cond, gbvp_t *cond_a, double *delta) {

    as->randomseed = 0;
    as->seed = 123L;
    as->Nopt = 5;
    as->gtol = 1e-5;
    as->tol = 1e-5;
    as->n0 = 16;

//    anns_solve_growreduce(as);
//    anns_solve_restarts(as);
    method_restarts(as);

    as_delete(as);
}
