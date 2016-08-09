#include "curvelinear.h"

static void ai_init(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
static void ai_newpoints(anns_instance_t *, double *bounds);
static void ai_initann(anns_instance_t *, double *);
static double ai_mse(anns_instance_t *, double *);

#define C_a 1.0
#define C_R (1./sqrt(2.))
#define C_angle (0.25*M_PI)

ANNS_CND(A) {return u_xx[0][0] + u_xx[0][1];}
ANNS_CND(G1) {return u[0] - v[0];}
ANNS_CND(G2) {return u[0] - v[0];}
ANNS_CND(G3) {return u[0] - v[0];}

ANNS_CND_G(A_a) {return u_axx[0] + u_axx[1];}
ANNS_CND_G(G1_a) {return u_a;}
ANNS_CND_G(G2_a) {return u_a;}
ANNS_CND_G(G3_a) {return u_a;}

ANNS_FUN(G1_f) {return C_a * (0.5 - 2*x[1]*x[1]);}
ANNS_FUN(G2_f) {return C_a * x[0] * x[0];}
ANNS_FUN(G3_f) {return C_a * (0.25 - x[1] * x[1]);}


void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_newpoints(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 1, 0, G1_f, NULL);
    ai_eval_v(ai, 2, 0, G2_f, NULL);
    ai_eval_v(ai, 3, 0, G3_f, NULL);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;
    double b = 1/sqrt(ai->nets[0]->nc), r, angle;

    for (i = 0; i < ai->nets[0]->lenc; i+=4) {

        x[i] = uniform(-2, 2.);   // Веса
//            x[i++] = uniform(0.5, 2.);   // Ширины
        x[i+1] = b;

        r = uniform(0.25, C_a+0.25);
        angle = uniform(0, C_angle);

        x[i+2] = r * cos(angle);
        x[i+3] = r * sin(angle);
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    double s, tmp, point[2], r, angle, x2, y2, b;
    int i, p, M = 500;

    s = 0;
    for (p = 0; p < M; p++) {

        angle = uniform(0, C_angle);

        b = C_R * (cos(C_angle)/cos(angle) );

        r = uniform(b, C_R);

        point[0] = r * cos(angle);
        point[1] = r * sin(angle);
        x2 = point[0] * point[0];
        y2 = point[1] * point[1];

        tmp = anns_eval(point, x, ai) - C_a * (x2 - y2);
        s += tmp * tmp;
    }

    return sqrt(s / M);
}




void curvelinear_ellipse_quad(void) {
    int data[] = {
        PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 3, 1,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 2, 4, 1, 0,
        100, 15,  0, /* Уравнение переноса */
        10, 5,   1, /* Верхняя граница */
        10, 5,   1, /* Левая граница */
        10, 5,   1, /* Правая граница */
        NORMALIZED | RADIAL | GAUSSIAN, 32, 2, 0,  2, 0, 0, 0};

    double bounds[] = {
        GP_UNIFORM_CUT_SECTOR, 0, 0, C_R, 0, C_angle,
        GP_UNIFORM_ARC, 0, 0, C_R, 0, C_angle,
        GP_UNIFORM_SOLID_CUBOID, 0.5, 1./sqrt(2), 0, 0,
        GP_UNIFORM_SOLID_CUBOID, 0.5, 0.5, 0, 0.5,
    };
//        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, 0, M_PI,
//        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, M_PI, M_PI,
//        GP_UNIFORM_SOLID_CUBOID, 0, 0, 0, M_PI,
//        GP_UNIFORM_SOLID_CUBOID, M_PI, M_PI, 0, M_PI,
//        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, 0, 0};

    fbvp_t cond[] = {A, G1, G2, G3,};
    gbvp_t cond_a[] = {A_a, G1_a, G2_a, G3_a,};
    double delta[] = {1./30, 1./10, 1./10, 1./10,};

    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints, ai_initann, ai_mse, NULL, data, bounds, cond, cond_a, NULL);

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

#undef C_a
#undef C_R
#undef C_angle
