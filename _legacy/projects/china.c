#include "china.h"

static void ai_init(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
static void ai_newpoints(anns_instance_t *, double *bounds);
static void ai_initann(anns_instance_t *, double *);
static double ai_mse(anns_instance_t *, double *);

ANNS_CND(Asines) {return u_xx[0][0] + u_xx[0][1] - v[0];}
ANNS_CND(Bsines) {return u[0];}

ANNS_CND_G(Asines_a) {return u_axx[0] + u_axx[1];}
ANNS_CND_G(Bsines_a) {return u_a;}

ANNS_FUN(Asines_f) {return sin(M_PI*x[0]) * sin(M_PI*x[1]);}

void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_newpoints(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, Asines_f, NULL);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i = 0, j, t, n, s;
    double b = M_PI/sqrt(ai->nets[0]->nc), d;

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

#define M_POINTS 5
#define STEP (1.0/(M_POINTS-1))

double ai_mse(anns_instance_t *ai, double *x) {
    double s, tmp, point[3], exact_sum = 0;
    int i, p, M = 10000, j;

    s = 0;
    for (i = 1; i < M_POINTS-1; i++) {
        for (j = 1; j < M_POINTS-1; j++) {
            point[0] = STEP * i;
            point[1] = STEP * j;

            tmp = sin(M_PI*point[0]) * sin(M_PI*point[1]) / (2 * M_PI * M_PI);
//            printf("%f ", tmp);
            exact_sum += tmp * tmp;

            tmp += anns_eval(point, x, ai);
            s += tmp * tmp;
        }
    }

//    for (p = 0; p < M; p++) {
//        for (j = 0; j < ai->dim; j++) {
//            point[j] = uniform(0, 1);
//        }
//        tmp = sin(M_PI*point[0]) * sin(M_PI*point[1]) / (2 * M_PI * M_PI);
//        exact_sum += tmp * tmp;
//
//        tmp += anns_eval(point, x, ai);
//        s += tmp * tmp;
//    }
printf("\n%g(%g) %g(%g) %g\n", s, sqrt(s), exact_sum, sqrt(exact_sum), sqrt(s) / sqrt(exact_sum));
//    return sqrt(s) / sqrt(exact_sum);
return sqrt(s/(M_POINTS*M_POINTS)) ;
}

void china_sines(void) {
    int ni = 5, nb = 5;

    int data[] = {
        PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 1, 2,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 2, 5, 1, 0,
        ni*ni, ni*ni,  1, /* Уравнение переноса */
        nb, nb,   0, /* Верхняя граница */
        nb, nb,   0, /* Левая граница */
        nb, nb,   0, /* Правая граница */
        nb, nb,   0, /* Нижняя граница */
        NORMALIZED | RADIAL | GAUSSIAN, 16, 2, 0,  2, 0, 0, 0, 0,};

//    double bounds[] = {
//        GP_UNIFORM_SOLID_CUBOID, 0, 1, 0, 1,
//        GP_UNIFORM_SOLID_CUBOID, 0, 1, 1, 1,
//        GP_UNIFORM_SOLID_CUBOID, 0, 0, 0, 1,
//        GP_UNIFORM_SOLID_CUBOID, 1, 1, 0, 1,
//        GP_UNIFORM_SOLID_CUBOID, 0, 1, 0, 0};
    double bounds[] = {
        GP_GRID_2D, ni, 0, 1,   ni, 0, 1,
        GP_GRID_2D, nb, 0, 1,   1, 1, 1,
        GP_GRID_2D, 1, 0, 0,   nb, 0, 1,
        GP_GRID_2D, 1, 1, 1,   nb, 0, 1,
        GP_GRID_2D, nb, 0, 1,   1, 0, 0};

    fbvp_t cond[] = {Asines, Bsines, Bsines, Bsines, Bsines};
    gbvp_t cond_a[] = {Asines_a, Bsines_a, Bsines_a, Bsines_a, Bsines_a};
    double delta[] = {1e-1, 1, 1, 1, 1};


    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints, ai_initann, ai_mse, NULL, data, bounds, cond, cond_a, delta);

    am->randomseed = 0;
    am->seed = 123111111L;
    am->gtol = 1e-6;
    am->tol = 1e-5;
    am->n0 = 8;

    method_restarts(am);

    as_delete(am);
}
