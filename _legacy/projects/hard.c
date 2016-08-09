#include "dir_ell_hard.h"

static void ai_newpoints_sf(anns_instance_t *, double *bounds);
static double ai_mse_sf(double *x, int n, void *instance);
static void ai_newpoints_spe(anns_instance_t *ai, double *bounds);
static void ai_initann(anns_instance_t *ai, double *x);
static void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta);

void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_newpoints_sf(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
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
//
//ANNS_CND(Asf) {return u_xx[0][0] + u_xx[0][1];}
//ANNS_CND(Bsf_l) {return u_x[0][0] -  0;}
////ANNS_CND(Bsf_r) {return u_x[0][0] - 1 + x[0];}
////ANNS_CND(Bsf_r) {return u_x[0][0] - 1;}
//ANNS_CND(Bsf_r) {return u_x[0][0]  -  0;}
////ANNS_CND(Bsf_t) {return u_x[0][1] - x[0] < 0.5 ? 0 : x[0] > 0.5 ? 1.0 : 0.5;}
////ANNS_CND(Bsf_t) {return u_x[0][1] - (x[0] < 0.5 ? 2*x[0] : (2 - 2*x[0]));}
//ANNS_CND(Bsf_t) {return u_x[0][1] - 1;}
//ANNS_CND(Bsf_b) {return u_x[0][1];}
//
//ANNS_CND_G(Asf_a) {return u_axx[0] + u_axx[1];}
//ANNS_CND_G(Bsf_l_a) {return u_ax[0];}
//ANNS_CND_G(Bsf_r_a) {return u_ax[0];}
//ANNS_CND_G(Bsf_b_a) {return u_ax[1];}
//ANNS_CND_G(Bsf_t_a) {return u_ax[1];}

ANNS_CND(Asf) {return u_xx[0][0] + u_xx[0][1];}
ANNS_CND(Bsf_l) {return u[0] - 0;}
ANNS_CND(Bsf_r) {return u[0] - x[1];}
//ANNS_CND(Bsf_t) {return u[0] - x[0] < 0.5 ? 0 : x[0] > 0.5 ? 1.0 : 0.5;}
//ANNS_CND(Bsf_t) {return u[0] - sqrt(x[0]);}
//ANNS_CND(Bsf_t) {return u[0] - ((x[0] < 0.5) ? 0.0 : 1.0);}
ANNS_CND(Bsf_t) {return u[0] - ((x[0] < 0.5) ? 2*x[0] : 1.0);}
ANNS_CND(Bsf_b) {return u[0] - 0;}

ANNS_CND_G(Asf_a) {return u_axx[0] + u_axx[1];}
ANNS_CND_G(Bsf_l_a) {return u_a;}
ANNS_CND_G(Bsf_r_a) {return u_a;}
ANNS_CND_G(Bsf_b_a) {return u_a;}
ANNS_CND_G(Bsf_t_a) {return u_a;}

double ai_mse_sf(double *x, int n, void *instance) {

    return -1;
}

void laplace_stepflux(void) {
    int data[] = {
        PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 6, 1,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 2, 5, 1, 0,
        40, 15,  0, /* Уравнение переноса */
        10, 5,   0, /* Левая граница */
        10, 5,   0, /* Правая граница */
        10, 5,   0, /* Нижняя граница */
        10, 5,   0, /* Верхняя граница */
        NORMALIZED | RADIAL | GAUSSIAN, 16, 2, 0,  2, 1, 1, 1, 1,};

    double bounds[] = {
        GP_UNIFORM_SOLID_CUBOID, 0, 1, 0, 1,
        GP_UNIFORM_SOLID_CUBOID, 0, 0, 0, 1, // Левая
        GP_UNIFORM_SOLID_CUBOID, 1, 1, 0, 1, // Правая
        GP_UNIFORM_SOLID_CUBOID, 0, 1, 0, 0, // Нижняя
        GP_UNIFORM_SOLID_CUBOID, 0, 1, 1, 1}; // Верхняя



    fbvp_t cond[] = {Asf, Bsf_l, Bsf_r, Bsf_b, Bsf_t};
    gbvp_t cond_a[] = {Asf_a, Bsf_l_a, Bsf_r_a, Bsf_b_a, Bsf_t_a};
    double delta[] = {1., 1, 1, 1, 1};

    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints_sf, ai_initann, ai_mse_sf, NULL, data, bounds, cond, cond_a, delta);

    am->randomseed = 0;
//    am->seed = 12311L;
    am->seed = 1231111111L;
    am->Nopt = 5;
    am->gtol = 1e-5;
    am->tol = 1e-5;
    am->nofileout = 0;

    method_restarts(am);

    as_delete(am);
}
