#include "evolutionary.h"

//#define C_x (2*M_PI)
#define C_x 1.0
#define C_a 1.0

static void ai_init(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
static void ai_newpoints_1d(anns_instance_t *, double *bounds);
static void ai_initann_1d(anns_instance_t *, double *);
static double ai_mse(anns_instance_t *, double *);

ANNS_FUN(A1d_i_f) {return sin(C_x*x[0]);}

ANNS_CND(A1d_i) {return u[0] - v[0];}
ANNS_CND_G(A1d_i_a) {return u_a;}

ANNS_CND(A1d) {anns_instance_t *ai = instance; return .5*(ai->step)*C_a*u_xx[0][0] - u[0] + v[0];}
ANNS_CND(B1d) {return u[0];}

ANNS_CND_G(A1d_a) {anns_instance_t *ai = instance; return .5*(ai->step)*C_a*u_axx[0] - u_a;}
ANNS_CND_G(B1d_a) {return u_a;}


static void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

static void ai_newpoints_1d(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, A1d_i_f, NULL);
}


static void ai_initann_1d(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;
    double b = 1/sqrt(ai->nets[0]->nc);

    for (i = 0; i < ai->nets[0]->lenc; i+=3) {
        x[i] = uniform(-2, 2.);   // Веса
        x[i+1] = b;
        x[i+2] = uniform(-0.05, 2*M_PI / C_x + 0.05);
    }
}

static void ai_update_1d(anns_instance_t *ai, double *x, int i, double tau) {
    int p;
    anns_cache_t *ac = ai->cache;
    // Сохраняем сеть i-1 шага в контейнер
    ai_put(ai, x, i-1);

    // Обновляем предысторию (v[0]) для конечно-разностной схемы (явно-неявная, Кранк-Никлсон)
    for (p = 0; p < ai->conds[0]->m; p++) {
        ai->nets[0]->ann_val[2](x, 0, p, 0, 0, ai);
        ac->v[0][p][0] = 0.5*tau*C_a*ac->u_xx[0][0] + ac->u[0];
    }
}

void evolutionary_1d(double Tm, int fdsteps) {
    int na = 30;

    int data[] = {
        PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 1, 1,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 1, 2, 1, fdsteps,
        na, 10,  1, /* Уравнение переноса */
        2, 2,   0, /* Границы */
        NORMALIZED | RADIAL | GAUSSIAN,2, 1, 0,  2, 0};

    double bounds[] = {
        GP_GRID_1D, na, 0, 2*M_PI/C_x,
        GP_POINTS, 0, 2*M_PI/C_x,
    };

    fbvp_t cond[] = {A1d, B1d};
    gbvp_t cond_a[] = {A1d_a, B1d_a};

    fbvp_t cond_i[] = {A1d_i, B1d};
    gbvp_t cond_i_a[] = {A1d_i_a, B1d_a};

    anns_methoddata_t *am = am_new3(ai_init, ai_newpoints_1d, ai_initann_1d, NULL, ai_update_1d,
                                    data, bounds,
                                    cond, cond_a, NULL,
                                    cond_i, cond_i_a, 0, Tm, fdsteps);

    am->randomseed = 0;
    am->seed = 123L;
    am->Nopt = 5;
    am->gtol = 1e-5;
    am->tol = 1e-5;
    am->n0 = 4;
    am->nofileout = 0;

    method_hybrid(am);

    am_delete(am);
}

/*
    Двумер
*/
ANNS_FUN(A2d_i_f) {return cos(2*x[0])*sinh(x[1]);}

ANNS_CND(A2d_i) {return u[0] - v[0];}
ANNS_CND_G(A2d_i_a) {return u_a;}

ANNS_FUN(B2d_left_f) {double *t = instance; return sinh(x[1]) * exp(-3*C_a* *t);}
ANNS_FUN(B2d_right_f) {double *t = instance; return -2*sinh(x[1]) * exp(-3*C_a* *t);}
ANNS_FUN(B2d_bottom_f) {double *t = instance; return cos(2*x[0]) * exp(-3*C_a* *t);}
ANNS_FUN(B2d_top_f) {double *t = instance; return .75 * cos(2*x[0]) * exp(-3*C_a* *t);}

ANNS_CND(A2d) {anns_instance_t *ai = instance; return .5*(ai->step)*C_a*(u_xx[0][0] + u_xx[0][1]) - u[0] + v[0];}
ANNS_CND(B2d_left) {return u[0] - v[0];}
ANNS_CND(B2d_right) {return u_x[0][0] - v[0];}
ANNS_CND(B2d_bottom) {return u_x[0][1] - v[0];}
ANNS_CND(B2d_top) {return u[0] - v[0];}

ANNS_CND_G(A2d_a) {anns_instance_t *ai = instance; return .5*(ai->step)*C_a*(u_axx[0] + u_axx[1]) - u_a;}
ANNS_CND_G(B2d_left_a) {return u_a;}
ANNS_CND_G(B2d_right_a) {return u_ax[0];}
ANNS_CND_G(B2d_bottom_a) {return u_ax[1];}
ANNS_CND_G(B2d_top_a) {return u_a;}

static void ai_update_2d(anns_instance_t *ai, double *x, int i, double tau) {
    int p;
    double t = i * tau;
    anns_cache_t *ac = ai->cache;
    // Сохраняем сеть i-1 шага в контейнер
    ai_put(ai, x, i-1);

    // Обновляем предысторию (v[0]) для конечно-разностной схемы (явно-неявная, Кранк-Никлсон)
    for (p = 0; p < ai->conds[0]->m; p++) {
        ai->nets[0]->ann_val[2](x, 0, p, 0, 0, ai);
        ac->v[0][p][0] = .5*tau*C_a*(ac->u_xx[0][0] + ac->u_xx[0][1])+ ac->u[0];
    }

    ai_eval_v(ai, 1, 0, B2d_left_f, &t);
    ai_eval_v(ai, 2, 0, B2d_right_f, &t);
    ai_eval_v(ai, 3, 0, B2d_bottom_f, &t);
    ai_eval_v(ai, 4, 0, B2d_top_f, &t);
}

static void ai_newpoints_2d(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, A2d_i_f, NULL);
    ai_eval_v(ai, 1, 0, A2d_i_f, NULL);
    ai_eval_v(ai, 2, 0, A2d_i_f, NULL);
    ai_eval_v(ai, 3, 0, A2d_i_f, NULL);
    ai_eval_v(ai, 4, 0, A2d_i_f, NULL);
}

static void ai_initann_2d(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;
    double b = 1/sqrt(ai->nets[0]->nc);

    for (i = 0; i < ai->nets[0]->lenc;) {
        x[i++] = uniform(-2, 2.);   // Веса
        x[i++] = b;
        x[i++] = uniform(-0.05, M_PI/4+ 0.05);
        x[i++] = uniform(-0.05, log(2) + 0.05);
    }
    an_setwidth_nearest(ai->nets[0], 1.5);
}


void evolutionary_2d(double Tm, int fdsteps) {
    int na = 6, nb = 10, ni = 10;

    int data[] = {
        PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 2, 1,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 2, 5, 1, fdsteps,
        na*na, 40, 1, /* Уравнение переноса */
        nb, 10, 1, /* Граница u(0, y)    */
        nb, 10, 1, /* Граница u(pi/4, y) */
        nb, 10, 1, /* Граница u(x, 0)    */
        nb, 10, 1, /* Граница u(x, ln2)  */

        NORMALIZED | RADIAL | GAUSSIAN, 8, 2, 0,  2, 0, 1, 1, 0};

    double bounds[] = {
        GP_GRID_2D, na, 0, M_PI/4,  na, 0, log(2),
        GP_GRID_2D, 1,  0, 0,       nb, 0, log(2),
        GP_GRID_2D, 1,M_PI/4,M_PI/4,nb, 0, log(2),
        GP_GRID_2D, nb, 0, M_PI/4,  1, 0, 0,
        GP_GRID_2D, nb, 0, M_PI/4,  1, log(2), log(2),
    };

    fbvp_t cond[] = {A2d, B2d_left, B2d_right, B2d_bottom, B2d_top};
    gbvp_t cond_a[] = {A2d_a, B2d_left_a, B2d_right_a, B2d_bottom_a, B2d_top_a};

    fbvp_t cond_i[] = {A2d_i, A2d_i, A2d_i, A2d_i, A2d_i};
    gbvp_t cond_i_a[] = {A2d_i_a, A2d_i_a, A2d_i_a,A2d_i_a, A2d_i_a};

    anns_methoddata_t *am = am_new3(ai_init, ai_newpoints_2d, ai_initann_2d, NULL, ai_update_2d,
                                    data, bounds,
                                    cond, cond_a, NULL,
                                    cond_i, cond_i_a, 0, Tm, fdsteps);

    am->randomseed = 0;
    am->seed = 123L;
    am->Nopt = 5;
    am->gtol = 1e-5;
    am->tol = 1e-5;
    am->n0 = 4;
    am->nofileout = 0;

    method_hybrid(am);

    am_delete(am);
}


/*

    Уравнение Бюргерса

*/

static void ai_newpoints_1d_nl(anns_instance_t *ai, double *bounds);
static void ai_update_1d_nonlinear(anns_instance_t *ai, double *x, int i, double tau);


ANNS_FUN(A1d_nl_i_f) {return 2.*M_PI*sin(M_PI*x[0]) / (2. + cos(M_PI*x[0]));}

ANNS_CND(A1d_nl_i) {return u[0] - v[0];}
ANNS_CND_G(A1d_nl_i_a) {return u_a;}

ANNS_CND(A1d_nl) {anns_instance_t *ai = instance; return u[0] + .5*(ai->step)*(u[0]*u_x[0][0] - u_xx[0][0]) + v[0];}
ANNS_CND(B1d_nl) {return u[0];}

ANNS_CND_G(A1d_nl_a) {anns_instance_t *ai = instance; return u_a + .5*(ai->step)*(u_a*u_x[0][0] + u[0]*u_ax[0] - u_axx[0]);}
ANNS_CND_G(B1d_nl_a) {return u_a;}


void ai_newpoints_1d_nl(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, A1d_nl_i_f, NULL);
}


void ai_update_1d_nonlinear(anns_instance_t *ai, double *x, int i, double tau) {
    int p;
    anns_cache_t *ac = ai->cache;

    // Сохраняем сеть i-1 шага в контейнер
    ai_put(ai, x, i-1);

    // Обновляем предысторию (v[0]) для конечно-разностной схемы (явно-неявная, Кранк-Никлсон)
    for (p = 0; p < ai->conds[0]->m; p++) {
        ai->nets[0]->ann_val[2](x, 0, p, 0, 0, ai);
        ac->v[0][p][0] = -ac->u[0] + .5*(ai->step)*(ac->u[0]*ac->u_x[0][0] - ac->u_xx[0][0]);
    }
}


void evolutionary_1d_nl(double Tm, int fdsteps) {
    int data[] = {
        PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 3, 0,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 1, 2, 1, fdsteps,
        20, 10,  1, /* Уравнение переноса */
        2, 2,   0, /* Границы */
        NORMALIZED | RADIAL | GAUSSIAN, 4, 1, 0,  2, 0};

    double bounds[] = {
        GP_GRID_1D, 20, 0, 1,
        GP_POINTS, 0, 1,
    };

    fbvp_t cond[] = {A1d_nl, B1d_nl};
    gbvp_t cond_a[] = {A1d_nl_a, B1d_nl_a};

    fbvp_t cond_i[] = {A1d_nl_i, B1d_nl};
    gbvp_t cond_i_a[] = {A1d_nl_i_a, B1d_nl_a};

    anns_methoddata_t *am = am_new3(ai_init, ai_newpoints_1d_nl, ai_initann_1d, NULL, ai_update_1d_nonlinear,
                                    data, bounds,
                                    cond, cond_a, NULL,
                                    cond_i, cond_i_a, 0, Tm, fdsteps);

    am->randomseed = 0;
    am->seed = 1231L;
    am->Nopt = 5;
    am->gtol = 1e-6;
    am->tol = 1e-5;
    am->n0 = 4;
    am->nofileout = 0;

    method_hybrid(am);

    am_delete(am);
}

#undef C_x
