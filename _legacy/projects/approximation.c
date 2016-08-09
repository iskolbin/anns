#include "approximation.h"

#define C_k 27.79

ANNS_CND(A_appr) {return u[0] - v[0];}
ANNS_CND_G(A_appr_a) {return u_a;}
ANNS_FUN(A_appr_bl_f) {
    double ek = exp(C_k), ek_ = exp(-C_k), ekx = exp(C_k * x[0]), ekx_ = exp(-C_k * x[0]);
    return ((ek - 1) * ekx_ + (1 - ek_) * ekx ) / (ek - ek_);
}

static void ai_init(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
static void ai_initann(anns_instance_t *, double *);
static void ai_newpoints_bl(anns_instance_t *, double *bounds);
static double ai_mse_bl(double *, int, void *);

void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_newpoints_bl(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, A_appr_bl_f, NULL);
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


double ai_mse_bl(double *x, int n, void *instance) {
    anns_instance_t *ai = instance;
    double s=0, tmp, point[1];
    int i, p, Np = 101;
    double dx = 1./(Np-1);
    double ek = exp(C_k), ek_ = exp(-C_k), ekx, ekx_;

    s = 0;
    for (p = 0; p < Np; p++) {
        point[0] =  dx * p;

        tmp = 0;
        tmp += anns_eval(point, x, ai);

        ekx = exp(C_k * point[0]);
        ekx_ = exp(-C_k * point[0]);

        tmp -= ((ek - 1) * ekx_ + (1 - ek_) * ekx ) / (ek - ek_);

        s += tmp * tmp;
    }

    return sqrt(s / (Np-1) );;
}



void approx_bl(void) {
    int np = 101;
    int data[] = {
        PJ_DIRECT, PJ_APPROXIMATION, 4, 3,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 1, 1, 1, 0,
        np, np, 1, /* Уравнение переноса */
        NORMALIZED | RADIAL | GAUSSIAN, 6, 1, 0,  0};

    double bounds[] = {
        GP_GRID_1D, np, 0, 1};

    fbvp_t cond[] = {A_appr};
    gbvp_t cond_a[] = {A_appr_a};

    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints_bl, ai_initann, ai_mse_bl, NULL, data, bounds, cond, cond_a, NULL);

    am->randomseed = 0;
//    am->seed = 12311L;
    am->seed = 12311111L;
    am->Nopt = 5;
    am->gtol = 1e-6;
    am->tol = 1e-5;
//    am->nofileout = 0;

    method_restarts(am);

    as_delete(am);
}




#define C_spec_xmin -2.0
#define C_spec_xmax 3.0
#define C_spec_npts 201

ANNS_FUN(A_appr_spec_f) {
    return x[0] < 0.25 ? 0 :
        x[0] < 0.5 ? x[0] - 0.25 :
            x[0] < 0.75 ?  x[0] - 0.25 -32*(x[0]-0.625)*(x[0]-0.625) + 0.75 :
                x[0] < 1 ? x[0] - 0.25:
                    0.75;
}

double ai_mse_approx_spec(double *x, int n, void *instance) {
    double s=0, tmp, point[1];
    int p, Np = C_spec_npts;
    double dx = (C_spec_xmax - C_spec_xmin)/(Np-1);

    s = 0;
    for (p = 0; p < Np; p++) {
        point[0] =  C_spec_xmin + dx * p;
        tmp = anns_eval(point, x, instance) - A_appr_spec_f(point,1,NULL);
        s += tmp * tmp;
    }

    return sqrt(s / (Np-1) );;
}

static void ai_newpoints_spec(anns_instance_t *, double *bounds);
static double ai_mse_spec(anns_instance_t *, double *);
static void ai_initann_spec(anns_instance_t *, double *);

void ai_newpoints_spec(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, A_appr_spec_f, NULL);
}

void ai_initann_spec(anns_instance_t *ai, double *x) {
    int i = 0, j, t, n, s;
    double b = M_PI/sqrt(ai->nets[0]->nc), d;

    for (t = 0; t < ai->nnets; t++) {
        for (i = ai->nets[t]->I; i < ai->nets[t]->lenc;) {
            x[i++] = uniform(-2., 2.);         // Веса
            x[i++] = b;
            for (j = 0; j < ai->nets[t]->dim; j++) {
                x[i++] = uniform(-0.05 + C_spec_xmin, C_spec_xmax + 0.05);       // Центры
            }
        }

        an_setwidth_nearest(ai->nets[t], 1.5);
    }
}

// Линейная зависимость с импульсом
void approx_spec(void) {
    int np = C_spec_npts;
    int data[] = {
        PJ_DIRECT, PJ_APPROXIMATION, 1, 2,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 1, 1, 1, 0,
        np, np, 1, /* Уравнение переноса */
        CLASSIC | RADIAL | GAUSSIAN, 8, 1, 0,  0};

    double bounds[] = {
        GP_GRID_1D, np, C_spec_xmin, C_spec_xmax};

    fbvp_t cond[] = {A_appr};
    gbvp_t cond_a[] = {A_appr_a};

    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints_spec, ai_initann, ai_mse_approx_spec, NULL, data, bounds, cond, cond_a, NULL);

    am->randomseed = 0;
    am->seed = 123111L;
//    am->seed = 12311111L;
    am->Nopt = 5;
    am->gtol = 1e-4;
    am->tol = 1e-4;
    am->nofileout = 0;

    method_restarts(am);

    as_delete(am);
}





#define C_spec2_xmin 0.0
#define C_spec2_xmax 1.0
#define C_spec2_npts 101

ANNS_FUN(A_appr_spec2_f) {
    return x[0] < 1./3. ? 0 :
        x[0] < 2./3. ? 3*x[0] - 1:
            1.;
}



static void ai_newpoints_spec2(anns_instance_t *, double *bounds);
static double ai_mse_spec2(anns_instance_t *, double *);
static void ai_initann_spec2(anns_instance_t *, double *);
double ai_mse_approx_spec2(double *x, int n, void *instance) {
    double s=0, tmp, point[1];
    int p, Np = C_spec2_npts;
    double dx = (C_spec2_xmax - C_spec2_xmin)/(Np-1);

    s = 0;
    for (p = 0; p < Np; p++) {
        point[0] =  C_spec2_xmin + dx * p;
        tmp = anns_eval(point, x, instance) - A_appr_spec2_f(point,1,NULL);
        s += tmp * tmp;
    }

    return sqrt(s / (Np-1) );;
}
void ai_newpoints_spec2(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, A_appr_spec2_f, NULL);
}

void ai_initann_spec2(anns_instance_t *ai, double *x) {
    int i = 0, j, t, n, s;
    double b = M_PI/sqrt(ai->nets[0]->nc), d;

    for (t = 0; t < ai->nnets; t++) {
        for (i = ai->nets[t]->I; i < ai->nets[t]->lenc;) {
            x[i++] = uniform(-2., 2.);         // Веса
            x[i++] = b;
            for (j = 0; j < ai->nets[t]->dim; j++) {
                x[i++] = uniform(-0.05 + C_spec2_xmin, C_spec2_xmax + 0.05);       // Центры
            }
        }

        an_setwidth_nearest(ai->nets[t], 1.5);
    }
}

// Линейная зависимость
void approx_spec2(void) {
    int np = C_spec2_npts;
    int data[] = {
        PJ_DIRECT, PJ_APPROXIMATION, 2, 2,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 1, 1, 1, 0,
        np, np, 1, /* Уравнение переноса */
        NORMALIZED | RADIAL | GAUSSIAN, 8, 1, 0,  0};

    double bounds[] = {
        GP_GRID_1D, np, C_spec2_xmin, C_spec2_xmax};

    fbvp_t cond[] = {A_appr};
    gbvp_t cond_a[] = {A_appr_a};

    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints_spec2, ai_initann, ai_mse_approx_spec2, NULL, data, bounds, cond, cond_a, NULL);

    am->randomseed = 0;
    am->seed = 123111L;
//    am->seed = 12311111L;
    am->Nopt = 5;
    am->gtol = 1e-5;
    am->tol = 1e-5;
    am->nofileout = 0;

    method_restarts(am);

    as_delete(am);
}



#define C_spec3_xmin 0.0
#define C_spec3_xmax 1.0
#define C_spec3_npts 501

ANNS_FUN(A_appr_spec3_f) {
    return x[0] < 0.5 ? 0 : 1;
}


static void ai_newpoints_spec3(anns_instance_t *, double *bounds);
static double ai_mse_spec3(anns_instance_t *, double *);
static void ai_initann_spec3(anns_instance_t *, double *);

double ai_mse_approx_spec3(double *x, int n, void *instance) {
    double s=0, tmp, point[1];
    int p, Np = C_spec3_npts;
    double dx = (C_spec3_xmax - C_spec3_xmin)/(Np-1);

    s = 0;
    for (p = 0; p < Np; p++) {
        point[0] =  C_spec3_xmin + dx * p;
        tmp = anns_eval(point, x, instance) - A_appr_spec3_f(point,1,NULL);
        s += tmp * tmp;
    }

    return sqrt(s / (Np-1) );;
}


void ai_newpoints_spec3(anns_instance_t *ai, double *bounds) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, A_appr_spec3_f, NULL);
}

void ai_initann_spec3(anns_instance_t *ai, double *x) {
    int i = 0, j, t, n, s;
    double b = M_PI/sqrt(ai->nets[0]->nc), d;

    for (t = 0; t < ai->nnets; t++) {
        for (i = ai->nets[t]->I; i < ai->nets[t]->lenc;) {
            x[i++] = uniform(-2., 2.);         // Веса
            x[i++] = b;
            for (j = 0; j < ai->nets[t]->dim; j++) {
                x[i++] = uniform(-0.05 + C_spec3_xmin, C_spec3_xmax + 0.05);       // Центры
            }
        }

        an_setwidth_nearest(ai->nets[t], 1.5);
    }
}

// Ступенька
void approx_spec3(void) {
    int np = C_spec3_npts;
    int data[] = {
        PJ_DIRECT, PJ_APPROXIMATION, 3, 2,
        FUNCTIONAL_SQUARE, DELTA_INV_ALL, 1, 1, 1, 0,
        np, np, 1, /* Уравнение переноса */
        CLASSIC | RADIAL | GAUSSIAN, 8, 1, 0,  0};

    double bounds[] = {
        GP_GRID_1D, np, C_spec3_xmin, C_spec3_xmax};

    fbvp_t cond[] = {A_appr};
    gbvp_t cond_a[] = {A_appr_a};

    anns_methoddata_t *am = am_new2(ai_init, ai_newpoints_spec3, ai_initann, ai_mse_approx_spec3, NULL, data, bounds, cond, cond_a, NULL);

    am->randomseed = 0;
    am->seed = 1231111L;
//    am->seed = 12311111L;
    am->Nopt = 5;
    am->gtol = 1e-5;
    am->tol = 1e-5;
    am->nofileout = 0;

    method_restarts(am);

    as_delete(am);
}
#undef C_k
