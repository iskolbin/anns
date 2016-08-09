#ifndef ANNS_TASK
#define ANNS_TASK

#define C_alpha 1.0
#define C_gamma 0.1
#define C_Nnu   1.0

int data[] = {
    PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 8, 1,

    FUNCTIONAL_SQUARE, 1, 3, 1, 0,

    20, 0,
    1, 0,
    1, 0,

    NORMALIZED | RADIAL | GAUSSIAN, 4, 1, 0,  2, 0, 0,
};


double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID,  0.01, 1.0, // A
    GP_POINTS, 0.0,
    GP_POINTS, 1.0,
};


ANNS_CND(A) {return u_xx[0][0] + 2*u_x[0][0]/x[0] + C_alpha*exp(u[0]/(1+u[0]/C_gamma));}
ANNS_CND(B) {return u_x[0][0];}
ANNS_CND(C) {return C_Nnu*u[0] + u_x[0][0];}

ANNS_CND_G(A_a) {
    double K = 1./(1.+u[0]/C_gamma);
    return u_axx[0] + 2*u_ax[0]/x[0] + C_alpha*exp(K*u[0])*u_a*K*(1-u[0]*K/C_gamma);}

ANNS_CND_G(B_a) {return u_ax[0];}
ANNS_CND_G(C_a) {return C_Nnu*u_a + u_ax[0];}

fbvp2_t cond[] = {A, B, C};
gbvp2_t cond_a[] = {
    A_a, B_a, C_a,
};

double delta[] = {1, 1, 1};

void ai_init(anns_instance_t *ai) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_reload(anns_instance_t *ai) {
    ai_genpoints(ai, bounds);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (t = 0; t < ai->Nt; t++) {
        for (i = ai->I[t]; i < ai->I[t+1]; ) {
            x[i++] = uniform(-2., 2.);         // Веса
            x[i++] = uniform(0.25, 1.);        // Ширины
            x[i++] = uniform(-0.25, 1.25);       // Центры
        }
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    return -1;
}

#endif // ANNS_TASK
