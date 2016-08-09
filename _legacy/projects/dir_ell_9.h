#ifndef ANNS_TASK
#define ANNS_TASK

#define C_alpha  (M_PI/30)
//#define C_beta   0.10014
#define C_beta   1.0
#define C_rho    0.5
#define C_theta  0.05

int data[] = {
    PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 9, 0,

    FUNCTIOAL_SQUARE, 1, 3, 1, 0,

   20, 2,
    1, 0,
    1, 0,

    NORMALIZED | RADIAL | GAUSSIAN, 4, 1, 0,  2, 0, 1,
};


double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID,  0.0, 1.0, // A
    GP_POINTS, 0.0,
    GP_POINTS, 1.0,
};


ANNS_CND(A) {return u_xx[0][0] + v[0]*u_x[0][0] + v[1]*u[0]*u[0]*u[0]*u[0];}
ANNS_CND(B) {return u[0]-1;}
ANNS_CND(C) {return u_x[0][0];}

ANNS_CND_G(A_a) {return u_axx[0] + v[0]*u_ax[0] + 4*v[1]*u_a*u[0]*u[0]*u[0];}

ANNS_CND_G(B_a) {return u_a;}
ANNS_CND_G(C_a) {return u_ax[0];}

ANNS_FUN(A_f1) {
    double K1 = 1/(x[0] + C_rho);
    double K2 = -tan(C_alpha) / ((1-x[0])*tan(C_alpha) + C_theta);
    return K1 + K2;}
ANNS_FUN(A_f2) {return -C_beta / ((1-x[0])*tan(C_alpha) + C_theta);}

fbvp_t cond[] = {A, B, C};
gbvp_t cond_a[] = {
    A_a, B_a, C_a,
};

double delta[] = {1e-1, 1, 1};

void ai_init(anns_instance_t *ai) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_reload(anns_instance_t *ai) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, A_f1, NULL);
    ai_eval_v(ai, 0, 1, A_f2, NULL);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (t = 0; t < ai->nnets; t++) {
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
