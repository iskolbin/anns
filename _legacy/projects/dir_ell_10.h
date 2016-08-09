#ifndef ANNS_TASK
#define ANNS_TASK

#define C_rho    62.
#define C_omega1 10.
#define C_omega2 30.
#define C_r1     5.
#define C_r2     6.
#define C_scale  1e-0

#define C_alpha  ((C_omega1*C_r1*C_r1 - C_omega2*C_r2*C_r2)/(C_r2*C_r2 - C_r1*C_r1))
#define C_beta   (C_r1*C_r1*C_r2*C_r2 * (C_omega1 - C_omega2) / (C_r2*C_r2 - C_r1*C_r1) )

int data[] = {
    PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 10, 0,

    FUNCTIOAL_SQUARE, 1, 4, 2, 0,

    20, 0,
    20, 0,
    1, 0,
    1, 0,

    NORMALIZED | RADIAL | GAUSSIAN, 4, 1, 0,  1, -1, -1, -1,
    NORMALIZED | RADIAL | GAUSSIAN, 4, 1, 0,  0, 2, 1, 1,
};


double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID,  C_r1, C_r2, //
    GP_UNIFORM_SOLID_CUBOID,  C_r1, C_r2,
    GP_POINTS, C_r1,
    GP_POINTS, C_r2,
};


ANNS_CND(A1) {return u[1] - u_x[0][0];}
ANNS_CND(A2) {
    double a1 = 3*u_x[1][0] + x[0]*u_xx[1][0] - u[1]/x[0];
    double a2 = (u[1] + x[0]*u_x[1][0]) * (u[1] + x[0]*u_x[1][0]) /(2*x[0]*u[1]);
    double c = 0.5 / sqrt(x[0]*u[1]);

    return sqrt(C_scale) * c * (a1 - a2);}

ANNS_CND(GLeft) {return C_scale * (u_x[1][0] - C_rho * C_r1 * C_omega1*C_omega1);}
ANNS_CND(GRight) {return C_scale * (u_x[1][0] - C_rho * C_r2 * C_omega2*C_omega2);}

ANNS_CND_G(A1_a_0) {return u_ax[0];}
ANNS_CND_G(A1_a_1) {return u_a;}

ANNS_CND_G(A2_a) {
    double a1 = x[0]*u_axx[0] - u_a/x[0] + 3*u_ax[0];

    double c = x[0]*u_x[1][0] + u[1];
    double d = 2*sqrt(x[0]*u[1]);

    double a2 = c*c*u_a / (2*x[0]*u[1]*u[1]);
    double a3 = c * (x[0]*u_ax[0] + u_a) / (x[0] * u[1]);

    double b1 = x[0]*u_xx[1][0] + 3*u_x[1][0] - u[1]/x[0];
    double b2 = c*c / (2*x[0]*u[1]);

    return sqrt(C_scale) * ( (a1 + a2 + a3) / d  +  x[0]*u_a * (b1 + b2) / (2*x[0]*u[1]*d) );
}


ANNS_CND_G(GLeft_a) {return u_ax[0];}
ANNS_CND_G(GRight_a) {return u_ax[0];}

fbvp_t cond[] = {A1, A2, GLeft, GRight};
gbvp_t cond_a[] = {
    A1_a_0, NULL, NULL, NULL,
    A1_a_1, A2_a, GLeft_a, GRight_a,
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
