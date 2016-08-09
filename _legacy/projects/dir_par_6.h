#ifndef ANNS_TASK
#define ANNS_TASK

int data[]={
    PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 6, 0,

    FUNCTIONAL_SQUARE, 3, 6, 1, 0,

    250, 0, /* A_omega */

    25, 0, /* A_initial */

    25, 0, /* A_Gamma1 */
    25, 0, /* A_Gamma2 */
    25, 0, /* A_gamma_star */
    25, 1, /* A_gamma */

    NORMALIZED | RADIAL | GAUSSIAN, 32, 3, 0,  2, 0, 1, 1, 1, 1,
};

#define C_l1 1.0
#define C_l2 0.5
#define C_T  0.25

double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID, 0, C_l1,    0, C_l2,    0, C_T,

    GP_UNIFORM_SOLID_CUBOID, 0, C_l1,    0, C_l2,    0, 0,

    GP_UNIFORM_SOLID_CUBOID, 0, 0,       0, C_l2,    0, C_T,
    GP_UNIFORM_SOLID_CUBOID, C_l1, C_l1, 0, C_l2,    0, C_T,
    GP_UNIFORM_SOLID_CUBOID, 0, C_l1,    0, 0,       0, C_T,
    GP_UNIFORM_SOLID_CUBOID, 0, C_l1,    C_l2, C_l2, 0, C_T,
};

ANNS_CND(A_omega) {return u_xx[0][0] + u_xx[0][1] - u_x[0][2];}
ANNS_CND(A_Gamma) {return u_x[0][0];}
ANNS_CND(A_gamma) {return u_x[0][1] - v[0];}
ANNS_CND(A_gamma_star) {return u_x[0][1];}
ANNS_CND(A_initial) {return u[0];}


ANNS_CND_G(A_omega_a) {return u_axx[0] + u_axx[1] - u_ax[2];}
ANNS_CND_G(A_Gamma_a) {return u_ax[0];}
ANNS_CND_G(A_gamma_a) {return u_ax[1];}
ANNS_CND_G(A_initial_a) {return u_a;}


ANNS_FUN(F_mu_1) {return x[0]*x[2]/C_T;}


fbvp2_t cond[] = {A_omega, A_initial, A_Gamma, A_Gamma, A_gamma_star, A_gamma};
gbvp2_t cond_a[] = {A_omega_a, A_initial_a, A_Gamma_a, A_Gamma_a, A_gamma_a, A_gamma_a};

void ai_init(anns_instance_t *ai) {
    ai_conditions(ai, cond, cond_a, NULL);
    ai_loadfuncs(ai);
}

void ai_reload(anns_instance_t *ai) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 5, 0, F_mu_1, NULL);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (t = 0; t < ai->Nt; t++) {
        for (i = ai->I[t]; i < ai->I[t+1]; ) {
            x[i++] = uniform(-4., 4.);         // Веса
            x[i++] = uniform(0.1, 2.);        // Ширины
            x[i++] = uniform(-0.25, C_l1 + 0.25);
            x[i++] = uniform(-0.25, C_l2 + 0.25);      // Центры
            x[i++] = uniform(-0.25, C_T + 0.25);      // Центры
        }
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
//    double s, tmp, point[2];
//    int i, p, Np = 10000;
//
//    s = 0;
//    for (p = 0; p < Np; p++) {
//        point[0] = uniform(0, C_l1);
//        point[1] = uniform(0, C_l2);
//
//        tmp = 0;
//        tmp += anns_eval(point, x, ai);
//        tmp -= exp(-x[1])*sin(x[0]);
//
//        s += tmp * tmp;
//    }
//
//    return sqrt(s / (Np-1) );
    return -1;
}

#endif // ANNS_TASK
