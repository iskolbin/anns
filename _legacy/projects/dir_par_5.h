#ifndef ANNS_TASK
#define ANNS_TASK

int data[]={
    PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 5, 0,

    FUNCTIONAL_SQUARE, 2, 4, 1, 0,

    20, 0, /* A (Heat equation) */
    15, 1, /* Initial condition */
    15, 0, /* Bleft */
    15, 0, /* Bright */

    NORMALIZED | RADIAL | GAUSSIAN, 16, 2, 0,  2, 0, 0, 0, 0
};

#define X_max (M_PI+M_PI)
#define T_max (0.1*(4*M_PI*M_PI))

double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID, 0, X_max,    0, T_max,
    GP_UNIFORM_SOLID_CUBOID, 0, X_max,    0, 0,
    GP_UNIFORM_SOLID_CUBOID, 0, 0,        0, T_max,
    GP_UNIFORM_SOLID_CUBOID, X_max, X_max,0, T_max
};


ANNS_CND(A) {return u_xx[0][0] - u_x[0][1];}
ANNS_CND(I) {return u[0] - v[0];}
ANNS_CND(B) {return u[0];}


ANNS_CND_G(A_a) {return u_axx[0] - u_ax[1];}
ANNS_CND_G(I_a) {return u_a;}
ANNS_CND_G(B_a) {return u_a;}

ANNS_FUN(I_f) {return sin(x[1]);}

fbvp2_t cond[] = {A, I, B, B};
gbvp2_t cond_a[] = {A_a, I_a, B_a, B_a};

void ai_init(anns_instance_t *ai) {
    ai_conditions(ai, cond, cond_a, NULL);
    ai_loadfuncs(ai);
}

void ai_reload(anns_instance_t *ai) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 1, 0, I_f, NULL);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (t = 0; t < ai->Nt; t++) {
        for (i = ai->I[t]; i < ai->I[t+1]; ) {
            x[i++] = uniform(-4., 4.);         // Веса
            x[i++] = uniform(0.1, 2.);        // Ширины
            x[i++] = uniform(-0.25, X_max + 0.25);    // Центры
            x[i++] = uniform(-0.25, T_max + 0.25);      // Центры
        }
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    double s, tmp, point[2];
    int i, p, Np = 10000;

    s = 0;
    for (p = 0; p < Np; p++) {
        point[0] = uniform(0, X_max);
        point[1] = uniform(0, T_max);

        tmp = 0;
        tmp += anns_eval(point, x, ai);
        tmp -= exp(-x[1])*sin(x[0]);

        s += tmp * tmp;
    }

    return sqrt(s / (Np-1) );
}

#endif // ANNS_TASK
