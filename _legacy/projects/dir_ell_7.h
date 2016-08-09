#ifndef ANNS_TASK
#define ANNS_TASK

int data[] = {
    PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 7, 1,

    FUNCTIONAL_SQUARE, 1, 3, 1, 0,

    20, 2,
    1, 0,
    1, 0,

    NORMALIZED | RADIAL | GAUSSIAN, 8, 1, 0,  2, 1, 0,
};


double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID,  0.0, 1.0, // A
    GP_POINTS, 0.0,
    GP_POINTS, 1.0,
};


ANNS_CND(A) {return u[0]*u[0]*u[0]*u_xx[0][0] + sin(u[0]*u_x[0][0]) + 1;}
ANNS_CND(B) {return u_x[0][0];}
ANNS_CND(C) {return u[0];}

ANNS_CND_G(A_a) {return 3*u[0]*u[0]*u_a*u_xx[0][0] + u[0]*u[0]*u[0]*u_axx[0] + cos(u[0]*u_x[0][0])*(u_a*u_x[0][0] + u[0]*u_ax[0]);}
ANNS_CND_G(B_a) {return u_ax[0];}
ANNS_CND_G(C_a) {return u_a;}

fbvp2_t cond[] = {A, B, C};
gbvp2_t cond_a[] = {
    A_a, B_a, C_a,
};

double delta[] = {1e-2, 1, 1};

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
    double s, tmp, point[1];
    int i, p;

    point[0] = 0;
    tmp = anns_eval(point, x, ai) - 0.90168;

    return sqrt(tmp*tmp);
}

#endif // ANNS_TASK
