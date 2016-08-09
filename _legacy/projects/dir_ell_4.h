#ifndef ANNS_TASK
#define ANNS_TASK

int data[] = {
    PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 4, 1,

    FUNCTIONAL_SQUARE, 2, 2, 1, 0,

    100, 1,
    10, 0,

    NORMALIZED | RADIAL | GAUSSIAN, 16, 2, 0,  2, 0, 0, 0, 0,
};


double bounds[] = {
    GP_UNIFORM_SOLID_SPHERE,  0.0, 1.0, 0.0, 0.0, // A
    GP_UNIFORM_SOLID_SPHERE,  1.0, 1.0, 0.0, 0.0,
};


ANNS_CND(A) {return u_xx[0][0] + u_xx[0][1] - v[0];}
ANNS_CND(B) {return u[0];}

ANNS_CND_G(A_a) {return u_axx[0] + u_axx[1];}
ANNS_CND_G(B_a) {return u_a;}

ANNS_FUN(A_f) {
    double xc = x[0] - 0.4, yc = x[1] - 0;
    return (xc*xc+yc*yc <= 0.4) ? 1 : 0;}

fbvp2_t cond[] = {A, B, B, B, B};
gbvp2_t cond_a[] = {
    A_a, B_a, B_a, B_a, B_a,
    A_a, B_a, B_a, B_a, B_a,
};

void ai_init(anns_instance_t *ai) {
    ai_conditions(ai, cond, cond_a, NULL);
    ai_loadfuncs(ai);
}

void ai_reload(anns_instance_t *ai) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, A_f, NULL);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (t = 0; t < ai->Nt; t++) {
        for (i = ai->I[t]; i < ai->I[t+1]; ) {
            x[i++] = uniform(-4., 4.);         // Веса
            x[i++] = uniform(0.25, 1.);        // Ширины
            x[i++] = uniform(-0.25, 1.25);       // Центры
            x[i++] = uniform(-0.25, 1.25);
        }
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    return -1;
}

#endif // ANNS_TASK
