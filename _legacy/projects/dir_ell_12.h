#ifndef ANNS_TASK
#define ANNS_TASK

int data[] = {
    PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 12, 0,

    FUNCTIONAL_SQUARE, 1, 3, 1, 0,

    50, 1,

    1, 0,
    1, 0,

    NORMALIZED | RADIAL | GAUSSIAN, 4, 1, 0,  2, 0, 0
};

#define C_k 0.05

double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID, 0.0, 1.0, // A

    GP_POINTS, 0.0,
    GP_POINTS, 1.0,
};


ANNS_CND(A) {return C_k * u_xx[0][0] + v[0];}
ANNS_CND(B) {return u[0];}

ANNS_CND_G(A_a) {return C_k * u_axx[0];}
ANNS_CND_G(B_a) {return u_a;}

ANNS_FUN(A_f) {return (x[0] <= 0.5) ? x[0] : 0.5;}

fbvp2_t cond[] = {A, B, B};
gbvp2_t cond_a[] = {A_a, B_a, B_a};

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
            x[i++] = uniform(-4., 4.);   // Веса
            x[i++] = uniform(0.5, 2.);   // Ширины
            x[i++] = uniform(-0.5, 1.5); // Центры
        }
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    return -1;
}

#endif // ANNS_TASK
