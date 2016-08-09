#ifndef ANNS_TASK
#define ANNS_TASK

int data[] = {
    PJ_DIRECT, PJ_APPROXIMATION, 2, 0,

    FUNCTIONAL_SQUARE, 1, 3, 2, 0,

    11, 1,
    11, 1,
    1, 0,

    NORMALIZED | RADIAL | GAUSSIAN, 4, 1, 0,  0, -1, 0,
    NORMALIZED | RADIAL | GAUSSIAN, 4, 1, 0,  -1, 0, 0,
};


double bounds[] = {
    GP_GRID_1D, 11, 0.0, 0.5,
    GP_GRID_1D, 11, 0.5, 1.0,
    GP_POINTS, 0.5,
};

ANNS_CND(A) {return u[0] - v[0];}
ANNS_CND(B) {return u[1] - v[0];}
ANNS_CND(C) {return u[1] - u[0];}

ANNS_CND_G(A_a) {return u_a;}
ANNS_CND_G(B_a) {return u_a;}
ANNS_CND_G(C_a_1) {return u_a;}
ANNS_CND_G(C_a_2) {return -u_a;}

ANNS_FUN(A_f) {return (x[0]<0.5) ? 2*x[0] : 2 - 2*x[0];}

fbvp2_t cond[] = {A, B, C};
gbvp2_t cond_a[] = {
    A_a, NULL, C_a_1,
    NULL, B_a, C_a_2};

void ai_init(anns_instance_t *ai) {
    ai_conditions(ai, cond, cond_a, NULL);
    ai_loadfuncs(ai);
}

void ai_reload(anns_instance_t *ai) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, A_f, NULL);
    ai_eval_v(ai, 1, 0, A_f, NULL);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (t = 0; t < ai->Nt; t++) {
        for (i = ai->I[t]; i < ai->I[t+1]; ) {
            x[i++] = uniform(-2., 2.);         // Веса
            x[i++] = uniform(0.25, 1.25);        // Ширины
            x[i++] = uniform(-0.25, 1.25);       // Центры
        }
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    double s, tmp, point[1];
    int i, p, Np = 10001;
    double dx = 1./(Np-1);

    s = 0;
    for (p = 0; p < Np; p++) {
        point[0] =  dx * p;

        tmp = (point[0] < 0.5) ? anns_evalsingle(0, point, x, ai) : anns_evalsingle(1, point, x, ai);
        tmp -= A_f(point, 1, NULL);

        s += tmp * tmp;
    }

    return sqrt(s / (Np-1) );;
}

#endif // ANNS_TASK
