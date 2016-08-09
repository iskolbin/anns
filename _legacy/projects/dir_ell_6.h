#ifndef ANNS_TASK
#define ANNS_TASK

#define C_k 0.5
#define C_N 0.1

int data[] = {
    PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 6, 1,

    FUNCTIONAL_SQUARE, 1, 3, 1, 0,

    20, 2,
    1, 0,
    1, 0,

    NORMALIZED | RADIAL | GAUSSIAN, 4, 1, 0,  2, 0, 0,
};


double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID,  C_k, 1.0, // A
    GP_POINTS, C_k,
    GP_POINTS, 1.0,
};


ANNS_CND(A) {return v[0]*u_x[0][0] + u_xx[0][0] + v[1];}
ANNS_CND(B) {return u[0];}
ANNS_CND(C) {return u[0] - 1;}

ANNS_CND_G(A_a) {return v[0]*u_ax[0] + u_axx[0];}
ANNS_CND_G(B_a) {return u_a;}
ANNS_CND_G(C_a) {return u_a;}

ANNS_FUN(A_f1) {return 1/x[0];}
ANNS_FUN(A_f2) {return 4*C_N/(x[0]*x[0]*x[0]*x[0]);}

fbvp2_t cond[] = {A, B, C};
gbvp2_t cond_a[] = {
    A_a, B_a, C_a,
};

double delta[] = {1e-4, 1, 1};

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

    for (t = 0; t < ai->Nt; t++) {
        for (i = ai->I[t]; i < ai->I[t+1]; ) {
            x[i++] = uniform(-2., 2.);         // Веса
            x[i++] = uniform(0.25, 1.);        // Ширины
            x[i++] = uniform(-0.25 + C_k, 1.25);       // Центры
        }
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    double s, tmp, point[1];
    int i, p, Np = 10000;

    s = 0;
    for (p = 0; p < Np; p++) {
        point[0] = uniform(C_k, 1.0);

        tmp = 0;
        tmp += anns_eval(point, x, ai);
        tmp -= (C_N + 1 - C_N/(point[0]*point[0])) - (C_N + 1 - C_N/(C_k*C_k))*log(point[0])/log(C_k);;

        s += tmp * tmp;
    }

    return sqrt(s / (Np-1) );
}

#endif // ANNS_TASK
