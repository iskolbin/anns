#ifndef ANNS_TASK
#define ANNS_TASK

#define C_k 27.79

int data[] = {
    PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 11, 0,

    FUNCTIONAL_SQUARE, 1, 3, 1, 0,

    20, 0,
    1, 0,
    1, 0,

    NORMALIZED | RADIAL | GAUSSIAN, 8, 1, 0,  2, 0, 0
};


double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID,  0.0, 1.0,
    GP_UNIFORM_SOLID_CUBOID,  0.0, 0.0,
    GP_UNIFORM_SOLID_CUBOID,  1.0, 1.0,
};

ANNS_CND(A) {return u_xx[0][0] - C_k*C_k*u[0];}
ANNS_CND(B) {return u[0]-1;}

ANNS_CND_G(A_a) {return u_axx[0] - C_k*C_k*u_a;}
ANNS_CND_G(B_a) {return u_a;}

fbvp2_t cond[] = {A, B, B};
gbvp2_t cond_a[] = {A_a, B_a, B_a};
double delta[] = {1e-6, 1, 1};

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
            x[i++] = uniform(0.25, 1.25);        // Ширины
            x[i++] = uniform(-0.25, 1.25);       // Центры
        }
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    double s, tmp, point[1];
    int i, p, Np = 10001;
    double dx = 1./(Np-1);
    double ek = exp(C_k), ek_ = exp(-C_k), ekx, ekx_;

    s = 0;
    for (p = 0; p < Np; p++) {
        point[0] =  dx * p;

        tmp = 0;
        tmp += anns_eval(point, x, ai);

        ekx = exp(C_k * point[0]);
        ekx_ = exp(-C_k * point[0]);

        tmp -= ((ek - 1) * ekx_ + (1 - ek_) * ekx ) / (ek - ek_);

        s += tmp * tmp;
    }

    return sqrt(s / (Np-1) );;
}

#endif // ANNS_TASK

