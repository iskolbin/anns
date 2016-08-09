#ifndef ANNS_TASK
#define ANNS_TASK

int data[] = {
    PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 2, 0,

    FUNCTIOAL_SQUARE, 3, 7, 1, 0,

    40, 1,

    10, 1,
    10, 1,
    10, 1,
    10, 1,
    10, 1,
    10, 1,

    NORMALIZED | RADIAL | GAUSSIAN, 4, 3, 0,  2, 0, 0, 0, 0, 0, 0
};


double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID, 0.0, 1.0,  0.0, 1.0,  0.0, 1.0, // A

    GP_UNIFORM_SOLID_CUBOID, 0.0, 1.0,  0.0, 1.0,  0.0, 0.0, // B z=0
    GP_UNIFORM_SOLID_CUBOID, 0.0, 1.0,  0.0, 1.0,  1.0, 1.0, // B z=1
    GP_UNIFORM_SOLID_CUBOID, 0.0, 1.0,  0.0, 0.0,  0.0, 1.0, // B y=0
    GP_UNIFORM_SOLID_CUBOID, 0.0, 1.0,  1.0, 1.0,  0.0, 1.0, // B y=1
    GP_UNIFORM_SOLID_CUBOID, 0.0, 0.0,  0.0, 1.0,  0.0, 1.0, // B x=0
    GP_UNIFORM_SOLID_CUBOID, 1.0, 1.0,  0.0, 1.0,  0.0, 1.0, // B x=1
};


ANNS_CND(A) {return u_xx[0][0] + u_xx[0][1] + u_xx[0][2] - sin(u[0]*u[0]) - v[0];}
ANNS_CND(B) {return u[0]- v[0];}

ANNS_CND_G(A_a) {return u_axx[0] + u_axx[1] + u_axx[2] - 2*u[0]*u_a*cos(u[0]*u[0]);}
ANNS_CND_G(B_a) {return u_a;}

ANNS_FUN(A_f) {double v = sin(x[0] + x[1] + x[2]); return -3*v - sin(v*v);}
ANNS_FUN(B_f) {return sin(x[0] + x[1] + x[2]);}

fbvp_t cond[] = {A, B, B, B, B, B, B};
gbvp_t cond_a[] = {A_a, B_a, B_a, B_a, B_a, B_a, B_a};
double delta[] = {1./40, .1, .1, .1, .1, .1, .1};

void ai_init(anns_instance_t *ai) {
    ai_conditions(ai, cond, cond_a, NULL);
    ai_loadfuncs(ai);
}

void ai_reload(anns_instance_t *ai) {
    int i;
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, A_f, NULL);
    for (i = 1; i < 7; i++) {
        ai_eval_v(ai, i, 0, B_f, NULL);
    }
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (t = 0; t < ai->nnets; t++) {
        for (i = ai->I[t]; i < ai->I[t+1]; ) {
            x[i++] = uniform(-4., 4.);   // Веса
            x[i++] = uniform(0.5, 2.);   // Ширины
            x[i++] = uniform(-0.5, 1.5); // Центры
            x[i++] = uniform(-0.5, 1.5);
            x[i++] = uniform(-0.5, 1.5);
        }
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    double s, tmp, point[3];
    int i, p, M = 10000;

    s = 0;
    for (p = 0; p < M; p++) {
        point[0] = uniform(0.0, 1.0);
        point[1] = uniform(0.0, 1.0);
        point[2] = uniform(0.0, 1.0);

        tmp = 0;
        tmp += anns_eval(point, x, ai);
        tmp += -sin(point[0] + point[1] + point[2]);

        s += tmp * tmp;
    }

    return sqrt(s / (M-1) );
}

#endif // ANNS_TASK
