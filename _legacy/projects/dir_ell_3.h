#ifndef ANNS_TASK
#define ANNS_TASK

#define NP_I 50
#define NP_B 20

int data[] = {
    PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 3, 0,

    FUNCTIOAL_SQUARE, 3, 7, 1, 0,

    NP_I, 1,
    NP_B, 0,
    NP_B, 0,
    NP_B, 0,
    NP_B, 0,
    NP_B, 0,
    NP_B, 0,

    NORMALIZED | RADIAL | GAUSSIAN, 64, 3, 0,   2, 0, 0, 0, 0, 0, 0,
};


double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID, 0.0, M_PI,  0.0, M_PI, 0.0, M_PI,// A

    GP_UNIFORM_SOLID_CUBOID, 0.0, M_PI,  0.0, M_PI,  0.0, 0.0,   // x, y, z=0
    GP_UNIFORM_SOLID_CUBOID, 0.0, M_PI,  0.0, M_PI,  M_PI, M_PI, // x, y, z=PI
    GP_UNIFORM_SOLID_CUBOID, 0.0, M_PI,  0.0, 0.0,   0.0, M_PI,  // x, y=0, z
    GP_UNIFORM_SOLID_CUBOID, 0.0, M_PI,  M_PI, M_PI, 0.0, M_PI,  // x, y=PI, z
    GP_UNIFORM_SOLID_CUBOID, 0.0, 0.0,   0.0, M_PI,  0.0, M_PI,  // x=0, y, z
    GP_UNIFORM_SOLID_CUBOID, M_PI, M_PI, 0.0, M_PI,  0.0, M_PI,  // x=PI, y, z
};


ANNS_CND(A) {return u_xx[0][0] + u_xx[0][1] + u_xx[0][2] - v[0];}
ANNS_CND(B) {return u[0];}

ANNS_CND_G(A_a) {return u_axx[0] + u_axx[1] + u_axx[2];}
ANNS_CND_G(B_a) {return u_a;}

ANNS_FUN(A_f) {return sin(x[0])*sin(x[1])*sin(x[2]);}

fbvp_t cond[] = {A, B, B, B, B, B, B};
gbvp_t cond_a[] = {A_a, B_a, B_a, B_a, B_a, B_a, B_a};
double delta[] = {1./NP_I, 1./NP_B, 1./NP_B, 1./NP_B, 1./NP_B, 1./NP_B, 1./NP_B};

void ai_init(anns_instance_t *ai) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_reload(anns_instance_t *ai) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, A_f, NULL);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (t = 0; t < ai->nnets; t++) {
        for (i = ai->I[t]; i < ai->I[t+1]; ) {
            x[i++] = uniform(-4., 4.);         // Веса
            x[i++] = M_PI/(2*pow(16, 1./3));
//            x[i++] = uniform(0.5, 4.);        // Ширины
            for (j = 0; j < ai->nd[t]; j++) {
                x[i++] = uniform(-M_PI/2,M_PI+M_PI/2);       // Центры
            }
            printf("%g %g %g %g %g\n", x[i-5], x[i-4], x[i-3], x[i-2], x[i-1]);
        }
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    double s, tmp, point[3];
    int i, p, M = 10000;

    s = 0;
    for (p = 0; p < M; p++) {
        point[0] = uniform(0, M_PI);
        point[1] = uniform(0, M_PI);
        point[2] = uniform(0, M_PI);

        tmp = 0;
        tmp += anns_eval(point, x, ai);
        tmp += sin(point[0])*sin(point[1])*sin(point[2]) / 3.0;

        s += tmp * tmp;
    }

    return sqrt(s / (M-1) );
}


#endif // ANNS_TASK
