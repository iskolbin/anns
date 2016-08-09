#ifndef ANNS_TASK
#define ANNS_TASK

int data[] = {
    PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 1, 1,

    FUNCTIONAL_SQUARE, 2, 5, 2, 0,

    20, 1,
    10, 0,
    10, 0,
    10, 0,
    10, 0,

    CLASSIC | RADIAL | GAUSSIAN, 8, 2, 0,  2, 0, 0, 0, 0,
    CLASSIC | RADIAL | GAUSSIAN, 8, 2, 0,  2, 0, 0, 0, 0,
};


double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID, 0.0, M_PI,  0.0, M_PI, // A

    GP_UNIFORM_SOLID_CUBOID, 0.0, M_PI,  M_PI, M_PI, // Btop
    GP_UNIFORM_SOLID_CUBOID, 0.0, 0.0,   0.0, M_PI,  // Bleft
    GP_UNIFORM_SOLID_CUBOID, M_PI, M_PI, 0.0, M_PI,  // Bright
    GP_UNIFORM_SOLID_CUBOID, 0.0, M_PI,  0.0, 0.0    // Bbottom
};


ANNS_CND(A) {return u_xx[0][0] + u_xx[0][1] +  u_xx[1][0] + u_xx[1][1] - v[0];}
ANNS_CND(B) {return u[0] + u[1];}

ANNS_CND_G(A_a) {return u_axx[0] + u_axx[1];}
ANNS_CND_G(B_a) {return u_a;}

ANNS_FUN(A_f) {return sin(x[0])*sin(x[1]);}

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
            x[i++] = uniform(-4., 4.);         // ����
            x[i++] = uniform(0.5, 4.);        // ������
            x[i++] = uniform(-M_PI/2,M_PI+M_PI/2);       // ������
            x[i++] = uniform(-M_PI/2,M_PI+M_PI/2);
        }
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    double s, tmp, lambda[32], point[2];
    int i, p, Np = 10000;

    s = 0;
    for (p = 0; p < Np; p++) {
        point[0] = uniform(0, M_PI);
        point[1] = uniform(0, M_PI);

        tmp = 0;
        tmp += anns_eval(point, x, ai);
        tmp += 0.5*sin(point[0])*sin(point[1]);

        s += tmp * tmp;
    }

    return sqrt(s / (Np-1) );
}

#endif // ANNS_TASK
