#ifndef ANNS_TASK
#define ANNS_TASK

#define PROBLEM_SENSORS 32
#define PROBLEM_ERROR   1
#define G_M 100

int data[] = {
    PJ_INVERSE_UNKNOWN_SOURCE, PJ_ELLIPTIC_PROBLEM, 2, 1,

    FUNCTIONAL_SQUARE, 2, 5, 2, 0,
//    FUNCTIONAL_SQUARE, 2, 4, 1, 0,
     100, 1,   // Уравнение переноса
//    100, 0,   // Уравнение переноса
    40, 0,   // x = 0
    40, 0,   // x = pi
    20, 0,   // y = 0
    PROBLEM_SENSORS, 1,   // y = 1
//    2, 0,

    NORMALIZED | RADIAL | GAUSSIAN, 24, 2, 0,   2, 0,  0,  0,  0,
//    -1,
//    NORMALIZED | RADIAL | GAUSSIAN, 16, 2, 0,   2,  0,  0,  0,  0,
    NORMALIZED | RADIAL | GAUSSIAN, 8,  1, 0,   0, -1, -1, -1, -1,
//    0,
//    NORMALIZED | RADIAL | GAUSSIAN, 16, 2, 0,   0, -1, -1, -1, -1,
};

double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID, 0.0, M_PI,   0.0, 2.5,  // Уравнение переноса
    GP_UNIFORM_SOLID_CUBOID, 0.0, 0.0,    0.0,2.5,  // x = 0
    GP_UNIFORM_SOLID_CUBOID, M_PI, M_PI,  0.0, 2.5,  // x = pi

    GP_UNIFORM_SOLID_CUBOID, 0.0, M_PI,   0.0, 0.0,        // y = 0
    GP_GRID_2D,              PROBLEM_SENSORS, 0.0, M_PI,   1, 1.0, 1.0,  // y = 1

//    GP_POINTS, 0.0, 0.0, M_PI, 0.0,
};

ANNS_CND(A)  {return u_xx[0][0] + u_xx[0][1] + u[1];}
ANNS_CND(B) {return u[0];}
ANNS_CND(Btop)  {return u[0] - v[0];}
//ANNS_CND(C) {return u[1];}

ANNS_CND_G(A_a_1)  {return u_axx[0] + u_axx[1];}
ANNS_CND_G(A_a_2)  {return u_a;}
ANNS_CND_G(B_a) {return u_a;}
//ANNS_CND_G(C_a) {return u_a;}

ANNS_FUN(g) {int i;
    double s = 0;
    for (i = 1; i < G_M; i++) {
        s += (1 - exp(-i)) * exp(-i) * sin(i*x[0]) / (i * i);
    }
    return s;}

ANNS_FUN(f) {int i;
    double s = 0;
    for (i = 1; i < G_M; i++) {
        s += exp(-i)*sin(i*x[0]);
    }
    return s;
}

fbvp2_t cond[] = {
    A,
    B, B, B, Btop,
//    C
    };

gbvp2_t cond_a[] = {
     A_a_1,
     B_a, B_a, B_a, B_a,
//      NULL,

    A_a_2, NULL, NULL, NULL, NULL,
//     C_a,
    };

double delta[] = {
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

void ai_init(anns_instance_t *ai) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_reload(anns_instance_t *ai) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 4, 0, g, NULL);
//    ai_eval_v(ai, 0, 0, f, NULL);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (i = ai->I[0]; i < ai->I[1]; ) {
        x[i++] = uniform(0.0,0.0);               // Веса
        x[i++] = uniform(0.05, 1.05);  // Ширины
        x[i++] = uniform(-0.05*M_PI, 1.05*M_PI); // Центры
        x[i++] = uniform(-0.05, 1.05); // Центры
    }
//
    for (i = ai->I[1]; i < ai->I[2]; ) {
        x[i++] = uniform(0.0,1.0);                // Веса
        x[i++] = uniform(0.2, 0.8);  // Ширины
        x[i++] = uniform(-0.05*M_PI, 1.05*M_PI); // Центры
//        x[i++] = uniform(-0.05, 1.05); // Центры
    }
}

double ai_mse(anns_instance_t *bi, double *x) {
    double hx = 0.001;
    int i, n = (int)(M_PI/hx);
    double point[] = {0.0, 0.0};
    double mse = 0, tmp;
//    point[1] = 1.0;
//    for (point[0] = 0.0; point[0] < M_PI; point[0] += hx) {
    for (i = 0; i < n; i++) {
        point[0] = uniform(0.0, M_PI);
        point[1] = uniform(1.0, 1.0);

        bi->ann_val[0][2](point, x, 0, 0, bi);
        tmp = bi->cache->u_xx[0][0] + bi->cache->u_xx[0][1] + f(point, 1, NULL);
//        tmp = bi->cache->u[0] -g(point, 1, NULL);
//        tmp = f(point, 1, NULL) - anns_evalsingle(1, point, x, bi);
        mse += tmp*tmp;
    }

    return sqrt(mse/(n-1));
//return 0;
}

#endif // ANNS_TASK
