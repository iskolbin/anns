#ifndef ANNS_TASK
#define ANNS_TASK

#define PROBLEM_SENSORS 32
#define G_M 100
//0.001
//double usensors[] = {0.000000,0.032557,0.062038,0.091117,0.118739,0.142629,0.167971,0.185000,0.200149,0.214470,0.223786,0.231251,0.234336,0.234491,0.233302,0.230187,0.224006,0.216805,0.206542,0.197106,0.185331,0.171418,0.156656,0.143179,0.125768,0.111135,0.092316,0.075018,0.055980,0.038895,0.018337,0.000000};
//0.01
//double usensors[] = {0.000000,0.023216,0.064338,0.097800,0.101116,0.130217,0.193204,0.183348,0.215233,0.210985,0.217102,0.215902,0.226992,0.233631,0.222099,0.228217,0.239201,0.226999,0.202151,0.187968,0.182759,0.169952,0.187589,0.145638,0.109991,0.102155,0.086553,0.085859,0.058607,0.044572,0.021358,0.000000};
//0.1
double usensors[] = {0.000000,0.039384,0.174621,0.113007,0.025761,0.144354,0.280506,0.119385,0.108032,0.210638,0.207091,0.226350,0.155892,0.273700,0.369279,0.479947,0.269660,0.343766,0.163797,0.199944,0.252433,0.238628,0.167695,0.234982,0.171119,0.133204,0.159144,0.095233,0.150588,0.127139,-0.123750,0.000000};

int data[] = {
    PJ_INVERSE_UNKNOWN_SOURCE, PJ_ELLIPTIC_PROBLEM, 4, 0,

    FUNCTIONAL_SQUARE, 2, 5, 2, 0,

    100, 1,   // Уравнение переноса
    40, 0,   // x = 0
    40, 0,   // x = pi
    20, 0,   // y = 0
    PROBLEM_SENSORS, 1,   // y = 1

    NORMALIZED | RADIAL | GAUSSIAN, 24, 2, 0,   2, 0,  0,  0,  0,
    NORMALIZED | RADIAL | GAUSSIAN, 8,  1, 0,   0, -1, -1, -1, -1,
};

double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID, 0.0, M_PI,   0.0, 2.5,  // Уравнение переноса
    GP_UNIFORM_SOLID_CUBOID, 0.0, 0.0,    0.0,2.5,  // x = 0
    GP_UNIFORM_SOLID_CUBOID, M_PI, M_PI,  0.0, 2.5,  // x = pi

    GP_UNIFORM_SOLID_CUBOID, 0.0, M_PI,   0.0, 0.0,        // y = 0
    GP_GRID_2D,              PROBLEM_SENSORS, 0.0, M_PI,   1, 1.0, 1.0,  // y = 1


};

ANNS_CND(A)  {return u_xx[0][0] + u_xx[0][1] + u[1];}
ANNS_CND(B) {return u[0];}
ANNS_CND(Btop)  {return u[0] - v[0];}

ANNS_CND_G(A_a_1)  {return u_axx[0] + u_axx[1];}
ANNS_CND_G(A_a_2)  {return u_a;}
ANNS_CND_G(B_a) {return u_a;}

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
    };

gbvp2_t cond_a[] = {
     A_a_1,
     B_a, B_a, B_a, B_a,

    A_a_2, NULL, NULL, NULL, NULL,
    };

double delta[] = {
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

void ai_init(anns_instance_t *ai) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_reload(anns_instance_t *ai) {
    ai_genpoints(ai, bounds);
    ai_load_v(ai, 4, 0, usensors);
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
    }
}

double ai_mse(anns_instance_t *bi, double *x) {
    double hx = 0.001;
    int i, n = (int)(M_PI/hx);
    double point[] = {0.0, 0.0};
    double mse = 0, tmp;

    for (i = 0; i < n; i++) {
        point[0] = uniform(0.0, M_PI);
        point[1] = uniform(1.0, 1.0);

        bi->ann_val[0][2](point, x, 0, 0, bi);
        tmp = bi->cache->u_xx[0][0] + bi->cache->u_xx[0][1] + f(point, 1, NULL);
        mse += tmp*tmp;
    }

    return sqrt(mse/(n-1));
}

#endif // ANNS_TASK

