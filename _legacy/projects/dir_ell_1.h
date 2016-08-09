#ifndef ANNS_TASK
#define ANNS_TASK

#define N_POINTS_OMEGA 60
#define N_POINTS_GAMMA 15

int data[]={PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 1, 0,
FUNCTIOAL_LINEAR, 2, 5, 1, 0,
N_POINTS_OMEGA, 1, /* A (Transfer equation) */
N_POINTS_GAMMA, 0, /* Btop */
N_POINTS_GAMMA, 0, /* Bleft */
N_POINTS_GAMMA, 0, /* Bright */
N_POINTS_GAMMA, 0, /* Bbottom */
NORMALIZED|RADIAL|GAUSSIAN, 16, 2, 0,  2, 0, 0, 0, 0,};
double bounds[] = {
GP_UNIFORM_SOLID_CUBOID,0,M_PI,0,M_PI,
GP_UNIFORM_SOLID_CUBOID,0,M_PI,M_PI,M_PI,
GP_UNIFORM_SOLID_CUBOID,0,0,0,M_PI,
GP_UNIFORM_SOLID_CUBOID,M_PI,M_PI,0,M_PI,
GP_UNIFORM_SOLID_CUBOID,0,M_PI,0,0};

ANNS_CND(A) {return (u_xx[0][0] + u_xx[0][1] - v[0])*(u_xx[0][0] + u_xx[0][1] - v[0]);}
ANNS_CND(B) {return u[0]*u[0];}

ANNS_CND_G(A_a) {return 2*(u_xx[0][0] + u_xx[0][1] - v[0])*(u_axx[0] + u_axx[1]);}
ANNS_CND_G(B_a) {return 2*u_a*u[0];}

//ANNS_CND(A) {return u_xx[0][0] + u_xx[0][1] - v[0];}
//ANNS_CND(B) {return u[0];}
//
//ANNS_CND_G(A_a) {return u_axx[0] + u_axx[1];}
//ANNS_CND_G(B_a) {return u_a;}

ANNS_FUN(A_f) {return sin(x[0])*sin(x[1]);}

fbvp_t cond[] = {A, B, B, B, B};
gbvp_t cond_a[] = {A_a, B_a, B_a, B_a, B_a};
double delta[] = {1./N_POINTS_OMEGA, 1./N_POINTS_GAMMA, 1./N_POINTS_GAMMA, 1./N_POINTS_GAMMA, 1./N_POINTS_GAMMA};


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
            x[i++] = uniform(-3., 3.);         // Веса
//            x[i++] = uniform(0.5, 4.);        // Ширины
            x[i++] = M_PI/sqrt(ai->nn[0]);
            x[i++] = uniform(-M_PI/4,M_PI+M_PI/4);       // Центры
            x[i++] = uniform(-M_PI/4,M_PI+M_PI/4);       // Центры
        }
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    double s, tmp, point[2];
    int i, p, M = 10000;

    s = 0;
    for (p = 0; p < M; p++) {
        point[0] = uniform(0, M_PI);
        point[1] = uniform(0, M_PI);

        tmp = 0;
        tmp += anns_eval(point, x, ai);
        tmp += 0.5*sin(point[0])*sin(point[1]);

        s += tmp * tmp;
    }

    return sqrt(s / (M-1) );
}

#endif // ANNS_TASK
