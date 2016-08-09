#ifndef ANNS_TASK
#define ANNS_TASK

double U_b[] = {0.0, 3.0000000000000001e-006, 2.1999999999999999e-005, 8.3999999999999995e-005, 0.000232, 0.000522, 0.0010150000000000001, 0.001776, 0.0028649999999999999, 0.0043369999999999997, 0.0062389999999999998, 0.0086079999999999993, 0.011474, 0.014859000000000001, 0.018778, 0.023238999999999999, 0.028244999999999999, 0.033796, 0.039886999999999999, 0.046509000000000002, 0.053651999999999998, 0.061302000000000002, 0.069445000000000007, 0.078064999999999996, 0.087142999999999998, 0.096660999999999997, 0.106599, 0.116938, 0.12765599999999999, 0.13873099999999999, 0.150143, 0.16186800000000001, 0.17388500000000001, 0.186171, 0.19870399999999999, 0.21145900000000001, 0.224416, 0.23755000000000001, 0.25083899999999998, 0.26425999999999999, 0.27778999999999998, 0.29140700000000003, 0.305087, 0.31880900000000001, 0.33255099999999999, 0.34628900000000001, 0.36000199999999999, 0.37366899999999997, 0.387268, 0.40077699999999999, 0.41417599999999999, 0.42744399999999999, 0.44056000000000001, 0.45350600000000002, 0.46626000000000001, 0.47880400000000001, 0.491118, 0.50318399999999996, 0.514984, 0.52649999999999997, 0.53771400000000003, 0.54861000000000004, 0.559172, 0.56938299999999997, 0.57922700000000005, 0.58869099999999996, 0.59775800000000001, 0.60641599999999996, 0.61465099999999995, 0.62244999999999995, 0.62980100000000006, 0.63669100000000001, 0.64310999999999996, 0.64904700000000004, 0.65449199999999996, 0.65943600000000002, 0.66386900000000004, 0.66778400000000004, 0.67117199999999999, 0.67402799999999996, 0.67634399999999995, 0.67811399999999999, 0.67933399999999999, 0.67999900000000002, 0.68010599999999999, 0.67964999999999998, 0.67863099999999998, 0.67704399999999998, 0.67488999999999999, 0.67216799999999999, 0.66887700000000005, 0.66501900000000003, 0.66059400000000001, 0.65560499999999999, 0.65005400000000002, 0.64394399999999996, 0.63727900000000004, 0.63006399999999996, 0.62230300000000005, 0.61400200000000005};
#define N_UB 100

//#include "../bvp/bvp_ann.h"

// x[0] - x
// x[1] - t

int data[] = {
    PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 4, 0,

    FUNCTIONAL_SQUARE, 2, 5, 1, 0,

    200, 1, // Уравнение переноса
    20, 0,  // Изолированная сторона
    0, 1,  // Нагреваемая сторона
    50, 0,  // Начальное условие

    N_UB, 1, // Измеренная температура на изолированной стороне

    NORMALIZED | RADIAL | GAUSSIAN, 32, 2, 0,  2, 1, 0, 0, 0,
};

#define SOURCE_SWITCH 0

double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID, 0.0, 1.0,  0.0, 1.0, // Уравнение переноса
    GP_UNIFORM_SOLID_CUBOID, 0.0, 0.0,  0.0, 1.0, // Изолированная сторона
    GP_UNIFORM_SOLID_CUBOID, 1.0, 1.0,  0.0, 1.0, // Нагреваемая сторона
    GP_UNIFORM_SOLID_CUBOID, 0.0, 1.0,  0.0, 0.0, // Начальное условие

    GP_GRID_2D, 1, 0.0, 0.0, N_UB, 0.0, 1.0,
};


#if SOURCE_SWITCH == 0
ANNS_FUN(F_q) {return sin(M_PI*x[1]);}
#elif SOURCE_SWITCH == 1
//ANNS_FUN(F_q) {double x_ = x[0]/X_m; return 0.7 + 0.3*sin(2.5*M_PI*x_);}
#else
//ANNS_FUN(F_q) {double x_ = x[0]/X_m; return ((x_ < 0.5) ? 1 : 0.5);}
#endif

ANNS_FUN(F_u) {return 10*x[0]*(1 - x[0]);}
//ANNS_FUN(F_u) {return x[0]*(1 - x[0]);}
//ANNS_FUN(F_u) {return 1.0;}

ANNS_CND(Ctransfer)   {return u_xx[0][0] - v[0]*u_x[0][1];}
ANNS_CND(Cinsulated)  {return u_x[0][0];}
ANNS_CND(Cheated)     {return u[0] - v[0];}
ANNS_CND(Cleft)       {return u[0];}
ANNS_CND(Csensors)    {return u[0] - v[0];}

ANNS_CND_G(Ctransfer_a)  {return u_axx[0] -  v[0]*u_ax[1];}
ANNS_CND_G(Cinsulated_a) {return u_ax[0];}
ANNS_CND_G(Cheated_a)    {return u_a;}
ANNS_CND_G(Cleft_a)      {return u_a;}
ANNS_CND_G(Csensors_a)   {return u_a;}

fbvp2_t cond[] = {Ctransfer, Cinsulated, Cheated, Cleft, Csensors};
gbvp2_t cond_a[] = {Ctransfer_a, Cinsulated_a, Cheated_a, Cleft_a, Csensors_a};
double delta[] = {1.0, 1.0, 1.0, 10.0, 1.0};

void ai_init(anns_instance_t *ai) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_reload(anns_instance_t *ai) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, F_u, NULL);
    ai_eval_v(ai, 2, 0, F_q, NULL);
    ai_load_v(ai, 4, 0, U_b);
}


//void init(bvp_instance_t *bi, double *x) {
//    bi_load_conditions(bi, cond, cond_a, delta);
//    bi_load_netfuncs(bi);
//}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (t = 0; t < ai->Nt; t++) {
        for (i = ai->I[t]; i < ai->I[t+1]; ) {
            x[i++] = uniform(-4., 4.);         // Веса
            x[i++] = uniform(0.5, 4.);        // Ширины
            x[i++] = uniform(-M_PI/2,M_PI+M_PI/2);       // Центры
            x[i++] = uniform(-M_PI/2,M_PI+M_PI/2);       // Центры
        }
    }
}
//
//void update(bvp_instance_t *bi, double *x) {
//    int i, j, t, n, s;
//
//    bi_generate_points(bi, bounds);
//    bi_eval_v(bi, 0, 0, F_u);
//    bi_eval_v(bi, 2, 0, F_q);
//
//    for (t = 0; t < bi->Nt; t++) {
//        for (i = bi->I[t]; i < bi->I[t+1]; ) {
//            x[i++] = uniform(-4.0, 4.);         // Веса
//            x[i++] = uniform(0.05, 1.25);      // Ширины
//            x[i++] = uniform(-0.05, 1.05);     // x
//            x[i++] = uniform(-0.05, 1.05);     // y
//        }
//    }
//}

//void init_ann(double *x, bvp_instance_t *bi) {
//
//}

//double eval_mse(double *x, bvp_instance_t *bi) {
//    return 0.0;
//}
double ai_mse(anns_instance_t *ai, double *x) {
    return 0;
}

#endif // ANNS_TASK
