#ifndef ANNS_TASK
#define ANNS_TASK

#define FD_STEPS_EVAL 161
#define FD_STEPS 161

int data[] = {
    PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 1, 1,

    FUNCTIOAL_SQUARE, 1, 3, 1, FD_STEPS,

    20, 2,  // Уравнение переноса (v[0] - предвычисленные значения f(x), v[1] - предыстория)
    1, 0,   // Изолированная сторона
    1, 1,   // Нагреваемая сторона (v[0] - мощность нагрева)

    NORMALIZED | RADIAL | GAUSSIAN, 5, 1, 0,   2, 1, 1,
};

#define X_m 5.0
#define T_m 10.0
#define U_m 1.0

#define TAU (T_m/(FD_STEPS-1))

#define X_1 (X_m*0.9)

#define C_alpha ( (T_m/(X_m*X_m)) * 202.33576642336 )
#define C_beta  ( (U_m/X_m) * 0.24912184549463 )

double bounds[] = {
    GP_GRID_1D, 20, 0.0, X_m,   // Уравнение переноса
    GP_UNIFORM_SOLID_CUBOID, 0.0, 0.0,     // Изолированная сторона
    GP_UNIFORM_SOLID_CUBOID, X_m, X_m, // Нагреваемая сторона
};

//ANNS_FUN(F_q) {t = x[1]/T_m; return C_beta * ((t <= 0.5) ? (3.0 + 14.0*t) : (10 - 14.0*(t-0.5)));}
//double F_q(double t) {t /= T_m; return C_beta * (7.0 + 3.0*sin(2.5*M_PI*t));}
//ANNS_FUN(F_q) {t = x[1]/T_m; return C_beta * ((t <= 0.5) ? 10.0 : 5.0);}

ANNS_FUN(F_u) {double x_ = x[0]/X_m; return C_alpha * x_ * (1-x_);}

ANNS_CND(A) {return 0.5*TAU*u_xx[0][0] - v[0]*u[0] + v[1];}
ANNS_CND(B) {return u_x[0][0];}
ANNS_CND(C) {return u_x[0][0] - v[0];}

ANNS_CND_G(A_a) {return 0.5*TAU*u_axx[0] - v[0]*u_a;}
ANNS_CND_G(B_a) {return u_ax[0];}
ANNS_CND_G(C_a) {return u_ax[0];}

fbvp_t cond[] = {A, B, C};
gbvp_t cond_a[] = {A_a, B_a, C_a};

double delta[] = {1.0, 1.0, 1.0};

void ai_init(anns_instance_t *ai) {
    ai_conditions(ai, cond, cond_a, NULL);
    ai_loadfuncs(ai);
}

void ai_reload(anns_instance_t *ai) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, F_u, NULL);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (i = ai->I[0]; i < ai->I[1]; ) {
        x[i++] = uniform(U_m, U_m);          // Веса (в начальный момент t = 0)
        x[i++] = uniform(0.05 * X_m, 1.05 * X_m);        // Ширины
        x[i++] = uniform(0.5* X_m, 1.05 * X_m);  // Центры
    }
}

void ai_update(anns_instance_t *ai, double *x, int i) {
    int p;
    anns_cache_t *ac = ai->cache;

    // Сохраняем сеть i-1 шага в контейнер
    ai_put(ai, x, i-1);

    // Обновляем предысторию (v[1]) для конечно-разностной схемы (явно-неявная, Кранк-Никлсон)
    for (p = 0; p < ai->np[0]; p++) {
        ai->ann_val[0][2](ac->points[0][p], x, 0, 0, ai);
        ac->v[0][p][1] = 0.5*TAU*ac->u_xx[0][0] + ac->v[0][p][0]*ac->u[0];
    }

    // Обновляем поток
    ac->v[2][0][0] = F_q(i * TAU);
}

double ai_mse(anns_instance_t *bi, double *x) {
    return 0.0;
}

#endif // ANNS_TASK
