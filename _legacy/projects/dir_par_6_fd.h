#ifndef ANNS_TASK
#define ANNS_TASK

#define C_l1 1.0
#define C_l2 0.5
#define C_T  1.0

#define FD_STEPS_EVAL 101
#define FD_STEPS 101

#define C_tau (C_T/(FD_STEPS_EVAL-1))

int data[] = {
    PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 6, 1,

    FUNCTIONAL_SQUARE, 2, 5, 1, FD_STEPS,

    500, 1, /* A_omega */

    25, 0, /* A_Gamma1 */
    25, 0, /* A_Gamma2 */
    25, 0, /* A_gamma_star */
    25, 1, /* A_gamma */

    NORMALIZED | RADIAL | GAUSSIAN, 16, 2, 0,   2, 1, 1, 1, 1,
};

double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID, 0, C_l1,    0, C_l2,

    GP_UNIFORM_SOLID_CUBOID, 0, 0,       0, C_l2,
    GP_UNIFORM_SOLID_CUBOID, C_l1, C_l1, 0, C_l2,
    GP_UNIFORM_SOLID_CUBOID, 0, C_l1,    0, 0,
    GP_UNIFORM_SOLID_CUBOID, 0, C_l1,    C_l2, C_l2,
};

ANNS_CND(A_omega) {return 0.5*C_tau*(u_xx[0][0] + u_xx[0][1]) - u[0] + v[0];}
ANNS_CND(A_Gamma) {return u_x[0][0];}
ANNS_CND(A_gamma) {return u_x[0][1] - v[0];}
ANNS_CND(A_gamma_star) {return u_x[0][1];}


ANNS_CND_G(A_omega_a) {return 0.5*C_tau*(u_axx[0] + u_axx[1]) - u_a;}
ANNS_CND_G(A_Gamma_a) {return u_ax[0];}
ANNS_CND_G(A_gamma_a) {return u_ax[1];}

//ANNS_FUN(F_mu_1) {return x[0]*x[2]/C_T;}

//double delta[] = {1.0, 1.0, 1.0};
fbvp2_t cond[] = {A_omega, A_Gamma, A_Gamma, A_gamma_star, A_gamma};
gbvp2_t cond_a[] = {A_omega_a, A_Gamma_a, A_Gamma_a, A_gamma_a, A_gamma_a};

void ai_init(anns_instance_t *ai) {
    int p;
    ai_conditions(ai, cond, cond_a, NULL);
    ai_loadfuncs(ai);
    ai_genpoints(ai, bounds);
    for (p = 0; p < ai->np[4]; p++) {
        ai->cache->v[4][p][0] = 0;
    }
}

void ai_reload(anns_instance_t *ai) {
    int p;

}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (i = ai->I[0]; i < ai->I[1]; ) {
        x[i++] = 0.0;          // Веса (в начальный момент t = 0)
        x[i++] = uniform(0.25, 1.25);        // Ширины
//        x[i++] = 0.5/sqrt(ai->nn[0]);
        x[i++] = uniform(-0.25, C_l1 + 0.25);  // Центры
        x[i++] = uniform(-0.25, C_l2 + 0.25);
    }
}

void ai_update(anns_instance_t *ai, double *x, int i) {
    int p;
    anns_cache_t *ac = ai->cache;

    // Сохраняем сеть i-1 шага в контейнер
    ai_put(ai, x, i-1);

    // Обновляем предысторию (v[0]) для конечно-разностной схемы (явно-неявная, Кранк-Никлсон)
    for (p = 0; p < ai->np[0]; p++) {
        ai->ann_val[0][2](ac->points[0][p], x, 0, 0, ai);
        ac->v[0][p][0] = 0.5*C_tau*(ac->u_xx[0][0] + ac->u_xx[0][1]) + ac->u[0];
    }

    // Обновляем поток
    for (p = 0; p < ai->np[4]; p++) {
        ac->v[4][p][0] = ai->cache->points[4][p][0]*C_tau*i/C_T;
    }
}

double ai_mse(anns_instance_t *bi, double *x) {
    return -1;
}

#endif // ANNS_TASK
