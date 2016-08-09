#include "anns_numeric.h"

// Первая производная: (u(x+dx) - u(x-dx)) / 2dx + O(h^2)
// Вторая производная: (u(x+dx) - 2u(x) + u(x-dx)) / dxdx + O(h^2)
// Смешанная вторая производная: (u(x+dx,y+dy)-u(x+dx,y-dy)-u(x-dy,y+dy)+u(x-dx,y-dy)) / 4dxdy + O(h^2)

void anns_numeric_val_d1(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
//    anns_instance_t *ai = anns_instance;
//    anns_cache_t *ac = ai->cache;
//    anns_func_t ann_val = ai->ann_val[t][0];
//    int j;
//
//    ann_val = ai->ann_val[t][0];
//    // Вычисляем первую производную
//    for (j = 0; j < ai->nd[t]; j++) {
//        // x - dx
//        point[j] -= ANNS_NUMERIC_STEP;
//        ann_val(point, x, t, mixed_derivatives, ai);
//        ac->u_x[t][j] = -ac->u[t];
//
//        // x + dx
//        point[j] += 2 * ANNS_NUMERIC_STEP;
//        ann_val(point, x, t, mixed_derivatives, ai);
//        ac->u_x[t][j] += ac->u[t];
//
//        ac->u_x[t][j] /= 2 * ANNS_NUMERIC_STEP;
//
//        point[j] -= ANNS_NUMERIC_STEP; // Возвращаем точку к исходному виду
//    }
//
//    ann_val(point, x, t, mixed_derivatives, anns_instance);
}

void anns_numeric_val_d2(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
//    anns_instance_t *ai = anns_instance;
//    anns_cache_t *ac = ai->cache;
//    anns_func_t ann_val = ai->ann_val[t][0];
//    int j, k;
//
//    ann_val = ai->ann_val[t][0];
//    // Вычисляем первую и вторую производную
//    for (j = 0; j < ai->nd[t]; j++) {
//        ann_val(point, x, t, mixed_derivatives, ai);
//        ac->u_xx[t][j] = -2*ac->u[t];
//
//        // x - dx
//        point[j] -= ANNS_NUMERIC_STEP;
//        ann_val(point, x, t, mixed_derivatives, ai);
//        ac->u_x[t][j] = -ac->u[t];
//        ac->u_xx[t][j] += ac->u[t];
//        if (mixed_derivatives) {
//            for (k = 0; k < j; k++) {
//                // x - dx, y - dy
//                point[k] -= ANNS_NUMERIC_STEP;
//                ann_val(point, x, t, mixed_derivatives, ai);
//                ac->u_xy[t][j][k] = ac->u[t];
//
//                // x - dx, y + dy
//                point[k] += 2*ANNS_NUMERIC_STEP;
//                ann_val(point, x, t, mixed_derivatives, ai);
//                ac->u_xy[t][j][k] -= ac->u[t];
//
//                point[k] -= ANNS_NUMERIC_STEP; // Возвращаем точку к исходному виду
//            }
//        }
//
//        // x + dx
//        point[j] += 2 * ANNS_NUMERIC_STEP;
//        ann_val(point, x, t, mixed_derivatives, ai);
//        ac->u_x[t][j] += ac->u[t];
//        ac->u_xx[t][j] += ac->u[t];
//
//        ac->u_x[t][j] /= 2 * ANNS_NUMERIC_STEP;
//        ac->u_xx[t][j] /= (ANNS_NUMERIC_STEP*ANNS_NUMERIC_STEP);
//
//        if (mixed_derivatives) {
//            for (k = 0; k < j; k++) {
//                // x + dx, y - dy
//                point[k] -= ANNS_NUMERIC_STEP;
//                ann_val(point, x, t, mixed_derivatives, ai);
//                ac->u_xy[t][j][k] -= ac->u[t];
//
//                // x + dx, y + dy
//                point[k] += 2*ANNS_NUMERIC_STEP;
//                ann_val(point, x, t, mixed_derivatives, ai);
//                ac->u_xy[t][j][k] += ac->u[t];
//
//                point[k] -= ANNS_NUMERIC_STEP; // Возвращаем точку к исходному виду
//
//                ac->u_xy[t][j][k] /= (4*ANNS_NUMERIC_STEP*ANNS_NUMERIC_STEP);
//            }
//        }
//
//        point[j] -= ANNS_NUMERIC_STEP; // Возвращаем точку к исходному виду
//    }
//
//    ann_val(point, x, t, mixed_derivatives, anns_instance);
}

void anns_numeric_valgrad(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
//    anns_instance_t *ai = anns_instance;
//    anns_cache_t *ac = ai->cache;
//    anns_func_t ann_val = ai->ann_val[t][0];
//    int i;
//
//    // x + dx
//    for (i = 0; i < ai->na[t]; i++) {
//        x[i + ai->I[t]] += ANNS_NUMERIC_STEP;
//        ann_val(point, x, t, mixed_derivatives, ai);
//        ac->u_a[t][i] = ac->u[t];
//        x[i + ai->I[t]] -= ANNS_NUMERIC_STEP;
//    }
//
//    // x - dx
//    for (i = 0; i < ai->na[t]; i++) {
//        x[i + ai->I[t]] -= ANNS_NUMERIC_STEP;
//        ann_val(point, x, t, mixed_derivatives, ai);
//        ac->u_a[t][i] -= ac->u[t];
//        ac->u_a[t][i] /= (2*ANNS_NUMERIC_STEP);
//        x[i + ai->I[t]] += ANNS_NUMERIC_STEP;
//    }
//
//    ann_val(point, x, t, mixed_derivatives, ai);
}

void anns_numeric_valgrad_d1(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
//    anns_instance_t *ai = anns_instance;
//    anns_cache_t *ac = ai->cache;
//    int i, j;
//
//    // x + dx
//    for (i = 0; i < ai->na[t]; i++) {
//        x[i + ai->I[t]] += ANNS_NUMERIC_STEP;
//        anns_numeric_val_d1(point, x, t, mixed_derivatives, ai);
//        x[i + ai->I[t]] -= ANNS_NUMERIC_STEP;
//        for (j = 0; j < ai->nd[t]; j++) {
//            ac->u_ax[t][i][j] = ac->u_x[t][j];
//        }
//    }
//
//    // x - dx
//    for (i = 0; i < ai->na[t]; i++) {
//        x[i + ai->I[t]] -= ANNS_NUMERIC_STEP;
//        anns_numeric_val_d1(point, x, t, mixed_derivatives, ai);
//        x[i + ai->I[t]] += ANNS_NUMERIC_STEP;
//        for (j = 0; j < ai->nd[t]; j++) {
//            ac->u_ax[t][i][j] -= ac->u_x[t][j];
//            ac->u_ax[t][i][j] /= (2*ANNS_NUMERIC_STEP);
//        }
//    }
//
//    anns_numeric_valgrad(point, x, t, mixed_derivatives, ai);
}

void anns_numeric_valgrad_d2(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
//    anns_instance_t *ai = anns_instance;
//    anns_cache_t *ac = ai->cache;
//    anns_func_t ann_val = ai->ann_val[t][0];
//    int i, j, k;
//
//    // x + dx
//    for (i = 0; i < ai->na[t]; i++) {
//        x[i + ai->I[t]] += ANNS_NUMERIC_STEP;
//        anns_numeric_val_d2(point, x, t, mixed_derivatives, ai);
//        x[i + ai->I[t]] -= ANNS_NUMERIC_STEP;
//        for (j = 0; j < ai->nd[t]; j++) {
//            ac->u_axx[t][i][j] = ac->u_xx[t][j];
//            if (mixed_derivatives) {
//                for (k = 0; k < j; k++) {
//                    ac->u_axy[t][i][j][k] = ac->u_xy[t][j][k];
//                }
//            }
//        }
//    }
//    // x - dx
//    for (i = 0; i < ai->na[t]; i++) {
//        x[i + ai->I[t]] -= ANNS_NUMERIC_STEP;
//        anns_numeric_val_d2(point, x, t, mixed_derivatives, ai);
//        x[i + ai->I[t]] += ANNS_NUMERIC_STEP;
//        for (j = 0; j < ai->nd[t]; j++) {
//            ac->u_axx[t][i][j] -= ac->u_xx[t][j];
//            ac->u_axx[t][i][j] /= (2*ANNS_NUMERIC_STEP);
//            if (mixed_derivatives) {
//                for (k = 0; k < j; k++) {
//                    ac->u_axy[t][i][j][k] -= ac->u_xy[t][j][k];
//                    ac->u_axy[t][i][j][k] /=(2*ANNS_NUMERIC_STEP);
//                }
//            }
//        }
//    }
//
//    anns_numeric_valgrad_d1(point, x, t, mixed_derivatives, ai);
}

void anns_numeric_val_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {
}

void anns_numeric_val_d1_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {

}

void anns_numeric_val_d2_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {

}

void anns_numeric_valgrad_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {

}

void anns_numeric_valgrad_d1_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {

}

void anns_numeric_valgrad_d2_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance) {

}

double anns_numeric_eval(double *z, double *y, int t, void *anns_instance) {
    return 0;
}

static anns_func_t valgrad_functions[] = {anns_numeric_valgrad, anns_numeric_valgrad_d1, anns_numeric_valgrad_d2, };
static anns_func_t val_wonly_functions[] = {anns_numeric_val_wonly, anns_numeric_val_d1_wonly, anns_numeric_val_d2_wonly, };
static anns_func_t valgrad_wonly_functions[] = {anns_numeric_valgrad_wonly, anns_numeric_valgrad_d1_wonly, anns_numeric_valgrad_d2_wonly, };

void anns_numeric_loadfuncs(anns_net_t *net) {
    anns_func_t *val_functions = calloc(3, sizeof *val_functions);

    val_functions[0] = net->ann_val[0];
    val_functions[1] = anns_numeric_val_d1;
    val_functions[2] = anns_numeric_val_d2;

    net->ann_eval = anns_numeric_eval;
    net->ann_valgrad = valgrad_functions;
    net->ann_val_wonly = val_wonly_functions;
    net->ann_valgrad_wonly = valgrad_wonly_functions;
}

//void anns_numeric_loadfuncs(int t, void *anns_instance) {
//    anns_instance_t *ai = anns_instance;
//    anns_func_t *val_functions = calloc(3, sizeof *val_functions);
//
//    val_functions[0] = ai->ann_val[t][0];
//    val_functions[1] = anns_numeric_val_d1;
//    val_functions[2] = anns_numeric_val_d2;
//
//    ai->ann_val[t] = val_functions;
//    ai->ann_valgrad[t] = valgrad_functions;
//}

