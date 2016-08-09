#ifndef ANNS_NUMERIC_H_INCLUDED
#define ANNS_NUMERIC_H_INCLUDED

#include "anns_instance.h"
#define ANNS_NUMERIC_STEP (1e-3) /* Шаг численного дифференциирования */

// Численное дифференциирование
// Для загрузки функций anns_numeric_loadfuncs. Перед этим надо загрузить непосредственно функцию
// (например с помощью nrbf_gauss_loadfuncs).
// Для производных по пространству используются центральные формулы, по параметру - правые.
// Применять стоит только в отладочных целях, ибо производительность крайне мала.
// Шаг надо выбирать очень разумно! Уже при 1е-5 потеря точности и вторая производная считается неправильно!

void anns_numeric_val_d1(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void anns_numeric_val_d2(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);

void anns_numeric_valgrad(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void anns_numeric_valgrad_d1(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void anns_numeric_valgrad_d2(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);

void anns_numeric_val_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void anns_numeric_val_d1_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void anns_numeric_val_d2_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);

void anns_numeric_valgrad_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void anns_numeric_valgrad_d1_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void anns_numeric_valgrad_d2_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);

void anns_numeric_loadfuncs(anns_net_t *net);
#define anns_numeric_freefuncs(t,ai) (free(ai->ann_val[t]))

#endif // ANNS_NUMERIC_H_INCLUDED
