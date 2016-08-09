#ifndef FMIN_LINESEARCH_H_INCLUDED
#define FMIN_LINESEARCH_H_INCLUDED

#include "fmin_vars.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#define FMIN_LINESEARCH_FAIL (DBL_MAX)

typedef struct {
    double c1; // Для правила Армихо
    double c2; // Для условий Вольфа
    double sigma; // Для сжатия и расширения интервала

    double delta1;
    double delta2;

    double tol;    // Точность (для точных методов)
    int max_steps; // Максимальное число шагов
} fmin_linesearch_args_t;


typedef double (*fmin_linesearch_t)(fmin_vars_t *vars, fmin_linesearch_args_t *args);

double fmin_linesearch_goldenratio(fmin_vars_t *vars, fmin_linesearch_args_t *args);
double fmin_linesearch_wolfe(fmin_vars_t *vars, fmin_linesearch_args_t *args);
double fmin_linesearch_cubic(double a, double fa, double fpa, double b, double fb, double c, double fc);
double fmin_linesearch_quad(double a, double fa, double fpa, double b, double fb);
double fmin_linesearch_zoom(double a_lo, double a_hi, double phi_lo, double phi_hi, double derphi_lo, fmin_vars_t *vars, fmin_linesearch_args_t *args);

#endif // FMIN_LINESEARCH_H_INCLUDED
