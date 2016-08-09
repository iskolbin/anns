#ifndef AURO_OPTIMIZE_FMIN_CG_DESCENT_H_INCLUDED
#define AURO_OPTIMIZE_FMIN_CG_DESCENT_H_INCLUDED

/**
 * CG_DESCENT optimization.
 *
 * Function prototype
 *
 * int cg_descent(
 *      double *x,
 *      int n,
 *      cg_parameter  *UParm,
 *      cg_stats *Stat,
 *      double grad_tol,
 *      double (*value)(double *, int, void *),
 *      void (*grad)(double *, double *, int, void *),
 *      double (*valgrad)(double *, double *, int, void *),
 *      void *instance)
 *
 * @param  x        input: starting guess, output: the solution
 * @param  n        problem dimension
 * @param  UParm    user parameters, NULL = use default parameters
 * @param  Stat     structure with statistics (can be NULL)
 * @param  grad_tol StopRule = 1: |g|_infty <= max (grad_tol,
                    StopFac*initial |g|_infty) [default]
                    topRule = 0: |g|_infty <= grad_tol(1+|f|)
 * @param  value    f = value (x, n, instance)
 * @param  grad     grad (g, x, n)
 * @param  valgrad  f = valgrad (g, x, n, instance), NULL = compute value & gradient using value & grad
 * @param  Work     either size 4n work array or NULL
 * @param  instance additional variable passed to value, grad, valgrad functions (NULL if not needed)
 * @retval int      -2 (function value became nan)
 *                  -1 (starting function value is nan)
 *                  0 (convergence tolerance satisfied)
 *                  1 (change in func <= feps*|f|)
 *                  2 (total iterations exceeded maxit)
 *                  3 (slope always negative in line search)
 *                  4 (number of line search iterations exceeds nline)
 *                  5 (search direction not a descent direction)
 *                  6 (excessive updating of eps)
 *                  7 (Wolfe conditions never satisfied)
 *                  8 (debugger is on and the function value increases)
 *                  9 (out of memory)
 *                  10 (function becomes nan)
 *                  11 (no cost or gradient improvement in 2n + Parm->nslow iterations)
 *
*/

#include "cg_descent/cg_user.h"

#endif // AURO_OPTIMIZE_FMIN_CG_DESCENT_H_INCLUDED
