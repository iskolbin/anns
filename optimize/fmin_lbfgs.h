#ifndef FMIN_LBFGS_H_INCLUDED
#define FMIN_LBFGS_H_INCLUDED

/**
 * L-BFGS optimization.
 *
 * Function prototype
 *
 * int lbfgs(
 *      int n,
 *      lbfgsfloatval_t *x,
 *      lbfgsfloatval_t *ptr_fx,
 *      lbfgs_evaluate_t proc_evaluate,
 *      lbfgs_progress_t proc_progress,
 *      void *instance,
 *      lbfgs_parameter_t *param)
 *
 *  @param  n           The number of variables.
 *  @param  x           The array of variables. A client program can set
 *                      default values for the optimization and receive the
 *                      optimization result through this array. This array
 *                      must be allocated by ::lbfgs_malloc function
 *                      for libLBFGS built with SSE/SSE2 optimization routine
 *                      enabled. The library built without SSE/SSE2
 *                      optimization does not have such a requirement.
 *  @param  ptr_fx      The pointer to the variable that receives the final
 *                      value of the objective function for the variables.
 *                      This argument can be set to \c NULL if the final
 *                      value of the objective function is unnecessary.
 *  @param  proc_evaluate   The callback function to provide function and
 *                          gradient evaluations given a current values of
 *                          variables. A client program must implement a
 *                          callback function compatible with \ref
 *                          lbfgs_evaluate_t and pass the pointer to the
 *                          callback function.
 *  @param  proc_progress   The callback function to receive the progress
 *                          (the number of iterations, the current value of
 *                          the objective function) of the minimization
 *                          process. This argument can be set to \c NULL if
 *                          a progress report is unnecessary.
 *  @param  instance    A user data for the client program. The callback
 *                      functions will receive the value of this argument.
 *  @param  param       The pointer to a structure representing parameters for
 *                      L-BFGS optimization. A client program can set this
 *                      parameter to \c NULL to use the default parameters.
 *                      Call lbfgs_parameter_init() function to fill a
 *                      structure with the default values.
 *  @retval int         The status code. This function returns zero if the
 *                      minimization process terminates without an error. A
 *                      non-zero value indicates an error.
 *
*/

#include "libbfgs/lbfgs.h"
#include "libbfgs/lbfgs_geterror.h"

#endif // FMIN_LBFGS_H_INCLUDED
