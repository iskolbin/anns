#ifndef FMIN_CG_DESCENT_H_
#define FMIN_CG_DESCENT_H_

#include "cg_user.h"

#include <math.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#define PRIVATE static
#define ZERO 0.0
#define ONE 1.0
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

typedef struct cg_com_struct /* common variables */
{
    /* parameters computed by the code */
    INT              n ; /* problem dimension, saved for reference */
    INT             nf ; /* number of function evaluations */
    INT             ng ; /* number of gradient evaluations */
    int         QuadOK ; /* T (quadratic step successful) */
    int       UseCubic ; /* T (use cubic step) F (use secant step) */
    int           neps ; /* number of time eps updated */
    int       PertRule ; /* T => estimated error in function value is eps*Ck,
                            F => estimated error in function value is eps */
    int          QuadF ; /* T => function appears to be quadratic */
    double   SmallCost ; /* |f| <= SmallCost => set PertRule = F */
    double       alpha ; /* stepsize along search direction */
    double           f ; /* function value for step alpha */
    double          df ; /* function derivative for step alpha */
    double       fpert ; /* perturbation is eps*|f| if PertRule is T */
    double         eps ; /* current value of eps */
    double         tol ; /* computing tolerance */
    double          f0 ; /* old function value */
    double         df0 ; /* old derivative */
    double          Ck ; /* average cost as given by the rule:
                            Qk = Qdecay*Qk + 1, Ck += (fabs (f) - Ck)/Qk */
    double    wolfe_hi ; /* upper bound for slope in Wolfe test */
    double    wolfe_lo ; /* lower bound for slope in Wolfe test */
    double   awolfe_hi ; /* upper bound for slope, approximate Wolfe test */
    int         AWolfe ; /* F (use Wolfe line search)
                                T (use approximate Wolfe line search)
                                do not change user's AWolfe, this value can be
                                changed based on AWolfeFac */
    int          Wolfe ; /* T (means code reached the Wolfe part of cg_line */
    double         rho ; /* either Parm->rho or Parm->An_rho */
    double          *x ; /* current iterate */
    double      *xtemp ; /* x + alpha*d */
    double          *d ; /* current search direction */
    double          *g ; /* gradient at x */
    double      *gtemp ; /* gradient at x + alpha*d */
    double   (*cg_value) (double *, INT, void *) ; /* f = cg_value (x, n) */
    void      (*cg_grad) (double *, double *, INT, void *) ; /* cg_grad (g, x, n) */
    double (*cg_valgrad) (double *, double *, INT, void *) ; /* f = cg_valgrad (g,x,n)*/
    cg_parameter *Parm ; /* user parameters */

    void *instance;
} cg_com ;

/* prototypes */

PRIVATE int cg_Wolfe
(
    double   alpha, /* stepsize */
    double       f, /* function value associated with stepsize alpha */
    double    dphi, /* derivative value associated with stepsize alpha */
    cg_com    *Com  /* cg com */
) ;

PRIVATE int cg_tol
(
    double     gnorm, /* gradient sup-norm */
    cg_com    *Com    /* cg com */
) ;

PRIVATE int cg_line
(
    cg_com   *Com  /* cg com structure */
) ;

PRIVATE int cg_contract
(
    double    *A, /* left side of bracketing interval */
    double   *fA, /* function value at a */
    double   *dA, /* derivative at a */
    double    *B, /* right side of bracketing interval */
    double   *fB, /* function value at b */
    double   *dB, /* derivative at b */
    cg_com  *Com  /* cg com structure */
) ;

PRIVATE int cg_evaluate
(
    char    *what, /* fg = evaluate func and grad, g = grad only,f = func only*/
    char     *An, /* y means check function/derivative values for An */
    cg_com   *Com
) ;

PRIVATE double cg_cubic
(
    double  a,
    double fa, /* function value at a */
    double da, /* derivative at a */
    double  b,
    double fb, /* function value at b */
    double db  /* derivative at b */
) ;

PRIVATE double cg_dot
(
    double *x, /* first vector */
    double *y, /* second vector */
    INT     n  /* length of vectors */
) ;

PRIVATE void cg_copy
(
    double *y, /* output of copy */
    double *x, /* input of copy */
    int     n  /* length of vectors */
) ;

PRIVATE void cg_step
(
    double *xtemp, /*output vector */
    double     *x, /* initial vector */
    double     *d, /* search direction */
    double  alpha, /* stepsize */
    INT         n   /* length of the vectors */
) ;

PRIVATE void cg_printParms
(
    cg_parameter  *Parm
) ;

#endif
