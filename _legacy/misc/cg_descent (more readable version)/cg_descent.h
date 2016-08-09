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

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

enum {CG_EVAL_FG, CG_EVAL_G, CG_EVAL_F};

typedef struct cg_com_struct /* common variables */
{
    /* parameters computed by the code */
    int              n ; /* problem dimension, saved for reference */
    int             nf ; /* number of function evaluations */
    int             ng ; /* number of gradient evaluations */
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
    double         rho ; /* either Parm->rho or Parm->nan_rho */
    double          *x ; /* current iterate */
    double      *xtemp ; /* x + alpha*d */
    double          *d ; /* current search direction */
    double          *g ; /* gradient at x */
    double      *gtemp ; /* gradient at x + alpha*d */
    double   (*cg_value) (double *, int, void *) ; /* f = cg_value (x, n) */
    void      (*cg_grad) (double *, double *, int, void *) ; /* cg_grad (g, x, n) */
    double (*cg_valgrad) (double *, double *, int, void *) ; /* f = cg_valgrad (g,x,n)*/
    cg_parameter *Parm ; /* user parameters */

    void *instance;
} cg_com ;

/* prototypes */

static int cg_Wolfe
(
    double   alpha, /* stepsize */
    double       f, /* function value associated with stepsize alpha */
    double    dphi, /* derivative value associated with stepsize alpha */
    cg_com    *Com  /* cg com */
) ;

static int cg_tol
(
    double     gnorm, /* gradient sup-norm */
    cg_com    *Com    /* cg com */
) ;

static int cg_line
(
    cg_com   *Com  /* cg com structure */
) ;

static int cg_contract
(
    double    *A, /* left side of bracketing interval */
    double   *fA, /* function value at a */
    double   *dA, /* derivative at a */
    double    *B, /* right side of bracketing interval */
    double   *fB, /* function value at b */
    double   *dB, /* derivative at b */
    cg_com  *Com  /* cg com structure */
) ;


static int cg_evaluate
(
    int    what, /* fg = evaluate func and grad, g = grad only,f = func only*/
    int     nan, /* y means check function/derivative values for nan */
    cg_com   *Com
) ;

static double cg_cubic
(
    double  a,
    double fa, /* function value at a */
    double da, /* derivative at a */
    double  b,
    double fb, /* function value at b */
    double db  /* derivative at b */
) ;

#define cg_dot(dest,x,y,n) { \
    int i; \
    dest = 0; \
    for (i = 0; i < (n); i++) \
        (dest) += (x)[i]*(y)[i]; }

#define cg_copy(dest,src,n) { \
    int j; \
    for (j = 0; j < (n); j++) \
        (dest)[j] = (src)[j]; }

#define cg_step(xtemp,x,d,alpha,n) { \
    int i; \
    for (i = 0; i < (n); i++) \
        (xtemp)[i] = (x)[i] + (alpha)*(d)[i]; }

#define cg_secant(a,b,da,db) \
    (-(da) < (db) ? (a) - ((a)-(b))*((da)/((da)-(db))) : \
        (da) != (b) ? (b) - ((a)-(b))*((db)/((da)-(db))) : -1)

static void cg_printParms
(
    cg_parameter  *Parm
) ;

#endif
