/* =========================================================================
   ============================ CG_DESCENT =================================
   =========================================================================
       ________________________________________________________________
      |      A conjugate gradient method with guaranteed descent       |
      |             C-code Version 1.1  (October 6, 2005)              |
      |                    Version 1.2  (November 14, 2005)            |
      |                    Version 2.0  (September 23, 2007)           |
      |                    Version 3.0  (May 18, 2008)                 |
      |                    Version 4.0  (March 28, 2011)               |
      |                    Version 4.1  (April 8, 2011)                |
      |                    Version 4.2  (April 14, 2011)               |
      |                    Version 5.0  (May 1, 2011)                  |
      |                                                                |
      |           William W. Hager    and   Hongchao Zhang             |
      |          hager@math.ufl.edu       hzhang@math.ufl.edu          |
      |                   Department of Mathematics                    |
      |                     University of Florida                      |
      |                 Gainesville, Florida 32611 USA                 |
      |                      352-392-0281 x 244                        |
      |                                                                |
      |                 Copyright by William W. Hager                  |
      |                                                                |
      |          http://www.math.ufl.edu/~hager/papers/CG              |
      |________________________________________________________________|
       ________________________________________________________________
      |This program is free software; you can redistribute it and/or   |
      |modify it under the terms of the GNU General Public License as  |
      |published by the Free Software Foundation; either version 2 of  |
      |the License, or (at your option) any later version.             |
      |This program is distributed in the hope that it will be useful, |
      |but WITHOUT ANY WARRANTY; without even the implied warranty of  |
      |MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   |
      |GNU General Public License for more details.                    |
      |                                                                |
      |You should have received a copy of the GNU General Public       |
      |License along with this program; if not, write to the Free      |
      |Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, |
      |MA  02110-1301  USA                                             |
      |________________________________________________________________|

      References:
      1. W. W. Hager and H. Zhang, A new conjugate gradient method
         with guaranteed descent and an efficient line search,
         SIAM Journal on Optimization, 16 (2005), 170-192.
      2. W. W. Hager and H. Zhang, Algorithm 851: CG_DESCENT,
         A conjugate gradient method with guaranteed descent,
         ACM Transactions on Mathematical Software, 32 (2006), 113-137.
      3. W. W. Hager and H. Zhang, A survey of nonlinear conjugate gradient
         methods, Pacific Journal of Optimization, 2 (2006), pp. 35-58. */

#include "cg_descent.h"

int cg_descent /*  return:
                      -2 (function value became nan)
                      -1 (starting function value is nan)
                       0 (convergence tolerance satisfied)
                       1 (change in func <= feps*|f|)
                       2 (total iterations exceeded maxit)
                       3 (slope always negative in line search)
                       4 (number of line search iterations exceeds nline)
                       5 (search direction not a descent direction)
                       6 (excessive updating of eps)
                       7 (Wolfe conditions never satisfied)
                       8 (debugger is on and the function value increases)
                       9 (out of memory)
                      10 (function becomes nan)
                      11 (no cost or gradient improvement in
                          2n + Parm->nslow iterations)*/
(
    double            *x, /* input: starting guess, output: the solution */
    int                n, /* problem dimension */
    cg_stats       *Stat, /* structure with statistics (can be NULL) */
    cg_parameter  *UParm, /* user parameters, NULL = use default parameters */
    double      grad_tol, /* StopRule = 1: |g|_infty <= max (grad_tol,
                                           StopFac*initial |g|_infty) [default]
                             StopRule = 0: |g|_infty <= grad_tol(1+|f|) */
    double      (*value) (double *, int, void *),  /* f = value (x, n) */
    void         (*grad) (double *, double *, int, void *), /* grad (g, x, n) */
    double    (*valgrad) (double *, double *, int, void *), /* f = valgrad (g, x, n),
                          NULL = compute value & gradient using value & grad */
    double         *Work,  /* either size 4n work array or NULL */
    void        *instance  /* additional variable passed to value, grad, valgrad functions (NULL if not needed) */
)
{
    int     n5, i, iter, IterRestart, maxit, nrestart ;
    int     nslow, slowlimit, IterQuad, status, PrintLevel, QuadF, StopRule ;
    double  delta2, Qk, Ck, fbest, gbest,
            f, ftemp, gnorm, xnorm, gnorm2, dnorm2, denom,
            t, dphi, dphi0, alpha,
            yk, ykyk, ykgk, dkyk, beta, QuadTrust, tol,
           *d, *g, *xtemp, *gtemp, *work, tmp ;
    cg_parameter *Parm, ParmStruc ;
    cg_com Com ;

    /* initialize the parameters */
    if ( UParm == NULL ) {
        Parm = &ParmStruc ;
        cg_default (Parm) ;
    }
    else Parm = UParm ;
    PrintLevel = Parm->PrintLevel ;
    Com.Parm = Parm ;
    Com.eps = Parm->eps ;
    Com.PertRule = Parm->PertRule ;
    Com.Wolfe = FALSE ; /* initially Wolfe line search not performed */
    QuadF = FALSE ;     /* initially function assumed to be nonquadratic */

    if ( Parm->PrintParms ) cg_printParms (Parm) ;

    /* allocate work arrays */
    if ( Work == NULL ) {
        work = malloc (4*n*sizeof (double)) ;
    } else {
        work = Work ;
    }

    if ( work == NULL ){
        printf ("Insufficient memory for specified problem dimension %d\n", n) ;
        status = 9 ;
        return (status) ;
    }

    Com.x = x ;
    Com.d = d = work ;
    Com.g = g = d+n ;
    Com.xtemp = xtemp = g+n ;
    Com.gtemp = gtemp = xtemp+n ;
    Com.n = n ;          /* problem dimension */
    Com.nf = 0 ;   /* number of function evaluations */
    Com.ng = 0 ;   /* number of gradient evaluations */
    Com.neps = 0 ;       /* number of times eps updated */
    Com.AWolfe = Parm->AWolfe ; /* do not touch user's AWolfe */
    Com.cg_value = value ;
    Com.cg_grad = grad ;
    Com.cg_valgrad = valgrad ;
    Com.instance = instance;
    StopRule = Parm->StopRule ;

    /* the conjugate gradient algorithm is restarted every nrestart iteration */
    nrestart = (int) n*Parm->restart_fac;

    /* abort when number of iterations reaches maxit */
    if ( Parm->maxit_fac == INF ) {
        maxit = INT_INF ;
    } else {
        maxit = (int) n*Parm->maxit_fac ;
    }

    f = 0;
    fbest = INF ;
    gbest = INF ;
    nslow = 0 ;
    slowlimit = 2*n + Parm->nslow ;

    Ck = 0 ;
    Qk = 0 ;

    /* initial function and gradient evaluations, initial direction */
    Com.alpha = 0 ;
    cg_evaluate (CG_EVAL_FG, 0, &Com) ;
    f = Com.f ;
    Com.f0 = f + f ;
    Com.SmallCost = fabs(f) * Parm->SmallCost ;

    xnorm = 0 ;

    for (i = 0; i < n; i++) {
        tmp = fabs(x[i]);
        if ( xnorm < tmp ) {
            xnorm = tmp ;
        }
    }

    gnorm = 0 ;
    gnorm2 = 0 ;

    for (i = 0; i < n; i++) {
        t = g[i] ;
        gnorm2 += t*t ;
        tmp = fabs(t);
        if ( gnorm < tmp ) {
            gnorm = tmp ;
        }
        d[i] = -t ;
    }

    dnorm2 = gnorm2 ;

    /* check if the starting function value is nan */
    if ( f != f ) {
        status = -1 ;
        goto Exit ;
    }

    if ( Parm->StopRule ) {
        tol = MAX (gnorm*Parm->StopFac, grad_tol) ;
    } else {
        tol = grad_tol ;
    }

    Com.tol = tol ;

    if ( PrintLevel >= 1 ) {
        printf ("iter: %5i f: %13.6e df: %13.6e gnorm: %13.6e\n",
        (int) 0, f, -gnorm2, gnorm) ;
    }

    if ( cg_tol (gnorm, &Com) ) {
        iter = 0 ;
        status = 0 ;
        goto Exit ;
    }

    dphi0 = -gnorm2 ;
    delta2 = 2*Parm->delta - 1 ;
    alpha = Parm->step ;
    if ( alpha == 0 ) {
        alpha = Parm->psi0*xnorm/gnorm ;
        if ( xnorm == 0 ) {
            if ( f != 0 ) {
                alpha = Parm->psi0*fabs (f)/gnorm2 ;
            } else {
                alpha = 1 ;
            }
        }
    }
    IterRestart = 0 ;  /* counts number of iterations since last restart */
    IterQuad = 0 ;     /* counts number of iterations that function change
                          is close to that of a quadratic */

    /* Start the conjugate gradient iteration.
       alpha starts as old step, ends as final step for current iteration
       f is function value for alpha = 0
       QuadOK = TRUE means that a quadratic step was taken */

    for (iter = 1; iter <= maxit; iter++) {
        Com.QuadOK = FALSE ;

        alpha = Parm->psi2*alpha ;

        if ( f != 0 ) {
            t = fabs ((f-Com.f0)/f) ;
        } else {
            t = 1 ;
        }

        Com.UseCubic = TRUE ;
        if ( (t < Parm->CubicCutOff) || !Parm->UseCubic ) {
            Com.UseCubic = FALSE ;
        }

        if ( Parm->QuadStep ) {
            /* test if quadratic interpolation step should be tried */
            if ( ((t > Parm->QuadCutOff)&&(fabs(f) >= Com.SmallCost)) || QuadF ) {
                Com.alpha = Parm->psi1*alpha ;
                if ( QuadF ) {
                    status = cg_evaluate (CG_EVAL_G, 1, &Com) ;
                    if ( status ) {
                        goto Exit ;
                    }
                    if ( Com.df > dphi0 ) {
                        alpha = -dphi0/((Com.df-dphi0)/Com.alpha) ;
                        Com.QuadOK = TRUE ;
                    }
                } else {
                    status = cg_evaluate (CG_EVAL_F, 1, &Com) ;
                    if ( status ) {
                        goto Exit ;
                    }
                    ftemp = Com.f ;
                    denom = 2*(((ftemp-f)/Com.alpha)-dphi0) ;
                    if ( denom > 0 ) {
                        t = -dphi0*Com.alpha/denom ;
                        alpha = t ;
                        if ( (ftemp < f) || QuadF ) {
                            Com.QuadOK = TRUE ;
                        } else { /* safeguard */
                            alpha = MAX (t, Com.alpha*Parm->QuadSafe) ;
                        }
                    }
                }
                if ( PrintLevel >= 1 ) {
                    if ( denom <= 0 ) {
                        printf ("Quad step fails (denom = %14.6e)\n", denom);
                    } else if ( Com.QuadOK ) {
                        printf ("Quad step %14.6e OK\n", alpha);
                    } else printf ("Quad step %14.6e done, but not OK\n", alpha) ;
                }
            } else if ( PrintLevel >= 1 ) {
                printf ("No quad step (chg: %14.6e, cut: %10.2e)\n",
                         t, Parm->QuadCutOff) ;
            }
        }
        Com.f0 = f ;                          /* f0 saved as prior value */
        Com.df0 = dphi0 ;

        /* parameters in Wolfe and approximate Wolfe conditions, and in update*/

        Qk = Parm->Qdecay*Qk + 1 ;
        Ck = Ck + (fabs (f) - Ck)/Qk ;        /* average cost magnitude */

        if ( Com.PertRule ) {
            Com.fpert = f + Com.eps*fabs (f) ;
        } else {
            Com.fpert = f + Com.eps ;
        }

        Com.wolfe_hi = Parm->delta*dphi0 ;
        Com.wolfe_lo = Parm->sigma*dphi0 ;
        Com.awolfe_hi = delta2*dphi0 ;
        Com.alpha = alpha ;

        /* perform line search */
        status = cg_line(&Com) ;

        /*try approximate Wolfe line search if ordinary Wolfe fails */
        if ( (status > 0) && !Com.AWolfe ) {
            if ( PrintLevel >= 1 ) {
                 printf ("\nWOLFE LINE SEARCH FAILS\n") ;
            }
            if ( status != 3 ) {
                Com.AWolfe = TRUE ;
                status = cg_line (&Com) ;
            }
        }

        alpha = Com.alpha ;
        f = Com.f ;
        dphi = Com.df ;

        if ( status ) {
            goto Exit ;
        }

        if ( (f != f) || (dphi != dphi) ) {
            status = 10 ;
            goto Exit ;
        }

        /* Test for convergence to within machine epsilon
           [set feps to zero to remove this test] */

        if ( -alpha*dphi0 <= Parm->feps*fabs(f) ) {
            status = 1 ;
            goto Exit ;
        }

        /* test how close the cost function changes are to that of a quadratic
           QuadTrust = 0 means the function change matches that of a quadratic*/
        t = alpha*(dphi+dphi0) ;
        if ( fabs (t) <= Parm->qeps*MIN (Ck, 1) ) {
            QuadTrust = 0 ;
        } else {
            QuadTrust = fabs((2.0*(f-Com.f0)/t)-1) ;
        }

        if ( QuadTrust <= Parm->qrule) {
            IterQuad++ ;
        } else {
            IterQuad = 0 ;
        }

        if ( IterQuad == Parm->qrestart ) {
            QuadF = TRUE ;
        }

        IterRestart++ ;

        if ( !Com.AWolfe ) {
            if ( fabs (f-Com.f0) < Parm->AWolfeFac*Ck ) {
                Com.AWolfe = TRUE ;
                if ( Com.Wolfe ) IterRestart = nrestart ;
            }
        }

        /* test if the CG algorithm should be restarted */
        if ( (IterRestart == nrestart) || ((IterQuad == Parm->qrestart) && (IterQuad != IterRestart)) ) {
            IterRestart = 0 ;
            IterQuad = 0 ;
            /* search direction d = -g */
            if ( PrintLevel >= 1 ) {
                printf ("RESTART CG\n") ;
            }

            gnorm = 0 ;
            gnorm2 = 0 ;
            cg_copy (x, xtemp, n) ;
            for (i = 0; i < n; i++)
            {
                t = gtemp [i] ;
                tmp = fabs(t);
                if ( gnorm < tmp ) {
                    gnorm = tmp ;
                }
                gnorm2 += t*t ;
                g [i] = t ;
                d [i] = -t ;
            }

            if ( cg_tol (gnorm, &Com) ) {
                status = 0 ;
                goto Exit ;
            }
            dphi0 = -gnorm2 ;
            dnorm2 = gnorm2 ;
        } else { /* compute beta, yk2, gnorm, gnorm2, dnorm2, update x and g */
            cg_copy (x, xtemp, n) ;

            gnorm = 0 ;
            ykyk = 0 ;
            ykgk = 0 ;

            for (i = 0; i < n; i++) {
                t = gtemp[i] ;
                tmp = fabs(t);
                if ( gnorm < tmp ) {
                    gnorm = tmp ;
                }
                yk = t - g[i] ;
                g[i] = t ;
                ykgk += yk*t ;
                ykyk += yk*yk ;
            }

            if ( cg_tol (gnorm, &Com) ) {
                status = 0 ;
                goto Exit ;
            }

            dkyk = dphi - dphi0 ;

            if ( Parm->AdaptiveBeta ) {
                t = 2 - 1/(0.1*QuadTrust + 1) ;
            } else {
                t = Parm->theta ;
            }

            beta = (ykgk - t*dphi*ykyk/dkyk)/dkyk ;

            /* faster: initialize dnorm2 = gnorm2 at start, then
                       dnorm2 = gnorm2 + beta**2*dnorm2 - 2.*beta*dphi
                       gnorm2 = ||g_{k+1}||^2
                       dnorm2 = ||d_{k+1}||^2
                       dpi = g_{k+1}' d_k */

            /* lower bound for beta is BetaLower*d_k'g_k/ ||d_k||^2 */
            beta = MAX (beta, Parm->BetaLower*dphi0/dnorm2) ;

            /* update search direction d = -g + beta*dold */
            gnorm2 = 0 ;
            dnorm2 = 0 ;

            for (i = 0; i < n; i++)
            {
                t = g [i] ;
                gnorm2 += t*t ;
                t = -t + beta*d [i] ;
                d [i] = t ;
                dnorm2 += t*t ;
            }

            dphi0 = -gnorm2 + beta*dphi ;

            if ( Parm->debug ) { /* Check the dphi0 = d'g */
                t = 0 ;
                for (i = 0; i < n; i++)  t = t + d [i]*g [i] ;
                if ( fabs(t-dphi0) > Parm->debugtol*fabs(dphi0) ) {
                    printf("Warning, dphi0 != d'g!\n");
                    printf("dphi0:%13.6e, d'g:%13.6e\n",dphi0, t) ;
                }
            }
        }

        /* test for slow convergence */
        if ( (f < fbest) || (gnorm2 < gbest) ) {
            nslow = 0 ;
            if ( f < fbest ) {
                fbest = f ;
            }
            if ( gnorm2 < gbest ) {
                gbest = gnorm2 ;
            }
        } else {
            nslow++ ;

        }
        if ( nslow > slowlimit ) {
            status = 11 ;
            goto Exit ;
        }

        if ( PrintLevel >= 1 ) {
            printf ("\niter: %5i f = %13.6e gnorm = %13.6e\n",(int) iter, f, gnorm) ;
        }

        if ( Parm->debug ) {
            if ( f > Com.f0 + Parm->debugtol*Ck ) {
                status = 8 ;
                goto Exit ;
            }
        }

        if ( dphi0 > 0 ) {
           status = 5 ;
           goto Exit ;
        }
    }
    status = 2 ;

Exit:
    if ( Stat != NULL ) {
        Stat->f = f ;
        Stat->gnorm = gnorm ;
        Stat->nfunc = Com.nf ;
        Stat->ngrad = Com.ng ;
        Stat->iter = iter ;
    }
    if ( status > 2 ) {
        gnorm = 0 ;
        for (i = 0; i < n; i++)
        {
            x [i] = xtemp [i] ;
            g [i] = gtemp [i] ;
            t = fabs (g [i]) ;
            gnorm = MAX (gnorm, t) ;
        }
        if ( Stat != NULL ) {
            Stat->gnorm = gnorm ;
        }
    }
    if ( Parm->PrintFinal || PrintLevel >= 1 ) {
        const char mess1 [] = "Possible causes of this error message:" ;
        const char mess2 [] = "   - your tolerance may be too strict: "
                              "grad_tol = " ;
        const char mess3 [] = "Line search fails" ;
        const char mess4 [] = "   - your gradient routine has an error" ;
        const char mess5 [] = "   - the parameter epsilon in cg_descent_c.parm "
                              "is too small" ;
        printf ("\nTermination status: %i\n", status) ;
        if ( status == -2 ) {
            printf ("At iteration %10.0f function value became nan\n",
                    (double) iter) ;
        } else if ( status == -1 ) {
            printf ("Objective function value is nan at starting point\n") ;
        } else if ( status == 0 ) {
            printf ("Convergence tolerance for gradient satisfied\n") ;
        } else if ( status == 1 ) {
            printf ("Terminating since change in function value "
                    "<= feps*|f|\n") ;
        } else if ( status == 2 ) {
            printf ("Number of iterations exceed specified limit\n") ;
            printf ("Iterations: %10.0f maxit: %10.0f\n",
                    (double) iter, (double) maxit) ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, grad_tol) ;
        } else if ( status == 3 ) {
            printf ("Slope always negative in line search\n") ;
            printf ("%s\n", mess1) ;
            printf ("   - your cost function has an error\n") ;
            printf ("%s\n", mess4) ;
        } else if ( status == 4 ) {
            printf ("Line search fails, too many iterations\n") ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, grad_tol) ;
        } else if ( status == 5 ) {
            printf ("Search direction not a descent direction\n") ;
        } else if ( status == 6 ) { /* line search fails */
            printf ("%s\n", mess3) ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, grad_tol) ;
            printf ("%s\n", mess4) ;
            printf ("%s\n", mess5) ;
        } else if ( status == 7 ) { /* line search fails */
            printf ("%s\n", mess3) ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, grad_tol) ;
        } else if ( status == 8 ) { /* line search fails */
            printf ("%s\n", mess3) ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, grad_tol) ;
            printf ("%s\n", mess4) ;
            printf ("%s\n", mess5) ;
        } else if ( status == 8 ) {
            printf ("Debugger is on, function value does not improve\n") ;
            printf ("new value: %25.16e old value: %25.16e\n", f, Com.f0) ;
        } else if ( status == 9 ) {
            printf ("Insufficient memory\n") ;
        } else if ( status == 10 ) {
            printf ("Function become nan\n") ;
        } else if ( status == 11 ) {
            printf ("%i iterations without stict improvement in cost or gradient\n", nslow) ;
        }

        printf ("maximum norm for gradient: %13.6e\n", gnorm) ;
        printf ("function value:            %13.6e\n\n", f) ;
        printf ("cg  iterations:          %10.0f\n", (double) iter) ;
        printf ("function evaluations:    %10.0f\n", (double) Com.nf) ;
        printf ("gradient evaluations:    %10.0f\n", (double) Com.ng) ;
        printf ("===================================\n\n") ;
    }

    if ( Work == NULL ) {
        free (work) ;
    }

    return (status) ;
}

/* =========================================================================
   ==== cg_Wolfe ===========================================================
   =========================================================================
   Check whether the Wolfe or the approximate Wolfe conditions are satisfied
   ========================================================================= */
static int cg_Wolfe
(
    double   alpha, /* stepsize */
    double       f, /* function value associated with stepsize alpha */
    double    dphi, /* derivative value associated with stepsize alpha */
    cg_com    *Com  /* cg com */
) {
    if ( dphi >= Com->wolfe_lo ) {
        /* test original Wolfe conditions */
        if ( f - Com->f0 <= alpha*Com->wolfe_hi ) {
            if ( Com->Parm->PrintLevel >= 2 ) {
                printf ("Wolfe conditions hold\n") ;
            }
            return (1) ;
        }
        /* test approximate Wolfe conditions */
        else if ( Com->AWolfe ) {
            if ( (f <= Com->fpert) && (dphi <= Com->awolfe_hi) ) {
                if ( Com->Parm->PrintLevel >= 2 ) {
                    printf ("Approximate Wolfe conditions hold\n") ;
                }
                return (1) ;
            }
        }
    }

    return (0) ;
}

/* =========================================================================
   ==== cg_tol =============================================================
   =========================================================================
   Check for convergence
   ========================================================================= */
static int cg_tol
(
    double     gnorm, /* gradient sup-norm */
    cg_com    *Com    /* cg com */
) {
    /* StopRule = T => |grad|_infty <=max (tol, |grad|_infty*StopFact)
                  F => |grad|_infty <= tol*(1+|f|)) */
    if ( Com->Parm->StopRule ) {
        if ( gnorm <= Com->tol ) {
            return 1 ;
        }
    } else if ( gnorm <= Com->tol*(1 + fabs(Com->f)) ) {
        return (1) ;
    }

    return (0) ;
}

/* =========================================================================
   ==== cg_line ============================================================
   =========================================================================
   Approximate Wolfe line search routine
   Return:
      -2 (function nan)
       0 (convergence tolerance satisfied)
       3 (slope always negative in line search)
       4 (number line search iterations exceed nline)
       6 (excessive updating of eps)
       7 (Wolfe conditions never satisfied)
   ========================================================================= */
static int cg_line
(
    cg_com   *Com /* cg com structure */
) {
    int AWolfe, iter, ngrow, PrintLevel, qb, qb0, status, toggle ;
    double alpha, a, a1, a2, b, B, da, db, d0, d1, d2, dB, df, f, fa, fb, fB,
           a0, b0, da0, db0, fa0, fb0, width, rho, tmp ;
    char *s1, *s2, *fmt1, *fmt2 ;
    cg_parameter *Parm ;

    AWolfe = Com->AWolfe ;
    Parm = Com->Parm ;
    PrintLevel = Parm->PrintLevel ;
    if ( PrintLevel >= 1 ) {
        if ( AWolfe ) {
            printf ("Approximate Wolfe line search\n") ;
            printf ("=============================\n") ;
        } else {
            printf ("Wolfe line search\n") ;
            printf ("=================\n") ;
        }
    }

    /* evaluate function or gradient at Com->alpha (starting guess) */
    b = Com->alpha ;
    if ( Com->QuadOK ) {
        status = cg_evaluate (CG_EVAL_FG, 1, Com) ;
        fb = Com->f ;
        if ( !AWolfe ) fb -= b*Com->wolfe_hi ;
        qb = TRUE ; /* function value at b known */
    } else {
        status = cg_evaluate (CG_EVAL_G, 1, Com) ;
        qb = FALSE ;
    }
    if ( status ) {
        return (status) ; /* function is nan */
    }

    if ( AWolfe ) {
        db = Com->df ;
        d0 = da = Com->df0 ;
    } else {
        db = Com->df - Com->wolfe_hi ;
        d0 = da = Com->df0 - Com->wolfe_hi ;
    }
    a = 0 ;
    a1 = 0 ;
    d1 = d0 ;
    fa = Com->f0 ;
    if ( PrintLevel >= 1 ) {
        fmt1 = "%9s %2s a: %13.6e b: %13.6e fa: %13.6e fb: %13.6e da: %13.6e db: %13.6e\n" ;
        fmt2 = "%9s %2s a: %13.6e b: %13.6e fa: %13.6e fb:  x.xxxxxxxxxx da: %13.6e db: %13.6e\n" ;
        s2 = Com->QuadOK ? "OK" : "" ;
        if ( qb ) {
            printf (fmt1, "start    ", s2, a, b, fa, fb, da, db);
        } else {
            printf(fmt2, "start    ", s2, a, b, fa, da, db) ;
        }
    }

    /* if a quadratic interpolation step performed, check Wolfe conditions */
    if ( (Com->QuadOK) && (Com->f <= Com->f0) ) {
        if ( cg_Wolfe (b, Com->f, Com->df, Com) ) {
            return (0) ;
        }
    }

    /* if a Wolfe line search and the Wolfe conditions have not be satisfied */
    if ( !AWolfe ) {
        Com->Wolfe = TRUE ;
    }

    /*Find initial interval [a,b] such that
      da <= 0, db >= 0, fa <= fpert = [(f0 + eps*fabs (f0)) or (f0 + eps)] */
    rho = Com->rho ;
    ngrow = 1 ;

    while ( db < 0 ) {
        if ( !qb ) {
            cg_evaluate (CG_EVAL_F, 0, Com) ;
            if ( AWolfe ) {
                fb = Com->f ;
            } else {
                fb = Com->f - b*Com->wolfe_hi ;
            }
            qb = TRUE ;
        } if ( fb > Com->fpert ) { /* contract interval [a, b] */
            status = cg_contract (&a, &fa, &da, &b, &fb, &db, Com) ;

            if ( status == 0 ) {
                return (0) ;   /* Wolfe conditions hold */
            }

            if ( status == -2 ) {
                goto Line ; /* db >= 0 */
            }

            if ( Com->neps > Parm->neps ) {
                return (6) ;
            }
        }

        /* expansion phase */
        ngrow++ ;

        if ( ngrow > Parm->nexpand ) {
            return (3) ;
        }

        /* update interval (a replaced by b) */
        a = b ;
        fa = fb ;
        da = db ;

        /* store old values of a and corresponding derivative */
        d2 = d1 ;
        d1 = da ;
        a2 = a1 ;
        a1 = a ;

        if ( (ngrow == 3) || (ngrow == 6) ) {
            if ( d1 > d2 ) {
                tmp = (d1-d2)/(a1-a2) ;
                if ( tmp >= (d2-d0)/a2 ) {
                    /* convex derivative, secant overestimates minimizer */
                    b = a1 - d1/tmp ;
                } else {
                    /* concave derivative, secant underestimates minimizer */
                    b = a1 - Parm->SecantAmp*d1/tmp ;
                }
                /* safeguard growth */
                b = MIN (b, Parm->ExpandSafe*a1) ;
            } else {
                rho *= Parm->RhoGrow ;
                b = rho*b ;
            }
        } else {
            b = rho*b ;
        }

        Com->alpha = b ;
        cg_evaluate (CG_EVAL_G, 0, Com) ;
        qb = FALSE ;

        if ( AWolfe ) {
            db = Com->df ;
        } else {
            db = Com->df - Com->wolfe_hi ;
        }

        if ( PrintLevel >= 2 ) {
            s2 = Com->QuadOK ? "OK" : "";
            printf (fmt2, "expand   ", s2, a, b, fa, da, db) ;
        }
    }

    /* we now have fa <= fpert, da >= 0, db <= 0 */
Line:
    toggle = 0 ;
    width = b - a ;
    qb0 = FALSE ;

    for (iter = 0; iter < Parm->nline; iter++) {
        /* determine the next iterate */
        if ( (toggle == 0) || ((toggle == 2) && ((b-a) <= width)) ) {
            Com->QuadOK = TRUE ;
            if ( Com->UseCubic && qb ) {
                s1 = "cubic    " ;
                alpha = cg_cubic (a, fa, da, b, fb, db) ;
                if ( alpha < 0 ) { /* use secant method */
                    s1 = "secant   " ;
                    alpha = cg_secant(a,b,da,db);
                }
            } else {
                s1 = "secant   " ;
                alpha = cg_secant(a,b,da,db);
            }
            width = Parm->gamma*(b - a) ;
        } else if ( toggle == 1 ) { /* iteration based on smallest value*/
            Com->QuadOK = TRUE ;
            if ( Com->UseCubic ) {
                s1 = "cubic    " ;
                if ( Com->alpha == a ) { /* a is most recent iterate */
                    alpha = cg_cubic (a0, fa0, da0, a, fa, da) ;
                } else if ( qb0 ) { /* b is most recent iterate */
                    alpha = cg_cubic (b, fb, db, b0, fb0, db0) ;
                } else {
                    alpha = -1 ;
                }

                /* if alpha no good, use cubic between a and b */
                if ( (alpha <= a) || (alpha >= b) ) {
                    if ( qb ) {
                        alpha = cg_cubic (a, fa, da, b, fb, db) ;
                    } else {
                        alpha = -1 ;
                    }
                }

                /* if alpha still no good, use secant method */
                if ( alpha < 0 ) {
                    s1 = "secant   " ;
                    alpha = cg_secant(a,b,da,db);
                }
            }
            else  { /* ( use secant ) */
                s1 = "secant   " ;
                if ( (Com->alpha == a) && (da > da0) ) { /* use a0 if possible */
                    alpha = a - (a-a0)*(da/(da-da0)) ;
                } else if ( db < db0 ) { /* use b0 if possible */
                    alpha = b - (b-b0)*(db/(db-db0)) ;
                } else { /* secant based on a and b */
                    alpha = cg_secant(a,b,da,db);
                }

                if ( (alpha <= a) || (alpha >= b) ) {
                    alpha = cg_secant(a,b,da,db);
                }
            }
        } else {
            alpha = 0.5 * (a+b) ; /* use bisection if b-a decays slowly */
            s1 = "bisection" ;
            Com->QuadOK = FALSE ;
        }

        if ( (alpha <= a) || (alpha >= b) ) {
            alpha = 0.5 * (a+b) ;
            s1 = "bisection" ;
            if ( (alpha == a) || (alpha == b) ) {
                return (7) ;
            }
            Com->QuadOK = FALSE ; /* bisection was used */
        }

        if ( toggle == 0 ) { /* save values for next iteration */
            a0 = a ;
            b0 = b ;
            da0 = da ;
            db0 = db ;
            fa0 = fa ;
            if ( qb ) {
                fb0 = fb ;
                qb0 = TRUE ;
            }
        }

        toggle = (toggle + 1) % 3;

        Com->alpha = alpha ;
        cg_evaluate (CG_EVAL_FG, 0, Com) ;
        Com->alpha = alpha ;
        f = Com->f ;
        df = Com->df ;

        if ( Com->QuadOK ) {
            if ( cg_Wolfe (alpha, f, df, Com) ) {
                if ( PrintLevel >= 2 ) {
                    printf ("             a: %13.6e f: %13.6e df: %13.6e %1s\n", alpha, f, df, s1) ;
                }
                return (0) ;
            }
        }

        if ( !AWolfe ) {
            f -= alpha*Com->wolfe_hi ;
            df -= Com->wolfe_hi ;
        }

        if ( df >= 0 ) {
            b = alpha ;
            fb = f ;
            db = df ;
            qb = TRUE ;
        } else if ( f <= Com->fpert ) {
            a = alpha ;
            da = df ;
            fa = f ;
        } else {
            B = b ;
            if ( qb ) {
                fB = fb ;
            }
            dB = db ;
            b = alpha ;
            fb = f ;
            db = df ;
            /* contract interval [a, alpha] */
            status = cg_contract (&a, &fa, &da, &b, &fb, &db, Com) ;

            if ( status == 0 ) {
                return (0) ;
            }

            if ( status == -1 ) { /* eps reduced, use [a, b] = [alpha, b] */
                if ( Com->neps > Parm->neps ) {
                    return (6) ;
                }
                a = b ;
                fa = fb ;
                da = db ;
                b = B ;

                if ( qb ) {
                    fb = fB ;
                }

                db = dB ;
            } else {
                qb = TRUE ;
            }
        }
        if ( PrintLevel >= 2 ) {
            s2 = Com->QuadOK ? "OK" : "";
            if ( !qb ) {
                printf (fmt2, s1, s2, a, b, fa, da, db) ;
            } else {
                printf (fmt1, s1, s2, a, b, fa, fb, da, db) ;
            }
        }
    }

    return (4) ;
}

/* =========================================================================
   ==== cg_contract ========================================================
   =========================================================================
   The input for this routine is an interval [a, b] with the property that
   fa <= fpert, da >= 0, db >= 0, and fb >= fpert. The returned status is

   0  if the Wolfe conditions are satisfied
  -1  if a new value for eps is generated with the property that for the
      corresponding fpert, we have fb <= fpert
  -2  if a subinterval, also denoted [a, b], is generated with the property
      that fa <= fpert, da >= 0, and db <= 0

   NOTE: The input arguments are unchanged when status = -1
   ========================================================================= */
static int cg_contract
(
    double    *A, /* left side of bracketing interval */
    double   *fA, /* function value at a */
    double   *dA, /* derivative at a */
    double    *B, /* right side of bracketing interval */
    double   *fB, /* function value at b */
    double   *dB, /* derivative at b */
    cg_com  *Com  /* cg com structure */
) {
    int AWolfe, iter, PrintLevel, toggle ;
    double a, alpha, b, old, da, db, df, d1, dold, f, fa, fb, f1, fold, t, width ;
    char *s ;
    cg_parameter *Parm ;

    AWolfe = Com->AWolfe ;
    Parm = Com->Parm ;
    PrintLevel = Parm->PrintLevel ;
    a = *A ;
    fa = *fA ;
    da = *dA ;
    b = *B ;
    fb = *fB ;
    db = *dB ;
    f1 = fb ;
    d1 = db ;
    toggle = 0 ;
    width = 0 ;

    for (iter = 0; iter < Parm->nshrink; iter++) {
        if ( (toggle == 0) || ((toggle == 2) && ((b-a) <= width)) ) {
            /* cubic based on bracketing interval */
            alpha = cg_cubic (a, fa, da, b, fb, db) ;
            toggle = 0 ;
            width = Parm->gamma*(b-a) ;
            if ( iter ) {
                Com->QuadOK = TRUE ; /* at least 2 cubic iterations */
            }
        } else if ( toggle == 1 ) {
            Com->QuadOK = TRUE ;
            /* cubic based on most recent iterate and smallest value */
            if ( old < a ) { /* a is most recent iterate */
                alpha = cg_cubic (a, fa, da, old, fold, dold) ;
            } else { /* b is most recent iterate */
                alpha = cg_cubic (a, fa, da, b, fb, db) ;
            }
        } else {
            alpha = 0.5*(a+b) ; /* use bisection if b-a decays slowly */
            Com->QuadOK = FALSE ;
        }

        if ( (alpha <= a) || (alpha >= b) ) {
            alpha = 0.5*(a+b) ;
            Com->QuadOK = FALSE ; /* bisection was used */
        }

        toggle = (toggle + 1) % 3;

        Com->alpha = alpha ;
        cg_evaluate (CG_EVAL_FG, 0, Com) ;
        f = Com->f ;
        df = Com->df ;

        if ( Com->QuadOK ) {
            if ( cg_Wolfe (alpha, f, df, Com) ) {
                return (0) ;
            }
        }

        if ( !AWolfe ) {
            f -= alpha*Com->wolfe_hi ;
            df -= Com->wolfe_hi ;
        }

        if ( df >= 0 ) {
            *B = alpha ;
            *fB = f ;
            *dB = df ;
            *A = a ;
            *fA = fa ;
            *dA = da ;
            return (-2) ;
        }

        if ( f <= Com->fpert ) { /* update a using alpha */
            old = a ;
            a = alpha ;
            fold = fa ;
            fa = f ;
            dold = da ;
            da = df ;
        } else { /* update b using alpha */
            old = b ;
            b = alpha ;
            fb = f ;
            db = df ;
        }

        if ( PrintLevel >= 2 ) {
            s = Com->QuadOK ? "OK" : "";
            printf ("contract  %2s a: %13.6e b: %13.6e fa: %13.6e fb: %13.6e da: %13.6e db: %13.6e\n", s, a, b, fa, fb, da, db) ;
        }
    }

    /* see if the cost is small enough to change the PertRule */
    if ( fabs (fb) <= Com->SmallCost ) {
        Com->PertRule = FALSE ;
    }

    /* increase eps if slope is negative after Parm->nshrink iterations */
    t = Com->f0 ;

    if ( Com->PertRule ) {
        if ( t != 0 ) {
            Com->eps = Parm->egrow*(f1-t)/fabs (t) ;
            Com->fpert = t + fabs (t)*Com->eps ;
        } else {
            Com->fpert = 2.*f1 ;
        }
    } else {
        Com->eps = Parm->egrow*(f1-t) ;
        Com->fpert = t + Com->eps ;
    }

    if ( PrintLevel >= 1 ) {
        printf ("--increase eps: %e fpert: %e\n", Com->eps, Com->fpert) ;
    }

    Com->neps++ ;

    return (-1) ;
}

/* =========================================================================
   ==== cg_fg_evaluate =====================================================
   Evaluate the function and/or gradient.  Also, possibly check if either is nan
   and if so, then reduce the stepsize. Only used at the start of an iteration.
   Return:
      -2 (function nan)
       0 (successful evaluation)
   =========================================================================*/

static int cg_evaluate
(
    int    what, /* fg = evaluate func and grad, g = grad only,f = func only*/
    int     nan, /* y means check function/derivative values for nan */
    cg_com   *Com
) {
    int n ;
    int i ;
    double alpha, *d, *gtemp, *x, *xtemp ;
    cg_parameter *Parm ;

    Parm = Com->Parm ;
    n = Com->n ;
    x = Com->x ;
    d = Com->d ;
    xtemp = Com->xtemp ;
    gtemp = Com->gtemp ;
    alpha = Com->alpha ;
    /* check to see if values are nan */
    if ( nan ) {
        if ( what == CG_EVAL_F ) { /* compute function */
            cg_step (xtemp, x, d, alpha, n) ;
            /* provisional function value */
            Com->f = Com->cg_value (xtemp, n, Com->instance) ;
            Com->nf++ ;

            /* reduce stepsize if function value is nan */
            if ( (Com->f != Com->f) || (Com->f >= INF) ) {
                for (i = 0; i < Parm->nexpand; i++) {
                    alpha *= Parm->nan_decay ;
                    cg_step (xtemp, x, d, alpha, n) ;
                    Com->f = Com->cg_value (xtemp, n, Com->instance) ;
                    Com->nf++ ;
                    if ( (Com->f == Com->f) && (Com->f < INF) ) {
                        break ;
                    }
                }
                if ( i == Parm->nexpand ) {
                    return (-2) ;
                }
            }
            Com->alpha = alpha ;
        } else if ( what == CG_EVAL_G ) { /* compute gradient */
            cg_step (xtemp, x, d, alpha, n) ;

            Com->cg_grad (gtemp, xtemp, n, Com->instance) ;
            Com->ng++ ;
            cg_dot(Com->df, gtemp, d, n) ;

            /* reduce stepsize if derivative is nan */
            if ( (Com->df != Com->df) || (Com->df >= INF) ) {
                for (i = 0; i < Parm->nexpand; i++) {
                    alpha *= Parm->nan_decay ;
                    cg_step (xtemp, x, d, alpha, n) ;
                    Com->cg_grad (gtemp, xtemp, n, Com->instance) ;
                    Com->ng++ ;
                    cg_dot(Com->df, gtemp, d, n) ;
                    if ( (Com->df == Com->df) && (Com->df < INF) ) {
                        break ;
                    }
                }
                if ( i == Parm->nexpand ) {
                    return (-2) ;
                }
                Com->rho = Parm->nan_rho ;
            }
            else {
                Com->rho = Parm->rho ;
            }

            Com->alpha = alpha ;
        } else { /* compute function and gradient */
            cg_step (xtemp, x, d, alpha, n) ;
            if ( Com->cg_valgrad != NULL ) {
                Com->f = Com->cg_valgrad (gtemp, xtemp, n, Com->instance) ;
            } else {
                Com->cg_grad (gtemp, xtemp, n, Com->instance) ;
                Com->f = Com->cg_value (xtemp, n, Com->instance) ;
            }

            cg_dot(Com->df, gtemp, d, n) ;
            Com->nf++ ;
            Com->ng++ ;

            /* reduce stepsize if derivative is nan */
            if ( (Com->df != Com->df) || (Com->f != Com->f) ) {
                for (i = 0; i < Parm->nexpand; i++) {
                    alpha *= Parm->nan_decay ;
                    cg_step (xtemp, x, d, alpha, n) ;
                    if ( Com->cg_valgrad != NULL ) {
                        Com->f = Com->cg_valgrad (gtemp, xtemp, n, Com->instance) ;
                    } else {
                        Com->cg_grad (gtemp, xtemp, n, Com->instance) ;
                        Com->f = Com->cg_value (xtemp, n, Com->instance) ;
                    }
                    cg_dot(Com->df, gtemp, d, n) ;
                    Com->nf++ ;
                    Com->ng++ ;
                    if ( (Com->df == Com->df) && (Com->f == Com->f) ) {
                        break ;
                    }
                }
                if ( i == Parm->nexpand ) {
                    return (-2) ;
                }
                Com->rho = Parm->nan_rho ;
            } else {
                Com->rho = Parm->rho ;
            }

            Com->alpha = alpha ;
        }
    } else { /* evaluate without nan checking */
        if ( what == CG_EVAL_FG ) { /* compute function and gradient */
            if ( alpha == 0 ) { /* evaluate at x */
                if ( Com->cg_valgrad != NULL ) {
                    Com->f = Com->cg_valgrad (Com->g, x, n, Com->instance) ;
                } else {
                    Com->cg_grad (Com->g, x, n, Com->instance) ;
                    Com->f = Com->cg_value (x, n, Com->instance) ;
                }
            } else {
                cg_step (xtemp, x, d, alpha, n) ;
                if ( Com->cg_valgrad != NULL ) {
                    Com->f = Com->cg_valgrad (gtemp, xtemp, n, Com->instance) ;
                } else {
                    Com->cg_grad (gtemp, xtemp, n, Com->instance) ;
                    Com->f = Com->cg_value (xtemp, n, Com->instance) ;
                }
                cg_dot(Com->df, gtemp, d, n) ;
            }
            Com->nf++ ;
            Com->ng++ ;
        } else if ( what == CG_EVAL_F ) { /* compute function */
            cg_step (xtemp, x, d, alpha, n) ;
            Com->f = Com->cg_value (xtemp, n, Com->instance) ;
            Com->nf++ ;
        } else {
            cg_step (xtemp, x, d, alpha, n) ;
            Com->cg_grad (gtemp, xtemp, n, Com->instance) ;
            cg_dot (Com->df, gtemp, d, n) ;
            Com->ng++ ;
        }
    }

    return (0) ;
}

/* =========================================================================
   ==== cg_cubic ===========================================================
   =========================================================================
   Compute the minimizer of a Hermite cubic. If the computed minimizer
   outside [a, b], return -1 (it is assumed that a >= 0).
   ========================================================================= */
static double cg_cubic
(
    double  a,
    double fa, /* function value at a */
    double da, /* derivative at a */
    double  b,
    double fb, /* function value at b */
    double db  /* derivative at b */
){
    double c, d1, d2, delta, t, v, w ;

    delta = b - a ;

    if ( delta == 0 ) {
        c = a ;
    } else {
        v = da + db - 3*(fb-fa)/delta ;
        t = v*v - da*db ;

        if ( t < 0 ) { /* complex roots, use secant method */
             if ( fabs (da) < fabs (db) ) {
                 c = a - (a-b)*(da/(da-db)) ;
             } else if ( da != db ) {
                 c = b - (a-b)*(db/(da-db)) ;
             } else {
                 c = -1 ;
             }
        } else {

            w = delta > 0 ? sqrt(t) : -sqrt(t);

            d1 = da + v - w ;
            d2 = db + v + w ;

            if ( (d1 == 0) && (d2 == 0) ) {
                c = -1 ;
            } else {
                if ( fabs (d1) >= fabs (d2) ) {
                    c = a + delta*da/d1 ;
                } else {
                    c = b - delta*db/d2 ;
                }
            }
        }
    }

    return (c) ;
}

/* =========================================================================
   === cg_default ==========================================================
   =========================================================================
   Set default conjugate gradient parameter values. If the parameter argument
   of cg_descent is NULL, this routine is called by cg_descent automatically.
   If the user wishes to set parameter values, then the cg_parameter structure
   should be allocated in the main program. The user could call cg_default
   to initialize the structure, and then individual elements in the structure
   could be changed, before passing the structure to cg_descent.
   =========================================================================*/
void cg_default
(
    cg_parameter   *Parm
) {
    /* T => print final function value
       F => no printout of final function value */
    Parm->PrintFinal = TRUE ;

   /* Level 0 = no printing, ... , Level 3 = maximum printing */
    Parm->PrintLevel = 0 ;

    /* T => print parameters values
       F => do not display parmeter values */
    Parm->PrintParms = FALSE ;

    /* T => use approximate Wolfe line search
       F => use ordinary Wolfe line search, switch to approximate Wolfe when
                |f_k+1-f_k| < AWolfeFac*C_k, C_k = average size of cost */
    Parm->AWolfe = FALSE ;
    Parm->AWolfeFac = 1.e-3 ;

    /* factor in [0, 1] used to compute average cost magnitude C_k as follows:
       Q_k = 1 + (Qdecay)Q_k-1, Q_0 = 0,  C_k = C_k-1 + (|f_k| - C_k-1)/Q_k */
    Parm->Qdecay = .7 ;

    /* terminate after 2*n + nslow iterations without strict improvement in
       either function value or gradient */
    Parm->nslow = 1000 ;

    /* Stop Rules:
       T => ||grad||_infty <= max(grad_tol, initial |grad|_infty*StopFact)
       F => ||grad||_infty <= grad_tol*(1 + |f_k|) */
    Parm->StopRule = TRUE ;
    Parm->StopFac = 0.e-12 ;

    /* T => estimated error in function value is eps*Ck,
       F => estimated error in function value is eps */
    Parm->PertRule = TRUE ;
    Parm->eps = 1.e-6 ;

    /* factor by which eps grows when line search fails during contraction */
    Parm->egrow = 10. ;

    /* T => attempt quadratic interpolation in line search when
                |f_k+1 - f_k|/|f_k| > QuadCutOff
       F => no quadratic interpolation step */
    Parm->QuadStep = TRUE ;
    Parm->QuadCutOff = 1.e-12 ;

    /* maximum factor by which a quad step can reduce the step size */
    Parm->QuadSafe = 1.e-3 ;

    /* T => when possible, use a cubic step in the line search */
    Parm->UseCubic = TRUE ;

    /* use cubic step when |f_k+1 - f_k|/|f_k| > CubicCutOff */
    Parm->CubicCutOff = 1.e-12 ;

    /* |f| < SmallCost*starting cost => skip QuadStep and set PertRule = FALSE*/
    Parm->SmallCost = 1.e-30 ;

    /* T => check that f_k+1 - f_k <= debugtol*C_k
       F => no checking of function values */
    Parm->debug = FALSE ;
    Parm->debugtol = 1.e-10 ;

    /* if step is nonzero, it is the initial step of the initial line search */
    Parm->step = 0 ;

    /* abort cg after maxit_fac*n iterations */
    Parm->maxit_fac = INF ;

    /* maximum number of times the bracketing interval grows during expansion */
    Parm->nexpand = (int) 50 ;

    /* maximum factor secant step increases stepsize in expansion phase */
    Parm->ExpandSafe = 200. ;

    /* factor by which secant step is amplified during expansion phase
       where minimizer is bracketed */
    Parm->SecantAmp = 1.05 ;

    /* factor by which rho grows during expansion phase where minimizer is
       bracketed */
    Parm->RhoGrow = 2.0 ;

    /* maximum number of times that eps is updated */
    Parm->neps = (int) 5 ;

    /* maximum number of times the bracketing interval shrinks */
    Parm->nshrink = (int) 10 ;

    /* maximum number of secant iterations in line search is nline */
    Parm->nline = (int) 50 ;

    /* conjugate gradient method restarts after (n*restart_fac) iterations */
    Parm->restart_fac = 6.0 ;

    /* stop when -alpha*dphi0 (estimated change in function value) <= feps*|f|*/
    Parm->feps = 0 ;

    /* after encountering nan, growth factor when searching for
       a bracketing interval */
    Parm->nan_rho = 1.3 ;

    /* after encountering nan, decay factor for stepsize */
    Parm->nan_decay = 0.1 ;

    /* Wolfe line search parameter, range [0, .5]
       phi (a) - phi (0) <= delta phi'(0) */
    Parm->delta = .1 ;

    /* Wolfe line search parameter, range [delta, 1]
       phi' (a) >= sigma phi' (0) */
    Parm->sigma = .9 ;

    /* decay factor for bracket interval width in line search, range (0, 1) */
    Parm->gamma = .66 ;

    /* growth factor in search for initial bracket interval */
    Parm->rho = 5. ;

    /* starting guess for line search =
         psi0 ||x_0||_infty over ||g_0||_infty if x_0 != 0
         psi0 |f(x_0)|/||g_0||_2               otherwise */
    Parm->psi0 = .01 ;      /* factor used in starting guess for iteration 1 */

    /* for a QuadStep, function evalutated at psi1*previous step */
    Parm->psi1 = .1 ;

    /* when starting a new cg iteration, our initial guess for the line
       search stepsize is psi2*previous step */
    Parm->psi2 = 2. ;

    /* choose theta adaptively if AdaptiveBeta = T */
    Parm->AdaptiveBeta = FALSE ;

    /* lower bound for beta is BetaLower*d_k'g_k/ ||d_k||^2 */
    Parm->BetaLower = 0.4 ;

    /* value of the parameter theta in the cg_descent update formula:
       W. W. Hager and H. Zhang, A survey of nonlinear conjugate gradient
       methods, Pacific Journal of Optimization, 2 (2006), pp. 35-58. */
    Parm->theta = 1.0 ;

    /* parameter used in cost error estimate for quadratic restart criterion */
    Parm->qeps = 1.e-12 ;

    /* number of iterations the function is nearly quadratic before a restart */
    Parm->qrestart = 3 ;

    /* treat cost as quadratic if
       |1 - (cost change)/(quadratic cost change)| <= qrule */
    Parm->qrule = 1.e-8 ;
}

/* =========================================================================
   ==== cg_printParms ======================================================
   =========================================================================
   Print the contents of the cg_parameter structure
   ========================================================================= */
static void cg_printParms
(
    cg_parameter  *Parm
) {
    printf ("PARAMETERS:\n") ;
    printf ("\n") ;
    printf ("Wolfe line search parameter ..................... delta: %e\n",
             Parm->delta) ;
    printf ("Wolfe line search parameter ..................... sigma: %e\n",
             Parm->sigma) ;
    printf ("decay factor for bracketing interval ............ gamma: %e\n",
             Parm->gamma) ;
    printf ("growth factor for bracket interval ................ rho: %e\n",
             Parm->rho) ;
    printf ("growth factor for bracket interval after nan .. nan_rho: %e\n",
             Parm->nan_rho) ;
    printf ("decay factor for stepsize after nan ......... nan_decay: %e\n",
             Parm->nan_decay) ;
    printf ("parameter in lower bound for beta ........... BetaLower: %e\n",
             Parm->BetaLower) ;
    printf ("parameter describing cg_descent family .......... theta: %e\n",
             Parm->theta) ;
    printf ("perturbation parameter for function value ......... eps: %e\n",
             Parm->eps) ;
    printf ("factor by which eps grows if necessary .......... egrow: %e\n",
             Parm->egrow) ;
    printf ("factor for computing average cost .............. Qdecay: %e\n",
             Parm->Qdecay) ;
    printf ("relative change in cost to stop quadstep ... QuadCutOff: %e\n",
             Parm->QuadCutOff) ;
    printf ("maximum factor quadstep reduces stepsize ..... QuadSafe: %e\n",
             Parm->QuadSafe) ;
    printf ("skip quadstep if |f| <= SmallCost*start cost  SmallCost: %e\n",
             Parm->SmallCost) ;
    printf ("relative change in cost to stop cubic step  CubicCutOff: %e\n",
             Parm->CubicCutOff) ;
    printf ("terminate if no improvement over nslow iter ..... nslow: %i\n",
             Parm->nslow) ;
    printf ("factor multiplying gradient in stop condition . StopFac: %e\n",
             Parm->StopFac) ;
    printf ("cost change factor, approx Wolfe transition . AWolfeFac: %e\n",
             Parm->AWolfeFac) ;
    printf ("restart cg every restart_fac*n iterations . restart_fac: %e\n",
             Parm->restart_fac) ;
    printf ("cost error in quadratic restart is qeps*cost ..... qeps: %e\n",
             Parm->qeps) ;
    printf ("number of quadratic iterations before restart  qrestart: %i\n",
             Parm->qrestart) ;
    printf ("parameter used to decide if cost is quadratic ... qrule: %e\n",
             Parm->qrule) ;
    printf ("stop when cost change <= feps*|f| ................ feps: %e\n",
             Parm->feps) ;
    printf ("starting guess parameter in first iteration ...... psi0: %e\n",
             Parm->psi0) ;
    printf ("starting step in first iteration if nonzero ...... step: %e\n",
             Parm->step) ;
    printf ("factor multiply starting guess in quad step ...... psi1: %e\n",
             Parm->psi1) ;
    printf ("initial guess factor for general iteration ....... psi2: %e\n",
             Parm->psi2) ;
    printf ("max iterations is n*maxit_fac ............... maxit_fac: %e\n",
             Parm->maxit_fac) ;
    printf ("max number of contracts in the line search .... nshrink: %i\n",
             Parm->nshrink) ;
    printf ("max expansions in line search ................. nexpand: %i\n",
             Parm->nexpand) ;
    printf ("maximum growth of secant step in expansion . ExpandSafe: %e\n",
             Parm->ExpandSafe) ;
    printf ("growth factor for secant step during expand . SecantAmp: %e\n",
             Parm->SecantAmp) ;
    printf ("growth factor for rho during expansion phase .. RhoGrow: %e\n",
             Parm->RhoGrow) ;
    printf ("max number of times that eps is updated .......... neps: %i\n",
             Parm->neps) ;
    printf ("max number of iterations in line search ......... nline: %i\n",
             Parm->nline) ;
    printf ("print level (0 = none, 2 = maximum) ........ PrintLevel: %i\n",
             Parm->PrintLevel) ;
    printf ("Logical parameters:\n") ;
    if ( Parm->PertRule )
        printf ("    Error estimate for function value is eps*Ck\n") ;
    else
        printf ("    Error estimate for function value is eps\n") ;
    if ( Parm->QuadStep )
        printf ("    Use quadratic interpolation step\n") ;
    else
        printf ("    No quadratic interpolation step\n") ;
    if ( Parm->UseCubic)
        printf ("    Use cubic interpolation step when possible\n") ;
    else
        printf ("    Avoid cubic interpolation steps\n") ;
    if ( Parm->AdaptiveBeta )
        printf ("    Adaptively adjust direction update parameter beta\n") ;
    else
        printf ("    Use fixed parameter theta in direction update\n") ;
    if ( Parm->PrintFinal )
        printf ("    Print final cost and statistics\n") ;
    else
        printf ("    Do not print final cost and statistics\n") ;
    if ( Parm->PrintParms )
        printf ("    Print the parameter structure\n") ;
    else
        printf ("    Do not print parameter structure\n") ;
    if ( Parm->AWolfe)
        printf ("    Approximate Wolfe line search\n") ;
    else
        printf ("    Wolfe line search") ;
        if ( Parm->AWolfeFac > 0. )
            printf (" ... switching to approximate Wolfe\n") ;
        else
            printf ("\n") ;
    if ( Parm->StopRule )
        printf ("    Stopping condition uses initial grad tolerance\n") ;
    else
        printf ("    Stopping condition weighted by absolute cost\n") ;
    if ( Parm->debug)
        printf ("    Check for decay of cost, debugger is on\n") ;
    else
        printf ("    Do not check for decay of cost, debugger is off\n") ;
}

/*
Version 1.2 Change:
  1. The variable dpsi needs to be included in the argument list for
     subroutine cg_updateW (update of a Wolfe line search)

Version 2.0 Changes:
     The user interface was redesigned. The parameters no longer need to
     be read from a file. For compatibility with earlier versions of the
     code, we include the routine cg_readParms to read parameters.
     In the simplest case, the user can use NULL for the
     parameter argument of cg_descent, and the code sets the default
     parameter values. If the user wishes to modify the parameters, call
     cg_default in the main program to initialize a cg_parameter
     structure. Individual elements of the structure could be modified.
     The header file cg_user.h contains the structures and prototypes
     that the user may need to reference or modify, while cg_descent.h
     contains header elements that only cg_descent will access.  Note
     that the arguments of cg_descent have changed.

Version 3.0 Changes:
     Major overhaul

Version 4.0 Changes:
  1. Set theta = 1.0 by default in the cg_descent rule for beta_k
  2. Increase the default value of restart_fac to 6 (a value larger than 1 is
     more efficient when the problem dimension is small)
  3. Restart the CG iteration if the objective function is nearly quadratic
     for several iterations (qrestart). This type of restart is described
     in the paper "A nonlinear conjugate gradient algorithm with an
     optimal property and an improved Wolfe line search by Yu-Hong Dai and
     Cai-Xia Kou.
  4. New lower bound for beta: BetaLower*d_k'g_k/ ||d_k||^2
  5. Evaluation of the objective function and gradient is now handled by
     the routine cg_evaluate.

Version 4.1 Changes:
  1. Change cg_tol to be consistent with corresponding routine in asa_cg
  2. Compute dnorm2 when d is evaluated and make loops consistent with asa_cg

Version 4.2 Changes:
  1. Modify the line search so that when there are too many contractions,
     the code will increase eps and switch to expansion of the search interval.
     This fixes some cases where the code terminates when eps is too small.
     When the estimated error in the cost function is too small, the algorithm
     could fail in cases where the slope is negative at both ends of the
     search interval and the objective function value on the right side of the
     interval is larger than the value at the left side (because the true
     objective function value on the right is not greater than the value on
     the left).
  2. Fix bug in cg_lineW

Version 5.0 Changes:
     Revise the line search routines to exploit steps based on the
     minimizer of a Hermite interpolating cubic. Combine the approximate
     and the ordinary Wolfe line search into a single routine.
     Include safeguarded extrapolation during the expansion phase
     of the line search. Employ a quadratic interpolation step even
     when the requirement ftemp < f for a quadstep is not satisfied.
*/
