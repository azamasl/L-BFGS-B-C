/* Tues, Feb 17 2015 
 * Stephen Becker, stephen.becker@colorado.edu
 * */

#ifndef lbfgsb_h
#define lbfgsb_h



 /* You could have to modify these
 * Noticed that on windows, long int is 32-bit
 * while on linux and mac long int is 64-bit.
  *Use long long to force 64-bit if you want */
typedef long int integer;
typedef long int ftnlen;
typedef long int logical;
#define TRUE_ (1)
#define FALSE_ (0)


#include <math.h>
#include <stdio.h>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
/* mex.h includes matrix.h, which uses abs() */

/* Visual C++ compiler often includes file that have already
 * defined min/max, so do a check first */
/* The below cause conflicts with other versions of abs, max, min.
   Replaced with fabs, fmin, fmax provided by math.h */
/* #ifndef abs */
/* #define abs(x) ((x) >= 0 ? (x) : -(x)) */
/* #endif */
/* #ifndef min */
/* #define min(a,b) ((a) <= (b) ? (a) : (b)) */
/* #endif */
/* #ifndef max */
/* #define max(a,b) ((a) >= (b) ? (a) : (b)) */
/* #endif */




/* Modified L-BFGS-B to use integers instead
 * of strings, for testing the "task". Way
 * simpler this way. For ease-of-use, use
 * these aliases. For each class of task,
 * make sure numbering is contiguous
 * so that I have easy tests for class
 * membership */
#define START 1 
#define NEW_X 2
#define ABNORMAL 3 /* message: ABNORMAL_TERMINATION_IN_LNSRCH. */
#define RESTART 4 /* message: RESTART_FROM_LNSRCH. */

#define FG      10
#define FG_END  15
#define IS_FG(x) ( ((x)>=FG) ?  ( ((x)<=FG_END) ? 1 : 0 ) : 0 )
#define FG_LN   11
#define FG_LNSRCH FG_LN
#define FG_ST   12
#define FG_START FG_ST

#define CONVERGENCE 20
#define CONVERGENCE_END  25
#define IS_CONVERGED(x) ( ((x)>=CONVERGENCE) ?  ( ((x)<=CONVERGENCE_END) ? 1 : 0 ) : 0 )
#define CONV_GRAD   21 /* message: CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL. */
#define CONV_F      22 /* message: CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH. */

#define STOP  30  
#define STOP_END 40
#define IS_STOP(x) ( ((x)>=STOP) ?  ( ((x)<=STOP_END) ? 1 : 0 ) : 0 )
#define STOP_CPU  31 /* message: STOP: CPU EXCEEDING THE TIME LIMIT. */
#define STOP_ITER 32 /* message: STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIM.  */
#define STOP_GRAD 33 /* message: STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL. */

#define WARNING 100
#define WARNING_END 110
#define IS_WARNING(x) ( ((x)>=WARNING) ?  ( ((x)<=WARNING_END) ? 1 : 0 ) : 0 )
#define WARNING_ROUND 101  /* WARNING: ROUNDING ERRORS PREVENT PROGRESS */
#define WARNING_XTOL  102  /* WARNING: XTOL TEST SATISIED */
#define WARNING_STPMAX 103 /* WARNING: STP = STPMAX */
#define WARNING_STPMIN 104 /* WARNING: STP = STPMIN */

#define ERROR 200
#define ERROR_END 240
#define IS_ERROR(x) ( ((x)>=ERROR) ?  ( ((x)<=ERROR_END) ? 1 : 0 ) : 0 )
/* More specific conditions below */
#define ERROR_SMALLSTP 201 /* message: ERROR: STP .LT. STPMIN  */
#define ERROR_LARGESTP 202 /* message: ERROR: STP .GT. STPMAX  */
#define ERROR_INITIAL  203 /* message: ERROR: INITIAL G .GE. ZERO */
#define ERROR_FTOL     204 /* message: ERROR: FTOL .LT. ZERO   */
#define ERROR_GTOL     205 /* message: ERROR: GTOL .LT. ZERO   */
#define ERROR_XTOL     206 /* message: ERROR: XTOL .LT. ZERO   */
#define ERROR_STP0     207 /* message: ERROR: STPMIN .LT. ZERO */
#define ERROR_STP1     208 /* message: ERROR: STPMAX .LT. STPMIN */
#define ERROR_N0       209 /* ERROR: N .LE. 0 */
#define ERROR_M0       210 /* ERROR: M .LE. 0 */
#define ERROR_FACTR    211 /* ERROR: FACTR .LT. 0 */
#define ERROR_NBD      212 /* ERROR: INVALID NBD */
#define ERROR_FEAS     213 /* ERROR: NO FEASIBLE SOLUTION */


/* and "word" was a char that was one fo these: */
#define WORD_DEFAULT 0 /* aka "---".  */
#define WORD_CON 1 /*  the subspace minimization converged. */
#define WORD_BND 2 /* the subspace minimization stopped at a bound. */
#define WORD_TNT 3 /*  the truncated Newton step has been used. */


/* If we are linking with fortran code,
 * then use gfortran to compile, not gcc,
 * and it will expect function names to have
 * underscores.
 * With gcc, if we execute the following,
 * it will complain _daxpy_ is undefined ... */
#if defined(_WIN32) || defined(__hpux) || !defined(__GFORTRAN__)
#define FORTRAN_WRAPPER(x) x
#else
#define FORTRAN_WRAPPER(x) x ## _
/* if we're not WIN32 or HPUX, then if we are linking
 * with gfortran (instead of gcc or g++), we need to mangle
 * names appropriately */
#endif


/* First decision: use the included BLAS code
 * (un-optimized version taken from Netlib;
 * this is the "reference" implementation, or
 * "Ref" for short; this is your best option
 * since these routines should not be a bottleneck
 * in your computation);
 * or, use your own BLAS library, such 
 * as Intel MKL, ATLAS, GotoBLAS, etc.
 * (For example, you could use libmwblas
 * or libmwblas_compat32, which are included
 * with Mathworks and usually based on Intel MKL).
 *
 * The reason not to use your own BLAS library
 * is that (1) it may be faster but it won't be
 * a bottleneck, and (2) some BLAS libraries
 * use 32-bit/4byte integers, and others use 64-bit/
 * 8byte integers, and you pass by reference,
 * so if you get it wrong, you crash.
 * In Matlab, useing -lmwblascompat32
 *   uses names like ddot32
 *
 * In short, use the Ref version unless you feel
 * lucky.
 * */
#if !defined( _USE_OPTIMIZED_BLAS )
/* make alias to our reference implementation */
#define daxpy FORTRAN_WRAPPER(daxpyRef)
#define dcopy FORTRAN_WRAPPER(dcopyRef)
#define ddot  FORTRAN_WRAPPER(ddotRef)
#define dscal FORTRAN_WRAPPER(dscalRef)
#else
#if defined( _USE_32_BIT_BLAS )
#define daxpy FORTRAN_WRAPPER(daxpy32)
#define dcopy FORTRAN_WRAPPER(dcopy32)
#define ddot  FORTRAN_WRAPPER(ddot32)
#define dscal FORTRAN_WRAPPER(dscal32)
#else
#define daxpy FORTRAN_WRAPPER(daxpy)
#define dcopy FORTRAN_WRAPPER(dcopy)
#define ddot  FORTRAN_WRAPPER(ddot)
#define dscal FORTRAN_WRAPPER(dscal)
#endif
#endif


/* Fortran style i/o file was via pointer to integer
 * but with C, want pointer to FILE */
typedef FILE* fileType; 
/* typedef integer* fileType; */


/* Some linesearch parameters. The default
 * values are the same as the FORTRAN code,
 * but you could conceivably change these
 * at compile time; they are used in
 * dcsrch()  */
#ifndef FTOL
#define FTOL .001
#endif
#ifndef GTOL
#define GTOL .9
#endif
#ifndef XTOL
#define XTOL .1
#endif
#ifndef STEPMIN
#define STEPMIN 0.
#endif


/* If we want machine precision in a nice fashion, do this: */
#include <float.h>
#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2e-16
#endif







#ifdef __cplusplus
    extern "C" {
#endif



/* math.h */
double sqrt(double);

extern double ddot(integer *, double *, integer *, double *, 
        integer *);

extern  int daxpy(integer *, double *, double *, 
        integer *, double *, integer *);

extern  int dscal(integer *, double *, double *, 
        integer *);

extern  int dcopy(integer *, double *, integer *, 
	    double *, integer *);

#define setulb FORTRAN_WRAPPER(setulb)
extern int setulb(integer *n, integer *m, double *x, 
	double *l, double *u, integer *nbd, double *f, double 
	*g, double *factr, double *pgtol, double *wa, integer *
	iwa, integer *task, integer *iprint, integer *csave, logical *lsave, 
	integer *isave, double *dsave,
    double      (*valgrad) (double *, double *, integer) /* f = valgrad (g,x,n)*/
    );

#define mainlb FORTRAN_WRAPPER(mainlb)
extern int mainlb(integer *n, integer *m, double *x, 
        double *l, double *u, integer *nbd, double *f, double 
        *g, double *factr, double *pgtol, double *ws, double *
        wy, double *sy, double *ss, double *wt, double *wn, 
        double *snd, double *z__, double *r__, double *d__, 
        double *t, double *xp, double *wa, integer *index, 
        integer *iwhere, integer *indx2, integer *task, integer *iprint, 
        integer *csave, logical *lsave, integer *isave, double *dsave,
        double      (*valgrad) (double *, double *, integer) /* f = valgrad (g,x,n)*/
        );


#define freev FORTRAN_WRAPPER(freev) 
extern  int freev(integer *, integer *, integer *, 
        integer *, integer *, integer *, integer *, logical *, logical *, 
        logical *, integer *, integer *);

#define timer FORTRAN_WRAPPER(timer) 
extern  int timer(double *);
#define formk FORTRAN_WRAPPER(formk) 
extern int formk(integer *, 
        integer *, integer *, integer *, integer *, integer *, integer *, 
        logical *, double *, double *, integer *, double *, 
        double *, double *, double *, integer *, integer *, 
        integer *);
#define formt FORTRAN_WRAPPER(formt) 
extern  int formt(integer *, double *, double *, 
        double *, integer *, double *, integer *);
#define subsm FORTRAN_WRAPPER(subsm) 
extern int subsm(integer *, integer *, integer *, integer *, double *, double *, 
        integer *, double *, double *, double *, double *,
        double *, double *, double *, double *, integer *
        , integer *, integer *, double *, double *, integer *, 
        integer *);
#define prn1lb FORTRAN_WRAPPER(prn1lb) 
extern int prn1lb(integer *n, integer *m, double *l, 
        double *u, double *x, integer *iprint, fileType itfile, double *epsmch); 

#define prn2lb FORTRAN_WRAPPER(prn2lb) 
extern int prn2lb(integer *n, double *x, double *f, double *g, 
        integer *iprint, fileType itfile, integer *iter, integer *nfgv, integer *nact, double 
        * sbgnrm, integer *nseg, integer * word, integer * iword, integer * iback, double * stp, 
        double * xstep, ftnlen);

#define prn3lb FORTRAN_WRAPPER(prn3lb) 
extern int prn3lb(integer *n, double *x, double *f, integer *
	task, integer *iprint, integer *info, fileType itfile, integer *iter, 
	integer *nfgv, integer *nintol, integer *nskip, integer *nact, 
	double *sbgnrm, double *time, integer *nseg, integer *word, 
	integer *iback, double *stp, double *xstep, integer *k, 
	double *cachyt, double *sbtime, double *lnscht, ftnlen 
	task_len, ftnlen word_len);

#define errclb FORTRAN_WRAPPER(errclb) 
extern int errclb(integer *n, integer *m, double *factr, 
        double *l, double *u, integer *nbd, integer *task, integer *info,
        integer *k, ftnlen task_len);

#define active FORTRAN_WRAPPER(active) 
extern  int active(integer *, double *, double *,
        integer *, double *, integer *, integer *, logical *, 
        logical *, logical *);
#define cauchy FORTRAN_WRAPPER(cauchy) 
extern int cauchy(integer *, double *, 
        double *, double *, integer *, double *, integer *, 
        integer *, double *, double *, double *, integer *, 
        double *, double *, double *, double *, 
        double *, integer *, integer *, double *, double *, 
        double *, double *, integer *, integer *, double *, 
        integer *, double *);
#define cmprlb FORTRAN_WRAPPER(cmprlb) 
extern  int cmprlb(integer *, integer *, double *, 
        double *, double *, double *, double *, 
        double *, double *, double *, double *, integer *,
        double *, integer *, integer *, integer *, logical *, 
        integer *);
#define matupd FORTRAN_WRAPPER(matupd) 
extern  int matupd(integer *, integer *, double *, 
        double *, double *, double *, double *, 
        double *, integer *, integer *, integer *, integer *, 
        double *, double *, double *, double *, 
        double *);
#define lnsrlb FORTRAN_WRAPPER(lnsrlb) 
extern int lnsrlb(integer *n, double *l, double *u, 
        integer *nbd, double *x, double *f, double *fold, 
        double *gd, double *gdold, double *g, double *d__, 
        double *r__, double *t, double *z__, double *stp, 
        double *dnorm, double *dtd, double *xstep, double *
        stpmx, integer *iter, integer *ifun, integer *iback, integer *nfgv, 
        integer *info, integer *task, logical *boxed, logical *cnstnd, integer *
        csave, integer *isave, integer *iprint, double *dsave,
        integer max_ls,
        double      (*valgrad) (double *, double *, integer) );

#define projgr FORTRAN_WRAPPER(projgr) 
extern  int projgr(integer *, double *, double *,
        integer *, double *, double *, double *);

/* in linesearch.c */
#define dcsrch FORTRAN_WRAPPER(dcsrch) 
extern int dcsrch(double *f, double *g, double *stp, 
        double *ftol, double *gtol, double *xtol, double *
        stpmin, double *stpmax, integer *task, integer *isave, double *
        dsave);/* ftnlen task_len);*/
#define dcstep FORTRAN_WRAPPER(dcstep) 
extern  int dcstep(double *, double *,
        double *, double *, double *, double *,
        double *, double *, double *, logical *, double *,
        double *);

#define dpofa FORTRAN_WRAPPER(dpofa)
extern int dpofa(double *, integer *, integer *, 
		integer *);
#define dtrsl FORTRAN_WRAPPER(dtrsl)
extern int  dtrsl(double *, integer *, integer *, 
		double *, integer *, integer *);

double evaluate(double*g, double *x, integer n );


#define INT long int
#define INT_INF LONG_MAX
#define INF DBL_MAX

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef NULL
#define NULL 0
#endif
/*============================================================================
   cg_parameter is a structure containing parameters used in cg_descent
   cg_default assigns default values to these parameters */
typedef struct cg_parameter_struct /* user controlled parameters */
{
/*============================================================================
      parameters that the user may wish to modify
  ----------------------------------------------------------------------------*/
  /* T => print final statistics
     F => no printout of statistics */
  int PrintFinal ;

  /* Level 0  = no printing), ... , Level 3 = maximum printing */
  int PrintLevel ;

  /* T => print parameters values
     F => do not display parmeter values */
  int PrintParms ;

  /* T => use LBFGS
     F => only use L-BFGS when memory >= n */
  int LBFGS ;

  /* number of vectors stored in memory */
  int memory ;

  /* SubCheck and SubSkip control the frequency with which the subspace
     condition is checked. It is checked for SubCheck*mem iterations and
     if not satisfied, then it is skipped for Subskip*mem iterations
     and Subskip is doubled. Whenever the subspace condition is statisfied,
     SubSkip is returned to its original value. */
  int SubCheck ;
  int SubSkip ;

  /* when relative distance from current gradient to subspace <= eta0,
     enter subspace if subspace dimension = mem */
  double eta0 ;

  /* when relative distance from current gradient to subspace >= eta1,
     leave subspace */
  double eta1 ;

  /* when relative distance from current direction to subspace <= eta2,
     always enter subspace (invariant space) */
  double eta2 ;

  /* T => use approximate Wolfe line search
     F => use ordinary Wolfe line search, switch to approximate Wolfe when
              |f_k+1-f_k| < AWolfeFac*C_k, C_k = average size of cost  */
  int    AWolfe ;
  double AWolfeFac ;

  /* factor in [0, 1] used to compute average cost magnitude C_k as follows:
     Q_k = 1 + (Qdecay)Q_k-1, Q_0 = 0,  C_k = C_k-1 + (|f_k| - C_k-1)/Q_k */
  double Qdecay ;

  /* terminate after nslow iterations without strict improvement in
     either function value or gradient */
  int nslow ;

  /* Stop Rules:
     T => ||proj_grad||_infty <= max(grad_tol,initial ||grad||_infty*StopFact)
     F => ||proj_grad||_infty <= grad_tol*(1 + |f_k|) */
  int    StopRule ;
  double StopFac ;

  /* T => estimated error in function value is eps*Ck,
     F => estimated error in function value is eps */
  int    PertRule ;
  double eps ;

  /* factor by which eps grows when line search fails during contraction */
  double egrow ;

  /* T => attempt quadratic interpolation in line search when
              |f_k+1 - f_k|/f_k <= QuadCutoff
     F => no quadratic interpolation step */
  int    QuadStep ;
  double QuadCutOff ;

  /* maximum factor by which a quad step can reduce the step size */
  double QuadSafe ;

  /* T => when possible, use a cubic step in the line search */
  int UseCubic ;

  /* use cubic step when |f_k+1 - f_k|/|f_k| > CubicCutOff */
  double CubicCutOff ;

  /* |f| < SmallCost*starting cost => skip QuadStep and set PertRule = FALSE*/
  double SmallCost ;

  /* T => check that f_k+1 - f_k <= debugtol*C_k
     F => no checking of function values */
  int    debug ;
  double debugtol ;

  /* if step is nonzero, it is the initial step of the initial line search */
  double step ;

  /* abort cg after maxit iterations */
  INT maxit ;

  /* maximum number of times the bracketing interval grows during expansion */
  int ntries ;

  /* maximum factor secant step increases stepsize in expansion phase */
  double ExpandSafe ;

  /* factor by which secant step is amplified during expansion phase
     where minimizer is bracketed */
  double SecantAmp ;

  /* factor by which rho grows during expansion phase where minimizer is
     bracketed */
  double RhoGrow ;

  /* maximum number of times that eps is updated */
  int neps ;

  /* maximum number of times the bracketing interval shrinks */
  int nshrink ;

  /* maximum number of iterations in line search */
  int nline ;

  /* conjugate gradient method restarts after (n*restart_fac) iterations */
  double restart_fac ;

  /* stop when -alpha*dphi0 (estimated change in function value) <= feps*|f|*/
  double feps ;

  /* after encountering nan, growth factor when searching for
     a bracketing interval */
  double nan_rho ;

  /* after encountering nan, decay factor for stepsize */
  double nan_decay ;

/*============================================================================
       technical parameters which the user probably should not touch
  ----------------------------------------------------------------------------*/
  double           delta ; /* Wolfe line search parameter */
  double           sigma ; /* Wolfe line search parameter */
  double           gamma ; /* decay factor for bracket interval width */
  double             rho ; /* growth factor when searching for initial
                                bracketing interval */
  double            psi0 ; /* factor used in starting guess for iteration 1 */
  double          psi_lo ; /* in performing a QuadStep, we evaluate at point
                                betweeen [psi_lo, psi_hi]*psi2*previous step */
  double          psi_hi ;
  double            psi1 ; /* for approximate quadratic, use gradient at
                                psi1*psi2*previous step for initial stepsize */
  double            psi2 ; /* when starting a new cg iteration, our initial
                                guess for the line search stepsize is
                                psi2*previous step */
  int       AdaptiveBeta ; /* T => choose beta adaptively, F => use theta */
  double       BetaLower ; /* lower bound factor for beta */
  double           theta ; /* parameter describing the cg_descent family */
  double            qeps ; /* parameter in cost error for quadratic restart
                                criterion */
  double           qrule ; /* parameter used to decide if cost is quadratic */
  int           qrestart ; /* number of iterations the function should be
                                nearly quadratic before a restart */
} cg_parameter ;

typedef struct cg_stats_struct /* statistics returned to user */
{
  double               f ; /*function value at solution */
  double           gnorm ; /* max abs component of gradient */
  INT               iter ; /* number of iterations */
  INT            IterSub ; /* number of subspace iterations */
  INT             NumSub ; /* total number subspaces */
  INT              nfunc ; /* number of function evaluations */
  INT              ngrad ; /* number of gradient evaluations */
} cg_stats ;

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
  double         rho ; /* either Parm->rho or Parm->nan_rho */
  double    alphaold ; /* previous value for stepsize alpha */
  double          *x ; /* current iterate */
  double      *xtemp ; /* x + alpha*d */
  double          *d ; /* current search direction */
  double          *g ; /* gradient at x */
  double      *gtemp ; /* gradient at x + alpha*d */
  double   (*cg_value) (double *, INT) ; /* f = cg_value (x, n) */
  void      (*cg_grad) (double *, double *, INT) ; /* cg_grad (g, x, n) */
  double (*cg_valgrad) (double *, double *, INT) ; /* f = cg_valgrad (g,x,n)*/
  cg_parameter *Parm ; /* user parameters */
} cg_com ;







#ifdef __cplusplus
    }   /* extern "C" */
#endif

#endif /* lbfgsb_h */
