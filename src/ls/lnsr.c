#include "../lbfgsb.h"
static integer c__1 = 1;


int lnsrlb(integer *n, double *l, double *u,
           integer *nbd, double *x, double *f, double *fold,
           double *gd, double *gdold, double *g, double *d__,
           double *r__, double *t, double *z__, double *stp,
           double *dnorm, double *dtd, double *xstep, double *
stpmx, integer *iter, integer *ifun, integer *iback, integer *nfgv,
           integer *info, integer *task, logical *boxed, logical *cnstnd, integer *
ctask, integer *isave, integer *iprint, double *dsave,
           integer max_ls,double      (*eval) (double *, double *, integer)
){
  /*
  **********

  Subroutine lnsrlb

  This subroutine calls subroutine dcsrch from the Minpack2 library
    to perform the line search.  Subroutine dscrch is safeguarded so
    that all trial points lie within the feasible region.

  Subprograms called:

    Minpack2 Library ... dcsrch.

    Linpack ... dtrsl, ddot.


                        *  *  *

  NEOS, November 1994. (Latest revision June 1996.)
  Optimization Technology Center.
  Argonne National Laboratory and Northwestern University.
  Written by
                     Ciyou Zhu
  in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.


  **********
  */


  /* Table of constant values */
  static double c_b14 = FTOL;
  static double c_b15 = GTOL;
  static double c_b16 = XTOL;
  static double c_b17 = STEPMIN;
  /* System generated locals */
  integer i__1;
  double d__1;


  /* Local variables */
  static integer i__;
  static double a1, a2;

  /* Parameter adjustments */
  --z__;
  --t;
  --r__;
  --d__;
  --g;
  --x;
  --nbd;
  --u;
  --l;
  --isave;
  --dsave;


  *dtd = ddot(n, &d__[1], &c__1, &d__[1], &c__1);
  *dnorm = sqrt(*dtd);
  /*     Determine the maximum step length. */
  *stpmx = 1e10;
  if (*cnstnd) {
    if (*iter == 0) {
      *stpmx = 1.;
    } else {
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        a1 = d__[i__];
        if (nbd[i__] != 0) {
          if (a1 < 0. && nbd[i__] <= 2) {
            a2 = l[i__] - x[i__];
            if (a2 >= 0.) {
              *stpmx = 0.;
            } else if (a1 * *stpmx < a2) {
              *stpmx = a2 / a1;
            }
          } else if (a1 > 0. && nbd[i__] >= 2) {
            a2 = u[i__] - x[i__];
            if (a2 <= 0.) {
              *stpmx = 0.;
            } else if (a1 * *stpmx > a2) {
              *stpmx = a2 / a1;
            }
          }
        }
        /* L43: */
      }
    }
  }
  if (*iter == 0 && ! (*boxed)) {
    /* Computing MIN */
    d__1 = 1. / *dnorm;
    *stp = fmin(d__1,*stpmx);
  } else {
    *stp = 1.;
  }
  dcopy(n, &x[1], &c__1, &t[1], &c__1);
  dcopy(n, &g[1], &c__1, &r__[1], &c__1);
  *fold = *f;
  *ifun = 0;
  *iback = 0;
  *ctask = START;


  //L556:
  /* Function Body */
 integer ls_iter = 0;
 while (ls_iter < max_ls ) { //
   ls_iter++;
   *gd = ddot(n, &g[1], &c__1, &d__[1], &c__1);
   if (*ifun == 0) {
     *gdold = *gd;
     if (*gd >= 0.) {
       if (*iprint >= 0) {
         printf("ascend direction in projection gd = %.2e\n", *gd);
       }
       *info = -4;
       return 0;
     }
   }
   //printf("            initial stp = %5e \n", *stp);
   dcsrch(f, gd, stp, &c_b14, &c_b15, &c_b16, &c_b17, stpmx, ctask, &isave[1], &dsave[1]);/* (ftnlen)60);*/
   *xstep = *stp * *dnorm;
   //printf("inside ls, task = %5ld \n", *task);
   //printf("inside ls, ctask = %5ld \n", ctask);

   if (!(IS_WARNING(*ctask)) && !(IS_CONVERGED(*ctask))) {
     *task = FG_LNSRCH;
     ++(*ifun);
     ++(*nfgv);
     *iback = *ifun - 1;
     if (*stp == 1.) {
       dcopy(n, &z__[1], &c__1, &x[1], &c__1);
     } else {
       i__1 = *n;
       for (i__ = 1; i__ <= i__1; ++i__) {
         x[i__] = *stp * d__[i__] + t[i__];
         /* L41: */
       }
     }
     *f = eval(&g[1], &x[1], *n);
   } else {
     //printf("nfgv = %5ld \n", *nfgv);
     //printf("final stp = %5e \n", *stp);
     *task = NEW_X;
     return 0;
   }

 }
  return 0;
} /* lnsrlb */

/* ======================= The end of lnsrlb ============================= */