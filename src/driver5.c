#include "lbfgsb.h"

double evaluate( double    *g, double    *x, integer    n){
  double ex, f, t ;
  f = (double) 0 ;
  for (int i = 1; i <= n; i++){
    t = sqrt (i) ;
    ex = exp (x [i-1]) ;
    f += ex - t*x [i-1] ;
    g [i-1] = ex -  t ;
  }
  return f ;
}



/* Main program */
//int MAIN__(void)
int main(void){
  /* System generated locals */
  integer i1;
  double d1, d2;
  /* Local variables */
  static double f, g[1024];
  static integer i;
  static double l[1024];
  static integer m, n;
  static double u[1024], x[1024], t1, t2, wa[43251];
  static integer nbd[1024], iwa[3072];
/*     static char task[60]; */
  static integer taskValue;
  static integer *task=&taskValue; /* must initialize !! */
/*      http://stackoverflow.com/a/11278093/269192 */
  static double factr;
/*     static char csave[60]; */
  static integer csaveValue;
  static integer *csave=&csaveValue;
  static double dsave[29];
  static integer isave[44];
  static logical lsave[4];
  static double pgtol;
  static integer iprint;

/*     We wish to have output at every iteration. */
  iprint = 1;
/*     iprint = 101; */
/*     We specify the tolerances in the stopping criteria. */
  factr = 1;//1e7;
  pgtol = 1e-9;

  n = 100;
  m = 5;

  for (i = 1; i <= n; i++) {
    nbd[i - 1] = 0;
    x[i - 1] = 1.;
  }

  /*     We start the iteration by initializing task. */

  *task = (integer)START;
/*     s_copy(task, "START", (ftnlen)60, (ftnlen)5); */
  /*        ------- the beginning of the loop ---------- */
  //L111:
  /*     This is the call to the L-BFGS-B code. */
  setulb(&n, &m, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task, &iprint, csave, lsave, isave, dsave,evaluate);// evaluate);

//  if ( IS_FG(*task) ) {
//    //f = evaluate(g, x , n);
//    f = myvalgrad(g, x, n);
//    goto L111;
//  }
//
//  if ( *task==NEW_X ) {
//    goto L111;
//  }
  return 0;
} /* MAIN__ */

