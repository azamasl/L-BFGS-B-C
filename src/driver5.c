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


int main(void){
  /* System generated locals */
  integer i1;
  double d1, d2;
  /* Local variables */
  static double f, g[1024];
  static integer i;
  static double l[1024];
  static integer m, n;
  static double u[1024], x[1024], t1, t2;
  static integer nbd[1024];

  static double pgtol;
  static double factr;
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

  setulb(&n, &m, x, l, u, nbd, &f, g, &factr, &pgtol, &iprint,evaluate);
  return 0;
} /* MAIN__ */

