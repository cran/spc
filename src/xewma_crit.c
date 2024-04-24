#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define cfar 6

double *vector (long n);
extern double rho0;

double xe_crit(int ctyp, double l, double L0, double zr, double hs, double m0, int ltyp, int N, double c0);
int xe2fr_crit(double l, double L0, double cinf, int N, int nmax, double *cn);

void xewma_crit(int *ctyp, double *l, double *L0, double *zr, double *hs, double *mu0, int *ltyp, int *r, double *c0, double *h)
{ int i, n;
  double *cn;
  cn  = vector(*ctyp);
  if ( *ltyp == cfar ) {
    n = xe2fr_crit(*l, *L0, *hs, *r, *ctyp, cn); /* hs = cinf, h = cn, ctyp = nmax */
    for (i=0; i<*ctyp; i++) h[i] = cn[i];
  } else {
    *h = xe_crit(*ctyp, *l, *L0, *zr, *hs, *mu0, *ltyp, *r, *c0);
  }
}
