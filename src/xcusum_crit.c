#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

extern double rho0;

double xc_crit(int ctyp, double k, double L0, double hs, double m0, int N);

void xcusum_crit
( int *ctyp, double *k, double *L0, double *hs, double *mu0, int *r, double *h)
{
  *h = xc_crit(*ctyp,*k,*L0,*hs,*mu0,*r);
}
