#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

extern double rho0;

double xe_crit(int ctyp, double l, double L0, double zr, double hs, 
               double m0, int ltyp, int N);

void xewma_crit
( int *ctyp, double *l, double *L0, double *zr, double *hs, 
  double *mu0, int *ltyp, int *r, double *h)
{
  *h = xe_crit(*ctyp,*l,*L0,*zr,*hs,*mu0,*ltyp,*r);
}
