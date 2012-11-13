#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define cusum1 0
#define cusum2 1

double *vector (long n);
double xc1_sf(double k, double h, double hs, double mu, int N, int nmax, double *p0);

void xcusum_sf(int *ctyp, double *k, double *h, double *hs, double *mu, int *r, int *n, double *sf)
{ int result=0, i;
  double *p0;
 p0  = vector(*n);
 if (*ctyp==cusum1) result = xc1_sf(*k, *h, *hs, *mu, *r, *n, p0);
 if ( result != 0 ) warning("trouble with xc1_sf called from xcusum_sf [package spc]");
 for (i=0; i<*n; i++) sf[i] = p0[i];
}
