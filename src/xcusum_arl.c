#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define cusum1 0
#define cusum2 1
#define cusumC 2

extern double rho0;

double xc1_iglarl(double k, double h, double hs, double mu, int N);
double xc2_iglarl(double k, double h, double hs, double mu, int N);
double xcC_iglarl(double k, double h, double hs, double mu, int N);

void xcusum_arl
( int *ctyp, double *k, double *h, double *hs, double *mu, int *r, double *arl)
{
 if (*ctyp==cusum1) *arl = xc1_iglarl(*k,*h,*hs,*mu,*r);
 if (*ctyp==cusum2) *arl = xc2_iglarl(*k,*h,*hs,*mu,*r); 
 if (*ctyp==cusumC) *arl = xcC_iglarl(*k,*h,*hs,*mu,*r);
}
