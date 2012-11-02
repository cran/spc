#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define cusum1 0
#define cusum2 1
#define cusumC 2

#define igl 0
#define mc 1

extern double rho0;

double xc1_iglarl(double k, double h, double hs, double mu, int N);
double xc1_arlm(double k, double h, double hs, int q, double mu0, double mu1, int N, int nmax);
double xc2_iglarl(double k, double h, double hs, double mu, int N);
double xc2_be_arl(double k, double h, double hs1, double hs2, double mu, int N);
double xcC_iglarl(double k, double h, double hs, double mu, int N);

void xcusum_arl
( int *ctyp, double *k, double *h, double *hs, double *mu, int *q, int *r, int *method, double *arl)
{ int nmax=100000;
  double lhs;

 if ( *ctyp == cusum1 && *q==1 ) *arl = xc1_iglarl(*k,*h,*hs,*mu,*r);
 if ( *ctyp == cusum1 && *q>1 )  *arl = xc1_arlm(*k, *h, *hs, *q, 0., *mu, *r, nmax);
 if ( *ctyp == cusum2 ) {
   if ( *method == igl ) *arl = xc2_iglarl(*k,*h,*hs,*mu,*r);
   lhs = - *hs;
   if ( *method == mc )  *arl = xc2_be_arl(*k,*h,*hs,lhs,*mu,*r); 
 }
 if ( *ctyp == cusumC ) *arl = xcC_iglarl(*k,*h,*hs,*mu,*r);
}
