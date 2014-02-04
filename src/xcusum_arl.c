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

double *vector (long n);
double xc1_iglarl(double k, double h, double hs, double mu, int N);
double xc1_arlm(double k, double h, double hs, int q, double mu0, double mu1, int N, int nmax);
double xc1_arlm_hom(double k, double h, double hs, int q, double mu0, double mu1, int N, double *ced);
double xc2_iglarl(double k, double h, double hs, double mu, int N);
double xc2_be_arl(double k, double h, double hs1, double hs2, double mu, int N);
double xcC_iglarl(double k, double h, double hs, double mu, int N);

void xcusum_arl
( int *ctyp, double *k, double *h, double *hs, double *mu, int *q, int *r, int *method, double *arl)
{ int i, /*nmax=100000,*/ result=0;
  double lhs, *ced, arl1=-1.;
 ced  = vector(*q);

 if ( *ctyp == cusum1 && *q==1 ) arl1 = xc1_iglarl(*k,*h,*hs,*mu,*r);
 if ( *ctyp == cusum1 && *q>1 ) result = xc1_arlm_hom(*k, *h, *hs, *q, 0., *mu, *r, ced);
 /* *arl = xc1_arlm(*k, *h, *hs, *q, 0., *mu, *r, nmax); */
 if ( *ctyp == cusum2 ) {
   if ( *method == igl ) arl1 = xc2_iglarl(*k,*h,*hs,*mu,*r);
   lhs = - *hs;
   if ( *method == mc )  arl1 = xc2_be_arl(*k,*h,*hs,lhs,*mu,*r); 
 }
 if ( *ctyp == cusumC ) arl1 = xcC_iglarl(*k,*h,*hs,*mu,*r);
 
 if ( result != 0 ) warning("trouble in xgrsr_arl [package spc]");
 
 if (  *q > 1 ) for (i=0; i<*q; i++) arl[i] = ced[i]; else *arl = arl1;
}
