#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define cusum1 0
#define cusum2 1
#define cusumC 2

extern double rho0;

double xc1_iglad (double k, double h, double mu0, double mu1, int N);
double xc2_iglad (double k, double h, double mu0, double mu1, int N);
double xc2_igladR(double k, double h, double mu0, double mu1, int r);
double xcC_iglad (double k, double h, double mu0, double mu1, int N);

void xcusum_ad
( int *ctyp, double *k, double *h, double *mu0, double *mu1,
  int *r, double *ad)
{ 
 if (*ctyp==cusum1) *ad = xc1_iglad(*k,*h,*mu0,*mu1,*r);
 if (*ctyp==cusum2 && *r>0) *ad = xc2_iglad(*k,*h,*mu0,*mu1,*r);
 if (*ctyp==cusum2 && *r<0) *ad = xc2_igladR(*k,*h,*mu0,*mu1,-*r); 
 if (*ctyp==cusumC) *ad = xcC_iglad(*k,*h,*mu0,*mu1,*r);
} 
