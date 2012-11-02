#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define grsr1 0
#define grsr2 1

extern double rho0;

double xsr1_iglad(double k, double h, double zr, double mu0, double mu1, int N);

void xgrsr_ad(int *ctyp, double *k, double *h, double *mu0, double *mu1, double *zr, int *r, double *ad)
{ 
 if (*ctyp==grsr1) *ad = xsr1_iglad(*k, *h, *zr, *mu0, *mu1, *r); 
} 
