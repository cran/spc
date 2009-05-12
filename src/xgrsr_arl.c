#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define grsr1 0
#define grsr2 1

extern double rho0;

double xsr1_iglarl(double k, double h, double zr, double hs, double mu, int N);

void xgrsr_arl
( int *ctyp, double *k, double *h, double *zr, double *hs, double *mu, int *r, double *arl)
{
 if (*ctyp==grsr1) *arl = xsr1_iglarl(*k, *h, *zr, *hs, *mu, *r); 
}
