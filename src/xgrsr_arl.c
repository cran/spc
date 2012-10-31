#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define grsr1 0
#define grsr2 1

extern double rho0;

double xsr1_iglarl(double k, double h, double zr, double hs, double mu, int N);
double xsr1_arlm(double k, double h, double zr, double hs, int q, double mu0, double mu1, int N, int nmax);

void xgrsr_arl(int *ctyp, double *k, double *h, double *zr, double *hs, double *mu, int *q, int *r, double *arl)
{ int nmax=100000;
 if (*ctyp==grsr1 && *q==1) *arl = xsr1_iglarl(*k, *h, *zr, *hs, *mu, *r); 
 if (*ctyp==grsr1 && *q>1)  *arl = xsr1_arlm(*k, *h, *zr, *hs, *q, 0., *mu, *r, nmax);
}
