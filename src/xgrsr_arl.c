#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define grsr1 0
#define grsr2 1

extern double rho0;

double *vector (long n);
double xsr1_iglarl(double k, double h, double zr, double hs, double mu, int N, int MPT);
double xsr1_arlm(double k, double h, double zr, double hs, int q, double mu0, double mu1, int N, int nmax, int MPT);
double xsr1_arlm_hom(double k, double h, double zr, double hs, int q, double mu0, double mu1, int N, int MPT, double *ced);

void xgrsr_arl(int *ctyp, double *k, double *h, double *zr, double *hs, double *mu, int *q, int *r, int *MPT, double *arl)
{ int i, /*nmax=100000,*/ result=0;
  double *ced, arl1=-1.;
 ced  = vector(*q);
 
 if ( *ctyp==grsr1 && *q==1 ) arl1 = xsr1_iglarl(*k, *h, *zr, *hs, *mu, *r, *MPT); 
 if ( *ctyp==grsr1 && *q>1  ) result = xsr1_arlm_hom(*k, *h, *zr, *hs, *q, 0., *mu, *r, *MPT, ced);
 /* *arl = xsr1_arlm(*k, *h, *zr, *hs, *q, 0., *mu, *r, nmax, *MPT);*/ 
 
 if ( result != 0 ) warning("trouble in xgrsr_arl [package spc]");
 
 if (  *q > 1 ) for (i=0; i<*q; i++) arl[i] = ced[i]; else *arl = arl1;
}
