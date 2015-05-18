#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewma1 0
#define ewma2 1
#define fix 0
#define vacl 1

double *vector (long n);

/*
double xe1_sf (double l, double c, double zr, double hs, double mu, int N, int nmax, double *p0);
double xe1_sfm(double l, double c, double zr, double hs, int q, double mu0, double mu1, int mode, int N, int nmax, double *p0);
*/

double xte2_sf(double l, double c, double hs, int df, double mu, int N, int nmax, double *p0, int subst);

double xte2_sfm(double l, double c, double hs, int df, int q, double mu0, double mu1, int mode, int N, int nmax, double *p0, int subst);


void xtewma_sf(int *ctyp, double *l, double *c, double *zr, double *hs, int *df, double *mu, int *ltyp, int *r, int *ntyp, int *q, int *n, double *sf)
{ int result=0, i;
  double *p0;
 p0  = vector(*n);
 
 if ( *ctyp==ewma2 && *ltyp==fix && *q==1 ) {
   result = xte2_sf(*l, *c, *hs, *df, *mu, *r, *n, p0, *ntyp);
 }  
 
 if ( *ctyp==ewma2 && ( ( *ltyp==fix && *q>1 ) || ( *ltyp>fix) ) ) {
   result = xte2_sfm(*l, *c, *hs, *df, *q, 0., *mu, *ltyp, *r, *n, p0, *ntyp);
 }

 if ( result != 0 ) warning("trouble in xtewma_sf [package spc]");

 for (i=0; i<*n; i++) sf[i] = p0[i];
}
