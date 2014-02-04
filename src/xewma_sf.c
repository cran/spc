#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewma1 0
#define ewma2 1
#define fix 0
#define vacl 1
#define fir 2
#define both 3
#define steiner 4
#define stat 5
#define test 6

double *vector (long n);

double xe1_sf (double l, double c, double zr, double hs, double mu, int N, int nmax, double *p0);
double xe1_sfm(double l, double c, double zr, double hs, int q, double mu0, double mu1, int mode, int N, int nmax, double *p0);

double xe2_sf (double l, double c, double hs, double mu, int N, int nmax, double *p0);
double xe2_sfm(double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax, double *p0);

void xewma_sf(int *ctyp, double *l, double *c, double *zr, double *hs, double *mu, int *ltyp, int *r, int *q, int *n, double *sf)
{ int result=0, i;
  double *p0;
 p0  = vector(*n);
 
 if ( *ctyp==ewma1 && *ltyp==fix && *q==1 )  result = xe1_sf (*l, *c, *zr, *hs, *mu, *r, *n, p0);
 if ( *ctyp==ewma1 && *ltyp==fix && *q>1 )   result = xe1_sfm(*l, *c, *zr, *hs, *q, 0., *mu, *ltyp, *r, *n, p0);
 if ( *ctyp==ewma1 && *ltyp>fix )            result = xe1_sfm(*l, *c, *zr, *hs, *q, 0., *mu, *ltyp, *r, *n, p0);
 
 if ( *ctyp==ewma2 && *ltyp==fix && *q==1 )  result = xe2_sf (*l, *c, *hs, *mu, *r, *n, p0);
 if ( *ctyp==ewma2 && *ltyp==fix && *q>1 )   result = xe2_sfm(*l, *c, *hs, *q, 0., *mu, *ltyp, *r, *n, p0);
 if ( *ctyp==ewma2 && *ltyp>fix )            result = xe2_sfm(*l, *c, *hs, *q, 0., *mu, *ltyp, *r, *n, p0);   

 if ( result != 0 ) warning("trouble in xewma_sf [package spc]");

 for (i=0; i<*n; i++) sf[i] = p0[i];
}
