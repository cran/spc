#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewmaU 0
#define ewmaUR 1
#define ewma2 2
#define ewmaLR 3

double *vector (long n);
double seU_sf(double l, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0);
double se2_sf(double l, double cl, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0);
double seUR_sf(double l, double cl, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0);
double seLR_sf(double l, double cl, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0);

void sewma_sf
( int *ctyp, 
      double *l, double *cl, double *cu, double *hs, int *N,
      double *sigma, 
      int *df, int *qm, int *n, double *sf)
{ int result=0, i;
  double *p0;
 p0  = vector(*n);

 if ( *ctyp == ewmaU )
   result = seU_sf(*l, *cu, *hs, *sigma, *df, *N, *n, *qm, p0);
 
 if ( *ctyp == ewmaUR )
   result = seUR_sf(*l, *cl, *cu, *hs, *sigma, *df, *N, *n, *qm, p0);
 
 if ( *ctyp == ewma2 )
   result = se2_sf(*l, *cl, *cu, *hs, *sigma, *df, *N, *n, *qm, p0);
 
 if ( *ctyp == ewmaLR )
   result = seLR_sf(*l, *cl, *cu, *hs, *sigma, *df, *N, *n, *qm, p0);

 if ( result != 0 ) warning("trouble in sewma_sf [package spc]");

 for (i=0; i<*n; i++) sf[i] = p0[i];
 
 Free(p0);
}
