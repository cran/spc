#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewmaU 0
#define ewmaUR 1
#define ewma2 2
#define ewmaLR 3

double *vector (long n);

double seU_sf_prerun_SIGMA_deluxe(double l, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0);
double seU_sf_prerun_SIGMA(double l, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0);

double seUR_sf_prerun_SIGMA_deluxe(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0);
double seUR_sf_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0);

double se2_sf_prerun_SIGMA_deluxe(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0);
double se2_sf_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0);

double seLR_sf_prerun_SIGMA_deluxe(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0);
double seLR_sf_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0);

void sewma_sf_prerun
( int *ctyp, double *l, double *cl, double *cu, double *hs, double *sigma, int *df1, int *qm1, int *n,
      int *df2, int *qm2, double *truncate, int *tail_approx, double *sf)
{ int result=0, i;
  double *p0;
 p0  = vector(*n);

 if ( *ctyp == ewmaU ) {
   if ( *tail_approx ) result = seU_sf_prerun_SIGMA_deluxe(*l, *cu, *hs, *sigma, *df1, *df2, *n, *qm1, *qm2, *truncate, p0);
   else result = seU_sf_prerun_SIGMA(*l, *cu, *hs, *sigma, *df1, *df2, *n, *qm1, *qm2, *truncate, p0);
 }
   
 if ( *ctyp == ewmaUR ) {
   if ( *tail_approx ) result = seUR_sf_prerun_SIGMA_deluxe(*l, *cl, *cu, *hs, *sigma, *df1, *df2, *n, *qm1, *qm2, *truncate, p0);
   else result = seUR_sf_prerun_SIGMA(*l, *cl, *cu, *hs, *sigma, *df1, *df2, *n, *qm1, *qm2, *truncate, p0);
 }
 
  if ( *ctyp == ewma2 ) {
   if ( *tail_approx ) result = se2_sf_prerun_SIGMA_deluxe(*l, *cl, *cu, *hs, *sigma, *df1, *df2, *n, *qm1, *qm2, *truncate, p0);
   else result = se2_sf_prerun_SIGMA(*l, *cl, *cu, *hs, *sigma, *df1, *df2, *n, *qm1, *qm2, *truncate, p0);
 }
 
  if ( *ctyp == ewmaLR ) {
   if ( *tail_approx ) result = seLR_sf_prerun_SIGMA_deluxe(*l, *cl, *cu, *hs, *sigma, *df1, *df2, *n, *qm1, *qm2, *truncate, p0);
   else result = seLR_sf_prerun_SIGMA(*l, *cl, *cu, *hs, *sigma, *df1, *df2, *n, *qm1, *qm2, *truncate, p0);
 }
 
 if ( result != 0 ) warning("trouble in sewma_sf_prerun [package spc]");

 for (i=0; i<*n; i++) sf[i] = p0[i];
}
