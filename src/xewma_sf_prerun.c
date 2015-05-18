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

#define MU 0
#define SIGMA 1
#define BOTH 2

double *vector (long n);

double xe2_sf_prerun_MU_deluxe(double l, double c, double hs, double mu, int pn, int nmax, int qm, double truncate, double BOUND, double *p0);
double xe2_sf_prerun_MU(double l, double c, double hs, double mu, int pn, int nmax, int qm, double truncate, double *p0);

double xe2_sf_prerun_SIGMA_deluxe(double l, double c, double hs, double mu, int pn, int nmax, int qm, double truncate, double BOUND, double *p0);
double xe2_sf_prerun_SIGMA(double l, double c, double hs, double mu, int pn, int nmax, int qm, double truncate, double *p0);

double xe2_sf_prerun_BOTH_deluxe(double l, double c, double hs, double mu, int pn, int df, int nmax, int qm1, int qm2, double truncate, double BOUND, double *p0);
double xe2_sf_prerun_BOTH(double l, double c, double hs, double mu, int pn, int df, int nmax, int qm1, int qm2, double truncate, double *p0);


double xe2_sfm_prerun_MU_deluxe(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate, double BOUND, double *p0);
double xe2_sfm_prerun_MU(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate, double *p0);

double xe2_sfm_prerun_SIGMA_deluxe(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate, double BOUND, double *p0);
double xe2_sfm_prerun_SIGMA(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate, double *p0);

double xe2_sfm_prerun_BOTH_deluxe(double l, double c, double hs, int q, double mu0, double mu1, int pn, int df, int mode, int nmax, int qm1, int qm2, double truncate, double BOUND, double *p0);
double xe2_sfm_prerun_BOTH(double l, double c, double hs, int q, double mu0, double mu1, int pn, int df, int mode, int nmax, int qm1, int qm2, double truncate, double *p0);


void xewma_sf_prerun
( int *ctyp, double *l, double *c, double *zr, double *hs, double *mu, int *ltyp, int *q, int *n, int *size, int *df, int *mode,
  int *qm1, int *qm2, double *truncate, int *tail_approx, double *bound, double *sf)

{ int i, result=0;
  double *p0;
 p0  = vector(*n); 

 if ( *mode == MU ) {
   if ( *ctyp==ewma2 && *ltyp==fix && *q==1 ) {
     if ( *tail_approx ) result = xe2_sf_prerun_MU_deluxe(*l, *c, *hs, *mu, *size, *n, *qm1, *truncate, *bound, p0);
     else result =  xe2_sf_prerun_MU(*l, *c, *hs, *mu, *size, *n, *qm1, *truncate, p0);
   }
   if ( *ctyp==ewma2 && *ltyp==fix && *q>1 )  {
     if ( *tail_approx ) result = xe2_sfm_prerun_MU_deluxe(*l, *c, *hs, *q, 0., *mu, *size, *ltyp, *n, *qm1, *truncate, *bound, p0);
     else result = xe2_sfm_prerun_MU(*l, *c, *hs, *q, 0., *mu, *size, *ltyp, *n, *qm1, *truncate, p0);
   }
   if ( *ctyp==ewma2 && *ltyp>fix ) {
     if ( *tail_approx ) result = xe2_sfm_prerun_MU_deluxe(*l, *c, *hs, *q, 0., *mu, *size, *ltyp, *n, *qm1, *truncate, *bound, p0);
     else result = xe2_sfm_prerun_MU(*l, *c, *hs, *q, 0., *mu, *size, *ltyp, *n, *qm1, *truncate, p0);
   }
 }
 
 if ( *mode == SIGMA ) {
   if ( *ctyp==ewma2 && *ltyp==fix && *q==1 ) {
     if ( *tail_approx ) result = xe2_sf_prerun_SIGMA_deluxe(*l, *c, *hs, *mu, *size, *n, *qm2, *truncate, *bound, p0);
     else result =  xe2_sf_prerun_SIGMA(*l, *c, *hs, *mu, *size, *n, *qm2, *truncate, p0);
   }
   if ( *ctyp==ewma2 && *ltyp==fix && *q>1 )  {
     if ( *tail_approx ) result = xe2_sfm_prerun_SIGMA_deluxe(*l, *c, *hs, *q, 0., *mu, *size, *ltyp, *n, *qm2, *truncate, *bound, p0);
     else result = xe2_sfm_prerun_SIGMA(*l, *c, *hs, *q, 0., *mu, *size, *ltyp, *n, *qm2, *truncate, p0);
   }
   if ( *ctyp==ewma2 && *ltyp>fix ) {
     if ( *tail_approx ) result = xe2_sfm_prerun_SIGMA_deluxe(*l, *c, *hs, *q, 0., *mu, *size, *ltyp, *n, *qm2, *truncate, *bound, p0);
     else result = xe2_sfm_prerun_SIGMA(*l, *c, *hs, *q, 0., *mu, *size, *ltyp, *n, *qm2, *truncate, p0);
   }
 }
 
 if ( *mode == BOTH ) {
   if ( *ctyp==ewma2 && *ltyp==fix && *q==1 ) {
     if ( *tail_approx ) result = xe2_sf_prerun_BOTH_deluxe(*l, *c, *hs, *mu, *size, *df, *n, *qm1, *qm2, *truncate, *bound, p0);
     else result =  xe2_sf_prerun_BOTH(*l, *c, *hs, *mu, *size, *df, *n, *qm1, *qm2, *truncate, p0);     
   }
   if ( *ctyp==ewma2 && *ltyp==fix && *q>1 ) {
     if ( *tail_approx ) result = xe2_sfm_prerun_BOTH_deluxe(*l, *c, *hs, *q, 0., *mu, *size, *df, *ltyp, *n, *qm1, *qm2, *truncate, *bound, p0);
     else result = xe2_sfm_prerun_BOTH(*l, *c, *hs, *q, 0., *mu, *size, *df, *ltyp, *n, *qm1, *qm2, *truncate, p0);
   }
   if ( *ctyp==ewma2 && *ltyp>fix ) {
     if ( *tail_approx ) result = xe2_sfm_prerun_BOTH_deluxe(*l, *c, *hs, *q, 0., *mu, *size, *df, *ltyp, *n, *qm1, *qm2, *truncate, *bound, p0);
     else result = xe2_sfm_prerun_BOTH(*l, *c, *hs, *q, 0., *mu, *size, *df, *ltyp, *n, *qm1, *qm2, *truncate, p0);
   }
 }
 
 if ( result != 0 ) warning("\nSomething bad happened!\n\n");
 
 for (i=0; i<*n; i++) sf[i] = p0[i];
}
