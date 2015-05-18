#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewma1 0
#define ewma2 1
#define fix 0
#define vacl 1

#define conditional 0
#define cyclical 1

extern double rho0;

/*
double xe2_iglad (double l, double c, double mu0, double mu1, int N);
double xe2_igladc(double l, double c, double mu0, double mu1, double z0, int N);
double xe2_arlm  (double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax);*/

double xte2_iglad(double l, double c, int df, double mu0, double mu1, int N, int subst);
double xte2_igladc(double l, double c, int df, double mu0, double mu1, double z0, int N, int subst);
double xte2_arlm (double l, double c, double hs, int df, int q, double mu0, double mu1, int mode, int N, int nmax, int subst);

void xtewma_ad(int *ctyp, double *l, double *c, double *zr, int *df, double *mu0, double *mu1, double *z0, int *ltyp, int *styp, int *r, int *ntyp, double *ad)
{ int nmax=1000000;
 if ( *styp==conditional ) {
   if ( *ctyp==ewma2 && *ltyp==fix ) *ad = xte2_iglad(*l,*c,*df,*mu0,*mu1,*r,*ntyp);
   if ( *ctyp==ewma2 && *ltyp>fix )  *ad =  xte2_arlm(*l,*c,0.,*df,200,*mu0,*mu1,*ltyp,*r,nmax,*ntyp);
 } else {
   if ( *ctyp==ewma2 && *ltyp==fix ) *ad = xte2_igladc(*l,*c,*df,*mu0,*mu1,*z0,*r,*ntyp);
 }
}
