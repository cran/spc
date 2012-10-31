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
#define fink 6

extern double rho0;

double xe1_iglad(double l, double c, double zr, double mu0, double mu1, int N);
double xe1_arlm(double l, double c, double zr, double hs, int q, double mu0, double mu1, int mode, int N, int nmax);
double xe2_iglad(double l, double c, double mu0, double mu1, int N);
double xe2_arlm(double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax);

void xewma_ad(int *ctyp, double *l, double *c, double *zr, double *mu0, double *mu1, int *ltyp, int *r, double *ad)
{
 if (*ctyp==ewma1 && *ltyp==fix) *ad = xe1_iglad(*l,*c,*zr,*mu0,*mu1,*r);
 if (*ctyp==ewma1 && *ltyp>fix)  *ad = xe1_arlm(*l,*c,*zr,0.,200,*mu0,*mu1,*ltyp,*r,10000);
 if (*ctyp==ewma2 && *ltyp==fix) *ad = xe2_iglad(*l,*c,*mu0,*mu1,*r);
 if (*ctyp==ewma2 && *ltyp>fix)  *ad = xe2_arlm(*l,*c,0.,200,*mu0,*mu1,*ltyp,*r,10000);
}
