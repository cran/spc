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
#define elimit 7
#define waldmann 8
#define collocation 9

extern double rho0;

double *vector (long n);
double xe1_iglarl(double l, double c, double zr, double hs, double mu, int N);
double xe1_arlm(double l, double c, double zr, double hs, int q, double mu0, double mu1,int mode, int N, int nmax);
double xe1_arlm_hom(double l, double c, double zr, double hs, int q, double mu0, double mu1, int N, double *ced);
double xlimit1_arlm(double c, double zr, int q, double mu0, double mu1, int N, int nmax);
double xe1_Warl(double l, double c, double zr, double hs, double mu, int N, int nmax);

double xe2_iglarl(double l, double c, double hs, double mu, int N);
double xe2_Warl(double l, double c, double hs, double mu, int N, int nmax);
double xe2_Carl(double l, double c, double hs, double mu, int N, int qm);
double xe2_arlm(double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax);
double xe2_arlm_hom(double l, double c, double hs, int q, double mu0, double mu1, int N, double *ced);

void xewma_arl(int *ctyp, double *l, double *c, double *zr, double *hs, double *mu, int *ltyp, int *r, int *q, double *arl)
{ int nmax=100000, i, result=0;
  double *ced, arl1=-1.;
 ced  = vector(*q);

 if (*ctyp==ewma1 && *ltyp==fix && *q==1)
                     arl1 = xe1_iglarl(*l,*c,*zr,*hs,*mu,*r);
 if (*ctyp==ewma1 && *ltyp==fix && *q>1)
                     result = xe1_arlm_hom(*l, *c, *zr, *hs, *q, 0., *mu, *r, ced);
                     /* *arl = xe1_arlm(*l,*c,*zr,*hs,*q,0.,*mu,*ltyp,*r,nmax);*/
 if (*ctyp==ewma1 && *ltyp>fix && *ltyp<elimit)
                     arl1 = xe1_arlm(*l,*c,*zr,*hs,*q,0.,*mu,*ltyp,*r,nmax);
 if (*ctyp==ewma1 && *ltyp==elimit)
                     arl1 = xlimit1_arlm(*c,*zr,*q,0.,*mu,*r,nmax);
 if (*ctyp==ewma1 && *ltyp==waldmann)
                     arl1 = xe1_Warl(*l,*c,*zr,*hs,*mu,*r,nmax);

 if (*ctyp==ewma2 && *ltyp==fix && *q==1)
                     arl1 = xe2_iglarl(*l,*c,*hs,*mu,*r);
 if (*ctyp==ewma2 && *ltyp==fix && *q>1)
                     result = xe2_arlm_hom(*l, *c, *hs, *q, 0., *mu, *r, ced);                 
                     /* arl1 = xe2_arlm(*l,*c,*hs,*q,0.,*mu,*ltyp,*r,nmax);*/
 if (*ctyp==ewma2 && *ltyp>fix && *ltyp<waldmann)
                     arl1 = xe2_arlm(*l,*c,*hs,*q,0.,*mu,*ltyp,*r,nmax);
 if (*ctyp==ewma2 && *ltyp==waldmann)
                     arl1 = xe2_Warl(*l,*c,*hs,*mu,*r,nmax);
 if (*ctyp==ewma2 && *ltyp==collocation)
                     arl1 = xe2_Carl(*l,*c,*hs,*mu,*r,50);
 
 if ( result != 0 ) warning("trouble in xewma_arl [package spc]");
 
 if ( *ltyp==fix && *q>1 ) for (i=0; i<*q; i++) arl[i] = ced[i]; else *arl = arl1;
}
