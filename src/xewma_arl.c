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

double xe1_iglarl(double l, double c, double zr, double hs, double mu, int N);
double xe1_arlm(double l, double c, double zr, double hs, int q, double mu0, double mu1,int mode, int N, int nmax);
double xlimit1_arlm(double c, double zr, int q, double mu0, double mu1, int N, int nmax);
double xe1_Warl(double l, double c, double zr, double hs, double mu, int N, int nmax);

double xe2_iglarl(double l, double c, double hs, double mu, int N);
double xe2_Warl(double l, double c, double hs, double mu, int N, int nmax);
double xe2_Carl(double l, double c, double hs, double mu, int N, int qm);
double xe2_arlm(double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax);

void xewma_arl(int *ctyp, double *l, double *c, double *zr, double *hs, double *mu, int *ltyp, int *r, int *q, double *arl)
{ int nmax=100000;

 if (*ctyp==ewma1 && *ltyp==fix && *q==1)
                     *arl = xe1_iglarl(*l,*c,*zr,*hs,*mu,*r);
 if (*ctyp==ewma1 && *ltyp==fix && *q>1)
                     *arl = xe1_arlm(*l,*c,*zr,*hs,*q,0.,*mu,*ltyp,*r,nmax);
 if (*ctyp==ewma1 && *ltyp>fix && *ltyp<elimit)
                     *arl = xe1_arlm(*l,*c,*zr,*hs,*q,0.,*mu,*ltyp,*r,nmax);
 if (*ctyp==ewma1 && *ltyp==elimit)
                   *arl = xlimit1_arlm(*c,*zr,*q,0.,*mu,*r,nmax);
 if (*ctyp==ewma1 && *ltyp==waldmann)
                   *arl = xe1_Warl(*l,*c,*zr,*hs,*mu,*r,nmax);

 if (*ctyp==ewma2 && *ltyp==fix && *q==1)
                     *arl = xe2_iglarl(*l,*c,*hs,*mu,*r);
 if (*ctyp==ewma2 && *ltyp==fix && *q>1)
                     *arl = xe2_arlm(*l,*c,*hs,*q,0.,*mu,*ltyp,*r,nmax);
 if (*ctyp==ewma2 && *ltyp>fix && *ltyp<waldmann)
                     *arl = xe2_arlm(*l,*c,*hs,*q,0.,*mu,*ltyp,*r,nmax);
 if (*ctyp==ewma2 && *ltyp==waldmann)
                     *arl = xe2_Warl(*l,*c,*hs,*mu,*r,nmax);
 if (*ctyp==ewma2 && *ltyp==collocation)
                     *arl = xe2_Carl(*l,*c,*hs,*mu,*r,50);
}
