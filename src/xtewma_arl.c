#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewma1 0
#define ewma2 1
#define fix 0
#define vacl 1

extern double rho0;

double *vector (long n);

double xte2_iglarl(double l, double c, double hs, int df, double mu, int N, int subst);
double xte2_arlm(double l, double c, double hs, int df, int q, double mu0, double mu1, int mode, int N, int nmax, int subst); 
double xte2_arlm_hom(double l, double c, double hs, int df, int q, double mu0, double mu1, int N, double *ced, int subst);

void xtewma_arl(int *ctyp, double *l, double *c, double *zr, double *hs, int *df, double *mu, int *ltyp, int *r, int *ntyp, int *q, double *arl)
{ int nmax=100000, i, result=0;
  double *ced, arl1=-1.;
  
 ced  = vector(*q);

 if (*ctyp==ewma2 && *ltyp==fix && *q==1) arl1 = xte2_iglarl(*l,*c,*hs,*df,*mu,*r,*ntyp);
 
 if (*ctyp==ewma2 && *ltyp==fix &&  *q>1) result = xte2_arlm_hom(*l,*c,*hs,*df,*q,0.,*mu,*r,ced,*ntyp);
 
 if (*ctyp==ewma2 && *ltyp>fix )          arl1 = xte2_arlm(*l,*c,*hs,*df,*q,0.,*mu,*ltyp,*r,nmax,*ntyp);
 
 if ( result != 0 ) warning("trouble in xtewma_arl [package spc]");
 
 if ( *ltyp==fix && *q>1 ) for (i=0; i<*q; i++) arl[i] = ced[i]; else *arl = arl1;
}
