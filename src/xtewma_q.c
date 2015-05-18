#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewma1 0
#define ewma2 1
#define fix 0
#define vacl 1

/*double xe1_Wq(double l, double c, double p, double zr, double hs, double mu, int N, int nmax);
double xe1_Wqm(double l, double c, double p, double zr, double hs, int q, double mu0, double mu1, int mode, int N, int nmax);*/

double xte2_Wq(double l, double c, double p, double hs, int df, double mu, int N, int nmax, int subst);

double xte2_Wqm(double l, double c, double p, double hs, int df, int q, double mu0, double mu1, int mode, int N, int nmax, int subst);

void xtewma_q(int *ctyp, double *l, double *c, double *p, double *zr, double *hs, int *df, double *mu, int *ltyp, int *r, int *ntyp, int *q, double *tq)
{ int nmax=1000000;
  
 if ( *ctyp==ewma2 && *ltyp==fix && *q==1 )  *tq = xte2_Wq(*l, *c, *p, *hs, *df, *mu, *r, nmax, *ntyp);
 if ( *ctyp==ewma2 && *ltyp==fix && *q>1 )   *tq = xte2_Wqm(*l, *c, *p, *hs, *df, *q, 0., *mu, *ltyp, *r, nmax, *ntyp);
 if ( *ctyp==ewma2 && *ltyp>fix )            *tq = xte2_Wqm(*l, *c, *p, *hs, *df, *q, 0., *mu, *ltyp, *r, nmax, *ntyp);
}
