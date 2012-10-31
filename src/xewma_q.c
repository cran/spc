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

double xe1_Wq(double l, double c, double p, double zr, double hs, double mu, int N, int nmax);
double xe2_Wq(double l, double c, double p, double hs, double mu, int N, int nmax);
double xe2_Wqm(double l, double c, double p, double hs, int q, double mu0, double mu1, int mode, int N, int nmax);

void xewma_q(int *ctyp, double *l, double *c, double *p, double *zr, double *hs, double *mu, int *ltyp, int *r, int *q, double *tq)
{
 if ( *ctyp==ewma1 )  *tq = xe1_Wq(*l, *c, *p, *zr, *hs, *mu, *r, 100);
 if ( *ctyp==ewma2 && *ltyp==fix && *q==1 )  *tq = xe2_Wq(*l, *c, *p, *hs, *mu, *r, 10000);
 if ( *ctyp==ewma2 && *ltyp==fix && *q>1 )   *tq = xe2_Wqm(*l, *c, *p, *hs, *q, 0., *mu, *ltyp, *r, 10000);
 if ( *ctyp==ewma2 && *ltyp>fix )            *tq = xe2_Wqm(*l, *c, *p, *hs, *q, 0., *mu, *ltyp, *r, 10000);
}
