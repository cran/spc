#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewmaU 0
#define ewmaUR 1
#define ewma2 2
#define ewmaLR 3

double seU_Wq(double l, double cu, double p, double hs, double sigma, int df, int N, int nmax, int qm);
double se2_Wq(double l, double cl, double cu, double p, double hs, double sigma, int df, int N, int nmax, int qm);
double seUR_Wq(double l, double cl, double cu, double p, double hs, double sigma, int df, int N, int nmax, int qm);
double seLR_Wq(double l, double cl, double cu, double p, double hs, double sigma, int df, int N, int nmax, int qm);

void sewma_q(int *ctyp, double *l, double *cl, double *cu, double *p, double *hs, int *N, double *sigma, int *df, int *qm, double *tq)
{ int nmax=100000;
 
 if ( *ctyp == ewmaU )  *tq = seU_Wq(*l, *cu, *p, *hs, *sigma, *df, *N, nmax, *qm);
 if ( *ctyp == ewma2 )  *tq = se2_Wq(*l, *cl, *cu, *p, *hs, *sigma, *df, *N, nmax, *qm);
 if ( *ctyp == ewmaUR ) *tq = seUR_Wq(*l, *cl, *cu, *p, *hs, *sigma, *df, *N, nmax, *qm);
 if ( *ctyp == ewmaLR ) *tq = seLR_Wq(*l, *cl, *cu, *p, *hs, *sigma, *df, *N, nmax, *qm);

}
