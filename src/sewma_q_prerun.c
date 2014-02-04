#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewmaU 0
#define ewmaUR 1
#define ewma2 2
#define ewmaLR 3

double seU_Wq_prerun_SIGMA_deluxe(double l, double cu, double p, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate);
double seUR_Wq_prerun_SIGMA_deluxe(double l, double cl, double cu, double p, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate);
double seLR_Wq_prerun_SIGMA_deluxe(double l, double cl, double cu, double p, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate);
double se2_Wq_prerun_SIGMA_deluxe(double l, double cl, double cu, double p, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate);

void sewma_q_prerun
( int *ctyp, double *l, double *cl, double *cu, double *p, double *hs, double *sigma, int *df1, int *r, int *qm1, int *df2, int *qm2, double *truncate, double *tq)
{ int nmax=100000;
 
 if ( *ctyp == ewmaU )  *tq =  seU_Wq_prerun_SIGMA_deluxe(*l, *cu, *p, *hs, *sigma, *df1, *df2, nmax, *qm1, *qm2, *truncate);
 if ( *ctyp == ewma2 )  *tq =  se2_Wq_prerun_SIGMA_deluxe(*l, *cl, *cu, *p, *hs, *sigma, *df1, *df2, nmax, *qm1, *qm2, *truncate);
 if ( *ctyp == ewmaUR ) *tq = seUR_Wq_prerun_SIGMA_deluxe(*l, *cl, *cu, *p, *hs, *sigma, *df1, *df2, nmax, *qm1, *qm2, *truncate);
 if ( *ctyp == ewmaLR ) *tq = seLR_Wq_prerun_SIGMA_deluxe(*l, *cl, *cu, *p, *hs, *sigma, *df1, *df2, nmax, *qm1, *qm2, *truncate);

}
