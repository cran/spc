#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewmaU 0
#define ewmaUR 1
#define ewma2 2
#define ewmaLR 3

double  seU_iglarl_prerun_SIGMA(double l, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate);
double seUR_iglarl_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate);
double  se2_iglarl_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate);
double seLR_iglarl_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate);

void sewma_arl_prerun
( int *ctyp, double *l, double *cl, double *cu, double *hs, double *sigma, int *df1, int *N, int *qm1, int *df2, int *qm2, double *truncate, double *arl)
{ 
  *arl = -1.;
  if ( *ctyp==ewmaU )  *arl =  seU_iglarl_prerun_SIGMA(*l, *cu, *hs, *sigma, *df1, *df2, *N, *qm1, *qm2, *truncate);
  if ( *ctyp==ewma2 )  *arl =  se2_iglarl_prerun_SIGMA(*l, *cl, *cu, *hs, *sigma, *df1, *df2, *N, *qm1, *qm2, *truncate);
  if ( *ctyp==ewmaUR ) *arl = seUR_iglarl_prerun_SIGMA(*l, *cl, *cu, *hs, *sigma, *df1, *df2, *N, *qm1, *qm2, *truncate);
  if ( *ctyp==ewmaLR ) *arl = seLR_iglarl_prerun_SIGMA(*l, *cl, *cu, *hs, *sigma, *df1, *df2, *N, *qm1, *qm2, *truncate);
}
