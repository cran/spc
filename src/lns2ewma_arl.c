#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewmaU 0
#define ewmaL 1
#define ewma2 2

double lns2ewmaU_arl_igl(double l, double cl, double cu, double hs, double sigma, int df, int N);
double lns2ewma2_arl_igl(double l, double cl, double cu, double hs, double sigma, int df, int N);

void lns2ewma_arl
( int *ctyp, double *l, double *cl, double *cu, double *hs, double *sigma, int *df, int *r, double *arl)
{
 *arl = -1.;
 if ( *ctyp==ewmaU ) *arl = lns2ewmaU_arl_igl(*l, *cl, *cu, *hs, *sigma, *df, *r);
 /*if ( *ctyp==ewmaL ) *arl = lns2ewmaL_arl_igl(*l, *cl, *cu, *hs, *sigma, *df, *r);*/
 if ( *ctyp==ewma2 ) *arl = lns2ewma2_arl_igl(*l, *cl, *cu, *hs, *sigma, *df, *r);
}
