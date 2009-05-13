#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewmaU 0
#define ewmaUR 1
#define ewma2 2
#define ewmaLR 3

double seU_iglarl(double l, double cu, double hs, double sigma, 
                  int df, int N, int qm, int s_squared);
double se2_iglarl(double l, double cl, double cu, double hs, double sigma, 
                  int df, int N, int qm);
double seUR_iglarl(double l, double cl, double cu, double hs, double sigma,
                   int df, int N, int qm);

void sewma_arl
( int *ctyp, double *l, double *cl, double *cu, double *hs, 
      double *sigma, int *df, int *r, int *qm, int *s_squared, double *arl)
{
 *arl = -1.;        
 if (*ctyp==ewmaU) *arl = seU_iglarl(*l,*cu,*hs,*sigma,*df,*r,*qm,*s_squared);
 if (*ctyp==ewma2) *arl = se2_iglarl(*l,*cl,*cu,*hs,*sigma,*df,*r,*qm);
 if (*ctyp==ewmaUR) *arl = seUR_iglarl(*l,*cl,*cu,*hs,*sigma,*df,*r,*qm);
}
