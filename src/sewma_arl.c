#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewmaU 0
#define ewmaUR 1
#define ewma2 2
#define ewmaLR 3

double seU_iglarl(double l, double cu, double hs, double sigma, int df,
                  int N, int qm);

void sewma_arl
( int *ctyp, double *l, double *cl, double *cu, double *hs, 
      double *sigma, int *df, int *r, int *qm, double *arl)
{        
 if (*ctyp==ewmaU) *arl = seU_iglarl(*l,*cu,*hs,*sigma,*df,*r,*qm);
 else *arl = -1.;
}
