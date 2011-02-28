#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double seU_iglarl_RES
  (double l, double cu, double hs, double sigma, int df, int N, int qm, double alpha, double mu);

void s_res_ewma_arl
( double *alpha, int *n, int *ctyp, double *l, double *cu, double *hs,
  double *sigma, double *mu, int *r, int *qm, double *arl)
{ 
 *arl = -1.;
 *arl = seU_iglarl_RES(*l,*cu,*hs,*sigma,*n,*r,*qm,*alpha,*mu);
}
