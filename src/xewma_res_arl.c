#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double xe2_iglarl_RES(double l, double c, double hs, double mu, int N, double alpha, int df);

void x_res_ewma_arl(double *alpha, int *n, int *ctyp, double *l, double *c, double *hs, double *mu, int *r, double *arl)
{ 
 *arl = -1.;
 *arl = xe2_iglarl_RES(*l,*c,*hs,*mu,*r,*alpha,*n);
}
