#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double tewma_arl(double lambda, int k, int lk, int uk, double z0, double mu);
double tewma_arl_R(double lambda, int k, int lk, int uk, double gl, double gu, double z0, double mu);

void tewma_arl_wowR
(int *rando, double *lambda, int *k, int *lk, int *uk, double *gl, double *gu, double *z0, double *mu, double *arl)
{ 
 *arl = -1.;
 if ( *rando==0 ) *arl = tewma_arl(*lambda, *k, *lk, *uk, *z0, *mu); 
 if ( *rando==1 ) *arl = tewma_arl_R(*lambda, *k, *lk, *uk, *gl, *gu, *z0, *mu);
}

