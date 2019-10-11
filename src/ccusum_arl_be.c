#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define cusumU 0
#define cusumL 1
#define cusum2 2

double ccusum_U_arl(double mu, int km, int hm, int m, int i0);
double ccusum_U_arl_rando(double mu, int km, int hm, int m, double gamma, int i0);

double ccusum_L_arl(double mu, int km, int hm, int m, int i0);
double ccusum_L_arl_rando(double mu, int km, int hm, int m, double gamma, int i0);

double ccusum_2_arl(double mu, int km1, int hm1, int m1, int i01, int km2, int hm2, int m2, int i02);
double ccusum_2_arl_rando(double mu, int km1, int hm1, int m1, double gamma1, int i01, int km2, int hm2, int m2, double gamma2, int i02);

void ccusum_arl_be
(int *ctyp, int *rando, double *mu,
 int *km, int *hm, int *m, int *i0, double *gamma,
 int *km2, int *hm2, int *m2, int *i02, double *gamma2,
 double *arl)
{ 
 *arl = -1.;
 if ( *ctyp==cusumU && *rando==0 )  *arl = ccusum_U_arl(*mu, *km, *hm, *m, *i0);
 if ( *ctyp==cusumU && *rando==1 )  *arl = ccusum_U_arl_rando(*mu, *km, *hm, *m, *gamma, *i0);
 if ( *ctyp==cusumL && *rando==0 )  *arl = ccusum_L_arl(*mu, *km, *hm, *m, *i0);
 if ( *ctyp==cusumL && *rando==1 )  *arl = ccusum_L_arl_rando(*mu, *km, *hm, *m, *gamma, *i0);
 if ( *ctyp==cusum2 && *rando==0 )  *arl = ccusum_2_arl(*mu, *km, *hm, *m, *i0, *km2, *hm2, *m2, *i02);
 if ( *ctyp==cusum2 && *rando==1 )  *arl = ccusum_2_arl_rando(*mu, *km, *hm, *m, *gamma, *i0, *km2, *hm2, *m2, *gamma2, *i02);
}
