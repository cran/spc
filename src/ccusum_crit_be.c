#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define cusumU 0
#define cusumL 1
#define cusum2 2

int ccusum_U_crit(double A, double mu0, int km, int m, int i0);
int ccusum_U_rando_crit(double A, double mu0, int km, int m, int i0, int *hm, double *gamma);

int ccusum_L_crit(double A, double mu0, int km, int m, int i0);
int ccusum_L_rando_crit(double A, double mu0, int km, int m, int i0, int *hm, double *gamma);

void ccusum_crit_be
(int *ctyp, int *rando, double *mu0, int *km, double *A, int *m, int *i0, double *c_values)
{ double gamma=0.;
  int result=0, hm=0;
  
 if ( *ctyp==cusumU && *rando==0 )  hm = ccusum_U_crit(*A, *mu0, *km, *m, *i0);
 if ( *ctyp==cusumU && *rando==1 )  result = ccusum_U_rando_crit(*A, *mu0, *km, *m, *i0, &hm, &gamma);
 
 if ( *ctyp==cusumL && *rando==0 )  hm = ccusum_L_crit(*A, *mu0, *km, *m, *i0);
 if ( *ctyp==cusumL && *rando==1 )  result = ccusum_L_rando_crit(*A, *mu0, *km, *m, *i0, &hm, &gamma);
 
 if ( result != 0 ) warning("something went wrong with ccusum_*_rando_crit");
 
 c_values[0] = (double)hm;
 c_values[1] = gamma;
}
