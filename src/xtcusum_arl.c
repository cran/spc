#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define cusum1 0
#define cusum2 1

double xtc1_iglarl(double k, double h, double hs, int df, double mu, int N, int subst);
double xtc2_iglarl(double k, double h, double hs, int df, double mu, int N, int subst);

void xtcusum_arl
( int *ctyp, double *k, double *h, double *hs, int *df, double *mu, int *r, int *ntyp, double *arl)
{ 
 if ( *ctyp == cusum1 ) *arl = xtc1_iglarl(*k, *h, *hs, *df, *mu, *r, *ntyp);
 if ( *ctyp == cusum2 ) *arl = xtc2_iglarl(*k, *h, *hs, *df, *mu, *r, *ntyp);
}
