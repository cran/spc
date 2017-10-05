#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define cusumU 0
#define cusumL 1
#define cusum2 2

double scs_U_iglarl_v1(double refk, double h, double hs, double cS, double sigma, int df, int N, int qm);

void scusum_s_arl
( int *ctyp, double *k, double *h, double *hs, double *cS, double *sigma, int *df, double *k2, double *h2, double *hs2, int *r, int *qm, int *version, double *arl)
{
 *arl = -1.; 
 if ( *ctyp==cusumU )  {
   if ( *version==1 ) *arl = scs_U_iglarl_v1(*k, *h, *hs, *cS, *sigma, *df, *r, *qm);
   if ( *version==2 ) *arl = scs_U_iglarl_v1(*k, *h, *hs, *cS, *sigma, *df, *r, *qm);
 }
}
