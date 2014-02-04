#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define cusumU 0
#define cusumL 1
#define cusum2 2

double scU_iglarl_v1(double refk, double h, double hs, double sigma, int df, int N, int qm);
double scU_iglarl_v2(double refk, double h, double hs, double sigma, int df, int N, int qm);
double scL_iglarl_v2(double refk, double h, double hs, double sigma, int df, int N, int qm);
double sc2_iglarl_v2(double refkl, double refku, double hl, double hu, double hsl, double hsu, double sigma, int df, int N, int qm);

void scusum_arl
( int *ctyp, double *k, double *h, double *hs, double *sigma, int *df, double *k2, double *h2, double *hs2, int *r, int *qm, int *version, double *arl)
{
 *arl = -1.; 
 if ( *ctyp==cusumU )  {
   if ( *version==1 ) *arl = scU_iglarl_v1(*k, *h, *hs, *sigma, *df, *r, *qm);
   if ( *version==2 ) *arl = scU_iglarl_v2(*k, *h, *hs, *sigma, *df, *r, *qm);
 }
 if ( *ctyp==cusumL )  {
   /*if ( *version==1 ) *arl = scL_iglarl_v1(*k, *h, *hs, *sigma, *df, *r, *qm);*/
   if ( *version==2 ) *arl = scL_iglarl_v2(*k, *h, *hs, *sigma, *df, *r, *qm);
 }
 if ( *ctyp==cusum2 )  {
   /*if ( *version==1 ) *arl = sc2_iglarl_v1(*k2, *k, *h2, *h, *hs2, *hs, *sigma, *df, *r, *qm);*/
   if ( *version==2 ) *arl = sc2_iglarl_v2(*k2, *k, *h2, *h, *hs2, *hs, *sigma, *df, *r, *qm);
 }
}
