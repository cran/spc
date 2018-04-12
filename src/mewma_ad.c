#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define GL 0
#define CO 1
#define RA 2
#define CC 3
#define MC 4
#define SR 5
#define CO2 6
#define GL2 7
#define GL3 8
#define GL4 9
#define GL5 10
#define CO3 11
#define CO4 12

#define nGL1 13
#define nGL2 14
#define nGL3 15
#define nGL4 16
#define nGL5 17

double mxewma_ad (double lambda, double ce, int p, double delta, int N, int qm2, int psi_type, double hs, int qtype, int qm0, int qm1);
double mxewma_ad_new(double lambda, double ce, int p, double delta, int N, int psi_type, double hs, int qtype);
double mxewma_ad_e(double lambda, double ce, int p, double delta, int psi_type, int N);

void mewma_ad(double *l, double *c, int *p, double *delta, int *r, int *qm2, int *ptype, double *hs, int *qtype, int *qm0, int *qm1, double *ad)
{ 
  if ( *qtype == MC ) *ad = mxewma_ad_e(*l, *c, *p, *delta, *ptype, *r);
  else {  
    if ( *qtype < nGL1 ) *ad = mxewma_ad(*l, *c, *p, *delta, *r, *qm2, *ptype, *hs, *qtype, *qm0, *qm1);
    else *ad = mxewma_ad_new(*l, *c, *p, *delta, *r, *ptype, *hs, *qtype);
  }
}
