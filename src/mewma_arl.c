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

double mxewma_arl_0a(double lambda, double ce, int p, double hs, int N);
double mxewma_arl_0a2(double lambda, double ce, int p, double hs, int N);
double mxewma_arl_0b(double lambda, double ce, int p, double hs, int N, int qm);
double mxewma_arl_0c(double lambda, double ce, int p, double hs, int N);
double mxewma_arl_0d(double lambda, double ce, int p, double hs, int N);
double mxewma_arl_0e(double lambda, double ce, int p, double hs, int N);
double mxewma_arl_0f(double lambda, double ce, int p, double hs, int N);

double mxewma_arl_1a(double lambda, double ce, int p, double delta, double hs, int N);
double mxewma_arl_1a2(double lambda, double ce, int p, double delta, double hs, int N);
double mxewma_arl_1a3(double lambda, double ce, int p, double delta, double hs, int N);
double mxewma_arl_1a4(double lambda, double ce, int p, double delta, double hs, int N);
double mxewma_arl_1a5(double lambda, double ce, int p, double delta, double hs, int N);

double mxewma_arl_1b(double lambda, double ce, int p, double delta, double hs, int N, int qm0, int qm1);
double mxewma_arl_1b2(double lambda, double ce, int p, double delta, double hs, int N, int qm0, int qm1);
double mxewma_arl_1b3(double lambda, double ce, int p, double delta, double hs, int N, int qm0, int qm1);
double mxewma_arl_1b4(double lambda, double ce, int p, double delta, double hs, int N, int qm0, int qm1);
double mxewma_arl_1c(double lambda, double ce, int p, double delta, double hs, int N);
double mxewma_arl_1d(double lambda, double ce, int p, double delta, double hs, int N);
double mxewma_arl_1e(double lambda, double ce, int p, double delta, double hs, int N);
double mxewma_arl_1f(double lambda, double ce, int p, double delta, double hs, int N);

void mewma_arl(double *l, double *c, int *p, double *delta, double *hs, int *r, int *qtype, int *qm0, int *qm1, double *arl)
{
 if ( fabs(*delta)<1e-10 ) {
   if ( *qtype == GL )  *arl =  mxewma_arl_0a(*l, *c, *p, *hs, *r);
   if ( *qtype == GL2 ) *arl =  mxewma_arl_0a2(*l, *c, *p, *hs, *r);
   if ( *qtype == CO )  *arl =  mxewma_arl_0b(*l, *c, *p, *hs, *r, *qm0);   
   if ( *qtype == RA )  *arl =  mxewma_arl_0c(*l, *c, *p, *hs, *r);
   if ( *qtype == CC )  *arl =  mxewma_arl_0d(*l, *c, *p, *hs, *r);
   if ( *qtype == MC )  *arl =  mxewma_arl_0e(*l, *c, *p, *hs, *r);
   if ( *qtype == SR )  *arl =  mxewma_arl_0f(*l, *c, *p, *hs, *r);
 }
 else {
   if ( *qtype == GL )  *arl = mxewma_arl_1a(*l, *c, *p, *delta, *hs, *r);
   if ( *qtype == GL2 ) *arl = mxewma_arl_1a2(*l, *c, *p, *delta, *hs, *r);
   if ( *qtype == GL3 ) *arl = mxewma_arl_1a3(*l, *c, *p, *delta, *hs, *r);   
   if ( *qtype == GL4 ) *arl = mxewma_arl_1a4(*l, *c, *p, *delta, *hs, *r);
   if ( *qtype == GL5 ) *arl = mxewma_arl_1a5(*l, *c, *p, *delta, *hs, *r);
   
   if ( *qtype == CO )  *arl = mxewma_arl_1b(*l, *c, *p, *delta, *hs, *r, *qm0, *qm1);
   if ( *qtype == CO2 ) *arl = mxewma_arl_1b2(*l, *c, *p, *delta, *hs, *r, *qm0, *qm1);
   if ( *qtype == CO3 ) *arl = mxewma_arl_1b3(*l, *c, *p, *delta, *hs, *r, *qm0, *qm1);
   if ( *qtype == CO4 ) *arl = mxewma_arl_1b4(*l, *c, *p, *delta, *hs, *r, *qm0, *qm1);
   if ( *qtype == RA ) *arl = mxewma_arl_1c(*l, *c, *p, *delta, *hs, *r);
   if ( *qtype == CC ) *arl = mxewma_arl_1d(*l, *c, *p, *delta, *hs, *r);
   if ( *qtype == MC ) *arl = mxewma_arl_1e(*l, *c, *p, *delta, *hs, *r);
   if ( *qtype == SR ) *arl = mxewma_arl_1f(*l, *c, *p, *delta, *hs, *r);
 }
}
