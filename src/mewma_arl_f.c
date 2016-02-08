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

double mxewma_arl_f_0a(double lambda, double ce, int p, int N, double *ARL, double *w, double *z);
double mxewma_arl_f_0a2(double lambda, double ce, int p, int N, double *ARL, double *w, double *z);
double mxewma_arl_f_0b(double lambda, double ce, int p, int N, int qm, double *ARL);
double mxewma_arl_f_0c(double lambda, double ce, int p, int N, double *ARL, double *w, double *z);
double mxewma_arl_f_0d(double lambda, double ce, int p, int N, double *ARL, double *w, double *z);
double mxewma_arl_f_0e(double lambda, double ce, int p, int N, double *ARL, double *z);
double mxewma_arl_f_0f(double lambda, double ce, int p, int N, double *ARL, double *w, double *z);

double mxewma_arl_f_1a (double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1); /* GL class */
double mxewma_arl_f_1a2(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1); /* GL mod */
double mxewma_arl_f_1a3(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1); /* GL again mod sin, default for 2 and 4 */
double mxewma_arl_f_1a4(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1); /* GL again mod tan */
double mxewma_arl_f_1a5(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1); /* GL again mod sinh, default for all other p */

double mxewma_arl_f_1b (double lambda, double ce, int p, double delta, int N, int qm0, int qm1, double *g); /* collocation with two halfs in the same step + sin() */
double mxewma_arl_f_1b3(double lambda, double ce, int p, double delta, int N, int qm0, int qm1, double *g); /* collocation with two halfs in the same step */
double mxewma_arl_f_1b2(double lambda, double ce, int p, double delta, int N, int qm0, int qm1, double *g); /* collocation with shrinked supports of the outer integral */
double mxewma_arl_f_1b4(double lambda, double ce, int p, double delta, int N, int qm0, int qm1, double *g); /* collocation with two halfs in the same step + sinh() instead of sin() */

double mxewma_arl_f_1c (double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1); /* GL + Radau (Rigdon) */
double mxewma_arl_f_1d (double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1); /* Clenshaw-Curtis */
double mxewma_arl_f_1e (double lambda, double ce, int p, double delta, int N, double *g, int *dQ); /* Markov Chain (Runger/Prabhu) */
double mxewma_arl_f_1f (double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z); /* Simpson rule */

double *vector (long n);

void mewma_arl_f(double *l, double *c, int *p, double *delta, int *r, int *qtype, int *qm0, int *qm1, double *zeug)
{ double *ARL, *w, *z, *w1, *z1, zahl=0.;
  int i, j, r2, dQ;
  
 if ( fabs(*delta)<1e-10 ) {
   ARL = vector(*r);
   w   = vector(*r);
   z   = vector(*r); 
 
   for (i = 0; i < *r; i++) { w[i] = -1.; z[i] = 0.; } /* init */
   
   if ( *qtype == GL )  zahl =  mxewma_arl_f_0a (*l, *c, *p, *r, ARL, w, z);
   if ( *qtype == GL2 ) zahl =  mxewma_arl_f_0a2(*l, *c, *p, *r, ARL, w, z);
   if ( *qtype == CO )  zahl =  mxewma_arl_f_0b (*l, *c, *p, *r, *qm0, ARL);   
   if ( *qtype == RA )  zahl =  mxewma_arl_f_0c (*l, *c, *p, *r, ARL, w, z);
   if ( *qtype == CC )  zahl =  mxewma_arl_f_0d (*l, *c, *p, *r, ARL, w, z);
   if ( *qtype == MC )  zahl =  mxewma_arl_f_0e (*l, *c, *p, *r, ARL, z);
   if ( *qtype == SR )  zahl =  mxewma_arl_f_0f (*l, *c, *p, *r, ARL, w, z);
   
   for (i = 0; i < *r; i++) {
     zeug[i] = ARL[i];
     zeug[i + *r] = w[i];
     zeug[i + *r + *r] = z[i];
   }
   
   Free(z);
   Free(w);
   Free(ARL);
 }
 else {
   r2 = (*r) * (*r);
   ARL = vector(r2);
   w   = vector(*r);   
   z   = vector(*r);
   w1  = vector(*r);
   z1  = vector(*r);
   
   if ( *qtype == GL )  zahl = mxewma_arl_f_1a (*l, *c, *p, *delta, *r, ARL, w, z, w1, z1);
   if ( *qtype == GL2 ) zahl = mxewma_arl_f_1a2(*l, *c, *p, *delta, *r, ARL, w, z, w1, z1);
   if ( *qtype == GL3 ) zahl = mxewma_arl_f_1a3(*l, *c, *p, *delta, *r, ARL, w, z, w1, z1);   
   if ( *qtype == GL4 ) zahl = mxewma_arl_f_1a4(*l, *c, *p, *delta, *r, ARL, w, z, w1, z1);
   if ( *qtype == GL5 ) zahl = mxewma_arl_f_1a5(*l, *c, *p, *delta, *r, ARL, w, z, w1, z1);

   if ( *qtype == CO )  zahl = mxewma_arl_f_1b (*l, *c, *p, *delta, *r, *qm0, *qm1, ARL);
   if ( *qtype == CO2 ) zahl = mxewma_arl_f_1b2(*l, *c, *p, *delta, *r, *qm0, *qm1, ARL);
   if ( *qtype == CO3 ) zahl = mxewma_arl_f_1b3(*l, *c, *p, *delta, *r, *qm0, *qm1, ARL);
   if ( *qtype == CO4 ) zahl = mxewma_arl_f_1b4(*l, *c, *p, *delta, *r, *qm0, *qm1, ARL);
   
   if ( *qtype == RA )  zahl = mxewma_arl_f_1c(*l, *c, *p, *delta, *r, ARL, w, z, w1, z1);
   if ( *qtype == CC )  zahl = mxewma_arl_f_1d(*l, *c, *p, *delta, *r, ARL, w, z, w1, z1);
   if ( *qtype == MC )  zahl = mxewma_arl_f_1e(*l, *c, *p, *delta, *r, zeug, &dQ);
   if ( *qtype == SR )  zahl = mxewma_arl_f_1f(*l, *c, *p, *delta, *r, ARL, w, z, w1, z1);
   
   if ( *qtype != MC ) {
     for (i = 0; i < *r; i++) {
       for (j = 0; j < *r; j++) zeug[ i*(*r) + j ] = ARL[ i*(*r) + j ];
       zeug[i + r2] = w[i];
       zeug[i + r2 + *r] = z[i];
       zeug[i + r2 + 2*(*r)] = w1[i];
       zeug[i + r2 + 3*(*r)] = z1[i];
     }
   } /*else {
     printf("\n\ndQ = %d\n\n", dQ);
     for (i=0; i < dQ; i++) zeug[i] = ARL[i]; 
   } */ 
   
   Free(z1);
   Free(w1);
   Free(z);
   Free(w);
   Free(ARL);
 }
 
 if ( fabs(zahl) > 1e-9 ) warning("trouble in mewma_arl_f [package spc]");
}
