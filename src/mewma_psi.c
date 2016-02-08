#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define cond 0
#define cycl 1

double *vector (long n);

double mxewma_psi (double lambda, double ce, int p, int N, double *PSI, double *w, double *z);
double mxewma_psiS(double lambda, double ce, int p, double hs, int N, double *PSI, double *w, double *z);

void mewma_psi(double *l, double *c, int *p, int *type, double *hs, int *r, double *zeug)
{ double *PSI, *w, *z, zahl=0.;
  int i;
  
 PSI = vector(*r);
 w   = vector(*r);
 z   = vector(*r);
 
 if ( *type == cond ) zahl = mxewma_psi (*l, *c, *p, *r, PSI, w, z);
 if ( *type == cycl ) zahl = mxewma_psiS(*l, *c, *p, *hs, *r, PSI, w, z);
 
 zeug[0] = zahl;
 for (i = 1; i <= *r; i++) {
   zeug[i] = PSI[i-1];
   zeug[i + *r] = w[i-1];
   zeug[i + *r + *r] = z[i-1];
 }
 
 Free(z);
 Free(w);
 Free(PSI);
}
