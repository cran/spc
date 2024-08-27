#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewma1 0
#define ewma2 1
#define fix 0
#define vacl 1
#define fir 2
#define both 3
#define steiner 4
#define stat 5
#define fink 6
#define elimit 7
#define waldmann 8
#define collocation 9

double xe2_iglarl_f(double l, double c, double mu, int N, double *g, double *w, double *z);

double *vector (long n);

void xewma_arl_f(int *ctyp, double *l, double *c, double *zr, double *mu, int *ltyp, int *r, double *zeug)
{ double *ARL, *w, *z, zahl=0.;
  int i;

 ARL = vector(*r);
 w   = vector(*r);
 z   = vector(*r); 
 
 for (i = 0; i < *r; i++) { w[i] = -1.; z[i] = 0.; ARL[i] = 0.; } /* init */
 
 if ( *ctyp==ewma2 && *ltyp==fix ) zahl =  xe2_iglarl_f(*l, *c, *mu, *r, ARL, w, z);
 
 for (i = 0; i < *r; i++) {
   zeug[i] = ARL[i];
   zeug[i + *r] = w[i];
   zeug[i + *r + *r] = z[i];
 }
   
 free(z);
 free(w);
 free(ARL);
 
 if ( fabs(zahl) > 1e-9 ) warning("trouble in xewma_arl [package spc]"); 
}  
