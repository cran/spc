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
#define sven 5
#define fink 6
#define waldmann 7
#define collocation 8

extern double rho0;

double xe1_Wq(double l, double c, double p, double zr, double hs, double mu, int N, int nmax);
double xe2_Wq(double l, double c, double p, double hs, double mu, int N, int nmax);

void xewma_q
( int *ctyp, double *l, double *c, double *p, double *zr, double *hs,
  double *mu, int *r, double *q)
{
 if (*ctyp==ewma1)
   *q = xe1_Wq(*l, *c, *p, *zr, *hs, *mu, *r, 100);
 if (*ctyp==ewma2)
   *q = xe2_Wq(*l, *c, *p, *hs, *mu, *r, 10000);
}
