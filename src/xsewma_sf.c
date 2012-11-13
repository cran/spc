#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewmaU 0
#define ewmaUR 1
#define ewma2 2

double *vector (long n);
double xseU_sf(double lx, double ls, double cx, double cs, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm, double *p0);
double xse2_sf(double lx, double ls, double cx, double csl, double csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm, double *p0);

void xsewma_sf
( int *ctyp, 
      double *lx, double *cx, double *hsx, int *Nx,
      double *ls, double *csl, double *csu, double *hss, int *Ns,
      double *mu, double *sigma, 
      int *df, int *qm, int *n, double *sf)
{ int result=0, i;
  double *p0;
 p0  = vector(*n);

 if ( *ctyp == ewmaU )
   result = xseU_sf(*lx, *ls, *cx, *csu, *hsx, *hss, *mu, *sigma, *df, *Nx, *Ns, *n, *qm, p0);
 if ( *ctyp == ewma2 )
   result = xse2_sf(*lx, *ls, *cx, *csl, *csu, *hsx, *hss, *mu, *sigma, *df, *Nx, *Ns, *n, *qm, p0);

 if ( result != 0 ) warning("trouble in xsewma_sf [package spc]");

 for (i=0; i<*n; i++) sf[i] = p0[i];
}
