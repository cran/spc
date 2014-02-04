#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewmaU 0
#define ewma2 1

double xseU_Wq(double lx, double ls, double cx, double cs, double p, double hsx, double hss,
	       double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm);
double xse2_Wq(double lx, double ls, double cx, double csl, double csu, double p, double hsx, double hss,
	       double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm);

void xsewma_q
( int *ctyp,
  double *alpha,
  double *lx, double *cx, double *hsx, int *Nx,
  double *ls, double *csl, double *csu, double *hss, int *Ns,
  double *mu, double *sigma, 
  int *df, int *qm, double *tq)
{ int nmax=100000;
 *tq = -1.;
 if ( *ctyp == ewmaU )  *tq = xseU_Wq(*lx, *ls, *cx, *csu, *alpha, *hsx, *hss, *mu, *sigma, *df, *Nx, *Ns, nmax, *qm);
 if ( *ctyp == ewma2 )  *tq = xse2_Wq(*lx, *ls, *cx, *csl, *csu, *alpha, *hsx, *hss, *mu, *sigma, *df, *Nx, *Ns, nmax, *qm);
}
