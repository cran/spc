#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double xseU_arl_RES
  (double lx, double ls, double cx, double cs, double hsx, double hss,
   double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm, double alpha);

void xsewma_res_arl
( double *alpha, int *n, int *ctyp, 
  double *lx, double *cx, double *hsx, int *Nx,
  double *ls, double *csu, double *hss, int *Ns,
  double *mu, double *sigma, 
  int *qm, double *arl) 
{
 *arl = -1.;        
 *arl = xseU_arl_RES(*lx,*ls,*cx,*csu,*hsx,*hss,*mu,*sigma,*n,*Nx,*Ns,10000,*qm,*alpha);
}
