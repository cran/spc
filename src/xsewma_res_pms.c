#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double xseU_mu_before_sigma_RES
  (double lx, double ls, double cx, double cs, double hsx, double hss,
   double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm, double alpha, int vice_versa);

void xsewma_res_pms
( double *alpha, int *n, int *ctyp, 
  double *lx, double *cx, double *hsx, int *Nx,
  double *ls, double *csu, double *hss, int *Ns,
  double *mu, double *sigma, 
  int *qm, int *vice_versa, double *pms) 
{
 *pms = -1.;        
 *pms = xseU_mu_before_sigma_RES(*lx,*ls,*cx,*csu,*hsx,*hss,*mu,*sigma,*n,*Nx,*Ns,10000,*qm,*alpha,*vice_versa);
}
