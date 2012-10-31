#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewmaU 0
#define ewmaUR 1
#define ewma2 2

double xseU_arl
  (double lx, double ls, double cx, double cs, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm);
double xse2_arl
  (double lx, double ls, double cx, double csl, double csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm);

void xsewma_arl
( int *ctyp, 
      double *lx, double *cx, double *hsx, int *Nx,
      double *ls, double *csl, double *csu, double *hss, int *Ns,
      double *mu, double *sigma, 
      int *df, int *qm, int *s_squared, double *arl) 
{
 *arl = -1.;        
 if (*ctyp==ewmaU) 
   *arl = xseU_arl(*lx,*ls,*cx,*csu,*hsx,*hss,*mu,*sigma,*df,*Nx,*Ns,10000,*qm);
 if (*ctyp==ewma2) 
   *arl = xse2_arl(*lx,*ls,*cx,*csl,*csu,*hsx,*hss,*mu,*sigma,*df,*Nx,*Ns,10000,*qm);
}
