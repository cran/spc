#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewmaU 0
#define ewmaUR 1
#define ewma2 2
#define ewmaLR 3
#define fixed 0
#define unbiased 1

int xseU_crit
  (double lx, double ls, double L0, double *cx, double *cs,
   double hsx, double hss,
   double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm);
int xse2lu_crit
  (double lx, double ls, double L0, double *cx, double csl, double *csu, 
   double hsx, double hss,
   double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm);
int xse2fu_crit
  (double lx, double ls, double L0, double *cx, double *csl, double csu, 
   double hsx, double hss,
   double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm);
int xse2_crit
  (double lx, double ls, double L0, double *cx, double *csl, double *csu, 
   double hsx, double hss,
   double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm);

void xsewma_crit
( int *ctyp, int *ltyp,
      double *lx, double *ls, double *L0, double *cu0, double *hsx, double *hss,
      double *mu, double *sigma, int *df,
      int *Nx, int *Ns, int *qm, double *c_values)
{ int result;
  double cx, cl, cu;
 cx = -1.;
 cl = 0.;
 cu = -1.;

 if (*ctyp==ewmaU) 
   result = xseU_crit(*lx,*ls,*L0,&cx,&cu,*hsx,*hss,*mu,*sigma,*df,*Nx,*Ns,10000,*qm);
 if (*ctyp==ewma2) {
   if (*ltyp==fixed) {
     result = xse2fu_crit(*lx,*ls,*L0,&cx,&cl,*cu0,*hsx,*hss,*mu,*sigma,*df,*Nx,*Ns,10000,*qm);
     cu = *cu0;
   }
   if (*ltyp==unbiased)
     result = xse2_crit(*lx,*ls,*L0,&cx,&cl,&cu,*hsx,*hss,*mu,*sigma,*df,*Nx,*Ns,10000,*qm);
 }

 c_values[0] = cx;
 c_values[1] = cl;
 c_values[2] = cu;
}
