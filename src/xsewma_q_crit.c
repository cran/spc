#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewmaU 0
#define ewma2 1

#define fixed 0
#define unbiased 1

int xseU_q_crit(double lx, double ls, int L0, double alpha, double *cx, double *cs, double hsx, double hss, double mu, double sigma, int df,
                int Nx, int Ns, int qm, double c_error, double a_error);
int xse2fu_q_crit(double lx, double ls, int L0, double alpha, double *cx, double *csl, double csu, double hsx, double hss, double mu, double sigma, int df,
                  int Nx, int Ns, int qm, double c_error, double a_error);
int xse2_q_crit(double lx, double ls, int L0, double alpha, double *cx, double *csl, double *csu, double hsx, double hss, double mu, double sigma, int df,
                int Nx, int Ns, int qm, double c_error, double a_error);
                
void xsewma_q_crit
( int *ctyp, int *ltyp,
  double *lx, double *ls, double *L0, double *alpha, double *cu0, double *hsx, double *hss,
  double *mu, double *sigma, int *df,
  int *Nx, int *Ns, int *qm,
  double *c_error, double *a_error, double *c_values)
  
{ int result=0;
  double cx=-1., cl=0., cu=-1.;

 if ( *ctyp==ewmaU ) result = xseU_q_crit(*lx, *ls, *L0, *alpha, &cx, &cu, *hsx, *hss, *mu, *sigma, *df, *Nx, *Ns, *qm, *c_error, *a_error);
 
 if ( *ctyp==ewma2 ) {
   if ( *ltyp==fixed ) {
     result = xse2fu_q_crit(*lx, *ls, *L0, *alpha, &cx, &cl, *cu0, *hsx, *hss, *mu, *sigma, *df, *Nx, *Ns, *qm, *c_error, *a_error);
     cu = *cu0;
   }
   if ( *ltyp==unbiased )
     result = xse2_q_crit(*lx, *ls, *L0, *alpha, &cx, &cl, &cu, *hsx, *hss, *mu, *sigma, *df, *Nx, *Ns, *qm, *c_error, *a_error);
 }
 
 if ( result != 0 ) warning("trouble with xsewma_q_crit [package spc]");
 
 c_values[0] = cx;
 c_values[1] = cl;
 c_values[2] = cu;
 
}
