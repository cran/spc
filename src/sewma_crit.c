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

double seU_crit(double l, double L0, double hs, double sigma, int df, int N, int qm, int s_squared);
double seUR_crit(double l, double L0, double cl, double hs, double sigma, int df, int N, int qm, int s_squared);
double seLR_crit(double l, double L0, double cu, double hs, double sigma, int df, int N, int qm, int s_squared);
double se2fu_crit(double l, double L0, double cu, double hs, double sigma, int df, int N, int qm);
int se2_crit(double l, double L0, double *cl, double *cu, double hs, double sigma, int df, int N, int qm);

void sewma_crit
( int *ctyp, int *ltyp, double *l, double *L0, double *cl0, double *cu0, double *hs, double *sigma, int *df, int *r, int *qm, int *s_squared, double *c_values)
{ int result;
  double cl=0., cu=1.;

 if (*ctyp==ewmaU)  cu = seU_crit(*l,*L0,*hs,*sigma,*df,*r,*qm,*s_squared);
 if (*ctyp==ewmaUR) cu = seUR_crit(*l,*L0,*cl0,*hs,*sigma,*df,*r,*qm,*s_squared);
 if (*ctyp==ewmaLR) cl = seLR_crit(*l,*L0,*cu0,*hs,*sigma,*df,*r,*qm,*s_squared);
 if (*ctyp==ewma2) {
   if (*ltyp==fixed) {
     cl = se2fu_crit(*l,*L0,*cu0,*hs,*sigma,*df,*r,*qm);
     cu = *cu0;
   }
   if (*ltyp==unbiased) result = se2_crit(*l,*L0,&cl,&cu,*hs,*sigma,*df,*r,*qm);
 }

 c_values[0] = cl;
 c_values[1] = cu;
}
