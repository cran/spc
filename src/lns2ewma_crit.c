#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewmaU 0
#define ewmaL 1
#define ewma2 2

#define fixed 0
#define unbiased 1
#define eqtails 2
#define sym 3

double lns2ewmaU_crit(double l, double L0, double cl, double hs, double sigma, int df, int N);

double lns2ewma2_crit_cufix(double l, double cu, double L0, double hs, double sigma, int df, int N);
double lns2ewma2_crit_sym(double l, double L0, double hs, double sigma, int df, int N);
int lns2ewma2_crit_unbiased(double l, double L0, double *cl, double *cu, double hs, double sigma, int df, int N);

void lns2ewma_crit
( int *ctyp, int *ltyp, double *l, double *L0, double *cl0, double *cu0, double *hs, double *sigma, int *df, int *r, double *c_values)
{ int result=0;
  double cl=0., cu=1., ddf=1., mitte=0.;
  
 ddf = (double)*df;
 mitte = -1./ddf - 1./3./ddf/ddf + 2./15./ddf/ddf/ddf/ddf; 

 if ( *ctyp==ewmaU ) cu = lns2ewmaU_crit(*l, *L0, *cl0, *hs, *sigma, *df, *r);
 /*if ( *ctyp==ewmaL ) cl = lns2ewmaL_crit(*l, *L0, *cu0, *hs, *sigma, *df, *r);*/
 
 if ( *ctyp==ewma2 ) {
   if (*ltyp==fixed) {
     cl = lns2ewma2_crit_cufix(*l, *cu0, *L0, *hs, *sigma, *df, *r);
     cu = *cu0;
   }
   if ( *ltyp==unbiased ) result = lns2ewma2_crit_unbiased(*l, *L0, &cl, &cu, *hs, *sigma, *df, *r);
   /*if ( *ltyp==eqtails )  result = lns2ewma2_crit_eqtails(*l, *L0, &cl, &cu, *hs, *sigma, *df, *r);*/
   if ( *ltyp==sym )     { cl =  lns2ewma2_crit_sym(*l, *L0, *hs, *sigma, *df, *r); cu = 2*mitte - cl; }
 }
 if ( result != 0 ) warning("trouble with lns2ewma2_crit_unbiased called from lns2ewma_crit [package spc]");
 c_values[0] = cl;
 c_values[1] = cu;
}
