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
#define eqtails 2
#define sym 3

double seU_crit(double l, double L0, double hs, double sigma, int df, int N, int qm);
double seUR_crit(double l, double L0, double cl, double hs, double sigma, int df, int N, int qm);
double seLR_crit(double l, double L0, double cu, double hs, double sigma, int df, int N, int qm);
double se2fu_crit(double l, double L0, double cu, double hs, double sigma, int df, int N, int qm);

int se2_crit_unbiased(double l, double L0, double *cl, double *cu, double hs, double sigma, int df, int N, int qm);
int se2_crit_eqtails(double l, double L0, double *cl, double *cu, double hs, double sigma, int df, double ur, int N, int qm);
double se2_crit_sym(double l, double L0, double hs, double sigma, int df, int N, int qm);

double stdeU_crit(double l, double L0, double hs, double sigma, int df, int N, int qm);
double stdeUR_crit(double l, double L0, double cl, double hs, double sigma, int df, int N, int qm);
double stdeLR_crit(double l, double L0, double cu, double hs, double sigma, int df, int N, int qm);
double stde2fu_crit(double l, double L0, double cu, double hs, double sigma, int df, int N, int qm);

int stde2_crit_unbiased(double l, double L0, double *cl, double *cu, double hs, double sigma, int df, int N, int qm);
int stde2_crit_eqtails(double l, double L0, double *cl, double *cu, double hs, double sigma, int df, double ur, int N, int qm);
double stde2_crit_sym(double l, double L0, double hs, double sigma, int df, int N, int qm);

double c_four(double ddf);

void sewma_crit
( int *ctyp, int *ltyp, double *l, double *L0, double *cl0, double *cu0, double *hs, double *sigma, int *df, int *r, int *qm, double *ur, int *s_squared, double *c_values)
{ int result=0;
  double cl=0., cu=1., mitte=1.;

 if ( *s_squared==1 ) {
   if (*ctyp==ewmaU)  cu = seU_crit(*l,*L0,*hs,*sigma,*df,*r,*qm);
   if (*ctyp==ewmaUR) cu = seUR_crit(*l,*L0,*cl0,*hs,*sigma,*df,*r,*qm);
   if (*ctyp==ewmaLR) cl = seLR_crit(*l,*L0,*cu0,*hs,*sigma,*df,*r,*qm);
   if (*ctyp==ewma2) {
     if (*ltyp==fixed) {
       cl = se2fu_crit(*l,*L0,*cu0,*hs,*sigma,*df,*r,*qm);
       cu = *cu0;
     }
     if (*ltyp==unbiased) result = se2_crit_unbiased(*l, *L0, &cl, &cu, *hs, *sigma, *df, *r, *qm);
     if (*ltyp==eqtails)  result =  se2_crit_eqtails(*l, *L0, &cl, &cu, *hs, *sigma, *df, *ur, *r, *qm);
     if (*ltyp==sym)     {    cu =  se2_crit_sym(*l, *L0, *hs, *sigma, *df, *r, *qm); cl = 2. - cu; }
   }
 } else {
   mitte = c_four((double)*df);
   if ( *ctyp==ewmaU )  cu = stdeU_crit(*l,*L0,*hs,*sigma,*df,*r,*qm);
   if ( *ctyp==ewmaUR ) cu = stdeUR_crit(*l,*L0,*cl0,*hs,*sigma,*df,*r,*qm);
   if ( *ctyp==ewmaLR ) cl = stdeLR_crit(*l,*L0,*cu0,*hs,*sigma,*df,*r,*qm);
   if ( *ctyp==ewma2 ) {
     if ( *ltyp==fixed ) {
       cl = stde2fu_crit(*l,*L0,*cu0,*hs,*sigma,*df,*r,*qm);
       cu = *cu0;
     }
     if ( *ltyp==unbiased ) result = stde2_crit_unbiased(*l, *L0, &cl, &cu, *hs, *sigma, *df, *r, *qm);
     if ( *ltyp==eqtails )  result =  stde2_crit_eqtails(*l, *L0, &cl, &cu, *hs, *sigma, *df, *ur, *r, *qm);
     if ( *ltyp==sym )     {    cu =  stde2_crit_sym(*l, *L0, *hs, *sigma, *df, *r, *qm); cl = 2.*mitte - cu; }
   }
 }
 
 if ( result != 0 ) warning("trouble with se2_crit called from sewma_crit [package spc]");
 
 c_values[0] = cl;
 c_values[1] = cu;
}
