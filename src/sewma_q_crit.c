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
#define classic 2

double seU_q_crit(double l, int L0, double alpha, double hs, double sigma, int df, int N, int qm, double c_error, double a_error);
/*double se2lu_q_crit(double l, int L0, double alpha, double cl, double hs, double sigma, int df, int N, int qm, double c_error, double a_error);*/
double se2fu_q_crit(double l, int L0, double alpha, double cu, double hs, double sigma, int df, int N, int qm, double c_error, double a_error);
int se2_q_crit(double l, int L0, double alpha, double *cl, double *cu, double hs, double sigma, int df, int N, int qm, double c_error, double a_error);
int se2_q_crit_class(double l, int L0, double alpha, double *cl, double *cu, double hs, double sigma, int df, double ur, int N, int qm, double c_error, double a_error);
double seUR_q_crit(double l, int L0, double alpha, double cl, double hs, double sigma, int df, int N, int qm, double c_error, double a_error);
double seLR_q_crit(double l, int L0, double alpha, double cu, double hs, double sigma, int df, int N, int qm, double c_error, double a_error);

void sewma_q_crit
( int *ctyp, int *ltyp, double *l, int *L0, double *alpha, double *cl0, double *cu0, double *hs, double *sigma, int *df, int *r, int *qm, double *ur,
  double *c_error, double *a_error, double *c_values)
{ int result=0;
  double cl=0., cu=1.;

 if ( *ctyp==ewmaU )  { cu = seU_q_crit(*l, *L0, *alpha, *hs, *sigma, *df, *r, *qm, *c_error, *a_error); cl = 0.; }
 if ( *ctyp==ewmaUR ) { cu = seUR_q_crit(*l, *L0, *alpha, *cl0, *hs, *sigma, *df, *r, *qm, *c_error, *a_error); cl = *cl0; }
 if ( *ctyp==ewmaLR ) { cl = seLR_q_crit(*l, *L0, *alpha, *cu0, *hs, *sigma, *df, *r, *qm, *c_error, *a_error); cu = *cu0; }
 if ( *ctyp==ewma2 ) {
   if ( *ltyp==fixed ) {
     cl = se2fu_q_crit(*l, *L0, *alpha, *cu0, *hs, *sigma, *df, *r, *qm, *c_error, *a_error);
     cu = *cu0;
   }
   if ( *ltyp==unbiased ) result = se2_q_crit(*l, *L0, *alpha, &cl, &cu, *hs, *sigma, *df, *r, *qm, *c_error, *a_error);
   if ( *ltyp==classic )  result = se2_q_crit_class(*l, *L0, *alpha, &cl, &cu, *hs, *sigma, *df, *ur, *r, *qm, *c_error, *a_error);
 }
 
 if ( result != 0 ) warning("trouble with se2_crit called from sewma_q_crit [package spc]");
 
 c_values[0] = cl;
 c_values[1] = cu;
}
