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

double seU_q_crit_prerun_SIGMA(double l, int L0, double alpha, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2,
			       double truncate, int tail_approx, double c_error, double a_error);
double se2fu_q_crit_prerun_SIGMA(double l, int L0, double alpha, double cu, double s0, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2,
				 double truncate, int tail_approx, double c_error, double a_error);
int se2_q_crit_prerun_SIGMA(double l, int L0, double alpha, double *cl, double *cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2,
			   double truncate, int tail_approx, double c_error, double a_error);
double seUR_q_crit_prerun_SIGMA(double l, int L0, double alpha, double cl, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2,
				double truncate, int tail_approx, double c_error, double a_error);
double seLR_q_crit_prerun_SIGMA(double l, int L0, double alpha, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2,
				double truncate, int tail_approx, double c_error, double a_error);

void sewma_q_crit_prerun
( int *ctyp, int *ltyp, double *l, int *L0, double *alpha, double *cl0, double *cu0, double *hs, double *sigma, int *df1, int *r, int *qm1,
  int *df2, int *qm2, double *truncate, int *tail_approx, double *c_error, double *a_error, double *c_values)
{ int result=0;
  double cl=0., cu=1.;

 if ( *ctyp==ewmaU )
   { cu = seU_q_crit_prerun_SIGMA(*l, *L0, *alpha, *hs, *sigma, *df1, *df2, *r, *qm1, *qm2, *truncate, *tail_approx, *c_error, *a_error); cl = 0.; }
   
 if ( *ctyp==ewmaUR )
   { cu = seUR_q_crit_prerun_SIGMA(*l, *L0, *alpha, *cl0, *hs, *sigma, *df1, *df2, *r, *qm1, *qm2, *truncate, *tail_approx, *c_error, *a_error); cl = *cl0; }
 if ( *ctyp==ewmaLR )
   { cl = seLR_q_crit_prerun_SIGMA(*l, *L0, *alpha, *cu0, *hs, *sigma, *df1, *df2, *r, *qm1, *qm2, *truncate, *tail_approx, *c_error, *a_error); cu = *cu0; }
 if ( *ctyp==ewma2 ) {
   if ( *ltyp==fixed ) {
     cl = se2fu_q_crit_prerun_SIGMA(*l, *L0, *alpha, *cu0, -1., *hs, *sigma, *df1, *df2, *r, *qm1, *qm2, *truncate, *tail_approx, *c_error, *a_error);
     cu = *cu0;
   }
   if ( *ltyp==unbiased )
     result = se2_q_crit_prerun_SIGMA(*l, *L0, *alpha, &cl, &cu, *hs, *sigma, *df1, *df2, *r, *qm1, *qm2, *truncate, *tail_approx, *c_error, *a_error);
 }
 
 if ( result != 0 ) warning("trouble with se2_crit called from sewma_q_crit_prerun [package spc]");
 
 c_values[0] = cl;
 c_values[1] = cu;
}
