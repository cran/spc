#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewma1 0
#define ewma2 1
#define fix 0
#define vacl 1
#define fir 2
#define both 3
#define steiner 4
#define stat 5

#define MU 0
#define SIGMA 1
#define BOTH 2

double xe2_Wq_prerun_MU_deluxe(double l, double c, double p, double hs, double mu, int pn, int nmax, int qm, double truncate);
double xe2_Wq_prerun_SIGMA_deluxe(double l, double c, double p, double hs, double mu, int pn, int nmax, int qm, double truncate);
double xe2_Wq_prerun_BOTH_deluxe(double l, double c, double p, double hs, double mu, int pn, int df, int nmax, int qm1, int qm2, double truncate);

double xe2_Wqm_prerun_MU_deluxe(double l, double c, double p, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate);
double xe2_Wqm_prerun_SIGMA_deluxe(double l, double c, double p, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate);
double xe2_Wqm_prerun_BOTH_deluxe(double l, double c, double p, double hs, int q, double mu0, double mu1, int pn, int df, int mode, int nmax, int qm1, int qm2, double truncate);

void xewma_q_prerun
( int *ctyp, double *l, double *c, double *p, double *zr, double *hs, double *mu, int *ltyp, int *q, int *size, int *df, int *mode, int *qm1, int *qm2, double *truncate, double *tq)
{ int nmax=1000000;

 if (*mode == MU) {
   if ( *ctyp==ewma2 && *ltyp==fix && *q==1 )  *tq =  xe2_Wq_prerun_MU_deluxe(*l, *c, *p, *hs, *mu, *size, nmax, *qm1, *truncate);
   if ( *ctyp==ewma2 && *ltyp==fix && *q>1 )   *tq = xe2_Wqm_prerun_MU_deluxe(*l, *c, *p, *hs, *q, 0., *mu, *size, *ltyp, nmax, *qm1, *truncate);
   if ( *ctyp==ewma2 && *ltyp>fix )            *tq = xe2_Wqm_prerun_MU_deluxe(*l, *c, *p, *hs, *q, 0., *mu, *size, *ltyp, nmax, *qm1, *truncate);
 }
 if (*mode == SIGMA) {
   if ( *ctyp==ewma2 && *ltyp==fix && *q==1 )  *tq =  xe2_Wq_prerun_SIGMA_deluxe(*l, *c, *p, *hs, *mu, *size, nmax, *qm2, *truncate);
   if ( *ctyp==ewma2 && *ltyp==fix && *q>1 )   *tq = xe2_Wqm_prerun_SIGMA_deluxe(*l, *c, *p, *hs, *q, 0., *mu, *size, *ltyp, nmax, *qm2, *truncate);
   if ( *ctyp==ewma2 && *ltyp>fix )            *tq = xe2_Wqm_prerun_SIGMA_deluxe(*l, *c, *p, *hs, *q, 0., *mu, *size, *ltyp, nmax, *qm2, *truncate);
 }
 if (*mode == BOTH) {
   if ( *ctyp==ewma2 && *ltyp==fix && *q==1 )  *tq =  xe2_Wq_prerun_BOTH_deluxe(*l, *c, *p, *hs, *mu, *size, *df, nmax, *qm1, *qm2, *truncate);
   if ( *ctyp==ewma2 && *ltyp==fix && *q>1 )   *tq = xe2_Wqm_prerun_BOTH_deluxe(*l, *c, *p, *hs, *q, 0., *mu, *size, *df, *ltyp, nmax, *qm1, *qm2, *truncate);
   if ( *ctyp==ewma2 && *ltyp>fix )            *tq = xe2_Wqm_prerun_BOTH_deluxe(*l, *c, *p, *hs, *q, 0., *mu, *size, *df, *ltyp, nmax, *qm1, *qm2, *truncate);
 }
}
