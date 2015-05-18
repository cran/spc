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

double xe2_iglarl_prerun_MU(double l, double c, double hs, double mu, int pn, int qm, double truncate);
double xe2_iglarl_prerun_SIGMA(double l, double c, double hs, double mu, int pn, int qm, double truncate);
double xe2_iglarl_prerun_BOTH(double l, double c, double hs, double mu, int pn, int df, int qm1, int qm2, double truncate);

double xe2_arlm_prerun_MU(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate);
double xe2_arlm_prerun_SIGMA(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate);
double xe2_arlm_prerun_BOTH(double l, double c, double hs, int q, double mu0, double mu1, int pn, int df, int mode, int nmax, int qm1, int qm2, double truncate);

void xewma_arl_prerun
( int *ctyp, double *l, double *c, double *zr, double *hs, double *mu, int *ltyp, int *q, int *size, int *df, int *mode,
  int *qm1, int *qm2, double *truncate, double *arl)
{ int nmax = 100000;

 if ( *mode == MU ) {
   if (*ctyp==ewma2 && *ltyp==fix && *q==1) *arl = xe2_iglarl_prerun_MU(*l, *c, *hs, *mu, *size, *qm1, *truncate);
   if (*ctyp==ewma2 && *ltyp==fix && *q>1)  *arl = xe2_arlm_prerun_MU(*l, *c, *hs, *q, 0., *mu, *size, *ltyp, nmax, *qm1, *truncate);
   if (*ctyp==ewma2 && *ltyp>fix)           *arl = xe2_arlm_prerun_MU(*l, *c, *hs, *q, 0., *mu, *size, *ltyp, nmax, *qm1, *truncate);
 }
 if ( *mode == SIGMA ) {
   if (*ctyp==ewma2 && *ltyp==fix && *q==1) *arl = xe2_iglarl_prerun_SIGMA(*l, *c, *hs, *mu, *size, *qm2, *truncate);
   if (*ctyp==ewma2 && *ltyp==fix && *q>1)  *arl = xe2_arlm_prerun_SIGMA(*l, *c, *hs, *q, 0., *mu, *size, *ltyp, nmax, *qm2, *truncate);
   if (*ctyp==ewma2 && *ltyp>fix)           *arl = xe2_arlm_prerun_SIGMA(*l, *c, *hs, *q, 0., *mu, *size, *ltyp, nmax, *qm2, *truncate);
 }
  if ( *mode == BOTH ) {
   if (*ctyp==ewma2 && *ltyp==fix && *q==1) *arl = xe2_iglarl_prerun_BOTH(*l, *c, *hs, *mu, *size, *df, *qm1, *qm2, *truncate);
   if (*ctyp==ewma2 && *ltyp==fix && *q>1)  *arl = xe2_arlm_prerun_BOTH(*l, *c, *hs, *q, 0., *mu, *size, *df, *ltyp, nmax, *qm1, *qm2, *truncate);
   if (*ctyp==ewma2 && *ltyp>fix)           *arl = xe2_arlm_prerun_BOTH(*l, *c, *hs, *q, 0., *mu, *size, *df, *ltyp, nmax, *qm1, *qm2, *truncate);
 }
}
