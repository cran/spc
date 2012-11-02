#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define cusum1 0
#define cusum2 1
#define cusumC 2

#define Gan 0
#define Knoth 1

extern double rho0;

double xc1_iglarl_drift(double k, double h, double hs, double delta, int m, int N, int with0);
double xc1_iglarl_drift_wo_m(double k, double h, double hs, double delta, int *m, int N, int with0);
double xc1_iglarlm_drift(double k, double h, double hs, int q, double delta, int N, int nmax, int with0);

void xDcusum_arl
( int *ctyp, double *k, double *h, double *hs, double *delta, int *m, int *r, int *with0, int *mode, int *q, double *arl)
{
 if (*ctyp==cusum1 && *m>0)  *arl = xc1_iglarl_drift(*k, *h, *hs, *delta, *m, *r, *with0);
 if (*ctyp==cusum1 && *m==0 && *mode==Gan)   *arl = xc1_iglarl_drift_wo_m(*k, *h, *hs, *delta, m, *r, *with0);
 if (*ctyp==cusum1 && *m==0 && *mode==Knoth) *arl = xc1_iglarlm_drift(*k, *h, *hs, *q, *delta, *r, 10000, *with0);
}
