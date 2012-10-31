#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define grsr1 0
#define grsr2 1

#define Gan 0
#define Knoth 1

extern double rho0;

double xsr1_iglarl_drift(double k, double h, double zr, double hs, double delta, int m, int N, int with0);
double xsr1_iglarl_drift_wo_m(double k, double h, double zr, double hs, double delta, int *m, int N, int with0);
double xsr1_iglarlm_drift(double k, double h, double zr, double hs, int q, double delta, int N, int nmax, int with0);

void xDgrsr_arl
( double *k, double *h, double *zr, double *hs, double *delta, int *m, int *r, int *with0, int *mode, int *q, double *arl)
{
 if (*m>0)  *arl = xsr1_iglarl_drift(*k, *h, *zr, *hs, *delta, *m, *r, *with0);
 if (*m==0 && *mode==Gan) *arl = xsr1_iglarl_drift_wo_m(*k, *h, *zr, *hs, *delta, m, *r, *with0);
 if (*m==0 && *mode==Knoth) *arl = xsr1_iglarlm_drift(*k, *h, *zr, *hs, *q, *delta, *r, 10000, *with0);
}
