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
#define sven 5
#define fink 6
#define waldmann 7
#define collocation 8

#define Gan 0
#define Knoth 1
#define Waldm 2

extern double rho0;

double xe1_iglarl_drift(double l, double c, double zr, double hs, double delta, int m, int N, int with0);
double xe1_iglarl_drift_wo_m(double l, double c, double zr, double hs, double delta, int *m, int N, int with0);
double xe1_iglarlm_drift(double l, double c, double zr, double hs, int q, double delta, int N, int nmax, int with0);
double xe2_iglarl_drift(double l, double c, double hs, double delta, int m, int N, int with0);
double xe2_iglarl_drift_wo_m(double l, double c, double hs, double delta, int *m, int N, int with0);
double xe2_iglarlm_drift(double l, double c, double hs, int q, double delta, int N, int nmax, int with0);
double xe2_Warl_drift(double l, double c, double hs, double delta, int N, int nmax, int with0);

void xDewma_arl
( int *ctyp, double *l, double *c, double *zr, double *hs, double *delta, int *ltyp, int *m, int *r, int *with0, int *mode, int *q, double *arl)
{
 if (*ctyp==ewma1 && *m>0)
   *arl = xe1_iglarl_drift(*l,*c,*zr,*hs,*delta,*m,*r,*with0);
 if (*ctyp==ewma1 && *m==0 && *mode==Gan)
   *arl = xe1_iglarl_drift_wo_m(*l,*c,*zr,*hs,*delta,m,*r,*with0);
 if (*ctyp==ewma1 && *m==0 && *mode==Knoth)
   *arl = xe1_iglarlm_drift(*l,*c,*zr,*hs,*q,*delta,*r,10000,*with0);

 if (*ctyp==ewma2 && *m>0)
   *arl = xe2_iglarl_drift(*l,*c,*hs,*delta,*m,*r,*with0);
 if (*ctyp==ewma2 && *m==0 && *mode==Gan)
   *arl = xe2_iglarl_drift_wo_m(*l,*c,*hs,*delta,m,*r,*with0);
 if (*ctyp==ewma2 && *m==0 && *mode==Knoth)
   *arl = xe2_iglarlm_drift(*l,*c,*hs,*q,*delta,*r,10000,*with0);
 if (*ctyp==ewma2 && *m==0 && *mode==Waldm)
   *arl = xe2_Warl_drift(*l,*c,*hs,*delta,*r,10000,*with0);
}
