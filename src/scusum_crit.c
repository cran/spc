#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define cusumU 0
#define cusumL 1
#define cusum2 2

double scU_crit(double refk, double L0, double hs, double sigma, int df, int N, int qm);
double scL_crit(double refk, double L0, double hs, double sigma, int df, int N, int qm);
int sc2_crit_unbiased(double refkl, double refku, double L0, double *hl, double *hu, double hsl, double hsu, double sigma, int df, int N, int qm);

void scusum_crit(int *ctyp, double *k, double *L0, double *hs, double *sigma, int *df, int *ltyp, double *k2, double *hs2, int *r, int *qm, double *h)
{ int result=0;
  double hl=0., hu=0.;
  
  if ( *ctyp==cusumU ) *h = scU_crit(*k, *L0, *hs, *sigma, *df, *r, *qm);
  if ( *ctyp==cusumL ) *h = scL_crit(*k, *L0, *hs, *sigma, *df, *r, *qm);
  if ( *ctyp==cusum2 ) {
    result = sc2_crit_unbiased(*k2, *k, *L0, &hl, &hu, *hs2, *hs, *sigma, *df, *r, *qm);
    if ( result != 0 ) warning("trouble with sc2_crit_unbiased called from scusum_crit [package spc]");
    h[0] = hl;
    h[1] = hu;
  }
}
