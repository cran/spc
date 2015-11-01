#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double mxewma_ad (double lambda, double ce, int p, double delta, int N, int qm2, int psi_type, double hs, int qtype, int qm0, int qm1);

void mewma_ad(double *l, double *c, int *p, double *delta, int *r, int *qm2, int *ptype, double *hs, int *qtype, int *qm0, int *qm1, double *ad)
{ 
  *ad = mxewma_ad(*l, *c, *p, *delta, *r, *qm2, *ptype, *hs, *qtype, *qm0, *qm1);
}
