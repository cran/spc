#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double mxewma_crit(double lambda, double L0, int p, double hs, int N);

void mewma_crit(double *l, double *L0, int *p, double *hs, int *r, double *h)
{
  *h = mxewma_crit(*l, *L0, *p, *hs, *r);
}
