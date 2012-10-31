#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define cusum1 0
#define cusum2 1

double xc1_Wq(double k, double h, double p, double hs, double mu, int N, int nmax);

void xcusum_q(int *ctyp, double *k, double *h, double *p, double *hs, double *mu, int *r, double *q)
{
 if (*ctyp==cusum1) *q = xc1_Wq(*k, *h, *p, *hs, *mu, *r, 10000);
}
