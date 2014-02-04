#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define grsr1 0
#define grsr2 1

extern double rho0;

double xsr1_crit(double k, double L0, double zr, double hs, double m0, int N, int MPT);

void xgrsr_crit(double *k, double *L0, double *zr, double *hs, double *mu0, int *r, int *MPT, double *h)
{
  *h = xsr1_crit(*k, *L0, *zr, *hs, *mu0, *r, *MPT);
}
