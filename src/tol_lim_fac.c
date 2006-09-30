#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define WW 0
#define exact 1

double kww(int n, double p, double a);
double tl_factor (int n, double p, double a, int m);

void tol_lim_fac
( int *n, double *p, double *a, int *mtype, int *m, double *tlf )
{
  if (*mtype==WW) *tlf = kww(*n,*p,*a);
  else *tlf = tl_factor(*n,*p,*a,*m);        
}
