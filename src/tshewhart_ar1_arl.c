#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double t_shewhart_ar1_arl(double alpha, double cS, double delta, int df, int N1, int N2, int N3, double INF, int subst);

void tshewhart_ar1_arl(double *alpha, double *cS, double *delta, int *df, int *N1, int *N2, int *N3, double *INFI, int *subst, double *arl)
{ 
 *arl = t_shewhart_ar1_arl(*alpha, *cS, *delta, *df, *N1, *N2, *N3, *INFI, *subst);
}
