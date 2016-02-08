#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double x_shewhart_ar1_arl(double alpha, double cS, double delta, int N1, int N2);

void xshewhart_ar1_arl(double *alpha, double *cS, double *delta, int *N1, int *N2, double *arl)
{ 
 *arl = x_shewhart_ar1_arl(*alpha, *cS, *delta, *N1, *N2);
}
