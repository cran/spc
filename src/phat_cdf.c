#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double cdf_phat(double p, double mu, double sigma, int n, double LSL, double USL);
double cdf_phat2(double p, double mu, double sigma, int n, double LSL, double USL, int nodes);

void phat_cdf
(double *x, int *n, double *mu, double *sigma, int *ctyp, double *LSL, double *USL, int *nodes, double *cdf)
{ 
 *cdf = -1.;
 if ( *ctyp == 0 ) *cdf = cdf_phat(*x, *mu, *sigma, *n, *LSL, *USL);
 if ( *ctyp == 1 ) *cdf = cdf_phat2(*x, *mu, *sigma, *n, *LSL, *USL, *nodes);
}
