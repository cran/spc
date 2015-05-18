#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double qf_phat(double p0, double mu, double sigma, int n, double LSL, double USL);
double qf_phat2(double p0, double mu, double sigma, int n, double LSL, double USL, int nodes);

void phat_qf
(double *x, int *n, double *mu, double *sigma, int *ctyp, double *LSL, double *USL, int *nodes, double *qf)
{ 
 *qf = -1.;
 if ( *ctyp == 0 ) *qf = qf_phat(*x, *mu, *sigma, *n, *LSL, *USL);
 if ( *ctyp == 1 ) *qf = qf_phat2(*x, *mu, *sigma, *n, *LSL, *USL, *nodes);
}
