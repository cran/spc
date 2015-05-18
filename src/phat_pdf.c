#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double pdf_phat (double p, double mu, double sigma, int n, double LSL, double USL);
double pdf_phat2(double p, double mu, double sigma, int n, double LSL, double USL, int nodes);

void phat_pdf
(double *x, int *n, double *mu, double *sigma, int *ctyp, double *LSL, double *USL, int *nodes, double *pdf)
{ 
 *pdf = -1.;
 if ( *ctyp == 0 ) *pdf = pdf_phat (*x, *mu, *sigma, *n, *LSL, *USL);
 if ( *ctyp == 1 ) *pdf = pdf_phat2(*x, *mu, *sigma, *n, *LSL, *USL, *nodes);
}
