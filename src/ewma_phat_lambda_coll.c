#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double ewma_phat_lambda(double L0, double mu, double sigma, double max_l, double min_l, int n, double z0, double LSL, double USL, int qm);
double ewma_phat_lambda2(double L0, double mu, double sigma, double max_l, double min_l, int n, double z0, double LSL, double USL, int qm, int M);

void ewma_phat_lambda_coll
(double *L0, double *mu, double *sigma, int *ctyp, double *max_l, double *min_l, int *n, double *z0, double *LSL, double *USL, int *qm, double *lambda)
{ int M=4; 
 *lambda = -1.;
 if ( *ctyp == 0 ) *lambda = ewma_phat_lambda(*L0, *mu, *sigma, *max_l, *min_l, *n, *z0, *LSL, *USL, *qm);
 if ( *ctyp == 1 ) *lambda = ewma_phat_lambda2(*L0, *mu, *sigma, *max_l, *min_l, *n, *z0, *LSL, *USL, *qm, M);
}
