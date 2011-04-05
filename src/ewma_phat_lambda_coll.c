#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double ewma_phat_lambda(double L0, double mu, double sigma, double max_l, double min_l, int n, double z0, double LSL, double USL, int qm);

void ewma_phat_lambda_coll
(double *L0, double *mu, double *sigma, double *max_l, double *min_l, int *n, double *z0, double *LSL, double *USL, int *qm, double *lambda)
{ 
 *lambda = -1.;
 *lambda = ewma_phat_lambda(*L0, *mu, *sigma, *max_l, *min_l, *n, *z0, *LSL, *USL, *qm);
}
