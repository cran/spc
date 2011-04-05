#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double ewma_phat_crit(double lambda, double L0, double mu, double sigma, int n, double z0, double LSL, double USL, int N, int qm);

void ewma_phat_crit_coll
(double *lambda, double *L0, double *mu, double *sigma, int *n, double *z0, double *LSL, double *USL, int *N, int *qm, double *ucl)
{ 
 *ucl = -1.;
 *ucl = ewma_phat_crit(*lambda, *L0, *mu, *sigma, *n, *z0, *LSL, *USL, *N, *qm);
}
