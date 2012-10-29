#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double ewma_phat_arl(double lambda, double ucl, double mu, double sigma, int n, double z0, double LSL, double USL, int N, int qm);

void ewma_phat_arl_coll
(double *lambda, double *ucl, double *mu, double *sigma, int *n, double *z0, double *LSL, double *USL, int *N, int *qm, double *arl)
{ 
 *arl = -1.;
 *arl = ewma_phat_arl(*lambda, *ucl, *mu, *sigma, *n, *z0, *LSL, *USL, *N, *qm);
}
