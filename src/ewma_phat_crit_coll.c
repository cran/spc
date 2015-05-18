#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double ewma_phat_crit(double lambda, double L0, double mu, double sigma, int n, double z0, double LSL, double USL, int N, int qm);
double ewma_phat_crit2(double lambda, double L0, double mu, double sigma, int n, double z0, double LSL, double USL, int N, int qm, int M);

void ewma_phat_crit_coll
(double *lambda, double *L0, double *mu, double *sigma, int *n, double *z0, int *ctyp, double *LSL, double *USL, int *N, int *qm, double *ucl)
{ int M=4;
 *ucl = -1.; 
 if ( *ctyp == 0 ) *ucl = ewma_phat_crit(*lambda, *L0, *mu, *sigma, *n, *z0, *LSL, *USL, *N, *qm);
 if ( *ctyp == 1 ) *ucl = ewma_phat_crit2(*lambda, *L0, *mu, *sigma, *n, *z0, *LSL, *USL, *N, *qm, M);
}
