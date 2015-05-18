#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double ewma_phat_arl(double lambda, double ucl, double mu, double sigma, int n, double z0, double LSL, double USL, int N, int qm);
double ewma_phat_arl_be(double lambda, double ucl, double mu, double sigma, int n, double z0, double LSL, double USL, int N);

double ewma_phat_arl2(double lambda, double ucl, double mu, double sigma, int n, double z0, double LSL, double USL, int N, int qm, int M);
double ewma_phat_arl2_be(double lambda, double ucl, double mu, double sigma, int n, double z0, double LSL, double USL, int N);


void ewma_phat_arl_coll
(double *lambda, double *ucl, double *mu, double *sigma, int *n, double *z0, int *ctyp, double *LSL, double *USL, int *N, int *qm, int *ntyp, double *arl)
{ int M=4;
 *arl = -1.;
 if ( *ctyp == 0 ) {
   if ( *ntyp == 0 ) *arl = ewma_phat_arl(*lambda, *ucl, *mu, *sigma, *n, *z0, *LSL, *USL, *N, *qm);
   if ( *ntyp == 1 ) *arl = ewma_phat_arl_be(*lambda, *ucl, *mu, *sigma, *n, *z0, *LSL, *USL, *N);
 }
 if ( *ctyp == 1 ) {
   if ( *ntyp == 0 ) *arl = ewma_phat_arl2(*lambda, *ucl, *mu, *sigma, *n, *z0, *LSL, *USL, *N, *qm, M);
   if ( *ntyp == 1 ) *arl = ewma_phat_arl2_be(*lambda, *ucl, *mu, *sigma, *n, *z0, *LSL, *USL, *N);
 }
}
