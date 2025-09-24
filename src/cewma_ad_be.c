#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define ewmaU 0
#define ewmaL 1
#define ewma2 2

#define classic 0
#define transfer 1

double cewma_2_ad(double lambda, double AL, double AU, double mu0, double mu, int N);

double cewma_2_ad_new(double lambda, double AL, double AU, double mu0, double mu, int N, int OLD);

void cewma_ad_be(int* ctyp, int* mcdesign, int* rando, double* lambda, double* AL, double* AU, double* gL, double* gU, double* mu0, double* mu, int* N, int* OLD, double* ad)
{
    *ad = -1.;

    if (*ctyp == ewma2 && *mcdesign == classic) {
        if (*rando == 0)
            *ad = cewma_2_ad(*lambda, *AL, *AU, *mu0, *mu, *N);
    }

    if (*ctyp == ewma2 && *mcdesign == transfer) {
        if (*rando == 0)
            *ad = cewma_2_ad_new(*lambda, *AL, *AU, *mu0, *mu, *N, *OLD);
    }
}
