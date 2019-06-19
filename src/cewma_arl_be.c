#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewmaU 0
#define ewmaL 1
#define ewma2 2

#define classic 0
#define transfer 1

double cewma_U_arl(double lambda, double AU, double mu0, double z0, double mu, int N);
double cewma_L_arl(double lambda, double AL, double AU, double mu0, double z0, double mu, int N);
double cewma_2_arl(double lambda, double AL, double AU, double mu0, double z0, double mu, int N);
double cewma_2_arl_rando(double lambda, double AL, double AU, double gammaL, double gammaU, double mu0, double z0, double mu, int N);

double cewma_U_arl_new(double lambda, double AU, double mu0, double z0, double mu, int N);
double cewma_L_arl_new(double lambda, double AL, double AU, double mu0, double z0, double mu, int N);
double cewma_2_arl_new(double lambda, double AL, double AU, double mu0, double z0, double mu, int N);
double cewma_2_arl_rando_new(double lambda, double AL, double AU, double gammaL, double gammaU, double mu0, double z0, double mu, int N);

void cewma_arl_be
(int *ctyp, int *mcdesign, int *rando, double *lambda, double *AL, double *AU, double *gL, double *gU, double *mu0, double *z0, double *mu, int *N, double *arl)
{ 
 *arl = -1.;
 if ( *ctyp==ewmaU && *mcdesign==classic )  *arl = cewma_U_arl(*lambda, *AU, *mu0, *z0, *mu, *N);
 if ( *ctyp==ewmaU && *mcdesign==transfer ) *arl = cewma_U_arl_new(*lambda, *AU, *mu0, *z0, *mu, *N);
 
 if ( *ctyp==ewmaL && *mcdesign==classic  ) *arl = cewma_L_arl(*lambda, *AL, *AU, *mu0, *z0, *mu, *N);
 if ( *ctyp==ewmaL && *mcdesign==transfer ) *arl = cewma_L_arl_new(*lambda, *AL, *AU, *mu0, *z0, *mu, *N);
 
 if ( *ctyp==ewma2 && *mcdesign==classic )  {
    if ( *rando==0 ) *arl = cewma_2_arl(*lambda, *AL, *AU, *mu0, *z0, *mu, *N);
    if ( *rando==1 ) *arl = cewma_2_arl_rando(*lambda, *AL, *AU, *gL, *gU, *mu0, *z0, *mu, *N);
 }
 if ( *ctyp==ewma2 && *mcdesign==transfer ) {
    if ( *rando==0 ) *arl = cewma_2_arl_new(*lambda, *AL, *AU, *mu0, *z0, *mu, *N);
    if ( *rando==1 ) *arl = cewma_2_arl_rando_new(*lambda, *AL, *AU, *gL, *gU, *mu0, *z0, *mu, *N);
 }
}
