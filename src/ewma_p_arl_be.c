#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewmaU 0
#define ewmaL 1
#define ewma2 2

double ewma_pU_arl(double lambda, double ucl, int n, double p, double z0, int d_res, int round_mode, int mid_mode);
double ewma_pL_arl(double lambda, double lcl, int n, double p, double z0, int d_res, int round_mode, int mid_mode);
double ewma_p2_arl(double lambda, double lcl, double ucl, int n, double p, double z0, int d_res, int round_mode, int mid_mode);

void ewma_p_arl_be
(int *ctyp, double *lambda, double *lcl, double *ucl, int *n, double *p, double *z0, int *d_res, int *round_mode, int *mid_mode, double *arl)
{ 
 *arl = -1.;
 if ( *ctyp==ewmaU ) *arl = ewma_pU_arl(*lambda, *ucl, *n, *p, *z0, *d_res, *round_mode, *mid_mode);
 if ( *ctyp==ewmaL ) *arl = ewma_pL_arl(*lambda, *lcl, *n, *p, *z0, *d_res, *round_mode, *mid_mode);
 if ( *ctyp==ewma2 ) *arl = ewma_p2_arl(*lambda, *lcl, *ucl, *n, *p, *z0, *d_res, *round_mode, *mid_mode);
}

