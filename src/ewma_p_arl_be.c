#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double ewma_p_arl(double lambda, double ucl, int n, double p, double z0, int d_res, int round_mode, int mid_mode);

void ewma_p_arl_be
(double *lambda, double *ucl, int *n, double *p, double *z0, int *d_res, int *round_mode, int *mid_mode, double *arl)
{ 
 *arl = -1.;
 *arl = ewma_p_arl(*lambda, *ucl, *n, *p, *z0, *d_res, *round_mode, *mid_mode);
}
