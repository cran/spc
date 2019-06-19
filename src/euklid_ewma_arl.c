#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

double eewma_arl(int gX, int gY, int kL, int kU, double mu, double y0, int r0);

void euklid_ewma_arl
(int *gX, int *gY, int *kL, int *kU, double *mu, double *y0, int *r0, double *arl)
{ 
 *arl = -1.;
 *arl = eewma_arl(*gX, *gY, *kL, *kU, *mu, *y0, *r0);
}
