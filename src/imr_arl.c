#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define upper 0
#define two 1

double imr_arl_case01(double M, double R, double mu, double sigma, int N, int qm);
double imr_arl_case02(double M, double R, double mu, double sigma, int N, int qm);
double imr2_arl(double M, double Rl, double Ru, double mu, double sigma, int N, int qm);
double imr2_arl_case03(double M, double Rl, double mu, double sigma, int N, int qm);

void imr_arl
(double *M, double *Rl, double *Ru, double *mu, double *sigma, int *vtyp, int *N, int *qm, double *arl)
{ 
 *arl = -1.;
 if ( *vtyp == upper ) {
   if ( *Ru >= *M ) *arl = imr_arl_case01(*M, *Ru, *mu, *sigma, *N, *qm); else *arl = imr_arl_case02(*M, *Ru, *mu, *sigma, *N, *qm);
 } else {
   if ( *Ru >= *M * 2. ) *arl = imr2_arl_case03(*M, *Rl, *mu, *sigma, *N, *qm); else *arl =  imr2_arl(*M, *Rl, *Ru, *mu, *sigma, *N, *qm); 
 }    
}
