#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define upper 0
#define lower 1
#define two 2

#define sym 0
#define unb 1

#define classic 0
#define transfer 1

double cewma_U_crit(double lambda, double L0, double mu0, double z0, int N, int jmax);
double cewma_L_crit(double lambda, double L0, double AU, double mu0, double z0, int N, int jmax);
double cewma_2_crit_sym(double lambda, double L0, double mu0, double z0, int N, int jmax);
int cewma_2_crit_unb(double lambda, double L0, double mu0, double z0, int N, int jmax, double* AL, double* AU);

double cewma_U_crit_new(double lambda, double L0, double mu0, double z0, int N, int jmax, int OLD);
double cewma_L_crit_new(double lambda, double L0, double AU, double mu0, double z0, int N, int jmax, int OLD);
double cewma_2_crit_sym_new(double lambda, double L0, double mu0, double z0, int N, int jmax, int OLD);
int cewma_2_crit_unb_new(double lambda, double L0, double mu0, double z0, int N, int jmax, double* AL, double* AU, int OLD);
int cewma_2_crit_unb_rando_new(double lambda, double L0, double mu0, double z0, int N, int jmax, double* AL, double* AU, double* gL, double* gU, int OLD);

void cewma_crit_be(int* ctyp, int* design, int* mcdesign, int* rando, double* lambda, double* L0, double* AU, double* mu0, double* z0, int* N, int* jmax, int* OLD, double* AA)
{
    int result = 0;
    double lAL, lAU, lgL, lgU;
    *AA = -1.;
    if (*ctyp == upper) {
        if (*mcdesign == classic)
            *AA = cewma_U_crit(*lambda, *L0, *mu0, *z0, *N, *jmax);
        if (*mcdesign == transfer)
            *AA = cewma_U_crit_new(*lambda, *L0, *mu0, *z0, *N, *jmax, *OLD);
    }
    if (*ctyp == lower) {
        if (*mcdesign == classic)
            *AA = cewma_L_crit(*lambda, *L0, *AU, *mu0, *z0, *N, *jmax);
        if (*mcdesign == transfer)
            *AA = cewma_L_crit_new(*lambda, *L0, *AU, *mu0, *z0, *N, *jmax, *OLD);
    }
    if (*ctyp == two) {
        if (*design == sym) {
            if (*mcdesign == classic)
                *AA = cewma_2_crit_sym(*lambda, *L0, *mu0, *z0, *N, *jmax);
            if (*mcdesign == transfer)
                *AA = cewma_2_crit_sym_new(*lambda, *L0, *mu0, *z0, *N, *jmax, *OLD);
        }
        if (*design == unb) {
            if (*rando == 0) {
                if (*mcdesign == classic)
                    result = cewma_2_crit_unb(*lambda, *L0, *mu0, *z0, *N, *jmax, &lAL, &lAU);
                if (*mcdesign == transfer)
                    result = cewma_2_crit_unb_new(*lambda, *L0, *mu0, *z0, *N, *jmax, &lAL, &lAU, *OLD);
                AA[0] = lAL;
                AA[1] = lAU;
            }
            if (*rando == 1) {
                /*if ( *mcdesign==classic )  result = cewma_2_crit_unb(*lambda, *L0, *mu0, *z0, *N, *jmax, &lAL, &lAU);*/
                if (*mcdesign == transfer)
                    result = cewma_2_crit_unb_rando_new(*lambda, *L0, *mu0, *z0, *N, *jmax, &lAL, &lAU, &lgL, &lgU, *OLD);
                AA[0] = lAL;
                AA[1] = lAU;
                AA[2] = lgL;
                AA[3] = lgU;
            }
        }
    }
    if (result != 0)
        warning("something went wrong with cewma_2_crit_unb_*");
}
