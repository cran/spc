#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h> 

#define LOG 0
#define TAIL 1

#define cusum1 0
#define cusum2 1
#define cusumC 2
#define ewma1 0
#define ewma2 1
#define fix 0
#define vacl 1
#define fir 2
#define both 3
#define steiner 4
#define stat 5
#define fink 6

#define FINALeps 1e-12
#define lmEPS 1e-4

#define IDENTITY 0
#define SIN 1
#define SINH 2
#define TAN 3

/*** export ***/


/* CUSUM */

double xc_crit(int ctyp, double k, double L0, double hs, double m0, int N);

/* one-sided CUSUM */
double xc1_iglarl(double k, double h, double hs, double mu, int N);
double xc1_be_arl(double k, double h, double hs, double mu, int N);
double xc1_beL_arl(double k, double h, double hs, double mu, int N);
double xc1_beT_arl(double k, double h, double hs, double mu, int N);
double xc1_iglad (double k, double h, double mu0, double mu1, int N);
double xc1_iglarl_drift(double k, double h, double hs, double delta, int m, int N, int with0);
double xc1_iglarl_drift_wo_m(double k, double h, double hs, double delta, int *m, int N, int with0);
double xc1_iglarlm_drift(double k, double h, double hs, int q, double delta, int N, int nmax, int with0);

double xtc1_iglarl(double k, double h, double hs, int df, double mu, int N, int subst);

double xc1_Wq(double k, double h, double p, double hs, double mu, int N, int nmax);
double xc1_sf(double k, double h, double hs, double mu, int N, int nmax, double *p0);
double xc1_arlm(double k, double h, double hs, int q, double mu0, double mu1, int N, int nmax);
double xc1_arlm_hom(double k, double h, double hs, int q, double mu0, double mu1, int N, double *ced);

/* classical two-sided (2 charts) CUSUM */
double xc2_iglarl(double k, double h, double hs, double mu, int N);
double xc2_be_arl(double k, double h, double hs1, double hs2, double mu, int N);
double xc2_be_arlm(double k, double h, double hs1, double hs2, int q, double mu0, double mu1, int N, double *ced);
double xc2_iglad (double k, double h, double mu0, double mu1, int N);
double xc2_iglarl_drift(double k, double h, double hs, double delta, int m, int N, int drift0); /* it is not accurate */
double xc2_iglarl_drift_wo_m(double k, double h, double hs, double delta, int *m, int N, int drift0); /* it is not accurate */

double xtc2_iglarl(double k, double h, double hs, int df, double mu, int N, int subst);

/* Crosier's two-sided CUSUM */
double xcC_iglarl(double k, double h, double hs, double mu, int N);
double xcC_iglad (double k, double h, double mu0, double mu1, int N);


/* variance charts */
double scU_iglarl_v1(double refk, double h, double hs, double sigma, int df, int N, int qm);
double scU_iglarl_v2(double refk, double h, double hs, double sigma, int df, int N, int qm);
double scL_iglarl_v2(double refk, double h, double hs, double sigma, int df, int N, int qm);
double sc2_iglarl_v2(double refkl, double refku, double hl, double hu, double hsl, double hsu, double sigma, int df, int N, int qm);

double scU_crit(double refk, double L0, double hs, double sigma, int df, int N, int qm);
double scL_crit(double refk, double L0, double hs, double sigma, int df, int N, int qm);
double scU_fl_crit(double refkl, double refku, double hl, double L0, double hsl, double hsu, double sigma, int df, int N, int qm);
double scL_fu_crit(double refkl, double refku, double hu, double L0, double hsl, double hsu, double sigma, int df, int N, int qm);
int sc2_crit_unbiased(double refkl, double refku, double L0, double *hl, double *hu, double hsl, double hsu, double sigma, int df, int N, int qm);
/*double  sc2_eqtails(double refkl, double refku, double L0, double *hl, double *hu, double hsl, double hsu, double sigma, int df, int N, int qm);*/

/* CUSUM-Shewhart combo */
double scs_U_iglarl_v1(double refk, double h, double hs, double cS, double sigma, int df, int N, int qm);


/* Shiryaev-Roberts (only the one-sided version is implemented) */

double xsr1_crit(double k, double L0, double zr, double hs, double m0, int N, int MPT);

double xsr1_iglarl(double k, double h, double zr, double hs, double mu, int N, int MPT);
double xsr1_iglad(double k, double h, double zr, double mu0, double mu1, int N, int MPT);
double xsr1_arlm(double k, double h, double zr, double hs, int q, double mu0, double mu1, int N, int nmax, int MPT);
double xsr1_arlm_hom(double k, double h, double zr, double hs, int q, double mu0, double mu1, int N, int MPT, double *ced);
double xsr1_iglarl_drift(double k, double h, double zr, double hs, double delta, int m, int N, int with0);
double xsr1_iglarl_drift_wo_m(double k, double h, double zr, double hs, double delta, int *m, int N, int with0);
double xsr1_iglarlm_drift(double k, double h, double zr, double hs, int q, double delta, int N, int nmax, int with0);


/* EWMA */

double xe_crit(int ctyp, double l, double L0, double zr, double hs, double m0, int ltyp, int N, double c0);
double xe_q_crit(int ctyp, double l, int L0, double alpha, double zr, double hs, double m0, int ltyp, int N, double c_error, double a_error);

/* one-sided EWMA */
double xe1_iglarl(double l, double c, double zr, double hs, double mu, int N);
double xe1_iglad (double l, double c, double zr, double mu0, double mu1, int N);
double xe1_arlm(double l, double c, double zr, double hs, int q, double mu0, double mu1, int mode, int N, int nmax);
double xe1_arlm_hom(double l, double c, double zr, double hs, int q, double mu0, double mu1, int N, double *ced);
double xe1_Warl(double l, double c, double zr, double hs, double mu, int N, int nmax);
double xe1_Wq(double l, double c, double p, double zr, double hs, double mu, int N, int nmax);
double xe1_sf(double l, double c, double zr, double hs, double mu, int N, int nmax, double *p0);
double xe1_sfm(double l, double c, double zr, double hs, int q, double mu0, double mu1, int mode, int N, int nmax, double *p0);
double xe1_Wqm(double l, double c, double p, double zr, double hs, int q, double mu0, double mu1, int mode, int N, int nmax);
double xe1_iglarl_drift(double l, double c, double zr, double hs, double delta, int m, int N, int with0);
double xe1_iglarl_drift_wo_m(double l, double c, double zr, double hs, double delta, int *m, int N, int with0);
double xe1_iglarlm_drift(double l, double c, double zr, double hs, int q, double delta, int N, int nmax, int with0);

double xlimit1_arlm(double c, double zr, int q, double mu0, double mu1, int N, int nmax);

/* two-sided EWMA */
double xe2_iglarl(double l, double c, double hs, double mu, int N);
double xe2_iglarl_f(double l, double c, double mu, int N, double *g, double *w, double *z);
double xe2_iglad (double l, double c, double mu0, double mu1, int N);
double xe2_igladc(double l, double c, double mu0, double mu1, double z0, int N);
double xe2_arlm(double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax);
double xe2_arlmc(double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax);
int xe2_arlm_special(double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax, double *pair);
double xe2_arlm_hom(double l, double c, double hs, int q, double mu0, double mu1, int N, double *ced);
double xe2_Wq(double l, double c, double p, double hs, double mu, int N, int nmax);
double xe2_sf(double l, double c, double hs, double mu, int N, int nmax, double *p0);
double xe2_sfm(double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax, double *p0);
double xe2_Wqm(double l, double c, double p, double hs, int q, double mu0, double mu1, int mode, int N, int nmax);

double xe2_Warl(double l, double c, double hs, double mu, int N, int nmax); /* Waldmann's ARL procedure */
double xe2_Carl(double l, double c, double hs, double mu, int N, int qm); /* collocation */

double xe2_iglarl_drift(double l, double c, double hs, double delta, int m, int N, int with0);
double xe2_iglarl_drift_wo_m(double l, double c, double hs, double delta, int *m, int N, int with0);
double xe2_iglarlm_drift(double l, double c, double hs, int q, double delta, int N, int nmax, int with0);
double xe2_Warl_drift(double l, double c, double hs, double delta, int N, int nmax, int with0);

/* functions based on Srivastava & Wu (1997) */
double xe2_SrWu_crit(double l, double L0);
double xe2_SrWu_arl(double l, double c, double mu);
double xe2_SrWu_arl_full(double l, double c, double mu);
double xe2_SrWu_lambda(double delta, double L0);

/* t distribution */
double xte2_iglarl(double l, double c, double hs, int df, double mu, int N, int subst);
double xte2_iglad (double l, double c, int df, double mu0, double mu1, int N, int subst);
double xte2_igladc(double l, double c, int df, double mu0, double mu1, double z0, int N, int subst);
double xte2_arlm(double l, double c, double hs, int df, int q, double mu0, double mu1, int mode, int N, int nmax, int subst);
double xte2_arlm_hom(double l, double c, double hs, int df, int q, double mu0, double mu1, int N, double *ced, int subst);
double xte2_Wq(double l, double c, double p, double hs, int df, double mu, int N, int nmax, int subst);
double xte2_sf(double l, double c, double hs, int df, double mu, int N, int nmax, double *p0, int subst);
double xte2_sfm(double l, double c, double hs, int df, int q, double mu0, double mu1, int mode, int N, int nmax, double *p0, int subst);
double xte2_Wqm(double l, double c, double p, double hs, int df, int q, double mu0, double mu1, int mode, int N, int nmax, int subst);

double xte1_iglarl(double l, double c, double zr, double hs, int df, double mu, int N, int subst);


/* incorporate pre-run uncertainty */
double xe2_iglarl_prerun_MU(double l, double c, double hs, double mu, int pn, int qm, double truncate);
double xe2_iglarl_prerun_SIGMA(double l, double c, double hs, double mu, int pn, int qm, double truncate);
double xe2_iglarl_prerun_BOTH(double l, double c, double hs, double mu, int pn, int df, int qm1, int qm2, double truncate);

double xe2_arlm_prerun_MU(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate);
double xe2_arlm_prerun_SIGMA(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate);
double xe2_arlm_prerun_BOTH(double l, double c, double hs, int q, double mu0, double mu1, int pn, int df, int mode, int nmax, int qm1, int qm2, double truncate);

double xe2_sf_deluxe(double l, double c, double hs, double mu, int N, int nmax, double BOUND, double *p0, int *nstop, double *rho);
double xe2_sf_prerun_MU_deluxe(double l, double c, double hs, double mu, int pn, int nmax, int qm, double truncate, double BOUND, double *p0);
double xe2_sf_prerun_MU(double l, double c, double hs, double mu, int pn, int nmax, int qm, double truncate, double *p0);
double xe2_sf_prerun_SIGMA_deluxe(double l, double c, double hs, double mu, int pn, int nmax, int qm, double truncate, double BOUND, double *p0);
double xe2_sf_prerun_SIGMA(double l, double c, double hs, double mu, int pn, int nmax, int qm, double truncate, double *p0);
double xe2_sf_prerun_BOTH_deluxe(double l, double c, double hs, double mu, int pn, int df, int nmax, int qm1, int qm2, double truncate, double BOUND, double *p0);
double xe2_sf_prerun_BOTH(double l, double c, double hs, double mu, int pn, int df, int nmax, int qm1, int qm2, double truncate, double *p0);

double xe2_sfm_simple(double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax, double *p0);
double xe2_sfm_deluxe(double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax, double BOUND, double *p0, int *nstop, double *rho);
double xe2_sfm_prerun_MU_deluxe(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate, double BOUND, double *p0);
double xe2_sfm_prerun_MU(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate, double *p0);
double xe2_sfm_prerun_SIGMA_deluxe(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate, double BOUND, double *p0);
double xe2_sfm_prerun_SIGMA(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate, double *p0);
double xe2_sfm_prerun_BOTH_deluxe(double l, double c, double hs, int q, double mu0, double mu1, int pn, int df, int mode, int nmax, int qm1, int qm2, double truncate, double BOUND, double *p0);
double xe2_sfm_prerun_BOTH(double l, double c, double hs, int q, double mu0, double mu1, int pn, int df, int mode, int nmax, int qm1, int qm2, double truncate, double *p0);

double xe2_Wq_prerun_MU_deluxe(double l, double c, double p, double hs, double mu, int pn, int nmax, int qm, double truncate, double BOUND);
double xe2_Wq_prerun_SIGMA_deluxe(double l, double c, double p, double hs, double mu, int pn, int nmax, int qm, double truncate, double BOUND);
double xe2_Wq_prerun_BOTH_deluxe(double l, double c, double p, double hs, double mu, int pn, int df, int nmax, int qm1, int qm2, double truncate, double BOUND);

double xe2_Wqm_prerun_MU_deluxe(double l, double c, double p, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate, double BOUND);
double xe2_Wqm_prerun_SIGMA_deluxe(double l, double c, double p, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate, double BOUND);
double xe2_Wqm_prerun_BOTH_deluxe(double l, double c, double p, double hs, int q, double mu0, double mu1, int pn, int df, int mode, int nmax, int qm1, int qm2, double truncate, double BOUND);


/* EWMA residual charts */

double xe2_iglarl_RES(double l, double c, double hs, double mu, int N, double alpha, int df);
double seU_iglarl_RES(double l, double cu, double hs, double sigma, int df, int N, int qm, double alpha, double mu);  
double xseU_arl_RES(double lx, double ls, double cx, double cs, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm, double alpha);   
double xseU_mu_before_sigma_RES(double lx, double ls, double cx, double cs, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm, double alpha, int vice_versa);


/* modified Shewhart charts for dependent data */
double x_shewhart_ar1_arl(double alpha, double cS, double mu, int N1, int N2);
double t_shewhart_ar1_arl(double alpha, double cS, double delta, int df, int N1, int N2, int N3, double INF, int subst);
   

/* variance charts */
double seU_iglarl(double l, double cu, double hs, double sigma, int df, int N, int qm);
double se2_iglarl(double l, double cl, double cu, double hs, double sigma, int df, int N, int qm);
double seUR_iglarl(double l, double cl, double cu, double hs, double sigma, int df, int N, int qm);
double seLR_iglarl(double l, double cl, double cu, double hs, double sigma, int df, int N, int qm);

double stdeU_iglarl(double l, double cu, double hs, double sigma, int df, int N, int qm);
double stde2_iglarl(double l, double cl, double cu, double hs, double sigma, int df, int N, int qm);
double stdeUR_iglarl(double l, double cl, double cu, double hs, double sigma, int df, int N, int qm);
double stdeLR_iglarl(double l, double cl, double cu, double hs, double sigma, int df, int N, int qm);

double lns2ewmaU_arl_igl(double l, double cl, double cu, double hs, double sigma, int df, int N);
double lns2ewma2_arl_igl(double l, double cl, double cu, double hs, double sigma, int df, int N);

double seU_crit(double l, double L0, double hs, double sigma, int df, int N, int qm);
double se2lu_crit(double l, double L0, double cl, double hs, double sigma, int df, int N, int qm);
double se2fu_crit(double l, double L0, double cu, double hs, double sigma, int df, int N, int qm);
int se2_crit_unbiased(double l, double L0, double *cl, double *cu, double hs, double sigma, int df, int N, int qm);
int se2_crit_eqtails(double l, double L0, double *cl, double *cu, double hs, double sigma, int df, double ur, int N, int qm);
double se2_crit_sym(double l, double L0, double hs, double sigma, int df, int N, int qm);
double seUR_crit(double l, double L0, double cl, double hs, double sigma, int df, int N, int qm);
double seLR_crit(double l, double L0, double cu, double hs, double sigma, int df, int N, int qm);

double stdeU_crit(double l, double L0, double hs, double sigma, int df, int N, int qm);
double stde2lu_crit(double l, double L0, double cl, double hs, double sigma, int df, int N, int qm);
double stde2fu_crit(double l, double L0, double cu, double hs, double sigma, int df, int N, int qm);
int stde2_crit_unbiased(double l, double L0, double *cl, double *cu, double hs, double sigma, int df, int N, int qm);
int stde2_crit_eqtails(double l, double L0, double *cl, double *cu, double hs, double sigma, int df, double ur, int N, int qm);
double stde2_crit_sym(double l, double L0, double hs, double sigma, int df, int N, int qm);
double stdeUR_crit(double l, double L0, double cl, double hs, double sigma, int df, int N, int qm);
double stdeLR_crit(double l, double L0, double cu, double hs, double sigma, int df, int N, int qm);

double lns2ewmaU_crit(double l, double L0, double cl, double hs, double sigma, int df, int N);
double lns2ewma2_crit_cufix(double l, double cu, double L0, double hs, double sigma, int df, int N);
double lns2ewma2_crit_sym(double l, double L0, double hs, double sigma, int df, int N);
int lns2ewma2_crit_unbiased(double l, double L0, double *cl, double *cu, double hs, double sigma, int df, int N);

double seU_sf(double l, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0);
double seU_sf_deluxe(double l, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0, int *nstop, double *rho);
double se2_sf(double l, double cl, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0);
double se2_sf_deluxe(double l, double cl, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0, int *nstop, double *rho);
double seUR_sf(double l, double cl, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0);
double seUR_sf_deluxe(double l, double cl, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0, int *nstop, double *rho);
double seLR_sf(double l, double cl, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0);
double seLR_sf_deluxe(double l, double cl, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0, int *nstop, double *rho);

double seU_q_crit(double l, int L0, double alpha, double hs, double sigma, int df, int N, int qm, double c_error, double a_error);
double se2lu_q_crit(double l, int L0, double alpha, double cl, double hs, double sigma, int df, int N, int qm, double c_error, double a_error);
double se2fu_q_crit(double l, int L0, double alpha, double cu, double hs, double sigma, int df, int N, int qm, double c_error, double a_error);
int se2_q_crit(double l, int L0, double alpha, double *cl, double *cu, double hs, double sigma, int df, int N, int qm, double c_error, double a_error);
int se2_q_crit_class(double l, int L0, double alpha, double *cl, double *cu, double hs, double sigma, int df, double ur, int N, int qm, double c_error, double a_error);
double seUR_q_crit(double l, int L0, double alpha, double cl, double hs, double sigma, int df, int N, int qm, double c_error, double a_error);
double seLR_q_crit(double l, int L0, double alpha, double cu, double hs, double sigma, int df, int N, int qm, double c_error, double a_error);

double seU_Wq(double l, double cu, double p, double hs, double sigma, int df, int N, int nmax, int qm);
double se2_Wq(double l, double cl, double cu, double p, double hs, double sigma, int df, int N, int nmax, int qm);
double seUR_Wq(double l, double cl, double cu, double p, double hs, double sigma, int df, int N, int nmax, int qm);
double seLR_Wq(double l, double cl, double cu, double p, double hs, double sigma, int df, int N, int nmax, int qm);


/* MEWMA: Rigdon (1995a,b) */
double mxewma_arl_0a(double lambda, double ce, int p, double hs, int N); /* GL class */
double mxewma_arl_0a2(double lambda, double ce, int p, double hs, int N); /* GL mod */
double mxewma_arl_0b(double lambda, double ce, int p, double hs, int N, int qm); /* collocation */
double mxewma_arl_0c(double lambda, double ce, int p, double hs, int N); /* Radau (Rigdon) */
double mxewma_arl_0d(double lambda, double ce, int p, double hs, int N); /* Clenshaw-Curtis */
double mxewma_arl_0e(double lambda, double ce, int p, double hs, int N); /* Markov chain (Runger/Prabhu) */
double mxewma_arl_0f(double lambda, double ce, int p, double hs, int N); /* Simpson rule (poor performance) */

double mxewma_arl_f_0a(double lambda, double ce, int p, int N, double *ARL, double *w, double *z);
double mxewma_arl_f_0a2(double lambda, double ce, int p, int N, double *ARL, double *w, double *z);
double mxewma_arl_f_0b(double lambda, double ce, int p, int N, int qm, double *ARL);
double mxewma_arl_f_0c(double lambda, double ce, int p, int N, double *ARL, double *w, double *z);
double mxewma_arl_f_0d(double lambda, double ce, int p, int N, double *ARL, double *w, double *z);
double mxewma_arl_f_0e(double lambda, double ce, int p, int N, double *ARL, double *z);
double mxewma_arl_f_0f(double lambda, double ce, int p, int N, double *ARL, double *w, double *z);

double mxewma_arl_1a (double lambda, double ce, int p, double delta, double hs, int N); /* GL class */
double mxewma_arl_1a2(double lambda, double ce, int p, double delta, double hs, int N); /* GL mod */
double mxewma_arl_1a3(double lambda, double ce, int p, double delta, double hs, int N); /* GL again mod sin, default for 2 and 4 */
double mxewma_arl_1a4(double lambda, double ce, int p, double delta, double hs, int N); /* GL again mod tan */
double mxewma_arl_1a5(double lambda, double ce, int p, double delta, double hs, int N); /* GL again mod sinh, default for all other p */

double mxewma_arl_1q (double lambda, double ce, int p, double delta, int N); /* GL, changed integration order */
double mxewma_arl_1r (double lambda, double ce, int p, double delta, int N);
double mxewma_arl_1s (double lambda, double ce, int p, double delta, int N);
double mxewma_arl_1t (double lambda, double ce, int p, double delta, int N);
double mxewma_arl_1u (double lambda, double ce, int p, double delta, int N);

double mxewma_arl_f_1a (double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1); /* GL class */
double mxewma_arl_f_1a2(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1); /* GL mod */
double mxewma_arl_f_1a3(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1); /* GL again mod sin, default for 2 and 4 */
double mxewma_arl_f_1a4(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1); /* GL again mod tan */
double mxewma_arl_f_1a5(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1); /* GL again mod sinh, default for all other p */

double mxewma_arl_f_1q (double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1); /* GL, changed integration order */
double mxewma_arl_f_1r (double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1);
double mxewma_arl_f_1s (double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1);
double mxewma_arl_f_1t (double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1);
double mxewma_arl_f_1u (double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1);

double mxewma_arl_1b(double lambda, double ce, int p, double delta, double hs, int N, int qm0, int qm1); /* collocation */
double mxewma_arl_1b2(double lambda, double ce, int p, double delta, double hs, int N, int qm0, int qm1); /* collocation, trimmed support of outer integral */
double mxewma_arl_1b3(double lambda, double ce, int p, double delta, double hs, int N, int qm0, int qm1); /* collocation, tan instead of sin */
double mxewma_arl_1b4(double lambda, double ce, int p, double delta, double hs, int N, int qm0, int qm1); /* collocation, sinh instead of sin */

double mxewma_arl_1c(double lambda, double ce, int p, double delta, double hs, int N); /* Radau (Rigdon) */
double mxewma_arl_1d(double lambda, double ce, int p, double delta, double hs, int N); /* Clenshaw-Curtis */
double mxewma_arl_1e(double lambda, double ce, int p, double delta, double hs, int N); /* Markov chain (Runger/Prabhu) */
double mxewma_arl_1f(double lambda, double ce, int p, double delta, double hs, int N); /* Simpson rule (poor performance) */

double mxewma_arl_f_1b (double lambda, double ce, int p, double delta, int N, int qm0, int qm1, double *g); /* collocation with two halfs in the same step + sin() */
double mxewma_arl_f_1b3(double lambda, double ce, int p, double delta, int N, int qm0, int qm1, double *g); /* collocation with two halfs in the same step */
double mxewma_arl_f_1b2(double lambda, double ce, int p, double delta, int N, int qm0, int qm1, double *g); /* collocation with shrinked supports of the outer integral */
double mxewma_arl_f_1b4(double lambda, double ce, int p, double delta, int N, int qm0, int qm1, double *g); /* collocation with two halfs in the same step + sinh() instead of sin() */

double mxewma_arl_f_1c (double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1); /* GL + Radau (Rigdon) */
double mxewma_arl_f_1d (double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1); /* Clenshaw-Curtis */
double mxewma_arl_f_1e (double lambda, double ce, int p, double delta, int N, double *g, int *dQ); /* Markov Chain (Runger/Prabhu) */
double mxewma_arl_f_1f (double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z); /* Simpson rule */

double mxewma_crit(double lambda, double L0, int p, double hs, int N);

double mxewma_psi (double lambda, double ce, int p, int N, double *PSI, double *w, double *z);
double mxewma_psiS(double lambda, double ce, int p, double hs, int N, double *PSI, double *w, double *z);

/* Markov chain (Runger/Prabhu) */
double mxewma_psi0_e(double lambda, double ce, int p, int N, double *PSI);
double mxewma_psi1_e(double lambda, double ce, int p, int N, double *PSI);
double mxewma_psiS0_e(double lambda, double ce, int p, int N, double *PSI);
double mxewma_psiS1_e(double lambda, double ce, int p, int N, double *PSI);

double mxewma_ad(double lambda, double ce, int p, double delta, int N, int qm2, int psi_type, double hs, int qtype, int qm0, int qm1);
double mxewma_ad_new(double lambda, double ce, int p, double delta, int N, int psi_type, double hs, int qtype);

/* Markov chain (Runger/Prabhu) */
double mxewma_ad_e(double lambda, double ce, int p, double delta, int psi_type, int N);


/* incorporate pre-run uncertainty */
double seU_sf_prerun_SIGMA_deluxe(double l, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0);
double seU_sf_prerun_SIGMA(double l, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0);
double seUR_sf_prerun_SIGMA_deluxe(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0);
double seUR_sf_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0);
double se2_sf_prerun_SIGMA_deluxe(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0);
double se2_sf_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0);
double seLR_sf_prerun_SIGMA_deluxe(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0);
double seLR_sf_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0);

double  seU_iglarl_prerun_SIGMA(double l, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate);
double seUR_iglarl_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate);
double  se2_iglarl_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate);
double seLR_iglarl_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate);

double seU_q_crit_prerun_SIGMA(double l, int L0, double alpha, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate, int tail_approx, double c_error, double a_error);
double se2lu_q_crit_prerun_SIGMA(double l, int L0, double alpha, double cl, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate, int tail_approx, double c_error, double a_error);
double se2fu_q_crit_prerun_SIGMA(double l, int L0, double alpha, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate, int tail_approx, double c_error, double a_error);
int se2_q_crit_prerun_SIGMA(double l, int L0, double alpha, double *cl, double *cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate, int tail_approx, double c_error, double a_error);
double seUR_q_crit_prerun_SIGMA(double l, int L0, double alpha, double cl, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate, int tail_approx, double c_error, double a_error);
double seLR_q_crit_prerun_SIGMA(double l, int L0, double alpha, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate, int tail_approx, double c_error, double a_error);

double seU_Wq_prerun_SIGMA_deluxe(double l, double cu, double p, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate);
double seUR_Wq_prerun_SIGMA_deluxe(double l, double cl, double cu, double p, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate);
double seLR_Wq_prerun_SIGMA_deluxe(double l, double cl, double cu, double p, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate);

double se2_Wq_prerun_SIGMA_deluxe(double l, double cl, double cu, double p, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate);

double seU_crit_prerun_SIGMA(double l, double L0, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate);
double se2lu_crit_prerun_SIGMA(double l, double L0, double cl, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate);
double se2fu_crit_prerun_SIGMA(double l, double L0, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate);
int se2_crit_prerun_SIGMA(double l, double L0, double *cl, double *cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate);
double seUR_crit_prerun_SIGMA(double l, double L0, double cl, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate);
double seLR_crit_prerun_SIGMA(double l, double L0, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate);


/* simultaneous EWMA charts */
double xseU_arl(double lx, double ls, double cx, double cs, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm);
double xse2_arl(double lx, double ls, double cx, double csl, double csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm);

int xseU_crit(double lx, double ls, double L0, double *cx, double *cs, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm);
int xse2lu_crit(double lx, double ls, double L0, double *cx, double csl, double *csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm);
int xse2fu_crit(double lx, double ls, double L0, double *cx, double *csl, double csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm);
int xse2_crit(double lx, double ls, double L0, double *cx, double *csl, double *csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm);

double xseU_sf(double lx, double ls, double cx, double cs, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm, double *p0);
double xseU_sf_deluxe(double lx, double ls, double cx, double cs, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm, double *p0, int *nstop, double *rho);
double xse2_sf(double lx, double ls, double cx, double csl, double csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm, double *p0);
double xse2_sf_deluxe(double lx, double ls, double cx, double csl, double csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm, double *p0, int *nstop, double *rho);

int xseU_q_crit(double lx, double ls, int L0, double alpha, double *cx, double *cs, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int qm, double c_error, double a_error);
int xse2fu_q_crit(double lx, double ls, int L0, double alpha, double *cx, double *csl, double csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int qm, double c_error, double a_error);
int xse2_q_crit(double lx, double ls, int L0, double alpha, double *cx, double *csl, double *csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int qm, double c_error, double a_error);

double xseU_Wq(double lx, double ls, double cx, double cs, double p, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm);
double xse2_Wq(double lx, double ls, double cx, double csl, double csu, double p, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm);


/* EWMA p under sampling by variables */

double WK_h(double mu, double sigma, double LSL, double USL);
double wk_h_mu(double mu, double sigma, double LSL, double USL);
double wk_h_sigma(double mu, double sigma, double LSL, double USL);
double WK_h_invers_mu(double p, double sigma, double LSL, double USL);
double WK_h_invers_sigma(double p, double mu, double LSL, double USL);

double wk_alpha(double p, double sigma, int n, double LSL, double USL);

double cdf_phat(double p, double mu, double sigma, int n, double LSL, double USL);
double pdf_phat(double p, double mu, double sigma, int n, double LSL, double USL);
double qf_phat(double p0, double mu, double sigma, int n, double LSL, double USL);

double wk_cdf_i(double y, double p, double mu, double sigma, int n, double LSL, double USL);
double wk_pdf_i(double y, double p, double mu, double sigma, int n, double LSL, double USL);

double cdf_phat2(double p, double mu, double sigma, int n, double LSL, double USL, int nodes);
double pdf_phat2(double p, double mu, double sigma, int n, double LSL, double USL, int nodes);
double qf_phat2(double p0, double mu, double sigma, int n, double LSL, double USL, int nodes);

double ewma_phat_arl (double lambda, double ucl, double mu, double sigma, int n, double z0, double LSL, double USL, int N, int qm);
double ewma_phat_arl_be(double lambda, double ucl, double mu, double sigma, int n, double z0, double LSL, double USL, int N);
double ewma_phat_crit(double lambda, double  L0, double mu, double sigma, int n, double z0, double LSL, double USL, int N, int qm);
double ewma_phat_lambda(double L0, double mu, double sigma, double max_l, double min_l, int n, double z0, double LSL, double USL, int qm);

double ewma_phat_arl2 (double lambda, double ucl, double mu, double sigma, int n, double z0, double LSL, double USL, int N, int qm, int M);
double ewma_phat_arl2_be(double lambda, double ucl, double mu, double sigma, int n, double z0, double LSL, double USL, int N);
double ewma_phat_crit2(double lambda, double  L0, double mu, double sigma, int n, double z0, double LSL, double USL, int N, int qm, int M);
double ewma_phat_lambda2(double L0, double mu, double sigma, double max_l, double min_l, int n, double z0, double LSL, double USL, int qm, int M);


/* attribute EWMA p (X follows binomial distribution) */

double ewma_pU_arl(double lambda, double ucl, int n, double p, double z0, int d_res, int round_mode, int mid_mode);
double ewma_pL_arl(double lambda, double lcl, int n, double p, double z0, int d_res, int round_mode, int mid_mode);
double ewma_p2_arl(double lambda, double lcl, double ucl, int n, double p, double z0, int d_res, int round_mode, int mid_mode);

/* attribute EWMA p (X follows Poisson distribution) */

double cewma_U_arl(double lambda, double AU, double mu0, double z0, double mu, int N);
double cewma_L_arl(double lambda, double AL, double AU, double mu0, double z0, double mu, int N);
double cewma_2_arl(double lambda, double AL, double AU, double mu0, double z0, double mu, int N);
double cewma_2_arl_rando(double lambda, double AL, double AU, double gammaL, double gammaU, double mu0, double z0, double mu, int N);
double cewma_U_crit(double lambda, double L0, double mu0, double z0, int N, int jmax);
double cewma_L_crit(double lambda, double L0, double AU, double mu0, double z0, int N, int jmax);
double cewma_2_crit_sym(double lambda, double L0, double mu0, double z0, int N, int jmax);
double cewma_2_crit_AL(double lambda, double L0, double AU, double mu0, double z0, int N, int jmax);
double cewma_2_crit_AU(double lambda, double L0, double AL, double mu0, double z0, int N, int jmax);
int cewma_2_crit_unb(double lambda, double L0, double mu0, double z0, int N, int jmax, double *AL, double *AU);

double cewma_U_arl_new(double lambda, double AU, double mu0, double z0, double mu, int N);
double cewma_L_arl_new(double lambda, double AL, double AU, double mu0, double z0, double mu, int N);
double cewma_2_arl_new(double lambda, double AL, double AU, double mu0, double z0, double mu, int N);
double cewma_2_arl_rando_new(double lambda, double AL, double AU, double gammaL, double gammaU, double mu0, double z0, double mu, int N);
double cewma_U_crit_new(double lambda, double L0, double mu0, double z0, int N, int jmax);
double cewma_L_crit_new(double lambda, double L0, double AU, double mu0, double z0, int N, int jmax);
double cewma_2_crit_sym_new(double lambda, double L0, double mu0, double z0, int N, int jmax);
double cewma_2_crit_AL_new(double lambda, double L0, double AU, double mu0, double z0, int N, int jmax);
double cewma_2_crit_AU_new(double lambda, double L0, double AL, double mu0, double z0, int N, int jmax);
int cewma_2_crit_unb_new(double lambda, double L0, double mu0, double z0, int N, int jmax, double *AL, double *AU);
int cewma_2_crit_unb_rando_new(double lambda, double L0, double mu0, double z0, int N, int jmax, double *AL, double *AU, double *gL, double *gU);


/* TEWMA (thinning operation -- X follows Poisson distribution */

double tewma_arl(double lambda, int k, int lk, int uk, double z0, double mu);
double tewma_arl_R(double lambda, int k, int lk, int uk, double gl, double gu, double z0, double mu);

/*  Rakitzis / Castagliola / Maravelakis (2015), A new memory-type monitoring technique for count data, doi 10.1016/j.cie.2015.03.021 */

double eewma_arl(int gX, int gY, int kL, int kU, double mu, double y0, int r0);


/* HIER */
/* attribute CUSUM (X follows Poisson distribution) */

double ccusum_U_arl(double mu, int km, int hm, int m, int i0);
double ccusum_U_arl_rando(double mu, int km, int hm, int m, double gamma, int i0);
int ccusum_U_crit(double A, double mu0, int km, int m, int i0);
int ccusum_U_rando_crit(double A, double mu0, int km, int m, int i0, int *hm, double *gamma);

double ccusum_L_arl(double mu, int km, int hm, int m, int i0);
double ccusum_L_arl_rando(double mu, int km, int hm, int m, double gamma, int i0);
int ccusum_L_crit(double A, double mu0, int km, int m, int i0);
int ccusum_L_rando_crit(double A, double mu0, int km, int m, int i0, int *hm, double *gamma); 

double ccusum_2_arl(double mu, int km1, int hm1, int m1, int i01, int km2, int hm2, int m2, int i02);
double ccusum_2_arl_rando(double mu, int km1, int hm1, int m1, double gamma1, int i01, int km2, int hm2, int m2, double gamma2, int i02);


/* tolerance intervals */

double kww(int n, double q, double a);
double tl_factor(int n, double q, double a, int m);


/* internal functions etc. */

int qm_for_l_and_c(double l, double c);
int choose_N_for_seU(double lambda);
int choose_N_for_se2(double lambda, double cl, double cu);

void gausslegendre(int n, double x1, double x2, double *x, double *w);
void radau(int n, double x1, double x2, double *x, double *w);

int LU_decompose(double *a, int *ps, int n);
void LU_solve(double *a, double *b, int n);
void LU_solve2(double *a, double *b, int *ps, int n);

void pmethod(int n, double *p, int *status, double *lambda, double x_[], int *noofit);

int *ivector(long n);
double *vector (long n);
double *matrix(long m, long n);

double phi(double x, double mu);
double PHI(double x, double mu);
double qPHI(double p);

double chi(double s, int df);
double CHI(double s, int df);
double qCHI(double p, int df);
double nchi(double s, int df, double ncp);
double nCHI(double s, int df, double ncp);
double nqCHI(double p, int df, double ncp);

double pdf_t(double x, int df);
double cdf_t(double x, int df);
double  qf_t(double x, int df);

double pdf_tn(double x, int df, double ncp);
double cdf_tn(double x, int df, double ncp);
double  qf_tn(double x, int df, double ncp);

double cdf_binom(double q, int n, double p);
double qf_binom(double q, int n, double p);
double pdf_binom(double x, int n, double p);

double cdf_pois(double q, double lambda);
double qf_pois(double q, double lambda);
double pdf_pois(double x, double lambda);

double Tn(double z, int n);    /* Chebyshev polynomials */
double iTn(double z, int n);   /* indefinite integrals of Chebyshev polynomials */
double dTn(double z, int n);   /* derivatives of Chebyshev polynomials */

double rho0;


/* ------------------- functions and procedures ------------- */

int *ivector(long n)
{
  return (int *) Calloc( n, int ); 
}

double *vector(long n)
{
  return (double *) Calloc( n, double );
}

double *matrix(long m, long n)
{
  return (double *) Calloc( m*n, double );
}

/* normal density (pdf) */

double phi(double x, double mu)
{
 return dnorm(x,mu,1.,LOG);
}

/* normal cumulative distribution function (cdf) */

double PHI(double x, double mu)
{
 return pnorm(x,mu,1.,TAIL,LOG);
}

/* qf of normal rv */

double qPHI(double p)
{
 return qnorm(p,0.,1.,TAIL,LOG);
}

/* pdf of chisquare rv */

double chi(double s, int df)
{
 return dchisq(s,(double)df,LOG);
}


/* pdf of non-central chisquare rv */

double nchi(double s, int df, double ncp)
{
 return dnchisq(s,(double)df,ncp,LOG);
}

/* cdf of chisquare rv */

double CHI(double s, int df)
{
 return pchisq(s,(double)df,TAIL,LOG);
}

/* cdf of non-central chisquare rv */

double nCHI(double s, int df, double ncp)
{
 return pnchisq(s,(double)df,ncp,TAIL,LOG);
}

/* qf of chisquare rv */

double qCHI(double p, int df)
{
 return qchisq(p,(double)df,TAIL,LOG);
}

/* qf of non-central chisquare rv */

double nqCHI(double p, int df, double ncp)
{
 return qnchisq(p,(double)df,ncp,TAIL,LOG);
}


/* pdf of t distribution */

double pdf_t(double x, int df)
{
 return dt(x,(double)df,LOG);
}

/* cdf of t distribution */

double cdf_t(double x, int df)
{
 return pt(x,(double)df,TAIL,LOG);
}

/* quantile function of t distribution */

double qf_t(double x, int df)
{
 return qt(x,(double)df,TAIL,LOG);
}


/* pdf of non-central t distribution */

double pdf_tn(double x, int df, double ncp)
{
 return dnt(x,(double)df,ncp,LOG);
}

/* cdf of non-central t distribution */

double cdf_tn(double x, int df, double ncp)
{
 return pnt(x,(double)df,ncp,TAIL,LOG);
}

/* quantile function of non-central t distribution */

double qf_tn(double x, int df, double ncp)
{
 return qnt(x,(double)df,ncp,TAIL,LOG);
}


/* cdf of binomial rv */
double cdf_binom(double q, int n, double p)
{
  return pbinom(q,(double)n,p,TAIL,LOG);
}  

/* qf of binomial rv */
double qf_binom(double q, int n, double p)
{
  return qbinom(q,(double)n,p,TAIL,LOG);
} 

/* pdf of binomial rv */
double pdf_binom(double x, int n, double p)
{
  return dbinom(x,(double)n,p,LOG);
}  


/* cdf of Poisson rv */
double cdf_pois(double q, double lambda)
{
  return ppois(q,lambda,TAIL,LOG);
}

/* qf of Poisson rv */
double qf_pois(double q, double lambda)
{
  return qpois(q,lambda,TAIL,LOG);
} 

/* pdf of Poisson rv */
double pdf_pois(double x, double lambda)
{
  return dpois(x,lambda,LOG);
}


/* expectation of log-gamma */
double E_log_gamma(double ddf)
{
  return log(2./ddf) + digamma(ddf/2.); 
}

/* variance of log-gamma */
double V_log_gamma(double ddf)
{
  return trigamma(ddf/2.); 
}

/* expectation of S (chi) */
double c_four(double ddf)
{
  return sqrt( 2./ddf ) * gammafn( (ddf+1)/2. ) / gammafn( ddf/2. );
}


/* lapack routine */
void solve(int *n, double *a, double *b)
{ int nrhs=1, lda, ldb, *ipiv, info=0;
  lda = *n;
  ldb = *n; 
  ipiv =  ivector(*n);
  F77_NAME(dgesv)(n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  Free(ipiv);
}    


/* abscissae and weights of Gauss-Legendre quadrature */

#define GLeps 3e-11

void gausslegendre(int n, double x1, double x2, double *x, double *w)
/*
   The following algorithm is based on ideas of Knut Petras
   (see http://www-public.tu-bs.de:8080/~petras/).

   The nodes are derived by means of the Newton method.
   Afterwards, the weights are obtained by utilizing
   (regarding the connection between the Christoffel function
    and the weight, which is also called Christoffel number)

   w_i = w(x_i) = 2 / sum_j=0^n ( (2j+1) * (P_j(x_i))^2 )

   which is more stable than to rely on the usual

   w_i = 2 / (1-x_i^2)/(P_n^'(x_i))^2.

   Note that the Newton method is stopped as soon as the distance
   between two successive iterates is smaller than GLeps, plus
   one extra step.

   By comparing with results in Yakimiw (1996)
   we may conclude that the code behaves very well and even better.
*/
{ double xw, xmid, z0, z1, diff, p0, p1, p2=0., a;
  int i, j, m, stop, odd;

 m = (n+1)/2;
 odd = n%2 == 1;
 xmid = .5*(x2+x1);   /* interval centre */
 xw = .5*(x2-x1);     /* half interval length */

 for (i=0;i<m;i++) {
  if (odd && i==m-1)
    z1 = 0.;
  else {
    z0 = -cos( PI*(i+.75)/(n+.5) );  /* initial guess */
    stop = 0;
    diff = 1;
    while (stop<2) {
      p0 = 1.;
      p1 = z0;
      for (j=1;j<n;j++) { /* iterate to get the nth Legendre polynomial at z0 */
        p2 = ( (2.*j+1.)*z0*p1 - j*p0 )/(j+1.);
        p0 = p1;
        p1 = p2;
      }
      z1 = z0 + (1.-z0*z0)*p2/n/(z0*p2-p0); /* Newton method update, where   */
      diff = fabs(z1-z0);                   /* derivative is based on P_n(x) */
      z0 = z1;                              /* and P_n-1(x)                  */
      if (diff<GLeps) stop++; /* stop as soon diff is small enough      */
    }                         /* (for 2 times -> kind of overiterating) */
  }

  x[i]     = xmid + xw*z1;
  x[n-1-i] = xmid - xw*z1; /* nodes on interval (x1,x2) */

  p0 = 1.;
  p1 = z1;
  a = 1. + 3.*z1*z1;
  for (j=1;j<n;j++) {
    p2 = ( (2.*j+1.)*z1*p1 - j*p0 )/(j+1.);
    p0 = p1;
    p1 = p2;
    a += p1*p1*(2.*j+3.);
  } /* Christoffel function based approach which is more stable */

  w[i]     = 2./a * xw;
  w[n-1-i] = w[i];       /* weights for interval (x1,x2) */
 }
}

#undef GLeps


/* helper functions */


double r8_epsilon ( )
/*
  Purpose:  R8_EPSILON returns the R8 roundoff unit.
  
  Discussion:
  
    The roundoff unit is a number R which is a power of 2 with the property 
    that, to the precision of the computer's arithmetic,
      1 < 1 + R
    but
      1 = ( 1 + R / 2 )
  
  Licensing:  This code is distributed under the GNU LGPL license.
  
  Modified:  01 July 2004
  
  Author:  John Burkardt
  
  Parameters: Output, double R8_EPSILON, the double precision round-off unit.
*/
{ double r;

  r = 1.0;
  while ( 1.0 < ( double ) ( 1.0 + r )  )
  {
    r = r / 2.0;
  }
  return ( 2.0 * r );
}


double r8_max ( double x, double y )
/*
  Purpose:  R8_MAX returns the maximum of two R8's.
  
  Licensing:  This code is distributed under the GNU LGPL license.
  
  Modified:  18 August 2004
  
  Author:  John Burkardt
  
  Parameters:  Input, double X, Y, the quantities to compare.
               Output, double R8_MAX, the maximum of X and Y.
*/
{ double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}


double r8_abs ( double x )
/*
  Purpose:  R8_ABS returns the absolute value of an R8.
  
  Licensing:  This code is distributed under the GNU LGPL license.
  
  Modified:  14 November 2006
  
  Author:  John Burkardt
  
  Parameters:  Input, double X, the quantity whose absolute value is desired.
               Output, double R8_ABS, the absolute value of X.
*/
{ double value;

  if ( 0.0 <= x )
  {
    value = x;
  }
  else
  {
    value = -x;
  }
  return value;
}


void radau(int n, double x1, double x2, double *x, double *w)

/******************************************************************************/
/*
  Purpose:
  
    RADAU_COMPUTE computes a Radau quadrature rule.
  
  Discussion:
  
    The Radau rule is distinguished by the fact that the left endpoint
    (-1) is always an abscissa.
    
    The integral:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
  
    The quadrature rule will integrate exactly all polynomials up to
    X**(2*NORDER-2).
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    28 August 2007
  
  Author:
  
    Original MATLAB version by Greg von Winckel.
    C version by John Burkardt.
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
    Claudio Canuto, Yousuff Hussaini, Alfio Quarteroni, Thomas Zang,
    Spectral Methods in Fluid Dynamics,
    Springer, 1993,
    ISNB13: 978-3540522058,
    LC: QA377.S676.
  
    Francis Hildebrand,
    Section 8.11,
    Introduction to Numerical Analysis,
    Dover, 1987,
    ISBN13: 978-0486653631,
    LC: QA300.H5.
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3,
    LC: QA47.M315.
  
  Parameters:
  
    Input, int N, the order.
    N must be at least 1.
  
    Output, double X[N], the abscissas.

    Output, double W[N], the weights.
*/
{
  int i;
  int iterate;
  int iterate_max = 25;
  int j;
  double pi = 3.141592653589793;
  double temp;
  double test;
  double tolerance;
  double xw, xmid;

  xmid = .5*(x2+x1);   /* interval centre */
  xw   = .5*(x2-x1);   /* half interval length */ 

  tolerance = 100.0 * r8_epsilon ( );
/*
  Initial estimate for the abscissas is the Chebyshev-Gauss-Radau nodes.
*/
  for ( i = 0; i < n; i++ )
  {
    x[i] = - cos ( 2.0 * pi * ( double ) (         i )
                            / ( double ) ( 2 * n - 1 ) );
  }
  double xold[n];
  double p[n*(n+1)];
  iterate = 0;

  do
  {
    for ( i = 0; i < n; i++ )
    {
      xold[i] = x[i];
    }

    temp = 1.0;
    for ( j = 0; j < n + 1; j++ )
    {
      p[0+j*n] = temp;
      temp = -temp;
    }

    for ( i = 1; i < n; i++ )
    {
      p[i+0*n] = 1.0;
    }
    for ( i = 1; i < n; i++ )
    {
      p[i+1*n] = x[i];
    }

    for ( j = 2; j <= n; j++ )
    {
      for ( i = 1; i < n; i++ )
      {
        p[i+j*n] = ( ( double ) ( 2 * j - 1 ) * x[i] * p[i+(j-1)*n]
                   + ( double ) (   - j + 1 ) *        p[i+(j-2)*n] )
                   / ( double ) (     j     );
      }
    }
    for ( i = 1; i < n; i++ )
    {
      x[i] = xold[i] - ( ( 1.0 - xold[i] ) / ( double ) ( n ) )
        * ( p[i+(n-1)*n] + p[i+n*n] ) / ( p[i+(n-1)*n] - p[i+n*n] );
    }
    test = 0.0;
    for ( i = 0; i < n; i++ )
    {
      test = r8_max ( test, r8_abs ( x[i] - xold[i] ) );
    }
    iterate = iterate + 1;
  } while ( tolerance < test && iterate < iterate_max );

  w[0] = xw * 2.0 / ( double ) ( n * n );
  x[0] = x1;
  for ( i = 1; i < n; i++ )  {
    w[i] = xw * ( 1.0 - x[i] ) / pow ( ( double ) ( n ) * p[i+(n-1)*n], 2 );
    x[i] = xw*x[i] + xmid;
  } 

  return;
}
/******************************************************************************/


void matvec(int n, double *p, double *z, double y_[])
{ int i, j;
 for (i=0;i<n;i++) {
  y_[i] = 0.;
  for (j=0;j<n;j++)
   y_[i] += p[i*n+j] * z[j];
 }
}


/* power method */
#define convgd          0
#define limit           1
#define epsilon         1e-12
#define maxits          100000

void pmethod(int n, double *p, int *status, double *lambda, 
             double x_[], int *noofit)
{ int count, i, newi, oldi;
  double newmu, oldmu, *z, *y_;
  void matvec();

 z  = vector(n);
 y_ = vector(n);

 for (i=1;i<n;i++) z[i] = 0.; z[0] = 1.;

 newmu = 0.; newi = 0;
 count = 0; *status = limit;

 while ( (count<maxits) && (*status==limit) ) {
  count++;
  matvec(n, p, z, y_);
  oldmu = newmu; oldi = newi; newmu = 0.;

  for (i=0;i<n;i++)
   if ( fabs(y_[i])>fabs(newmu) ) { newmu = y_[i]; newi = i; }

  for (i=0;i<n;i++) z[i] = y_[i] / newmu;

  if ( fabs(newmu-oldmu)<=epsilon && newi==oldi ) *status = convgd;
 }

 for (i=0;i<n;i++) x_[i] = z[i];

 if (*status == convgd) { *lambda = newmu;  *noofit = count; }
 else { *noofit = maxits; }
}


/* Brownian Motion ARL approximations for CUSUM */
double BM_xc_arl(double k, double h, double mu)
{ double Delta, b, arl, offset=1.166;  
/*  offset examples
    0     -- Bagshaw/Johnson (1975)
    1.2   -- Reynolds (1975)
    1.166 -- Siegmund (1985)
*/
  Delta = mu - k;
  b = h + offset;
  if ( fabs(Delta) > 1e-10 ) arl = ( exp(-2.*Delta*b) + 2.*Delta*b - 1. )/2./Delta/Delta;
  else arl = b*b;
  return arl;
}


double BM_xc_crit(double k, double L0, double m0)
{ double c1, c2, c3, L1=0., L2=0., L3=0., dc;

  c2 = 0.;
  do {
    c2 += .1;
    L2 = BM_xc_arl(k, c2, m0);
  } while ( L2<L0 );

  c1 = c2 - .1;
  L1 = BM_xc_arl(k, c1, m0);

  do {
    if ( fabs(L2-L1) > 1e-10 ) {
      c3 = c1 + (L0-L1)/(L2-L1) * (c2-c1);
      L3 = BM_xc_arl(k, c3, m0);
      dc = c3-c2; c1 = c2; L1 = L2; c2 = c3; L2 = L3;
    } else {
      dc = 1e-12;
      c3 = c2;
    }
  } while ( (fabs(L0-L3)>1e-6) && (fabs(dc)>1e-9) );
  return c3;
}   



/* ************************************************************************* */
/*      zero-state and steady-state ARl and critical value routines          */

double xc_crit(int ctyp, double k, double L0, double hs, double m0, int N)
{ double c1, c2, c3, L1=0., L2=0., L3=0., dc, k_bm;

 if ( ctyp==cusumC || fabs(hs)>1e-9 ) {
   c2 = 0.;
   do {
     c2 += .5;
     if (ctyp==cusum1) L2 = xc1_iglarl ( k,c2,hs,m0,N );
     if (ctyp==cusum2) L2 = xc2_iglarl ( k,c2,hs,m0,N );
     if (ctyp==cusumC) L2 = xcC_iglarl ( k,c2,hs,m0,N );
   } while (L2<L0);

   c1 = c2 - .5;
   
   if (ctyp==cusum1) L1 = xc1_iglarl ( k,c1,hs,m0,N );
   if (ctyp==cusum2) L1 = xc2_iglarl ( k,c1,hs,m0,N );
   if (ctyp==cusumC) L1 = xcC_iglarl ( k,c1,hs,m0,N );
 } else {
   k_bm = k;
   /*if ( fabs(m0 - k) < 1e-3 ) k_bm = 1e-3;*/
   if ( ctyp==cusum1 ) {
     c2 = BM_xc_crit(k_bm, L0, m0);
   } else {
     c2 = BM_xc_crit(k_bm, 2.*L0, m0);
   }
   c1 = c2 - .2;
   if ( ctyp==cusum1 ) {
     L1 = xc1_iglarl ( k,c1,hs,m0,N );
     L2 = xc1_iglarl ( k,c2,hs,m0,N );
   } else {
     L1 = xc2_iglarl ( k,c1,hs,m0,N );
     L2 = xc2_iglarl ( k,c2,hs,m0,N );
   }
 }
   
 do {
   c3 = c1 + (L0-L1)/(L2-L1) * (c2-c1);
   if (ctyp==cusum1) L3 = xc1_iglarl ( k,c3,hs,m0,N );
   if (ctyp==cusum2) L3 = xc2_iglarl ( k,c3,hs,m0,N );
   if (ctyp==cusumC) L3 = xcC_iglarl ( k,c3,hs,m0,N );
   dc = c3-c2; c1 = c2; L1 = L2; c2 = c3; L2 = L3;
 } while ( (fabs(L0-L3)>1e-6) && (fabs(dc)>1e-9) );
 return c3;
}


double xsr1_crit(double k, double L0, double zr, double hs, double m0, int N, int MPT)
{ double c1, c2, c3, L1, L2, L3, dc;

 c2 = 0.;
 do {
   c2 += .5;
   L2 = xsr1_iglarl(k, c2, zr, hs, m0, N, MPT);
 } while ( L2<L0 );

 c1 = c2 - .5;
 L1 = xsr1_iglarl(k, c1, zr, hs, m0, N, MPT);

 do {
   c3 = c1 + (L0-L1)/(L2-L1) * (c2-c1);
   L3 = xsr1_iglarl(k, c3, zr, hs, m0, N, MPT);
   dc = c3-c2; c1 = c2; L1 = L2; c2 = c3; L2 = L3;
 } while ( (fabs(L0-L3)>1e-6) && (fabs(dc)>1e-9) );
 return c3;
}


double xe_crit(int ctyp, double l, double L0, double zr, double hs, double m0, int ltyp, int N, double c0)
{ double c1, c2, c3, L1=0., L2=0., L3=0., dc, norm, L2old=0., c2old=0.;
  int nmax=100000;

 if ( (ctyp==ewma1 && c0 < zr) || (ctyp==ewma2 && c0 < 0.) ) c2 = 1.; else c2 = c0;

 do {
   if ( ctyp==ewma1 ) {
     if ( ltyp==fix && hs>=0. ) L2 = xe1_iglarl ( l,c2,zr,hs,m0,N );
     if ( ltyp==fix && hs<0.  ) L2 = xe1_iglarl ( l,c2,zr,c2/2,m0,N );
     if ( ltyp>fix )            L2 = xe1_arlm ( l,c2,zr,hs,1,m0,m0,ltyp,N,nmax );
   }
   if ( ctyp==ewma2 ) {
     if ( ltyp==fix ) L2 = xe2_iglarl ( l,c2,hs,m0,N );
     if ( ltyp>fix ) {
       if ( hs<0. && ltyp==fir  ) L2 = xe2_arlm ( l,c2,c2/2.,1,m0,m0,ltyp,N,nmax );
       if ( hs<0. && ltyp==both ) L2 = xe2_arlm ( l,c2,c2/2.*sqrt(l*(2.-l)),1,m0,m0,ltyp,N,nmax );
       if ( hs>=0. )              L2 = xe2_arlm ( l,c2,hs,1,m0,m0,ltyp,N,nmax );
     }
   }  
   if ( L2 < 1. ) c2 -= .1;
 } while ( L2 < 1. && c2 > .00001 );
 
 if ( L2 < 1. ) error("invalid ARL value");
 if ( L2 > L0 ) { norm = -.1; } else { norm = .5; }
 if ( L2 < 1. + 1e-12 ) { c2 = 0.; norm = .1; }
 if ( (ctyp==ewma1 && c0 > zr) || (ctyp==ewma2 && c0 > 0.) ) norm /= 10.;

 do {
   L2old = L2;
   c2old = c2;
   c2 += norm;
   do {
     if ( ctyp==ewma1 ) {
       if ( ltyp==fix && hs>=0. ) L2 = xe1_iglarl ( l,c2,zr,hs,m0,N );
       if ( ltyp==fix && hs<0.  ) L2 = xe1_iglarl ( l,c2,zr,c2/2,m0,N );
       if ( ltyp>fix )            L2 = xe1_arlm ( l,c2,zr,hs,1,m0,m0,ltyp,N,nmax );
     }
     if ( ctyp==ewma2 ) {
       if ( ltyp==fix ) L2 = xe2_iglarl ( l,c2,hs,m0,N );
       if ( ltyp>fix  ) {
         if ( hs<0. && ltyp==fir  ) L2 = xe2_arlm ( l,c2,c2/2.,1,m0,m0,ltyp,N,nmax );
         if ( hs<0. && ltyp==both ) L2 = xe2_arlm ( l,c2,c2/2.*sqrt(l*(2.-l)),1,m0,m0,ltyp,N,nmax );
         if ( hs>=0. )              L2 = xe2_arlm ( l,c2,hs,1,m0,m0,ltyp,N,nmax );
       }
     }
     if ( L2 < 1. ) { norm /= 2.; c2 -= norm; }
     if ( c2 <= 1e-9 && fabs(L2-L2old)>100. ) norm = -.001;     
   } while ( L2 < 1. );
 } while ( ((L2 < L0 && norm>0.) || (L2 > L0 && norm<0.)) && (fabs(norm)>1e-8) );

 c1 = c2old;
 L1 = L2old;

 do {
   c3 = c1 + (L0-L1)/(L2-L1) * (c2-c1);
   norm = .5;
   do {
     if ( ctyp==ewma1 ){
       if ( ltyp==fix && hs>=0. ) L3 = xe1_iglarl ( l,c3,zr,hs,m0,N );
       if ( ltyp==fix && hs<0.  ) L3 = xe1_iglarl ( l,c3,zr,c3/2,m0,N );
       if ( ltyp>fix )            L3 = xe1_arlm ( l,c3,zr,hs,1,m0,m0,ltyp,N,nmax );
     }
     if ( ctyp==ewma2 ) {
       if ( ltyp==fix ) L3 = xe2_iglarl ( l,c3,hs,m0,N );
       if ( ltyp>fix ) {
         if ( hs<0. && ltyp==fir  ) L3 = xe2_arlm ( l,c3,c3/2.,1,m0,m0,ltyp,N,nmax );
         if ( hs<0. && ltyp==both ) L3 = xe2_arlm ( l,c3,c3/2.*sqrt(l*(2.-l)),1,m0,m0,ltyp,N,nmax );
         if ( hs>=0. )              L3 = xe2_arlm ( l,c3,hs,1,m0,m0,ltyp,N,nmax );
       }
     }
     if ( L3 < 1. ) {
       c3 = c1 + norm*(L0-L1)/(L2-L1) * (c2-c1);
       norm /= 2.;
     }
   } while ( (L3 < 1.) &&  (fabs(norm)>1e-8) );
   dc = c3-c2; c1 = c2; L1 = L2; c2 = c3; L2 = L3;
   if ( L3 < 1. ) error("invalid ARL value");
 } while ( (fabs(L0-L3)>1e-6) && (fabs(dc)>1e-9) );
 if ( fabs(L0-L3)>1e-6 ) warning("did not converge");
 return c3;
}


double xe_q_crit(int ctyp, double l, int L0, double alpha, double zr, double hs, double m0, int ltyp, int N, double c_error, double a_error)
{ double c1=0., c2=0., c3=0., p1=1., p2=1., p3=1., dc, *SF;
  int result=1;

 SF  = vector(L0);

 c2 = 0.; p2 = 1.;
 do {
   p1 = p2;
   c2 += .5;
   if ( ctyp==ewma1 && ltyp==fix ) result = xe1_sf(l, c2, zr, hs, m0, N, L0, SF);
   if ( ctyp==ewma1 && ltyp>fix  ) error("not implemented yet for one-sided EWMA and varying limits");
   if ( ctyp==ewma2 && ltyp==fix ) result = xe2_sf(l, c2, hs, m0, N, L0, SF);
   if ( ctyp==ewma2 && ltyp>fix  ) result = xe2_sfm(l, c2, hs, 1, m0, m0, ltyp, N, L0, SF);
   if ( result != 0 ) warning("trouble in xe_q_crit [package spc]");
   p2 = 1. - SF[L0-1];
 } while ( p2 > alpha );
 c1 = c2 - .5;
 
 do {
   c3 = c1 + ( alpha - p1 )/( p2 - p1 ) * ( c2 - c1 );
   if ( ctyp==ewma1 && ltyp==fix ) result = xe1_sf(l, c3, zr, hs, m0, N, L0, SF);
   if ( ctyp==ewma1 && ltyp>fix  ) error("not implemented yet for one-sided EWMA and varying limits");
   if ( ctyp==ewma2 && ltyp==fix ) result = xe2_sf(l, c3, hs, m0, N, L0, SF);
   if ( ctyp==ewma2 && ltyp>fix  ) result = xe2_sfm(l, c3, hs, 1, m0, m0, ltyp, N, L0, SF);
   if ( result != 0 ) warning("trouble in xe_q_crit [package spc]");
   p3 = 1. - SF[L0-1];

   dc = c3 - c2; c1 = c2; p1 = p2; c2 = c3; p2 = p3;
 } while ( fabs(alpha - p3)>a_error && fabs(dc)>c_error );

 Free(SF);
 
 return c3;
}


double xc1_iglarl(double k, double h, double hs, double mu, int N)
{ double *a, *g, *w, *z, arl;
  int i, j, NN;

 NN = N + 1;
 a = matrix(NN,NN);
 g = vector(NN);
 w = vector(N);
 z = vector(N);

 gausslegendre(N,0.,h,z,w);

 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) a[i*NN+j] = -w[j]*phi(z[j]+k-z[i],mu);
   ++a[i*NN+i];
   a[i*NN+N] = -PHI(k-z[i],mu);
 }
 for (j=0;j<N;j++)
    a[N*NN+j] = -w[j]*phi(z[j]+k,mu);
 a[N*NN+N] = 1. - PHI(k,mu);

 for (j=0;j<NN;j++) g[j] = 1.;
 LU_solve(a,g,NN);

 arl = 1. + PHI(k-hs,mu)*g[N];
 for (j=0;j<N;j++)
   arl += w[j]*phi(z[j]+k-hs,mu) * g[j];

 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}


double xc1_be_arl(double k, double h, double hs, double mu, int N)
{ double *a, *g, arl, z1, z2, w;
  int i, j;

 a = matrix(N,N);
 g = vector(N);  
 
  w = 2.*h/(2.*N - 1.);

 for (i=0; i<N; i++) 
   for (j=0; j<N; j++) {
     z1 = (j-i)*w - w/2.  + k;
     if (j==0) z1 = -10000.;
     z2 = (j-i)*w + w/2.  + k;
     a[i*N + j] = -PHI(z2, mu) + PHI(z1, mu);
     if ( i==j ) a[i*N + i]++;
   }
  
 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);
 
 i = (int) floor(hs/w + .5);
 arl = g[i];

 Free(a);
 Free(g);

 return arl;
}


double xc1_beL_arl(double k, double h, double hs, double mu, int N)
{ double *a, *g, arl, z1, z2, w;
  int i, j;

 a = matrix(N,N);
 g = vector(N);  
 
  w = 2.*h/(2.*N - 1.);

 for (i=0; i<N; i++) 
   for (j=0; j<N; j++) {
     z1 = (j-i)*w - w/2.  + k;
     if (j==0) z1 = -10000.;
     z2 = (j-i)*w + w/2.  + k;
     a[j*N + i] = -PHI(z2, mu) + PHI(z1, mu);
     if ( i==j ) a[i*N + i]++;
   }
  
 for (j=0; j<N; j++) g[j] = 1.;
 solve(&N, a, g);
 
 i = (int) floor(hs/w + .5);
 arl = g[i];

 Free(a);
 Free(g);

 return arl;
}


double xc1_beT_arl(double k, double h, double hs, double mu, int N)
{ double *a, *b1, *b2, *x, *y, *z, *phi, *psi, *g, w, al, ga, et, de, be, arl;
  int i, j, N1;

 N1 = N - 1; 
 
 a   = vector(2*N-1);
 b1  = vector(N);
 b2  = vector(N);
 x   = vector(N);
 y   = vector(N);
 z   = vector(N);
 phi = vector(N);
 psi = vector(N);
 g   = vector(N);

 w = 2.*h / (2.*N-1.);
 
 for (i=0; i<2*N-1; i++) a[i] = - ( PHI( -w*(i-N1) + w/2. + k, mu) - PHI( -w*(i-N1) - w/2. + k, mu) );
 a[N1] += 1.; 

 for (i=0; i<N; i++) {
   b1[i] = 1.;
   b2[i] = PHI( -(double)i*w - w/2. + k, mu);
 }
 
 x[0] = 1./a[N1];
 y[0] = 1./a[N1];
 phi[0] = b1[0]/a[N1];
 psi[0] = b2[0]/a[N1];
 
 for (i=1; i<N; i++) {
   al = 0.;
   for (j=0; j<i; j++) al += a[N1 + i - j] * x[j];
   ga = 0.;
   for (j=0; j<i; j++) ga += a[N1 - 1 - j] * y[j];
   et = -b1[i];
   for (j=0; j<i; j++) et += a[N1 + i - j] * phi[j];
   de = -b2[i];
   for (j=0; j<i; j++) de += a[N1 + i - j] * psi[j];

   be = 1. - al*ga;
   
   z[0] = -ga*x[0] / be;
   for (j=1; j<i; j++) z[j] = ( y[j-1] - ga*x[j] ) / be;
   z[i] = y[i-1] / be;   

   x[0] = x[0] / be;
   for (j=1; j<i; j++) x[j] = ( x[j] - al*y[j-1] ) / be;
   x[i] = -al*y[i-1] / be;
   
   for (j=0; j<=i; j++) y[j] = z[j];
       
   for (j=0; j<i; j++) {
     phi[j] = ( phi[j] - et*z[j] );
     psi[j] = ( psi[j] - de*z[j] );
   }
   phi[i] = -et*z[i];
   psi[i] = -de*z[i];
 }

 be = phi[0] / ( 1. - psi[0] );
 for (i=0; i<N; i++) g[i] = phi[i] + psi[i] * be;

 arl = 1. + PHI(-hs + w/2 + k, mu) * g[0];
 for (i=1; i<N; i++) arl += ( PHI(-hs + (double)i*w + w/2 + k, mu) - PHI(-hs + (double)i*w - w/2 + k, mu) ) * g[i];

 Free(g);
 Free(psi);
 Free(phi);
 Free(z);
 Free(y);
 Free(x);
 Free(b2);
 Free(b1);
 Free(a);

 return arl;
}


double xtc1_iglarl(double k, double h, double hs, int df, double mu, int N, int subst)
{ double *a, *g, *w, *z, arl, norm=1., arg=0., korr=1.;
  int i, j, NN;

 NN = N + 1;
 a = matrix(NN,NN);
 g = vector(NN);
 w = vector(N);
 z = vector(N);
 
 switch ( subst ) {
   case IDENTITY: gausslegendre(N, 0,      h, z, w); norm = 1.; break;
   case SIN:      gausslegendre(N, 0., PI/2., z, w); norm = 1.; break;
   case SINH:     gausslegendre(N, 0.,    1., z, w); norm = sinh(1.); break;
   case TAN:      gausslegendre(N, 0., PI/4., z, w); norm = 1.; break;
 }     
 
 h /= norm; 

 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) {
     switch ( subst ) {
       case IDENTITY: arg =       z[j]   + k -         z[i]; korr = 1.; break;
       case SIN:      arg = h*sin(z[j])  + k -  h*sin(z[i]); korr = h*cos(z[j]); break;
       case SINH:     arg = h*sinh(z[j]) + k - h*sinh(z[i]); korr = h*cosh(z[j]); break;
       case TAN:      arg = h*tan(z[j])  + k -  h*tan(z[i]); korr = h/( cos(z[j])*cos(z[j]) ); break;
     }
     a[i*NN+j] = -w[j] * pdf_t( arg - mu, df) * korr;
   }  
   ++a[i*NN+i];
   switch ( subst ) {
       case IDENTITY: arg = k -         z[i]; break;
       case SIN:      arg = k -  h*sin(z[i]); break;
       case SINH:     arg = k - h*sinh(z[i]); break;
       case TAN:      arg = k -  h*tan(z[i]); break;
   }
   a[i*NN+N] = - cdf_t(arg - mu, df);
 }
 
 for (j=0;j<N;j++) {
   switch ( subst ) {
     case IDENTITY: arg =       z[j]   + k; korr = 1.; break;
     case SIN:      arg = h*sin(z[j])  + k; korr = h*cos(z[j]); break;
     case SINH:     arg = h*sinh(z[j]) + k; korr = h*cosh(z[j]); break;
     case TAN:      arg = h*tan(z[j])  + k; korr = h/( cos(z[j])*cos(z[j]) ); break;
   }
   a[N*NN+j] = -w[j] * pdf_t( arg - mu, df) * korr;
 }  
 a[N*NN+N] = 1. - cdf_t(k - mu, df);

 for (j=0;j<NN;j++) g[j] = 1.;
 LU_solve(a, g, NN);  
 
 switch ( subst ) {
     case IDENTITY: arg = k -         hs; korr = 1.; break;
     case SIN:      arg = k -  h*sin(hs); korr = h*cos(z[j]); break;
     case SINH:     arg = k - h*sinh(hs); korr = h*cosh(z[j]); break;
     case TAN:      arg = k -  h*tan(hs); korr = h/( cos(z[j])*cos(z[j]) ); break;
 }
 arl = 1. + cdf_t(k - hs - mu, df) * g[N]; 
 for (j=0;j<N;j++) {
   switch ( subst ) {
     case IDENTITY: arg =       z[j]   + k -         hs; korr = 1.; break;
     case SIN:      arg = h*sin(z[j])  + k -  h*sin(hs); korr = h*cos(z[j]); break;
     case SINH:     arg = h*sinh(z[j]) + k - h*sinh(hs); korr = h*cosh(z[j]); break;
     case TAN:      arg = h*tan(z[j])  + k -  h*tan(hs); korr = h/( cos(z[j])*cos(z[j]) ); break;
   }
   arl += w[j] * pdf_t( arg - mu, df) * korr * g[j];
 }

 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}


double xc1_Wq(double k, double h, double p, double hs, double mu, int N, int nmax)
{ double *Pn, *w, *z, *p0, *atom, ratio, q_minus=0., q_plus=0., mn_minus=1., mn_plus=0., enumerator=0., Wq=0.;
  int i, j, n;

 w  = vector(N);
 z  = vector(N);
 Pn = matrix(nmax,N);
 p0 = vector(nmax);
 atom = vector(nmax);

 gausslegendre(N,0,h,z,w);

 for (n=1;n<=nmax;n++) {

   if (n==1) {
     for (i=0;i<N;i++)
       Pn[i] = PHI( -z[i]+h+k, mu);
     atom[0] = PHI( h+k, mu);
   } else {
     for (i=0;i<N;i++) {
       Pn[(n-1)*N+i] = PHI( -z[i]+k, mu) * atom[n-2];
       for (j=0;j<N;j++) Pn[(n-1)*N+i] += w[j] * phi( z[j]-z[i]+k, mu) * Pn[(n-2)*N+j];
     }
     atom[n-1] = PHI( k, mu) * atom[n-2];
     for (j=0;j<N;j++) atom[n-1] += w[j] * phi( z[j]+k, mu) * Pn[(n-2)*N+j];
   }

   if (n==1)
     p0[0] = PHI( h-hs+k, mu);
   else {
     p0[n-1] = PHI( -hs+k, mu) * atom[n-2];
     for (j=0;j<N;j++) p0[n-1] += w[j] * phi( z[j]-hs+k, mu) * Pn[(n-2)*N+j];
   }

   if ( p0[n-1] < 1.-p ) {
     Wq = (double)n;
     n = nmax+1;
   } else {     
     mn_minus = 1.; mn_plus = 0.;
     if ( n>1 ) {
       for (i=0; i<N; i++) {
         if (Pn[(n-2)*N+i]==0)
           if (Pn[(n-1)*N+i]==0) ratio = 0.;
           else ratio = 1.;
         else ratio = Pn[(n-1)*N+i]/Pn[(n-2)*N+i];
        if ( ratio<mn_minus ) mn_minus = ratio;
        if ( ratio>mn_plus ) mn_plus = ratio;
       }
       enumerator = log( (1.-p)/p0[n-1] );
       q_minus = (double)n + enumerator/log(mn_minus);
       q_plus  = (double)n + enumerator/log(mn_plus);
       if ( fabs( ceil(q_plus) - ceil(q_minus) ) < .5 ) {
	 Wq = ceil(q_plus);
	 n = nmax +1;
       }
     } /* n > 1 */  
   } /* p0[n-1] >= 1.-p */
 } /* n=1; n<=nmax; n++ */ 

 Free(p0);
 Free(Pn);
 Free(z);
 Free(w);
 Free(atom);

 return Wq;
}


double xc1_sf(double k, double h, double hs, double mu, int N, int nmax, double *p0)
{ double *Pn, *w, *z, *atom;
  int i, j, n;

 w  = vector(N);
 z  = vector(N);
 Pn = matrix(nmax,N);
 atom = vector(nmax);

 gausslegendre(N,0,h,z,w);

 for (n=1;n<=nmax;n++) {
   if (n==1) {
     for (i=0;i<N;i++)
       Pn[i] = PHI( -z[i]+h+k, mu);
     atom[0] = PHI( h+k, mu);
   } else {
     for (i=0;i<N;i++) {
       Pn[(n-1)*N+i] = PHI( -z[i]+k, mu) * atom[n-2];
       for (j=0;j<N;j++) Pn[(n-1)*N+i] += w[j] * phi( z[j]-z[i]+k, mu) * Pn[(n-2)*N+j];
     }
     atom[n-1] = PHI( k, mu) * atom[n-2];
     for (j=0;j<N;j++) atom[n-1] += w[j] * phi( z[j]+k, mu) * Pn[(n-2)*N+j];
   }

   if (n==1)
     p0[0] = PHI( h-hs+k, mu);
   else {
     p0[n-1] = PHI( -hs+k, mu) * atom[n-2];
     for (j=0;j<N;j++) p0[n-1] += w[j] * phi( z[j]-hs+k, mu) * Pn[(n-2)*N+j];
   }
 }

 Free(Pn);
 Free(z);
 Free(w);
 Free(atom);
 
 return 0;
}


double xc1_arlm(double k, double h, double hs, int q, double mu0, double mu1, int N, int nmax)
{ double *p0, *fn, *w, *z, arl0, rho, arl_minus=0., arl, arl_plus=0., mn_minus=0., mn_plus=0., ratio=0.;
  int i, j, n, NN;

 NN = N + 1;
 w   = vector(NN);
 z   = vector(NN);
 fn  = matrix(nmax, NN);
 p0  = vector(nmax);

 gausslegendre(N, 0., h, z, w);

 /* in-control, i. e. n<=q-1 */
 for (n=1;n<=q-1;n++) {
  /* determine f_n, n=1,2,...,q-1 */
  if (n==1) {
    for (i=0;i<N;i++) fn[0*NN+i] = phi( z[i]+k-hs, mu0);
    fn[0*NN+N] = PHI( k-hs, mu0);
  } else {
    for (i=0;i<N;i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * phi(z[i] + k, mu0);
      for (j=0;j<N;j++) {
        fn[(n-1)*NN+i] += w[j] * fn[(n-2)*NN+j] * phi(z[i] + k - z[j], mu0);
      }
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI(k, mu0);
    for (j=0;j<N;j++) fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI(k - z[j], mu0);
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];
 }

 /* out-of-control, i.e. t>=q */
 arl0 = 1.; rho = 0.;

 for (n=q;n<=nmax;n++) {
  /* determine f_n, n=q,q+1,... */

  if (n==1) {
    for (i=0;i<N;i++) fn[0*NN+i] = phi( z[i]+k-hs, mu1);
    fn[0*NN+N] = PHI( k-hs, mu1);
  } else {
    for (i=0;i<N;i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * phi(z[i] + k, mu1);
      for (j=0;j<N;j++) {
        fn[(n-1)*NN+i] += w[j] * fn[(n-2)*NN+j] * phi(z[i] + k - z[j], mu1);
      }
      if (n==q && q>1) fn[(n-1)*NN+i] /= p0[q-2];
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI(k, mu1);
    for (j=0;j<N;j++) fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI(k - z[j], mu1);
    if (n==q && q>1) fn[(n-1)*NN+N] /= p0[q-2];
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];
  
  /* computation of m_n+1^- and m_n+1^+, n=m-1,m,... */
  mn_minus = 1.; mn_plus = 0.;
  if (n>q) {
   for (i=0;i<NN;i++) {
    if (fn[(n-2)*NN+i]==0)
     if (fn[(n-1)*NN+i]==0) ratio = 0.; else ratio = 1.;
    else ratio = fn[(n-1)*NN+i]/fn[(n-2)*NN+i];
    if ( ratio<mn_minus ) mn_minus = ratio;
    if ( ratio>mn_plus ) mn_plus = ratio;
   }
  }

  if (n>q) rho = p0[n-1]/p0[n-2];

  /* computation of ARL, ARL^-, and ARL^+ */
  arl = arl0 + p0[n-1]/(1.-rho);
  if (mn_minus<1.) arl_minus = arl0 + p0[n-1]/(1.-mn_minus);
  else             arl_minus = -1.;
  if (mn_plus<1.) arl_plus = arl0 + p0[n-1]/(1.-mn_plus);
  else            arl_plus = -1.;
  arl0 += p0[n-1];

  if ( fabs((arl_plus-arl_minus)) < 1e-7 ) n = nmax+1;
 }

 arl = (arl_plus+arl_minus)/2.; rho0 = rho;

 Free(w);
 Free(z);
 Free(fn);
 Free(p0);

 return arl;
}


double xc1_arlm_hom(double k, double h, double hs, int q, double mu0, double mu1, int N, double *ced)
{ double *fn, *w, *z, *a, *arl, norm;
  int i, j, n, NN;
  
  NN = N + 1;
  w   = vector(NN);
  z   = vector(NN);
  fn  = matrix(q+1, NN);
  a   = matrix(NN,NN);
  arl = vector(NN);
 
  gausslegendre(N, 0., h, z, w); 
 
  /* ARL vector */
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) a[i*NN+j] = -w[j] * phi( z[j]+k-z[i], mu1); 
    ++a[i*NN+i];
    a[i*NN+N] = - PHI(k-z[i], mu1);
  }
  for (j=0;j<N;j++) a[N*NN+j] = -w[j] * phi( z[j]+k, mu1);
  a[N*NN+N] = 1. - PHI(k, mu1);

  for (j=0;j<NN;j++) arl[j] = 1.;
  LU_solve(a,arl,NN);
 
  /* q == 1 */
  ced[0] = 1. + PHI( k-hs, mu1) * arl[N];
  for (j=0; j<N; j++) ced[0] +=  w[j] * phi( z[j]+k-hs, mu1) * arl[j];


  /* density sequence for q > 1 */
  for (n=1; n<=q-1; n++) {
    if (n==1) {
      for (i=0; i<N; i++) fn[0*NN+i] = phi( z[i]+k-hs, mu0);
      fn[0*NN+N] = PHI( k-hs, mu0);
    } else {
      for (i=0; i<N; i++) {
        fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * phi( z[i] + k, mu0);
        for (j=0; j<N; j++) fn[(n-1)*NN+i] += w[j] * fn[(n-2)*NN+j] * phi( z[i] + k - z[j], mu0);
      }    
      fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI(k, mu0);
      for (j=0;j<N;j++) fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI(k - z[j], mu0);
    }
  
    ced[n] = fn[(n-1)*NN+N] * arl[N];
    norm = fn[(n-1)*NN+N];
    for (j=0; j<N; j++) {
      ced[n] += w[j] * fn[(n-1)*NN+j] * arl[j]; 
      norm += w[j] * fn[(n-1)*NN+j];
    }
    ced[n] /= norm;
  }  
  
  Free(w);
  Free(z);
  Free(fn);
  Free(a);
  Free(arl);

  return 0;
}


double xc1_iglarl_drift(double k, double h, double hs, double delta, int m, int N, int with0)
{ double *a, *g, *w, *z, arl, *MUs, *ARLs;
  int i, j, NN, m_;

 NN = N + 1;
 a = matrix(NN, NN);
 g = vector(NN);
 w = vector(NN);
 z = vector(NN);
 ARLs = vector(NN);
 MUs  = vector(m+1);

 gausslegendre(N, 0., h, z, w);

 if ( with0 ) {
   for (i=0;i<=m;i++) MUs[i] = (double)i * delta;
 } else {
   for (i=0;i<=m;i++) MUs[i] = (double)(i+1.) * delta;
 }

 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) a[i*NN+j] = -w[j] * phi( z[j]+k-z[i], MUs[m]);
   ++a[i*NN+i];
   a[i*NN+N] = -PHI( k-z[i], MUs[m]);
 }
 for (j=0;j<N;j++)
    a[N*NN+j] = -w[j] * phi( z[j]+k, MUs[m]);
 a[N*NN+N] = 1. - PHI(k, MUs[m]);

 for (j=0;j<NN;j++) g[j] = 1.;
 LU_solve(a, g, NN);

 for (m_=0;m_<m;m_++) {
   for (i=0;i<=N;i++) {
     ARLs[i] = 1. + PHI( k-z[i], MUs[m-m_]) * g[N];
     for (j=0;j<N;j++) { 
       ARLs[i] += w[j] * phi( z[j]+k-z[i], MUs[m-m_]) * g[j];
     }
   }
   for (j=0;j<=N;j++) g[j] = ARLs[j];
 }

 arl = 1. + PHI( k-hs, MUs[0]) * ARLs[N];
 for (j=0;j<N;j++) arl += w[j] * phi( z[j]+k-hs, MUs[0]) * ARLs[j];

 Free(a);
 Free(g);
 Free(w);
 Free(z);
 Free(ARLs);
 Free(MUs);

 return arl;
}


double xc1_iglarl_drift_wo_m(double k, double h, double hs, double delta, int *m, int N, int with0)
{ int m_;
  double arl1, arl2, eps=1e-6;
 m_ = 4;
 arl1 = xc1_iglarl_drift(k, h, hs, delta, m_, N, with0);
 arl2 = arl1 + 2.*eps;
 while ( fabs(arl2-arl1)>eps && (double)m_<1e4 ) {
   m_ = (int)round(1.5 * m_);
   arl1 = xc1_iglarl_drift(k, h, hs, delta, m_, N, with0);
   arl2 = xc1_iglarl_drift(k, h, hs, delta, m_+1, N, with0);
 }
 *m = m_;
 return arl1;
}


double xc1_iglarlm_drift(double k, double h, double hs, int q, double delta, int N, int nmax, int with0)
{ double *p0, *fn, *w, *z, arl0, rho, MEAN=0.,
         arl_minus=0., arl, arl_plus=0., mn_minus=0., mn_plus=0., nn, ratio=0.;
  int i, j, n, NN;

 NN = N + 1;
 w   = vector(NN);
 z   = vector(NN);
 fn  = matrix(nmax, NN);
 p0  = vector(nmax);

 gausslegendre(N, 0, h, z, w);

 /* in-control, i. e. n<=q-1 */
 MEAN = 0.;

 for (n=1;n<=q-1;n++) {
  nn = (double) n;

  /* determine f_n, n=1,2,...,q-1 */
  if (n==1) {
    for (i=0;i<N;i++) fn[0*NN+i] = phi( z[i]+k-hs, MEAN);
    fn[0*NN+N] = PHI( k-hs, MEAN);
  } else {
    for (i=0;i<N;i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * phi(z[i] + k, MEAN);
      for (j=0;j<N;j++) {
        fn[(n-1)*NN+i] += w[j] * fn[(n-2)*NN+j] * phi(z[i] + k - z[j], MEAN);
      }
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI(k, MEAN);
    for (j=0;j<N;j++) fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI(k - z[j], MEAN);
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];
 }

 /* out-of-control, i.e. t>=q */
 arl0 = 1.; rho = 0.;

 for (n=q;n<=nmax;n++) {
  nn = (double) n;

  /* determine f_n, n=q,q+1,... */
  if ( with0 ) {
    MEAN = (nn-(double)q) * delta;
  } else {
    MEAN = (nn-(double)q+1.) * delta;
  }

  if (n==1) {
    for (i=0;i<N;i++) fn[0*NN+i] = phi( z[i]+k-hs, MEAN);
    fn[0*NN+N] = PHI( k-hs, MEAN);
  } else {
    for (i=0;i<N;i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * phi(z[i] + k, MEAN);
      for (j=0;j<N;j++) {
        fn[(n-1)*NN+i] += w[j] * fn[(n-2)*NN+j] * phi(z[i] + k - z[j], MEAN);
      }
      if (n==q && q>1) fn[(n-1)*NN+i] /= p0[q-2];
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI(k, MEAN);
    for (j=0;j<N;j++) fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI(k - z[j], MEAN);
    if (n==q && q>1) fn[(n-1)*NN+N] /= p0[q-2];
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];

  /* computation of m_n+1^- and m_n+1^+, n=m-1,m,... */
  mn_minus = 1.; mn_plus = 0.;
  if (n>q) {
   for (i=0;i<NN;i++) {
    if (fn[(n-2)*NN+i]==0)
     if (fn[(n-1)*NN+i]==0) ratio = 0.; else ratio = 1.;
    else ratio = fn[(n-1)*NN+i]/fn[(n-2)*NN+i];
    if ( ratio<mn_minus ) mn_minus = ratio;
    if ( ratio>mn_plus ) mn_plus = ratio;
   }
  }

  if (n>q) rho = p0[n-1]/p0[n-2];

  /* computation of ARL, ARL^-, and ARL^+ */
  arl = arl0 + p0[n-1]/(1.-rho);
  if (mn_minus<1.) arl_minus = arl0 + p0[n-1]/(1.-mn_minus);
  else             arl_minus = -1.;
  if (mn_plus<1.) arl_plus = arl0 + p0[n-1]/(1.-mn_plus);
  else            arl_plus = -1.;
  arl0 += p0[n-1];

  if ( fabs((arl_plus-arl_minus)) < 1e-7 ) n = nmax+1;
 }

 arl = (arl_plus+arl_minus)/2.; rho0 = rho;

 Free(w);
 Free(z);
 Free(fn);
 Free(p0);

 return arl;
}


double xc2_iglarl(double k, double h, double hs, double mu, int N)
{ double arl1, arl2, arl3, arl4, arl;

/* relation between 1- and 2-sided CUSUM schemes due to Lucas/Crosier 1982,
   Technometrics 24, 199-205;
   only for headstart hs smaller than h/2 + k !!
*/

 arl1 = xc1_iglarl(k,h,0.,mu,N);
 arl2 = xc1_iglarl(k,h,hs,mu,N);
 arl3 = xc1_iglarl(k,h,0.,-mu,N);
 arl4 = xc1_iglarl(k,h,hs,-mu,N);
 arl = ( arl2*arl3 + arl1*arl4 - arl1*arl3 ) / ( arl1 + arl3 );
 return arl;
}


double xtc2_iglarl(double k, double h, double hs, int df, double mu, int N, int subst)
{ double arl1, arl2, arl3, arl4, arl;

/* relation between 1- and 2-sided CUSUM schemes due to Lucas/Crosier 1982,
   Technometrics 24, 199-205;
   only for headstart hs smaller than h/2 + k !!
*/

 arl1 = xtc1_iglarl(k, h, 0., df,  mu, N, subst);
 arl2 = xtc1_iglarl(k, h, hs, df,  mu, N, subst);
 arl3 = xtc1_iglarl(k, h, 0., df, -mu, N, subst);
 arl4 = xtc1_iglarl(k, h, hs, df, -mu, N, subst);
 arl = ( arl2*arl3 + arl1*arl4 - arl1*arl3 ) / ( arl1 + arl3 );
 return arl;
}


double xc2_be_arl (double k, double h, double hs1, double hs2, double mu, int N)
{ double *a, *g, arl, z1, z2, z11, z12, z21, z22, w;
  int i1, i2, j1, j2, NN, N3;

/* two-dimensional Markov chain approximation */

 NN = N*N; N3 = NN*N;
 a = matrix(NN,NN);
 g = vector(NN);

 w = 2.*h/(2.*N - 1.);

 for (i1=0;i1<N;i1++) 
   for (j1=0;j1<N;j1++) 
     for (i2=0;i2<N;i2++)
       for (j2=0;j2<N;j2++) {
         z11 = (i2-i1)*w - w/2.  + k; if (i2==0) z11 = -10000.;
         z12 = (i2-i1)*w + w/2.  + k;
         z21 = -2.*k - (j2-j1)*w - w/2.  + k;
         z22 = -2.*k - (j2-j1)*w + w/2.  + k; if (j2==0) z22 = 10000.;
         if ( z11 < z21 ) z1 = z21; else z1 = z11;
         if ( z12 < z22 ) z2 = z12; else z2 = z22;
         if ( z1 > z2 ) a[i1*N3+j1*NN+i2*N+j2] = 0.;
         else a[i1*N3+j1*NN+i2*N+j2] = -PHI(z2, mu) + PHI(z1, mu);
         if ( i1==i2 && j1==j2 ) a[i1*N3+j1*NN+i2*N+j2]++;
       }

 for (j1=0;j1<NN;j1++) g[j1] = 1.;
 LU_solve(a, g, NN);
 
 i1 = (int) ceil(hs1/w - .5);
 i2 = (int) ceil(hs2/w - .5);
 arl = g[i1*N + i2];

 Free(a);
 Free(g);

 return arl;
}


double xc2_iglarl_drift(double k, double h, double hs, double delta, int m, int N, int drift0)
{ double arl1, arl2, arl3, arl4, arl;

/* relation between 1- and 2-sided CUSUM schemes due to Lucas/Crosier 1982,
   Technometrics 24, 199-205;
   only for headstart hs smaller than h/2 + k !!
*/

 arl1 = xc1_iglarl_drift(k, h, 0., delta, m, N, drift0);
 arl2 = xc1_iglarl_drift(k, h, hs, delta, m, N, drift0);
 arl3 = xc1_iglarl_drift(k, h, 0., -delta, m, N, drift0);
 arl4 = xc1_iglarl_drift(k, h, hs, -delta, m, N, drift0);
 arl = ( arl2*arl3 + arl1*arl4 - arl1*arl3 ) / ( arl1 + arl3 );
 return arl;
}


double xc2_iglarl_drift_wo_m(double k, double h, double hs, double delta, int *m, int N, int drift0)
{ int m_;
  double arl1, arl2, eps=1e-6;
 m_ = 4;
 arl1 = xc2_iglarl_drift(k, h, hs, delta, m_, N, drift0);
 arl2 = arl1 + 2.*eps;
 while ( fabs(arl2-arl1)>eps && (double)m_<1e4 ) {
   m_ = (int)round(1.5 * m_);
   arl1 = xc2_iglarl_drift(k, h, hs, delta, m_, N, drift0);
   arl2 = xc2_iglarl_drift(k, h, hs, delta, m_+1, N, drift0);
 }
 *m = m_;
 return arl1;
}


double xcC_iglarl (double k, double h, double hs, double mu, int N)
{ double *a, *g, *w, *z, arl;
  int i, j, NN;

 NN = 2*N + 1;
 a = matrix(NN,NN);
 g = vector(NN);
 w = vector(NN);
 z = vector(NN);

 gausslegendre(N,0.,h,z,w); 

 for (i=0;i<N;i++) { /* upper */
   for (j=0;j<N;j++)    a[i*NN+j] = -w[j]  *phi( z[j]  +k-z[i],mu);
   for (j=N;j<NN-1;j++) a[i*NN+j] = -w[j-N]*phi(-z[j-N]-k-z[i],mu);
   ++a[i*NN+i];
   a[i*NN+NN-1] = - ( PHI(k-z[i],mu) - PHI(-k-z[i],mu) );
 }

 for (i=N;i<NN-1;i++) { /* lower */
   for (j=0;j<N;j++)    a[i*NN+j] = -w[j]  *phi( z[j]  +k+z[i-N],mu);
   for (j=N;j<NN-1;j++) a[i*NN+j] = -w[j-N]*phi(-z[j-N]-k+z[i-N],mu);
   ++a[i*NN+i];
   a[i*NN+NN-1] = - ( PHI(k+z[i-N],mu) - PHI(-k+z[i-N],mu) );
 }

 /* "fat" zero  */
 for (j=0;j<N;j++)
    a[(NN-1)*NN+j] = -w[j]  *phi( z[j]  +k,mu);
 for (j=N;j<NN-1;j++)
    a[(NN-1)*NN+j] = -w[j-N]*phi(-z[j-N]-k,mu);
 a[(NN-1)*NN+NN-1] = 1. - ( PHI(k,mu) - PHI(-k,mu) );

 for (j=0;j<NN;j++) g[j] = 1.;
 LU_solve(a,g,NN);

 arl = 1. + ( PHI(k-hs,mu) - PHI(-k-hs,mu) )*g[NN-1];
 for (j=0;j<N;j++)
   arl += w[j]  *phi( z[j]  +k-hs,mu) * g[j];
 for (j=N;j<NN-1;j++)
   arl += w[j-N]*phi(-z[j-N]-k+hs,mu) * g[j];

 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}


/* variance CUSUM charts */
/* double s2cusumU_arl_igl */
double scU_iglarl_v1(double refk, double h, double hs, double sigma, int df, int N, int qm)
{ double *a, *g, *w, *z, *b, arl, Hij, xl, xu, za, dN, ddf, s2, alpha, *zch;
  int i, j, k, M, Ntilde, NN, ii, jj, ihs, qi, qj;
 
 M = ceil( h/refk );
 ihs = ceil( hs/refk );
 if ( ihs<=0 ) ihs = 1;
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde;
 alpha = ddf/2./s2;
 
 a = matrix(NN,NN);
 g = vector(NN);
 b = vector(M+1);
 w = vector(qm);
 z = vector(qm);
 zch = matrix(M,Ntilde);

 /* interval borders b_i */
 b[0] = 0.;
 for (i=1; i<M; i++) b[i] = (double)(i)*refk;
 b[M] = h;

 /* Chebyshev nodes on [b_1,b_2],[b_2,b_3],...,[b_M,hu] */
 for (i=1; i<=M; i++)
   for (j=1; j<=Ntilde; j++)
     zch[(i-1)*Ntilde + j-1] = b[i-1] + (b[i]-b[i-1])/2.*(1.+cos(PI*(2.*(Ntilde-j+1.)-1.)/2./dN));
   
 for (i=1; i<=M; i++)
   for (j=1; j<=Ntilde ;j++) {
     qi = (i-1)*Ntilde + j-1;
     za = zch[(i-1)*Ntilde + j-1] - refk;

     for (ii=1; ii<=M; ii++) {
       if ( za>b[ii-1] ) xl = za; else xl = b[ii-1];
       xu = b[ii];
       if ( df!=2 && b[ii]>za ) { xl = sqrt(xl-za); xu = sqrt(xu-za); }

       for (jj=1; jj<=Ntilde; jj++) {
	 qj = (ii-1)*Ntilde + jj-1;
         if ( b[ii]<za ) a[qi*NN + qj] = 0.;
         else {
	   gausslegendre(qm, xl, xu, z, w);
           Hij = 0.;
           for (k=0; k<qm; k++)
             if ( df==2 )
               Hij += w[k]*Tn((2.*z[k]-b[ii]-b[ii-1])/(b[ii]-b[ii-1]),jj-1)*exp(-z[k]/s2);
             else
               Hij += w[k]*Tn((2.*(za+z[k]*z[k])-b[ii]-b[ii-1])/(b[ii]-b[ii-1]),jj-1)
                 * 2. * pow(z[k], ddf-1.) * exp(-alpha*z[k]*z[k]);
           if ( df==2 ) Hij *= exp(za/s2)/s2;
           else Hij *= pow(alpha,ddf/2.)/gammafn(ddf/2.);
           a[qi*NN + qj] = -Hij;
         }
       }
     }

     for (jj=1; jj<=Ntilde; jj++)
       a[qi*NN + jj-1] -= CHI(-ddf/s2*za, df) * Tn(-1.,jj-1);

     for (jj=1; jj<=Ntilde; jj++)
       a[qi*NN + (i-1)*Ntilde + jj-1] += Tn((2.*zch[(i-1)*Ntilde + j-1]-b[i]-b[i-1])/(b[i]-b[i-1]),jj-1);
   }

 for (j=0;j<NN;j++) g[j] = 1.;

 LU_solve(a, g, NN);

 arl = 0.;
 for (j=1; j<=Ntilde; j++)
   arl += g[(ihs-1)*Ntilde + j-1] * Tn((2.*hs-b[ihs]-b[ihs-1])/(b[ihs]-b[ihs-1]),j-1);

 Free(zch);
 Free(z);
 Free(w);
 Free(b);
 Free(g);
 Free(a);
 
 return arl;
}


double scs_U_iglarl_v1(double refk, double h, double hs, double cS, double sigma, int df, int N, int qm)
{ double *a, *g, *w, *z, *b, *b1, *b2, arl, Hij, xl, xu, za, zb, dN, ddf, s2, alpha, *zch, eps;
  int i, j, i1, i2, k, M, M1, M2, Ntilde, NN, ii, jj, ihs, qi, qj;
 
 eps = cS - refk;
 M2  = ceil( h/eps );  
 M1  = ceil( h/refk ); 
  
 s2    = sigma*sigma;
 ddf   = (double)df;
 alpha = ddf/2./s2;
 
 M = M1 + M2 - 1;
 Ntilde = ceil( (double)N/(double)M );
 dN = (double)Ntilde;
 NN = M*Ntilde;
 
 /*printf("\n\nM1 = %d,\tM2 = %d,\tM = %d\n\n", M1, M2, M);*/
 
 a = matrix(NN,NN);
 g = vector(NN);
 b1 = vector(M1+1);
 b2 = vector(M2+1);
 b  = vector(M+1);
 w = vector(qm);
 z = vector(qm);
 zch = matrix(M,Ntilde);

 /* interval borders support */
 for (i=1; i<M1; i++) b1[i] = (double)(i)*refk;
 b1[M1] = h;
 
 /* interval borders Shewhart */
 for (i=1; i<M2; i++) b2[i] = h - ( (double)M2 - (double)(i) )*eps;
 b2[M2] = h;
 
 /* merge */
 b[0] = 0.;
 i1 = 1; i2 = 1;
 for (i=1; i<M; i++) {
   if ( b1[i1] < b2[i2] ) {
     b[i] = b1[i1];
     i1++;
   } else {
     b[i] = b2[i2];
     i2++;
   }  
 }    
 b[M] = h;    
 
 ihs = M;
 for ( i=2; i<=M; i++ ) {
   if ( b[i-1] > hs ) {
     ihs = i-1; i = M+1;  
   }
 }
 
 /* Chebyshev nodes on [b_1,b_2],[b_2,b_3],...,[b_M,hu] */
 for (i=1; i<=M; i++)
   for (j=1; j<=Ntilde; j++)
     zch[(i-1)*Ntilde + j-1] = b[i-1] + (b[i]-b[i-1])/2.*(1.+cos(PI*(2.*(Ntilde-j+1.)-1.)/2./dN));
   
 for (i=1; i<=M; i++)
   for (j=1; j<=Ntilde ;j++) {
     qi = (i-1)*Ntilde + j-1;
     za = zch[(i-1)*Ntilde + j-1] - refk;
     zb = zch[(i-1)*Ntilde + j-1] + eps;

     for (ii=1; ii<=M; ii++) {
       if ( za>b[ii-1] ) xl = za; else xl = b[ii-1];
       if ( zb<b[ii] ) xu = zb; else xu = b[ii];
       if ( df!=2 && xu>za ) { xl = sqrt(xl-za); xu = sqrt(xu-za); }

       for (jj=1; jj<=Ntilde; jj++) {
	 qj = (ii-1)*Ntilde + jj-1;
         if ( xu<xl ) a[qi*NN + qj] = 0.;
         else {
	   gausslegendre(qm, xl, xu, z, w);
           Hij = 0.;
           for (k=0; k<qm; k++)
             if ( df==2 )
               Hij += w[k]*Tn((2.*z[k]-b[ii]-b[ii-1])/(b[ii]-b[ii-1]),jj-1)*exp(-z[k]/s2);
             else
               Hij += w[k]*Tn((2.*(za+z[k]*z[k])-b[ii]-b[ii-1])/(b[ii]-b[ii-1]),jj-1)
                 * 2. * pow(z[k], ddf-1.) * exp(-alpha*z[k]*z[k]);
           if ( df==2 ) Hij *= exp(za/s2)/s2;
           else Hij *= pow(alpha,ddf/2.)/gammafn(ddf/2.);
           a[qi*NN + qj] = -Hij;
         }
       }
     }

     for (jj=1; jj<=Ntilde; jj++)
       a[qi*NN + jj-1] -= CHI(-ddf/s2*za, df) * Tn(-1.,jj-1);

     for (jj=1; jj<=Ntilde; jj++)
       a[qi*NN + (i-1)*Ntilde + jj-1] += Tn((2.*zch[(i-1)*Ntilde + j-1]-b[i]-b[i-1])/(b[i]-b[i-1]),jj-1);
   }

 for (j=0;j<NN;j++) g[j] = 1.;

 LU_solve(a, g, NN);

 arl = 0.;
 for (j=1; j<=Ntilde; j++)
   arl += g[(ihs-1)*Ntilde + j-1] * Tn((2.*hs-b[ihs]-b[ihs-1])/(b[ihs]-b[ihs-1]),j-1);

 Free(zch);
 Free(z);
 Free(w);
 Free(b1);
 Free(b2);
 Free(b);
 Free(g);
 Free(a);
 
 return arl;
}


/* double Cs2arlGCRpw */
double scU_iglarl_v2(double refk, double h, double hs, double sigma, int df, int N, int qm)
{ double *a, *g, *w, *z, arl, Hij, xl, za, dN, ddf, s2, *t, t0, t1, th, x0, x1, dummy;
  int i, j, k, M, Ntilde, NN, ii, jj, it, qi, qj;
 
 M = ceil( h/refk );
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde;
 
 a = matrix(NN,NN);
 g = vector(NN);
 w = vector(qm);
 z = vector(qm);
 t = vector(NN);

 /* Chebyshev Gauss-Lobatto nodes */ 
 for(i=1; i<=M; i++) {
   t0 = (double)(i-1.)*refk;
   t1 = t0 + refk;
   if ( t1>h ) t1 = h;
   for (j=1; j<Ntilde; j++) {
     th = cos( PI/(dN-1.) * (dN-j-1.) );
     t[(i-1)*(Ntilde-1)+j] = t0 + (th+1.)/2.*(t1-t0);
   }
 }
 t[0] = 0.;

 for (i=1; i<=M; i++)
   for (j=1; j<=Ntilde ;j++) {
     qi = (i-1)*Ntilde + j-1;  it = (i-1)*(Ntilde-1) + j-1;
     
     za = t[it] - refk;
     if ( za<0. ) xl = 0.; else xl = za;

     for (ii=1; ii<i-1; ii++)
       for (jj=1; jj<=Ntilde; jj++) {
         qj = (ii-1)*Ntilde + jj-1;
         a[qi*NN + qj] = 0.;
       } /* ii = 1 .. i-2, jj = 1 .. Ntilde */       
      
     if ( i>1 ) {
       ii = i-1;
       t0 = (double)(ii-1.)*refk;
       t1 = t0 + refk;
       if ( t1>h ) t1 = h;
       if ( t0<xl ) x0 = xl; else x0 = t0;
       if ( df!=2 ) {
         if ( x0-za>1e-10 ) x0 = sqrt(x0-za); else x0 = 0.;
         if ( t1-za>1e-10 ) x1 = sqrt(t1-za); else x1 = 0.; }
       else x1 = t1;
       
       for (jj=1; jj<=Ntilde; jj++) {
         qj = (ii-1)*Ntilde + jj-1;

         if ( j==1 ) a[qi*NN + qj] = - Tn((2.*t[it]-t0-t1)/(t1-t0),jj-1);
         else {
           if ( fabs(x1-x0)>1e-12 ) {
             gausslegendre(qm, x0, x1, z, w);
             Hij = 0.;
             for (k=0; k<qm; k++) {
               if ( df==2 )
                 Hij += w[k] * Tn((2.*z[k]-t0-t1)/(t1-t0),jj-1) * ddf/s2*chi(ddf/s2*(z[k]-za),df);
               else
                 Hij += w[k] * Tn((2.*(z[k]*z[k]+za)-t0-t1)/(t1-t0),jj-1)*2.*pow(z[k],ddf-1.)*exp(-ddf*z[k]*z[k]/2./s2); 
             } /* k = 0 .. qm-1 */
             if ( df!=2 ) Hij /= gammafn(ddf/2.) * pow(2.*s2/ddf,ddf/2.); 
             a[qi*NN + qj] = - Hij;
           }
           else a[qi*NN + qj] = 0.;
         } /* j != 1*/ 
       } /* jj = 1 .. Ntilde */
     } /* i > 1 */
      
     for (ii=i; ii<=M; ii++) {
       t0 = (double)(ii-1.)*refk;
       t1 = t0 + refk;
       if ( t1>h ) t1 = h;
       if ( t0<xl ) x0 = xl; else x0 = t0;
       if ( df!=2 ) {
         if ( x0-za>1e-10 ) x0 = sqrt(x0-za); else x0 = 0.;
         if ( t1-za>1e-10 ) x1 = sqrt(t1-za); else x1 = 0.; }
       else x1 = t1;
	
       if ( i>1 && j==1 && ii==i ) {
         for (jj=1; jj<=Ntilde; jj++) {
           qj = (ii-1)*Ntilde + jj-1;
           a[qi*NN + qj] = Tn((2.*t[it]-t0-t1)/(t1-t0),jj-1);
         } /* jj = 1 .. Ntilde */
       } /* i>1 && j==1 && ii==i */

       if ( i>1 && j==1 && ii>i ) {
         for (jj=1; jj<=Ntilde; jj++) {
           qj = (ii-1)*Ntilde + jj-1;
           a[qi*NN + qj] = 0.;
         } /* jj = 1 .. Ntilde */
       } /* i>1 && j==1 && ii>i */

       if ( i==1 || j>1 ) {
         for (jj=1; jj<=Ntilde; jj++) {
           qj = (ii-1)*Ntilde + jj-1;
           gausslegendre(qm, x0, x1, z, w);
           Hij = 0.;
           for (k=0; k<qm; k++)
             if ( df==2 )
               Hij += w[k] * Tn((2.*z[k]-t0-t1)/(t1-t0),jj-1) * ddf/s2*chi(ddf/s2*(z[k]-za),df);
             else
               Hij += w[k] * Tn((2.*(z[k]*z[k]+za)-t0-t1)/(t1-t0),jj-1)*2.*pow(z[k],ddf-1.)*exp(-ddf*z[k]*z[k]/2./s2);
           if ( df!=2 ) Hij /= gammafn(ddf/2.) * pow(2.*s2/ddf,ddf/2.);
           if ( ii==i ) a[qi*NN + qj] = Tn((2.*t[it]-t0-t1)/(t1-t0),jj-1) - Hij;
           else a[qi*NN + qj] = -Hij;
         } /* jj = 1 .. Ntilde */
       } /* i==1 || j>1 */
     } /* ii = i .. M */
        
     if ( i==1 ) {
       t0 = 0.;
       t1 = refk;
       if ( t1>h ) t1 = h;
       for (jj=1; jj<=Ntilde; jj++) {
         dummy = -za/s2;
         if ( dummy>0. ) {
           if ( df==1 ) dummy = 2.*PHI( sqrt(dummy), 0. ) - 1.;
           if ( df==2 ) dummy = 1. - exp( -dummy );
           if ( df>2  ) dummy = CHI( ddf*dummy, df);
         }
         else dummy = 0.;
         a[qi*NN + jj-1] -= dummy * Tn(-1.,jj-1);
       } /* jj = 1 .. Ntilde */
     } /* i==1 */ 
      
   } /* i = 1 .. M, j = 1 .. Ntilde */
    
 for (j=0; j<NN; j++) g[j] = 1.;
 for (j=1; j<M; j++)  g[Ntilde*j] = 0.;
 
 LU_solve(a, g, NN);

 arl = 0.;
 for (i=1; i<=M; i++) {
    t0 = (double)(i-1.)*refk;
    t1 = t0 + refk;
    if ( t1>h ) t1 = h;
    if ( t0<=hs && hs<t1 )
      for (j=1; j<=Ntilde; j++) {
        ii = (i-1)*Ntilde + j-1;
        arl += g[ii] * Tn((2.*hs-t0-t1)/(t1-t0),j-1);
      } /* j = 1 .. Ntilde */
 } /* i = 1 .. M */ 

 Free(t);
 Free(z);
 Free(w);
 Free(g);
 Free(a);
 
 return arl;
}


/* double lCs2arlGCRpw */
double scL_iglarl_v2(double refk, double h, double hs, double sigma, int df, int N, int qm)
{ double *a, *g, *w, *z, arl, Hij, xu, za, dN, ddf, s2, *t, t0, t1, th, x0, x1, dummy;
  int i, j, k, M, Ntilde, NN, ii, jj, it, qi, qj, imax;
 
 M = ceil( h/refk );
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde;
  
 a = matrix(NN,NN);
 g = vector(NN);
 w = vector(qm);
 z = vector(qm);
 t = vector(NN);

 /* Chebyshev Gauss-Lobatto nodes */ 
 for(i=1; i<=M; i++) {
   t0 = h - (double)(M-i+1.)*refk;
   t1 = t0 + refk;
   if ( t0<0. ) t0 = 0.;
   for (j=1; j<Ntilde; j++) {
     th = cos( PI/(dN-1.) * (dN-j-1.) );
     t[(i-1)*(Ntilde-1)+j] = t0 + (th+1.)/2.*(t1-t0);
   }
 }
 t[0] = 0.;

 for (i=1; i<=M; i++)
   for (j=1; j<=Ntilde ;j++) {
     qi = (i-1)*Ntilde + j-1;  it = (i-1)*(Ntilde-1) + j-1;
     
     za = t[it] + refk;
     if ( za<h ) xu = za; else xu = h;

     imax = i+1; if ( imax>M ) imax = M;
     for (ii=1; ii<=imax; ii++) {
       t0 = h - (double)(M-ii+1.)*refk;
       t1 = t0 + refk;
       if ( t0<0. ) t0 = 0.;
       if ( t1<xu ) x1 = t1; else x1 = xu;
       
       if ( df!=2 ) {
         if ( za-x1>1e-10 ) x0 = sqrt(za-x1); else x0 = 0.;
         if ( za-t0>1e-10 ) x1 = sqrt(za-t0); else x1 = 0.;
       }
       else x0 = t0;
       
       for (jj=1; jj<=Ntilde; jj++) {
         qj = (ii-1)*Ntilde + jj-1;

         if ( i>1 && j==1 ) { /* continuity condition */
           if ( ii==i-1 ) a[qi*NN + qj] = - Tn((2.*t[it]-t0-t1)/(t1-t0),jj-1);
           if ( ii==i )   a[qi*NN + qj] =   Tn((2.*t[it]-t0-t1)/(t1-t0),jj-1);
           if ( ii<i-1 || ii>i) a[qi*NN + qj] = 0.;
         } else {
           gausslegendre(qm, x0, x1, z, w);
           Hij = 0.;
           for (k=0; k<qm; k++)
             if (df==2) 
               Hij += w[k] * Tn((2.*z[k]-t0-t1)/(t1-t0),jj-1) * ddf/s2*chi(ddf/s2*(za-z[k]),df);
             else
               Hij += w[k] * Tn((2.*(za-z[k]*z[k])-t0-t1)/(t1-t0),jj-1) *2.*pow(z[k],ddf-1.)*exp(-ddf*z[k]*z[k]/2./s2);
           if ( df!=2 ) Hij /= gammafn(ddf/2.) * pow(2.*s2/ddf,ddf/2.);
           if ( ii==i ) a[qi*NN + qj] = Tn((2.*t[it]-t0-t1)/(t1-t0),jj-1) - Hij;
           else a[qi*NN + qj] = -Hij;
         } /* (! i>1 && j==1) */
       } /* jj = 1 .. Ntilde */
     } /* ii = 1 .. imax <= M */
     
     for (ii=i+2; ii<=M; ii++)
       for (jj=1; jj<=Ntilde; jj++) {
         qj = (ii-1)*N + jj-1;
         a[qi*NN + qj] = 0.;
       }

     if ( i==1 || j>1 ) {
       for ( jj=1; jj<=Ntilde; jj++) { /* ii = 1  -- atom */
          dummy = za/s2;
          if ( df==1 ) dummy = 2.*( 1. - PHI( sqrt(dummy), 0. ) );
          if ( df==2 ) dummy = exp( -dummy );
          if (  df>2 ) dummy = 1. - CHI( ddf*dummy, df);
          a[qi*NN + jj-1] -= dummy * Tn(-1.,jj-1);
        } /* jj = 1 .. Ntilde */
      } /* i==1 || j>1 */       
   } /* i = 1 .. M, j = 1 .. Ntilde */
    
 for (j=0; j<NN; j++) g[j] = 1.;
 for (j=1; j<M; j++)  g[Ntilde*j] = 0.;
 
 LU_solve(a, g, NN);

 arl = 0.;
 for (i=1; i<=M; i++) {
    t0 = h - (double)(M-i+1.)*refk;
    t1 = t0 + refk;
    if ( t0<0. ) t0 = 0.;
    if ( t0<=hs && hs<t1 )
      for (j=1; j<=Ntilde; j++) {
        ii = (i-1)*Ntilde + j-1;
        arl += g[ii] * Tn((2.*hs-t0-t1)/(t1-t0),j-1);
      } /* j = 1 .. Ntilde */
 } /* i = 1 .. M */ 

 Free(t);
 Free(z);
 Free(w);
 Free(g);
 Free(a);
 
 return arl;
}


double sc2_iglarl_v2(double refkl, double refku, double hl, double hu, double hsl, double hsu, double sigma, int df, int N, int qm)
{ double arl1, arl2, arl3, arl4, arl;

/* relation between 1- and 2-sided CUSUM schemes due to Lucas/Crosier 1982,
   Technometrics 24, 199-205;
   only for headstart hs smaller than h/2 + k !!
   Chang/Gan 1995 claim that it is valid also for 2-sided S^2 CUSUM
   (JQT 27(2), 109-119
*/

 arl1 = scU_iglarl_v2(refku, hu,  0., sigma, df, N, qm);
 arl2 = scU_iglarl_v2(refku, hu, hsu, sigma, df, N, qm);
 arl3 = scL_iglarl_v2(refkl, hl,  0., sigma, df, N, qm);
 arl4 = scL_iglarl_v2(refkl, hl, hsl, sigma, df, N, qm);
 arl = ( arl2*arl3 + arl1*arl4 - arl1*arl3 ) / ( arl1 + arl3 );
 return arl;
}


double scU_crit(double refk, double L0, double hs, double sigma, int df, int N, int qm)
{ double c1, c2, c3, L1=0., L2=0., L3=0., dc;

 c2 = 0.;
 L2 = 1.;
 do {
   c1 = c2;
   L1 = L2;
   c2 += 1.;
   L2 = scU_iglarl_v2(refk, c2, hs, sigma, df, N, qm);
 } while ( L2<L0 );
 
 do {
   c3 = c1 + (L0-L1)/(L2-L1) * (c2-c1);
   L3 = scU_iglarl_v2(refk, c3, hs, sigma, df, N, qm);
   dc = c3-c2; c1 = c2; L1 = L2; c2 = c3; L2 = L3;
 } while ( (fabs(L0-L3)>1e-6) && (fabs(dc)>1e-9) );
 return c3;
}


double scL_crit(double refk, double L0, double hs, double sigma, int df, int N, int qm)
{ double c1, c2, c3, L1=0., L2=0., L3=0., dc;

 c2 = 0.;
 L2 = 1.;
 do {
   c1 = c2;
   L1 = L2;
   c2 += 1;
   L2 = scL_iglarl_v2(refk, c2, hs, sigma, df, N, qm);
 } while ( L2<L0 );
 
 do {
   c3 = c1 + (L0-L1)/(L2-L1) * (c2-c1);
   L3 = scL_iglarl_v2(refk, c3, hs, sigma, df, N, qm);
   dc = c3-c2; c1 = c2; L1 = L2; c2 = c3; L2 = L3;
 } while ( (fabs(L0-L3)>1e-6) && (fabs(dc)>1e-9) );
 return c3;
}


double scL_fu_crit(double refkl, double refku, double hu, double L0, double hsl, double hsu, double sigma, int df, int N, int qm)
{ double c1, c2, c3, L1=0., L2=0., L3=0., dc;

 c2 = 0.;
 L2 = 1.;
 do {
   c1 = c2;
   L1 = L2;
   c2 += 1;
   L2 = sc2_iglarl_v2(refkl, refku, c2, hu, hsl, hsu, sigma, df, N, qm);
 } while ( L2<L0 );
 
 do {
   c3 = c1 + (L0-L1)/(L2-L1) * (c2-c1);
   L3 = sc2_iglarl_v2(refkl, refku, c3, hu, hsl, hsu, sigma, df, N, qm);
   dc = c3-c2; c1 = c2; L1 = L2; c2 = c3; L2 = L3;
 } while ( (fabs(L0-L3)>1e-6) && (fabs(dc)>1e-9) );
 return c3;
}


double scU_fl_crit(double refkl, double refku, double hl, double L0, double hsl, double hsu, double sigma, int df, int N, int qm)
{ double c1, c2, c3, L1=0., L2=0., L3=0., dc;

 c2 = 0.;
 L2 = 1.;
 do {
   c1 = c2;
   L1 = L2;
   c2 += 1;
   L2 = sc2_iglarl_v2(refkl, refku, hl, c2, hsl, hsu, sigma, df, N, qm);
 } while ( L2<L0 );
 
 do {
   c3 = c1 + (L0-L1)/(L2-L1) * (c2-c1);
   L3 = sc2_iglarl_v2(refkl, refku, hl, c3, hsl, hsu, sigma, df, N, qm);
   dc = c3-c2; c1 = c2; L1 = L2; c2 = c3; L2 = L3;
 } while ( (fabs(L0-L3)>1e-6) && (fabs(dc)>1e-9) );
 return c3;
}


int sc2_crit_unbiased(double refkl, double refku, double L0, double *hl, double *hu, double hsl, double hsu, double sigma, int df, int N, int qm)
{ double h1, h2, h3, dh, lh, sl1, sl2, sl3, Lm, Lp, step;

 step = .2/sqrt(df); 
 
 h1 = scU_crit(refku, 2.*L0, hsu, sigma, df, N, qm);
 lh = scL_crit(refkl, 2.*L0, hsl, sigma, df, N, qm);
 Lm = sc2_iglarl_v2(refkl, refku, lh, h1, hsl, hsu, sigma-lmEPS, df, N, qm);
 Lp = sc2_iglarl_v2(refkl, refku, lh, h1, hsl, hsu, sigma+lmEPS, df, N, qm);
 sl1 = (Lp-Lm)/(2.*lmEPS);
 
 h2 = h1;
 sl2 = sl1;
 do { 
   h1 = h2;
   sl1 = sl2;
   h2 = h1 + step;
   lh = scL_fu_crit(refkl, refku, h2, L0, hsl, hsu, sigma, df, N, qm);   
   Lm = sc2_iglarl_v2(refkl, refku, lh, h2, hsl, hsu, sigma-lmEPS, df, N, qm);
   Lp = sc2_iglarl_v2(refkl, refku, lh, h2, hsl, hsu, sigma+lmEPS, df, N, qm);   
   sl2 = (Lp-Lm)/(2.*lmEPS);
 } while ( sl2 < 0. );

 do {
   h3 = h1 - sl1/(sl2-sl1) * (h2-h1);
   lh = scL_fu_crit(refkl, refku, h3, L0, hsl, hsu, sigma, df, N, qm);   
   Lm = sc2_iglarl_v2(refkl, refku, lh, h3, hsl, hsu, sigma-lmEPS, df, N, qm);
   Lp = sc2_iglarl_v2(refkl, refku, lh, h3, hsl, hsu, sigma+lmEPS, df, N, qm);      
   sl3 = (Lp-Lm)/(2.*lmEPS);
   dh = h3-h2; h1 = h2; sl1 = sl2; h2 = h3; sl2 = sl3;
 } while ( fabs(sl3)>1e-7 && fabs(dh)>1e-9 );

 *hl = lh; *hu = h3;

 return 0;
}


/* MPT = Moustakides/Polunchenko/Tartakovsky */
double xsr1_iglarl(double k, double h, double zr, double hs, double mu, int N, int MPT)
{ double *a, *g, *w, *z, arl, adjust=1.;
  int i, j, NN;

 adjust = 1.;
 if ( MPT ) adjust = 2.*k;
  
 NN = N + 1;
 a = matrix(NN,NN);
 g = vector(NN);
 w = vector(NN);
 z = vector(NN);

 gausslegendre(N, zr, h, z, w);

 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) 
     a[i*NN+j] = -w[j] * phi( (z[j]-log(1.+exp(z[i])))/adjust + k, mu)/adjust;    
   ++a[i*NN+i];
   a[i*NN+N] = - PHI( (zr-log(1.+exp(z[i])))/adjust + k, mu);
 }
 for (j=0;j<N;j++)
    a[N*NN+j] = -w[j] * phi( (z[j]-log(1.+exp(zr)))/adjust + k, mu)/adjust;
 a[N*NN+N] = 1. - PHI( (zr-log(1.+exp(zr)))/adjust + k, mu);

 for (j=0;j<NN;j++) g[j] = 1.;
 LU_solve(a,g,NN);

 if (hs > h) {
   arl = 1. + PHI( zr/adjust + k, mu) * g[N];
   for (j=0;j<N;j++)
     arl += w[j] * phi( z[j]/adjust + k, mu)/adjust * g[j];
 } else {
   arl = 1. + PHI( (zr-log(1.+exp(hs)))/adjust + k, mu) * g[N];
   for (j=0;j<N;j++)
     arl += w[j] * phi( (z[j]-log(1.+exp(hs)))/adjust + k, mu)/adjust * g[j];
 }

 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}


double xsr1_arlm(double k, double h, double zr, double hs, int q, double mu0, double mu1, int N, int nmax, int MPT)
{ double *p0, *fn, *w, *z, arl0, rho, arl_minus=0., arl, arl_plus=0., mn_minus=0., mn_plus=0., ratio=0., adjust=1.;
  int i, j, n, NN;
 
 adjust = 1.;
 if ( MPT ) adjust = 2.*k; 

 NN = N + 1;
 w   = vector(NN);
 z   = vector(NN);
 fn  = matrix(nmax, NN);
 p0  = vector(nmax);

 gausslegendre(N, zr, h, z, w);

 /* in-control, i. e. n<=q-1 */
 for (n=1; n<=q-1; n++) {
  /* determine f_n, n=1,2,...,q-1 */
  if (n==1) {
    if ( hs > h ) {
      for (i=0; i<N; i++) fn[0*NN+i] = phi( z[i]/adjust + k, mu0)/adjust;
      fn[0*NN+N] = PHI( zr/adjust + k, mu0);
    } else {
      for (i=0; i<N; i++) fn[0*NN+i] = phi( (z[i]-log(1.+exp(hs)))/adjust + k, mu0)/adjust;
      fn[0*NN+N] = PHI( (zr-log(1.+exp(hs)))/adjust + k, mu0);
    }
  } else {
    for (i=0; i<N; i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * phi( (z[i]-log(1.+exp(zr)))/adjust + k, mu0)/adjust;
      for (j=0; j<N; j++) {
        fn[(n-1)*NN+i] += w[j] * fn[(n-2)*NN+j] * phi( (z[i]-log(1.+exp(z[j])))/adjust + k, mu0)/adjust;
      }
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI( (zr-log(1.+exp(zr)))/adjust + k, mu0);
    for (j=0; j<N; j++) fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI( (zr-log(1.+exp(z[j])))/adjust + k, mu0);
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0; i<N; i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];
 }

 /* out-of-control, i.e. t>=q */
 arl0 = 1.; rho = 0.;

 for (n=q; n<=nmax; n++) {
  if ( n==1 ) {
    if ( hs > h ) {
      for (i=0; i<N; i++) fn[0*NN+i] = phi( z[i]/adjust + k, mu1)/adjust;
      fn[0*NN+N] = PHI( zr/adjust + k, mu1);
    } else {
      for (i=0; i<N; i++) fn[0*NN+i] = phi( (z[i]-log(1.+exp(hs)))/adjust + k, mu1)/adjust;
      fn[0*NN+N] = PHI( (zr-log(1.+exp(hs)))/adjust + k, mu1);
    }
  } else {
    for (i=0; i<N; i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * phi( (z[i]-log(1.+exp(zr)))/adjust + k, mu1)/adjust;
      for (j=0; j<N; j++) {
        fn[(n-1)*NN+i] += w[j] * fn[(n-2)*NN+j] * phi( (z[i]-log(1.+exp(z[j])))/adjust + k, mu1)/adjust;
      }
      if ( n==q && q>1 ) fn[(n-1)*NN+i] /= p0[q-2];
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI( (zr-log(1.+exp(zr)))/adjust + k, mu1);
    for (j=0; j<N; j++) fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI( (zr-log(1.+exp(z[j])))/adjust + k, mu1);
    if ( n==q && q>1 ) fn[(n-1)*NN+N] /= p0[q-2];
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];

  /* computation of m_n+1^- and m_n+1^+, n=m-1,m,... */
  mn_minus = 1.; mn_plus = 0.;
  if ( n > q ) {
    for (i=0; i<NN; i++) {
      if (fn[(n-2)*NN+i]==0)
      if (fn[(n-1)*NN+i]==0) ratio = 0.; else ratio = 1.;
      else ratio = fn[(n-1)*NN+i]/fn[(n-2)*NN+i];
      if ( ratio<mn_minus ) mn_minus = ratio;
      if ( ratio>mn_plus ) mn_plus = ratio;
    }
    rho = p0[n-1]/p0[n-2];
  }

  /* computation of ARL, ARL^-, and ARL^+ */
  arl = arl0 + p0[n-1]/(1.-rho);
  if (mn_minus<1.) arl_minus = arl0 + p0[n-1]/(1.-mn_minus);
  else             arl_minus = -1.;
  if (mn_plus<1.) arl_plus = arl0 + p0[n-1]/(1.-mn_plus);
  else            arl_plus = -1.;
  arl0 += p0[n-1];

  if ( fabs((arl_plus-arl_minus)) < 1e-7 ) n = nmax+1;
 }

 arl = (arl_plus+arl_minus)/2.; rho0 = rho;

 Free(w);
 Free(z);
 Free(fn);
 Free(p0);

 return arl;
}


double xsr1_arlm_hom(double k, double h, double zr, double hs, int q, double mu0, double mu1, int N, int MPT, double *ced)
{ double *fn, *w, *z, *a, *arl, adjust=1., norm;
  int i, j, n, NN;
 
 adjust = 1.;
 if ( MPT ) adjust = 2.*k; 

 NN = N + 1;
 w   = vector(NN);
 z   = vector(NN);
 fn  = matrix(q+1, NN);
 a   = matrix(NN,NN);
 arl = vector(NN);
 
 gausslegendre(N, zr, h, z, w); 
 
 /* ARL vector */
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*NN+j] = -w[j] * phi( (z[j]-log(1.+exp(z[i])))/adjust + k, mu1)/adjust;
   ++a[i*NN+i];
   a[i*NN+N] = - PHI( (zr-log(1.+exp(z[i])))/adjust + k, mu1);
 }
 for (j=0; j<N; j++)
    a[N*NN+j] = -w[j] * phi( (z[j]-log(1.+exp(zr)))/adjust + k, mu1)/adjust;
 a[N*NN+N] = 1. - PHI( (zr-log(1.+exp(zr)))/adjust + k, mu1);

 for (j=0; j<NN; j++) arl[j] = 1.;
 LU_solve(a, arl, NN);

 /* q == 1 */
  if ( hs > h ) {
   ced[0] = 1. + PHI( zr/adjust + k, mu1) * arl[N];
   for (j=0; j<N; j++)
     ced[0] += w[j] * phi( z[j]/adjust + k, mu1)/adjust * arl[j];
 } else {
   ced[0] = 1. + PHI( (zr-log(1.+exp(hs)))/adjust + k, mu1) * arl[N];
   for (j=0; j<N; j++)
     ced[0] += w[j] * phi( (z[j]-log(1.+exp(hs)))/adjust + k, mu1)/adjust * arl[j];
 }
 
 /* density sequence for q > 1 */
 for (n=1; n<=q-1; n++) {
   if ( n == 1 ) {
     if ( hs > h ) {
       for (i=0; i<N; i++) fn[0*NN+i] = phi( z[i]/adjust + k, mu0)/adjust;
       fn[0*NN+N] = PHI( zr/adjust + k, mu0);
     } else {
       for (i=0; i<N; i++) fn[0*NN+i] = phi( (z[i]-log(1.+exp(hs)))/adjust + k, mu0)/adjust;
       fn[0*NN+N] = PHI( (zr-log(1.+exp(hs)))/adjust + k, mu0);
     }
   } else {
     for (i=0; i<N; i++) {
       fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * phi( (z[i]-log(1.+exp(zr)))/adjust + k, mu0)/adjust;
       for (j=0; j<N; j++) {
         fn[(n-1)*NN+i] += w[j] * fn[(n-2)*NN+j] * phi( (z[i]-log(1.+exp(z[j])))/adjust + k, mu0)/adjust;
       }
     }
     fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI( (zr-log(1.+exp(zr)))/adjust + k, mu0);
     for (j=0; j<N; j++) fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI( (zr-log(1.+exp(z[j])))/adjust + k, mu0);    
   }
  
   ced[n] = fn[(n-1)*NN+N] * arl[N];
   norm = fn[(n-1)*NN+N];
   for (j=0; j<N; j++) {
     ced[n] += w[j] * fn[(n-1)*NN+j] * arl[j]; 
     norm += w[j] * fn[(n-1)*NN+j];
   }
   ced[n] /= norm;
 }
 
 Free(w);
 Free(z);
 Free(fn);
 Free(a);
 Free(arl);

 return 0;
}


double xsr1_iglarl_drift(double k, double h, double zr, double hs, double delta, int m, int N, int with0)
{ double *a, *g, *w, *z, arl, *MUs, *ARLs;
  int i, j, NN, m_;

 NN = N + 1;
 a = matrix(NN, NN);
 g = vector(NN);
 w = vector(NN);
 z = vector(NN);
 ARLs = vector(NN);
 MUs  = vector(m+1); 

 gausslegendre(N, zr, h, z, w);

 if ( with0 ) {
   for (i=0;i<=m;i++) MUs[i] = (double)i * delta;
 } else {
   for (i=0;i<=m;i++) MUs[i] = (double)(i+1.) * delta;
 }

 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) a[i*NN+j] = -w[j] * phi( z[j]-log(1.+exp(z[i]))+k, MUs[m]);
   ++a[i*NN+i];
   a[i*NN+N] = - PHI( zr-log(1.+exp(z[i]))+k, MUs[m]);
 }
 for (j=0;j<N;j++)
    a[N*NN+j] = -w[j] * phi( z[j]-log(1.+exp(zr))+k, MUs[m]);
 a[N*NN+N] = 1. - PHI( zr-log(1.+exp(zr))+k, MUs[m]);

 for (j=0;j<NN;j++) g[j] = 1.;
 LU_solve(a, g, NN);

 for (m_=0;m_<m;m_++) {
   for (i=0;i<=N;i++) {
     ARLs[i] = 1. + PHI( zr-log(1.+exp(z[i]))+k, MUs[m-m_]) * g[N];
     for (j=0;j<N;j++) { 
       ARLs[i] += w[j] * phi( z[j]-log(1.+exp(z[i]))+k, MUs[m-m_]) * g[j];
     }
   }
   for (j=0;j<=N;j++) g[j] = ARLs[j];
 }

 if (hs > h) {
   arl = 1. + PHI( zr+k, MUs[0]) * ARLs[N];
   for (j=0;j<N;j++) arl += w[j] * phi( z[j]+k, MUs[0]) * ARLs[j];
 } else {
   arl = 1. + PHI( zr-log(1.+exp(hs))+k, MUs[0]) * ARLs[N];
   for (j=0;j<N;j++) arl += w[j] * phi( z[j]-log(1.+exp(hs))+k, MUs[0]) * ARLs[j];
 }

 Free(a);
 Free(g);
 Free(w);
 Free(z);
 Free(ARLs);
 Free(MUs);

 return arl;
}


double xsr1_iglarl_drift_wo_m(double k, double h, double zr, double hs, double delta, int *m, int N, int with0)
{ int m_;
  double arl1, arl2, eps=1e-6;
 m_ = 4;
 arl1 = xsr1_iglarl_drift(k, h, zr, hs, delta, m_, N, with0);
 arl2 = arl1 + 2.*eps;
 while ( fabs(arl2-arl1)>eps && (double)m_<1e4 ) {
   m_ = (int)round(1.5 * m_);
   arl1 = xsr1_iglarl_drift(k, h, zr, hs, delta, m_, N, with0);
   arl2 = xsr1_iglarl_drift(k, h, zr, hs, delta, m_+1, N, with0);
 }
 *m = m_;
 return arl1;
}


double xsr1_iglarlm_drift(double k, double h, double zr, double hs, int q, double delta, int N, int nmax, int with0)
{ double *p0, *fn, *w, *z, arl0, rho, MEAN=0.,
         arl_minus=0., arl, arl_plus=0., mn_minus=0., mn_plus=0., nn, ratio=0.;
  int i, j, n, NN;

 NN = N + 1;
 w   = vector(NN);
 z   = vector(NN);
 fn  = matrix(nmax, NN);
 p0  = vector(nmax);

 gausslegendre(N, zr, h, z, w);

 /* in-control, i. e. n<=q-1 */
 MEAN = 0.;

 for (n=1;n<=q-1;n++) {
  nn = (double) n;

  /* determine f_n, n=1,2,...,q-1 */
  if (n==1) {
    for (i=0;i<N;i++) fn[0*NN+i] = phi( z[i]-log(1.+exp(hs))+k, MEAN);
    fn[0*NN+N] = PHI( zr-log(1.+exp(hs))+k, MEAN);
  } else {
    for (i=0;i<N;i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * phi( z[i]-log(1.+exp(zr))+k, MEAN);
      for (j=0;j<N;j++) {
        fn[(n-1)*NN+i] += w[j] * fn[(n-2)*NN+j] * phi( z[i]-log(1.+exp(z[j]))+k, MEAN);
      }
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI( zr-log(1.+exp(zr))+k, MEAN);
    for (j=0;j<N;j++) fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI( zr-log(1.+exp(z[j]))+k, MEAN);
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];
 }

 /* out-of-control, i.e. t>=q */
 arl0 = 1.; rho = 0.;

 for (n=q;n<=nmax;n++) {
  nn = (double) n;

  /* determine f_n, n=q,q+1,... */
  if ( with0 ) {
    MEAN = (nn-(double)q) * delta;
  } else {
    MEAN = (nn-(double)q+1.) * delta;
  }

  if (n==1) {
    for (i=0;i<N;i++) fn[0*NN+i] = phi( z[i]-log(1.+exp(hs))+k, MEAN);
    fn[0*NN+N] = PHI( zr-log(1.+exp(hs))+k, MEAN);
  } else {
    for (i=0;i<N;i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * phi( z[i]-log(1.+exp(zr))+k, MEAN);
      for (j=0;j<N;j++) {
        fn[(n-1)*NN+i] += w[j] * fn[(n-2)*NN+j] * phi( z[i]-log(1.+exp(z[j]))+k, MEAN);
      }
      if (n==q && q>1) fn[(n-1)*NN+i] /= p0[q-2];
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI( zr-log(1.+exp(zr))+k, MEAN);
    for (j=0;j<N;j++) fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI( zr-log(1.+exp(z[j]))+k, MEAN);
    if (n==q && q>1) fn[(n-1)*NN+N] /= p0[q-2];
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];

  /* computation of m_n+1^- and m_n+1^+, n=m-1,m,... */
  mn_minus = 1.; mn_plus = 0.;
  if (n>q) {
   for (i=0;i<NN;i++) {
    if (fn[(n-2)*NN+i]==0)
     if (fn[(n-1)*NN+i]==0) ratio = 0.; else ratio = 1.;
    else ratio = fn[(n-1)*NN+i]/fn[(n-2)*NN+i];
    if ( ratio<mn_minus ) mn_minus = ratio;
    if ( ratio>mn_plus ) mn_plus = ratio;
   }
  }

  if (n>q) rho = p0[n-1]/p0[n-2];

  /* computation of ARL, ARL^-, and ARL^+ */
  arl = arl0 + p0[n-1]/(1.-rho);
  if (mn_minus<1.) arl_minus = arl0 + p0[n-1]/(1.-mn_minus);
  else             arl_minus = -1.;
  if (mn_plus<1.) arl_plus = arl0 + p0[n-1]/(1.-mn_plus);
  else            arl_plus = -1.;
  arl0 += p0[n-1];

  if ( fabs((arl_plus-arl_minus)) < 1e-7 ) n = nmax+1;
 }

 arl = (arl_plus+arl_minus)/2.; rho0 = rho;

 Free(w);
 Free(z);
 Free(fn);
 Free(p0);

 return arl;
}


double xsr1_iglad(double k, double h, double zr, double mu0, double mu1, int N, int MPT)
{ double *a, *w, *z, *arl, *psi, rho, ad, norm, adjust=1.;
  int i, j, status, noofit, NN;

 adjust = 1.;
 if ( MPT ) adjust = 2.*k;
   
 NN = N + 1;
 a = matrix(NN,NN);
 arl = vector(NN);
 psi = vector(NN);
 w = vector(NN);
 z = vector(NN);

 gausslegendre(N, zr, h, z, w);

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*NN+j] = -w[j] * phi( (z[j]-log(1.+exp(z[i])))/adjust + k, mu1)/adjust;
   ++a[i*NN+i];
   a[i*NN+N] = - PHI( (zr-log(1.+exp(z[i])))/adjust + k, mu1);
 }
 for (j=0; j<N; j++)
    a[N*NN+j] = -w[j] * phi( (z[j]-log(1.+exp(zr)))/adjust + k, mu1)/adjust;
 a[N*NN+N] = 1. - PHI( (zr-log(1.+exp(zr)))/adjust + k, mu1);

 for (j=0; j<NN; j++) arl[j] = 1.;
 LU_solve(a,arl,NN);

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*NN+j] = w[j] * phi( (z[i]-log(1.+exp(z[j])))/adjust + k, mu0)/adjust;
   a[i*NN+N] = phi( (z[i]-log(1.+exp(zr)))/adjust + k, mu0)/adjust;
 }
 for (j=0; j<N; j++)
    a[N*NN+j] = w[j] * PHI( (zr-log(1.+exp(z[j])))/adjust + k, mu0);
 a[N*NN+N] = PHI( (zr-log(1.+exp(zr)))/adjust + k, mu0);

 pmethod(NN, a, &status, &rho, psi, &noofit);

 ad = psi[N]*arl[N];
 norm = psi[N];
 for (j=0; j<N; j++) {
   ad += w[j] * arl[j] * psi[j];
   norm += w[j] * psi[j];
 }
 ad /= norm;
 rho0 = rho;

 Free(a);
 Free(arl);
 Free(psi);
 Free(w);
 Free(z);

 return ad;
}


/* functions based on Srivastava & Wu (1997)
   Evaluation of optimum weights and average run lengths in EWMA control schemes,
   Commun. Stat., Theory Methods 26, 1253-1267.
   DOI: 10.1080/03610929708831980
*/


double xe2_SrWu_crit(double l, double L0)
{ double a, c;
 a = 2. * log( l * L0 * sqrt(2/PI) );
 c = sqrt( a - log(a-1.) ) + (1.-l)/2.;
 return c;
}


double xe2_SrWu_arl(double l, double c, double mu)
{ double g, w, arl=-1.;
 g = c * sqrt( l/2./mu/mu );
 w = c + 1.166*sqrt( mu * l ) - sqrt( 2.*mu*mu/l );
 if ( g < 1. ) arl = -log( 1.- g) /l - g/4./(1.-g)/mu/mu + .75;
 if ( g > 1. && fabs(mu) > 1. ) arl = PHI(w,0.)/phi(w,0.)/l/w;
 return arl;
}


double xe2_SrWu_arl_full(double l, double c, double mu)
{ double eta, Lmu, alpha1, alpha2, h1, h2, f1, f2, arl=-1., *w, *z;
  int i, qm=50;
  
 mu = fabs(mu); 

 w = vector(qm);
 z = vector(qm); 
 
 Lmu = c + 1.16*sqrt(l*mu);
 
 eta = mu * sqrt(2./l);
 
 gausslegendre(qm, 0, Lmu, z, w);
 
 alpha1 = 0.; alpha2 = 0.;
 for (i=0; i<qm; i++) {
   alpha1 += w[i] / phi(z[i]+eta, 0.);
   alpha2 += w[i] / phi(z[i]-eta, 0.);
 }
 
 h1 = alpha1 / (alpha1 + alpha2);
 h2 = alpha2 / (alpha1 + alpha2);
 
 f1 = 0.; f2 = 0.;
 for (i=0; i<qm; i++) {
   f1 += w[i] * ( PHI(z[i]+eta, 0.) - PHI( eta, 0.) ) / phi(z[i]+eta, 0.);
   f2 += w[i] * ( PHI(z[i]-eta, 0.) - PHI(-eta, 0.) ) / phi(z[i]-eta, 0.);
 }
 
 arl = ( h1*f2 + h2*f1 )/l;
 
 Free(w);
 Free(z);
 
 return arl;
}


double xe2_SrWu_lambda(double delta, double L0)
{ double dstar, b, l;
 dstar = 0.5117;
 b = 2.*log( 2.*sqrt(2./PI)*dstar*delta*delta*L0 );
 l = 2*dstar*delta*delta/( b - log(b) );
 return l;
}



double xe1_iglarl(double l, double c, double zr, double hs, double mu, int N)
{ double *a, *g, *w, *z, arl;
  int i, j, NN;

 NN = N + 1;
 a = matrix(NN,NN);
 g = vector(NN);
 w = vector(NN);
 z = vector(NN);

 c  *= sqrt( l/(2.-l) );
 zr *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );

 gausslegendre(N,zr,c,z,w);

 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) a[i*NN+j] = -w[j]/l * phi((z[j]-(1.-l)*z[i])/l,mu);
   ++a[i*NN+i];
   a[i*NN+N] = - PHI((zr-(1.-l)*z[i])/l,mu);
 }
 for (j=0;j<N;j++)
    a[N*NN+j] = -w[j]/l * phi((z[j]-(1.-l)*zr)/l,mu);
 a[N*NN+N] = 1. - PHI(zr,mu);

 for (j=0;j<NN;j++) g[j] = 1.;
 LU_solve(a,g,NN);

 arl = 1. + PHI((zr-(1.-l)*hs)/l,mu) * g[N];
 for (j=0;j<N;j++)
   arl += w[j]/l * phi((z[j]-(1.-l)*hs)/l,mu) * g[j];

 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}


double xe2_iglarl(double l, double c, double hs, double mu, int N)
{ double *a, *g, *w, *z, arl;
  int i, j;

 a = matrix(N,N);
 g = vector(N);
 w = vector(N);
 z = vector(N);

 c  *= sqrt( l/(2.-l) ); 
 hs *= sqrt( l/(2.-l) );

 gausslegendre(N,-c,c,z,w);

 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) a[j*N+i] = -w[j]/l * phi( (z[j]-(1.-l)*z[i])/l,mu);
   ++a[i*N+i];
 }

 for (j=0;j<N;j++) g[j] = 1.;
 solve(&N, a, g);

 arl = 1.;
 for (j=0;j<N;j++)
   arl += w[j]/l * phi( (z[j]-(1.-l)*hs)/l,mu) * g[j];

 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}


double xe2_iglarl_f(double l, double c, double mu, int N, double *g, double *w, double *z)
{ double *a;
  int i, j;

 a = matrix(N,N);

 c  *= sqrt( l/(2.-l) ); 

 gausslegendre(N, -c, c, z, w);

 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) a[j*N+i] = -w[j]/l * phi( (z[j]-(1.-l)*z[i])/l,mu);
   ++a[i*N+i];
 }

 for (j=0;j<N;j++) g[j] = 1.;
 solve(&N, a, g); 

 Free(a);

 return 0.;
}


double xte2_iglarl(double l, double c, double hs, int df, double mu, int N, int subst)
{ double *a, *g, *w, *z, arl, norm=1., arg=0., korr=1.;
  int i, j;

 a = matrix(N,N);
 g = vector(N);
 w = vector(N);
 z = vector(N);

 c  *= sqrt( l/(2.-l) ); 
 hs *= sqrt( l/(2.-l) );

 switch ( subst ) {
   case IDENTITY: gausslegendre(N, -c, c, z, w); norm = 1.; break;
   case SIN:   gausslegendre(N, -PI/2., PI/2., z, w); norm = 1.; break;
   case SINH:  gausslegendre(N, -1., 1., z, w); norm = sinh(1.); break;
   case TAN:   gausslegendre(N, -PI/4., PI/4., z, w); norm = 1.; break;
 }     
 
 c /= norm;

 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) {
     switch ( subst ) {
       case IDENTITY: arg = z[j] - (1.-l)*z[i]; korr = 1.; break;
       case SIN:   arg = c*sin(z[j]) - (1.-l)*c*sin(z[i]); korr = c*cos(z[j]); break;
       case SINH:  arg = c*sinh(z[j]) - (1.-l)*c*sinh(z[i]); korr = c*cosh(z[j]); break;
       case TAN:   arg = c*tan(z[j]) - (1.-l)*c*tan(z[i]); korr = c/( cos(z[j])*cos(z[j]) ); break;
     }
     a[i*N+j] = -w[j]/l * pdf_t( arg/l - mu, df) * korr;
   }
   ++a[i*N+i];
 }

 for (j=0;j<N;j++) g[j] = 1.;
 LU_solve(a,g,N);

 arl = 1.;
 for (j=0;j<N;j++) {
   switch ( subst ) {
     case IDENTITY: arg = z[j] - (1.-l)*hs; korr = 1.; break;
     case SIN:   arg = c*sin(z[j]) - (1.-l)*hs; korr = c*cos(z[j]); break;
     case SINH:  arg = c*sinh(z[j]) - (1.-l)*hs; korr = c*cosh(z[j]); break;
     case TAN:   arg = c*tan(z[j]) - (1.-l)*hs; korr = c/( cos(z[j])*cos(z[j]) ); break;
   }
   arl += w[j]/l * pdf_t( arg/l - mu, df) * g[j] * korr;
 }

 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}


double xte1_iglarl(double l, double c, double zr, double hs, int df, double mu, int N, int subst)
{ double *a, *g, *w, *z, arl, norm=1., arg=0., korr=1., z1, z2;
  int i, j, NN;

 NN = N + 1;
 a = matrix(NN,NN);
 g = vector(NN);
 w = vector(NN);
 z = vector(NN);

 c  *= sqrt( l/(2.-l) ); 
 hs *= sqrt( l/(2.-l) );
 zr *= sqrt( l/(2.-l) );
 
 /* IDENTITY as default */
 z1 = zr;
 z2 = c;

 switch ( subst ) {
   case IDENTITY: z1 = zr; z2 = c; break;
   case SIN:      if ( zr >= -c ) { z1 = asin(zr/c);  z2 = PI/2.;     norm=c; } else { z1 = -PI/2.;     z2 = asin(c/fabs(zr));  norm=fabs(zr); } break;
   case SINH:     if ( zr >= -c ) { z1 = asinh(zr/c); z2 = asinh(1.); norm=c; } else { z1 = asinh(-1.); z2 = asinh(c/fabs(zr)); norm=fabs(zr); } break;
   case TAN:      if ( zr >= -c ) { z1 = atan(zr/c);  z2 = PI/4.;     norm=c; } else { z1 = -PI/4.;     z2 = atan(c/fabs(zr));  norm=fabs(zr); } break;
 }     

 gausslegendre(N, z1, z2, z, w); 

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) {
     switch ( subst ) {
       case IDENTITY: arg = z[j] - (1.-l)*z[i];                        korr = 1.; break;
       case SIN:      arg = norm * ( sin(z[j]) - (1.-l)*sin(z[i]) );   korr = norm*cos(z[j]); break;
       case SINH:     arg = norm * ( sinh(z[j]) - (1.-l)*sinh(z[i]) ); korr = norm*cosh(z[j]); break;
       case TAN:      arg = norm * ( tan(z[j]) - (1.-l)*tan(z[i]) );   korr = norm/( cos(z[j])*cos(z[j]) ); break;
     }
     a[i*NN+j] = -w[j]/l * pdf_t( arg/l - mu, df) * korr;
   }
   ++a[i*NN+i];
   switch ( subst ) {
     case IDENTITY: arg = zr - (1.-l)*z[i];         break;
     case SIN:      arg = zr - (1.-l)*norm*sin(z[i]);  break;
     case SINH:     arg = zr - (1.-l)*norm*sinh(z[i]); break;
     case TAN:      arg = zr - (1.-l)*norm*tan(z[i]);  break;
   }
   a[i*NN+N] = - cdf_t( arg/l - mu, df); 
 }

 for (j=0; j<N; j++) {
   switch ( subst ) {
     case IDENTITY: arg = z[j] - (1.-l)*zr;            korr = 1.; break;
     case SIN:      arg = norm*sin(z[j]) - (1.-l)*zr;  korr = norm*cos(z[j]); break;
     case SINH:     arg = norm*sinh(z[j]) - (1.-l)*zr; korr = norm*cosh(z[j]); break;
     case TAN:      arg = norm*tan(z[j]) - (1.-l)*zr;  korr = norm/( cos(z[j])*cos(z[j]) ); break;
   }
   a[N*NN+j] = -w[j]/l * pdf_t( arg/l - mu, df) * korr;
 }

 a[N*NN+N] = 1. - cdf_t( zr - mu, df); 

 for (j=0; j<NN; j++) g[j] = 1.;
 LU_solve(a, g, NN);

 arl = 1. + cdf_t( (zr-(1.-l)*hs)/l - mu, df) * g[N];
 for (j=0;j<N;j++) {
   switch ( subst ) {
     case IDENTITY: arg = z[j] - (1.-l)*hs;            korr = 1.; break;
     case SIN:      arg = norm*sin(z[j]) - (1.-l)*hs;  korr = norm*cos(z[j]); break;
     case SINH:     arg = norm*sinh(z[j]) - (1.-l)*hs; korr = norm*cosh(z[j]); break;
     case TAN:      arg = norm*tan(z[j]) - (1.-l)*hs;  korr = norm/( cos(z[j])*cos(z[j]) ); break;
   }
   arl += w[j]/l * pdf_t( arg/l - mu, df) * g[j] * korr;
 }

 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}


double xe1_iglarl_drift(double l, double c, double zr, double hs, double delta, int m, int N, int with0)
{ double *a, *g, *w, *z, arl, *MUs, *ARLs;
  int i, j, NN, m_;

 NN = N + 1;
 a = matrix(NN, NN);
 g = vector(NN);
 w = vector(NN);
 z = vector(NN);
 ARLs = vector(NN);
 MUs  = vector(m+1); 

 c  *= sqrt( l/(2.-l) );
 zr *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );

 gausslegendre(N, zr, c, z, w);

 if ( with0 ) {
   for (i=0;i<=m;i++) MUs[i] = (double)i * delta;
 } else {
   for (i=0;i<=m;i++) MUs[i] = (double)(i+1.) * delta;
 }

 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) a[i*NN+j] = -w[j]/l * phi( (z[j]-(1.-l)*z[i])/l, MUs[m]);
   ++a[i*NN+i];
   a[i*NN+N] = - PHI( (zr-(1.-l)*z[i])/l, MUs[m]);
 }
 for (j=0;j<N;j++)
    a[N*NN+j] = -w[j]/l * phi( (z[j]-(1.-l)*zr)/l, MUs[m]);
 a[N*NN+N] = 1. - PHI(zr, MUs[m]);

 for (j=0;j<NN;j++) g[j] = 1.;
 LU_solve(a, g, N);

 for (m_=0;m_<m;m_++) {
   for (i=0;i<N;i++) {
     ARLs[i] = 1. + PHI(zr, MUs[m-m_]) * g[N];
     for (j=0;j<=N;j++) { 
       ARLs[i] += w[j]/l * phi( (z[j]-(1.-l)*z[i])/l, MUs[m-m_]) * g[j];
     }
   }
   for (j=0;j<=N;j++) g[j] = ARLs[j];
 }

 arl = 1. + PHI( (zr-(1.-l)*hs)/l, MUs[0]) * ARLs[N];
 for (j=0;j<N;j++) arl += w[j]/l * phi( (z[j]-(1.-l)*hs)/l, MUs[0]) * ARLs[j];

 Free(a);
 Free(g);
 Free(w);
 Free(z);
 Free(ARLs);
 Free(MUs);

 return arl;
}


double xe1_iglarl_drift_wo_m(double l, double c, double zr, double hs, double delta, int *m, int N, int with0)
{ int m_;
  double arl1, arl2, eps=1e-6;
 m_ = 4;
 arl1 = xe1_iglarl_drift(l, c, zr, hs, delta, m_, N, with0);
 arl2 = arl1 + 2.*eps;
 while ( fabs(arl2-arl1)>eps && (double)m_<1e4 ) {
   m_ = (int)round(1.5 * m_);
   arl1 = xe1_iglarl_drift(l, c, zr, hs, delta, m_, N, with0);
   arl2 = xe1_iglarl_drift(l, c, zr, hs, delta, m_+1, N, with0);
 }
 *m = m_;
 return arl1;
}


double xe1_iglarlm_drift(double l, double c, double zr, double hs, int q, double delta, int N, int nmax, int with0)
{ double *p0, *fn, *w, *z, arl0, rho, MEAN=0.,
         arl_minus=0., arl, arl_plus=0., mn_minus=0., mn_plus=0., nn, ratio=0.;
  int i, j, n, NN;

 NN = N + 1;
 w   = vector(NN);
 z   = vector(NN);
 fn  = matrix(nmax, NN);
 p0  = vector(nmax);

 c  *= sqrt( l/(2.-l) );
 zr *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );

 gausslegendre(N, zr, c, z, w);

 /* in-control, i. e. n<=q-1 */
 MEAN = 0.;

 for (n=1;n<=q-1;n++) {
  nn = (double) n;

  /* determine f_n, n=1,2,...,q-1 */
  if (n==1) {
    for (i=0;i<N;i++) fn[0*NN+i] = phi( (z[i]-(1.-l)*hs)/l, MEAN)/l;
    fn[0*NN+N] = PHI( (zr-(1.-l)*hs)/l, MEAN);
  } else {
    for (i=0;i<N;i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * phi( (z[i]-(1.-l)*zr)/l, MEAN)/l;
      for (j=0;j<N;j++) {
        fn[(n-1)*NN+i] += w[j] * fn[(n-2)*NN+j] * phi( (z[i]-(1.-l)*z[j])/l, MEAN)/l;
      }
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI( zr, MEAN);
    for (j=0;j<N;j++) fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI( (zr-(1.-l)*z[j])/l, MEAN);
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];
 }

 /* out-of-control, i.e. t>=q */
 arl0 = 1.; rho = 0.;

 for (n=q;n<=nmax;n++) {
  nn = (double) n;

  /* determine f_n, n=q,q+1,... */
  if ( with0 ) {
    MEAN = (nn-(double)q) * delta;
  } else {
    MEAN = (nn-(double)q+1.) * delta;
  }

  if (n==1) {
    for (i=0;i<N;i++) fn[0*NN+i] = phi( (z[i]-(1.-l)*hs)/l, MEAN)/l;
    fn[0*NN+N] = PHI( (zr-(1.-l)*hs)/l, MEAN);
  } else {
    for (i=0;i<N;i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * phi( (z[i]-(1.-l)*zr)/l, MEAN)/l;
      for (j=0;j<N;j++) {
        fn[(n-1)*NN+i] += w[j] * fn[(n-2)*NN+j] * phi( (z[i]-(1.-l)*z[j])/l, MEAN)/l;
      }
      if (n==q && q>1) fn[(n-1)*NN+i] /= p0[q-2];
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI( zr, MEAN);
    for (j=0;j<N;j++) fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI( (zr-(1.-l)*z[j])/l, MEAN);
    if (n==q && q>1) fn[(n-1)*NN+N] /= p0[q-2];
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];

  /* computation of m_n+1^- and m_n+1^+, n=m-1,m,... */
  mn_minus = 1.; mn_plus = 0.;
  if (n>q) {
   for (i=0;i<NN;i++) {
    if (fn[(n-2)*NN+i]==0)
     if (fn[(n-1)*NN+i]==0) ratio = 0.; else ratio = 1.;
    else ratio = fn[(n-1)*NN+i]/fn[(n-2)*NN+i];
    if ( ratio<mn_minus ) mn_minus = ratio;
    if ( ratio>mn_plus ) mn_plus = ratio;
   }
  }

  if (n>q) rho = p0[n-1]/p0[n-2];

  /* computation of ARL, ARL^-, and ARL^+ */
  arl = arl0 + p0[n-1]/(1.-rho);
  if (mn_minus<1.) arl_minus = arl0 + p0[n-1]/(1.-mn_minus);
  else             arl_minus = -1.;
  if (mn_plus<1.) arl_plus = arl0 + p0[n-1]/(1.-mn_plus);
  else            arl_plus = -1.;
  arl0 += p0[n-1];

  if ( fabs((arl_plus-arl_minus)) < 1e-7 ) n = nmax+1;
 }

 arl = (arl_plus+arl_minus)/2.; rho0 = rho;

 Free(w);
 Free(z);
 Free(fn);
 Free(p0);

 return arl;
}


double xe2_iglarl_drift(double l, double c, double hs, double delta, int m, int N, int with0)
{ double *a, *g, *w, *z, arl, *MUs, *ARLs;
  int i, j, m_;

 a = matrix(N,N);
 g = vector(N);
 w = vector(N);
 z = vector(N);
 ARLs = vector(N);
 MUs  = vector(m+1); 

 c  *= sqrt( l/(2.-l) ); 
 hs *= sqrt( l/(2.-l) );

 gausslegendre(N, -c, c, z, w);

 if ( with0 ) {
   for (i=0;i<=m;i++) MUs[i] = (double)i * delta;
 } else {
   for (i=0;i<=m;i++) MUs[i] = (double)(i+1.) * delta;
 }

 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) a[i*N+j] = -w[j]/l * phi( (z[j]-(1.-l)*z[i])/l, MUs[m]);
   ++a[i*N+i];
 }
 for (j=0;j<N;j++) g[j] = 1.;
 LU_solve(a, g, N);

 for (m_=0;m_<m;m_++) {
   for (i=0;i<N;i++) {
     ARLs[i] = 1.;
     for (j=0;j<N;j++) {
       ARLs[i] += w[j]/l * phi( (z[j]-(1.-l)*z[i])/l, MUs[m-m_]) * g[j];
     }
   }
   for (j=0;j<N;j++) g[j] = ARLs[j];
 }

 arl = 1.;
 for (j=0;j<N;j++) arl += w[j]/l * phi( (z[j]-(1.-l)*hs)/l, MUs[0]) * ARLs[j];

 Free(a);
 Free(g);
 Free(w);
 Free(z);
 Free(ARLs);
 Free(MUs);

 return arl;
}


double xe2_iglarl_drift_wo_m(double l, double c, double hs, double delta, int *m, int N, int with0)
{ int m_;
  double arl1, arl2, eps=1e-6;
 m_ = 4;
 arl1 = xe2_iglarl_drift(l, c, hs, delta, m_, N, with0);
 arl2 = arl1 + 2.*eps;
 while ( fabs(arl2-arl1)>eps && (double)m_<1e4 ) {
   m_ = (int)round(1.5 * m_);
   arl1 = xe2_iglarl_drift(l, c, hs, delta, m_, N, with0);
   arl2 = xe2_iglarl_drift(l, c, hs, delta, m_+1, N, with0);
 }
 *m = m_;
 return arl1;
}


double xe2_iglarlm_drift(double l, double c, double hs, int q, double delta, int N, int nmax, int with0)
{ double *p0, *fn, *w, *z, arl0, rho, MEAN=0.,
         arl_minus=0., arl, arl_plus=0., mn_minus=0., mn_plus=0., nn, ratio=0.;
  int i, j, n;

 w   = vector(N);
 z   = vector(N);
 fn  = matrix(nmax, N);
 p0  = vector(nmax);

 c  *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );

 gausslegendre(N, -c, c, z, w);

 /* in-control, i. e. n<=q-1 */
 MEAN = 0.;

 for (n=1;n<=q-1;n++) {
  nn = (double) n;

  /* determine f_n, n=1,2,...,q-1 */
  if (n==1) {
    for (i=0;i<N;i++) fn[0*N+i] = phi( (z[i]-(1.-l)*hs)/l, MEAN)/l;
  } else {
    for (i=0;i<N;i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0;j<N;j++) {
        fn[(n-1)*N+i] += w[j] * fn[(n-2)*N+j] * phi( (z[i]-(1.-l)*z[j])/l, MEAN)/l;
      }
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*N+i];
 }

 /* out-of-control, i.e. t>=q */
 arl0 = 1.; rho = 0.;
 for (n=q;n<=nmax;n++) {
  nn = (double) n;

  /* determine f_n, n=q,q+1,... */
  if ( with0 ) {
    MEAN = (nn-(double)q) * delta;
  } else {
    MEAN = (nn-(double)q+1.) * delta;
  }

  if (n==1) {
    for (i=0;i<N;i++) fn[0*N+i] = phi( (z[i]-(1.-l)*hs)/l, MEAN)/l;
  } else {
    for (i=0;i<N;i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0;j<N;j++) {
        fn[(n-1)*N+i] += w[j] * fn[(n-2)*N+j] * phi( (z[i]-(1.-l)*z[j])/l, MEAN)/l;
      }
      if (n==q && q>1) fn[(n-1)*N+i] /= p0[q-2];
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*N+i];

  /* computation of m_n+1^- and m_n+1^+, n=m-1,m,... */
  mn_minus = 1.; mn_plus = 0.;
  if (n>q) {
   for (i=0;i<N;i++) {
    if (fn[(n-2)*N+i]==0)
     if (fn[(n-1)*N+i]==0) ratio = 0.; else ratio = 1.;
    else ratio = fn[(n-1)*N+i]/fn[(n-2)*N+i];
    if ( ratio<mn_minus ) mn_minus = ratio;
    if ( ratio>mn_plus ) mn_plus = ratio;
   }
  }

  if (n>q) rho = p0[n-1]/p0[n-2];

  /* computation of ARL, ARL^-, and ARL^+ */
  arl = arl0 + p0[n-1]/(1.-rho);
  if (mn_minus<1.) arl_minus = arl0 + p0[n-1]/(1.-mn_minus);
  else             arl_minus = -2.;
  if (mn_plus<1.) arl_plus = arl0 + p0[n-1]/(1.-mn_plus);
  else            arl_plus = -1.;
  arl0 += p0[n-1];

  if ( fabs((arl_plus-arl_minus)) < 1e-7 ) n = nmax+1;
 }

 arl = (arl_plus+arl_minus)/2.; rho0 = rho;

 Free(w);
 Free(z);
 Free(fn);
 Free(p0);

 return arl;
}


double xe2_Warl_drift(double l, double c, double hs, double delta, int N, int nmax, int with0)
{ double *Pn, *w, *z, *p0, MEAN, nn, ratio, arl_minus=0., arl0=1., arl_plus=0., mn_minus=1., mn_plus=0.;
  int i, j, n;

 c  *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );

 w  = vector(N);
 z  = vector(N);
 Pn = matrix(nmax,N);
 p0 = vector(nmax);

 gausslegendre(N,-c,c,z,w);

 arl0 = 1.;

 for (n=1;n<=nmax;n++) {
   nn = (double)n;
   if ( with0 ) {
     MEAN = (nn-1.) * delta;
   } else {
     MEAN = nn * delta;
   }

   if (n==1)
     for (i=0;i<N;i++)
       Pn[i] = PHI( (c-(1.-l)*z[i])/l, MEAN) - PHI( (-c-(1.-l)*z[i])/l, MEAN);
   else
     for (i=0;i<N;i++) {
       Pn[(n-1)*N+i] = 0.;
       for (j=0;j<N;j++)
         Pn[(n-1)*N+i] += w[j]/l*phi( (z[j]-(1.-l)*z[i])/l, MEAN)*Pn[(n-2)*N+j];
     }
   
   if (n==1)
     p0[0] = PHI( (c-(1.-l)*hs)/l, MEAN) - PHI( (-c-(1.-l)*hs)/l, MEAN);
   else {
     p0[n-1] = 0.;
     for (j=0;j<N;j++) p0[n-1] += w[j]/l * phi( (z[j]-(1.-l)*hs)/l, MEAN) * Pn[(n-2)*N+j];
   }

   mn_minus = 1.; mn_plus = 0.;
   if (n>1) {
     for (i=0;i<N;i++) {
       if (Pn[(n-2)*N+i]==0)
         if (Pn[(n-1)*N+i]==0) ratio = 0.;
         else ratio = 1.;
       else ratio = Pn[(n-1)*N+i]/Pn[(n-2)*N+i];
      if ( ratio<mn_minus ) mn_minus = ratio;
      if ( ratio>mn_plus ) mn_plus = ratio;
     }
   }

   if (0.<mn_minus && mn_minus<1.) arl_minus = arl0 + p0[n-1]/(1.-mn_minus);
   else             arl_minus = -2.;
   if (0.<mn_plus && mn_plus<1.) arl_plus = arl0 + p0[n-1]/(1.-mn_plus);
   else            arl_plus = -1.;
   arl0 += p0[n-1];

   if ( fabs( (arl_plus-arl_minus)/arl_minus )<FINALeps ) n = nmax+1;
 }

 Free(p0);
 Free(Pn);
 Free(z);
 Free(w);

 return (arl_plus+arl_minus)/2.;
}


double xe2_Warl(double l, double c, double hs, double mu, int N, int nmax)
{ double *Sm, *Pn, *w, *z, *p0, ratio, arl_minus=0., arl=1., arl_plus=0., mn_minus=1., mn_plus=0.;
  int i, j, n;

 c  *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );

 Sm = matrix(N,N);
 w  = vector(N);
 z  = vector(N);
 Pn = matrix(nmax,N);
 p0 = vector(nmax);

 gausslegendre(N,-c,c,z,w);

 for (i=0;i<N;i++)
   for (j=0;j<N;j++)
     Sm[i*N+j] = w[j]/l * phi( (z[j]-(1.-l)*z[i])/l, mu);

 arl = 1.;

 for (n=1;n<=nmax;n++) {

   if (n==1)
     for (i=0;i<N;i++)
       Pn[i] = PHI( (c-(1.-l)*z[i])/l, mu) - PHI( (-c-(1.-l)*z[i])/l, mu);
   else
     for (i=0;i<N;i++) {
       Pn[(n-1)*N+i] = 0.;
       for (j=0;j<N;j++)
         Pn[(n-1)*N+i] += Sm[i*N+j] * Pn[(n-2)*N+j];
     }

   if (n==1)
     p0[0] = PHI( (c-(1.-l)*hs)/l, mu) - PHI( (-c-(1.-l)*hs)/l, mu);
   else {
     p0[n-1] = 0.;
     for (j=0;j<N;j++) p0[n-1] += w[j]/l * phi( (z[j]-(1.-l)*hs)/l, mu) * Pn[(n-2)*N+j];
   }

   mn_minus = 1.; mn_plus = 0.;
   if (n>1) {
     for (i=0;i<N;i++) {
       if (Pn[(n-2)*N+i]==0)
         if (Pn[(n-1)*N+i]==0) ratio = 0.;
         else ratio = 1.;
       else ratio = Pn[(n-1)*N+i]/Pn[(n-2)*N+i];
      if ( ratio<mn_minus ) mn_minus = ratio;
      if ( ratio>mn_plus ) mn_plus = ratio;
     }

     arl_minus = arl + p0[n-1]/(1.-mn_minus);
     arl_plus = arl + p0[n-1]/(1.-mn_plus);
   }
   arl += p0[n-1];

   if ( fabs( (arl_plus-arl_minus)/arl_minus )<FINALeps ) n = nmax+1;
 }

 Free(p0);
 Free(Pn);
 Free(z);
 Free(w);
 Free(Sm);

 return (arl_plus+arl_minus)/2.;
}


double xe2_Wq(double l, double c, double p, double hs, double mu, int N, int nmax)
{ double *Sm, *Pn, *w, *z, *p0, ratio, q_minus=0., q_plus=0., mn_minus=1., mn_plus=0., enumerator=0., Wq=0.;
  int i, j, n;

 c  *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );

 Sm = matrix(N, N);
 w  = vector(N);
 z  = vector(N);
 Pn = matrix(nmax, N);
 p0 = vector(nmax);
 gausslegendre(N, -c, c, z, w);

 for (i=0; i<N; i++)
   for (j=0; j<N; j++) Sm[i*N+j] = w[j]/l * phi( (z[j]-(1.-l)*z[i])/l, mu);

 for (n=1; n<=nmax; n++) {

   if ( n==1 )
     for (i=0; i<N; i++) Pn[i] = PHI( (c-(1.-l)*z[i])/l, mu) - PHI( (-c-(1.-l)*z[i])/l, mu);
   else
     for (i=0; i<N; i++) {
       Pn[(n-1)*N+i] = 0.;
       for (j=0; j<N; j++) Pn[(n-1)*N+i] += Sm[i*N+j] * Pn[(n-2)*N+j];
     }

   if ( n==1 )
     p0[0] = PHI( (c-(1.-l)*hs)/l, mu) - PHI( (-c-(1.-l)*hs)/l, mu);
   else {
     p0[n-1] = 0.;
     for (j=0; j<N; j++) p0[n-1] += w[j]/l * phi( (z[j]-(1.-l)*hs)/l, mu) * Pn[(n-2)*N+j];
   }

   if ( p0[n-1] < 1.-p ) {
     Wq = (double)n;
     n = nmax+1;
   } else {     
     mn_minus = 1.; mn_plus = 0.;
     if ( n>1 ) {
       for (i=0; i<N; i++) {
         if (Pn[(n-2)*N+i]==0)
           if (Pn[(n-1)*N+i]==0) ratio = 0.;
           else ratio = 1.;
         else ratio = Pn[(n-1)*N+i]/Pn[(n-2)*N+i];
        if ( ratio<mn_minus ) mn_minus = ratio;
        if ( ratio>mn_plus ) mn_plus = ratio;
       }
       enumerator = log( (1.-p)/p0[n-1] );
       q_minus = (double)n + enumerator/log(mn_minus);
       q_plus  = (double)n + enumerator/log(mn_plus);
       if ( fabs( ceil(q_plus) - ceil(q_minus) ) < .5 ) {
	 Wq = ceil(q_plus);
	 n = nmax +1;
       }
     }  /* n > 1 */ 
   } /* p0[n-1] >= 1.-p */
 } /* n=1; n<=nmax; n++ */

 Free(p0);
 Free(Pn);
 Free(z);
 Free(w);
 Free(Sm);

 return Wq;
}


double xte2_Wq(double l, double c, double p, double hs, int df, double mu, int N, int nmax, int subst)
{ double *Sm, *Pn, *w, *z, *p0, ratio, q_minus=0., q_plus=0., mn_minus=1., mn_plus=0., enumerator=0., Wq=0., norm=1., arg=0., korr=1.;
  int i, j, n;
 
 c  *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );

 Sm = matrix(N, N);
 w  = vector(N);
 z  = vector(N);
 Pn = matrix(nmax, N);
 p0 = vector(nmax);
 
 switch ( subst ) {
   case IDENTITY: gausslegendre(N, -c, c, z, w);         norm = 1.; break;
   case SIN:   gausslegendre(N, -PI/2., PI/2., z, w); norm = 1.; break;
   case SINH:  gausslegendre(N, -1., 1., z, w);       norm = sinh(1.); break;
   case TAN:   gausslegendre(N, -PI/4., PI/4., z, w); norm = 1.; break;
 }     
 
 c /= norm;

 for (i=0; i<N; i++)
   for (j=0; j<N; j++) {
     switch ( subst ) {
       case IDENTITY: arg = z[j] - (1.-l)*z[i];                 korr = 1.; break;
       case SIN:   arg = c*sin(z[j]) - (1.-l)*c*sin(z[i]);   korr = c*cos(z[j]); break;
       case SINH:  arg = c*sinh(z[j]) - (1.-l)*c*sinh(z[i]); korr = c*cosh(z[j]); break;
       case TAN:   arg = c*tan(z[j]) - (1.-l)*c*tan(z[i]);   korr = c/( cos(z[j])*cos(z[j]) ); break;
     }
     Sm[i*N+j] = w[j]/l * pdf_t( arg/l - mu, df) * korr;
   }

/* for (n=1; n<=nmax; n++) {*/
 for (n=1; n<=100; n++) {

   if ( n==1 )
     for (i=0; i<N; i++) {
       switch ( subst ) {
         case IDENTITY: arg = z[i]; break;
         case SIN:   arg = c*sin(z[i]); break;
         case SINH:  arg = c*sinh(z[i]); break;
         case TAN:   arg = c*tan(z[i]); break;
       }
       Pn[i] = cdf_t( ( c*norm - (1.-l)*arg )/l - mu, df) - cdf_t( ( -c*norm - (1.-l)*arg )/l - mu, df);
     }
   else
     for (i=0; i<N; i++) {
       Pn[(n-1)*N+i] = 0.;
       for (j=0; j<N; j++) Pn[(n-1)*N+i] += Sm[i*N+j] * Pn[(n-2)*N+j];
     }

   if ( n==1 )
     p0[0] = cdf_t( ( c*norm - (1.-l)*hs )/l - mu, df) - cdf_t( ( -c*norm - (1.-l)*hs )/l - mu, df);
   else {
     p0[n-1] = 0.;
     for (j=0; j<N; j++) {
       switch ( subst ) {
         case IDENTITY: arg = z[j];         korr = 1.; break;
         case SIN:   arg = c*sin(z[j]);  korr = c*cos(z[j]); break;
         case SINH:  arg = c*sinh(z[j]); korr = c*cosh(z[j]); break;
         case TAN:   arg = c*tan(z[j]);  korr = c/( cos(z[j])*cos(z[j]) ); break;
       }
       p0[n-1] += w[j]/l * pdf_t( ( arg - (1.-l)*hs )/l - mu, df) * Pn[(n-2)*N+j] * korr;
     }
   }

   /*printf("%4d\t\t%.6f\n", n, p0[n-1]);*/
   
   if ( p0[n-1] < 1.-p ) {
     Wq = (double)n;
     n = nmax+1;
   } else {     
     mn_minus = 1.; mn_plus = 0.;
     if ( n>1 ) {
       for (i=0; i<N; i++) {
         if (Pn[(n-2)*N+i]==0)
           if (Pn[(n-1)*N+i]==0) ratio = 0.;
           else ratio = 1.;
         else ratio = Pn[(n-1)*N+i]/Pn[(n-2)*N+i];
        if ( ratio<mn_minus ) mn_minus = ratio;
        if ( ratio>mn_plus ) mn_plus = ratio;
       }
       enumerator = log( (1.-p)/p0[n-1] );
       q_minus = (double)n + enumerator/log(mn_minus);
       q_plus  = (double)n + enumerator/log(mn_plus);
       if ( fabs( ceil(q_plus) - ceil(q_minus) ) < .5 ) {
	 Wq = ceil(q_plus);
	 n = nmax +1;
       }
     }  /* n > 1 */ 
   } /* p0[n-1] >= 1.-p */
 } /* n=1; n<=nmax; n++ */

 Free(p0);
 Free(Pn);
 Free(z);
 Free(w);
 Free(Sm);

 return Wq;
}


double xe2_Wqm(double l, double c, double p, double hs, int q, double mu0, double mu1, int mode, int N, int nmax)
{ double *Smatrix, *p0, *fn, *w, *z, dn, rn, cn, rn0, cn0, delta=0.,
         q_minus=2., q_plus=3., mn_minus, mn_plus, nn, fSt, aSt, ratio, enumerator=0., nq, Wq=0.;
  int i, j, n;

 fSt = 0.5;
 aSt = ( -2./log10(1.-fSt) - 1.)/19.;

 c  *= sqrt( l/(2.-l) ); 
 hs *= sqrt( l/(2.-l) );
 if ( mode==fir || mode==both ) delta = 2.*hs;

 Smatrix = matrix(N, N);
 w       = vector(N);
 z       = vector(N);
 fn      = matrix(nmax, N);
 p0      = vector(nmax);

 gausslegendre(N, -c, c, z, w);

 rn = 1.; cn = 0.; rn0 = 1., cn0 = 0.;

 /* in-control, i. e. n<=q-1 */
 for (n=1; n<=q-1; n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=1,2,...,q-1 */
  switch ( mode ) {
    case vacl: rn = sqrt( 1. - pow(1.-l,2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }
 
  /* determine f_n, n=1,2,...,q-1 */
  if ( n==1 ) {
    for (i=0; i<N; i++)
      if ( mode==stat )
        fn[0*N+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)),mu0);
      else
        fn[0*N+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l,mu0);
  } 
  else {
    for (i=0; i<N; i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0; j<N; j++) {
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l*phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l,mu0);
      }
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0; i<N; i++) p0[n-1] += w[i] * fn[(n-1)*N+i];

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;
 }

 /* out-of-control, i.e. t>=q */
 for (n=q; n<=nmax; n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=q,q+1,... */
  switch ( mode ) {
    case vacl: rn = sqrt( 1. - pow(1.-l,2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }

  /* determine f_n, n=q,q+1,... */
  if ( n==1 ) {
    for (i=0; i<N; i++)
      if ( mode==stat )
        fn[0*N+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)),mu1);
      else
        fn[0*N+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l,mu1);
  }
  else {
    for (i=0; i<N; i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0; j<N; j++) 
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l*phi( (cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l,mu1);
      if ( n==q && q>1 ) fn[(n-1)*N+i] /= p0[q-2];
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0; i<N; i++) p0[n-1] += w[i] * fn[(n-1)*N+i];
  nq = (double)(n-q+1);

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;

  if ( p0[n-1] < 1.-p ) {
    Wq = nq;
    n = nmax+1;
  } else {     
    /* computation of m_n+1^- and m_n+1^+, n=m-1,m,... */
    mn_minus = 1.; mn_plus = 0.;
    if ( n > q ) {
      for (i=0; i<N; i++) {
        if (fn[(n-2)*N+i]==0)
          if (fn[(n-1)*N+i]==0) ratio = 0.; else ratio = 1.;
        else ratio = fn[(n-1)*N+i]/fn[(n-2)*N+i];
        if ( ratio<mn_minus ) mn_minus = ratio;
        if ( ratio>mn_plus ) mn_plus = ratio;
      }
      enumerator = log( (1.-p)/p0[n-1] );         
      q_minus = nq + enumerator/log(mn_minus);
      q_plus  = nq + enumerator/log(mn_plus);
      if ( fabs( ceil(q_plus) - ceil(q_minus) ) < .5 ) {
	 Wq = ceil(q_plus);
	 n = nmax +1;
      }   
    } /* n > q */  
  } /* p0[n-1] >= 1.-p */
 } /* n=q; n<=nmax; n++ */ 

 Free(Smatrix);
 Free(w);
 Free(z);
 Free(fn);
 Free(p0);

 return Wq;
}


double xte2_Wqm(double l, double c, double p, double hs, int df, int q, double mu0, double mu1, int mode, int N, int nmax, int subst)
{ double *Smatrix, *p0, *fn, *w, *z, dn, rn, cn, rn0, cn0, delta=0.,
         q_minus=2., q_plus=3., mn_minus, mn_plus, nn, fSt, aSt, ratio, enumerator=0., nq, Wq=0., norm=1., arg=0., korr=1.;
  int i, j, n;

 fSt = 0.5;
 aSt = ( -2./log10(1.-fSt) - 1.)/19.;

 c  *= sqrt( l/(2.-l) ); 
 hs *= sqrt( l/(2.-l) );
 if ( mode==fir || mode==both ) delta = 2.*hs;

 Smatrix = matrix(N, N);
 w       = vector(N);
 z       = vector(N);
 fn      = matrix(nmax, N);
 p0      = vector(nmax);

 switch ( subst ) {
   case IDENTITY: gausslegendre(N, -c, c, z, w);         norm = 1.; break;
   case SIN:   gausslegendre(N, -PI/2., PI/2., z, w); norm = 1.; break;
   case SINH:  gausslegendre(N, -1., 1., z, w);       norm = sinh(1.); break;
   case TAN:   gausslegendre(N, -PI/4., PI/4., z, w); norm = 1.; break;
 }     
 
 c /= norm;

 rn = 1.; cn = 0.; rn0 = 1., cn0 = 0.;

 /* in-control, i. e. n<=q-1 */
 for (n=1; n<=q-1; n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=1,2,...,q-1 */
  switch ( mode ) {
    case vacl: rn = sqrt( 1. - pow(1.-l,2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c*norm);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c*norm);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }
 
  /* determine f_n, n=1,2,...,q-1 */
  if ( n==1 ) {
    for (i=0; i<N; i++) {
      switch ( subst ) {
        case IDENTITY: arg = z[i]; break;
        case SIN:   arg = c*sin(z[i]); break;
        case SINH:  arg = c*sinh(z[i]); break;
        case TAN:   arg = c*tan(z[i]); break;
      }
      if ( mode==stat )
        fn[0*N+i] = 1./sqrt(l/(2.-l)) * pdf_t( ( cn+rn*arg )/sqrt(l/(2.-l)) - mu0, df);
      else
        fn[0*N+i] = rn/l * pdf_t( ( cn+rn*arg - (1.-l)*hs )/l - mu0, df);
    }
  } 
  else {
    for (i=0; i<N; i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0; j<N; j++) {
	switch ( subst ) {
          case IDENTITY: arg = cn+rn*z[i] - (1.-l)*(cn0+rn0*z[j]);                 korr = 1.; break;
          case SIN:   arg = cn+rn*c*sin(z[i]) - (1.-l)*(cn0+rn0*c*sin(z[j]));   korr = c*cos(z[j]); break;
          case SINH:  arg = cn+rn*c*sinh(z[i]) - (1.-l)*(cn0+rn0*c*sinh(z[j])); korr = c*cosh(z[j]); break;
          case TAN:   arg = cn+rn*c*tan(z[i]) - (1.-l)*(cn0+rn0*c*tan(z[j]));   korr = c/( cos(z[j])*cos(z[j]) ); break;
        }
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l * pdf_t( arg/l - mu0, df) * korr;
      }
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0; i<N; i++) {
    switch ( subst ) {
      case IDENTITY: korr = 1.; break;
      case SIN:   korr = c*cos(z[i]); break;
      case SINH:  korr = c*cosh(z[i]); break;
      case TAN:   korr = c/( cos(z[i])*cos(z[i]) ); break;
    }
    p0[n-1] += w[i] * fn[(n-1)*N+i] * korr;
  }

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;
 }

 /* out-of-control, i.e. t>=q */
 for (n=q; n<=nmax; n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=q,q+1,... */
  switch ( mode ) {
    case vacl: rn = sqrt( 1. - pow(1.-l,2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c*norm);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c*norm);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }

  /* determine f_n, n=q,q+1,... */
  if ( n==1 ) {
    for (i=0; i<N; i++) {
      switch ( subst ) {
        case IDENTITY: arg = z[i]; break;
        case SIN:   arg = c*sin(z[i]); break;
        case SINH:  arg = c*sinh(z[i]); break;
        case TAN:   arg = c*tan(z[i]); break;
      }
      if ( mode==stat )
        fn[0*N+i] = 1./sqrt(l/(2.-l)) * pdf_t( ( cn+rn*arg )/sqrt(l/(2.-l)) - mu1, df);
      else
        fn[0*N+i] = rn/l * pdf_t( ( cn+rn*arg - (1.-l)*hs )/l - mu1, df);
    }
  }
  else {
    for (i=0; i<N; i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0; j<N; j++) {
	switch ( subst ) {
          case IDENTITY: arg = cn+rn*z[i] - (1.-l)*(cn0+rn0*z[j]);                 korr = 1.; break;
          case SIN:   arg = cn+rn*c*sin(z[i]) - (1.-l)*(cn0+rn0*c*sin(z[j]));   korr = c*cos(z[j]); break;
          case SINH:  arg = cn+rn*c*sinh(z[i]) - (1.-l)*(cn0+rn0*c*sinh(z[j])); korr = c*cosh(z[j]); break;
          case TAN:   arg = cn+rn*c*tan(z[i]) - (1.-l)*(cn0+rn0*c*tan(z[j]));   korr = c/( cos(z[j])*cos(z[j]) ); break;
        }
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l * pdf_t( arg/l - mu1, df) * korr;
      }
      if ( n==q && q>1 ) fn[(n-1)*N+i] /= p0[q-2];
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0;i<N;i++) {
    switch ( subst ) {
      case IDENTITY: korr = 1.; break;
      case SIN:   korr = c*cos(z[i]); break;
      case SINH:  korr = c*cosh(z[i]); break;
      case TAN:   korr = c/( cos(z[i])*cos(z[i]) ); break;
    }
    p0[n-1] += w[i] * fn[(n-1)*N+i] * korr;
  }
  nq = (double)(n-q+1.);

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;

  if ( p0[n-1] < 1.-p ) {
    Wq = nq;
    n = nmax+1;
  } else {     
    /* computation of m_n+1^- and m_n+1^+, n=m-1,m,... */
    mn_minus = 1.; mn_plus = 0.;
    if ( n > q ) {
      for (i=0; i<N; i++) {
        if (fn[(n-2)*N+i]==0)
          if (fn[(n-1)*N+i]==0) ratio = 0.; else ratio = 1.;
        else ratio = fn[(n-1)*N+i]/fn[(n-2)*N+i];
        if ( ratio<mn_minus ) mn_minus = ratio;
        if ( ratio>mn_plus ) mn_plus = ratio;
      }
      enumerator = log( (1.-p)/p0[n-1] );         
      q_minus = nq + enumerator/log(mn_minus);
      q_plus  = nq + enumerator/log(mn_plus);
      if ( fabs( ceil(q_plus) - ceil(q_minus) ) < .5 ) {
	 Wq = ceil(q_plus);
	 n = nmax +1;
      }   
    } /* n > q */  
  } /* p0[n-1] >= 1.-p */
 } /* n=q; n<=nmax; n++ */ 

 Free(Smatrix);
 Free(w);
 Free(z);
 Free(fn);
 Free(p0);

 return Wq;
}


double xe1_Warl(double l, double c, double zr, double hs,
double mu, int N, int nmax)
{ double *Pn, *w, *z, *p0, *atom, ratio, arl_minus=0., arl=1., arl_plus=0., mn_minus=1., mn_plus=0.;
  int i, j, n;

 c  *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );
 zr *= sqrt( l/(2.-l) );

 w  = vector(N);
 z  = vector(N);
 Pn = matrix(nmax,N);
 p0 = vector(nmax);
 atom = vector(nmax);

 gausslegendre(N,zr,c,z,w);

 for (n=1;n<=nmax;n++) {

   if (n==1) {
     for (i=0;i<N;i++)
       Pn[i] = PHI( (c-(1.-l)*z[i])/l, mu);
     atom[0] = PHI( (c-(1.-l)*zr)/l, mu);
   } else {
     for (i=0;i<N;i++) {
       Pn[(n-1)*N+i] = PHI( (zr-(1.-l)*z[i])/l, mu) * atom[n-2];
       for (j=0;j<N;j++)
         Pn[(n-1)*N+i] += w[j]/l * phi( (z[j]-(1.-l)*z[i])/l, mu) * Pn[(n-2)*N+j];
     }
     atom[n-1] = PHI( zr, mu) * atom[n-2];
     for (j=0;j<N;j++)
       atom[n-1] += w[j]/l * phi( (z[j]-(1.-l)*zr)/l, mu) * Pn[(n-2)*N+j];
   }
  
   if (n==1)
     p0[0] = PHI( (c-(1.-l)*hs)/l, mu);
   else {
     p0[n-1] = PHI( (zr-(1.-l)*hs)/l, mu) * atom[n-2];
     for (j=0;j<N;j++)
       p0[n-1] += w[j]/l * phi( (z[j]-(1.-l)*hs)/l, mu) * Pn[(n-2)*N+j];
   }

   mn_minus = 1.; mn_plus = 0.;
   if (n>1) {
     mn_minus = atom[n-1]/atom[n-2];
     mn_plus = atom[n-1]/atom[n-2];
     for (i=0;i<N;i++) {
       if (Pn[(n-2)*N+i]==0)
         if (Pn[(n-1)*N+i]==0) ratio = 0.;
         else ratio = 1.;
       else ratio = Pn[(n-1)*N+i]/Pn[(n-2)*N+i];
      if ( ratio<mn_minus ) mn_minus = ratio;
      if ( ratio>mn_plus ) mn_plus = ratio;
     }

     arl_minus = arl + p0[n-1]/(1.-mn_minus);
     arl_plus = arl + p0[n-1]/(1.-mn_plus);
   }
   arl += p0[n-1];

   if ( fabs( (arl_plus-arl_minus)/arl_minus )<FINALeps ) n = nmax+1;
 }

 Free(p0);
 Free(Pn);
 Free(z);
 Free(w);
 Free(atom);

 return (arl_plus+arl_minus)/2.;
}


double xe1_Wq(double l, double c, double p, double zr, double hs, double mu, int N, int nmax)
{ double *Pn, *w, *z, *p0, *atom, ratio, q_minus=0., q_plus=0., mn_minus=1., mn_plus=0., enumerator=0., Wq=0.;
  int i, j, n;

 c  *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );
 zr *= sqrt( l/(2.-l) );

 w  = vector(N);
 z  = vector(N);
 Pn = matrix(nmax,N);
 p0 = vector(nmax);
 atom = vector(nmax);

 gausslegendre(N,zr,c,z,w);

 for (n=1;n<=nmax;n++) {

   if (n==1) {
     for (i=0;i<N;i++)
       Pn[i] = PHI( (c-(1.-l)*z[i])/l, mu);
     atom[0] = PHI( (c-(1.-l)*zr)/l, mu);
   } else {
     for (i=0;i<N;i++) {
       Pn[(n-1)*N+i] = PHI( (zr-(1.-l)*z[i])/l, mu) * atom[n-2];
       for (j=0;j<N;j++) Pn[(n-1)*N+i] += w[j]/l * phi( (z[j]-(1.-l)*z[i])/l, mu) * Pn[(n-2)*N+j];
     }
     atom[n-1] = PHI( zr, mu) * atom[n-2];
     for (j=0;j<N;j++) atom[n-1] += w[j]/l * phi( (z[j]-(1.-l)*zr)/l, mu) * Pn[(n-2)*N+j];
   }

   if (n==1)
     p0[0] = PHI( (c-(1.-l)*hs)/l, mu);
   else {
     p0[n-1] = PHI( (zr-(1.-l)*hs)/l, mu) * atom[n-2];
     for (j=0;j<N;j++) p0[n-1] += w[j]/l * phi( (z[j]-(1.-l)*hs)/l, mu) * Pn[(n-2)*N+j];
   }

   if ( p0[n-1] < 1.-p ) {
     Wq = (double)n;
     n = nmax+1;
   } else {
     mn_minus = 1.; mn_plus = 0.;
     if ( n>1 ) {
       mn_minus = atom[n-1]/atom[n-2];
       mn_plus = atom[n-1]/atom[n-2];
       for (i=0;i<N;i++) {
         if (Pn[(n-2)*N+i]==0)
           if (Pn[(n-1)*N+i]==0) ratio = 0.;
           else ratio = 1.;
         else ratio = Pn[(n-1)*N+i]/Pn[(n-2)*N+i];
        if ( ratio<mn_minus ) mn_minus = ratio;
        if ( ratio>mn_plus ) mn_plus = ratio;
       }
       enumerator = log( (1.-p)/p0[n-1] );
       q_minus = (double)n + enumerator/log(mn_minus);
       q_plus  = (double)n + enumerator/log(mn_plus);
       if ( fabs( ceil(q_plus) - ceil(q_minus) ) < .5 ) {
	 Wq = ceil(q_plus);
	 n = nmax +1;
       }
     } /* n > 1 */
   } /* p0[n-1] >= 1.-p */
 } /* n=1; n<=nmax; n++ */ 

 Free(p0);
 Free(Pn);
 Free(z);
 Free(w);
 Free(atom);

 return Wq;
}


double xe1_Wqm(double l, double c, double p, double zr, double hs, int q, double mu0, double mu1, int mode, int N, int nmax)
{ double *Smatrix, *p0, *fn, *w, *z, rn, cn, rn0, cn0, q_minus=2., q_plus=3., mn_minus, mn_plus, nn, ratio, enumerator=0., nq, Wq=0.;
  int i, j, n, NN;

 c  *= sqrt( l/(2.-l) );
 zr *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );

 NN = N + 1;
 Smatrix = matrix(NN, NN);
 w       = vector(NN);
 z       = vector(NN);
 fn      = matrix(nmax, NN);
 p0      = vector(nmax);

 gausslegendre(N, zr, c, z, w);

 rn = 1.; cn = 0.; rn0 = 1., cn0 = 0.;

 /* in-control, i. e. n<=q-1 */
 for (n=1;n<=q-1;n++) {
  nn = (double) n;

  /* determine r_n, n=1,2,...,q-1 */
  if ( mode==vacl ) {
    rn = sqrt( 1. - pow(1.-l, 2.*nn) );
  }

  /* determine f_n, n=1,2,...,q-1 */
  if (n==1) {
    for (i=0;i<N;i++) {
      if ( mode==stat ) {
        fn[0*NN+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)), mu0);
      }
      else {
        fn[0*NN+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l, mu0);
      }
    }
    if ( mode==stat ) {
      fn[0*NN+N] = PHI( (cn+rn*zr)/sqrt(l/(2.-l)), mu0);
    }
    else {
      fn[0*NN+N] = PHI( (cn+rn*zr-(1.-l)*hs)/l, mu0);
    }
  }
  else {
    for (i=0;i<N;i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * rn/l * phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*zr))/l, mu0);
      for (j=0;j<N;j++) {
        fn[(n-1)*NN+i] += w[j]*fn[(n-2)*NN+j] * rn/l
                   *phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l, mu0);
      }
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI( (cn+rn*zr-(1.-l)*(cn0+rn0*zr))/l, mu0);
    for (j=0;j<N;j++)
      fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI( (cn+rn*zr-(1.-l)*(cn0+rn0*z[j]))/l, mu0);
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;
 }

 /* out-of-control, i.e. t>=q */
 for (n=q;n<=nmax;n++) {
  nn = (double) n;

  /* determine r_n, n=1,2,...,q-1 */
  if ( mode==vacl ) {
    rn = sqrt( 1. - pow(1.-l, 2.*nn) );
  }

  /* determine f_n, n=q,q+1,... */
  if (n==1) {
    for (i=0;i<N;i++) {
      if ( mode==stat ) {
        fn[0*NN+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)), mu1);
      }
      else {
        fn[0*NN+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l, mu1);
      }
    }
    if ( mode==stat ) {
      fn[0*NN+N] = PHI( (cn+rn*zr)/sqrt(l/(2.-l)), mu1);
    }
    else {
      fn[0*NN+N] = PHI( (cn+rn*zr-(1.-l)*hs)/l, mu1);
    }
  }
  else {
    for (i=0;i<N;i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * rn/l * phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*zr))/l, mu1);
      for (j=0;j<N;j++) {
        fn[(n-1)*NN+i] += w[j]*fn[(n-2)*NN+j] * rn/l * phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l, mu1);
      }
      if ( n==q && q>1 ) fn[(n-1)*NN+i] /= p0[q-2];
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI( (cn+rn*zr-(1.-l)*(cn0+rn0*zr))/l, mu1);
    for (j=0;j<N;j++)
      fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI( (cn+rn*zr-(1.-l)*(cn0+rn0*z[j]))/l, mu1);
    if (n==q && q>1) fn[(n-1)*NN+N] /= p0[q-2];
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];
  nq = (double)(n-q+1);

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;

  if ( p0[n-1] < 1.-p ) {
    Wq = nq;
    n = nmax+1;
  } else {     
    /* computation of m_n+1^- and m_n+1^+, n=m-1,m,... */
    mn_minus = 1.; mn_plus = 0.;
    if ( n > q ) {
      for (i=0;i<NN;i++) {
        if (fn[(n-2)*NN+i]==0)
        if (fn[(n-1)*NN+i]==0) ratio = 0.; else ratio = 1.;
        else ratio = fn[(n-1)*NN+i]/fn[(n-2)*NN+i];
        if ( ratio<mn_minus ) mn_minus = ratio;
        if ( ratio>mn_plus ) mn_plus = ratio;
      }
      enumerator = log( (1.-p)/p0[n-1] );         
      q_minus = nq + enumerator/log(mn_minus);
      q_plus  = nq + enumerator/log(mn_plus);
      if ( fabs( ceil(q_plus) - ceil(q_minus) ) < .5 ) {
	 Wq = ceil(q_plus);
	 n = nmax +1;
      }   
    } /* n > q */  
  } /* p0[n-1] >= 1.-p */
 } /* n=q; n<=nmax; n++ */ 

 Free(Smatrix);
 Free(w);
 Free(z);
 Free(fn);
 Free(p0);

 return Wq;
}


double xe2_Carl(double l, double c, double hs, double mu, int N, int qm)
{ double *a, *g, *w, *z, arl, Hij, zi, lzi, dN;
  int i, j, k;

 c  *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );

 dN = (double)N;

 a = matrix(N,N);
 g = vector(N);
 w = vector(qm);
 z = vector(qm);

 gausslegendre(qm,-c,c,z,w);

 for (i=0;i<N;i++) {
   zi = c * cos( (2.*(i+1.)-1.)*PI/2./dN );
   lzi = (1.-l)*zi;

   a[i*N] = 1 - ( PHI( (c-lzi)/l, mu) - PHI( (-c-lzi)/l, mu) );

   for (j=1;j<N;j++) {
     Hij = 0.;
     for (k=0;k<qm;k++) 
       Hij += w[k]/l * Tn( z[k]/c, j) * phi( (z[k]-lzi)/l, mu);
     a[i*N+j] = Tn( zi/c, j) - Hij;
   }
 }

 for (j=0;j<N;j++) g[j] = 1.;
 LU_solve(a,g,N);

 arl = g[0];
 for (j=1;j<N;j++) arl += g[j] * Tn( hs/c, j);

 Free(z);
 Free(w);
 Free(g);
 Free(a);

 return arl;
}


/* Manuel's PMS stuff */


double xe2_iglarl_RES
(double l, double c, double hs, double mu, int N, double alpha, int df)
{ double *a, *g, *w, *z, arl, ddf;
  int i, j;  

/* residual preliminaries */
 ddf = (double)df;
 mu *= ( 1. + ddf*sqrt( (1.-alpha)/(1.+alpha) ) )/(ddf+1.);
  
 a = matrix(N,N);
 g = vector(N);
 w = vector(N);
 z = vector(N);

 c  *= sqrt( l/(2.-l) ); 
 hs *= sqrt( l/(2.-l) );

 gausslegendre(N,-c,c,z,w);
  
 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) a[i*N+j] = -w[j]/l * phi( (z[j]-(1.-l)*z[i])/l,mu);
   ++a[i*N+i];
 }

 for (j=0;j<N;j++) g[j] = 1.;
 LU_solve(a,g,N);
  
 arl = 1.;
 for (j=0;j<N;j++)
   arl += w[j]/l * phi( (z[j]-(1.-l)*hs)/l,mu) * g[j];
 
 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}


double seU_iglarl_RES
  (double l, double cu, double hs, double sigma, int df, int N, int qm, double alpha, double mu)
{ double *a, *g, *w, *z, arl, Hij, xi, xl, za, xu, dN, ddf, s2, v, ncp;
  int i, j, k;

 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)N;
 
 /* residual preliminaries */
 mu *= ( 1. + ddf*sqrt( (1.-alpha)/(1.+alpha) ) )/(ddf+1.);
 ncp = ddf/(ddf+1.)*mu*mu/s2*pow( 1.-sqrt((1.-alpha)/(1.+alpha)), 2.);

 a = matrix(N,N);
 g = vector(N);
 w = vector(qm);
 z = vector(qm);

 for (i=0;i<N;i++) {
   xi = cu/2.*(1.+cos(PI*(2.*(i+1.)-1.)/2./dN));

   za = (1.-l)*xi;
   xl = 0.; 
   xu = sqrt(cu-za); 

   gausslegendre(qm,xl,xu,z,w);

   v = (cu - za)/l;
   a[i*N] = 1. - nCHI( ddf/s2*v, df, ncp);
   
   for (j=1;j<N;j++) {
     Hij = 0.;
     for (k=0;k<qm;k++) {
       v = (z[k] - za) / l;
       Hij += w[k] * Tn( (2.*(z[k]*z[k]+za)-cu)/cu ,j)
              * 2. * z[k]/l * ddf/s2 * nchi( ddf/s2*z[k]*z[k]/l, df, ncp);
     }
     a[i*N+j] = Tn( (2.*xi-cu)/cu ,j) - Hij;
   }
 }

 for (j=0;j<N;j++) g[j] = 1.;
 LU_solve(a,g,N);

 arl = g[0];
 for (j=1;j<N;j++)
   arl += g[j] * Tn( (2.*hs-cu)/cu ,j);

 Free(z);
 Free(w);
 Free(g);
 Free(a);

 return arl;
}


double xseU_arl_RES
  (double lx, double ls, double cx, double cs, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm, double alpha)
{ double *Sx, *Pnx, *wx, *zx, *p0x, *p0,
         *S1s, *S2s, *Pns, *ws, *zs, *p0s, q, *zch, *rside,
         za=0., s2,
         arl_minus=0., arl, arl_plus=0., mn_minus=1., mn_plus=0.,
         mn_minusx, mn_minuss, mn_plusx, mn_pluss, ddf, xl, xu,
         oben, unten, ncp;
  int i, j, k, n, *ps;

 cx  *= sqrt( lx/(2.-lx) ); 
 hsx *= sqrt( lx/(2.-lx) );

 s2 = sigma*sigma;
 ddf = (double)df;
 
 /* residual preliminaries */
 ncp = ddf/(ddf+1.)/(ddf+1.)*mu*mu/s2*pow( 1.-sqrt((1.-alpha)/(1.+alpha)), 2.);
 mu *= ( 1. + ddf*sqrt( (1.-alpha)/(1.+alpha) ) )/(ddf+1.);

 Sx  = matrix(Nx,Nx);
 wx  = vector(Nx);
 zx  = vector(Nx);
 Pnx = matrix(nmax,Nx);
 p0x = vector(nmax);

 S1s = matrix(Ns,Ns);
 S2s = matrix(Ns,Ns);
 ps = ivector(Ns);
 zch = vector(Ns);
 rside = vector(Ns);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,Ns);
 p0s = vector(nmax);

 p0  = vector(nmax);

 gausslegendre(Nx,-cx,cx,zx,wx);

 for (i=0;i<Nx;i++) {
   za = (1.-lx)*zx[i];
   for (j=0;j<Nx;j++)
     Sx[i*Nx+j] = wx[j]/lx*phi( ((zx[j]-za)/lx-mu)/sigma, 0.)/sigma;
 }  

/* Chebyshev nodes on [0,cs] */
 for (i=0;i<Ns;i++) 
   zch[i] = cs/2.*(1.+cos(PI*(2.*(i+1.)-1.)/2./(double)Ns) );

/* P(L>1)(zch[i]) */
 for (i=0;i<Ns;i++)
   rside[i] = nCHI( ddf/s2*(cs-(1.-ls)*zch[i])/ls, df, ncp); 

 for (i=0;i<Ns;i++) {
   za = (1.-ls)*zch[i];
   xl = 0.; xu = sqrt(cs-za);
   gausslegendre(qm,xl,xu,zs,ws);
   for (j=0;j<Ns;j++) {
     S1s[i*Ns+j] = 0.;
     for (k=0;k<qm;k++)
       S1s[i*Ns+j] += ws[k]*Tn((2.*(zs[k]*zs[k]+za)-cs)/cs, j)
                      * 2.*zs[k]/ls * ddf/s2 * nchi( ddf/s2 * zs[k]*zs[k]/ls, df, ncp);
   }
 }

 for (i=0;i<Ns;i++)
   for (j=0;j<Ns;j++) S2s[i*Ns+j] = Tn( (2.*zch[i]-cs)/cs, j);

 LU_decompose(S2s,ps,Ns);

 arl = 1.;
 for (n=1;n<=nmax;n++) {

   if (n==1)
     for (i=0;i<Nx;i++)
       Pnx[i] = PHI( (( cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.) - 
                PHI( ((-cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.);
   else
     for (i=0;i<Nx;i++) {
       Pnx[(n-1)*Nx+i] = 0.;
       for (j=0;j<Nx;j++)
         Pnx[(n-1)*Nx+i] += Sx[i*Nx+j] * Pnx[(n-2)*Nx+j];
     }

   p0x[n-1] = 0.;
   if (n==1)
     p0x[0] = PHI( (( cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.) - 
              PHI( ((-cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.);
   else
     for (j=0;j<Nx;j++)
       p0x[n-1] += wx[j]/lx * phi( ((zx[j]-(1.-lx)*hsx)/lx-mu)/sigma, 0.)/sigma * Pnx[(n-2)*Nx+j];


   if (n==1)
     for (i=0;i<Ns;i++) {
       Pns[i] = 0.;
       for (j=0;j<Ns;j++)
         Pns[i] += 2./Ns * Tn( (2.*zch[j]-cs)/cs, i) * rside[j];
       if (i==0) Pns[i] /= 2.;
     }
   else {
     for (i=0;i<Ns;i++) {
       rside[i] = 0.;
       for (j=0;j<Ns;j++) rside[i] += S1s[i*Ns+j] * Pns[(n-2)*Ns+j];
     }
     LU_solve2(S2s,rside,ps,Ns);
     for (i=0;i<Ns;i++) Pns[(n-1)*Ns+i] = rside[i];
   }

   p0s[n-1] = 0.;  
   if (n==1)
     p0s[0] = nCHI(ddf/s2*(cs-(1.-ls)*hss)/ls, df, ncp);
   else
     for (j=0;j<Ns;j++)
       p0s[n-1] += Pns[(n-1)*Ns+j] * Tn( (2.*hss-cs)/cs, j);

   p0[n-1] = p0x[n-1] * p0s[n-1];

   mn_minusx = 1.; mn_plusx = 0.;
   mn_minuss = 1.; mn_pluss = 0.;
   if (n>1) {
     for (i=0;i<Nx;i++) {
       if (Pnx[(n-1)*Nx+i]==0)
         if (Pnx[(n-1)*Nx+i]==0) q = 0.;
         else q = 1.;
       else q = Pnx[(n-1)*Nx+i]/Pnx[(n-2)*Nx+i];
      if ( q<mn_minusx ) mn_minusx = q;
      if ( q>mn_plusx ) mn_plusx = q;
     }

     for (i=0;i<Ns;i++) {
       oben = 0.; unten = 0.;
       for (j=0;j<Ns;j++) {
         oben += Pns[(n-1)*Ns+j] * Tn( (2.*zch[i]-cs)/cs, j);
         unten+= Pns[(n-2)*Ns+j] * Tn( (2.*zch[i]-cs)/cs, j);
       }
       if (fabs(unten)<1e-16)
         if (fabs(oben)<1e-16) q = 0.;
         else q = 1.;
       else q = oben/unten;
      if ( q<mn_minuss ) mn_minuss = q;
      if ( q>mn_pluss ) mn_pluss = q;
     }

     mn_minus = mn_minusx * mn_minuss;
     mn_plus  = mn_plusx * mn_pluss;

     arl_minus = arl + p0[n-1]/(1.-mn_minus);
     arl_plus = arl + p0[n-1]/(1.-mn_plus);
   }
   arl += p0[n-1];

   if ( fabs( (arl_plus-arl_minus)/arl_minus )<FINALeps ) n = nmax+1;
 }

 Free(p0);

 Free(p0s);
 Free(Pns);
 Free(zs);
 Free(ws);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);
 
 Free(p0x);
 Free(Pnx);
 Free(zx);
 Free(wx);
 Free(Sx);

 return (arl_plus+arl_minus)/2.;
}


double xseU_mu_before_sigma_RES
  (double lx, double ls, double cx, double cs, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm, double alpha, int vice_versa)
{ double *Sx, *Pnx, *wx, *zx, *p0x,
         *S1s, *S2s, *Pns, *ws, *zs, *p0s, *zch, *rside,
         za=0., s2, mu_before_sigma=0., ddf, xl, xu, ncp;
  int i, j, k, n, *ps;
 
 cx  *= sqrt( lx/(2.-lx) );
 hsx *= sqrt( lx/(2.-lx) );

 s2 = sigma*sigma;
 ddf = (double)df;
 
 /* residual preliminaries */
 ncp = ddf/(ddf+1.)/(ddf+1.)*mu*mu/s2*pow( 1.-sqrt((1.-alpha)/(1.+alpha)), 2.);
 mu *= ( 1. + ddf*sqrt( (1.-alpha)/(1.+alpha) ) )/(ddf+1.);

 Sx  = matrix(Nx,Nx);
 wx  = vector(Nx);
 zx  = vector(Nx);
 Pnx = matrix(nmax,Nx);
 p0x = vector(nmax);

 S1s = matrix(Ns,Ns);
 S2s = matrix(Ns,Ns);
 ps = ivector(Ns);
 zch = vector(Ns);
 rside = vector(Ns);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,Ns);
 p0s = vector(nmax);

 gausslegendre(Nx,-cx,cx,zx,wx);

 for (i=0;i<Nx;i++) {
   za = (1.-lx)*zx[i];
   for (j=0;j<Nx;j++)
     Sx[i*Nx+j] = wx[j]/lx*phi( ((zx[j]-za)/lx-mu)/sigma, 0.)/sigma;
 }  

/* Chebyshev nodes on [0,cs] */
 for (i=0;i<Ns;i++) 
   zch[i] = cs/2.*(1.+cos(PI*(2.*(i+1.)-1.)/2./(double)Ns) );

/* P(L>1)(zch[i]) */
 for (i=0;i<Ns;i++)
   rside[i] = nCHI( ddf/s2*(cs-(1.-ls)*zch[i])/ls, df, ncp); 

 for (i=0;i<Ns;i++) {
   za = (1.-ls)*zch[i];
   xl = 0.; xu = sqrt(cs-za);
   gausslegendre(qm,xl,xu,zs,ws);
   for (j=0;j<Ns;j++) {
     S1s[i*Ns+j] = 0.;
     for (k=0;k<qm;k++)
       S1s[i*Ns+j] += ws[k]*Tn((2.*(zs[k]*zs[k]+za)-cs)/cs, j)
                      * 2.*zs[k]/ls * ddf/s2 * nchi( ddf/s2 * zs[k]*zs[k]/ls, df, ncp);
   }
 }

 for (i=0;i<Ns;i++)
   for (j=0;j<Ns;j++) S2s[i*Ns+j] = Tn( (2.*zch[i]-cs)/cs, j);

 LU_decompose(S2s,ps,Ns);

 mu_before_sigma = 0.;
 for (n=1;n<=nmax;n++) {

   if (n==1)
     for (i=0;i<Nx;i++)
       Pnx[i] = PHI( (( cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.) - 
                PHI( ((-cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.);
   else
     for (i=0;i<Nx;i++) {
       Pnx[(n-1)*Nx+i] = 0.;
       for (j=0;j<Nx;j++)
         Pnx[(n-1)*Nx+i] += Sx[i*Nx+j] * Pnx[(n-2)*Nx+j];
     }

   p0x[n-1] = 0.;
   if (n==1)
     p0x[0] = PHI( (( cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.) - 
              PHI( ((-cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.);
   else
     for (j=0;j<Nx;j++)
       p0x[n-1] += wx[j]/lx * phi( ((zx[j]-(1.-lx)*hsx)/lx-mu)/sigma, 0.)/sigma * Pnx[(n-2)*Nx+j];


   if (n==1)
     for (i=0;i<Ns;i++) {
       Pns[i] = 0.;
       for (j=0;j<Ns;j++)
         Pns[i] += 2./Ns * Tn( (2.*zch[j]-cs)/cs, i) * rside[j];
       if (i==0) Pns[i] /= 2.;
     }
   else {
     for (i=0;i<Ns;i++) {
       rside[i] = 0.;
       for (j=0;j<Ns;j++) rside[i] += S1s[i*Ns+j] * Pns[(n-2)*Ns+j];
     }
     LU_solve2(S2s,rside,ps,Ns);
     for (i=0;i<Ns;i++) Pns[(n-1)*Ns+i] = rside[i];
   }

   p0s[n-1] = 0.;  
   if (n==1)
     p0s[0] = nCHI(ddf/s2*(cs-(1.-ls)*hss)/ls, df, ncp);
   else
     for (j=0;j<Ns;j++)
       p0s[n-1] += Pns[(n-1)*Ns+j] * Tn( (2.*hss-cs)/cs, j);
 
   if ( vice_versa ) { /* S chart before X chart -- PMS IV */
     if (n>1)
       mu_before_sigma += ( p0s[n-2] - p0s[n-1] ) * p0x[n-1];
     else
       mu_before_sigma = ( 1. - p0s[n-1] ) * p0x[n-1];
     if ( p0s[n-1]<FINALeps ) n = nmax+1;
   } else { /* X chart before S chart -- PMS III */
     if (n>1)
       mu_before_sigma += ( p0x[n-2]-p0x[n-1] ) * p0s[n-1];
     else
       mu_before_sigma = ( 1.-p0x[n-1] ) * p0s[n-1];
     if ( p0x[n-1]<FINALeps ) n = nmax+1;
   }
 }

 Free(p0s);
 Free(Pns);
 Free(zs);
 Free(ws);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);
 
 Free(p0x);
 Free(Pnx);
 Free(zx);
 Free(wx);
 Free(Sx);

 return mu_before_sigma;
}


/* end of Manuel's stuff */


/* For Christian */

/* Shewhart charts for dependent data, Gaussian case */
double x_shewhart_ar1_arl(double alpha, double cS, double delta, int N1, int N2)
{ double *a, *g, *w1, *z1, *w2, *z2, arl, arl1, mdelta, l, kappa, cE;
  int i, j;

 a  = matrix(N1,N1);
 g  = vector(N1);
 w1 = vector(N1);
 z1 = vector(N1);
 w2 = vector(N2);
 z2 = vector(N2); 
 
 l = 1. - alpha;
 kappa = sqrt( (1. - alpha) / (1. + alpha) ); /* = sqrt( lambda / ( 2 - lambda ) ) */
 mdelta = kappa*delta;
 cE = kappa*cS;
 
 gausslegendre(N1, -cE, cE, z1, w1);

 for (i=0; i<N1; i++) {
   for (j=0; j<N1; j++) a[i*N1+j] = - w1[j]/l * phi( ( z1[j] - (1.-l)*z1[i] )/l, mdelta);
   ++a[i*N1 + i];
 }

 for (j=0; j<N1; j++) g[j] = 1.;
 LU_solve(a, g, N1);
 
 gausslegendre(N2, -cS, cS, z2, w2);
 
 arl = 1.;
 for (i=0; i<N2; i++) {
    arl1 = 1.; 
    for (j=0; j<N1; j++) arl1 += w1[j]/l * phi( ( z1[j] - (1.-l)*kappa*z2[i] )/l, mdelta) * g[j];
    arl += w2[i] * phi(z2[i], delta) * arl1;
 }

 Free(a);
 Free(g);
 Free(w1);
 Free(z1);
 Free(w2);
 Free(z2);

 return arl;
}


/* Shewhart charts for dependent data, Student t case */
double t_shewhart_ar1_arl(double alpha, double cS, double delta, int df, int N1, int N2, int N3, double INF, int subst)
{ double *a, *g, *w1, *z1, *w2, *z2, *a3, *w3, *z3, *psi, arl, arl1, mdelta, l, kappa, cE, k1, norm=1., arg, ddf, k2, hs, k3, rho, *pdf;
  int i, j, status, noofit;

 a  = matrix(N1,N1);
 g  = vector(N1);
 w1 = vector(N1);
 z1 = vector(N1);
 
 w2 = vector(N2);
 z2 = vector(N2);
 pdf = vector(N2);
 
 w3 = vector(N3);
 z3 = vector(N3);
 a3 = matrix(N3,N3);
 psi = vector(N3);
 
 l = 1. - alpha;
 kappa = sqrt( (1. - alpha) / (1. + alpha) ); /* = sqrt( lambda / ( 2 - lambda ) ) */
 mdelta = kappa*delta; 
 cE = kappa*cS;
 
 ddf = (double)df;
 k2 = sqrt( ddf/(ddf-2.) );
 k3 = 1. / sqrt( 1. - alpha*alpha );

 switch ( subst ) {
   case IDENTITY: gausslegendre(N1, -cE, cE, z1, w1); norm = 1.; break;
   case SIN:   gausslegendre(N1, -PI/2., PI/2., z1, w1); norm = 1.; break;
   case SINH:  gausslegendre(N1, -1., 1., z1, w1); norm = sinh(1.); break;
   case TAN:   gausslegendre(N1, -PI/4., PI/4., z1, w1); norm = 1.; break;
 }     
 
 cE /= norm;

 for (i=0; i<N1; i++) {
   for (j=0; j<N1; j++) {
     /* IDENTITY as default */  
     arg = z1[j] - (1.-l)*z1[i];
     k1 = 1.;       
     switch ( subst ) {
       case IDENTITY: arg = z1[j] - (1.-l)*z1[i]; k1 = 1.; break;
       case SIN:   arg = cE*sin(z1[j]) - (1.-l)*cE*sin(z1[i]); k1 = cE*cos(z1[j]); break;
       case SINH:  arg = cE*sinh(z1[j]) - (1.-l)*cE*sinh(z1[i]); k1 = cE*cosh(z1[j]); break;
       case TAN:   arg = cE*tan(z1[j]) - (1.-l)*cE*tan(z1[i]); k1 = cE/( cos(z1[j])*cos(z1[j]) ); break;
     }
     a[i*N1+j] = - w1[j]/l * k2*pdf_t( k2*(arg/l - mdelta), df) * k1;
   }
   ++a[i*N1 + i];
 }

 for (j=0; j<N1; j++) g[j] = 1.;
 LU_solve(a, g, N1);
 
/* stationary distribution */
 gausslegendre(N3, -INF, INF, z3, w3);
 for (i=0;i<N3;i++)
   for (j=0;j<N3;j++) a3[i*N3+j] = w3[j] * pdf_t( (z3[i] - alpha*z3[j] - (1.-alpha)*delta)*k2*k3, df) * k2*k3;

 pmethod(N3, a3, &status, &rho, psi, &noofit);
 
 norm = 0.;
 for (j=0; j<N3; j++) norm += w3[j] * psi[j];
 
 gausslegendre(N2, -cS, cS, z2, w2);
 
 for (i=0; i<N2; i++) {
   pdf[i] = 0.;
   for (j=0; j<N3; j++) pdf[i] += w3[j] * psi[j] * pdf_t( (z2[i] - alpha*z3[j] - (1.-alpha)*delta)*k2*k3, df) * k2*k3;
   pdf[i] /= norm;
 }
 
 arl = 1.;
 for (i=0; i<N2; i++) {
    arl1 = 1.;     
    for (j=0; j<N1; j++) {
      hs = kappa*z2[i];
      /* IDENTITY as default */
      arg = z1[j] - (1.-l)*hs;
      k1 = 1.;
      switch ( subst ) {
        case IDENTITY: arg = z1[j] - (1.-l)*hs; k1 = 1.; break;
        case SIN:   arg = cE*sin(z1[j]) - (1.-l)*hs; k1 = cE*cos(z1[j]); break;
        case SINH:  arg = cE*sinh(z1[j]) - (1.-l)*hs; k1 = cE*cosh(z1[j]); break;
        case TAN:   arg = cE*tan(z1[j]) - (1.-l)*hs; k1 = cE/( cos(z1[j])*cos(z1[j]) ); break;
      }
     arl1 += w1[j]/l * k2*pdf_t( k2*(arg/l - mdelta), df) * g[j] * k1;
   }
   arl += w2[i] * pdf[i] * arl1;
 }

 Free(a);
 Free(g);
 Free(w1);
 Free(z1);
 
 Free(w2);
 Free(z2);
 Free(pdf);
  
 Free(a3);
 Free(w3);
 Free(z3);
 Free(psi);
 
 return arl;
}


/* end of Christian's stuff */


double xc1_iglad (double k, double h, double mu0, double mu1, int N)
{ double *a, *w, *z, *arl, *psi, rho, ad, norm;
  int i, j, status, noofit, NN;

 NN = N + 1;
 a = matrix(NN,NN);
 arl = vector(NN);
 psi = vector(NN);
 w = vector(NN);
 z = vector(NN);

 gausslegendre(N,0.,h,z,w);

 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) a[i*NN+j] = -w[j]*phi(z[j]+k-z[i],mu1); 
   ++a[i*NN+i];
   a[i*NN+N] = - PHI(k-z[i],mu1);
 }
 for (j=0;j<N;j++)
    a[N*NN+j] = -w[j]*phi( z[j]+k,mu1);
 a[N*NN+N] = 1. - PHI(k,mu1);

 for (j=0;j<NN;j++) arl[j] = 1.;
 LU_solve(a,arl,NN);

 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) a[i*NN+j] = w[j]*phi(z[i]+k-z[j],mu0);
   a[i*NN+N] = phi(z[i]+k,mu0); 
 }
 for (j=0;j<N;j++)
    a[N*NN+j] = w[j] * PHI(k-z[j],mu0);
 a[N*NN+N] = PHI(k,mu0);

 pmethod(NN,a,&status,&rho,psi,&noofit);

 ad = psi[N]*arl[N];
 norm = psi[N];
 for (j=0;j<N;j++) {
   ad += w[j] * arl[j] * psi[j];
   norm += w[j] * psi[j];
 }
 ad /= norm;
 rho0 = rho;

 Free(a);
 Free(arl);
 Free(psi);
 Free(w);
 Free(z);

 return ad;
}


double xcC_iglad (double k, double h, double mu0, double mu1, int N)
{ double *a, *w, *z, *arl, *psi, rho, ad, norm;
  int i, j, status, noofit, NN;

 NN = 2*N + 1;
 a = matrix(NN,NN);
 arl = vector(NN);
 psi = vector(NN);
 w = vector(NN);
 z = vector(NN);

 gausslegendre(N,0.,h,z,w);

 for (i=0;i<N;i++) { /* upper */
   for (j=0;j<N;j++)    a[i*NN+j] = -w[j]  *phi( z[j]  +k-z[i],mu1);
   for (j=N;j<NN-1;j++) a[i*NN+j] = -w[j-N]*phi(-z[j-N]-k-z[i],mu1);
   ++a[i*NN+i];
   a[i*NN+NN-1] = - ( PHI(k-z[i],mu1) - PHI(-k-z[i],mu1) );
 }

 for (i=N;i<NN-1;i++) { /* lower */
   for (j=0;j<N;j++)    a[i*NN+j] = -w[j]  *phi( z[j]  +k+z[i-N],mu1);
   for (j=N;j<NN-1;j++) a[i*NN+j] = -w[j-N]*phi(-z[j-N]-k+z[i-N],mu1);
   ++a[i*NN+i];
   a[i*NN+NN-1] = - ( PHI(k+z[i-N],mu1) - PHI(-k+z[i-N],mu1) );
 }

 /* "fat" zero */
 for (j=0;j<N;j++)
    a[(NN-1)*NN+j] = -w[j] * phi( z[j]  +k,mu1);
 for (j=N;j<NN-1;j++)
    a[(NN-1)*NN+j] = -w[j-N]*phi(-z[j-N]-k,mu1);
 a[(NN-1)*NN+NN-1] = 1. - ( PHI(k,mu1) - PHI(-k,mu1) );

 for (j=0;j<NN;j++) arl[j] = 1.;
 LU_solve(a,arl,NN);

 for (i=0;i<N;i++) { /* upper */
   for (j=0;j<N;j++)    a[i*NN+j] = w[j]  *phi( z[i]+k-z[j]  ,mu0);
   for (j=N;j<NN-1;j++) a[i*NN+j] = w[j-N]*phi( z[i]+k+z[j-N],mu0);
   a[i*NN+NN-1] = phi(z[i]+k,mu0);
 }

 for (i=N;i<NN-1;i++) { /* lower */
   for (j=0;j<N;j++)    a[i*NN+j] = w[j]  *phi( -z[i-N]-k-z[j]  ,mu0);
   for (j=N;j<NN-1;j++) a[i*NN+j] = w[j-N]*phi( -z[i-N]-k+z[j-N],mu0);
   a[i*NN+NN-1] = phi(-z[i-N]-k,mu0);
 }

 /* "fat" zero */
 for (j=0;j<N;j++)
    a[(NN-1)*NN+j] = w[j]  * ( PHI(k-z[j]  ,mu0) - PHI(-k-z[j]  ,mu0) );
 for (j=N;j<NN-1;j++)
    a[(NN-1)*NN+j] = w[j-N]* ( PHI(k+z[j-N],mu0) - PHI(-k+z[j-N],mu0) );
 a[(NN-1)*NN+NN-1] = PHI(k,mu0) - PHI(-k,mu0);

 pmethod(NN,a,&status,&rho,psi,&noofit);

 ad = psi[NN-1]*arl[NN-1];
 norm = psi[NN-1];
 for (j=0;j<N;j++) {
   ad += w[j] * arl[j] * psi[j];
   norm += w[j] * psi[j];
 }
 for (j=N;j<NN-1;j++) {
   ad += w[j-N] * arl[j] * psi[j];
   norm += w[j-N] * psi[j];
 }
 ad /= norm;
 rho0 = rho;

 Free(a);
 Free(arl);
 Free(psi);
 Free(w);
 Free(z);

 return ad;
}


double xe1_iglad (double l, double c, double zr, double mu0, double mu1, int N)
{ double *a, *w, *z, *arl, *psi, rho, ad, norm;
  int i, j, status, noofit, NN;

 NN = N + 1;
 a = matrix(NN,NN);
 arl = vector(NN);
 psi = vector(NN);
 w = vector(NN);
 z = vector(NN);

 c  *= sqrt( l/(2.-l) );
 zr *= sqrt( l/(2.-l) );

 gausslegendre(N,zr,c,z,w);

 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) a[i*NN+j] = -w[j]/l * phi( (z[j]-(1.-l)*z[i])/l,mu1);
   ++a[i*NN+i];
   a[i*NN+N] = - PHI((zr-(1.-l)*z[i])/l,mu1);
 }
 for (j=0;j<N;j++)
    a[N*NN+j] = -w[j]/l * phi( (z[j]-(1.-l)*zr)/l,mu1);
 a[N*NN+N] = 1. - PHI(zr,mu1);

 for (j=0;j<NN;j++) arl[j] = 1.;
 LU_solve(a,arl,NN); 

 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) a[i*NN+j] = w[j]/l * phi((z[i]-(1.-l)*z[j])/l,mu0);
   a[i*NN+N] = 1./l * phi((z[i]-(1.-l)*zr)/l,mu0);
 }
 for (j=0;j<N;j++)
    a[N*NN+j] = w[j] * PHI((zr-(1.-l)*z[j])/l,mu0);
 a[N*NN+N] = PHI(zr,mu0);

 pmethod(NN,a,&status,&rho,psi,&noofit);

 ad = psi[N]*arl[N];
 norm = psi[N];
 for (j=0;j<N;j++) {
   ad += w[j] * arl[j] * psi[j];
   norm += w[j] * psi[j];
 }
 ad /= norm;
 rho0 = rho;

 Free(a);
 Free(arl);
 Free(psi);
 Free(w);
 Free(z);

 return ad;
}


double xe2_iglad (double l, double c, double mu0, double mu1, int N)
{ double *a, *w, *z, *arl, *psi, rho, ad, norm;
  int i, j, status, noofit;

 a = matrix(N,N);
 arl = vector(N);
 psi = vector(N);
 w = vector(N);
 z = vector(N);

 c *= sqrt( l/(2.-l) );

 gausslegendre(N,-c,c,z,w);

 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) a[i*N+j] = -w[j]/l * phi((z[j]-(1.-l)*z[i])/l,mu1);
   ++a[i*N+i];
 }

 for (j=0;j<N;j++) arl[j] = 1.;
 LU_solve(a,arl,N);

 for (i=0;i<N;i++)
   for (j=0;j<N;j++) a[i*N+j] = w[j]/l * phi((z[i]-(1.-l)*z[j])/l,mu0);

 pmethod(N,a,&status,&rho,psi,&noofit);

 ad = 0.; norm = 0.;
 for (j=0;j<N;j++) {
   ad += w[j] * arl[j] * psi[j];
   norm += w[j] * psi[j];
 }
 ad /= norm;
 rho0 = rho;

 Free(a);
 Free(arl);
 Free(psi);
 Free(w);
 Free(z);

 return ad;
}


double xte2_iglad (double l, double c, int df, double mu0, double mu1, int N, int subst)
{ double *a, *w, *z, *arl, *psi, rho, ad, nenner, norm=1., arg=0., korr=1.;
  int i, j, status, noofit;

 a = matrix(N,N);
 arl = vector(N);
 psi = vector(N);
 w = vector(N);
 z = vector(N);

 c *= sqrt( l/(2.-l) );

 switch ( subst ) {
   case IDENTITY: gausslegendre(N, -c, c, z, w);         norm = 1.; break;
   case SIN:      gausslegendre(N, -PI/2., PI/2., z, w); norm = 1.; break;
   case SINH:     gausslegendre(N, -1., 1., z, w);       norm = sinh(1.); break;
   case TAN:      gausslegendre(N, -PI/4., PI/4., z, w); norm = 1.; break;
 }     
 
 c /= norm;

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) {
     switch ( subst ) {
       case IDENTITY: arg = z[j] - (1.-l)*z[i];                 korr = 1.; break;
       case SIN:      arg = c*sin(z[j]) - (1.-l)*c*sin(z[i]);   korr = c*cos(z[j]); break;
       case SINH:     arg = c*sinh(z[j]) - (1.-l)*c*sinh(z[i]); korr = c*cosh(z[j]); break;
       case TAN:      arg = c*tan(z[j]) - (1.-l)*c*tan(z[i]);   korr = c/( cos(z[j])*cos(z[j]) ); break;
     }
     a[i*N+j] = -w[j]/l * pdf_t( arg/l - mu1, df) * korr;
   }
   ++a[i*N+i];
 }

 for (j=0;j<N;j++) arl[j] = 1.;
 LU_solve(a, arl, N);

 for (i=0; i<N; i++)
   for (j=0; j<N; j++) {
     switch ( subst ) {
       case IDENTITY: arg = z[i] - (1.-l)*z[j];                 korr = 1.; break;
       case SIN:      arg = c*sin(z[i]) - (1.-l)*c*sin(z[j]);   korr = c*cos(z[j]); break;
       case SINH:     arg = c*sinh(z[i]) - (1.-l)*c*sinh(z[j]); korr = c*cosh(z[j]); break;
       case TAN:      arg = c*tan(z[i]) - (1.-l)*c*tan(z[j]);   korr = c/( cos(z[j])*cos(z[j]) ); break;
     }
     a[i*N+j] = w[j]/l * pdf_t( arg/l - mu0, df) * korr;
   }

 pmethod(N, a, &status, &rho, psi, &noofit);

 ad = 0.; nenner = 0.;
 for (j=0; j<N; j++) {
   switch ( subst ) {
     case IDENTITY: korr = 1.; break;
     case SIN:      korr = c*cos(z[j]); break;
     case SINH:     korr = c*cosh(z[j]); break;
     case TAN:      korr = c/( cos(z[j])*cos(z[j]) ); break;
   }
   ad += w[j] * arl[j] * psi[j] * korr;
   nenner += w[j] * psi[j] * korr;
 }
 ad /= nenner;
 rho0 = rho;

 Free(a);
 Free(arl);
 Free(psi);
 Free(w);
 Free(z);

 return ad;
}


double xe2_igladc(double l, double c, double mu0, double mu1, double z0, int N)
{ double *a, *w, *z, *arl, *psi, rho, ad, norm, L0;
  int i, j, status, noofit, NN;

 NN = N + 1;
 a = matrix(NN,NN);
 arl = vector(N);
 psi = vector(NN);
 w = vector(N);
 z = vector(N);

 c *= sqrt( l/(2.-l) );
 z0 *= sqrt( l/(2.-l) );

 gausslegendre(N, -c, c, z, w);

 /* ooc vector */
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*N+j] = -w[j]/l * phi((z[j]-(1.-l)*z[i])/l, mu1);
   ++a[i*N+i];
 }

 for (j=0; j<N; j++) arl[j] = 1.;
 LU_solve(a, arl, N);
 
 /* ooc ARL at restart point */
 L0 = 1.;
 for (j=0; j<N; j++) L0 += w[j]/l * phi( (z[j]-(1.-l)*z0)/l, mu1) * arl[j];
 

 /* left eigenvector psi */
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*NN+j] = w[j]/l * phi( (z[i]-(1.-l)*z[j])/l, mu0);
   a[i*NN+N] = 1./l * phi( (z[i]-(1.-l)*z0)/l, mu0);
 }
 for (j=0;j<N;j++) a[N*NN+j] = w[j] * ( 1 - PHI( (c-(1.-l)*z[j])/l, mu0) + PHI( (-c-(1.-l)*z[j])/l, mu0) );
 a[N*NN+N] = 1 - PHI( (c-(1.-l)*z0)/l, mu0) + PHI( (-c-(1.-l)*z0)/l, mu0);

 pmethod(NN, a, &status, &rho, psi, &noofit);

 ad = L0 * psi[N];
 norm = psi[N];
 for (j=0; j<N; j++) {
   ad += w[j] * arl[j] * psi[j];
   norm += w[j] * psi[j];
 }
 ad /= norm;
 rho0 = rho;

 Free(a);
 Free(arl);
 Free(psi);
 Free(w);
 Free(z);

 return ad;
}


double xte2_igladc(double l, double c, int df, double mu0, double mu1, double z0, int N, int subst)
{ double *a, *w, *z, *arl, *psi, rho, ad, nenner, L0, norm=1., arg=0., korr=1.;
  int i, j, status, noofit, NN;

 NN = N + 1;
 a = matrix(NN,NN);
 arl = vector(N);
 psi = vector(NN);
 w = vector(N);
 z = vector(N);

 c *= sqrt( l/(2.-l) );
 z0 *= sqrt( l/(2.-l) );

 switch ( subst ) {
   case IDENTITY: gausslegendre(N, -c, c, z, w);         norm = 1.; break;
   case SIN:      gausslegendre(N, -PI/2., PI/2., z, w); norm = 1.; break;
   case SINH:     gausslegendre(N, -1., 1., z, w);       norm = sinh(1.); break;
   case TAN:      gausslegendre(N, -PI/4., PI/4., z, w); norm = 1.; break;
 }     
 
 c /= norm;

 /* ooc vector */
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) {
     switch ( subst ) {
       case IDENTITY: arg = z[j] - (1.-l)*z[i];                 korr = 1.; break;
       case SIN:      arg = c*sin(z[j]) - (1.-l)*c*sin(z[i]);   korr = c*cos(z[j]); break;
       case SINH:     arg = c*sinh(z[j]) - (1.-l)*c*sinh(z[i]); korr = c*cosh(z[j]); break;
       case TAN:      arg = c*tan(z[j]) - (1.-l)*c*tan(z[i]);   korr = c/( cos(z[j])*cos(z[j]) ); break;
     }
     a[i*N+j] = -w[j]/l * pdf_t( arg/l - mu1, df) * korr;
   }
   ++a[i*N+i];
 }

 for (j=0; j<N; j++) arl[j] = 1.;
 LU_solve(a, arl, N);
 
 
 /* ooc ARL at restart point */
 L0 = 1.;
 for (j=0; j<N; j++) {
   switch ( subst ) {
     case IDENTITY: arg = z[j];         korr = 1.; break;
     case SIN:      arg = c*sin(z[j]);  korr = c*cos(z[j]); break;
     case SINH:     arg = c*sinh(z[j]); korr = c*cosh(z[j]); break;
     case TAN:      arg = c*tan(z[j]);  korr = c/( cos(z[j])*cos(z[j]) ); break;
   }
   L0 += w[j]/l * pdf_t( ( arg - (1.-l)*z0 )/l - mu1, df) * korr * arl[j];
 } 
 
 /* left eigenvector psi */
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) {
     switch ( subst ) {
       case IDENTITY: arg = z[i] - (1.-l)*z[j];                 korr = 1.; break;
       case SIN:      arg = c*sin(z[i]) - (1.-l)*c*sin(z[j]);   korr = c*cos(z[j]); break;
       case SINH:     arg = c*sinh(z[i]) - (1.-l)*c*sinh(z[j]); korr = c*cosh(z[j]); break;
       case TAN:      arg = c*tan(z[i]) - (1.-l)*c*tan(z[j]);   korr = c/( cos(z[j])*cos(z[j]) ); break;
     }
     a[i*NN+j] = w[j]/l * pdf_t( arg/l - mu0, df) * korr;
   }
   switch ( subst ) {
     case IDENTITY: arg = z[i]; break;
     case SIN:      arg = c*sin(z[i]); break;
     case SINH:     arg = c*sinh(z[i]); break;
     case TAN:      arg = c*tan(z[i]); break;
   }
   a[i*NN+N] = 1./l * pdf_t( ( arg - (1.-l)*z0 )/l - mu0, df);
 }   
 for (j=0;j<N;j++) {
   switch ( subst ) {
     case IDENTITY: arg = z[j];         korr = 1.; break;
     case SIN:      arg = c*sin(z[j]);  korr = c*cos(z[j]); break;
     case SINH:     arg = c*sinh(z[j]); korr = c*cosh(z[j]); break;
     case TAN:      arg = c*tan(z[j]);  korr = c/( cos(z[j])*cos(z[j]) ); break;
   }
   a[N*NN+j] = w[j] * ( 1. - cdf_t( ( c*norm - (1.-l)*arg )/l - mu0, df) + cdf_t( ( -c*norm - (1.-l)*arg )/l - mu0, df) ) * korr;
 } 
 a[N*NN+N] = 1. - cdf_t( ( c*norm - (1.-l)*z0 )/l - mu0, df) + cdf_t( ( -c*norm - (1.-l)*z0 )/l - mu0, df);
 
 pmethod(NN, a, &status, &rho, psi, &noofit);

 ad = L0 * psi[N];
 nenner = psi[N];
 for (j=0; j<N; j++) {
   switch ( subst ) {
     case IDENTITY: korr = 1.; break;
     case SIN:      korr = c*cos(z[j]); break;
     case SINH:     korr = c*cosh(z[j]); break;
     case TAN:      korr = c/( cos(z[j])*cos(z[j]) ); break;
   }
   ad += w[j] * arl[j] * psi[j] * korr;
   nenner += w[j] * psi[j] * korr;
 }
 ad /= nenner;
 rho0 = rho;

 Free(a);
 Free(arl);
 Free(psi);
 Free(w);
 Free(z);

 return ad;
}


double xe1_sf(double l, double c, double zr, double hs, double mu, int N, int nmax, double *p0)
{ double *Pn, *w, *z, *atom;
  int i, j, n;

 c  *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );
 zr *= sqrt( l/(2.-l) );

 w  = vector(N);
 z  = vector(N);
 Pn = matrix(nmax,N);
 atom = vector(nmax);

 gausslegendre(N,zr,c,z,w);

 for (n=1;n<=nmax;n++) {
   if (n==1) {
     for (i=0;i<N;i++)
       Pn[i] = PHI( (c-(1.-l)*z[i])/l, mu);
     atom[0] = PHI( (c-(1.-l)*zr)/l, mu);
   } else {
     for (i=0;i<N;i++) {
       Pn[(n-1)*N+i] = PHI( (zr-(1.-l)*z[i])/l, mu) * atom[n-2];
       for (j=0;j<N;j++) Pn[(n-1)*N+i] += w[j]/l * phi( (z[j]-(1.-l)*z[i])/l, mu) * Pn[(n-2)*N+j];
     }
     atom[n-1] = PHI( zr, mu) * atom[n-2];
     for (j=0;j<N;j++) atom[n-1] += w[j]/l * phi( (z[j]-(1.-l)*zr)/l, mu) * Pn[(n-2)*N+j];
   }

   if (n==1)
     p0[0] = PHI( (c-(1.-l)*hs)/l, mu);
   else {
     p0[n-1] = PHI( (zr-(1.-l)*hs)/l, mu) * atom[n-2];
     for (j=0;j<N;j++) p0[n-1] += w[j]/l * phi( (z[j]-(1.-l)*hs)/l, mu) * Pn[(n-2)*N+j];
   }
 }

 Free(Pn);
 Free(z);
 Free(w);
 Free(atom);

 return 0;
}


double xe1_sfm(double l, double c, double zr, double hs, int q, double mu0, double mu1, int mode, int N, int nmax, double *p0)
{ double *Smatrix, *fn, *w, *z, rn, cn, rn0, cn0, nn;
  int i, j, n, NN;
 
 c  *= sqrt( l/(2.-l) );
 zr *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );

 NN = N + 1;
 Smatrix = matrix(NN, NN);
 w       = vector(NN);
 z       = vector(NN);
 fn      = matrix(nmax, NN);

 gausslegendre(N, zr, c, z, w);

 rn = 1.; cn = 0.; rn0 = 1., cn0 = 0.;

 /* in-control, i. e. n<=q-1 */
 for (n=1;n<=q-1;n++) {
  nn = (double) n;

  /* determine r_n, n=1,2,...,q-1 */
  if ( mode==vacl ) {
    rn = sqrt( 1. - pow(1.-l, 2.*nn) );
  }

  /* determine f_n, n=1,2,...,q-1 */
  if (n==1) {
    for (i=0;i<N;i++) {
      if (mode==stat) {
        fn[0*NN+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)), mu0);
      } else {
        fn[0*NN+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l, mu0);
      }
    }
    if (mode==stat) {
      fn[0*NN+N] = PHI( (cn+rn*zr)/sqrt(l/(2.-l)), mu0);
    } else {
      fn[0*NN+N] = PHI( (cn+rn*zr-(1.-l)*hs)/l, mu0);
    }
  } else {
    for (i=0;i<N;i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * rn/l * phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*zr))/l, mu0);
      for (j=0;j<N;j++) {
        fn[(n-1)*NN+i] += w[j]*fn[(n-2)*NN+j] * rn/l
                   *phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l, mu0);
      }
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI( (cn+rn*zr-(1.-l)*(cn0+rn0*zr))/l, mu0);
    for (j=0;j<N;j++)
      fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI( (cn+rn*zr-(1.-l)*(cn0+rn0*z[j]))/l, mu0);
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;
 }

 /* out-of-control, i.e. t>=q */
 
 for (n=q;n<=nmax;n++) {
  nn = (double) n;

  /* determine r_n, n=1,2,...,q-1 */
  if ( mode==vacl ) {
    rn = sqrt( 1. - pow(1.-l, 2.*nn) );
  }

  /* determine f_n, n=q,q+1,... */  
  if ( n==1 ) {
    for (i=0;i<N;i++) {
      if (mode==stat) {
        fn[0*NN+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)), mu1);
      } else {
        fn[0*NN+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l, mu1);
      }
    }
    if ( mode==stat ) {
      fn[0*NN+N] = PHI( (cn+rn*zr)/sqrt(l/(2.-l)), mu1);
    } else {
      fn[0*NN+N] = PHI( (cn+rn*zr-(1.-l)*hs)/l, mu1);
    }
  } else {
    for (i=0;i<N;i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * rn/l * phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*zr))/l, mu1);
      for (j=0;j<N;j++) {
        fn[(n-1)*NN+i] += w[j]*fn[(n-2)*NN+j] * rn/l * phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l, mu1);
      }
      if ( n==q && q>1 ) fn[(n-1)*NN+i] /= p0[q-2];
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI( (cn+rn*zr-(1.-l)*(cn0+rn0*zr))/l, mu1);
    for (j=0;j<N;j++)
      fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI( (cn+rn*zr-(1.-l)*(cn0+rn0*z[j]))/l, mu1);
    if ( n==q && q>1 ) fn[(n-1)*NN+N] /= p0[q-2];
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;
 }

 Free(Smatrix);
 Free(w);
 Free(z);
 Free(fn);

 return 0;
}


double xe2_sf(double l, double c, double hs, double mu, int N, int nmax, double *p0)
{ double *Sm, *Pn, *w, *z;
  int i, j, n;

 c  *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );

 Sm = matrix(N, N);
 w  = vector(N);
 z  = vector(N);
 Pn = matrix(nmax, N);

 gausslegendre(N, -c, c, z, w);

 for (i=0; i<N; i++)
   for (j=0; j<N; j++) Sm[i*N+j] = w[j]/l * phi( (z[j]-(1.-l)*z[i])/l, mu);

 for (n=1; n<=nmax; n++) {
   if ( n==1 )
     for (i=0; i<N; i++) Pn[i] = PHI( (c-(1.-l)*z[i])/l, mu) - PHI( (-c-(1.-l)*z[i])/l, mu);
   else
     for (i=0; i<N; i++) {
       Pn[(n-1)*N+i] = 0.;
       for (j=0; j<N; j++) Pn[(n-1)*N+i] += Sm[i*N+j] * Pn[(n-2)*N+j];
     }

   if ( n==1 )
     p0[0] = PHI( (c-(1.-l)*hs)/l, mu) - PHI( (-c-(1.-l)*hs)/l, mu);
   else {
     p0[n-1] = 0.;
     for (j=0; j<N; j++) p0[n-1] += w[j]/l * phi( (z[j]-(1.-l)*hs)/l, mu) * Pn[(n-2)*N+j];
   }   
 }

 Free(Pn);
 Free(z);
 Free(w);
 Free(Sm);

 return 0;
}


double xte2_sf(double l, double c, double hs, int df, double mu, int N, int nmax, double *p0, int subst)
{ double *Sm, *Pn, *w, *z, norm=1., arg=0., korr=1.;
  int i, j, n;

 c  *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );

 Sm = matrix(N, N);
 w  = vector(N);
 z  = vector(N);
 Pn = matrix(nmax, N);

 switch ( subst ) {
   case IDENTITY: gausslegendre(N, -c, c, z, w);         norm = 1.; break;
   case SIN:   gausslegendre(N, -PI/2., PI/2., z, w); norm = 1.; break;
   case SINH:  gausslegendre(N, -1., 1., z, w);       norm = sinh(1.); break;
   case TAN:   gausslegendre(N, -PI/4., PI/4., z, w); norm = 1.; break;
 }     
 
 c /= norm;

 for (i=0; i<N; i++)
   for (j=0; j<N; j++) {
     switch ( subst ) {
       case IDENTITY: arg = z[j] - (1.-l)*z[i];                 korr = 1.; break;
       case SIN:   arg = c*sin(z[j]) - (1.-l)*c*sin(z[i]);   korr = c*cos(z[j]); break;
       case SINH:  arg = c*sinh(z[j]) - (1.-l)*c*sinh(z[i]); korr = c*cosh(z[j]); break;
       case TAN:   arg = c*tan(z[j]) - (1.-l)*c*tan(z[i]);   korr = c/( cos(z[j])*cos(z[j]) ); break;
     }
     Sm[i*N+j] = w[j]/l * pdf_t( arg/l - mu, df) * korr;
   }

 for (n=1; n<=nmax; n++) {
   if ( n==1 )
     for (i=0; i<N; i++) {
       switch ( subst ) {
         case IDENTITY: arg = z[i]; break;
         case SIN:   arg = c*sin(z[i]); break;
         case SINH:  arg = c*sinh(z[i]); break;
         case TAN:   arg = c*tan(z[i]); break;
       }
       Pn[i] = cdf_t( ( c*norm - (1.-l)*arg )/l - mu, df) - cdf_t( ( -c*norm - (1.-l)*arg )/l - mu, df);
     }
   else
     for (i=0; i<N; i++) {
       Pn[(n-1)*N+i] = 0.;
       for (j=0; j<N; j++) Pn[(n-1)*N+i] += Sm[i*N+j] * Pn[(n-2)*N+j];
     }

   if ( n==1 )
     p0[0] = cdf_t( ( c*norm - (1.-l)*hs )/l - mu, df) - cdf_t( ( -c*norm - (1.-l)*hs )/l - mu, df);
   else {
     p0[n-1] = 0.;
     for (j=0; j<N; j++) {
       switch ( subst ) {
         case IDENTITY: arg = z[j];         korr = 1.; break;
         case SIN:   arg = c*sin(z[j]);  korr = c*cos(z[j]); break;
         case SINH:  arg = c*sinh(z[j]); korr = c*cosh(z[j]); break;
         case TAN:   arg = c*tan(z[j]);  korr = c/( cos(z[j])*cos(z[j]) ); break;
       }
       p0[n-1] += w[j]/l * pdf_t( ( arg - (1.-l)*hs )/l - mu, df) * Pn[(n-2)*N+j] * korr;
     }
   }   
 }

 Free(Pn);
 Free(z);
 Free(w);
 Free(Sm);

 return 0;
}


double xe2_sfm(double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax, double *p0)
{ double *Smatrix, *fn, *w, *z, dn, rn, cn, rn0, cn0, delta=0., nn, fSt, aSt;
  int i, j, n;

 fSt = 0.5;
 aSt = ( -2./log10(1.-fSt) - 1.)/19.;

 c  *= sqrt( l/(2.-l) ); 
 hs *= sqrt( l/(2.-l) );
 if ( mode==fir || mode==both ) delta = 2.*hs;

 Smatrix = matrix(N, N);
 w       = vector(N);
 z       = vector(N);
 fn      = matrix(nmax, N);

 gausslegendre(N, -c, c, z, w);

 rn = 1.; cn = 0.; rn0 = 1., cn0 = 0.;

 /* in-control, i. e. n<=q-1 */

 for (n=1; n<=q-1; n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=1,2,...,q-1 */
  switch ( mode ) {
    case vacl: rn = sqrt( 1. - pow(1.-l,2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }
 
  /* determine f_n, n=1,2,...,q-1 */
  if ( n==1 ) {
    for (i=0; i<N; i++)
      if ( mode==stat )
        fn[0*N+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)),mu0);
      else
        fn[0*N+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l,mu0);
  } 
  else {
    for (i=0; i<N; i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0; j<N; j++) {
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l*phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l,mu0);
      }
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0; i<N; i++) p0[n-1] += w[i] * fn[(n-1)*N+i];

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;
 }

 /* out-of-control, i.e. t>=q */

 for (n=q; n<=nmax; n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=q,q+1,... */
  switch ( mode ) {
    case vacl: rn = sqrt( 1. - pow(1.-l,2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }

  /* determine f_n, n=q,q+1,... */
  if ( n==1 ) {
    for (i=0; i<N; i++)
      if ( mode==stat )
        fn[0*N+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)),mu1);
      else
        fn[0*N+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l,mu1);
  }
  else {
    for (i=0; i<N; i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0; j<N; j++) 
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l*phi( (cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l,mu1);
      if ( n==q && q>1 ) fn[(n-1)*N+i] /= p0[q-2];
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0; i<N; i++) p0[n-1] += w[i] * fn[(n-1)*N+i];

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;
 }

 Free(Smatrix);
 Free(w);
 Free(z);
 Free(fn);

 return 0;
}


double xte2_sfm(double l, double c, double hs, int df, int q, double mu0, double mu1, int mode, int N, int nmax, double *p0, int subst)
{ double *Smatrix, *fn, *w, *z, dn, rn, cn, rn0, cn0, delta=0., nn, fSt, aSt, norm=1., arg=0., korr=1.;
  int i, j, n;

 fSt = 0.5;
 aSt = ( -2./log10(1.-fSt) - 1.)/19.;

 c  *= sqrt( l/(2.-l) ); 
 hs *= sqrt( l/(2.-l) );
 if ( mode==fir || mode==both ) delta = 2.*hs;

 Smatrix = matrix(N, N);
 w       = vector(N);
 z       = vector(N);
 fn      = matrix(nmax, N);

 switch ( subst ) {
   case IDENTITY: gausslegendre(N, -c, c, z, w);         norm = 1.; break;
   case SIN:   gausslegendre(N, -PI/2., PI/2., z, w); norm = 1.; break;
   case SINH:  gausslegendre(N, -1., 1., z, w);       norm = sinh(1.); break;
   case TAN:   gausslegendre(N, -PI/4., PI/4., z, w); norm = 1.; break;
 }     
 
 c /= norm;

 rn = 1.; cn = 0.; rn0 = 1., cn0 = 0.;

 /* in-control, i. e. n<=q-1 */

 for (n=1; n<=q-1; n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=1,2,...,q-1 */
  switch ( mode ) {
    case vacl: rn = sqrt( 1. - pow(1.-l,2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c*norm);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c*norm);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }
 
  /* determine f_n, n=1,2,...,q-1 */
  if ( n==1 ) {
    for (i=0; i<N; i++) {
      switch ( subst ) {
        case IDENTITY: arg = z[i]; break;
        case SIN:   arg = c*sin(z[i]); break;
        case SINH:  arg = c*sinh(z[i]); break;
        case TAN:   arg = c*tan(z[i]); break;
      }
      if ( mode==stat )
        fn[0*N+i] = 1./sqrt(l/(2.-l)) * pdf_t( ( cn+rn*arg )/sqrt(l/(2.-l)) - mu0, df);
      else
        fn[0*N+i] = rn/l * pdf_t( ( cn+rn*arg - (1.-l)*hs )/l - mu0, df);
    }
  } 
  else {
    for (i=0; i<N; i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0; j<N; j++) {
	switch ( subst ) {
          case IDENTITY: arg = cn+rn*z[i] - (1.-l)*(cn0+rn0*z[j]);                 korr = 1.; break;
          case SIN:   arg = cn+rn*c*sin(z[i]) - (1.-l)*(cn0+rn0*c*sin(z[j]));   korr = c*cos(z[j]); break;
          case SINH:  arg = cn+rn*c*sinh(z[i]) - (1.-l)*(cn0+rn0*c*sinh(z[j])); korr = c*cosh(z[j]); break;
          case TAN:   arg = cn+rn*c*tan(z[i]) - (1.-l)*(cn0+rn0*c*tan(z[j]));   korr = c/( cos(z[j])*cos(z[j]) ); break;
        }
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l * pdf_t( arg/l - mu0, df) * korr;
      }
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0; i<N; i++) {
    switch ( subst ) {
      case IDENTITY: korr = 1.; break;
      case SIN:   korr = c*cos(z[i]); break;
      case SINH:  korr = c*cosh(z[i]); break;
      case TAN:   korr = c/( cos(z[i])*cos(z[i]) ); break;
    }
    p0[n-1] += w[i] * fn[(n-1)*N+i] * korr;
  }

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;
 }

 /* out-of-control, i.e. t>=q */

 for (n=q; n<=nmax; n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=q,q+1,... */
  switch ( mode ) {
    case vacl: rn = sqrt( 1. - pow(1.-l,2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c*norm);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c*norm);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }

  /* determine f_n, n=q,q+1,... */
  if ( n==1 ) {
    for (i=0; i<N; i++) {
      switch ( subst ) {
        case IDENTITY: arg = z[i]; break;
        case SIN:   arg = c*sin(z[i]); break;
        case SINH:  arg = c*sinh(z[i]); break;
        case TAN:   arg = c*tan(z[i]); break;
      }
      if ( mode==stat )
        fn[0*N+i] = 1./sqrt(l/(2.-l)) * pdf_t( ( cn+rn*arg )/sqrt(l/(2.-l)) - mu1, df);
      else
        fn[0*N+i] = rn/l * pdf_t( ( cn+rn*arg - (1.-l)*hs )/l - mu1, df);
    }
  }
  else {
    for (i=0; i<N; i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0; j<N; j++) {
	switch ( subst ) {
          case IDENTITY: arg = cn+rn*z[i] - (1.-l)*(cn0+rn0*z[j]);                 korr = 1.; break;
          case SIN:   arg = cn+rn*c*sin(z[i]) - (1.-l)*(cn0+rn0*c*sin(z[j]));   korr = c*cos(z[j]); break;
          case SINH:  arg = cn+rn*c*sinh(z[i]) - (1.-l)*(cn0+rn0*c*sinh(z[j])); korr = c*cosh(z[j]); break;
          case TAN:   arg = cn+rn*c*tan(z[i]) - (1.-l)*(cn0+rn0*c*tan(z[j]));   korr = c/( cos(z[j])*cos(z[j]) ); break;
        }
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l * pdf_t( arg/l - mu1, df) * korr;
      }
      if ( n==q && q>1 ) fn[(n-1)*N+i] /= p0[q-2];
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0;i<N;i++) {
    switch ( subst ) {
      case IDENTITY: korr = 1.; break;
      case SIN:   korr = c*cos(z[i]); break;
      case SINH:  korr = c*cosh(z[i]); break;
      case TAN:   korr = c/( cos(z[i])*cos(z[i]) ); break;
    }
    p0[n-1] += w[i] * fn[(n-1)*N+i] * korr;
  }

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;
 }

 Free(Smatrix);
 Free(w);
 Free(z);
 Free(fn);

 return 0;
}


double xe1_arlm(double l, double c, double zr, double hs, int q, double mu0, double mu1,
                int mode, int N, int nmax)
{ double *Smatrix, *p0, *fn, *w, *z, arl0, rho, rn, cn, rn0, cn0,
         arl_minus=0, arl, arl_plus=0, mn_minus, mn_plus, nn, ratio;
  int i=0, j=0, n, NN;

 c  *= sqrt( l/(2.-l) );
 zr *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );

 NN = N + 1;
 Smatrix = matrix(NN, NN);
 w       = vector(NN);
 z       = vector(NN);
 fn      = matrix(nmax, NN);
 p0      = vector(nmax);

 gausslegendre(N, zr, c, z, w);

 rn = 1.; cn = 0.; rn0 = 1., cn0 = 0.;

 /* in-control, i. e. n<=q-1 */
 for (n=1;n<=q-1;n++) {
  nn = (double) n;

  /* determine r_n, n=1,2,...,q-1 */
  if ( mode==vacl ) {
    rn = sqrt( 1. - pow(1.-l, 2.*nn) );
  }

  /* determine f_n, n=1,2,...,q-1 */
  if (n==1) {
    for (i=0;i<N;i++) {
      if (mode==stat) {
        fn[0*NN+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)), mu0);
      }
      else {
        fn[0*NN+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l, mu0);
      }
    }
    if (mode==stat) {
      fn[0*NN+N] = PHI( (cn+rn*zr)/sqrt(l/(2.-l)), mu0);
    }
    else {
      fn[0*NN+N] = PHI( (cn+rn*zr-(1.-l)*hs)/l, mu0);
    }
  }
  else {
    for (i=0;i<N;i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * rn/l * phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*zr))/l, mu0);
      for (j=0;j<N;j++) {
        fn[(n-1)*NN+i] += w[j]*fn[(n-2)*NN+j] * rn/l
                   *phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l, mu0);
      }
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI( (cn+rn*zr-(1.-l)*(cn0+rn0*zr))/l, mu0);
    for (j=0;j<N;j++)
      fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI( (cn+rn*zr-(1.-l)*(cn0+rn0*z[j]))/l, mu0);
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;
 }

 /* out-of-control, i.e. t>=q */
 arl0 = 1.; rho = 0.;

 for (n=q;n<=nmax;n++) {
  nn = (double) n;

  /* determine r_n, n=1,2,...,q-1 */
  if ( mode==vacl ) {
    rn = sqrt( 1. - pow(1.-l, 2.*nn) );
  }

  /* determine f_n, n=q,q+1,... */
  if (n==1) {
    for (i=0;i<N;i++) {
      if (mode==stat) {
        fn[0*NN+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)), mu1);
      }
      else {
        fn[0*NN+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l, mu1);
      }
    }
    if (mode==stat) {
      fn[0*NN+N] = PHI( (cn+rn*zr)/sqrt(l/(2.-l)), mu1);
    }
    else {
      fn[0*NN+N] = PHI( (cn+rn*zr-(1.-l)*hs)/l, mu1);
    }
  }
  else {
    for (i=0;i<N;i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * rn/l * phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*zr))/l, mu1);
      for (j=0;j<N;j++) {
        fn[(n-1)*NN+i] += w[j]*fn[(n-2)*NN+j] * rn/l * phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l, mu1);
      }
      if (n==q && q>1) fn[(n-1)*NN+i] /= p0[q-2];
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI( (cn+rn*zr-(1.-l)*(cn0+rn0*zr))/l, mu1);
    for (j=0;j<N;j++)
      fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI( (cn+rn*zr-(1.-l)*(cn0+rn0*z[j]))/l, mu1);
    if (n==q && q>1) fn[(n-1)*NN+N] /= p0[q-2];
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;

  /* computation of m_n+1^- and m_n+1^+, n=m-1,m,... */
  mn_minus = 1.; mn_plus = 0.;
  if (n>q) {
   for (i=0;i<NN;i++) {
    if (fn[(n-2)*NN+i]==0)
      if (fn[(n-1)*NN+i]==0) ratio = 0.; else ratio = 1.;
    else ratio = fn[(n-1)*NN+i]/fn[(n-2)*NN+i];
    if ( ratio<mn_minus ) mn_minus = ratio;
    if ( ratio>mn_plus ) mn_plus = ratio;
   }
  }

  if (n>q) rho = p0[n-1]/p0[n-2];

  /* computation of ARL, ARL^-, and ARL^+ */
  arl = arl0 + p0[n-1]/(1.-rho);
  if (mn_minus<1.) arl_minus = arl0 + p0[n-1]/(1.-mn_minus);
  else             arl_minus = -1.;
  if (mn_plus<1.) arl_plus = arl0 + p0[n-1]/(1.-mn_plus);
  else            arl_plus = -1.;
  arl0 += p0[n-1];

  if ( fabs((arl_plus-arl_minus)) < 1e-7 ) n = nmax+1;
 }

 arl = (arl_plus+arl_minus)/2; rho0 = rho;

 Free(Smatrix);
 Free(w);
 Free(z);
 Free(fn);
 Free(p0);

 return arl;
}


double xe1_arlm_hom(double l, double c, double zr, double hs, int q, double mu0, double mu1, int N, double *ced)
{ double *fn, *w, *z, *a, *arl, norm;
  int i, j, n, NN;
  
  NN = N + 1;
  w   = vector(NN);
  z   = vector(NN);
  fn  = matrix(q+1, NN);
  a   = matrix(NN,NN);
  arl = vector(NN);
  
  c  *= sqrt( l/(2.-l) );
  zr *= sqrt( l/(2.-l) );
  hs *= sqrt( l/(2.-l) );

  gausslegendre(N, zr, c, z, w);

  /* ARL vector */
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) a[i*NN+j] = -w[j]/l * phi( (z[j]-(1.-l)*z[i])/l, mu1);
    ++a[i*NN+i];
    a[i*NN+N] = - PHI( (zr-(1.-l)*z[i])/l, mu1);
  }
  for (j=0; j<N; j++) a[N*NN+j] = -w[j]/l * phi( (z[j]-(1.-l)*zr)/l, mu1);
  a[N*NN+N] = 1. - PHI( zr, mu1);

  for (j=0; j<NN; j++) arl[j] = 1.;
  LU_solve(a, arl, NN);
  
  /* q == 1 */
  ced[0] = 1. + PHI( (zr-(1.-l)*hs)/l, mu1) * arl[N];
  for (j=0; j<N; j++) ced[0] +=  w[j]/l * phi( (z[j]-(1.-l)*hs)/l, mu1) * arl[j];

/* density sequence for q > 1 */
  for (n=1; n<=q-1; n++) {
    if (n==1) {      
      for (i=0;i<N;i++) fn[0*NN+i] = phi( (z[i]-(1.-l)*hs)/l, mu0)/l;
      fn[0*NN+N] = PHI( (zr-(1.-l)*hs)/l, mu0);  
    } else {      
      for (i=0; i<N; i++) {
        fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * phi( (z[i]-(1.-l)*zr)/l, mu0)/l;
        for (j=0; j<N; j++) fn[(n-1)*NN+i] += w[j] * fn[(n-2)*NN+j] * phi( (z[i]-(1.-l)*z[j])/l, mu0)/l;
      }
      fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI( zr, mu0);
      for (j=0; j<N; j++) fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI( (zr-(1.-l)*z[j])/l, mu0);      
    }
  
    ced[n] = fn[(n-1)*NN+N] * arl[N];
    norm = fn[(n-1)*NN+N];
    for (j=0; j<N; j++) {
      ced[n] += w[j] * fn[(n-1)*NN+j] * arl[j]; 
      norm += w[j] * fn[(n-1)*NN+j];
    }
    ced[n] /= norm;
  }  

  Free(w);
  Free(z);
  Free(fn);
  Free(a);
  Free(arl);

  return 0;
}


double xlimit1_arlm(double c, double zr, int q, double mu0, double mu1, int N, int nmax)
{ double *Smatrix, *p0, *fn, *w, *z, l1, l2, arl0, rho,
         arl_minus=0, arl, arl_plus=0, mn_minus, mn_plus, nn, ratio;
  int i=0, j=0, n, NN;

 NN = N + 1;
 Smatrix = matrix(NN, NN);
 w       = vector(NN);
 z       = vector(NN);
 fn      = matrix(nmax, NN);
 p0      = vector(nmax);

 gausslegendre(N, zr, c, z, w);

 /* in-control, i. e. n<=q-1 */
 for (n=1;n<=q-1;n++) {
  nn = (double) n;
  /* determine f_n, n=1,2,...,q-1 */
  if (n==1) {
    for (i=0;i<N;i++) fn[0*NN+i] = phi( z[i], mu0);
    fn[0*NN+N] = PHI( zr, mu0);
  }
  else {
    l1 = sqrt( (nn-1.)/nn );
    l2 = sqrt( 1./nn );
    for (i=0;i<N;i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * phi( (z[i]-l1*zr)/l2, mu0)/l2;
      for (j=0;j<N;j++) {
        fn[(n-1)*NN+i] += w[j]*fn[(n-2)*NN+j] * phi( (z[i]-l1*z[j])/l2, mu0)/l2;
      }
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI( (zr-l1*zr)/l2, mu0);
    for (j=0;j<N;j++)
      fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI( (zr-l1*z[j])/l2, mu0);
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];
 }

 /* out-of-control, i.e. t>=q */
 arl0 = 1.; rho = 0.;

 for (n=q;n<=nmax;n++) {
  nn = (double) n;

  /* determine f_n, n=q,q+1,... */
  if (n==1) {
    for (i=0;i<N;i++) fn[0*NN+i] = phi( z[i], mu1);
    fn[0*NN+N] = PHI( zr, mu1);
  }
  else {
    l1 = sqrt( (nn-1.)/nn );
    l2 = sqrt( 1./nn );
    for (i=0;i<N;i++) {
      fn[(n-1)*NN+i] = fn[(n-2)*NN+N] * phi( (z[i]-l1*zr)/l2, mu1)/l2;
      for (j=0;j<N;j++) {
        fn[(n-1)*NN+i] += w[j] * fn[(n-2)*NN+j] * phi( (z[i]-l1*z[j])/l2, mu1)/l2;
      }
      if (n==q && q>1) fn[(n-1)*NN+i] /= p0[q-2];
    }
    fn[(n-1)*NN+N] = fn[(n-2)*NN+N] * PHI( (zr-l1*zr)/l2, mu1);
    for (j=0;j<N;j++)
      fn[(n-1)*NN+N] += w[j] * fn[(n-2)*NN+j] * PHI( (zr-l1*z[j])/l2, mu1);
    if (n==q && q>1) fn[(n-1)*NN+N] /= p0[q-2];
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = fn[(n-1)*NN+N];
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*NN+i];

  /* computation of m_n+1^- and m_n+1^+, n=m-1,m,... */
  mn_minus = 1.; mn_plus = 0.;
  if (n>q) {
   for (i=0;i<NN;i++) {
    if (fn[(n-2)*NN+i]==0)
      if (fn[(n-1)*NN+i]==0) ratio = 0.; else ratio = 1.;
    else ratio = fn[(n-1)*NN+i]/fn[(n-2)*NN+i];
    if ( ratio<mn_minus ) mn_minus = ratio;
    if ( ratio>mn_plus ) mn_plus = ratio;
   }
  }

  if (n>q) rho = p0[n-1]/p0[n-2];

  /* computation of ARL, ARL^-, and ARL^+ */
  arl = arl0 + p0[n-1]/(1.-rho);
  if (mn_minus<1.) arl_minus = arl0 + p0[n-1]/(1.-mn_minus);
  else             arl_minus = -1.;
  if (mn_plus<1.) arl_plus = arl0 + p0[n-1]/(1.-mn_plus);
  else            arl_plus = -1.;
  arl0 += p0[n-1];

  if ( (p0[n-1]>p0[n-2] || rho>1.) && n>10 ) error("invalid ARL value");
  if ( fabs((arl_plus-arl_minus)) < 1e-5 ) n = nmax+1;
 }

 arl = (arl_plus+arl_minus)/2; rho0 = rho;

 Free(Smatrix);
 Free(w);
 Free(z);
 Free(fn);
 Free(p0);

 return arl;
}


double xe2_arlm(double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax)
{ double *Smatrix, *p0, *fn, *w, *z,
         arl0, rho, dn, rn, cn, rn0, cn0, delta=0.,
         arl_minus=0, arl, arl_plus=0, mn_minus, mn_plus, nn,
         fSt, aSt, ratio;
  int i, j, n;

 fSt = 0.5;
 aSt = ( -2./log10(1.-fSt) - 1.)/19.;

 c  *= sqrt( l/(2.-l) ); 
 hs *= sqrt( l/(2.-l) );
 if (mode==fir || mode==both) delta = 2.*hs;

 Smatrix = matrix(N,N);
 w       = vector(N);
 z       = vector(N);
 fn      = matrix(nmax,N);
 p0      = vector(nmax);

 gausslegendre(N,-c,c,z,w);

 rn = 1.; cn = 0.; rn0 = 1., cn0 = 0.;

 /* in-control, i. e. n<=q-1 */
 for (n=1;n<=q-1;n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=1,2,...,q-1 */
  switch (mode) {
    case vacl: rn = sqrt( 1. - pow(1.-l, 2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }
 
  /* determine f_n, n=1,2,...,q-1 */
  if (n==1) {
    for (i=0;i<N;i++)
      if (mode==stat)
        fn[0*N+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)),mu0);
      else
        fn[0*N+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l,mu0);
  } 
  else {
    for (i=0;i<N;i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0;j<N;j++) {
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l
                   *phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l,mu0);
      }
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*N+i];

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;
 }

 /* out-of-control, i.e. t>=q */
 arl0 = 1.; rho = 0.;

 for (n=q;n<=nmax;n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=q,q+1,... */
  switch (mode) {
    case vacl: rn = sqrt( 1. - pow(1.-l,2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }

  /* determine f_n, n=q,q+1,... */
  if (n==1) {
    for (i=0;i<N;i++)
      if (mode==stat)
        fn[0*N+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)),mu1);
      else
        fn[0*N+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l,mu1);
  }
  else {
    for (i=0;i<N;i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0;j<N;j++) 
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l
                   *phi( (cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l,mu1);
      if (n==q && q>1) fn[(n-1)*N+i] /= p0[q-2];
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*N+i];

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;

  /* computation of m_n+1^- and m_n+1^+, n=m-1,m,... */
  mn_minus = 1.; mn_plus = 0.;
  if (n>q) {
   for (i=0;i<N;i++) {
    if (fn[(n-2)*N+i]==0)
     if (fn[(n-1)*N+i]==0) ratio = 0.; else ratio = 1.;
    else ratio = fn[(n-1)*N+i]/fn[(n-2)*N+i];
    if ( ratio<mn_minus ) mn_minus = ratio;
    if ( ratio>mn_plus ) mn_plus = ratio;
   }
  }
 
  if (n>q) rho = p0[n-1]/p0[n-2];

  /* computation of ARL, ARL^-, and ARL^+ */
  arl = arl0 + p0[n-1]/(1.-rho);
  if (mn_minus<1.) arl_minus = arl0 + p0[n-1]/(1.-mn_minus);
  else             arl_minus = -1.;
  if (mn_plus<1.) arl_plus = arl0 + p0[n-1]/(1.-mn_plus);
  else            arl_plus = -1.;
  arl0 += p0[n-1]; 

  if ( fabs((arl_plus-arl_minus)) < 1e-7 ) n = nmax+1;
 }

 arl = (arl_plus+arl_minus)/2; rho0 = rho;

 Free(Smatrix);
 Free(w);
 Free(z);
 Free(fn);
 Free(p0);

 return arl;
}


double xe2_arlmc(double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax)
{ double *Smatrix, *p0, *Psi, *fn, *w, *z,
         arl0, rho, dn, rn, cn, rn0, cn0, delta=0.,
         arl_minus=0, arl, arl_plus=0, mn_minus, mn_plus, nn,
         fSt, aSt, ratio;
  int i, j, n;

 fSt = 0.5;
 aSt = ( -2./log10(1.-fSt) - 1.)/19.;

 c  *= sqrt( l/(2.-l) ); 
 hs *= sqrt( l/(2.-l) );
 if (mode==fir || mode==both) delta = 2.*hs;

 Smatrix = matrix(N,N);
 w       = vector(N);
 z       = vector(N);
 fn      = matrix(nmax,N);
 p0      = vector(nmax);
 Psi     = vector(nmax);

 gausslegendre(N,-c,c,z,w);

 rn = 1.; cn = 0.; rn0 = 1., cn0 = 0.;

 /* in-control, i. e. n<=q-1 */
 for (n=1;n<=q-1;n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=1,2,...,q-1 */
  switch (mode) {
    case vacl: rn = sqrt( 1. - pow(1.-l, 2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }
 
  /* determine f_n, n=1,2,...,q-1 */
  if (n==1) {
    for (i=0;i<N;i++)
      if (mode==stat)
        fn[0*N+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)),mu0);
      else
        fn[0*N+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l,mu0);
    Psi[0] = PHI( (cn-rn*c-(1.-l)*hs)/l, mu0) + 1. - PHI( (cn+rn*c-(1.-l)*hs)/l, mu0);  
  } 
  else {
    for (i=0;i<N;i++) {
      fn[(n-1)*N+i] = Psi[n-1] * rn/l*phi( (cn+rn*z[i]-(1.-l)*hs)/l, mu0);
      for (j=0;j<N;j++) {
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j] * rn/l*phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l,mu0);
      }
      Psi[n] = Psi[n-1] * ( PHI( (cn-rn*c-(1.-l)*hs)/l, mu0) + 1. - PHI( (cn+rn*c-(1.-l)*hs)/l, mu0) );
      for (j=0;j<N;j++) {
        Psi[n] += w[j] * fn[(n-2)*N+j] * ( PHI( (cn-rn*c-(1.-l)*(cn0+rn0*z[j]))/l, mu0) + 1. - PHI( (cn+rn*c-(1.-l)*(cn0+rn0*z[j]))/l, mu0) );
      }
    }
  }

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;
 }

 /* out-of-control, i.e. t>=q */
 arl0 = 1.; rho = 0.;

 for (n=q;n<=nmax;n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=q,q+1,... */
  switch (mode) {
    case vacl: rn = sqrt( 1. - pow(1.-l,2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }

  /* determine f_n, n=q,q+1,... */
  if (n==1) {
    for (i=0;i<N;i++)
      if (mode==stat)
        fn[0*N+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)),mu1);
      else
        fn[0*N+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l,mu1);
  }
  else {
    for (i=0;i<N;i++) {
      fn[(n-1)*N+i] = 0.;
      if (n==q && q>1) fn[(n-1)*N+i] = Psi[n-1] * rn/l*phi( (cn+rn*z[i]-(1.-l)*hs)/l, mu1);
      for (j=0;j<N;j++) 
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j] * rn/l*phi( (cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l, mu1);
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*N+i];

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;

  /* computation of m_n+1^- and m_n+1^+, n=m-1,m,... */
  mn_minus = 1.; mn_plus = 0.;
  if (n>q) {
   for (i=0;i<N;i++) {
    if (fn[(n-2)*N+i]==0)
     if (fn[(n-1)*N+i]==0) ratio = 0.; else ratio = 1.;
    else ratio = fn[(n-1)*N+i]/fn[(n-2)*N+i];
    if ( ratio<mn_minus ) mn_minus = ratio;
    if ( ratio>mn_plus ) mn_plus = ratio;
   }
  }
 
  if (n>q) rho = p0[n-1]/p0[n-2];

  /* computation of ARL, ARL^-, and ARL^+ */
  arl = arl0 + p0[n-1]/(1.-rho);
  if (mn_minus<1.) arl_minus = arl0 + p0[n-1]/(1.-mn_minus);
  else             arl_minus = -1.;
  if (mn_plus<1.) arl_plus = arl0 + p0[n-1]/(1.-mn_plus);
  else            arl_plus = -1.;
  arl0 += p0[n-1]; 

  if ( fabs((arl_plus-arl_minus)) < 1e-7 ) n = nmax+1;
 }

 arl = (arl_plus+arl_minus)/2; rho0 = rho; 

 Free(Smatrix);
 Free(w);
 Free(z);
 Free(fn);
 Free(p0);

 return arl;
}


int xe2_arlm_special(double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax, double *pair)
{ double *Smatrix, *p0, *fn, *w, *z,
         arl0, rho, dn, rn, cn, rn0, cn0, delta=0.,
         arl_minus=0, arl, arl_plus=0, mn_minus, mn_plus, nn,
         fSt, aSt, ratio;
  int i, j, n;

 fSt = 0.5;
 aSt = ( -2./log10(1.-fSt) - 1.)/19.;

 c  *= sqrt( l/(2.-l) ); 
 hs *= sqrt( l/(2.-l) );
 if (mode==fir || mode==both) delta = 2.*hs;

 Smatrix = matrix(N,N);
 w       = vector(N);
 z       = vector(N);
 fn      = matrix(nmax,N);
 p0      = vector(nmax);

 gausslegendre(N,-c,c,z,w);

 rn = 1.; cn = 0.; rn0 = 1., cn0 = 0.;

 /* in-control, i. e. n<=q-1 */
 for (n=1;n<=q-1;n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=1,2,...,q-1 */
  switch (mode) {
    case vacl: rn = sqrt( 1. - pow(1.-l, 2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }
 
  /* determine f_n, n=1,2,...,q-1 */
  if (n==1) {
    for (i=0;i<N;i++)
      if (mode==stat)
        fn[0*N+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)),mu0);
      else
        fn[0*N+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l,mu0);
  } 
  else {
    for (i=0;i<N;i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0;j<N;j++) {
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l
                   *phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l,mu0);
      }
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*N+i];

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;
 }

 /* out-of-control, i.e. t>=q */
 arl0 = 1.; rho = 0.;

 for (n=q;n<=nmax;n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=q,q+1,... */
  switch (mode) {
    case vacl: rn = sqrt( 1. - pow(1.-l,2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }

  /* determine f_n, n=q,q+1,... */
  if (n==1) {
    for (i=0;i<N;i++)
      if (mode==stat)
        fn[0*N+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)),mu1);
      else
        fn[0*N+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l,mu1);
  }
  else {
    for (i=0;i<N;i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0;j<N;j++) 
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l
                   *phi( (cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l,mu1);
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0;i<N;i++) p0[n-1] += w[i] * fn[(n-1)*N+i];

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;

  /* computation of m_n+1^- and m_n+1^+, n=m-1,m,... */
  mn_minus = 1.; mn_plus = 0.;
  if (n>q) {
   for (i=0;i<N;i++) {
    if (fn[(n-2)*N+i]==0)
     if (fn[(n-1)*N+i]==0) ratio = 0.; else ratio = 1.;
    else ratio = fn[(n-1)*N+i]/fn[(n-2)*N+i];
    if ( ratio<mn_minus ) mn_minus = ratio;
    if ( ratio>mn_plus ) mn_plus = ratio;
   }
  }
 
  if ( n > q ) rho = p0[n-1]/p0[n-2];

  /* computation of ARL, ARL^-, and ARL^+ */
  arl = arl0 + p0[n-1]/(1.-rho);
  if (mn_minus<1.) arl_minus = arl0 + p0[n-1]/(1.-mn_minus);
  else             arl_minus = -1.;
  if (mn_plus<1.) arl_plus = arl0 + p0[n-1]/(1.-mn_plus);
  else            arl_plus = -1.;
  arl0 += p0[n-1]; 

  if ( fabs((arl_plus-arl_minus)) < 1e-7 ) n = nmax+1;
 }

 arl = ( arl_plus + arl_minus )/2;
 
 pair[0] = 1.;
 if ( q > 1 ) pair[0] = p0[q-2];
 pair[1] = arl;

 Free(Smatrix);
 Free(w);
 Free(z);
 Free(fn);
 Free(p0);

 return 0;
}


double xte2_arlm(double l, double c, double hs, int df, int q, double mu0, double mu1, int mode, int N, int nmax, int subst)
{ double *Smatrix, *p0, *fn, *w, *z,
         arl0, rho, dn, rn, cn, rn0, cn0, delta=0.,
         arl_minus=0, arl, arl_plus=0, mn_minus, mn_plus, nn,
         fSt, aSt, ratio, norm=1., arg=0., korr=1.;
  int i, j, n;

 fSt = 0.5;
 aSt = ( -2./log10(1.-fSt) - 1.)/19.;

 c  *= sqrt( l/(2.-l) ); 
 hs *= sqrt( l/(2.-l) );
 if (mode==fir || mode==both) delta = 2.*hs;

 Smatrix = matrix(N,N);
 w       = vector(N);
 z       = vector(N);
 fn      = matrix(nmax,N);
 p0      = vector(nmax);

 switch ( subst ) {
   case IDENTITY: gausslegendre(N, -c, c, z, w);         norm = 1.; break;
   case SIN:   gausslegendre(N, -PI/2., PI/2., z, w); norm = 1.; break;
   case SINH:  gausslegendre(N, -1., 1., z, w);       norm = sinh(1.); break;
   case TAN:   gausslegendre(N, -PI/4., PI/4., z, w); norm = 1.; break;
 }     
 
 c /= norm;

 rn = 1.; cn = 0.; rn0 = 1., cn0 = 0.;

 /* in-control, i. e. n<=q-1 */
 for (n=1; n<=q-1; n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=1,2,...,q-1 */
  switch ( mode ) {
    case vacl: rn = sqrt( 1. - pow(1.-l, 2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c*norm);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c*norm);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }
 
  /* determine f_n, n=1,2,...,q-1 */
  if ( n==1 ) {
    for (i=0; i<N; i++) {
      switch ( subst ) {
        case IDENTITY: arg = z[i]; break;
        case SIN:   arg = c*sin(z[i]); break;
        case SINH:  arg = c*sinh(z[i]); break;
        case TAN:   arg = c*tan(z[i]); break;
      }
      if ( mode==stat )
        fn[0*N+i] = 1./sqrt(l/(2.-l)) * pdf_t( ( cn+rn*arg )/sqrt(l/(2.-l)) - mu0, df);
      else
        fn[0*N+i] = rn/l * pdf_t( ( cn+rn*arg - (1.-l)*hs )/l - mu0, df);
    }
  } 
  else {
    for (i=0;i<N;i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0;j<N;j++) {
	switch ( subst ) {
          case IDENTITY: arg = cn+rn*z[i] - (1.-l)*(cn0+rn0*z[j]);                 korr = 1.; break;
          case SIN:   arg = cn+rn*c*sin(z[i]) - (1.-l)*(cn0+rn0*c*sin(z[j]));   korr = c*cos(z[j]); break;
          case SINH:  arg = cn+rn*c*sinh(z[i]) - (1.-l)*(cn0+rn0*c*sinh(z[j])); korr = c*cosh(z[j]); break;
          case TAN:   arg = cn+rn*c*tan(z[i]) - (1.-l)*(cn0+rn0*c*tan(z[j]));   korr = c/( cos(z[j])*cos(z[j]) ); break;
        }
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l * pdf_t( arg/l - mu0, df) * korr;
      }
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0; i<N; i++) {
    switch ( subst ) {
      case IDENTITY: korr = 1.; break;
      case SIN:   korr = c*cos(z[i]); break;
      case SINH:  korr = c*cosh(z[i]); break;
      case TAN:   korr = c/( cos(z[i])*cos(z[i]) ); break;
    }
    p0[n-1] += w[i] * fn[(n-1)*N+i] * korr;
  }

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;
 }

 /* out-of-control, i.e. t>=q */
 arl0 = 1.; rho = 0.;

 for (n=q; n<=nmax; n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=q,q+1,... */
  switch ( mode ) {
    case vacl: rn = sqrt( 1. - pow(1.-l,2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c*norm);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c*norm);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }

  /* determine f_n, n=q,q+1,... */
  if ( n==1 ) {
    for (i=0; i<N; i++) {
      switch ( subst ) {
        case IDENTITY: arg = z[i]; break;
        case SIN:   arg = c*sin(z[i]); break;
        case SINH:  arg = c*sinh(z[i]); break;
        case TAN:   arg = c*tan(z[i]); break;
      }
      if ( mode==stat )
        fn[0*N+i] = 1./sqrt(l/(2.-l)) * pdf_t( ( cn+rn*arg )/sqrt(l/(2.-l)) - mu1, df);
      else
        fn[0*N+i] = rn/l * pdf_t( ( cn+rn*arg - (1.-l)*hs )/l - mu1, df);
    }
  }
  else {
    for (i=0; i<N; i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0; j<N; j++) {
	switch ( subst ) {
          case IDENTITY: arg = cn+rn*z[i] - (1.-l)*(cn0+rn0*z[j]);                 korr = 1.; break;
          case SIN:   arg = cn+rn*c*sin(z[i]) - (1.-l)*(cn0+rn0*c*sin(z[j]));   korr = c*cos(z[j]); break;
          case SINH:  arg = cn+rn*c*sinh(z[i]) - (1.-l)*(cn0+rn0*c*sinh(z[j])); korr = c*cosh(z[j]); break;
          case TAN:   arg = cn+rn*c*tan(z[i]) - (1.-l)*(cn0+rn0*c*tan(z[j]));   korr = c/( cos(z[j])*cos(z[j]) ); break;
        }
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l * pdf_t( arg/l - mu1, df) * korr;
      }
      if ( n==q && q>1 ) fn[(n-1)*N+i] /= p0[q-2];
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0;i<N;i++) {
    switch ( subst ) {
      case IDENTITY: korr = 1.; break;
      case SIN:   korr = c*cos(z[i]); break;
      case SINH:  korr = c*cosh(z[i]); break;
      case TAN:   korr = c/( cos(z[i])*cos(z[i]) ); break;
    }
    p0[n-1] += w[i] * fn[(n-1)*N+i] * korr;
  }

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;

  /* computation of m_n+1^- and m_n+1^+, n=m-1,m,... */
  mn_minus = 1.; mn_plus = 0.;
  if ( n>q ) {
   for (i=0; i<N; i++) {
    if (fn[(n-2)*N+i]==0)
     if (fn[(n-1)*N+i]==0) ratio = 0.; else ratio = 1.;
    else ratio = fn[(n-1)*N+i]/fn[(n-2)*N+i];
    if ( ratio<mn_minus ) mn_minus = ratio;
    if ( ratio>mn_plus ) mn_plus = ratio;
   }
  }
 
  if ( n>q ) rho = p0[n-1]/p0[n-2];

  /* computation of ARL, ARL^-, and ARL^+ */
  arl = arl0 + p0[n-1]/(1.-rho);
  if ( mn_minus<1. ) arl_minus = arl0 + p0[n-1]/(1.-mn_minus);
  else               arl_minus = -1.;
  if ( mn_plus<1. ) arl_plus = arl0 + p0[n-1]/(1.-mn_plus);
  else              arl_plus = -1.;
  arl0 += p0[n-1]; 

  if ( fabs((arl_plus-arl_minus)) < 1e-7 ) n = nmax+1;
 }

 arl = (arl_plus+arl_minus)/2; rho0 = rho;

 Free(Smatrix);
 Free(w);
 Free(z);
 Free(fn);
 Free(p0);

 return arl;
}


double xe2_arlm_hom(double l, double c, double hs, int q, double mu0, double mu1, int N, double *ced)
{ double *fn, *w, *z, *a, *arl, norm;
  int i, j, n;
  
  w   = vector(N);
  z   = vector(N);
  fn  = matrix(q+1, N);
  a   = matrix(N,N);
  arl = vector(N);
  
  c  *= sqrt( l/(2.-l) ); 
  hs *= sqrt( l/(2.-l) );

  gausslegendre(N, -c, c, z, w);

  for (i=0; i<N; i++) {
    for (j=0;j<N;j++) a[i*N+j] = -w[j]/l * phi( (z[j]-(1.-l)*z[i])/l, mu1);
    ++a[i*N+i];
  }

  for (j=0; j<N; j++) arl[j] = 1.;
  LU_solve(a, arl, N);
 
  /* q == 1 */
  ced[0] = 1.;
  for (j=0; j<N; j++) ced[0] +=  w[j]/l * phi( (z[j]-(1.-l)*hs)/l, mu1) * arl[j];
  
  /* density sequence for q > 1 */
  for (n=1; n<=q-1; n++) {    
    if ( n==1 ) {
      for (i=0; i<N; i++) fn[0*N+i] = phi( (z[i]-(1.-l)*hs)/l, mu0)/l;
    } else {
      for (i=0; i<N; i++) {
        fn[(n-1)*N+i] = 0.;
        for (j=0; j<N; j++) fn[(n-1)*N+i] += w[j] * fn[(n-2)*N+j] * phi((z[i]-(1.-l)*z[j])/l, mu0)/l;
      }
    }
 
    ced[n] = 0.;
    norm = 0.;
    for (j=0; j<N; j++) {
      ced[n] += w[j] * fn[(n-1)*N+j] * arl[j]; 
      norm += w[j] * fn[(n-1)*N+j];
    }
    ced[n] /= norm;
  }
 
  Free(w);
  Free(z);
  Free(fn);
  Free(a);
  Free(arl);

  return 0;
}


double xte2_arlm_hom(double l, double c, double hs, int df, int q, double mu0, double mu1, int N, double *ced, int subst)
{ double *fn, *w, *z, *a, *arl, nenner=1., norm=1., arg=0., korr=1.;
  int i, j, n;
  
  w   = vector(N);
  z   = vector(N);
  fn  = matrix(q+1, N);
  a   = matrix(N,N);
  arl = vector(N);
  
  c  *= sqrt( l/(2.-l) ); 
  hs *= sqrt( l/(2.-l) );

  switch ( subst ) {
    case IDENTITY: gausslegendre(N, -c, c, z, w); norm = 1.; break;
    case SIN:   gausslegendre(N, -PI/2., PI/2., z, w); norm = 1.; break;
    case SINH:  gausslegendre(N, -1., 1., z, w); norm = sinh(1.); break;
    case TAN:   gausslegendre(N, -PI/4., PI/4., z, w); norm = 1.; break;
  }     
 
  c /= norm;

  for (i=0;i<N;i++) {
    for (j=0;j<N;j++) {
      switch ( subst ) {
        case IDENTITY: arg = z[j] - (1.-l)*z[i]; korr = 1.; break;
        case SIN:   arg = c*sin(z[j]) - (1.-l)*c*sin(z[i]); korr = c*cos(z[j]); break;
        case SINH:  arg = c*sinh(z[j]) - (1.-l)*c*sinh(z[i]); korr = c*cosh(z[j]); break;
        case TAN:   arg = c*tan(z[j]) - (1.-l)*c*tan(z[i]); korr = c/( cos(z[j])*cos(z[j]) ); break;
      }
      a[i*N+j] = -w[j]/l * pdf_t( arg/l - mu1, df) * korr;
    }
    ++a[i*N+i];
  }

  for (j=0;j<N;j++) arl[j] = 1.;
  LU_solve(a,arl,N);
 
  /* q == 1 */
  ced[0] = 1.;
  for (j=0; j<N; j++) {
    switch ( subst ) {
      case IDENTITY: arg = z[j] - (1.-l)*hs; korr = 1.; break;
      case SIN:   arg = c*sin(z[j]) - (1.-l)*hs; korr = c*cos(z[j]); break;
      case SINH:  arg = c*sinh(z[j]) - (1.-l)*hs; korr = c*cosh(z[j]); break;
      case TAN:   arg = c*tan(z[j]) - (1.-l)*hs; korr = c/( cos(z[j])*cos(z[j]) ); break;
    }
    ced[0] += w[j]/l * pdf_t( arg/l - mu1, df) * arl[j] * korr;
  }
  
  /* density sequence for q > 1 */
  for (n=1; n<=q-1; n++) {    
    if ( n==1 ) {
      for (i=0; i<N; i++) {
	switch ( subst ) {
	  case IDENTITY: arg = z[i]; break;
          case SIN:   arg = c*sin(z[i]); break;
          case SINH:  arg = c*sinh(z[i]); break;
          case TAN:   arg = c*tan(z[i]); break;
        }
	fn[0*N+i] = pdf_t( ( arg - (1.-l)*hs )/l - mu0, df)/l;
      }
    } else {
      for (i=0; i<N; i++) {
        fn[(n-1)*N+i] = 0.;
        for (j=0; j<N; j++) {
	  switch ( subst ) {
            case IDENTITY: arg = z[i] - (1.-l)*z[j];                 korr = 1.; break;
            case SIN:   arg = c*sin(z[i]) - (1.-l)*c*sin(z[j]);   korr = c*cos(z[j]); break;
            case SINH:  arg = c*sinh(z[i]) - (1.-l)*c*sinh(z[j]); korr = c*cosh(z[j]); break;
            case TAN:   arg = c*tan(z[i]) - (1.-l)*c*tan(z[j]);   korr = c/( cos(z[j])*cos(z[j]) ); break;
          }
	  fn[(n-1)*N+i] += w[j] * fn[(n-2)*N+j] * pdf_t( arg/l - mu0, df)/l * korr;
	}
      }
    }
 
    ced[n] = 0.;
    nenner = 0.;
    for (j=0; j<N; j++) {
      switch ( subst ) {
        case IDENTITY: korr = 1.; break;
        case SIN:   korr = c*cos(z[j]); break;
        case SINH:  korr = c*cosh(z[j]); break;
        case TAN:   korr = c/( cos(z[j])*cos(z[j]) ); break;
      }
      ced[n] += w[j] * fn[(n-1)*N+j] * arl[j] * korr; 
      nenner += w[j] * fn[(n-1)*N+j] * korr;
    }
    ced[n] /= nenner;
  }
 
  Free(w);
  Free(z);
  Free(fn);
  Free(a);
  Free(arl);

  return 0;
}


int qm_for_l_and_c(double l, double c) {
  int qm=20;
  qm = (int)ceil( 3.141 * c / sqrt(l) );
  if ( qm < 20 ) qm = 20;
  /*if ( qm > 1000 ) qm = 1000;*/
  return qm;
}


/* routines for prerun impact on average ARL, QRL performance */

/* 1. ARL (fixed limits) */

double xe2_iglarl_prerun_MU(double l, double c, double hs, double mu, int pn, int qm, double truncate)
{ double *w, *z, b, result, dn, sdn;
  int i, Nlocal;

  w = vector(qm);
  z = vector(qm);
  dn = (double)pn;
  sdn = sqrt(dn);
  b = -qPHI(truncate/2.)/sdn;
  gausslegendre(qm, -b, b, z, w);
  Nlocal = qm_for_l_and_c(l, c);
  result = 0.;
  for (i=0; i<qm; i++) result += w[i] * sdn*phi( z[i]*sdn, 0. ) * xe2_iglarl(l, c, hs, z[i]+mu, Nlocal);
  Free(w);
  Free(z);
  return result;
}


double xe2_iglarl_prerun_SIGMA(double l, double c, double hs, double mu, int pn, int qm, double truncate)
{ double *w, *z, b1, b2, result, ddf;
  int i, Nlocal;
  
  w = vector(qm);
  z = vector(qm);
  ddf = (double)(pn-1);  
  b1 = sqrt(qCHI(     truncate/2., pn-1)/ddf);
  b2 = sqrt(qCHI(1. - truncate/2., pn-1)/ddf); 
  gausslegendre(qm, b1, b2, z, w);  
  result = 0.;
  for (i=0; i<qm; i++) {
    Nlocal = qm_for_l_and_c(l, z[i]*c);
    result += w[i] * 2.*ddf*z[i]*chi( ddf*z[i]*z[i], pn-1) * xe2_iglarl(l, z[i]*c, hs, mu, Nlocal);
  }
  Free(w);
  Free(z);
  return result;
}


double xe2_iglarl_prerun_BOTH(double l, double c, double hs, double mu, int pn, int df, int qm1, int qm2, double truncate)
{ double *w1, *z1, *w2, *z2, b, b1, b2, result, dn, sdn, ddf;
  int i, j, Nlocal;

  w1 = vector(qm1);
  z1 = vector(qm1);
  w2 = vector(qm2);
  z2 = vector(qm2);  
  dn = (double)pn;
  sdn = sqrt(dn);
  b = -qPHI(truncate/2.)/sdn;
  gausslegendre(qm1, -b, b, z1, w1);  
  ddf = (double)(df);  
  b1 = sqrt(qCHI(     truncate/2., df)/ddf);
  b2 = sqrt(qCHI(1. - truncate/2., df)/ddf);
  w2 = vector(qm2);
  z2 = vector(qm2);
  gausslegendre(qm2, b1, b2, z2, w2);    
  result = 0.;
  for (j=0; j<qm2; j++) {
    Nlocal = qm_for_l_and_c(l, z2[j]*c);
    for (i=0; i<qm1; i++)     
      result += w1[i]*sdn*phi( z1[i]*sdn, 0.) * w2[j]*2.*ddf*z2[j]*chi( ddf*z2[j]*z2[j], df) * xe2_iglarl(l, z2[j]*c, hs, z1[i]+mu, Nlocal);  
  }
  Free(w1);
  Free(z1);
  Free(w2);
  Free(z2);  
  return result;
}


/* 2. ARL (varying limits and conditional) */

double xe2_arlm_prerun_MU(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate)
{ double *w, *z, b, result1, result0, dn, sdn, *pair;
  int i, Nlocal, fahne;

  w = vector(qm);
  z = vector(qm);
  pair = vector(2);
  dn = (double)pn;
  sdn = sqrt(dn);
  b = -qPHI(truncate/2.)/sdn;
  gausslegendre(qm, -b, b, z, w);
  Nlocal = qm_for_l_and_c(l, c);  
  result1 = 0.;
  result0 = 0.;
  for (i=0; i<qm; i++) {
    fahne = xe2_arlm_special(l, c, hs, q, z[i]+mu0, z[i]+mu1, mode, Nlocal, nmax, pair);
    if ( fahne!= 0 ) warning("something happened with xe2_arlm_special");
    result1 += w[i] * sdn*phi( z[i]*sdn, 0.) * pair[1];
    result0 += w[i] * sdn*phi( z[i]*sdn, 0.) * pair[0];
  }
  result1 /= result0;
  Free(pair);
  Free(w);
  Free(z);
  return result1;
}


double xe2_arlm_prerun_SIGMA(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate)
{ double *w, *z, b1, b2, result1, result0, ddf, *pair;
  int i, Nlocal, fahne;
 
  w = vector(qm);
  z = vector(qm);
  pair = vector(2);
  ddf = (double)(pn-1);  
  b1 = sqrt(qCHI(     truncate/2., pn-1)/ddf);
  b2 = sqrt(qCHI(1. - truncate/2., pn-1)/ddf); 
  gausslegendre(qm, b1, b2, z, w);  
  result1 = 0.;
  result0 = 0.;
  for (i=0; i<qm; i++) {
    Nlocal = qm_for_l_and_c(l, z[i]*c);
    fahne = xe2_arlm_special(l, z[i]*c, hs, q, mu0, mu1, mode, Nlocal, nmax, pair);
    if ( fahne!= 0 ) warning("something happened with xe2_arlm_special");
    result1 += w[i] * 2.*ddf*z[i]*chi( ddf*z[i]*z[i], pn-1) * pair[1];
    result0 += w[i] * 2.*ddf*z[i]*chi( ddf*z[i]*z[i], pn-1) * pair[0];
  }
  result1 /= result0;
  Free(pair);  
  Free(w);
  Free(z);
  return result1;
}


double xe2_arlm_prerun_BOTH(double l, double c, double hs, int q, double mu0, double mu1, int pn, int df, int mode, int nmax, int qm1, int qm2, double truncate)
{ double *w1, *z1, *w2, *z2, b, b1, b2, result1, result0, dn, sdn, ddf, *pair;
  int i, j, Nlocal, fahne;

  w1 = vector(qm1);
  z1 = vector(qm1);
  w2 = vector(qm2);
  z2 = vector(qm2);  
  pair = vector(2);
  dn = (double)pn;
  sdn = sqrt(dn);
  b = -qPHI(truncate/2.)/sdn;  
  gausslegendre(qm1, -b, b, z1, w1);  
  ddf = (double)(df);  
  b1 = sqrt(qCHI(     truncate/2., df)/ddf);
  b2 = sqrt(qCHI(1. - truncate/2., df)/ddf);
  w2 = vector(qm2);
  z2 = vector(qm2);
  gausslegendre(qm2, b1, b2, z2, w2);    
  result1 = 0.;
  result0 = 0.;
  for (j=0; j<qm2; j++) {
    Nlocal = qm_for_l_and_c(l, z2[j]*c);
    for (i=0; i<qm1; i++) {
      fahne = xe2_arlm_special(l, z2[j]*c, hs, q, z1[i]+mu0, z1[i]+mu1, mode, Nlocal, nmax, pair);
      if ( fahne!= 0 ) warning("something happened with xe2_arlm_special");
      result1 += w1[i]*sdn*phi( z1[i]*sdn, 0.) * w2[j]*2.*ddf*z2[j]*chi( ddf*z2[j]*z2[j], df) * pair[1];
      result0 += w1[i]*sdn*phi( z1[i]*sdn, 0.) * w2[j]*2.*ddf*z2[j]*chi( ddf*z2[j]*z2[j], df) * pair[0];
    }
  }
  result1 /= result0;
  Free(pair);
  Free(w1);
  Free(z1);
  Free(w2);
  Free(z2);  
  return result1;
}



/* some helper functions */


double xe2_sf_deluxe(double l, double c, double hs, double mu, int N, int nmax, double BOUND, double *p0, int *nstop, double *rho)
{ double *Sm, *Pn, *w, *z, mn_minus=1., mn_plus=0., ratio;
  int i, j, n;

 c  *= sqrt( l/(2.-l) );
 hs *= sqrt( l/(2.-l) );
 Sm = matrix(N, N);
 w  = vector(N);
 z  = vector(N);
 Pn = matrix(nmax, N);
 gausslegendre(N, -c, c, z, w);
 
 *nstop = 0;

 for (i=0; i<N; i++)
   for (j=0; j<N; j++)
     Sm[i*N+j] = w[j]/l * phi( (z[j]-(1.-l)*z[i])/l, mu);

 for (n=1; n<=nmax; n++) {
   if ( n==1 )
     for (i=0; i<N; i++)
       Pn[i] = PHI( (c-(1.-l)*z[i])/l, mu) - PHI( (-c-(1.-l)*z[i])/l, mu);
   else
     for (i=0; i<N; i++) {
       Pn[(n-1)*N+i] = 0.;
       for (j=0; j<N; j++)
	 Pn[(n-1)*N+i] += Sm[i*N+j] * Pn[(n-2)*N+j];
     }

   if ( n==1 )
     p0[0] = PHI( (c-(1.-l)*hs)/l, mu) - PHI( (-c-(1.-l)*hs)/l, mu);
   else {
     p0[n-1] = 0.;
     for (j=0; j<N; j++)
       p0[n-1] += w[j]/l * phi( (z[j]-(1.-l)*hs)/l, mu) * Pn[(n-2)*N+j];
   }
   
   mn_minus = 1.; mn_plus = 0.;
   if ( n>1 ) {
     for (i=0;i<N;i++) {
       if (Pn[(n-2)*N+i]==0)
         if (Pn[(n-1)*N+i]==0) ratio = 0.;
         else ratio = 1.;
       else ratio = Pn[(n-1)*N+i]/Pn[(n-2)*N+i];
       if ( ratio<mn_minus ) mn_minus = ratio;
       if ( ratio>mn_plus ) mn_plus = ratio;
     }
     *rho = (mn_minus + mn_plus)/2.;
     if ( fabs(mn_plus - mn_minus) < BOUND ) {       
       *nstop = n;
       n = nmax + 1;
     }
   }
 }
 
 Free(Pn);
 Free(z);
 Free(w);
 Free(Sm);

 return 0;
}


double xe2_sfm_simple(double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax, double *p0)
{ double *Smatrix, *fn, *w, *z, dn, rn, cn, rn0, cn0, delta=0., nn, fSt, aSt;
  int i, j, n;

 fSt = 0.5;
 aSt = ( -2./log10(1.-fSt) - 1.)/19.;

 c  *= sqrt( l/(2.-l) ); 
 hs *= sqrt( l/(2.-l) );
 if ( mode==fir || mode==both ) delta = 2.*hs;

 Smatrix = matrix(N, N);
 w       = vector(N);
 z       = vector(N);
 fn      = matrix(nmax, N);

 gausslegendre(N, -c, c, z, w);

 rn = 1.; cn = 0.; rn0 = 1., cn0 = 0.;

 /* in-control, i. e. n<=q-1 */

 for (n=1; n<=q-1; n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=1,2,...,q-1 */
  switch ( mode ) {
    case vacl: rn = sqrt( 1. - pow(1.-l,2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }
 
  /* determine f_n, n=1,2,...,q-1 */
  if ( n==1 ) {
    for (i=0; i<N; i++)
      if ( mode==stat )
        fn[0*N+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)),mu0);
      else
        fn[0*N+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l,mu0);
  } 
  else {
    for (i=0; i<N; i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0; j<N; j++) {
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l*phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l,mu0);
      }
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0; i<N; i++) p0[n-1] += w[i] * fn[(n-1)*N+i];

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;
 }

 /* out-of-control, i.e. t>=q */

 for (n=q; n<=nmax; n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=q,q+1,... */
  switch ( mode ) {
    case vacl: rn = sqrt( 1. - pow(1.-l,2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }

  /* determine f_n, n=q,q+1,... */
  if ( n==1 ) {
    for (i=0; i<N; i++)
      if ( mode==stat )
        fn[0*N+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)),mu1);
      else
        fn[0*N+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l,mu1);
  }
  else {
    for (i=0; i<N; i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0; j<N; j++) 
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l*phi( (cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l,mu1);
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0; i<N; i++) p0[n-1] += w[i] * fn[(n-1)*N+i];

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;
 }

 Free(Smatrix);
 Free(w);
 Free(z);
 Free(fn);

 return 0;
}


double xe2_sfm_deluxe(double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax, double BOUND, double *p0, int *nstop, double *rho)
{ double *Smatrix, *fn, *w, *z, dn, rn, cn, rn0, cn0, delta=0., nn, fSt, aSt, mn_minus, mn_plus, ratio;
  int i, j, n;

 fSt = 0.5;
 aSt = ( -2./log10(1.-fSt) - 1.)/19.;
 c  *= sqrt( l/(2.-l) ); 
 hs *= sqrt( l/(2.-l) );
 if ( mode==fir || mode==both ) delta = 2.*hs;
 Smatrix = matrix(N, N);
 w       = vector(N);
 z       = vector(N);
 fn      = matrix(nmax, N);
 gausslegendre(N, -c, c, z, w);
 *nstop = 0;

 rn = 1.; cn = 0.; rn0 = 1., cn0 = 0.;

 /* in-control, i. e. n<=q-1 */

 for (n=1; n<=q-1; n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=1,2,...,q-1 */
  switch ( mode ) {
    case vacl: rn = sqrt( 1. - pow(1.-l,2.*nn) );
         break;
    case fir: dn = delta*pow(1.-l,nn);
              rn = 1. - dn/(2.*c);
              cn = dn/2.;
         break;
    case both: dn = delta*pow(1.-l,nn);
               rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c);
               cn = dn/2.;
         break;
    case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
         break;
  }
 
  /* determine f_n, n=1,2,...,q-1 */
  if ( n==1 ) {
    for (i=0; i<N; i++)
      if ( mode==stat )
        fn[0*N+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)), mu0);
      else
        fn[0*N+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l, mu0);
  } 
  else {
    for (i=0; i<N; i++) {
      fn[(n-1)*N+i] = 0.;
      for (j=0; j<N; j++) {
        fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l*phi((cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l, mu0);
      }
    }
  }

  /* determine P(L>n), n=1,2,...,q-1 */
  p0[n-1] = 0.;
  for (i=0; i<N; i++) p0[n-1] += w[i] * fn[(n-1)*N+i];

  /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
  cn0 = cn; rn0 = rn;
 }

 /* out-of-control, i.e. t>=q */

 for (n=q; n<=nmax; n++) {
   nn = (double) n;

   /* determine c_n and r_n, n=q,q+1,... */
   switch ( mode ) {
     case vacl: rn = sqrt( 1. - pow(1.-l,2.*nn) );
          break;
     case fir: dn = delta*pow(1.-l,nn);
               rn = 1. - dn/(2.*c);
               cn = dn/2.;
          break;
     case both: dn = delta*pow(1.-l,nn);
                rn = sqrt( 1. - pow(1.-l,2.*nn) ) - dn/(2.*c);
                cn = dn/2.;
          break;
     case steiner: rn = sqrt(1.-pow(1.-l,2.*nn))*(1.-pow(1.-fSt,1.+aSt*(nn-1.)));
          break;
   }

   /* determine f_n, n=q,q+1,... */
   if ( n==1 ) {
     for (i=0; i<N; i++)
       if ( mode==stat )
         fn[0*N+i] = 1./sqrt(l/(2.-l))*phi( (cn+rn*z[i])/sqrt(l/(2.-l)), mu1);
       else
         fn[0*N+i] = rn/l * phi( (cn+rn*z[i]-(1.-l)*hs)/l, mu1);
   }
   else {
     for (i=0; i<N; i++) {
       fn[(n-1)*N+i] = 0.;
       for (j=0; j<N; j++) 
         fn[(n-1)*N+i] += w[j]*fn[(n-2)*N+j]*rn/l*phi( (cn+rn*z[i]-(1.-l)*(cn0+rn0*z[j]))/l, mu1);
     }
   }

   /* determine P(L>n), n=1,2,...,q-1 */
   p0[n-1] = 0.;
   for (i=0; i<N; i++) p0[n-1] += w[i] * fn[(n-1)*N+i];
 
   /* weights and nodes w.r.t. O_n become w. a. n. w.r.t. O_n-1 */
   cn0 = cn; rn0 = rn;
   
   mn_minus = 1.; mn_plus = 0.;
   if ( n>q ) {
     for (i=0;i<N;i++) {
       if (fn[(n-2)*N+i]==0)
         if (fn[(n-1)*N+i]==0) ratio = 0.; else ratio = 1.;
       else ratio = fn[(n-1)*N+i]/fn[(n-2)*N+i];
       if ( ratio<mn_minus ) mn_minus = ratio;
       if ( ratio>mn_plus ) mn_plus = ratio;
     }
     *rho = (mn_minus + mn_plus)/2.;
     if ( fabs(mn_plus - mn_minus) < BOUND ) {
       *nstop = n;
       n = nmax + 1;
     }
   }
 }

 Free(Smatrix);
 Free(w);
 Free(z);
 Free(fn);

 return 0;
}



/* P(L>n) */


double xe2_sf_prerun_MU_deluxe(double l, double c, double hs, double mu, int pn, int nmax, int qm, double truncate, double BOUND, double *p0)
{ double *ww, *zz, b, dn, sdn, *SF, rho;
  int i, m, n, nstop, Nlocal;

 SF = vector(nmax); 
 ww = vector(qm);
 zz = vector(qm);
 
 dn = (double)pn;
 sdn = sqrt(dn);
 b = -qPHI(truncate/2.)/sdn;  
 gausslegendre(qm, -b, b, zz, ww); 
 for (i=0; i<qm; i++) ww[i] *= sdn*phi( zz[i]*sdn, 0. );
 
 for (n=0; n<nmax; n++) p0[n] = 0.;

 Nlocal = qm_for_l_and_c(l, c);
 for (i=0; i<qm; i++) {     
   m = xe2_sf_deluxe(l, c, hs, zz[i]+mu, Nlocal, nmax, BOUND, SF, &nstop, &rho);
   if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
   if ( nstop > 0 ) {
     for (n=0; n<nstop; n++)  p0[n] += ww[i] * SF[n];
     for (n=nstop; n<nmax; n++)  p0[n] += ww[i] * SF[nstop-1] * pow(rho, n-nstop+1);
   } else {
     for (n=0; n<nmax; n++)  p0[n] += ww[i] * SF[n];
   }
 }
 
 Free(ww);
 Free(zz); 
 Free(SF);
 
 return 0;
}


double xe2_sf_prerun_MU(double l, double c, double hs, double mu, int pn, int nmax, int qm, double truncate, double *p0)
{ double *ww, *zz, b, dn, sdn, *SF;
  int i, m, n, Nlocal;

 SF = vector(nmax); 
 ww = vector(qm);
 zz = vector(qm);
 
 dn = (double)pn;
 sdn = sqrt(dn);
 b = -qPHI(truncate/2.)/sdn;  
 gausslegendre(qm, -b, b, zz, ww); 
 for (i=0; i<qm; i++) ww[i] *= sdn*phi( zz[i]*sdn, 0. );

 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 Nlocal = qm_for_l_and_c(l, c);
 for (i=0; i<qm; i++) {     
   m = xe2_sf(l, c, hs, zz[i]+mu, Nlocal, nmax, SF);
   if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sf");
   for (n=0; n<nmax; n++)  p0[n] += ww[i] * SF[n];
 }
 
 Free(ww);
 Free(zz); 
 Free(SF);
 
 return 0;
}


double xe2_sfm_prerun_MU_deluxe(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate, double BOUND, double *p0)
{ double *ww, *zz, b, dn, sdn, *SF, rho;
  int i, m, n, nstop, Nlocal;

 SF = vector(nmax); 
 ww = vector(qm);
 zz = vector(qm);
 
 dn = (double)pn;
 sdn = sqrt(dn);
 b = -qPHI(truncate/2.)/sdn;  
 gausslegendre(qm, -b, b, zz, ww); 
 for (i=0; i<qm; i++) ww[i] *= sdn*phi( zz[i]*sdn, 0. );
 
 for (n=0; n<nmax; n++) p0[n] = 0.;

 Nlocal = qm_for_l_and_c(l, c);
 for (i=0; i<qm; i++) {     
   m = xe2_sfm_deluxe(l, c, hs, q, zz[i]+mu0, zz[i]+mu1, mode, Nlocal, nmax, BOUND, SF, &nstop, &rho);
   if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
   if ( nstop > 0 ) {
     for (n=0; n<nstop; n++)  p0[n] += ww[i] * SF[n];
     for (n=nstop; n<nmax; n++)  p0[n] += ww[i] * SF[nstop-1] * pow(rho, n-nstop+1);
   } else {
     for (n=0; n<nmax; n++)  p0[n] += ww[i] * SF[n];
   }
 }
 
 if ( q > 1 ) for (n=q-1; n<nmax; n++) p0[n] /= p0[q-2];
 
 Free(ww);
 Free(zz); 
 Free(SF);
 
 return 0;
}


double xe2_sfm_prerun_MU(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate, double *p0)
{ double *ww, *zz, b, dn, sdn, *SF;
  int i, m, n, Nlocal;

 SF = vector(nmax); 
 ww = vector(qm);
 zz = vector(qm);
 
 dn = (double)pn;
 sdn = sqrt(dn);
 b = -qPHI(truncate/2.)/sdn;  
 gausslegendre(qm, -b, b, zz, ww); 
 for (i=0; i<qm; i++) ww[i] *= sdn*phi( zz[i]*sdn, 0. );

 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 Nlocal = qm_for_l_and_c(l, c);
 for (i=0; i<qm; i++) {     
   m = xe2_sfm_simple(l, c, hs, q, zz[i]+mu0, zz[i]+mu1, mode, Nlocal, nmax, SF);
   if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sfm");
   for (n=0; n<nmax; n++)  p0[n] += ww[i] * SF[n];
 }
 
 if ( q > 1 ) for (n=q-1; n<nmax; n++) p0[n] /= p0[q-2];
 
 Free(ww);
 Free(zz); 
 Free(SF);
 
 return 0;
}


double xe2_sf_prerun_SIGMA_deluxe(double l, double c, double hs, double mu, int pn, int nmax, int qm, double truncate, double BOUND, double *p0)
{ double *ww, *zz, b1, b2, ddf, *SF, rho;
  int i, m, n, nstop, Nlocal;

 SF = vector(nmax);
 ww = vector(qm);
 zz = vector(qm);
 
 ddf = (double)(pn-1);  
 b1 = sqrt(qCHI(     truncate/2., pn-1)/ddf);
 b2 = sqrt(qCHI(1. - truncate/2., pn-1)/ddf); 
 gausslegendre(qm, b1, b2, zz, ww); 
 for (i=0; i<qm; i++) ww[i] *= 2.*ddf*zz[i] * chi( ddf*zz[i]*zz[i], pn-1);
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm; i++) {
    Nlocal = qm_for_l_and_c(l, zz[i]*c);
    m = xe2_sf_deluxe(l, zz[i]*c, hs, mu, Nlocal, nmax, BOUND, SF, &nstop, &rho);
    if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
    if ( nstop > 0 ) {
      for (n=0; n<nstop; n++)  p0[n] += ww[i] * SF[n];
      for (n=nstop; n<nmax; n++)  p0[n] += ww[i] * SF[nstop-1] * pow(rho, n-nstop+1);
    } else {
      for (n=0; n<nmax; n++)  p0[n] += ww[i] * SF[n];
    }
 }  
 
 Free(ww);
 Free(zz); 
 Free(SF);
 
 return 0;
}


double xe2_sf_prerun_SIGMA(double l, double c, double hs, double mu, int pn, int nmax, int qm, double truncate, double *p0)
{ double *ww, *zz, b1, b2, ddf, *SF;
  int i, m, n, Nlocal;

 SF = vector(nmax);
 ww = vector(qm);
 zz = vector(qm);
 
 ddf = (double)(pn-1);  
 b1 = sqrt(qCHI(     truncate/2., pn-1)/ddf);
 b2 = sqrt(qCHI(1. - truncate/2., pn-1)/ddf); 
 gausslegendre(qm, b1, b2, zz, ww); 
 for (i=0; i<qm; i++) ww[i] *= 2.*ddf*zz[i] * chi( ddf*zz[i]*zz[i], pn-1);
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm; i++) {
    Nlocal = qm_for_l_and_c(l, zz[i]*c);
    m = xe2_sf(l, zz[i]*c, hs, mu, Nlocal, nmax, SF);
    if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sf");
    for (n=0; n<nmax; n++)  p0[n] += ww[i] * SF[n];
 }  
 
 Free(ww);
 Free(zz); 
 Free(SF);
 
 return 0;
}


double xe2_sfm_prerun_SIGMA_deluxe(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate, double BOUND, double *p0)
{ double *ww, *zz, b1, b2, ddf, *SF, rho;
  int i, m, n, nstop, Nlocal;

 SF = vector(nmax);
 ww = vector(qm);
 zz = vector(qm);
 
 ddf = (double)(pn-1);  
 b1 = sqrt(qCHI(     truncate/2., pn-1)/ddf);
 b2 = sqrt(qCHI(1. - truncate/2., pn-1)/ddf); 
 gausslegendre(qm, b1, b2, zz, ww); 
 for (i=0; i<qm; i++) ww[i] *= 2.*ddf*zz[i] * chi( ddf*zz[i]*zz[i], pn-1);
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm; i++) {
    Nlocal = qm_for_l_and_c(l, zz[i]*c);
    m = xe2_sfm_deluxe(l, zz[i]*c, hs, q, mu0, mu1, mode, Nlocal, nmax, BOUND, SF, &nstop, &rho);
    if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
    if ( nstop > 0 ) {
      for (n=0; n<nstop; n++)  p0[n] += ww[i] * SF[n];
      for (n=nstop; n<nmax; n++)  p0[n] += ww[i] * SF[nstop-1] * pow(rho, n-nstop+1);
    } else {
      for (n=0; n<nmax; n++)  p0[n] += ww[i] * SF[n];
    }
 }  
 
 if ( q > 1 ) for (n=q-1; n<nmax; n++) p0[n] /= p0[q-2];
 
 Free(ww);
 Free(zz); 
 Free(SF);
 
 return 0;
}


double xe2_sfm_prerun_SIGMA(double l, double c, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate, double *p0)
{ double *ww, *zz, b1, b2, ddf, *SF;
  int i, m, n, Nlocal;

 SF = vector(nmax);
 ww = vector(qm);
 zz = vector(qm);
 
 ddf = (double)(pn-1);  
 b1 = sqrt(qCHI(     truncate/2., pn-1)/ddf);
 b2 = sqrt(qCHI(1. - truncate/2., pn-1)/ddf); 
 gausslegendre(qm, b1, b2, zz, ww); 
 for (i=0; i<qm; i++) ww[i] *= 2.*ddf*zz[i] * chi( ddf*zz[i]*zz[i], pn-1);
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm; i++) {
    Nlocal = qm_for_l_and_c(l, zz[i]*c);
    m = xe2_sfm_simple(l, zz[i]*c, hs, q, mu0, mu1, mode, Nlocal, nmax, SF);
    if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sfm");
    for (n=0; n<nmax; n++)  p0[n] += ww[i] * SF[n];
 }  
 
 if ( q > 1 ) for (n=q-1; n<nmax; n++) p0[n] /= p0[q-2];
 
 Free(ww);
 Free(zz); 
 Free(SF);
 
 return 0;
}


double xe2_sf_prerun_BOTH_deluxe(double l, double c, double hs, double mu, int pn, int df, int nmax, int qm1, int qm2, double truncate, double BOUND, double *p0)
{ double *ww1, *zz1, *ww2, *zz2, b, b1, b2, dn, sdn, ddf, *SF, rho;
  int i, j, m, n, nstop, Nlocal;

 SF = vector(nmax); 
 ww1 = vector(qm1);
 zz1 = vector(qm1); 
 ww2 = vector(qm2);
 zz2 = vector(qm2);

 dn = (double)pn;
 sdn = sqrt(dn);
 b = -qPHI(truncate/2.)/sdn;  
 gausslegendre(qm1, -b, b, zz1, ww1);
 for (i=0; i<qm1; i++) ww1[i] *= sdn * phi( zz1[i]*sdn, 0.);
 
 ddf = (double)(df);  
 b1 = sqrt(qCHI(     truncate/2., df)/ddf);
 b2 = sqrt(qCHI(1. - truncate/2., df)/ddf);
 gausslegendre(qm2, b1, b2, zz2, ww2); 
 for (j=0; j<qm2; j++) ww2[j] *= 2.*ddf*zz2[j] * chi( ddf*zz2[j]*zz2[j], df); 
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm1; i++) {
   for (j=0; j<qm2; j++) {
     Nlocal = qm_for_l_and_c(l, zz2[j]*c);
     m = xe2_sf_deluxe(l, zz2[j]*c, hs, zz1[i]+mu, Nlocal, nmax, BOUND, SF, &nstop, &rho);
     if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
     if ( nstop > 0 ) {
       for (n=0; n<nstop; n++)  p0[n] += ww1[i] * ww2[j] * SF[n];
       for (n=nstop; n<nmax; n++)  p0[n] += ww1[i] * ww2[j] * SF[nstop-1] * pow(rho, n-nstop+1);
     } else {
       for (n=0; n<nmax; n++)  p0[n] += ww1[i] * ww2[j] * SF[n];
     }
   }
 }
 
 Free(ww1);
 Free(zz1); 
 Free(ww2);
 Free(zz2);
 Free(SF);
 
 return 0;
}


double xe2_sf_prerun_BOTH(double l, double c, double hs, double mu, int pn, int df, int nmax, int qm1, int qm2, double truncate, double *p0)
{ double *ww1, *zz1, *ww2, *zz2, b, b1, b2, dn, sdn, ddf, *SF;
  int i, j, m, n, Nlocal;

 SF = vector(nmax); 
 ww1 = vector(qm1);
 zz1 = vector(qm1); 
 ww2 = vector(qm2);
 zz2 = vector(qm2);

 dn = (double)pn;
 sdn = sqrt(dn);
 b = -qPHI(truncate/2.)/sdn;  
 gausslegendre(qm1, -b, b, zz1, ww1);
 for (i=0; i<qm1; i++) ww1[i] *= sdn * phi( zz1[i]*sdn, 0.);
 
 ddf = (double)(df);  
 b1 = sqrt(qCHI(     truncate/2., df)/ddf);
 b2 = sqrt(qCHI(1. - truncate/2., df)/ddf);
 gausslegendre(qm2, b1, b2, zz2, ww2); 
 for (j=0; j<qm2; j++) ww2[j] *= 2.*ddf*zz2[j] * chi( ddf*zz2[j]*zz2[j], df); 
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm1; i++) {
   for (j=0; j<qm2; j++) {
     Nlocal = qm_for_l_and_c(l, zz2[j]*c);
     m = xe2_sf(l, zz2[j]*c, hs, zz1[i]+mu, Nlocal, nmax, SF);
     if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sf");
     for (n=0; n<nmax; n++)  p0[n] += ww1[i] * ww2[j] * SF[n];
   }
 }
 
 Free(ww1);
 Free(zz1); 
 Free(ww2);
 Free(zz2);
 Free(SF);
 
 return 0;
}



double xe2_sfm_prerun_BOTH_deluxe(double l, double c, double hs, int q, double mu0, double mu1, int pn, int df, int mode, int nmax, int qm1, int qm2, double truncate, double BOUND, double *p0)
{ double *ww1, *zz1, *ww2, *zz2, b, b1, b2, dn, sdn, ddf, *SF, rho;
  int i, j, m, n, nstop, Nlocal;

 SF = vector(nmax); 
 ww1 = vector(qm1);
 zz1 = vector(qm1); 
 ww2 = vector(qm2);
 zz2 = vector(qm2);

 dn = (double)pn;
 sdn = sqrt(dn);
 b = -qPHI(truncate/2.)/sdn;  
 gausslegendre(qm1, -b, b, zz1, ww1);
 for (i=0; i<qm1; i++) ww1[i] *= sdn * phi( zz1[i]*sdn, 0.);
 
 ddf = (double)(df);  
 b1 = sqrt(qCHI(     truncate/2., df)/ddf);
 b2 = sqrt(qCHI(1. - truncate/2., df)/ddf);
 gausslegendre(qm2, b1, b2, zz2, ww2); 
 for (j=0; j<qm2; j++) ww2[j] *= 2.*ddf*zz2[j] * chi( ddf*zz2[j]*zz2[j], df); 
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm1; i++) {
   for (j=0; j<qm2; j++) {
     Nlocal = qm_for_l_and_c(l, zz2[j]*c);
     m = xe2_sfm_deluxe(l, zz2[j]*c, hs, q, zz1[i]+mu0, zz1[i]+mu1, mode, Nlocal, nmax, BOUND, SF, &nstop, &rho);
     if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
     if ( nstop > 0 ) {
       for (n=0; n<nstop; n++)  p0[n] += ww1[i] * ww2[j] * SF[n];
       for (n=nstop; n<nmax; n++)  p0[n] += ww1[i] * ww2[j] * SF[nstop-1] * pow(rho, n-nstop+1);
     } else {
       for (n=0; n<nmax; n++)  p0[n] += ww1[i] * ww2[j] * SF[n];
     }
   }
 }
 
 if ( q > 1 ) for (n=q-1; n<nmax; n++) p0[n] /= p0[q-2];
 
 Free(ww1);
 Free(zz1); 
 Free(ww2);
 Free(zz2);
 Free(SF);
 
 return 0;
}


double xe2_sfm_prerun_BOTH(double l, double c, double hs, int q, double mu0, double mu1, int pn, int df, int mode, int nmax, int qm1, int qm2, double truncate, double *p0)
{ double *ww1, *zz1, *ww2, *zz2, b, b1, b2, dn, sdn, ddf, *SF;
  int i, j, m, n, Nlocal;

 SF = vector(nmax); 
 ww1 = vector(qm1);
 zz1 = vector(qm1); 
 ww2 = vector(qm2);
 zz2 = vector(qm2);

 dn = (double)pn;
 sdn = sqrt(dn);
 b = -qPHI(truncate/2.)/sdn;  
 gausslegendre(qm1, -b, b, zz1, ww1);
 for (i=0; i<qm1; i++) ww1[i] *= sdn * phi( zz1[i]*sdn, 0.);
 
 ddf = (double)(df);  
 b1 = sqrt(qCHI(     truncate/2., df)/ddf);
 b2 = sqrt(qCHI(1. - truncate/2., df)/ddf);
 gausslegendre(qm2, b1, b2, zz2, ww2); 
 for (j=0; j<qm2; j++) ww2[j] *= 2.*ddf*zz2[j] * chi( ddf*zz2[j]*zz2[j], df); 
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm1; i++) {
   for (j=0; j<qm2; j++) {
     Nlocal = qm_for_l_and_c(l, zz2[j]*c);
     m = xe2_sfm_simple(l, zz2[j]*c, hs, q, zz1[i]+mu0, zz1[i]+mu1, mode, Nlocal, nmax, SF);
     if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sfm");
     for (n=0; n<nmax; n++)  p0[n] += ww1[i] * ww2[j] * SF[n];
   }
 }
 
 if ( q > 1 ) for (n=q-1; n<nmax; n++) p0[n] /= p0[q-2];
 
 Free(ww1);
 Free(zz1); 
 Free(ww2);
 Free(zz2);
 Free(SF);
 
 return 0;
}


/* quantile function */


double xe2_Wq_prerun_MU_deluxe(double l, double c, double p, double hs, double mu, int pn, int nmax, int qm, double truncate, double BOUND)
{ double *ww, *zz, b, dn, sdn, *SF, *p0, rho, *rhomany, *SFlast, Lp=-1., sf_level_adj;
  int i, j, n, nstop, nstop_, nn, nsm, qnspecial=0, Nlocal;

 p0 = vector(nmax);
 SF = vector(nmax);
 
 rhomany = vector(qm);
 SFlast = vector(qm);
 ww = vector(qm);
 zz = vector(qm); 
 
 dn = (double)pn;
 sdn = sqrt(dn);
 b = -qPHI(truncate/2.)/sdn;  
 gausslegendre(qm, -b, b, zz, ww); 
 for (i=0; i<qm; i++) ww[i] *= sdn*phi( zz[i]*sdn, 0. );
 
 Nlocal = qm_for_l_and_c(l, c);
 
 qnspecial = (qm+1) / 2;
  
 j = xe2_sf_deluxe(l, c, hs, zz[qnspecial]+mu, Nlocal, nmax, BOUND, SF, &nsm, &rho);
 if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
 n = nsm;
 
 if ( nsm < 1 ) { /* did not converge yet -- should be the rare case */
   nn = nmax;
   warning("The geometric tail approximation might not work.");
 } else {
   nstop = nsm;
   
   j = xe2_sf_deluxe(l, c, hs, zz[qnspecial+1]+mu, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax ) {
       nstop = nstop_;
       i++;
       j = xe2_sf_deluxe(l, c, hs, zz[qnspecial+i]+mu, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
       if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }
   
   nstop = n;
   j = xe2_sf_deluxe(l, c, hs, zz[qnspecial-1]+mu, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax ) {
       nstop = nstop_;
       i++;
       j = xe2_sf_deluxe(l, c, hs, zz[qnspecial-i]+mu, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
       if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }  
   nn = nsm;
 }
  
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm; i++) {
   j = xe2_sf_deluxe(l, c, hs, zz[i]+mu, Nlocal, nn, BOUND, SF, &nstop, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
   if ( nstop < 1 ) {
     nstop = nn;
     warning("The geometric tail approximation might not work.");     
   }
   rhomany[i] = rho;        
   for (n=0; n<nstop; n++)  p0[n] += ww[i] * SF[n];
   if ( nstop < nn) {
     for (n=nstop; n<nn; n++) p0[n] += ww[i] * SF[nstop-1] * pow(rho, n-nstop+1);          
   }
   SFlast[i] = SF[nstop-1] * pow(rho, nn-nstop);
 } 
 
 sf_level_adj = 1. - p;
 if ( p0[nn-1] <= sf_level_adj ) { 
   n = nn-1;
   while (  p0[n] <= sf_level_adj &&  n > 0  ) n--;
   if ( p0[n] > 1.-p ) Lp = (double)( n + 2 ); else Lp = 1.;
 } else {
   for (n=nn; n<nmax; n++) {
      p0[n] = 0.;
      for (i=0; i<qm; i++) p0[n] += ww[i] * SFlast[i] * pow(rhomany[i], n-nn+1);
      if (  p0[n] <= sf_level_adj  ) {
        Lp = (double)( n + 1 );
        n = nmax+1;
      }
   }
 }
 
 Free(p0);
 Free(ww);
 Free(zz); 
 Free(SF);
 Free(SFlast);
 Free(rhomany);

 return Lp;
}


double xe2_Wqm_prerun_MU_deluxe(double l, double c, double p, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate, double BOUND)
{ double *ww, *zz, b, dn, sdn, *SF, *p0, rho, *rhomany, *SFlast, Lp=-1., sf_level_adj;
  int i, j, n, nstop, nstop_, nn, nsm, qnspecial=0, Nlocal;

 p0 = vector(nmax);
 SF = vector(nmax);
 
 rhomany = vector(qm);
 SFlast = vector(qm);
 ww = vector(qm);
 zz = vector(qm); 
 
 dn = (double)pn;
 sdn = sqrt(dn);
 b = -qPHI(truncate/2.)/sdn;  
 gausslegendre(qm, -b, b, zz, ww); 
 for (i=0; i<qm; i++) ww[i] *= sdn*phi( zz[i]*sdn, 0. );
 
 Nlocal = qm_for_l_and_c(l, c);
 
 qnspecial = (qm+1) / 2;
 
 j = xe2_sfm_deluxe(l, c, hs, q, zz[qnspecial]+mu0, zz[qnspecial]+mu1, mode, Nlocal, nmax, BOUND, SF, &nsm, &rho);
 if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
 n = nsm;
  
 if ( nsm < 1 ) { /* did not converge yet -- should be the rare case */
   nn = nmax;
   warning("The geometric tail approximation might not work.");
 } else {
   nstop = nsm;
   
   j = xe2_sfm_deluxe(l, c, hs, q, zz[qnspecial+1]+mu0, zz[qnspecial+1]+mu1, mode, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax) {
       nstop = nstop_;
       i++;
       j = xe2_sfm_deluxe(l, c, hs, q, zz[qnspecial+i]+mu0, zz[qnspecial+i]+mu1, mode, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
       if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }
   
   nstop = n;
   j = xe2_sfm_deluxe(l, c, hs, q, zz[qnspecial-1]+mu0, zz[qnspecial-1]+mu1, mode, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax ) {
       nstop = nstop_;
       i++;
       j = xe2_sfm_deluxe(l, c, hs, q, zz[qnspecial-i]+mu0, zz[qnspecial-i]+mu1, mode, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
       if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }  
   nn = nsm;
 } 
  
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm; i++) {
   j = xe2_sfm_deluxe(l, c, hs, q, zz[i]+mu0, zz[i]+mu1, mode, Nlocal, nn, BOUND, SF, &nstop, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
   if ( nstop < 1 ) {
     nstop = nn;
     warning("The geometric tail approximation might not work.");     
   }
   rhomany[i] = rho;        
   for (n=0; n<nstop; n++)  p0[n] += ww[i] * SF[n];
   if ( nstop < nn) {
     for (n=nstop; n<nn; n++) p0[n] += ww[i] * SF[nstop-1] * pow(rho, n-nstop+1);          
   }
   SFlast[i] = SF[nstop-1] * pow(rho, nn-nstop);
 } 

 sf_level_adj = 1. - p;
 if ( q > 1 ) sf_level_adj *= p0[q-2];
 if ( p0[nn-1] <= sf_level_adj ) { 
   n = nn-1;
   while (  p0[n] <= sf_level_adj &&  n > 0  ) n--;
   if ( p0[n] > sf_level_adj ) Lp = (double)( n + 2 - q + 1 ); else Lp = 1.;
 } else {
   for (n=nn; n<nmax; n++) {
      p0[n] = 0.;
      for (i=0; i<qm; i++) p0[n] += ww[i] * SFlast[i] * pow(rhomany[i], n-nn+1);
      if (  p0[n] <= sf_level_adj  ) {
        Lp = (double)( n + 1  -  q + 1 );
        n = nmax+1;
      }
   }
 }
 
 Free(p0);
 Free(ww);
 Free(zz); 
 Free(SF);
 Free(SFlast);
 Free(rhomany);

 return Lp;
}


double xe2_Wq_prerun_SIGMA_deluxe(double l, double c, double p, double hs, double mu, int pn, int nmax, int qm, double truncate, double BOUND)
{ double *ww, *zz, b1, b2, ddf, *SF, *p0, rho, *rhomany, *SFlast, Lp=-1., sf_level_adj;
  int i, j, n, nstop, nstop_, nsm, nn, qnspecial=0, Nlocal;

 p0 = vector(nmax);
 SF = vector(nmax);
 
 rhomany = vector(qm);
 SFlast = vector(qm);
 ww = vector(qm);
 zz = vector(qm); 
 
 ddf = (double)(pn-1);  
 b1 = sqrt(qCHI(     truncate/2., pn-1)/ddf);
 b2 = sqrt(qCHI(1. - truncate/2., pn-1)/ddf); 
 gausslegendre(qm, b1, b2, zz, ww); 
 for (i=0; i<qm; i++) ww[i] *= 2.*ddf*zz[i] * chi( ddf*zz[i]*zz[i], pn-1);
 
 /*qnspecial = qm-1;*/
 qnspecial = (qm+1) / 2;
 
 Nlocal = qm_for_l_and_c(l, zz[qnspecial]*c);
 j = xe2_sf_deluxe(l, zz[qnspecial]*c, hs, mu, Nlocal, nmax, BOUND, SF, &nsm, &rho);
 if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
 n = nsm;
 
 if ( nsm < 1 ) { /* did not converge yet -- should be the rare case */
   nn = nmax;
   warning("The geometric tail approximation might not work.");
 } else {
   nstop = nsm;
   
   Nlocal = qm_for_l_and_c(l, zz[qnspecial+1]*c);
   j = xe2_sf_deluxe(l, zz[qnspecial+1]*c, hs, mu, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax ) {
       nstop = nstop_;
       i++;
       Nlocal = qm_for_l_and_c(l, zz[qnspecial+i]*c);
       j = xe2_sf_deluxe(l, zz[qnspecial+i]*c, hs, mu, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
       if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }
 
   nstop = n;
   Nlocal = qm_for_l_and_c(l, zz[qnspecial-1]*c);
   j = xe2_sf_deluxe(l, zz[qnspecial-1]*c, hs, mu, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax ) {
       nstop = nstop_;
       i++;
       Nlocal = qm_for_l_and_c(l, zz[qnspecial-i]*c);
       j = xe2_sf_deluxe(l, zz[qnspecial-i]*c, hs, mu, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
       if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }  
   nn = nsm;
 }
  
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm; i++) {
   Nlocal = qm_for_l_and_c(l, zz[i]*c);
   j = xe2_sf_deluxe(l, zz[i]*c, hs, mu, Nlocal, nn, BOUND, SF, &nstop, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
   if ( nstop < 1 ) {
     nstop = nn;
     warning("The geometric tail approximation might not work.");     
   }
   rhomany[i] = rho;        
   for (n=0; n<nstop; n++)  p0[n] += ww[i] * SF[n];
   if ( nstop < nn) {
     for (n=nstop; n<nn; n++) p0[n] += ww[i] * SF[nstop-1] * pow(rho, n-nstop+1);          
   }
   SFlast[i] = SF[nstop-1] * pow(rho, nn-nstop);
 } 
 
 sf_level_adj = 1. - p;
 if ( p0[nn-1] <= sf_level_adj ) { 
   n = nn-1;
   while (  p0[n] <= sf_level_adj &&  n > 0  ) n--;
   if ( p0[n] > 1.-p ) Lp = (double)( n + 2 ); else Lp = 1.;
 } else {
   for (n=nn; n<nmax; n++) {
      p0[n] = 0.;
      for (i=0; i<qm; i++) p0[n] += ww[i] * SFlast[i] * pow(rhomany[i], n-nn+1);
      if (  p0[n] <= sf_level_adj  ) {
        Lp = (double)( n + 1 );
        n = nmax+1;
      }
   }
 }
 
 Free(p0);
 Free(ww);
 Free(zz); 
 Free(SF);
 Free(SFlast);
 Free(rhomany);

 return Lp;
}


double xe2_Wqm_prerun_SIGMA_deluxe(double l, double c, double p, double hs, int q, double mu0, double mu1, int pn, int mode, int nmax, int qm, double truncate, double BOUND)
{ double *ww, *zz, b1, b2, ddf, *SF, *p0, rho, *rhomany, *SFlast, Lp=-1., sf_level_adj;
  int i, j, n, nstop, nstop_, nsm, nn, qnspecial=0, Nlocal;

 p0 = vector(nmax);
 SF = vector(nmax);
 
 rhomany = vector(qm);
 SFlast = vector(qm);
 ww = vector(qm);
 zz = vector(qm); 
 
 ddf = (double)(pn-1);  
 b1 = sqrt(qCHI(     truncate/2., pn-1)/ddf);
 b2 = sqrt(qCHI(1. - truncate/2., pn-1)/ddf); 
 gausslegendre(qm, b1, b2, zz, ww); 
 for (i=0; i<qm; i++) ww[i] *= 2.*ddf*zz[i] * chi( ddf*zz[i]*zz[i], pn-1);
 
 /*qnspecial = qm-1;*/
 qnspecial = (qm+1) / 2;

 Nlocal = qm_for_l_and_c(l, zz[qnspecial]*c);
 j = xe2_sfm_deluxe(l, zz[qnspecial]*c, hs, q, mu0, mu1, mode, Nlocal, nmax, BOUND, SF, &nsm, &rho);
 if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
 n = nsm;
 
 if ( nsm < 1 ) { /* did not converge yet -- should be the rare case */
   nn = nmax;
   warning("The geometric tail approximation might not work.");
 } else {
   nstop = nsm;
   
   Nlocal = qm_for_l_and_c(l, zz[qnspecial+1]*c);
   j = xe2_sfm_deluxe(l, zz[qnspecial+1]*c, hs, q, mu0, mu1, mode, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax ) {
       nstop = nstop_;
       i++;
       Nlocal = qm_for_l_and_c(l, zz[qnspecial+i]*c);
       j = xe2_sfm_deluxe(l, zz[qnspecial+i]*c, hs, q, mu0, mu1, mode, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
       if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }
 
   nstop = n;
   Nlocal = qm_for_l_and_c(l, zz[qnspecial-1]*c);
   j = xe2_sfm_deluxe(l, zz[qnspecial-1]*c, hs, q, mu0, mu1, mode, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax ) {
       nstop = nstop_;
       i++;
       Nlocal = qm_for_l_and_c(l, zz[qnspecial-i]*c);
       j = xe2_sfm_deluxe(l, zz[qnspecial-i]*c, hs, q, mu0, mu1, mode, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
       if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }  
   nn = nsm;
 }
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm; i++) {
   Nlocal = qm_for_l_and_c(l, zz[i]*c);
   j = xe2_sfm_deluxe(l, zz[i]*c, hs, q, mu0, mu1, mode, Nlocal, nn, BOUND, SF, &nstop, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
   if ( nstop < 1 ) {
     nstop = nn;
     warning("The geometric tail approximation might not work.");     
   }
   rhomany[i] = rho;        
   for (n=0; n<nstop; n++)  p0[n] += ww[i] * SF[n];
   if ( nstop < nn) {
     for (n=nstop; n<nn; n++) p0[n] += ww[i] * SF[nstop-1] * pow(rho, n-nstop+1);          
   }
   SFlast[i] = SF[nstop-1] * pow(rho, nn-nstop);
 } 
 
 sf_level_adj = 1. - p;
 if ( q > 1 ) sf_level_adj *= p0[q-2];
 if ( p0[nn-1] <= sf_level_adj ) { 
   n = nn-1;
   while (  p0[n] <= sf_level_adj &&  n > 0  ) n--;
   if ( p0[n] > sf_level_adj ) Lp = (double)( n + 2 - q + 1 ); else Lp = 1.;
 } else {
   for (n=nn; n<nmax; n++) {
      p0[n] = 0.;
      for (i=0; i<qm; i++) p0[n] += ww[i] * SFlast[i] * pow(rhomany[i], n-nn+1);
      if (  p0[n] <= sf_level_adj  ) {
        Lp = (double)( n + 1  -  q + 1 );
        n = nmax+1;
      }
   }
 }
 
 Free(p0);
 Free(ww);
 Free(zz); 
 Free(SF);
 Free(SFlast);
 Free(rhomany);

 return Lp;
}


double xe2_Wq_prerun_BOTH_deluxe(double l, double c, double p, double hs, double mu, int pn, int df, int nmax, int qm1, int qm2, double truncate, double BOUND)  
{ double *ww1, *zz1, *ww2, *zz2, b, b1, b2, dn, sdn, ddf, *p0, *SF, rho, *rhomany, *SFlast, Lp=-1., sf_level_adj;
  int i, j, m, n, nstop, nstop_, nn, nsm, qnspecial1=0, qnspecial2=0, Nlocal;
 
 p0 = vector(nmax);  
 SF = vector(nmax);
 
 rhomany = vector(qm1*qm2);
 SFlast = vector(qm1*qm2);
 
 ww1 = vector(qm1);
 zz1 = vector(qm1); 
 ww2 = vector(qm2);
 zz2 = vector(qm2);
 
 dn = (double)pn;
 sdn = sqrt(dn);
 b = -qPHI(truncate/2.)/sdn;  
 gausslegendre(qm1, -b, b, zz1, ww1);
 for (i=0; i<qm1; i++) ww1[i] *= sdn * phi( zz1[i]*sdn, 0.);
 qnspecial1 = qm1 - 1;
 
 ddf = (double)(df);  
 b1 = sqrt(qCHI(     truncate/2., df)/ddf);
 b2 = sqrt(qCHI(1. - truncate/2., df)/ddf);
 gausslegendre(qm2, b1, b2, zz2, ww2); 
 for (j=0; j<qm2; j++) ww2[j] *= 2.*ddf*zz2[j] * chi( ddf*zz2[j]*zz2[j], df); 
 qnspecial2 = qm2 - 1;
 
 Nlocal = qm_for_l_and_c(l, zz2[qnspecial2]*c);
 m = xe2_sf_deluxe(l, zz2[qnspecial2]*c, hs, zz1[qnspecial1]+mu, Nlocal, nmax, BOUND, SF, &nsm, &rho);
 if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
 
 if ( nsm < 1 ) { /* did not converge yet -- should be the rare case */
   nn = nmax;
   warning("The geometric tail approximation might not work.");
 } else {
   nstop = nsm;
   
   m = xe2_sf_deluxe(l, zz2[qnspecial2]*c, hs, zz1[qnspecial1-1]+mu, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
   if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax ) {
       nstop = nstop_;
       i++;
       m = xe2_sf_deluxe(l, zz2[qnspecial2]*c, hs, zz1[qnspecial1-i]+mu, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
       if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }
   nn = nsm;
 }
   
 for (n=0; n<nmax; n++) p0[n] = 0.; 

 for (i=0; i<qm1; i++) {
   for (j=0; j<qm2; j++) {
     Nlocal = qm_for_l_and_c(l, zz2[j]*c);
     m = xe2_sf_deluxe(l, zz2[j]*c, hs, zz1[i]+mu, Nlocal, nn, BOUND, SF, &nstop, &rho);
     if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sf_deluxe");
     if ( nstop < 1 ) {
       nstop = nn;
       warning("The geometric tail approximation might not work.");     
     }
     rhomany[i + j*qm1] = rho;
     for (n=0; n<nstop; n++)  p0[n] += ww1[i] * ww2[j] * SF[n];
     if ( nn > nstop ) {
       for (n=nstop; n<nn; n++) p0[n] += ww1[i] * ww2[j] * SF[nstop-1] * pow(rho, n-nstop+1);
     }
     SFlast[i + j*qm1] = SF[nstop-1] * pow(rho, nn-nstop);
   }
 }
 
 sf_level_adj = 1. - p;
 if ( p0[nn-1] <= sf_level_adj ) { 
   n = nn-1;
   while (  p0[n] <= sf_level_adj &&  n > 0  ) n--;
   if ( p0[n] > 1.-p ) Lp = (double)( n + 2 ); else Lp = 1.;
 } else {
   for (n=nn; n<nmax; n++) {
      p0[n] = 0.;
      for (i=0; i<qm1; i++)
	for (j=0; j<qm2; j++)
          p0[n] += ww1[i] * ww2[j] * SFlast[i + j*qm1] * pow(rhomany[i + j*qm1], n-nn+1);
      if (  p0[n] <= sf_level_adj  ) {
        Lp = (double)( n + 1 );
        n = nmax+1;
      }
   }
 }
 
 Free(p0);
 Free(ww1);
 Free(zz1); 
 Free(ww2);
 Free(zz2);
 Free(SF);
 Free(SFlast);
 Free(rhomany);
 
 return Lp;
}


double xe2_Wqm_prerun_BOTH_deluxe(double l, double c, double p, double hs, int q, double mu0, double mu1, int pn, int df, int mode, int nmax, int qm1, int qm2, double truncate, double BOUND)  
{ double *ww1, *zz1, *ww2, *zz2, b, b1, b2, dn, sdn, ddf, *p0, *SF, rho, *rhomany, *SFlast, Lp=-1., sf_level_adj;
  int i, j, m, n, nstop, nstop_, nn, nsm, qnspecial1=0, qnspecial2=0, Nlocal;
 
 p0 = vector(nmax);  
 SF = vector(nmax);
 
 rhomany = vector(qm1*qm2);
 SFlast = vector(qm1*qm2);
 
 ww1 = vector(qm1);
 zz1 = vector(qm1); 
 ww2 = vector(qm2);
 zz2 = vector(qm2);
 
 dn = (double)pn;
 sdn = sqrt(dn);
 b = -qPHI(truncate/2.)/sdn;  
 gausslegendre(qm1, -b, b, zz1, ww1);
 for (i=0; i<qm1; i++) ww1[i] *= sdn * phi( zz1[i]*sdn, 0.);
 qnspecial1 = qm1 - 1;
 
 ddf = (double)(df);  
 b1 = sqrt(qCHI(     truncate/2., df)/ddf);
 b2 = sqrt(qCHI(1. - truncate/2., df)/ddf);
 gausslegendre(qm2, b1, b2, zz2, ww2); 
 for (j=0; j<qm2; j++) ww2[j] *= 2.*ddf*zz2[j] * chi( ddf*zz2[j]*zz2[j], df); 
 qnspecial2 = qm2 - 1;
 
 Nlocal = qm_for_l_and_c(l, zz2[qnspecial2]*c);
 m = xe2_sfm_deluxe(l, zz2[qnspecial2]*c, hs, q, zz1[qnspecial1]+mu0, zz1[qnspecial1]+mu1, mode, Nlocal, nmax, BOUND, SF, &nsm, &rho);
 if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
 
 if ( nsm < 1 ) { /* did not converge yet -- should be the rare case */
   nn = nmax;
   warning("The geometric tail approximation might not work.");
 } else {
   nstop = nsm;
   
   m = xe2_sfm_deluxe(l, zz2[qnspecial2]*c, hs, q, zz1[qnspecial1-1]+mu0, zz1[qnspecial1-1]+mu1, mode, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
   if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax ) {
       nstop = nstop_;
       i++;
       m = xe2_sfm_deluxe(l, zz2[qnspecial2]*c, hs, q, zz1[qnspecial1-i]+mu0, zz1[qnspecial1-i]+mu1, mode, Nlocal, nmax, BOUND, SF, &nstop_, &rho);
       if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }
   nn = nsm;
 }
   
 for (n=0; n<nmax; n++) p0[n] = 0.; 

 for (i=0; i<qm1; i++) {
   for (j=0; j<qm2; j++) {
     Nlocal = qm_for_l_and_c(l, zz2[j]*c);
     m = xe2_sfm_deluxe(l, zz2[j]*c, hs, q, zz1[i]+mu0, zz1[i]+mu1, mode, Nlocal, nn, BOUND, SF, &nstop, &rho);
     if ( m != 0 ) warning("trouble with internal [package spc] function xe2_sfm_deluxe");
     if ( nstop < 1 ) {
       nstop = nn;
       warning("The geometric tail approximation might not work.");     
     }
     rhomany[i + j*qm1] = rho;
     for (n=0; n<nstop; n++)  p0[n] += ww1[i] * ww2[j] * SF[n];
     if ( nn > nstop ) {
       for (n=nstop; n<nn; n++) p0[n] += ww1[i] * ww2[j] * SF[nstop-1] * pow(rho, n-nstop+1);
     }
     SFlast[i + j*qm1] = SF[nstop-1] * pow(rho, nn-nstop);
   }
 }
 
 sf_level_adj = 1. - p;
 if ( q > 1 ) sf_level_adj *= p0[q-2];
 if ( p0[nn-1] <= sf_level_adj ) { 
   n = nn-1;
   while (  p0[n] <= sf_level_adj &&  n > 0  ) n--;
   if ( p0[n] > sf_level_adj ) Lp = (double)( n + 2 - q + 1 ); else Lp = 1.;
 } else {
   for (n=nn; n<nmax; n++) {
      p0[n] = 0.;
      for (i=0; i<qm1; i++)
        for (j=0; j<qm2; j++)
          p0[n] += ww1[i] * ww2[j] * SFlast[i + j*qm1] * pow(rhomany[i + j*qm1], n-nn+1);
      if (  p0[n] <= sf_level_adj  ) {
        Lp = (double)( n + 1  -  q + 1 );
        n = nmax+1;
      }
   }
 }
 
 Free(p0);
 Free(ww1);
 Free(zz1); 
 Free(ww2);
 Free(zz2);
 Free(SF);
 Free(SFlast);
 Free(rhomany);
 
 return Lp;
}


double xc2_iglad (double k, double h, double mu0, double mu1, int N)
{ double *a, *arl, *psi, rho, ad, norm, z1, z2, z11, z12, z21, z22, w;
  int i1, i2, j1, j2, status, noofit, NN, N3;

 NN = N*N; N3 = NN*N;
 a = matrix(NN,NN);
 arl = vector(NN);
 psi = vector(NN);

 w = 2.*h/(2.*N - 1.);

 for (i1=0;i1<N;i1++) 
   for (j1=0;j1<N;j1++) 
     for (i2=0;i2<N;i2++)
       for (j2=0;j2<N;j2++) {
         z11 = (i2-i1)*w - w/2.  + k; if (i2==0) z11 = -10000.;
         z12 = (i2-i1)*w + w/2.  + k;
         z21 = -2.*k - (j2-j1)*w - w/2.  + k;
         z22 = -2.*k - (j2-j1)*w + w/2.  + k; if (j2==0) z22 = 10000.;
         if (z11<z21) z1 = z21; else z1 = z11;
         if (z12<z22) z2 = z12; else z2 = z22;
         if (z1>z2) a[i1*N3+j1*NN+i2*N+j2] = 0.;
         else a[i1*N3+j1*NN+i2*N+j2] = -PHI(z2,mu1) + PHI(z1,mu1);
         if (i1==i2 && j1==j2) a[i1*N3+j1*NN+i2*N+j2]++;
       }

 for (j1=0;j1<NN;j1++) arl[j1] = 1.;
 LU_solve(a,arl,NN);

 for (i1=0;i1<N;i1++)
   for (j1=0;j1<N;j1++)
     for (i2=0;i2<N;i2++)
       for (j2=0;j2<N;j2++) {
         z11 = (i2-i1)*w - w/2.  + k; if (i2==0) z11 = -10000.;
         z12 = (i2-i1)*w + w/2.  + k;
         z21 = -2.*k - (j2-j1)*w - w/2.  + k;
         z22 = -2.*k - (j2-j1)*w + w/2.  + k; if (j2==0) z22 = 10000.;
         if (z11<z21) z1 = z21; else z1 = z11;
         if (z12<z22) z2 = z12; else z2 = z22;
         if (z1>z2) a[i2*N3+j2*NN+i1*N+j1] = 0.;
         else a[i2*N3+j2*NN+i1*N+j1] = PHI(z2,mu0) - PHI(z1,mu0);
       }

 pmethod(NN,a,&status,&rho,psi,&noofit);

 ad = 0.;
 norm = 0.;
 for (i1=0;i1<N;i1++) 
   for (j1=0;j1<N;j1++) {
     ad += arl[i1*N+j1] * psi[i1*N+j1];
     norm += psi[i1*N+j1];
   }
 ad /= norm;
 rho0 = rho;

 Free(a);
 Free(arl);
 Free(psi);

 return ad;
}


double xc2_be_arlm(double k, double h, double hs1, double hs2, int q, double mu0, double mu1, int N, double *ced)
{ double *a, *arl, *fn, norm, z1, z2, z11, z12, z21, z22, w;
  int i1, i2, j1, j2, n, NN, N3;
  
 NN  = N*N;
 N3  = NN*N;
 a   = matrix(NN,NN);
 arl = vector(NN);
 fn  = matrix(q+1, NN);

 w = 2.*h/(2.*N - 1.);
 
 /* ARL vector */
 for (i1=0;i1<N;i1++) 
   for (j1=0;j1<N;j1++) 
     for (i2=0;i2<N;i2++)
       for (j2=0;j2<N;j2++) {
         z11 = (i2-i1)*w - w/2.  + k; if (i2==0) z11 = -10000.;
         z12 = (i2-i1)*w + w/2.  + k;
         z21 = -2.*k - (j2-j1)*w - w/2.  + k;
         z22 = -2.*k - (j2-j1)*w + w/2.  + k; if (j2==0) z22 = 10000.;
         if (z11<z21) z1 = z21; else z1 = z11;
         if (z12<z22) z2 = z12; else z2 = z22;
         if (z1>z2) a[i1*N3+j1*NN+i2*N+j2] = 0.;
         else a[i1*N3+j1*NN+i2*N+j2] = -PHI(z2,mu1) + PHI(z1,mu1);
         if (i1==i2 && j1==j2) a[i1*N3+j1*NN+i2*N+j2]++;
       }

 for (j1=0;j1<NN;j1++) arl[j1] = 1.;
 LU_solve(a,arl,NN); 
 
 /* q == 1 */
 i1 = (int) ceil(hs1/w - .5);
 j1 = (int) ceil(hs2/w - .5);
 ced[0] = arl[i1*N + j1];

 /* transition matrix */
 for (i1=0;i1<N;i1++)
   for (j1=0;j1<N;j1++)
     for (i2=0;i2<N;i2++)
       for (j2=0;j2<N;j2++) {
         z11 = (i2-i1)*w - w/2.  + k; if (i2==0) z11 = -10000.;
         z12 = (i2-i1)*w + w/2.  + k;
         z21 = -2.*k - (j2-j1)*w - w/2.  + k;
         z22 = -2.*k - (j2-j1)*w + w/2.  + k; if (j2==0) z22 = 10000.;
         if (z11<z21) z1 = z21; else z1 = z11;
         if (z12<z22) z2 = z12; else z2 = z22;
         if (z1>z2) a[i2*N3+j2*NN+i1*N+j1] = 0.;
         else a[i2*N3+j2*NN+i1*N+j1] = PHI(z2,mu0) - PHI(z1,mu0);
       }
 
 /* density sequence for q > 1 */
 for (n=1; n<=q-1; n++) {
   if (n==1) {
     for (i1=0; i1<N; i1++)
       for (j1=0; j1<N; j1++) fn[0*NN + i1*N+j1] = 0.; 
     i1 = (int) ceil(hs1/w - .5);
     j1 = (int) ceil(hs2/w - .5);
     fn[0*NN + i1*N+j1] = 1.;
     /*printf("\n-- 0 --\n\n");
     for (i1=0; i1<N; i1++) {
       for (j1=0; j1<N; j1++) printf("%.3f\t", fn[0*NN + i1*N+j1]);
       printf("\n");
     }
     printf("\n");*/
   }
   
   for (i1=0; i1<N; i1++)
     for (j1=0; j1<N; j1++) {
       fn[n*NN + i1*N+j1] = 0.;
       for (i2=0; i2<N; i2++)
         /*for (j2=0; j2<N; j2++) fn[n*NN + i1*N+j1] += a[i2*N3+j2*NN+i1*N+j1] * fn[(n-1)*NN + i2*N+j2];*/
         for (j2=0; j2<N; j2++) fn[n*NN + i1*N+j1] += a[i1*N3+j1*NN+i2*N+j2] * fn[(n-1)*NN + i2*N+j2];
     }   
  
   /*printf("\n-- %d --\n\n", n);
   for (i1=0; i1<N; i1++) {
     for (j1=0; j1<N; j1++) printf("%.3f\t", fn[n*NN + i1*N+j1]);
     printf("\n");
   }
   printf("\n");*/
  
   ced[n] = 0.;
   norm = 0.;
   for (i1=0; i1<N; i1++) 
     for (j1=0; j1<N; j1++) {
       ced[n] += arl[i1*N+j1] * fn[n*NN + i1*N+j1];
       norm   += fn[n*NN + i1*N+j1];
     }
   ced[n] /= norm;
 }  
  
 Free(fn);
 Free(a);
 Free(arl);

 return 0;
}


/* Richardson extrapolation */
double xc2_igladR (double k, double h, double mu0, double mu1, int r)
{ double *a, *b, ad;
  int i, j, N;

  a = matrix(r,r);
  b = vector(r);

  for (i=0;i<r;i++) {
    N = (int)pow(2.,(double)(i)+1.);
    b[i] = -xc2_iglad(k,h,mu0,mu1,N);
    a[i*r+0] = -1.;
    for (j=0;j<r;j++)
      if (i==0) a[i*r+j] = 1.;
      else a[i*r+j] = pow( 2, -(double)(j+1.)*(double)(i) );
  }

  LU_solve(a,b,r);

  ad = b[0];

  Free(a);
  Free(b);

  return ad;
}


/* variance control charts */

/* -------------- Chebyshev polynomials on [-1,1] ----------------- */

double Tn(double z, int n)
{ double result=1.;

 if ( fabs(z)<1-1e-12 ) {
   switch (n) {
     case 0: result = 1.; break;
     case 1: result = z; break;
     case 2: result = 2.*z*z-1.; break;
     case 3: result = 4.*z*z*z-3.*z; break;
     case 4: result = 8.*pow(z,4.)-8.*z*z+1.; break;
     case 5: result = 16.*pow(z,5.)-20.*z*z*z+5.*z; break;
   }
   if ( n > 5 ) result = cos( (double)(n)*acos(z) );
 }
 else { if ( z<0. && (n % 2 == 1) ) result = -1.; else result = 1.; }
 return result;
}


/* -------------- indefinite integrals of Chebyshev polynomials on [-1,1] ----------------- */ 
double iTn(double z, int n)
{ double result=1.;
 
 switch (n) {
   case 0: result = z; break;
   case 1: result = z*z/2.; break;
   case 2: result = 2.*z*z*z/3. - z; break;
 }
 if ( n > 2 ) result = ( Tn(z,n+1)/(n+1.) - Tn(z,n-1)/(n-1.) )/2.;
 return result;
}


/* -------------- derivatives of Chebyshev polynomials on [-1,1] ----------------- */ 
double dTn(double z, int n)
{ double result=1., dn;
 dn = (double)n;
 if ( fabs(z)<1-1e-12 ) {
   switch (n) {
     case 0: result = 0.; break;
     case 1: result = 1.; break;
     case 2: result = 4.*z; break;
     case 3: result = 12.*z*z-3.; break;
     case 4: result = 32.*z*z*z-16.*z; break;
     case 5: result = 80.*pow(z,4.)-60.*z*z+5.; break;
   }
   if ( n > 5 ) result = dn * ( Tn(z,n-1) - z*Tn(z,n) ) / (1.-z*z);
 }
 else { if ( z<0. && (n % 2 == 0) ) result = -dn*dn; else result = dn*dn; }
 return result;
}


double seU_iglarl(double l, double cu, double hs, double sigma, int df, int N, int qm)
{ double *a, *g, *w, *z, arl, Hij, xi, xl, za, xu, dN, ddf, s2, v;
  int i, j, k;

 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)N;

 a = matrix(N,N);
 g = vector(N);
 w = vector(qm);
 z = vector(qm);

 for (i=0;i<N;i++) {
   xi = cu/2.*(1.+cos(PI*(2.*(i+1.)-1.)/2./dN));

   za = (1.-l)*xi;

   xl = za; 
   xu = cu;
   if ( df!=2 ) { 
     xl = 0.; 
     xu = sqrt(cu-za); 
   }

   gausslegendre(qm,xl,xu,z,w);

   v = (cu - za) / l;
   if (df==2) a[i*N] = exp(-v/s2);
   else       a[i*N] = 1. - CHI( ddf/s2*v, df);

   for (j=1;j<N;j++) {
     Hij = 0.;
     for (k=0;k<qm;k++) {
       v = (z[k] - za) / l;
       if ( df==2 )
         Hij += w[k] * Tn( (2.*z[k]-cu)/cu, j) * exp(-v/s2);
       if ( df!=2 )
         Hij += w[k] * Tn( (2.*(z[k]*z[k]+za)-cu)/cu ,j) * 2. * pow(z[k], ddf-1.) * exp(-ddf*z[k]*z[k]/2./s2/l);
     }
     if (df==2) Hij /= s2*l;
     else       Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf,ddf/2.);
     a[i*N+j] = Tn( (2.*xi-cu)/cu ,j) - Hij;
   }
 }

 for (j=0;j<N;j++) g[j] = 1.;
 LU_solve(a,g,N);

 arl = g[0];
 for (j=1;j<N;j++)
   arl += g[j] * Tn( (2.*hs-cu)/cu ,j);

 Free(z);
 Free(w);
 Free(g);
 Free(a);

 return arl;
}


double stdeU_iglarl(double l, double cu, double hs, double sigma, int df, int N, int qm)
{ double *a, *g, *w, *z, arl, Hij, xi, xl, za, xu, dN, ddf, s2, v;
  int i, j, k;

 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)N;

 a = matrix(N,N);
 g = vector(N);
 w = vector(qm);
 z = vector(qm);

 for (i=0;i<N;i++) {
   xi = cu/2.*(1.+cos(PI*(2.*(i+1.)-1.)/2./dN));

   za = (1.-l)*xi;

   xl = za; 
   xu = cu;

   gausslegendre(qm,xl,xu,z,w);

   v = (cu - za) / l;
   a[i*N] = 1. - CHI( ddf/s2*v*v, df);

   for (j=1;j<N;j++) {
     Hij = 0.;
     for (k=0;k<qm;k++) {
       v = (z[k] - za) / l;
       Hij += w[k] * Tn( (2.*z[k]-cu)/cu ,j) * pow(v,ddf-1.)*exp(-ddf/2./s2*v*v);
     }
     Hij *= 2./l/gammafn(ddf/2.)/pow(2.*s2/ddf,ddf/2.);
     a[i*N+j] = Tn( (2.*xi-cu)/cu ,j) - Hij;
   }
 }

 for (j=0;j<N;j++) g[j] = 1.;
 LU_solve(a,g,N);

 arl = g[0];
 for (j=1;j<N;j++)
   arl += g[j] * Tn( (2.*hs-cu)/cu ,j);

 Free(z);
 Free(w);
 Free(g);
 Free(a);

 return arl;
}


double seU_sf(double l, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0)
{ double *S1s, *S2s, *Pns, *ws, *zs, *zch, *rside, za=0., s2, ddf, xl, xu;
  int i, j, k, n, *ps;

 s2 = sigma*sigma;
 ddf = (double)df;

 S1s = matrix(N,N);
 S2s = matrix(N,N);
 ps = ivector(N);
 zch = vector(N);
 rside = vector(N);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,N);

/* Chebyshev nodes on [0,cu] */
 for (i=0; i<N; i++) zch[i] = cu/2.*(1.+cos(PI*(2.*(i+1.)-1.)/2./(double)N) );

/* P(L>1)(zch[i]) */
 for (i=0; i<N; i++) rside[i] = CHI( ddf/s2*(cu-(1.-l)*zch[i])/l, df);

 for (i=0; i<N; i++) {
   za = (1.-l)*zch[i];
   if ( df==2 ) { xl = za; xu = cu; }
   else         { xl = 0.; xu = sqrt(cu-za); }
   gausslegendre(qm, xl, xu, zs, ws);
   for (j=0; j<N; j++) {
     S1s[i*N+j] = 0.;
     for (k=0; k<qm; k++)
       if ( df==2 )
         S1s[i*N+j] += ws[k]*Tn((2.*zs[k]-cu)/cu, j) * exp((za-zs[k])/s2/l); 
       else
         S1s[i*N+j] += ws[k]*Tn((2.*(zs[k]*zs[k]+za)-cu)/cu, j) * 2.*pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
     if ( df==2 ) S1s[i*N+j] /= s2*l;
     else         S1s[i*N+j] /= gammafn(ddf/2.) * pow(2.*s2*l/ddf,ddf/2.);
   }
 }

 for (i=0; i<N; i++)
   for (j=0; j<N; j++) S2s[i*N+j] = Tn( (2.*zch[i]-cu)/cu, j);

 LU_decompose(S2s, ps, N);

 for (n=1;n<=nmax;n++) {
   if (n==1)
     for (i=0; i<N; i++) {
       Pns[i] = 0.;
       for (j=0; j<N; j++)
         Pns[i] += 2./N * Tn( (2.*zch[j]-cu)/cu, i) * rside[j];
       if ( i==0 ) Pns[i] /= 2.;
     }
   else {
     for (i=0; i<N; i++) {
       rside[i] = 0.;
       for (j=0; j<N; j++) rside[i] += S1s[i*N+j] * Pns[(n-2)*N+j];
     }
     LU_solve2(S2s, rside, ps, N);
     for (i=0; i<N; i++) Pns[(n-1)*N+i] = rside[i];
   }

   p0[n-1] = 0.;
   if ( n==1 )
     p0[0] = CHI(ddf/s2*(cu-(1.-l)*hs)/l, df);
   else
     for (j=0; j<N; j++)
       p0[n-1] += Pns[(n-1)*N+j] * Tn( (2.*hs-cu)/cu, j);
 }

 Free(Pns);
 Free(zs);
 Free(ws);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);

 return 0;
}


double seU_sf_deluxe(double l, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0, int *nstop, double *rho)
{ double *S1s, *S2s, *Pns, *ws, *zs, *zch, *rside, za=0., s2, ddf, xl, xu, mn_minus=1., mn_plus=0., oben, unten, q;
  int i, j, k, n, *ps;

 s2 = sigma*sigma;
 ddf = (double)df;

 S1s = matrix(N,N);
 S2s = matrix(N,N);
 ps = ivector(N);
 zch = vector(N);
 rside = vector(N);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,N);
 
 *nstop = 0;

/* Chebyshev nodes on [0,cu] */
 for (i=0; i<N; i++) zch[i] = cu/2.*(1.+cos(PI*(2.*(i+1.)-1.)/2./(double)N) );

/* P(L>1)(zch[i]) */
 for (i=0; i<N; i++) rside[i] = CHI( ddf/s2*(cu-(1.-l)*zch[i])/l, df);

 for (i=0; i<N; i++) {
   za = (1.-l)*zch[i];
   if ( df==2 ) { xl = za; xu = cu; }
   else         { xl = 0.; xu = sqrt(cu-za); }
   gausslegendre(qm, xl, xu, zs, ws);
   for (j=0; j<N; j++) {
     S1s[i*N+j] = 0.;
     for (k=0; k<qm; k++)
       if ( df==2 )
         S1s[i*N+j] += ws[k]*Tn((2.*zs[k]-cu)/cu, j) * exp((za-zs[k])/s2/l); 
       else
         S1s[i*N+j] += ws[k]*Tn((2.*(zs[k]*zs[k]+za)-cu)/cu, j) * 2.*pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
     if ( df==2 ) S1s[i*N+j] /= s2*l;
     else         S1s[i*N+j] /= gammafn(ddf/2.) * pow(2.*s2*l/ddf,ddf/2.);
   }
 }

 for (i=0; i<N; i++)
   for (j=0; j<N; j++) S2s[i*N+j] = Tn( (2.*zch[i]-cu)/cu, j);

 LU_decompose(S2s, ps, N);

 for (n=1;n<=nmax;n++) {
   if (n==1)
     for (i=0; i<N; i++) {
       Pns[i] = 0.;
       for (j=0; j<N; j++)
         Pns[i] += 2./N * Tn( (2.*zch[j]-cu)/cu, i) * rside[j];
       if ( i==0 ) Pns[i] /= 2.;
     }
   else {
     for (i=0; i<N; i++) {
       rside[i] = 0.;
       for (j=0; j<N; j++) rside[i] += S1s[i*N+j] * Pns[(n-2)*N+j];
     }
     LU_solve2(S2s, rside, ps, N);
     for (i=0; i<N; i++) Pns[(n-1)*N+i] = rside[i];
   }

   p0[n-1] = 0.;
   if ( n==1 )
     p0[0] = CHI(ddf/s2*(cu-(1.-l)*hs)/l, df);
   else
     for (j=0; j<N; j++)
       p0[n-1] += Pns[(n-1)*N+j] * Tn( (2.*hs-cu)/cu, j);
     
   mn_minus = 1.; mn_plus = 0.;
   if ( n > 1 ) {
     for (i=0; i<N; i++) {
       oben = 0.; unten = 0.;
       for (j=0; j<N; j++) {
         oben  += Pns[(n-1)*N+j] * Tn( (2.*zch[i]-cu)/cu, j);
         unten += Pns[(n-2)*N+j] * Tn( (2.*zch[i]-cu)/cu, j);
       }
       if ( fabs(unten)<1e-16 )
         if ( fabs(oben)<1e-16 ) q = 0.;
         else q = 1.;
       else q = oben/unten;

       if ( q<mn_minus ) mn_minus = q;
       if ( q>mn_plus ) mn_plus = q;
     }
     *rho = (mn_minus + mn_plus)/2.;
     if ( fabs(mn_plus - mn_minus) < FINALeps ) {       
       *nstop = n;
       n = nmax + 1;
     }     
   } /* n > 1 */
 } /* n=1; n<=nmax; n++ */ 

 Free(Pns);
 Free(zs);
 Free(ws);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);

 return 0;
}


int choose_N_for_seU(double lambda)
{ int N=20;  

  N = 25;
  if ( 0.1  <= lambda && lambda < 0.2 )  N = 35;
  if ( 0.05 <= lambda && lambda < 0.1 )  N = 50;
  if ( 0.02 <= lambda && lambda < 0.05)  N = 70;
  if ( 0.01 <= lambda && lambda < 0.02)  N = 100;
  if (                   lambda < 0.01 ) N = 150;
  
  return N;  
}


int choose_N_for_se2(double lambda, double cl, double cu)
{ int N=20, M=1;

  M = ceil( ( log(cl) - log(cu) )/log( 1. - lambda ) );

  N = 5;
  if ( 0.1  <= lambda && lambda < 0.2 )  N = 10;
  if ( 0.05 <= lambda && lambda < 0.1 )  N = 20;
  if ( 0.02 <= lambda && lambda < 0.05)  N = 40;
  if ( 0.01 <= lambda && lambda < 0.02)  N = 60;
  if (                   lambda < 0.01 ) N = 90;
  N *= M;
  
  if ( N < 30 ) N = 30;
  if ( N > 200 ) N = 200;
  
  return N;  
}


double seU_sf_prerun_SIGMA_deluxe(double l, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0)
{ double *ww, *zz, b1, b2, ddf2, *SF, rho, s2;
  int i, m, n, nstop, Nlocal;

 Nlocal = choose_N_for_seU(l); 
 
 SF = vector(nmax);
 ww = vector(qm2);
 zz = vector(qm2);
 
 ddf2 = (double)(df2);  
 b1 = qCHI(     truncate/2., df2)/ddf2;
 b2 = qCHI(1. - truncate/2., df2)/ddf2; 
 gausslegendre(qm2, b1, b2, zz, ww); 
 for (i=0; i<qm2; i++) ww[i] *= ddf2 * chi( ddf2*zz[i], df2);
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm2; i++) {    
    s2 = zz[i];
    m = seU_sf_deluxe(l, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop, &rho);
    if ( m != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
    if ( nstop > 0  ) {
      for (n=0; n<nstop; n++)  p0[n] += ww[i] * SF[n];
      for (n=nstop; n<nmax; n++)  p0[n] += ww[i] * SF[nstop-1] * pow(rho, n-nstop+1);
    } else {
      for (n=0; n<nmax; n++)  p0[n] += ww[i] * SF[n];
    }
 }  
 
 Free(ww);
 Free(zz); 
 Free(SF); 
 
 return 0;
}


double seU_sf_prerun_SIGMA(double l, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0)
{ double *ww, *zz, b1, b2, ddf2, *SF, s2;
  int i, m, n, Nlocal;

 Nlocal = choose_N_for_seU(l); 
  
 SF = vector(nmax);
 ww = vector(qm2);
 zz = vector(qm2);
 
 ddf2 = (double)(df2);  
 b1 = qCHI(     truncate/2., df2)/ddf2;
 b2 = qCHI(1. - truncate/2., df2)/ddf2; 
 gausslegendre(qm2, b1, b2, zz, ww); 
 for (i=0; i<qm2; i++) ww[i] *= ddf2 * chi( ddf2*zz[i], df2);
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm2; i++) {
   s2 = zz[i];
   m = seU_sf(l, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF);
   if ( m != 0 ) warning("trouble with internal [package spc] function seU_sf");
   for (n=0; n<nmax; n++)  p0[n] += ww[i] * SF[n];
 }  
 
 Free(ww);
 Free(zz); 
 Free(SF);
 
 return 0;
}


double seU_Wq_prerun_SIGMA_deluxe(double l, double cu, double p, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate)
{ double *ww, *zz, b1, b2, ddf2, *SF, *p0, rho, *rhomany, *SFlast, Lp=-1., sf_level_adj, s2;
  int i, j, n, nstop, nstop_, nsm, nn, qnspecial=0, Nlocal;
  
 Nlocal = choose_N_for_seU(l); 
  
 p0 = vector(nmax);
 SF = vector(nmax); 
 rhomany = vector(qm2);
 SFlast = vector(qm2);
 ww = vector(qm2);
 zz = vector(qm2);
 
 ddf2 = (double)(df2);  
 b1 = qCHI(     truncate/2., df2)/ddf2;
 b2 = qCHI(1. - truncate/2., df2)/ddf2; 
 gausslegendre(qm2, b1, b2, zz, ww); 
 for (i=0; i<qm2; i++) ww[i] *= ddf2 * chi( ddf2*zz[i], df2);
 
 qnspecial = (qm2+1) / 2;
 
 s2 = zz[qnspecial];
 j = seU_sf_deluxe(l, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nsm, &rho);
 if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
 n = nsm;
 
 if ( nsm < 1 ) { /* did not converge yet -- should be the rare case */
   nn = nmax;
   warning("The geometric tail approximation might not work.");
 } else {
   nstop = nsm;

   s2 = zz[qnspecial+1];
   j = seU_sf_deluxe(l, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop_, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax ) {
       nstop = nstop_;
       i++;
       s2 = zz[qnspecial+i];
       j = seU_sf_deluxe(l, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop_, &rho);
       if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }
 
   nstop = n;
   s2 = zz[qnspecial-1];
   j = seU_sf_deluxe(l, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop_, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax ) {
       nstop = nstop_;
       i++;   
       s2 = zz[qnspecial-i];
       j = seU_sf_deluxe(l, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop_, &rho);
       if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }  
   nn = nsm;
 }
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm2; i++) {
   s2 = zz[i];
   j = seU_sf_deluxe(l, s2*cu, s2*hs, sigma, df1, Nlocal, nn, qm1, SF, &nstop, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
   if ( nstop < 1 ) {
     nstop = nn;
     warning("The geometric tail approximation might not work.");     
   }
   rhomany[i] = rho;        
   for (n=0; n<nstop; n++)  p0[n] += ww[i] * SF[n];
   if ( nstop < nn) {
     for (n=nstop; n<nn; n++) p0[n] += ww[i] * SF[nstop-1] * pow(rho, n-nstop+1);          
   }
   SFlast[i] = SF[nstop-1] * pow(rho, nn-nstop);
 }
 
 sf_level_adj = 1.-p;
 if ( p0[nn-1] <= sf_level_adj ) { 
   n = nn-1;
   while (  p0[n] <= sf_level_adj &&  n > 0  ) n--;
   if ( p0[n] > 1.-p ) Lp = (double)( n + 2 ); else Lp = 1.;
 } else {
   for (n=nn; n<nmax; n++) {
      p0[n] = 0.;
      for (i=0; i<qm2; i++) p0[n] += ww[i] * SFlast[i] * pow(rhomany[i], n-nn+1);
      if (  p0[n] <= sf_level_adj  ) {
        Lp = (double)( n + 1 );
        n = nmax+1;
      }
   }
 }
 
 Free(p0);
 Free(ww);
 Free(zz); 
 Free(SF);
 Free(SFlast);
 Free(rhomany);

 return Lp;
}



double seU_iglarl_prerun_SIGMA(double l, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate)
{ double *ww, *zz, b1, b2, result, ddf2, s2;
  int i;
  
  ww = vector(qm2);
  zz = vector(qm2);  
  ddf2 = (double)(df2);  
  b1 = qCHI(     truncate/2., df2)/ddf2;
  b2 = qCHI(1. - truncate/2., df2)/ddf2; 
  gausslegendre(qm2, b1, b2, zz, ww);  
  result = 0.;
  for (i=0; i<qm2; i++) {
    s2 = zz[i];
    result += ww[i] * ddf2 * chi( ddf2*s2, df2) * seU_iglarl(l, s2*cu, s2*hs, sigma, df1, N, qm1);
  }
  Free(ww);
  Free(zz);
  
  return result;
}


double seU_Wq(double l, double cu, double p, double hs, double sigma, int df, int N, int nmax, int qm)
{ double *S1s, *S2s, *Pns, *p0, *ws, *zs, *zch, *rside, za=0., s2, ddf, xl, xu, q_minus=0., q_plus=0., mn_minus=1., mn_plus=0., oben, unten, q, enumerator=0., Wq=0.;
  int i, j, k, n, *ps;
 
 s2 = sigma*sigma;
 ddf = (double)df;

 S1s = matrix(N,N);
 S2s = matrix(N,N);
 ps = ivector(N);
 zch = vector(N);
 rside = vector(N);
 ws  = vector(qm);
 zs  = vector(qm);
 p0 = vector(nmax);
 Pns = matrix(nmax,N);

/* Chebyshev nodes on [0,cu] */
 for (i=0; i<N; i++) zch[i] = cu/2.*(1.+cos(PI*(2.*(i+1.)-1.)/2./(double)N) );

/* P(L>1)(zch[i]) */
 for (i=0; i<N; i++) rside[i] = CHI( ddf/s2*(cu-(1.-l)*zch[i])/l, df);

 for (i=0; i<N; i++) {
   za = (1.-l)*zch[i];
   if ( df==2 ) { xl = za; xu = cu; }
   else         { xl = 0.; xu = sqrt(cu-za); }
   gausslegendre(qm, xl, xu, zs, ws);
   for (j=0; j<N; j++) {
     S1s[i*N+j] = 0.;
     for (k=0; k<qm; k++)
       if ( df==2 )
         S1s[i*N+j] += ws[k]*Tn((2.*zs[k]-cu)/cu, j) * exp((za-zs[k])/s2/l); 
       else
         S1s[i*N+j] += ws[k]*Tn((2.*(zs[k]*zs[k]+za)-cu)/cu, j) * 2.*pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
     if ( df==2 ) S1s[i*N+j] /= s2*l;
     else         S1s[i*N+j] /= gammafn(ddf/2.) * pow(2.*s2*l/ddf,ddf/2.);
   }
 }

 for (i=0; i<N; i++)
   for (j=0; j<N; j++) S2s[i*N+j] = Tn( (2.*zch[i]-cu)/cu, j);

 LU_decompose(S2s, ps, N);

 for (n=1;n<=nmax;n++) {
   if (n==1)
     for (i=0; i<N; i++) {
       Pns[i] = 0.;
       for (j=0; j<N; j++)
         Pns[i] += 2./N * Tn( (2.*zch[j]-cu)/cu, i) * rside[j];
       if ( i==0 ) Pns[i] /= 2.;
     }
   else {
     for (i=0; i<N; i++) {
       rside[i] = 0.;
       for (j=0; j<N; j++) rside[i] += S1s[i*N+j] * Pns[(n-2)*N+j];
     }
     LU_solve2(S2s, rside, ps, N);
     for (i=0; i<N; i++) Pns[(n-1)*N+i] = rside[i];
   }

   p0[n-1] = 0.;
   if ( n==1 )
     p0[0] = CHI(ddf/s2*(cu-(1.-l)*hs)/l, df);
   else
     for (j=0; j<N; j++)
       p0[n-1] += Pns[(n-1)*N+j] * Tn( (2.*hs-cu)/cu, j);
     
   if ( p0[n-1] < 1.-p ) {
     Wq = (double)n;
     n = nmax+1;
   } else {    
     mn_minus = 1.; mn_plus = 0.;
     if ( n > 1 ) {
       for (i=0; i<N; i++) {
         oben = 0.; unten = 0.;
         for (j=0; j<N; j++) {
           oben  += Pns[(n-1)*N+j] * Tn( (2.*zch[i]-cu)/cu, j);
           unten += Pns[(n-2)*N+j] * Tn( (2.*zch[i]-cu)/cu, j);
         }
         if ( fabs(unten)<1e-16 )
           if ( fabs(oben)<1e-16 ) q = 0.;
           else q = 1.;
         else q = oben/unten;
         if ( q<mn_minus ) mn_minus = q;
         if ( q>mn_plus ) mn_plus = q;
       }
       enumerator = log( (1.-p)/p0[n-1] );
       q_minus = (double)n + enumerator/log(mn_minus);
       q_plus  = (double)n + enumerator/log(mn_plus);
       /*if ( fabs( (q_plus-q_minus)/q_minus )<FINALeps ) n = nmax+1;*/
       if ( fabs( ceil(q_plus) - ceil(q_minus) ) < .5 ) {
	 Wq = ceil(q_plus);
	 n = nmax +1;
       }
     } /* n > 1 */
   } /* p0[n-1] >= 1.-p */
 } /* n=1; n<=nmax; n++ */ 

 Free(Pns);
 Free(p0);
 Free(zs);
 Free(ws);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);

 return Wq;
}


double seU_crit(double l, double L0, double hs, double sigma, int df, int N, int qm)
{ double s1, s2, s3, ds, L1, L2, L3, norm;
 
 norm = sqrt(df);
 s2 = hs - .15;
 L2 = 0.;
 do {
   s1 = s2;
   L1 = L2;
   s2 += .2/norm;
   L2 = seU_iglarl(l,s2,hs,sigma,df,N,qm);
 } while ( L2 < L0 );

 do {
   s3 = s1 + (L0-L1)/(L2-L1) * (s2-s1);
   L3 = seU_iglarl(l,s3,hs,sigma,df,N,qm);
   ds = s3-s2; s1 = s2; L1 = L2; s2 = s3; L2 = L3;
 } while ( fabs(L0-L3)>1e-7 && fabs(ds)>1e-9 );

 return s3;
}


double stdeU_crit(double l, double L0, double hs, double sigma, int df, int N, int qm)
{ double s1, s2, s3, ds, L1, L2, L3, norm;
 
 norm = sqrt(df);
 s2 = hs - .15;
 L2 = 0.;
 do {
   s1 = s2;
   L1 = L2;
   s2 += .2/norm;
   L2 = stdeU_iglarl(l,s2,hs,sigma,df,N,qm);
 } while ( L2 < L0 );

 do {
   s3 = s1 + (L0-L1)/(L2-L1) * (s2-s1);
   L3 = stdeU_iglarl(l,s3,hs,sigma,df,N,qm);
   ds = s3-s2; s1 = s2; L1 = L2; s2 = s3; L2 = L3;
 } while ( fabs(L0-L3)>1e-7 && fabs(ds)>1e-9 );

 return s3;
}


double seU_crit_prerun_SIGMA(double l, double L0, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate)
{ double s1, s2, s3, ds, L1=0., L2=0., L3=0.;

 s2 = hs;
 do {
   L1 = L2;
   s2 += .2;
   L2 = seU_iglarl_prerun_SIGMA(l, s2, hs, sigma, df1, df2, N, qm1, qm2, truncate);
 } while ( L2 < L0 );

 s1 = s2 - .2;

 do {
   s3 = s1 + (L0-L1)/(L2-L1) * (s2-s1);
   L3 = seU_iglarl_prerun_SIGMA(l, s3, hs, sigma, df1, df2, N, qm1, qm2, truncate);
   ds = s3-s2; s1 = s2; L1 = L2; s2 = s3; L2 = L3;
 } while ( fabs(L0-L3)>1e-6 && fabs(ds)>1e-9 );

 return s3;
}


double seU_q_crit(double l, int L0, double alpha, double hs, double sigma, int df, int N, int qm, double c_error, double a_error)
{ double s1, s2, s3, ds, p1, p2, p3, *SF;
  int result=1;

 SF  = vector(L0); 
  
 s2 = hs; p2 = 1.;
 do {
   p1 = p2;
   s2 += .2;
   result = seU_sf(l, s2, hs, sigma, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in seU_q_crit [package spc]");
   p2 = 1. - SF[L0-1];
 } while ( p2 > alpha );

 s1 = s2 - .2;

 do {
   s3 = s1 + (alpha - p1)/( p2 - p1 ) * (s2-s1);
   result = seU_sf(l, s3, hs, sigma, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in seU_q_crit [package spc]");
   p3 = 1. - SF[L0-1];
   ds = s3 - s2; s1 = s2; p1 = p2; s2 = s3; p2 = p3;
 } while ( fabs(alpha - p3)>a_error && fabs(ds)>c_error );

 Free(SF);
 
 return s3;
}


double seU_q_crit_prerun_SIGMA(double l, int L0, double alpha, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate, int tail_approx, double c_error, double a_error)
{ double s1, s2, s3, ds, p1, p2, p3, *SF;
  int result=1;  
  
 SF  = vector(L0); 
 
 s2 = seU_q_crit(l, L0, alpha, hs, sigma, df1, N, qm1, c_error, a_error); 
 if ( tail_approx ) result = seU_sf_prerun_SIGMA_deluxe(l, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
 else result = seU_sf_prerun_SIGMA(l, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
 if ( result != 0 ) warning("trouble in seU_q_crit_prerun_SIGMA [package spc]");
 p2 = 1. - SF[L0-1];

 if ( p2 > alpha ) {
   do {
     p1 = p2;
     s2 += .2;
     if ( tail_approx ) result = seU_sf_prerun_SIGMA_deluxe(l, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     else result = seU_sf_prerun_SIGMA(l, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     if ( result != 0 ) warning("trouble in seU_q_crit_prerun_SIGMA [package spc]");
     p2 = 1. - SF[L0-1];
   } while ( p2 > alpha );
   s1 = s2 - .2;
 } else {
   do {
     p1 = p2;
     s2 -= .2;
     if ( tail_approx ) result = seU_sf_prerun_SIGMA_deluxe(l, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     else result = seU_sf_prerun_SIGMA(l, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     if ( result != 0 ) warning("trouble in seU_q_crit_prerun_SIGMA [package spc]");
     p2 = 1. - SF[L0-1];
   } while ( p2 <= alpha && s2 > hs );
   s1 = s2 + .2;
 }
 
 do {
   s3 = s1 + (alpha - p1)/( p2 - p1 ) * (s2-s1);
   if ( tail_approx ) result = seU_sf_prerun_SIGMA_deluxe(l, s3, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
   else result = seU_sf_prerun_SIGMA(l, s3, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
   if ( result != 0 ) warning("trouble in seU_q_crit_prerun_SIGMA [package spc]");
   p3 = 1. - SF[L0-1];
   ds = s3 - s2; s1 = s2; p1 = p2; s2 = s3; p2 = p3;
 } while ( fabs(alpha - p3)>a_error && fabs(ds)>c_error );

 Free(SF);
 
 return s3;
}


double se2_iglarl(double l, double cl, double cu, double hs, double sigma, int df, int N, int qm)
{ double *a, *g, *w, *z, *t, h, arl, Hij, xl, za, dN, ddf, s2,
         t0, t1, x0, x1;
  int i, j, k, qi, qj, M, Ntilde, NN, ii, it, jj;

 M = ceil( (log(cl) - log(cu))/log(1.-l) );
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde - 1.;  

 a = matrix(NN, NN);
 g = vector(NN);
 t = vector(NN);
 w = vector(qm);
 z = vector(qm);

 for(i=0;i<M;i++) {
   t0 = cl/pow(1.-l,(double)(i));
   t1 = t0/(1.-l);
   if (t1>cu) t1 = cu;

   for (j=1;j<Ntilde;j++) { /* node_i,Ntilde-1 = node_i+1,0 */
     h = cos( PI/dN *(dN-j) );
     t[i*(Ntilde-1)+j] = t0 + (h+1.)/2.*(t1-t0);
     /* Chebyshev Gauss-Lobatto nodes on [t0,t1] */
   }
 }
 t[0] = cl;

 for (i=0;i<M;i++) {
   for (j=1;j<=Ntilde;j++) {
     ii = i*Ntilde + j-1;
     it = i*(Ntilde-1) + j-1;

     za = (1.-l)*t[it];
     if (za<cl) xl = cl; else xl = za;

     for (qi=0;qi<i-1;qi++)
       for (qj=1;qj<=Ntilde;qj++) {
         jj = qi*Ntilde + qj-1;
         a[ii*NN+jj] = 0.;
       }

     if (i>0) {
       qi = i-1;
       t0 = cl/pow(1.-l,(double)qi);
       t1 = t0/(1.-l);
       if (t1>cu) t1 = cu;
       if (t0<xl) x0 = xl; else x0 = t0;
       if (df==2)
         x1 = t1;
       else {
         if (x0-za>1e-10) x0 = sqrt(x0-za); else x0 = 0.;
         if (t1-za>1e-10) x1 = sqrt(t1-za); else x1 = 0.;
       }

       for (qj=1;qj<=Ntilde;qj++) {
         jj = qi*Ntilde + qj-1;

         if (j==1) a[ii*NN+jj] = - Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1);
         else {
           if (fabs(t1-x0)>1e-8) {
             gausslegendre(qm,x0,x1,z,w);
             Hij = 0.;
             for (k=0;k<qm;k++) {
               if (df==2)
                 Hij += w[k] * Tn( (2.*z[k]-t0-t1)/(t1-t0), qj-1) * 
                        exp((za-z[k])/s2/l);
               if (df!=2)
                 Hij += w[k] * Tn( (2.*(z[k]*z[k]+za)-t0-t1)/(t1-t0) ,qj-1) *
                        2. * pow(z[k], ddf-1.) * exp(-ddf*z[k]*z[k]/2./s2/l);
             }
             if (df==2) Hij /= s2*l;
             else       Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf,ddf/2.);
             a[ii*NN+jj] = -Hij;
           }
           else a[ii*NN+jj] = 0.;
         }
       }
     }

     for (qi=i;qi<M;qi++) {
       t0 = cl/pow(1.-l,(double)qi);
       t1 = t0/(1.-l);
       if (t1>cu) t1 = cu;
       if (t0<xl) x0 = xl; else x0 = t0;
       if (df==2)
         x1 = t1;
       else {
        if (x0-za>1e-10) x0 = sqrt(x0-za); else x0 = 0.;
        if (t1-za>1e-10) x1 = sqrt(t1-za); else x1 = 0.;
       }

       if (i>0 && j==1 && qi==i) {
         for (qj=1;qj<=Ntilde;qj++) {
           jj = qi*Ntilde + qj-1;
           a[ii*NN+jj] = Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1);
         }
       }

       if (i>0 && j==1 && qi>i) {
         for (qj=1;qj<=Ntilde;qj++) {
           jj = qi*Ntilde + qj-1;
           a[ii*NN+jj] = 0.;
         }
       }

       if (i==0 || j>1) {
         for (qj=1;qj<=Ntilde;qj++) {
           jj = qi*Ntilde + qj-1;
           gausslegendre(qm,x0,x1,z,w);
           Hij = 0.;
           for (k=0;k<qm;k++) {
             if (df==2)
               Hij += w[k] * Tn( (2.*z[k]-t0-t1)/(t1-t0), qj-1) * 
                      exp((za-z[k])/s2/l);
             if (df!=2)
               Hij += w[k] * Tn( (2.*(z[k]*z[k]+za)-t0-t1)/(t1-t0),qj-1) *
                      2. * pow(z[k], ddf-1.) * exp(-ddf*z[k]*z[k]/2./s2/l);
           }
           if (df==2) Hij /= s2*l;
           else       Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf,ddf/2.);
           if (qi==i) a[ii*NN+jj] = Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1) - 
                                        Hij;
           else a[ii*NN+jj] = -Hij;
         }
       }
     }
   }
 }

 for (j=0;j<NN;j++) g[j] = 1.;
 for (j=1;j<M;j++) g[Ntilde*j] = 0.;

 LU_solve(a,g,NN);

 arl = 0.;
 for (i=0;i<M;i++) {
   t0 = cl/pow(1.-l,(double)i);
   t1 = t0/(1.-l);
   if (t1>cu) t1 = cu;

   if (t0<=hs && hs<t1)
     for (j=1;j<=Ntilde;j++) {
        ii = i*Ntilde + j-1;
        arl += g[ii] * Tn((2.*hs-t0-t1)/(t1-t0),j-1);
     }
 }

 Free(z);
 Free(w);
 Free(t);
 Free(g);
 Free(a);

 return arl;
}


double stde2_iglarl(double l, double cl, double cu, double hs, double sigma, int df, int N, int qm)
{ double *a, *g, *w, *z, *t, h, arl, Hij, xl, za, dN, ddf, s2, t0, t1, x0, x1, v;
  int i, j, k, qi, qj, M, Ntilde, NN, ii, it, jj;

 M = ceil( (log(cl) - log(cu))/log(1.-l) );
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde - 1.;  

 a = matrix(NN, NN);
 g = vector(NN);
 t = vector(NN);
 w = vector(qm);
 z = vector(qm);

 for(i=0; i<M; i++) {
   t0 = cl/pow(1.-l,(double)(i));
   t1 = t0/(1.-l);
   if ( t1>cu ) t1 = cu;

   for (j=1; j<Ntilde; j++) { /* node_i,Ntilde-1 = node_i+1,0 */
     h = cos( PI/dN *(dN-j) );
     t[i*(Ntilde-1)+j] = t0 + (h+1.)/2.*(t1-t0); /* Chebyshev Gauss-Lobatto nodes on [t0,t1] */
   }
 }
 t[0] = cl;

 for (i=0; i<M; i++) {
   for (j=1; j<=Ntilde; j++) {
     ii = i*Ntilde + j-1;
     it = i*(Ntilde-1) + j-1;

     za = (1.-l)*t[it];
     if ( za<cl ) xl = cl; else xl = za;

     for (qi=0; qi<i-1; qi++)
       for (qj=1; qj<=Ntilde; qj++) {
         jj = qi*Ntilde + qj-1;
         a[ii*NN+jj] = 0.;
       }

     if ( i>0 ) {
       qi = i-1;
       t0 = cl/pow(1.-l,(double)qi);
       t1 = t0/(1.-l);
       if ( t1>cu ) t1 = cu;
       if ( t0<xl ) x0 = xl; else x0 = t0;
       x1 = t1;

       for (qj=1; qj<=Ntilde; qj++) {
         jj = qi*Ntilde + qj-1;

         if ( j==1 ) a[ii*NN+jj] = - Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1);
         else {
           if ( fabs(t1-x0)>1e-8 ) {
             gausslegendre(qm, x0, x1, z, w);
             Hij = 0.;
             for (k=0; k<qm; k++) {
	       v = (z[k] - za) / l;
               Hij += w[k] * Tn( (2.*z[k]-t0-t1)/(t1-t0), qj-1) * pow(v,ddf-1.)*exp(-ddf/2./s2*v*v);
             }
             Hij *= 2./l/gammafn(ddf/2.)/pow(2.*s2/ddf,ddf/2.);
             a[ii*NN+jj] = -Hij;
           }
           else a[ii*NN+jj] = 0.;
         }
       }
     }

     for (qi=i; qi<M; qi++) {
       t0 = cl/pow(1.-l,(double)qi);
       t1 = t0/(1.-l);
       if ( t1>cu ) t1 = cu;
       if ( t0<xl ) x0 = xl; else x0 = t0;
       x1 = t1;

       if ( i>0 && j==1 && qi==i ) {
         for (qj=1; qj<=Ntilde; qj++) {
           jj = qi*Ntilde + qj-1;
           a[ii*NN+jj] = Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1);
         }
       }

       if ( i>0 && j==1 && qi>i ) {
         for (qj=1;qj<=Ntilde;qj++) {
           jj = qi*Ntilde + qj-1;
           a[ii*NN+jj] = 0.;
         }
       }

       if ( i==0 || j>1 ) {
         for ( qj=1; qj<=Ntilde; qj++) {
           jj = qi*Ntilde + qj-1;
           gausslegendre(qm, x0, x1, z, w);
           Hij = 0.;
           for (k=0; k<qm; k++) {
	     v = (z[k] - za) / l;
             Hij += w[k] * Tn( (2.*z[k]-t0-t1)/(t1-t0), qj-1) * pow(v,ddf-1.)*exp(-ddf/2./s2*v*v);
           }
           Hij *= 2./l/gammafn(ddf/2.)/pow(2.*s2/ddf,ddf/2.);
           if ( qi==i ) a[ii*NN+jj] = Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1) - Hij;
           else a[ii*NN+jj] = -Hij;
         }
       }
     }
   }
 }

 for (j=0; j<NN; j++) g[j] = 1.;
 for (j=1; j<M; j++) g[Ntilde*j] = 0.;

 LU_solve(a,g,NN);

 arl = 0.;
 for (i=0; i<M; i++) {
   t0 = cl/pow(1.-l,(double)i);
   t1 = t0/(1.-l);
   if ( t1>cu ) t1 = cu;

   if ( t0<=hs && hs<t1 )
     for (j=1; j<=Ntilde; j++) {
        ii = i*Ntilde + j-1;
        arl += g[ii] * Tn((2.*hs-t0-t1)/(t1-t0),j-1);
     }
 }

 Free(z);
 Free(w);
 Free(t);
 Free(g);
 Free(a);

 return arl;
}


double se2_sf(double l, double cl, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0)
{ double *S1s, *S2s, *Pns, *ws, *zs, *zch, *rside, *b, za=0., s2, ddf, xl, xu, dN, Hij;
  int i, j, k, n, *ps, Ntilde, ihs, M, NN, ii, jj;

 M = ceil( (log(cl) - log(cu))/log(1.-l) );
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde;
 
 ihs = floor( (log(cl) - log(hs))/log(1.-l) );
 if ( ihs<0 ) ihs = 0;

 S1s = matrix(NN,NN);
 S2s = matrix(NN,NN);
 ps = ivector(NN);
 zch = matrix(M,Ntilde);
 rside = vector(NN);
 b   = vector(M+1);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,NN);

/* interval borders b_i = cl/(1-l)^i */
 for (i=0; i<M; i++) b[i] = cl/pow(1.-l, (double)(i));
 b[M] = cu;

 /* Chebyshev nodes on [b_0,b_1],[b_1,b_2],...,[b_M-1,cu] */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     zch[ i*Ntilde+j ] = b[i] + (b[i+1]-b[i])/2.*(1.+cos(PI*(2.*j+1.)/2./dN));
   }

 /* P(L>1)(zch[i,j]) */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     rside[ i*Ntilde+j ] =  CHI( ddf/s2*(cu-(1.-l)*zch[ i*Ntilde+j ])/l, df)
                          - CHI( ddf/s2*(cl-(1.-l)*zch[ i*Ntilde+j ])/l, df);
   }

 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     za = (1.-l)*zch[ i*Ntilde+j ];
     for (ii=0; ii<M; ii++)
       for (jj=0; jj<Ntilde; jj++) {
         if ( b[ii+1]<za ) S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = 0.;
         else {
           if ( za<b[ii] ) xl = b[ii]; else xl = za;
           xu = b[ii+1];
           if ( df!=2 ) {
             xl = sqrt(xl-za);
             xu = sqrt(xu-za);
           }
           gausslegendre(qm, xl, xu, zs, ws);
           Hij = 0.;
           for (k=0; k<qm; k++)
             if ( df==2 )
               Hij += ws[k]*Tn( (2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj) * exp(-zs[k]/s2/l);
             else
               Hij += ws[k]*Tn( (2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                      * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
           if ( df==2 ) Hij *= exp(za/s2/l)/s2/l;
           else         Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf, ddf/2.);
           S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = Hij;
         }
       }
   }

 for (i=0; i<NN; i++)
   for (j=0; j<NN; j++) S2s[i*NN+j] = 0.;

 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++)
     for (jj=0; jj<Ntilde; jj++)
       S2s[ (i*Ntilde+j)*NN + i*Ntilde+jj ] = Tn( (2.*zch[ i*Ntilde+j ]-b[i+1]-b[i])/(b[i+1]-b[i]), jj); 

 LU_decompose(S2s,ps,NN);

 for (n=1; n<=nmax; n++) {
   if ( n==1)
     for (i=0; i<M; i++)
       for (j=0; j<Ntilde; j++) {
         Pns[ i*Ntilde+j ] = 0.;
         for (jj=0; jj<Ntilde; jj++)
           Pns[ i*Ntilde+j ] += /*  usual Chebyshev approximation  */
             2./Ntilde * Tn( (2.*zch[i*Ntilde+jj]-b[i+1]-b[i])/(b[i+1]-b[i]), j) * rside[ i*Ntilde+jj ];
         if ( j==0 ) Pns[ i*Ntilde+j ] /= 2.;
       }
   else {
     for (i=0; i<NN; i++) {
       rside[i] = 0.;
       for (j=0; j<NN; j++) rside[i] += S1s[ i*NN+j ] * Pns[ (n-2)*NN+j ];
     }
     LU_solve2(S2s, rside, ps, NN);
     for (i=0; i<NN; i++) Pns[ (n-1)*NN+i ] = rside[i];
   }

   p0[n-1] = 0.;
   if ( n==1 )
     p0[0] =  CHI( ddf/s2*(cu-(1.-l)*hs)/l, df)
            - CHI( ddf/s2*(cl-(1.-l)*hs)/l, df);
   else
     for (j=0; j<Ntilde; j++)
       p0[n-1] += Pns[ (n-1)*NN + ihs*Ntilde+j ] * Tn( (2.*hs-b[ihs+1]-b[ihs])/(b[ihs+1]-b[ihs]), j);
 }

 Free(Pns);
 Free(zs);
 Free(ws);
 Free(b);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);

 return 0;
}


double se2_sf_deluxe(double l, double cl, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0, int *nstop, double *rho)
{ double *S1s, *S2s, *Pns, *ws, *zs, *zch, *rside, *b, za=0., s2, ddf, xl, xu, dN, Hij, mn_minus=1., mn_plus=0., oben, unten, q;
  int i, j, k, n, *ps, Ntilde, ihs, M, NN, ii, jj;  

 M = ceil( (log(cl)-log(cu))/log(1.-l) );
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde;
 
 ihs = floor( (log(cl) - log(hs))/log(1.-l) );
 if ( ihs<0 ) ihs = 0;

 S1s = matrix(NN,NN);
 S2s = matrix(NN,NN);
 ps = ivector(NN);
 zch = matrix(M,Ntilde);
 rside = vector(NN);
 b   = vector(M+1);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,NN);
 
 *nstop = 0;

/* interval borders b_i = cl/(1-l)^i */
 for (i=0; i<M; i++) b[i] = cl/pow(1.-l, (double)(i));
 b[M] = cu;

 /* Chebyshev nodes on [b_0,b_1],[b_1,b_2],...,[b_M-1,cu] */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     zch[ i*Ntilde+j ] = b[i] + (b[i+1]-b[i])/2.*(1.+cos(PI*(2.*j+1.)/2./dN));
   }

 /* P(L>1)(zch[i,j]) */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     rside[ i*Ntilde+j ] =  CHI( ddf/s2*(cu-(1.-l)*zch[ i*Ntilde+j ])/l, df)
                          - CHI( ddf/s2*(cl-(1.-l)*zch[ i*Ntilde+j ])/l, df);
   }

 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     za = (1.-l)*zch[ i*Ntilde+j ];
     for (ii=0; ii<M; ii++)
       for (jj=0; jj<Ntilde; jj++) {
         if ( b[ii+1]<za ) S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = 0.;
         else {
           if ( za<b[ii] ) xl = b[ii]; else xl = za;
           xu = b[ii+1];
           if ( df!=2 ) {
             xl = sqrt(xl-za);
             xu = sqrt(xu-za);
           }
           gausslegendre(qm, xl, xu, zs, ws);
           Hij = 0.;
           for (k=0; k<qm; k++)
             if ( df==2 )
               Hij += ws[k]*Tn( (2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj) * exp(-zs[k]/s2/l);
             else
               Hij += ws[k]*Tn( (2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                      * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
           if ( df==2 ) Hij *= exp(za/s2/l)/s2/l;
           else         Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf, ddf/2.);
           S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = Hij;
         }
       }
   }

 for (i=0; i<NN; i++)
   for (j=0; j<NN; j++) S2s[i*NN+j] = 0.;

 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++)
     for (jj=0; jj<Ntilde; jj++)
       S2s[ (i*Ntilde+j)*NN + i*Ntilde+jj ] = Tn( (2.*zch[ i*Ntilde+j ]-b[i+1]-b[i])/(b[i+1]-b[i]), jj); 

 LU_decompose(S2s,ps,NN);

 for (n=1; n<=nmax; n++) {
   if ( n==1)
     for (i=0; i<M; i++)
       for (j=0; j<Ntilde; j++) {
         Pns[ i*Ntilde+j ] = 0.;
         for (jj=0; jj<Ntilde; jj++)
           Pns[ i*Ntilde+j ] += /*  usual Chebyshev approximation  */
             2./Ntilde * Tn( (2.*zch[i*Ntilde+jj]-b[i+1]-b[i])/(b[i+1]-b[i]), j) * rside[ i*Ntilde+jj ];
         if ( j==0 ) Pns[ i*Ntilde+j ] /= 2.;
       }
   else {
     for (i=0; i<NN; i++) {
       rside[i] = 0.;
       for (j=0; j<NN; j++) rside[i] += S1s[ i*NN+j ] * Pns[ (n-2)*NN+j ];
     }
     LU_solve2(S2s, rside, ps, NN);
     for (i=0; i<NN; i++) Pns[ (n-1)*NN+i ] = rside[i];
   }

   p0[n-1] = 0.;
   if ( n==1 )
     p0[0] =  CHI( ddf/s2*(cu-(1.-l)*hs)/l, df)
            - CHI( ddf/s2*(cl-(1.-l)*hs)/l, df);
   else
     for (j=0; j<Ntilde; j++)
       p0[n-1] += Pns[ (n-1)*NN + ihs*Ntilde+j ] * Tn( (2.*hs-b[ihs+1]-b[ihs])/(b[ihs+1]-b[ihs]), j);
     
   mn_minus = 1.; mn_plus = 0.;
   if ( n > 1) {
     for (i=0; i<M; i++)
       for (j=0; j<Ntilde; j++) {
         oben = 0.;
         unten = 0.;
         for (jj=0; jj<Ntilde; jj++) {
           oben += Pns[ (n-1)*NN + i*Ntilde+jj ] * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
           unten+= Pns[ (n-2)*NN + i*Ntilde+jj ] * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
         }
         if ( fabs(unten)<1e-16 )
           if ( fabs(oben)<1e-16 ) q = 0.;
           else q = 1.;
         else q = oben/unten;
         if ( q<mn_minus ) mn_minus = q;
         if ( q>mn_plus ) mn_plus = q;
       }
     *rho = (mn_minus + mn_plus)/2.;
     if ( fabs(mn_plus - mn_minus) < FINALeps ) {       
       *nstop = n;
       n = nmax + 1;
     }     
   } /* n > 1 */
 } /* n=1; n<=nmax; n++ */ 


 Free(Pns);
 Free(zs);
 Free(ws);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);

 return 0;
}


double se2_sf_prerun_SIGMA_deluxe(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0)
{ double *ww, *zz, b1, b2, ddf2, *SF, rho, s2;
  int i, m, n, nstop, Nlocal;

 Nlocal = choose_N_for_se2(l, cl, cu); 
  
 SF = vector(nmax);
 ww = vector(qm2);
 zz = vector(qm2);
 
 ddf2 = (double)(df2);  
 b1 = qCHI(     truncate/2., df2)/ddf2;
 b2 = qCHI(1. - truncate/2., df2)/ddf2; 
 gausslegendre(qm2, b1, b2, zz, ww); 
 for (i=0; i<qm2; i++) ww[i] *= ddf2 * chi( ddf2*zz[i], df2);
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm2; i++) {    
    s2 = zz[i];
    m = se2_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop, &rho);
    if ( m != 0 ) warning("trouble with internal [package spc] function se2_sf_deluxe");
    if ( nstop > 0 ) {
      for (n=0; n<nstop; n++)  p0[n] += ww[i] * SF[n];
      for (n=nstop; n<nmax; n++)  p0[n] += ww[i] * SF[nstop-1] * pow(rho, n-nstop+1);
    } else {
      for (n=0; n<nmax; n++)  p0[n] += ww[i] * SF[n];
    }
 }  
 
 Free(ww);
 Free(zz); 
 Free(SF);
 
 return 0;
}


double se2_sf_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0)
{ double *ww, *zz, b1, b2, ddf2, *SF, s2;
  int i, m, n, Nlocal;

 Nlocal = choose_N_for_se2(l, cl, cu);
  
 SF = vector(nmax);
 ww = vector(qm2);
 zz = vector(qm2);
 
 ddf2 = (double)(df2);  
 b1 = qCHI(     truncate/2., df2)/ddf2;
 b2 = qCHI(1. - truncate/2., df2)/ddf2; 
 gausslegendre(qm2, b1, b2, zz, ww); 
 for (i=0; i<qm2; i++) ww[i] *= ddf2 * chi( ddf2*zz[i], df2);
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm2; i++) {
   s2 = zz[i];
   m = se2_sf(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF);
   if ( m != 0 ) warning("trouble with internal [package spc] function se2_sf");
   for (n=0; n<nmax; n++)  p0[n] += ww[i] * SF[n];
 }  
 
 Free(ww);
 Free(zz); 
 Free(SF);
 
 return 0;
}


double se2_Wq_prerun_SIGMA_deluxe(double l, double cl, double cu, double p, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate)
{ double *ww, *zz, b1, b2, ddf2, *SF, *p0, rho, *rhomany, *SFlast, Lp=-1., sf_level_adj, s2;
  int i, j, n, nstop, nstop_, nsm, nn, qnspecial=0, Nlocal;
  
 Nlocal = choose_N_for_se2(l, cl, cu); 
  
 p0 = vector(nmax);
 SF = vector(nmax); 
 rhomany = vector(qm2);
 SFlast = vector(qm2);
 ww = vector(qm2);
 zz = vector(qm2);
 
 ddf2 = (double)(df2);  
 b1 = qCHI(     truncate/2., df2)/ddf2;
 b2 = qCHI(1. - truncate/2., df2)/ddf2; 
 gausslegendre(qm2, b1, b2, zz, ww); 
 for (i=0; i<qm2; i++) ww[i] *= ddf2 * chi( ddf2*zz[i], df2);
 
 qnspecial = (qm2+1) / 2;
 
 s2 = zz[qnspecial];
 j = se2_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nsm, &rho);
 if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
 n = nsm;
 
 if ( nsm < 1 ) { /* did not converge yet -- should be the rare case */
   nn = nmax;
   warning("The geometric tail approximation might not work.");
 } else {
   nstop = nsm;

   s2 = zz[qnspecial+1];
   j = se2_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop_, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax ) {
       nstop = nstop_;
       i++;
       s2 = zz[qnspecial+i];
       j = se2_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop_, &rho);
       if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }
 
   nstop = n;
   s2 = zz[qnspecial-1];
   j = se2_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop_, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax ) {
       nstop = nstop_;
       i++;   
       s2 = zz[qnspecial-i];
       j = se2_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop_, &rho);
       if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }  
   nn = nsm;
 }
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm2; i++) {
   s2 = zz[i];
   j = se2_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nn, qm1, SF, &nstop, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
   if ( nstop < 1 ) {
     nstop = nn;
     warning("The geometric tail approximation might not work.");     
   }
   rhomany[i] = rho;        
   for (n=0; n<nstop; n++)  p0[n] += ww[i] * SF[n];
   if ( nstop < nn) {
     for (n=nstop; n<nn; n++) p0[n] += ww[i] * SF[nstop-1] * pow(rho, n-nstop+1);          
   }
   SFlast[i] = SF[nstop-1] * pow(rho, nn-nstop);
 }
 
 sf_level_adj = 1.-p;
 if ( p0[nn-1] <= sf_level_adj ) { 
   n = nn-1;
   while (  p0[n] <= sf_level_adj &&  n > 0  ) n--;
   if ( p0[n] > 1.-p ) Lp = (double)( n + 2 ); else Lp = 1.;
 } else {
   for (n=nn; n<nmax; n++) {
      p0[n] = 0.;
      for (i=0; i<qm2; i++) p0[n] += ww[i] * SFlast[i] * pow(rhomany[i], n-nn+1);
      if (  p0[n] <= sf_level_adj  ) {
        Lp = (double)( n + 1 );
        n = nmax+1;
      }
   }
 }
 
 Free(p0);
 Free(ww);
 Free(zz); 
 Free(SF);
 Free(SFlast);
 Free(rhomany);

 return Lp;
}


double se2_iglarl_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate)
{ double *ww, *zz, b1, b2, result, ddf2, s2;
  int i;
 
  ww = vector(qm2);
  zz = vector(qm2);  
  ddf2 = (double)(df2);  
  b1 = qCHI(     truncate/2., df2)/ddf2;
  b2 = qCHI(1. - truncate/2., df2)/ddf2; 
  gausslegendre(qm2, b1, b2, zz, ww);  
  result = 0.;
  for (i=0; i<qm2; i++) {
    s2 = zz[i];
    result += ww[i] * ddf2 * chi( ddf2*s2, df2) * se2_iglarl(l, s2*cl, s2*cu, s2*hs, sigma, df1, N, qm1);
  }
  Free(ww);
  Free(zz);
  
  return result;
}


double se2_Wq(double l, double cl, double cu, double p, double hs, double sigma, int df, int N, int nmax, int qm)
{ double *S1s, *S2s, *Pns, *p0, *ws, *zs, *zch, *rside, *b, za=0., s2, ddf, xl, xu, q_minus=0., q_plus=0., dN, Hij, mn_minus=1., mn_plus=0., oben, unten, q, enumerator=0., Wq=0.;
  int i, j, k, n, *ps, Ntilde, ihs, M, NN, ii, jj;  

 M = ceil( (log(cl)-log(cu))/log(1.-l) );
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde;
 
 ihs = floor( (log(cl) - log(hs))/log(1.-l) );
 if ( ihs<0 ) ihs = 0;

 S1s = matrix(NN,NN);
 S2s = matrix(NN,NN);
 ps = ivector(NN);
 zch = matrix(M,Ntilde);
 rside = vector(NN);
 b   = vector(M+1);
 ws  = vector(qm);
 zs  = vector(qm);
 p0 = vector(nmax);
 Pns = matrix(nmax,NN);

/* interval borders b_i = cl/(1-l)^i */
 for (i=0; i<M; i++) b[i] = cl/pow(1.-l, (double)(i));
 b[M] = cu;

 /* Chebyshev nodes on [b_0,b_1],[b_1,b_2],...,[b_M-1,cu] */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     zch[ i*Ntilde+j ] = b[i] + (b[i+1]-b[i])/2.*(1.+cos(PI*(2.*j+1.)/2./dN));
   }

 /* P(L>1)(zch[i,j]) */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     rside[ i*Ntilde+j ] =  CHI( ddf/s2*(cu-(1.-l)*zch[ i*Ntilde+j ])/l, df)
                          - CHI( ddf/s2*(cl-(1.-l)*zch[ i*Ntilde+j ])/l, df);
   }

 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     za = (1.-l)*zch[ i*Ntilde+j ];
     for (ii=0; ii<M; ii++)
       for (jj=0; jj<Ntilde; jj++) {
         if ( b[ii+1]<za ) S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = 0.;
         else {
           if ( za<b[ii] ) xl = b[ii]; else xl = za;
           xu = b[ii+1];
           if ( df!=2 ) {
             xl = sqrt(xl-za);
             xu = sqrt(xu-za);
           }
           gausslegendre(qm, xl, xu, zs, ws);
           Hij = 0.;
           for (k=0; k<qm; k++)
             if ( df==2 )
               Hij += ws[k]*Tn( (2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj) * exp(-zs[k]/s2/l);
             else
               Hij += ws[k]*Tn( (2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                      * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
           if ( df==2 ) Hij *= exp(za/s2/l)/s2/l;
           else         Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf, ddf/2.);
           S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = Hij;
         }
       }
   }

 for (i=0; i<NN; i++)
   for (j=0; j<NN; j++) S2s[i*NN+j] = 0.;

 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++)
     for (jj=0; jj<Ntilde; jj++)
       S2s[ (i*Ntilde+j)*NN + i*Ntilde+jj ] = Tn( (2.*zch[ i*Ntilde+j ]-b[i+1]-b[i])/(b[i+1]-b[i]), jj); 

 LU_decompose(S2s,ps,NN);

 for (n=1; n<=nmax; n++) {
   if ( n==1)
     for (i=0; i<M; i++)
       for (j=0; j<Ntilde; j++) {
         Pns[ i*Ntilde+j ] = 0.;
         for (jj=0; jj<Ntilde; jj++)
           Pns[ i*Ntilde+j ] += /*  usual Chebyshev approximation  */
             2./Ntilde * Tn( (2.*zch[i*Ntilde+jj]-b[i+1]-b[i])/(b[i+1]-b[i]), j) * rside[ i*Ntilde+jj ];
         if ( j==0 ) Pns[ i*Ntilde+j ] /= 2.;
       }
   else {
     for (i=0; i<NN; i++) {
       rside[i] = 0.;
       for (j=0; j<NN; j++) rside[i] += S1s[ i*NN+j ] * Pns[ (n-2)*NN+j ];
     }
     LU_solve2(S2s, rside, ps, NN);
     for (i=0; i<NN; i++) Pns[ (n-1)*NN+i ] = rside[i];
   }

   p0[n-1] = 0.;
   if ( n==1 )
     p0[0] =  CHI( ddf/s2*(cu-(1.-l)*hs)/l, df)
            - CHI( ddf/s2*(cl-(1.-l)*hs)/l, df);
   else
     for (j=0; j<Ntilde; j++)
       p0[n-1] += Pns[ (n-1)*NN + ihs*Ntilde+j ] * Tn( (2.*hs-b[ihs+1]-b[ihs])/(b[ihs+1]-b[ihs]), j);
     
   if ( p0[n-1] < 1.-p ) {
     Wq = (double)n;
     n = nmax+1;
   } else {     
     mn_minus = 1.; mn_plus = 0.;
     if ( n > 1) {
       for (i=0; i<M; i++)
         for (j=0; j<Ntilde; j++) {
           oben = 0.;
           unten = 0.;
           for (jj=0; jj<Ntilde; jj++) {
             oben += Pns[ (n-1)*NN + i*Ntilde+jj ] * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
             unten+= Pns[ (n-2)*NN + i*Ntilde+jj ] * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
           }
           if ( fabs(unten)<1e-16 )
             if ( fabs(oben)<1e-16 ) q = 0.;
             else q = 1.;
           else q = oben/unten;
           if ( q<mn_minus ) mn_minus = q;
           if ( q>mn_plus ) mn_plus = q;
         }
       enumerator = log( (1.-p)/p0[n-1] );
       q_minus = (double)n + enumerator/log(mn_minus);
       q_plus  = (double)n + enumerator/log(mn_plus); 
       /*if ( fabs( (q_plus-q_minus)/q_minus )<FINALeps ) n = nmax+1;*/
       if ( fabs( ceil(q_plus) - ceil(q_minus) ) < .5 ) {
	 Wq = ceil(q_plus);
	 n = nmax +1;
       }
     } /* n > 1 */
   } /* p0[n-1] >= 1.-p */
 } /* n=1; n<=nmax; n++ */ 


 Free(Pns);
 Free(p0);
 Free(zs);
 Free(ws);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);

 return Wq;
}


double se2lu_crit(double l, double L0, double cl, double hs, double sigma, int df, int N, int qm)
{ double s1, s2, s3, ds, L1, L2, L3;

 s2 = hs;
 do {
   s2 += .2;
   L2 = se2_iglarl(l,cl,s2,hs,sigma,df,N,qm);
 } while ( L2 < L0 );

 s1 = s2 - .2;
 L1 = se2_iglarl(l,cl,s1,hs,sigma,df,N,qm);

 do {
   s3 = s1 + (L0-L1)/(L2-L1) * (s2-s1);
   L3 = se2_iglarl(l,cl,s3,hs,sigma,df,N,qm);
   ds = s3-s2; s1 = s2; L1 = L2; s2 = s3; L2 = L3;
 } while ( fabs(L0-L3)>1e-6 && fabs(ds)>1e-9 );

 return s3;
}


double stde2lu_crit(double l, double L0, double cl, double hs, double sigma, int df, int N, int qm)
{ double s1, s2, s3, ds, L1, L2, L3;

 s2 = hs;
 L2 = 0.;
 do {
   s1 = s2;
   L1 = L2;
   s2 += .2;
   L2 = stde2_iglarl(l, cl, s2, hs, sigma, df, N, qm);
 } while ( L2 < L0 );

 do {
   s3 = s1 + (L0-L1)/(L2-L1) * (s2-s1);
   L3 = stde2_iglarl(l, cl, s3, hs, sigma, df, N, qm);
   ds = s3-s2; s1 = s2; L1 = L2; s2 = s3; L2 = L3;
 } while ( fabs(L0-L3)>1e-7 && fabs(ds)>1e-9 );

 return s3;
}


double se2lu_crit_prerun_SIGMA(double l, double L0, double cl, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate)
{ double s1, s2, s3, ds, L1=0., L2=0., L3=0.;

 s2 = hs;
 do {
   L1 = L2;
   s2 += .2;
   L2 = se2_iglarl_prerun_SIGMA(l, cl, s2, hs, sigma, df1, df2, N, qm1, qm2, truncate);
 } while ( L2 < L0 );

 s1 = s2 - .2;

 do {
   s3 = s1 + (L0-L1)/(L2-L1) * (s2-s1);
   L3 = se2_iglarl_prerun_SIGMA(l, cl, s3, hs, sigma, df1, df2, N, qm1, qm2, truncate);
   ds = s3 - s2; s1 = s2; L1 = L2; s2 = s3; L2 = L3;
 } while ( fabs(L0-L3)>1e-6 && fabs(ds)>1e-9 );

 return s3;
}


double se2lu_q_crit(double l, int L0, double alpha, double cl, double hs, double sigma, int df, int N, int qm, double c_error, double a_error)
{ double s1, s2, s3, ds, p1, p2, p3, *SF;
  int result=1;

 SF  = vector(L0); 
  
 s2 = hs; p2 = 1.;
 do {
   p1 = p2;
   s2 += .2;
   result = se2_sf(l, cl, s2, hs, sigma, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in se2lu_q_crit [package spc]");
   p2 = 1. - SF[L0-1];
 } while ( p2 > alpha );

 s1 = s2 - .2;

 do {
   s3 = s1 + (alpha - p1)/( p2 - p1 ) * (s2-s1);
   result = se2_sf(l, cl, s3, hs, sigma, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in se2lu_q_crit [package spc]");
   p3 = 1. - SF[L0-1];
   ds = s3 - s2; s1 = s2; p1 = p2; s2 = s3; p2 = p3;
 } while ( fabs(alpha - p3)>a_error && fabs(ds)>c_error );

 Free(SF);
 
 return s3;
}


double se2lu_q_crit_prerun_SIGMA(double l, int L0, double alpha, double cl, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate, int tail_approx, double c_error, double a_error)
{ double s1, s2, s3, ds, p1, p2, p3, *SF;
  int result=1;

 SF  = vector(L0);
 
 s2 = se2lu_q_crit(l, L0, alpha, cl, hs, sigma, df1, N, qm1, c_error, a_error);                     
 if ( tail_approx ) result = se2_sf_prerun_SIGMA_deluxe(l, cl, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
 else result = se2_sf_prerun_SIGMA(l, cl, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
 if ( result != 0 ) warning("trouble in se2lu_q_crit_prerun_SIGMA [package spc]");
 p2 = 1. - SF[L0-1];
 
 if ( p2 > alpha ) {
   do {
     p1 = p2;
     s2 += .2;
     if ( tail_approx ) result = se2_sf_prerun_SIGMA_deluxe(l, cl, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     else result = se2_sf_prerun_SIGMA(l, cl, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     if ( result != 0 ) warning("trouble in se2lu_q_crit_prerun_SIGMA [package spc]");
     p2 = 1. - SF[L0-1];
   } while ( p2 > alpha );
   s1 = s2 - .2;
 } else {
   do {
     p1 = p2;
     s2 -= .2;
     if ( tail_approx ) result = se2_sf_prerun_SIGMA_deluxe(l, cl, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     else result = se2_sf_prerun_SIGMA(l, cl, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     if ( result != 0 ) warning("trouble in se2lu_q_crit_prerun_SIGMA [package spc]");
     p2 = 1. - SF[L0-1];
   } while ( p2 <= alpha && s2 > hs );
   s1 = s2 + .2;
 }
 
 do {
   s3 = s1 + (alpha - p1)/( p2 - p1 ) * (s2-s1);
   if ( tail_approx ) result = se2_sf_prerun_SIGMA_deluxe(l, cl, s3, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
   else result = se2_sf_prerun_SIGMA(l, cl, s3, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
   if ( result != 0 ) warning("trouble in se2lu_q_crit_prerun_SIGMA [package spc]");
   p3 = 1. - SF[L0-1];
   ds = s3 - s2; s1 = s2; p1 = p2; s2 = s3; p2 = p3;
 } while ( fabs(alpha - p3)>a_error && fabs(ds)>c_error );

 Free(SF);
 
 return s3;
}


double se2fu_crit(double l, double L0, double cu, double hs, double sigma,  int df, int N, int qm)
{ double s1, s2, s3, ds, L1, L2, L3;

 s2 = 2. - cu;
 if ( s2 < 0.1 ) s2 = 0.1;
 L2 = se2_iglarl(l,s2,cu,hs,sigma,df,N,qm);
 if ( L2 < L0 ) {
   do {
     s1 = s2;
     s2 *= 0.8;
     L2 = se2_iglarl(l,s2,cu,hs,sigma,df,N,qm);
   } while ( L2 < L0 );
 } else {
   do {
     s1 = s2;
     s2 *= 1.2;
     L2 = se2_iglarl(l,s2,cu,hs,sigma,df,N,qm);
   } while ( L2 > L0 );
 }

 L1 = se2_iglarl(l,s1,cu,hs,sigma,df,N,qm);

 do {
   s3 = s1 + (L0-L1)/(L2-L1) * (s2-s1);
   L3 = se2_iglarl(l,s3,cu,hs,sigma,df,N,qm);
   ds = s3-s2; s1 = s2; L1 = L2; s2 = s3; L2 = L3;
 } while ( fabs(L0-L3)>1e-7 && fabs(ds)>1e-9 );
 
 return s3;
}


double stde2fu_crit(double l, double L0, double cu, double hs, double sigma,  int df, int N, int qm)
{ double s1, s2, s3, ds, L1, L2, L3, norm;

 norm = sqrt(df);
 s2 = 2. - cu;
 if ( s2 < 0.1 ) s2 = 0.1;
 L2 = stde2_iglarl(l, s2, cu, hs, sigma, df, N, qm);
 
 if ( L2 < L0 ) {
   do {
     s1 = s2;
     L1 = L2;
     s2 -= .2/norm;
     L2 = stde2_iglarl(l, s2, cu, hs, sigma, df, N, qm);
   } while ( L2 < L0 );
 } else {
   do {
     s1 = s2;
     L1 = L2;
     s2 += .2/norm;
     L2 = stde2_iglarl(l, s2, cu, hs, sigma, df, N, qm);
   } while ( L2 > L0 );
 }

 do {
   s3 = s1 + (L0-L1)/(L2-L1) * (s2-s1);
   L3 = stde2_iglarl(l, s3, cu, hs, sigma, df, N,qm);
   ds = s3-s2; s1 = s2; L1 = L2; s2 = s3; L2 = L3;
 } while ( fabs(L0-L3)>1e-7 && fabs(ds)>1e-9 );
 
 return s3;
}


double se2fu_crit_prerun_SIGMA(double l, double L0, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate)
{ double s1, s2, s3, ds, L1=0., L2=0., L3=0.;

 s2 = cu/2.;
 L2 = se2_iglarl_prerun_SIGMA(l, s2, cu, hs, sigma, df1, df2, N, qm1, qm2, truncate);
 if ( L2 < L0 ) {
   do {
     L1 = L2;
     s2 -= .1;
     L2 = se2_iglarl_prerun_SIGMA(l, s2, cu, hs, sigma, df1, df2, N, qm1, qm2, truncate);
   } while ( L2 < L0 && s2 > 0.);
   s1 = s2 + .1;
 } else {
   do {
     L1 = L2;
     s2 += .1;
     L2 = se2_iglarl_prerun_SIGMA(l, s2, cu, hs, sigma, df1, df2, N, qm1, qm2, truncate);
   } while ( L2 > L0 && s2 < hs );
   s1 = s2 - .1;
 }
 
 do {
   s3 = s1 + (L0-L1)/(L2-L1) * (s2-s1);
   L3 = se2_iglarl_prerun_SIGMA(l, s3, cu, hs, sigma, df1, df2, N, qm1, qm2, truncate);
   ds = s3-s2; s1 = s2; L1 = L2; s2 = s3; L2 = L3;
 } while ( fabs(L0-L3)>1e-6 && fabs(ds)>1e-9 );

 return s3;
}


double se2fu_q_crit(double l, int L0, double alpha, double cu, double hs, double sigma, int df, int N, int qm, double c_error, double a_error)
{ double s1, s2, s3, ds, p1, p2, p3, *SF;
  int result=1;

 SF  = vector(L0);  

 /*s2 = cu/2.; */
 s2 = hs/2.;
 result = se2_sf(l, s2, cu, hs, sigma, df, N, L0, qm, SF);
 if ( result != 0 ) warning("trouble in se2fu_q_crit [package spc]");
 p2 = 1. - SF[L0-1];
 
 if ( p2 < alpha ) {
   do {
     p1 = p2;
     s2 *= 1.2;
     result = se2_sf(l, s2, cu, hs, sigma, df, N, L0, qm, SF);
     if ( result != 0 ) warning("trouble in se2fu_q_crit [package spc]");
     p2 = 1. - SF[L0-1];
   } while ( p2 < alpha );
   s1 = s2 - .1;
 } else {
   do {
     p1 = p2;
     s2 /= 1.2;
     result = se2_sf(l, s2, cu, hs, sigma, df, N, L0, qm, SF);
     if ( result != 0 ) warning("trouble in se2fu_q_crit [package spc]");
     p2 = 1. - SF[L0-1];
   } while ( p2 >= alpha );
   s1 = s2 + .1;
 }

 do {
   s3 = s1 + (alpha - p1)/( p2 - p1 ) * (s2-s1);
   result = se2_sf(l, s3, cu, hs, sigma, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in se2fu_q_crit [package spc]");
   p3 = 1. - SF[L0-1];
   ds = s3 - s2; s1 = s2; p1 = p2; s2 = s3; p2 = p3;
 } while ( fabs(alpha - p3)>a_error && fabs(ds)>c_error );

 Free(SF);
 
 return s3;
}


double se2fu_q_crit_prerun_SIGMA(double l, int L0, double alpha, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate, int tail_approx, double c_error, double a_error)
{ double s1, s2, s3, ds, p1, p2, p3, *SF;
  int result=1, schritt=0, maxschritt=30;
  
 SF  = vector(L0);
 
 s2 = se2fu_q_crit(l, L0, alpha, cu, hs, sigma, df1, N, qm1, c_error, a_error);                     
 if ( tail_approx ) result = se2_sf_prerun_SIGMA_deluxe(l, s2, cu, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
 else result = se2_sf_prerun_SIGMA(l, s2, cu, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
 if ( result != 0 ) warning("trouble in se2fu_q_crit_prerun_SIGMA [package spc]");
 p2 = 1. - SF[L0-1];
 
 if ( p2 < alpha ) {
   do {
     p1 = p2;
     s1 = s2;
     s2 *= 1.1;
     if ( tail_approx ) result = se2_sf_prerun_SIGMA_deluxe(l, s2, cu, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     else result = se2_sf_prerun_SIGMA(l, s2, cu, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     if ( result != 0 ) warning("trouble in se2fu_q_crit_prerun_SIGMA [package spc]");
     p2 = 1. - SF[L0-1];
   } while ( p2 < alpha && s2 < hs );
 } else {
   do {
     p1 = p2;
     s1 = s2;
     s2 /= 1.1;
     if ( tail_approx ) result = se2_sf_prerun_SIGMA_deluxe(l, s2, cu, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     else result = se2_sf_prerun_SIGMA(l, s2, cu, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     if ( result != 0 ) warning("trouble in se2fu_q_crit_prerun_SIGMA [package spc]");
     p2 = 1. - SF[L0-1];
   } while ( p2 >= alpha && s2 > 0. );
 }
 
 schritt = 0;
 do {
   schritt++;
   s3 = s1 + (alpha - p1)/( p2 - p1 ) * (s2-s1);
   if ( tail_approx ) result = se2_sf_prerun_SIGMA_deluxe(l, s3, cu, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
   else result = se2_sf_prerun_SIGMA(l, s3, cu, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
   if ( result != 0 ) warning("trouble in se2fu_q_crit_prerun_SIGMA [package spc]");
   p3 = 1. - SF[L0-1];
   ds = s3 - s2; s1 = s2; p1 = p2; s2 = s3; p2 = p3;
 } while ( fabs(alpha - p3)>a_error && fabs(ds)>c_error && schritt<maxschritt );

 if ( schritt >= maxschritt ) warning("secant rule in se2fu_q_crit_prerun_SIGMA did not converge");
 
 Free(SF);
 
 return s3;
}


int se2_crit_prerun_SIGMA(double l, double L0, double *cl, double *cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate)
{ double s1, s2, s3, ds, sl1, sl2, sl3, csl, Lm, Lp;

 csl = hs/2.;
 s1 = se2lu_crit_prerun_SIGMA(l, L0, csl, hs, sigma, df1, df2, N, qm1, qm2, truncate); 
 Lm = se2_iglarl_prerun_SIGMA(l, csl, s1, hs, sigma-lmEPS, df1, df2, N, qm1, qm2, truncate);
 Lp = se2_iglarl_prerun_SIGMA(l, csl, s1, hs, sigma+lmEPS, df1, df2, N, qm1, qm2, truncate);  
 sl1 = (Lp-Lm)/(2.*lmEPS);
 
 s2 = s1 + .05;
 csl = se2fu_crit_prerun_SIGMA(l, L0, s2, hs, sigma, df1, df2, N, qm1, qm2, truncate); 
 Lm = se2_iglarl_prerun_SIGMA(l, csl, s2, hs, sigma-lmEPS, df1, df2, N, qm1, qm2, truncate);
 Lp = se2_iglarl_prerun_SIGMA(l, csl, s2, hs, sigma+lmEPS, df1, df2, N, qm1, qm2, truncate); 
 sl2 = (Lp-Lm)/(2.*lmEPS);

 do {
   s3 = s1 - sl1/(sl2-sl1) * (s2-s1);
   csl = se2fu_crit_prerun_SIGMA(l, L0, s3, hs, sigma, df1, df2, N, qm1, qm2, truncate);
   Lm = se2_iglarl_prerun_SIGMA(l, csl, s3, hs, sigma-lmEPS, df1, df2, N, qm1, qm2, truncate);
   Lp = se2_iglarl_prerun_SIGMA(l, csl, s3, hs, sigma+lmEPS, df1, df2, N, qm1, qm2, truncate);
   sl3 = (Lp-Lm)/(2.*lmEPS);
   ds = s3-s2; s1 = s2; sl1 = sl2; s2 = s3; sl2 = sl3;
 } while ( fabs(sl3)>1e-6 && fabs(ds)>1e-9 );

 *cl = csl; *cu = s3;

 return 0;
}


int se2_crit_unbiased(double l, double L0, double *cl, double *cu, double hs, double sigma, int df, int N, int qm)
{ double s1, s2, s3, ds, sl1, sl2, sl3, csl, Lm, Lp, step;

/* printf("\n\nse2_crit_unbiased\n\n");*/

 step = .1/sqrt(df); 
 s1 = seU_crit(l,L0,hs,sigma,df,N,qm);
 csl = 0.;
 Lm = seU_iglarl(l,s1,hs,sigma-lmEPS,df,N,qm);
 Lp = seU_iglarl(l,s1,hs,sigma+lmEPS,df,N,qm);
 sl1 = (Lp-Lm)/(2.*lmEPS);
/* printf("0 :: cl = %.4f,\tcu = %.12f,\tsl = %.6f\n", csl, s1, sl1);*/
 
 s2 = s1;
 sl2 = sl1;
 do { 
   s1 = s2;
   sl1 = sl2;
   s2 = s1 + step;
   csl = se2fu_crit(l,L0,s2,hs,sigma,df,N,qm);
   Lm = se2_iglarl(l,csl,s2,hs,sigma-lmEPS,df,N,qm);
   Lp = se2_iglarl(l,csl,s2,hs,sigma+lmEPS,df,N,qm);
   sl2 = (Lp-Lm)/(2.*lmEPS);
/*   printf("1 :: cl = %.4f,\tcu = %.12f,\tsl = %.6f\n", csl, s2, sl2);*/
 } while ( sl2 < 0. );

 do {
   s3 = s1 - sl1/(sl2-sl1) * (s2-s1);
   csl = se2fu_crit(l,L0,s3,hs,sigma,df,N,qm);
   Lm = se2_iglarl(l,csl,s3,hs,sigma-lmEPS,df,N,qm);
   Lp = se2_iglarl(l,csl,s3,hs,sigma+lmEPS,df,N,qm);
   sl3 = (Lp-Lm)/(2.*lmEPS);
/*   printf("2 :: cl = %.4f,\tcu = %.12f,\tsl = %.6f\n", csl, s3, sl3);*/
   ds = s3-s2; s1 = s2; sl1 = sl2; s2 = s3; sl2 = sl3;
 } while ( fabs(sl3)>1e-6 && fabs(ds)>1e-12 );

 *cl = csl; *cu = s3;

/* printf("\n\n");*/

 return 0;
}


int stde2_crit_unbiased(double l, double L0, double *cl, double *cu, double hs, double sigma, int df, int N, int qm)
{ double s1, s2, s3, ds, sl1, sl2, sl3, csl, Lm, Lp, step;

 step = .1/sqrt(df); 
 s1 = stdeU_crit(l, L0, hs, sigma, df, N, qm);
 csl = 0.;
 Lm = stdeU_iglarl(l, s1, hs, sigma-lmEPS, df, N, qm);
 Lp = stdeU_iglarl(l, s1, hs, sigma+lmEPS, df, N, qm);
 sl1 = (Lp-Lm)/(2.*lmEPS);
 
 s2 = s1;
 sl2 = sl1;
 do { 
   s1 = s2;
   sl1 = sl2;
   s2 = s1 + step;
   csl = stde2fu_crit(l, L0, s2, hs, sigma, df, N, qm);
   Lm = stde2_iglarl(l, csl, s2, hs, sigma-lmEPS, df, N, qm);
   Lp = stde2_iglarl(l, csl, s2, hs, sigma+lmEPS, df, N, qm);
   sl2 = (Lp-Lm)/(2.*lmEPS);
 } while ( sl2 < 0. );

 do {
   s3 = s1 - sl1/(sl2-sl1) * (s2-s1);
   csl = stde2fu_crit(l, L0, s3, hs, sigma, df, N, qm);
   Lm = stde2_iglarl(l, csl, s3, hs, sigma-lmEPS, df, N, qm);
   Lp = stde2_iglarl(l, csl, s3, hs, sigma+lmEPS, df, N, qm);
   sl3 = (Lp-Lm)/(2.*lmEPS);
   ds = s3-s2; s1 = s2; sl1 = sl2; s2 = s3; sl2 = sl3;
 } while ( fabs(sl3)>1e-7 && fabs(ds)>1e-9 );

 *cl = csl; *cu = s3;

 return 0;
}


int se2_crit_eqtails(double l, double L0, double *cl, double *cu, double hs, double sigma, int df, double ur, int N, int qm)
{ double u1, u2, du, l1, l2, dl, lARL1, lARL2, uARL1, uARL2, ARL22, ARL12, ARL21,
         f11, f22, f21, f12, d11, d22, d21, d12, nenner;
 
 l1 = seLR_crit(l, 2.*L0, ur, hs, sigma, df, N, qm);
 l2 = l1 * 0.9;;
 
 u1 = seU_crit(l, 2.*L0, hs, sigma, df, N, qm); 
 u2 = u1 * 1.1;
 
 /*ARL22 =  se2_iglarl(l, l1, u1, hs, sigma, df, N, qm);*/
 
 lARL2 = seLR_iglarl(l, l2, ur, hs, sigma, df, N, qm);
 uARL2 =  seU_iglarl(l, u2, hs, sigma, df, N, qm);
 ARL22 =  se2_iglarl(l, l2, u2, hs, sigma, df, N, qm);
 
 /*printf("(+)\tl1 = %.6f,\tu1 = %.6f\n", l1, u1);
 printf("(+)\tl2 = %.6f,\tu2 = %.6f,\tllARL2 = %.2f,\tuARL2 = %.2f,\tARL22 = %.2f\n\n", l2, u2, lARL2, uARL2, ARL22);*/
 
  do {
   lARL1 = seLR_iglarl(l, l1, ur, hs, sigma, df, N, qm);
   uARL1 =  seU_iglarl(l, u1, hs, sigma, df, N, qm);
   ARL12 =  se2_iglarl(l, l1, u2, hs, sigma, df, N, qm);
   ARL21 =  se2_iglarl(l, l2, u1, hs, sigma, df, N, qm);
   
   /*printf("(*)\tlARL1 = %.2f,\tuARL1 = %.2f,\tARL12 = %.2f,\tARL21 = %.2f\n", lARL1, uARL1, ARL12, ARL21);*/
  
   /* difference quotient */
   f11 = (ARL22 - ARL12)/(l2-l1); f12 = (ARL22 - ARL21)/(u2-u1);
   f21 = (lARL2 - lARL1)/(l2-l1); f22 = (uARL1 - uARL2)/(u2-u1);

   /* inverse of the difference quotient */
   nenner = f11*f22 - f12*f21;
   d11 =  f22/nenner;  d12 = -f12/nenner;
   d21 = -f21/nenner;  d22 =  f11/nenner;

   dl = d11*(ARL22-L0) + d12*(lARL2-uARL2);
   du = d21*(ARL22-L0) + d22*(lARL2-uARL2);

   l1 = l2;   u1 = u2;
   l2 -= dl;  u2 -= du;
   
   lARL2 = seLR_iglarl(l, l2, ur, hs, sigma, df, N, qm);
   uARL2 =  seU_iglarl(l, u2, hs, sigma, df, N, qm);
   ARL22 =  se2_iglarl(l, l2, u2, hs, sigma, df, N, qm);
   
   /*printf("(*)\tl2 = %.6f,\tu2 = %.6f,\tlARL2 = %.2f,\tuARL2 = %.2f,\tARL22 = %.2f\n\n", l2, u2, lARL2, uARL2, ARL22);*/
   
 } while (  (fabs(L0-ARL22)>1e-6 || fabs(lARL2-uARL2)>1e-6) && (fabs(l2-l1)>1e-9 || fabs(u2-u1)>1e-9)  );

 *cl = l2; *cu = u2;

 return 0;
}


int stde2_crit_eqtails(double l, double L0, double *cl, double *cu, double hs, double sigma, int df, double ur, int N, int qm)
{ double u1, u2, du, l1, l2, dl, lARL1, lARL2, uARL1, uARL2, ARL22, ARL12, ARL21,
         f11, f22, f21, f12, d11, d22, d21, d12, nenner;
 
 l1 = stdeLR_crit(l, 2.*L0, ur, hs, sigma, df, N, qm);
 l2 = l1 - .05;
 u1 =  stdeU_crit(l, 2.*L0, hs, sigma, df, N, qm); 
 u2 = u1 + .05;
 ARL22 =  stde2_iglarl(l, l1, u1, hs, sigma, df, N, qm);
 
 lARL2 = stdeLR_iglarl(l, l2, ur, hs, sigma, df, N, qm);
 uARL2 =  stdeU_iglarl(l, u2, hs, sigma, df, N, qm);
 ARL22 =  stde2_iglarl(l, l2, u2, hs, sigma, df, N, qm);
 
  do {
   lARL1 = stdeLR_iglarl(l, l1, ur, hs, sigma, df, N, qm);
   uARL1 =  stdeU_iglarl(l, u1, hs, sigma, df, N, qm);
   ARL12 =  stde2_iglarl(l, l1, u2, hs, sigma, df, N, qm);
   ARL21 =  stde2_iglarl(l, l2, u1, hs, sigma, df, N, qm);
  
   /* difference quotient */
   f11 = (ARL22 - ARL12)/(l2-l1); f12 = (ARL22 - ARL21)/(u2-u1);
   f21 = (lARL2 - lARL1)/(l2-l1); f22 = (uARL1 - uARL2)/(u2-u1);

   /* inverse of the difference quotient */
   nenner = f11*f22 - f12*f21;
   d11 =  f22/nenner;  d12 = -f12/nenner;
   d21 = -f21/nenner;  d22 =  f11/nenner;

   dl = d11*(ARL22-L0) + d12*(lARL2-uARL2);
   du = d21*(ARL22-L0) + d22*(lARL2-uARL2);

   l1 = l2;   u1 = u2;
   l2 -= dl;  u2 -= du;
   
   lARL2 = stdeLR_iglarl(l, l2, ur, hs, sigma, df, N, qm);
   uARL2 =  stdeU_iglarl(l, u2, hs, sigma, df, N, qm);
   ARL22 =  stde2_iglarl(l, l2, u2, hs, sigma, df, N, qm);
   
 } while (  (fabs(L0-ARL22)>1e-6 || fabs(lARL2-uARL2)>1e-6) && (fabs(l2-l1)>1e-9 || fabs(u2-u1)>1e-9)  );

 *cl = l2; *cu = u2;

 return 0;
}


double se2_crit_sym(double l, double L0, double hs, double sigma, int df, int N, int qm)
{ double cu1, cu2, cu3, cl1, cl2, cl3, L1, L2, L3, du, step;

  cu2 = seU_crit(l, L0, hs, sigma, df, N, qm);
  if ( cu2 < 2. ) {
    step = (2.-cu2)/10.;
    cu2 += step;
    cl2 = 2. - cu2;
    L2 = se2_iglarl(l, cl2, cu2, hs, sigma, df, N, qm);
    cu1 = cu2 + step;
    cl1 = 2. - cu1;
    L1 = se2_iglarl(l, cl1, cu1, hs, sigma, df, N, qm);

    do {
      cu3 = cu1 + (L0-L1)/(L2-L1) * (cu2-cu1);
      cl3 = 2. - cu3;
      L3 = se2_iglarl(l, cl3, cu3, hs, sigma, df, N, qm);
      du = cu3-cu2; cu1 = cu2; L1 = L2; cu2 = cu3; L2 = L3;
      if ( L3 < 1. ) error("invalid ARL value");
    } while ( (fabs(L0-L3)>1e-6) && (fabs(du)>1e-9) ); 
  } else {
    error("symmetric design not possible");
    cu3 = -1.;
  }
  return cu3;
}


double stde2_crit_sym(double l, double L0, double hs, double sigma, int df, int N, int qm)
{ double cu1, cu2, cu3, cl1, cl2, cl3, L1, L2, L3, du, step, mitte;

  mitte = c_four((double)df);
  cu2 = stdeU_crit(l, L0, hs, sigma, df, N, qm);
  if ( cu2 < 2. ) {
    step = (2.-cu2)/10.;
    cu2 += step;
    cl2 = 2.*mitte - cu2;
    L2 = stde2_iglarl(l, cl2, cu2, hs, sigma, df, N, qm);
    
    cu1 = cu2 + step;
    cl1 = 2.*mitte - cu1;
    L1 = stde2_iglarl(l, cl1, cu1, hs, sigma, df, N, qm);

    do {
      cu3 = cu1 + (L0-L1)/(L2-L1) * (cu2-cu1);
      cl3 = 2.*mitte - cu3;
      L3 = stde2_iglarl(l, cl3, cu3, hs, sigma, df, N, qm);
      du = cu3-cu2; cu1 = cu2; L1 = L2; cu2 = cu3; L2 = L3;
      if ( L3 < 1. ) error("invalid ARL value");
    } while ( (fabs(L0-L3)>1e-7) && (fabs(du)>1e-9) ); 
  } else {
    error("symmetric design not possible");
    cu3 = -1.;
  }
  return cu3;
}


int se2_q_crit(double l, int L0, double alpha, double *cl, double *cu, double hs, double sigma, int df, int N, int qm, double c_error, double a_error)
{ double s1, s2, s3, ds, sl1, sl2, sl3, csl, Pm, Pp, *SF;
  int result=1;

 SF  = vector(L0);

 s1 = seU_q_crit(l, L0, alpha, hs, sigma, df, N, qm, c_error, a_error);
 csl = 0.;
 result = seU_sf(l, s1, hs, sigma-lmEPS, df, N, L0, qm, SF);
 if ( result != 0 ) warning("trouble in se2_q_crit [package spc]");
 Pm = 1. - SF[L0-1];
 result = seU_sf(l, s1, hs, sigma+lmEPS, df, N, L0, qm, SF);
 if ( result != 0 ) warning("trouble in se2_q_crit [package spc]");
 Pp = 1. - SF[L0-1]; 
 sl1 = ( Pp - Pm )/(2.*lmEPS);

 s2 = s1 + .05;
 csl = se2fu_q_crit(l, L0, alpha, s2, hs, sigma, df, N, qm, c_error, a_error);
 result = se2_sf(l, csl, s2, hs, sigma-lmEPS, df, N, L0, qm, SF);
 if ( result != 0 ) warning("trouble in se2_q_crit [package spc]");
 Pm = 1. - SF[L0-1];
 result = se2_sf(l, csl, s2, hs, sigma+lmEPS, df, N, L0, qm, SF);
 if ( result != 0 ) warning("trouble in se2_q_crit [package spc]");
 Pp = 1. - SF[L0-1]; 
 sl2 = ( Pp - Pm )/(2.*lmEPS);

 do {
   s3 = s1 - sl1/(sl2-sl1) * (s2-s1);
   csl = se2fu_q_crit(l, L0, alpha, s3, hs, sigma, df, N, qm, c_error, a_error);
   result = se2_sf(l, csl, s3, hs, sigma-lmEPS, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in se2_q_crit [package spc]");
   Pm = 1. - SF[L0-1];
   result = se2_sf(l, csl, s3, hs, sigma+lmEPS, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in se2_q_crit [package spc]");
   Pp = 1. - SF[L0-1]; 
   sl3 = ( Pp - Pm )/(2.*lmEPS);
   ds = s3-s2; s1 = s2; sl1 = sl2; s2 = s3; sl2 = sl3;
 } while ( fabs(sl3)>a_error && fabs(ds)>c_error );

 *cl = csl; *cu = s3;
 
 Free(SF);

 return 0;
}


int se2_q_crit_class(double l, int L0, double alpha, double *cl, double *cu, double hs, double sigma, int df, double ur, int N, int qm, double c_error, double a_error)
{ double u1, u2, du, l1, l2, dl, lA1, lA2, uA1, uA2, A22, A12, A21,
         f11, f22, f21, f12, d11, d22, d21, d12, nenner, *SF;
  int result=1;
  
 SF  = vector(L0);
   
 l1 = seLR_q_crit(l, L0, alpha/2., ur, hs, sigma, df, N, qm, c_error, a_error);
 l2 = l1 - .05;
 u1 = seU_q_crit(l, L0, alpha/2., hs, sigma, df, N, qm, c_error, a_error);
 u2 = u1 + .05;
 
 result = seLR_sf(l, l2, ur, hs, sigma, df, N, L0, qm, SF);
 if ( result != 0 ) warning("trouble in se2_q_crit_class [package spc]");
 lA2 = 1. - SF[L0-1];
 result =  seU_sf(l,     u2, hs, sigma, df, N, L0, qm, SF);
 if ( result != 0 ) warning("trouble in se2_q_crit_class [package spc]");
 uA2 = 1. - SF[L0-1];
 result =  se2_sf(l, l2, u2, hs, sigma, df, N, L0, qm, SF);
 if ( result != 0 ) warning("trouble in se2_q_crit_class [package spc]");
 A22 = 1. - SF[L0-1];
 
  do {
   result = seLR_sf(l, l1, ur, hs, sigma, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in se2_q_crit_class [package spc]");
   lA1 = 1. - SF[L0-1];
   result =  seU_sf(l,     u1, hs, sigma, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in se2_q_crit_class [package spc]");
   uA1 = 1. - SF[L0-1];
   result =  se2_sf(l, l1, u2, hs, sigma, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in se2_q_crit_class [package spc]");
   A12 = 1. - SF[L0-1];
   result =  se2_sf(l, l2, u1, hs, sigma, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in se2_q_crit_class [package spc]");
   A21 = 1. - SF[L0-1];
  
   /* difference quotient */
   f11 = (A22 - A12)/(l2-l1); f12 = (A22 - A21)/(u2-u1);
   f21 = (lA2 - lA1)/(l2-l1); f22 = (uA1 - uA2)/(u2-u1);

   /* inverse of the difference quotient */
   nenner = f11*f22 - f12*f21;
   d11 =  f22/nenner;  d12 = -f12/nenner;
   d21 = -f21/nenner;  d22 =  f11/nenner;

   dl = d11*(A22-alpha) + d12*(lA2-uA2);
   du = d21*(A22-alpha) + d22*(lA2-uA2);

   l1 = l2;   u1 = u2;
   l2 -= dl;  u2 -= du;
   
   result = seLR_sf(l, l2, ur, hs, sigma, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in se2_q_crit_class [package spc]");
   lA2 = 1. - SF[L0-1];
   result =  seU_sf(l,     u2, hs, sigma, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in se2_q_crit_class [package spc]");
   uA2 = 1. - SF[L0-1];
   result =  se2_sf(l, l2, u2, hs, sigma, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in se2_q_crit_class [package spc]");
   A22 = 1. - SF[L0-1];
   
 } while (  (fabs(alpha-A22)>1e-9 || fabs(lA2-uA2)>1e-9) && (fabs(l2-l1)>1e-9 || fabs(u2-u1)>1e-9)  );

 *cl = l2; *cu = u2;
 
 Free(SF);

 return 0;
}


int se2_q_crit_prerun_SIGMA(double l, int L0, double alpha, double *cl, double *cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate, int tail_approx, double c_error, double a_error)
{ double s1, s2, s3, ds, sl1, sl2, sl3, csl, Pm, Pp, *SF;
  int result=1;

 /*printf("\n\nEs geht los!\n\n");*/
  
 SF  = vector(L0);

 s1 = seU_q_crit_prerun_SIGMA(l, L0, alpha, hs, sigma, df1, df2, N, qm1, qm2, truncate, tail_approx, c_error, a_error);
 csl = 0.;
 
 /*printf("Startwert berechnet, s1 = %.6f\n", s1);*/
 
 if ( tail_approx ) result = seU_sf_prerun_SIGMA_deluxe(l, s1, hs, sigma-lmEPS, df1, df2, L0, qm1, qm2, truncate, SF);
 else result = seU_sf_prerun_SIGMA(l, s1, hs, sigma-lmEPS, df1, df2, L0, qm1, qm2, truncate, SF);
 if ( result != 0 ) warning("trouble in se2_q_crit_prerun_SIGMA [package spc]");
 Pm = 1. - SF[L0-1];
 if ( tail_approx ) result = seU_sf_prerun_SIGMA_deluxe(l, s1, hs, sigma+lmEPS, df1, df2, L0, qm1, qm2, truncate, SF);
 else result = seU_sf_prerun_SIGMA(l, s1, hs, sigma+lmEPS, df1, df2, L0, qm1, qm2, truncate, SF);
 if ( result != 0 ) warning("trouble in se2_q_crit_prerun_SIGMA [package spc]");
 Pp = 1. - SF[L0-1]; 
 sl1 = ( Pp - Pm )/(2.*lmEPS);
 
 /*printf("slope = %.6f\n\n", sl1);*/

 if ( sl1 > 0 ) {
   do {
     s2 = s1;
     sl2 = sl1;
     s1 *= 1.1;
     csl = se2fu_q_crit_prerun_SIGMA(l, L0, alpha, s1, hs, sigma, df1, df2, N, qm1, qm2, truncate, tail_approx, c_error, a_error);
     if ( tail_approx ) result = se2_sf_prerun_SIGMA_deluxe(l, csl, s1, hs, sigma-lmEPS, df1, df2, L0, qm1, qm2, truncate, SF);
     else result = se2_sf_prerun_SIGMA(l, csl, s1, hs, sigma-lmEPS, df1, df2, L0, qm1, qm2, truncate, SF);
     if ( result != 0 ) warning("trouble in se2_q_crit_prerun_SIGMA [package spc]");
     Pm = 1. - SF[L0-1];
     if ( tail_approx ) result = se2_sf_prerun_SIGMA_deluxe(l, csl, s1, hs, sigma+lmEPS, df1, df2, L0, qm1, qm2, truncate, SF);
     else result = se2_sf_prerun_SIGMA(l, csl, s1, hs, sigma+lmEPS, df1, df2, L0, qm1, qm2, truncate, SF);
     if ( result != 0 ) warning("trouble in se2_q_crit_prerun_SIGMA [package spc]");
     Pp = 1. - SF[L0-1]; 
     sl1 = ( Pp - Pm )/(2.*lmEPS);
    /*printf("(i) s1 = %.6f, slope = %.6f\n", s1, sl1);*/
   } while ( sl1 > 0 );
 } else {
   do {
     s2 = s1;
     sl2 = sl1;
     s1 /= 1.1;
     csl = se2fu_q_crit_prerun_SIGMA(l, L0, alpha, s1, hs, sigma, df1, df2, N, qm1, qm2, truncate, tail_approx, c_error, a_error);
     if ( tail_approx ) result = se2_sf_prerun_SIGMA_deluxe(l, csl, s1, hs, sigma-lmEPS, df1, df2, L0, qm1, qm2, truncate, SF);
     else result = se2_sf_prerun_SIGMA(l, csl, s1, hs, sigma-lmEPS, df1, df2, L0, qm1, qm2, truncate, SF);
     if ( result != 0 ) warning("trouble in se2_q_crit_prerun_SIGMA [package spc]");
     Pm = 1. - SF[L0-1];
     if ( tail_approx ) result = se2_sf_prerun_SIGMA_deluxe(l, csl, s1, hs, sigma+lmEPS, df1, df2, L0, qm1, qm2, truncate, SF);
     else result = se2_sf_prerun_SIGMA(l, csl, s1, hs, sigma+lmEPS, df1, df2, L0, qm1, qm2, truncate, SF);
     if ( result != 0 ) warning("trouble in se2_q_crit_prerun_SIGMA [package spc]");
     Pp = 1. - SF[L0-1]; 
     sl1 = ( Pp - Pm )/(2.*lmEPS);
     /*printf("(ii) s1 = %.6f, slope = %.6f\n", s1, sl1);*/
   } while ( sl1 < 0 );
 }

 do {
   s3 = s1 - sl1/(sl2-sl1) * (s2-s1);
   csl = se2fu_q_crit_prerun_SIGMA(l, L0, alpha, s3, hs, sigma, df1, df2, N, qm1, qm2, truncate, tail_approx, c_error, a_error);
   if ( tail_approx ) result = se2_sf_prerun_SIGMA_deluxe(l, csl, s3, hs, sigma-lmEPS, df1, df2, L0, qm1, qm2, truncate, SF);
   else result = se2_sf_prerun_SIGMA(l, csl, s3, hs, sigma-lmEPS, df1, df2, L0, qm1, qm2, truncate, SF);
   if ( result != 0 ) warning("trouble in se2_q_crit_prerun_SIGMA [package spc]");
   Pm = 1. - SF[L0-1];
   if ( tail_approx ) result = se2_sf_prerun_SIGMA_deluxe(l, csl, s3, hs, sigma+lmEPS, df1, df2, L0, qm1, qm2, truncate, SF);
   else result = se2_sf_prerun_SIGMA(l, csl, s3, hs, sigma+lmEPS, df1, df2, L0, qm1, qm2, truncate, SF);
   if ( result != 0 ) warning("trouble in se2_q_crit_prerun_SIGMA [package spc]");
   Pp = 1. - SF[L0-1]; 
   sl3 = ( Pp - Pm )/(2.*lmEPS);
   /*printf("(iii) s3 = %.6f, slope = %.6f\n", s3, sl3);*/
   ds = s3-s2; s1 = s2; sl1 = sl2; s2 = s3; sl2 = sl3;
 } while ( fabs(sl3)>a_error && fabs(ds)>c_error );

 *cl = csl; *cu = s3;
 
 Free(SF);

 return 0;
}


double seUR_iglarl(double l, double cl, double cu, double hs, double sigma, int df, int N, int qm)
{ double *a, *g, *w, *z, *t, h, arl, Hij, xl, za, dN, ddf, s2,
         t0, t1, x0, x1, dummy;
  int i, j, k, qi, qj, M, Ntilde, NN, ii, it, jj;

 M = ceil( (log(cl)-log(cu))/log(1.-l) );
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde - 1.;

 a = matrix(NN,NN);
 g = vector(NN);
 t = vector(NN);
 w = vector(qm);
 z = vector(qm);

 for(i=0;i<M;i++) {
   t0 = cl/pow(1.-l,(double)(i));
   t1 = t0/(1.-l);
   if (t1>cu) t1 = cu;

   for (j=1;j<Ntilde;j++) { /* node_i,Ntilde-1 = node_i+1,0 */
     h = cos( PI/dN *(dN-j) );
     t[i*(Ntilde-1)+j] = t0 + (h+1.)/2.*(t1-t0);
     /* Chebyshev Gauss-Lobatto nodes on [t0,t1] */
   }
 }
 t[0] = cl;

 for (i=0;i<M;i++) {
   for (j=1;j<=Ntilde;j++) {
     ii = i*Ntilde + j-1;
     it = i*(Ntilde-1) + j-1;

     za = (1.-l)*t[it];
     if (za<cl) xl = cl; else xl = za;

     for (qi=0;qi<i-1;qi++)
       for (qj=1;qj<=Ntilde;qj++) {
         jj = qi*Ntilde + qj-1;
         a[ii*NN+jj] = 0.;
       }

     if (i>0) {
       qi = i-1;
       t0 = cl/pow(1.-l,(double)qi);
       t1 = t0/(1.-l);
       if (t1>cu) t1 = cu;
       if (t0<xl) x0 = xl; else x0 = t0;
       if (df==2)
         x1 = t1;
       else {
         if (x0-za>1e-10) x0 = sqrt(x0-za); else x0 = 0.;
         if (t1-za>1e-10) x1 = sqrt(t1-za); else x1 = 0.;
       }

       for (qj=1;qj<=Ntilde;qj++) {
         jj = qi*Ntilde + qj-1;

         if (j==1) a[ii*NN+jj] = - Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1);
         else {
           if (fabs(t1-x0)>1e-8) {
             gausslegendre(qm,x0,x1,z,w);
             Hij = 0.;
             for (k=0;k<qm;k++) {
               if (df==2)
                 Hij += w[k] * Tn( (2.*z[k]-t0-t1)/(t1-t0), qj-1) *
                        exp((za-z[k])/s2/l);
               if (df!=2)
                 Hij += w[k] * Tn( (2.*(z[k]*z[k]+za)-t0-t1)/(t1-t0) ,qj-1) *
                        2. * pow(z[k], ddf-1.) * exp(-ddf*z[k]*z[k]/2./s2/l);
             }
             if (df==2) Hij /= s2*l;
             else       Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf,ddf/2.);
             a[ii*NN+jj] = -Hij;
           }
           else a[ii*NN+jj] = 0.;
         }
       }
     }

     for (qi=i;qi<M;qi++) {
       t0 = cl/pow(1.-l,(double)qi);
       t1 = t0/(1.-l);
       if (t1>cu) t1 = cu;
       if (t0<xl) x0 = xl; else x0 = t0;
       if (df==2)
         x1 = t1;
       else {
        if (x0-za>1e-10) x0 = sqrt(x0-za); else x0 = 0.;
        if (t1-za>1e-10) x1 = sqrt(t1-za); else x1 = 0.;
       }

       if (i>0 && j==1 && qi==i) {
         for (qj=1;qj<=Ntilde;qj++) {
           jj = qi*Ntilde + qj-1;
           a[ii*NN+jj] = Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1);
         }
       }

       if (i>0 && j==1 && qi>i) {
         for (qj=1;qj<=Ntilde;qj++) {
           jj = qi*Ntilde + qj-1;
           a[ii*NN+jj] = 0.;
         }
       }

       if (i==0 || j>1) {
         for (qj=1;qj<=Ntilde;qj++) {
           jj = qi*Ntilde + qj-1;
           gausslegendre(qm,x0,x1,z,w);
           Hij = 0.;
           for (k=0;k<qm;k++) {
             if (df==2)
               Hij += w[k] * Tn( (2.*z[k]-t0-t1)/(t1-t0), qj-1) *
                      exp((za-z[k])/s2/l);
             if (df!=2)
               Hij += w[k] * Tn( (2.*(z[k]*z[k]+za)-t0-t1)/(t1-t0),qj-1) *
                      2. * pow(z[k], ddf-1.) * exp(-ddf*z[k]*z[k]/2./s2/l);
           }
           if (df==2) Hij /= s2*l;
           else       Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf,ddf/2.);
           if (qi==i) a[ii*NN+jj] = Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1) -
                                        Hij;
           else a[ii*NN+jj] = -Hij;
         }
       }
     }
     if (i==0) {
       t0 = cl;
       t1 = t0/(1.-l);
       if (t1>cu) t1 = cu;

       for (qj=1;qj<=Ntilde;qj++) {
         dummy = (cl-za)/l/s2;
         if (dummy>0.) {
           if (df==1) dummy = 2.*PHI( sqrt(dummy), 0. ) - 1.;
           if (df==2) dummy = 1. - exp( -dummy );
           if (df>2)  dummy = CHI( df*dummy, df);
         }
         else dummy = 0.;

         a[ii*NN+qj-1] -= dummy * Tn((2.*cl-t0-t1)/(t1-t0),qj-1);
       }
     }
   }
 }

 for (j=0;j<NN;j++) g[j] = 1.;
 for (j=1;j<M;j++) g[Ntilde*j] = 0.;

 LU_solve(a,g,NN);

 arl = 0.;
 for (i=0;i<M;i++) {
   t0 = cl/pow(1.-l,(double)i);
   t1 = t0/(1.-l);
   if (t1>cu) t1 = cu;

   if (t0<=hs && hs<t1)
     for (j=1;j<=Ntilde;j++) {
        ii = i*Ntilde + j-1;
        arl += g[ii] * Tn((2.*hs-t0-t1)/(t1-t0),j-1);
     }
 }

 Free(z);
 Free(w);
 Free(t);
 Free(g);
 Free(a);

 return arl;
}


double stdeUR_iglarl(double l, double cl, double cu, double hs, double sigma, int df, int N, int qm)
{ double *a, *g, *w, *z, *t, h, arl, Hij, xl, za, dN, ddf, s2, t0, t1, x0=0., x1, dummy, v;
  int i, j, k, qi, qj, M, Ntilde, NN, ii, it, jj;

 M = ceil( (log(cl)-log(cu))/log(1.-l) );
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde - 1.;

 a = matrix(NN,NN);
 g = vector(NN);
 t = vector(NN);
 w = vector(qm);
 z = vector(qm);

 for(i=0; i<M; i++) {
   t0 = cl/pow(1.-l,(double)(i));
   t1 = t0/(1.-l);
   if ( t1>cu ) t1 = cu;

   for (j=1; j<Ntilde; j++) { /* node_i,Ntilde-1 = node_i+1,0 */
     h = cos( PI/dN *(dN-j) );
     t[i*(Ntilde-1)+j] = t0 + (h+1.)/2.*(t1-t0); /* Chebyshev Gauss-Lobatto nodes on [t0,t1] */
   }
 }
 t[0] = cl;

 for (i=0; i<M; i++) {
   for (j=1; j<=Ntilde; j++) {
     ii = i*Ntilde + j-1;
     it = i*(Ntilde-1) + j-1;

     za = (1.-l)*t[it];
     if ( za<cl ) xl = cl; else xl = za;

     for (qi=0; qi<i-1; qi++)
       for (qj=1; qj<=Ntilde; qj++) {
         jj = qi*Ntilde + qj-1;
         a[ii*NN+jj] = 0.;
       }

     if ( i>0 ) {
       qi = i-1;
       t0 = cl/pow(1.-l,(double)qi);
       t1 = t0/(1.-l);
       if ( t1>cu ) t1 = cu;
       if ( t0<xl ) x0 = xl; else x0 = t0;
       x1 = t1;

       for (qj=1; qj<=Ntilde; qj++) {
         jj = qi*Ntilde + qj-1;

         if ( j==1 ) a[ii*NN+jj] = - Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1);
         else {
           if (fabs(t1-x0)>1e-8) {
             gausslegendre(qm, x0, x1, z, w);
             Hij = 0.;
             for (k=0; k<qm; k++) {
	       v = (z[k] - za) / l;
               Hij += w[k] * Tn( (2.*z[k]-t0-t1)/(t1-t0), qj-1) * pow(v,ddf-1.)*exp(-ddf/2./s2*v*v);
             }
             Hij *= 2./l/gammafn(ddf/2.)/pow(2.*s2/ddf,ddf/2.);
             a[ii*NN+jj] = -Hij;
           }
           else a[ii*NN+jj] = 0.;
         }
       }
     }

     for (qi=i; qi<M; qi++) {
       t0 = cl/pow(1.-l,(double)qi);
       t1 = t0/(1.-l);
       if ( t1>cu ) t1 = cu;
       if ( t0<xl ) x0 = xl; else x0 = t0; /* Hong Kong & Inez */
       x1 = t1;

       if ( i>0 && j==1 && qi==i ) {
         for (qj=1; qj<=Ntilde; qj++) {
           jj = qi*Ntilde + qj-1;
           a[ii*NN+jj] = Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1);
         }
       }

       if ( i>0 && j==1 && qi>i ) {
         for (qj=1;qj<=Ntilde;qj++) {
           jj = qi*Ntilde + qj-1;
           a[ii*NN+jj] = 0.;
         }
       }

       if ( i==0 || j>1 ) {
         for (qj=1; qj<=Ntilde; qj++) {
           jj = qi*Ntilde + qj-1;
           gausslegendre(qm, x0, x1, z, w);
           Hij = 0.;
           for (k=0;k<qm;k++) {
	     v = (z[k] - za) / l;
             Hij += w[k] * Tn( (2.*z[k]-t0-t1)/(t1-t0), qj-1) * pow(v,ddf-1.)*exp(-ddf/2./s2*v*v);
           }
           Hij *= 2./l/gammafn(ddf/2.)/pow(2.*s2/ddf,ddf/2.);
           if ( qi==i ) a[ii*NN+jj] = Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1) - Hij;
           else a[ii*NN+jj] = -Hij;
         }
       }
     }
     
     if ( i==0 ) {
       t0 = cl;
       t1 = t0/(1.-l);
       if ( t1>cu ) t1 = cu;

       for (qj=1; qj<=Ntilde; qj++) {
	 dummy = 0.;
         v = (cl-za)/l;
         if ( v>0. ) dummy = CHI(ddf/s2*v*v, df);
         a[ii*NN+qj-1] -= dummy * Tn((2.*cl-t0-t1)/(t1-t0),qj-1);
       }
     }
   }
 }

 for ( j=0; j<NN; j++) g[j] = 1.;
 for ( j=1; j<M; j++) g[Ntilde*j] = 0.;

 LU_solve(a, g, NN);

 arl = 0.;
 for (i=0; i<M; i++) {
   t0 = cl/pow(1.-l,(double)i);
   t1 = t0/(1.-l);
   if ( t1>cu ) t1 = cu;

   if ( t0<=hs && hs<t1 )
     for (j=1; j<=Ntilde; j++) {
        ii = i*Ntilde + j-1;
        arl += g[ii] * Tn((2.*hs-t0-t1)/(t1-t0),j-1);
     }
 }

 Free(z);
 Free(w);
 Free(t);
 Free(g);
 Free(a);

 return arl;
}


double seUR_sf(double l, double cl, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0)
{ double *S1s, *S2s, *Pns, *ws, *zs, *zch, *rside, *b, za=0., s2, ddf, xl, xu, dN, Hij, *S00, *p00, *VF0;
  int i, j, k, n, *ps, Ntilde, ihs, M, NN, ii, jj;
 
 M = ceil( (log(cl)-log(cu))/log(1.-l) );
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde;
 
 ihs = floor( (log(cl) - log(hs))/log(1.-l) );
 if ( ihs<0 ) ihs = 0;
 
 S1s = matrix(NN,NN);
 S2s = matrix(NN,NN);
 ps = ivector(NN);
 zch = matrix(M,Ntilde);
 rside = vector(NN+1);
 b   = vector(M+1);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,NN);
 
 S00 = vector(NN);
 p00 = vector(nmax);
 VF0 = vector(NN+1);

/* interval borders b_i = cl/(1-l)^i */
 for (i=0; i<M; i++) b[i] = cl/pow(1.-l, (double)(i));
 b[M] = cu;

 /* Chebyshev nodes on [b_0,b_1],[b_1,b_2],...,[b_M-1,cu] */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     zch[ i*Ntilde+j ] = b[i] + (b[i+1]-b[i])/2.*(1.+cos(PI*(2.*j+1.)/2./dN));
   }

 /* P(L>1)(zch[i,j]) */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     rside[ i*Ntilde+j ] =  CHI( ddf/s2*(cu-(1.-l)*zch[ i*Ntilde+j ])/l, df);
   }
 rside[NN] = CHI( ddf/s2*(cu-(1.-l)*cl)/l, df); /* reflexion at cl */ 

 /* P(zch[i,j] -> zreflect) */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++)
     VF0[ i*Ntilde+j ] = CHI( ddf/s2*(cl-(1.-l)*zch[ i*Ntilde+j ])/l, df) ;
 VF0[NN] = CHI( ddf/s2*cl, df);

 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     za = (1.-l)*zch[ i*Ntilde+j ];
     for (ii=0; ii<M; ii++)
       for (jj=0; jj<Ntilde; jj++) {
         if ( b[ii+1]<za ) S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = 0.;
         else {
           if ( za<b[ii] ) xl = b[ii]; else xl = za;
           xu = b[ii+1];
           if ( df!=2 ) {
             xl = sqrt(xl-za);
             xu = sqrt(xu-za);
           }
           gausslegendre(qm, xl, xu, zs, ws);
           Hij = 0.;
           for (k=0; k<qm; k++)
             if ( df==2 )
               Hij += ws[k]*Tn( (2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj) * exp(-zs[k]/s2/l);
             else
               Hij += ws[k]*Tn( (2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                      * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
           if ( df==2 ) Hij *= exp(za/s2/l)/s2/l;
           else         Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf, ddf/2.);
           S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = Hij;
         }
       }
   }

 za = (1.-l)*cl;
 for (ii=0; ii<M; ii++)
   for (jj=0; jj<Ntilde; jj++) {
     if ( b[ii+1]<za ) S00[ ii*Ntilde+jj ] = 0.;
     else {
       if ( za<b[ii] ) xl = b[ii]; else xl = za;
       xu = b[ii+1];
       if ( df!=2 ) {
         xl = sqrt(xl-za);
         xu = sqrt(xu-za);
       }
       gausslegendre(qm, xl, xu, zs, ws);
       Hij = 0.;
       for (k=0; k<qm; k++)
         if ( df==2 )
           Hij += ws[k]*Tn((2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj) * exp(-zs[k]/s2/l);
         else
           Hij += ws[k] * Tn((2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                  * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
       if ( df==2 ) Hij *= exp(za/s2/l)/s2/l;
       else         Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf, ddf/2.);
       S00[ ii*Ntilde+jj ] = Hij;
     }
   }
 
 for (i=0; i<NN; i++)
   for (j=0; j<NN; j++) S2s[i*NN+j] = 0.;

 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++)
     for (jj=0; jj<Ntilde; jj++)
       S2s[ (i*Ntilde+j)*NN + i*Ntilde+jj ] = Tn( (2.*zch[ i*Ntilde+j ]-b[i+1]-b[i])/(b[i+1]-b[i]), jj); 

 LU_decompose(S2s,ps,NN);

 for (n=1; n<=nmax; n++) {
   if ( n==1) {
     for (i=0; i<M; i++)
       for (j=0; j<Ntilde; j++) {
         Pns[ i*Ntilde+j ] = 0.;
         for (jj=0; jj<Ntilde; jj++)
           Pns[ i*Ntilde+j ] += /*  usual Chebyshev approximation  */
             2./Ntilde * Tn( (2.*zch[i*Ntilde+jj]-b[i+1]-b[i])/(b[i+1]-b[i]), j) * rside[ i*Ntilde+jj ];
         if ( j==0 ) Pns[ i*Ntilde+j ] /= 2.;
       }
     p00[0] = rside[NN];
   }
   else {
     for (i=0; i<NN; i++) {
       rside[i] = VF0[i] * p00[n-2];
       for (j=0; j<NN; j++) rside[i] += S1s[ i*NN+j ] * Pns[ (n-2)*NN+j ];
     }
     LU_solve2(S2s, rside, ps, NN);     
     for (i=0; i<NN; i++) Pns[ (n-1)*NN+i ] = rside[i];    
     p00[n-1] = VF0[NN] * p00[n-2];     
     for (i=0 ;i<NN; i++) p00[n-1] += S00[i] * Pns[ (n-2)*NN+i ];
   }

   p0[n-1] = 0.;
   if ( n==1 )
     p0[0] =  CHI( ddf/s2*(cu-(1.-l)*hs)/l, df);
   else
     for (j=0; j<Ntilde; j++)
       p0[n-1] += Pns[ (n-1)*NN + ihs*Ntilde+j ] * Tn( (2.*hs-b[ihs+1]-b[ihs])/(b[ihs+1]-b[ihs]), j);
 }
 
 Free(Pns);
 Free(zs);
 Free(ws);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);
 
 Free(S00);
 Free(p00);
 Free(VF0);

 return 0;
}


double seUR_sf_deluxe(double l, double cl, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0, int *nstop, double *rho)
{ double *S1s, *S2s, *Pns, *ws, *zs, *zch, *rside, *b, za=0., s2, ddf, xl, xu, dN, Hij, *S00, *p00, *VF0, mn_minus=1., mn_plus=0., oben, unten, q;
  int i, j, k, n, *ps, Ntilde, ihs, M, NN, ii, jj;

 M = ceil( (log(cl)-log(cu))/log(1.-l) );
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde;
 
 ihs = floor( (log(cl) - log(hs))/log(1.-l) );
 if ( ihs<0 ) ihs = 0;

 S1s = matrix(NN,NN);
 S2s = matrix(NN,NN);
 ps = ivector(NN);
 zch = matrix(M,Ntilde);
 rside = vector(NN+1);
 b   = vector(M+1);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,NN);
 
 S00 = vector(NN);
 p00 = vector(nmax);
 VF0 = vector(NN+1);
 
 *nstop = 0;

/* interval borders b_i = cl/(1-l)^i */
 for (i=0; i<M; i++) b[i] = cl/pow(1.-l, (double)(i));
 b[M] = cu;

 /* Chebyshev nodes on [b_0,b_1],[b_1,b_2],...,[b_M-1,cu] */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     zch[ i*Ntilde+j ] = b[i] + (b[i+1]-b[i])/2.*(1.+cos(PI*(2.*j+1.)/2./dN));
   }

 /* P(L>1)(zch[i,j]) */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     rside[ i*Ntilde+j ] =  CHI( ddf/s2*(cu-(1.-l)*zch[ i*Ntilde+j ])/l, df);
   }
 rside[NN] = CHI( ddf/s2*(cu-(1.-l)*cl)/l, df); /* reflexion at cl */
 
 /* P(zch[i,j] -> zreflect) */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++)
     VF0[ i*Ntilde+j ] = CHI( ddf/s2*(cl-(1.-l)*zch[ i*Ntilde+j ])/l, df) ;
 VF0[NN] = CHI( ddf/s2*cl, df);
 
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     za = (1.-l)*zch[ i*Ntilde+j ];
     for (ii=0; ii<M; ii++)
       for (jj=0; jj<Ntilde; jj++) {
         if ( b[ii+1]<za ) S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = 0.;
         else {
           if ( za<b[ii] ) xl = b[ii]; else xl = za;
           xu = b[ii+1];
           if ( df!=2 ) {
             xl = sqrt(xl-za);
             xu = sqrt(xu-za);
           }
           gausslegendre(qm, xl, xu, zs, ws);
           Hij = 0.;
           for (k=0; k<qm; k++)
             if ( df==2 )
               Hij += ws[k]*Tn( (2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj) * exp(-zs[k]/s2/l);
             else
               Hij += ws[k]*Tn( (2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                      * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
           if ( df==2 ) Hij *= exp(za/s2/l)/s2/l;
           else         Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf, ddf/2.);
           S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = Hij;
         }
       }
   }

 za = (1.-l)*cl;
 for (ii=0; ii<M; ii++)
   for (jj=0; jj<Ntilde; jj++) {
     if ( b[ii+1]<za ) S00[ ii*Ntilde+jj ] = 0.;
     else {
       if ( za<b[ii] ) xl = b[ii]; else xl = za;
       xu = b[ii+1];
       if ( df!=2 ) {
         xl = sqrt(xl-za);
         xu = sqrt(xu-za);
       }
       gausslegendre(qm, xl, xu, zs, ws);
       Hij = 0.;
       for (k=0; k<qm; k++)
         if ( df==2 )
           Hij += ws[k]*Tn((2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj) * exp(-zs[k]/s2/l);
         else
           Hij += ws[k] * Tn((2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                  * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
       if ( df==2 ) Hij *= exp(za/s2/l)/s2/l;
       else         Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf, ddf/2.);
       S00[ ii*Ntilde+jj ] = Hij;
     }
   }  
 
 for (i=0; i<NN; i++)
   for (j=0; j<NN; j++) S2s[i*NN+j] = 0.;

 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++)
     for (jj=0; jj<Ntilde; jj++)
       S2s[ (i*Ntilde+j)*NN + i*Ntilde+jj ] = Tn( (2.*zch[ i*Ntilde+j ]-b[i+1]-b[i])/(b[i+1]-b[i]), jj); 

 LU_decompose(S2s,ps,NN);

 for (n=1; n<=nmax; n++) {
   if ( n==1) {
     for (i=0; i<M; i++)
       for (j=0; j<Ntilde; j++) {
         Pns[ i*Ntilde+j ] = 0.;
         for (jj=0; jj<Ntilde; jj++)
           Pns[ i*Ntilde+j ] += /*  usual Chebyshev approximation  */
             2./Ntilde * Tn( (2.*zch[i*Ntilde+jj]-b[i+1]-b[i])/(b[i+1]-b[i]), j) * rside[ i*Ntilde+jj ];
         if ( j==0 ) Pns[ i*Ntilde+j ] /= 2.;
       }
     p00[0] = rside[NN];
   }
   else {
     for (i=0; i<NN; i++) {
       rside[i] = VF0[i] * p00[n-2];
       for (j=0; j<NN; j++) rside[i] += S1s[ i*NN+j ] * Pns[ (n-2)*NN+j ];
     }
     LU_solve2(S2s, rside, ps, NN);     
     for (i=0; i<NN; i++) Pns[ (n-1)*NN+i ] = rside[i];    
     p00[n-1] = VF0[NN] * p00[n-2];     
     for (i=0 ;i<NN; i++) p00[n-1] += S00[i] * Pns[ (n-2)*NN+i ];
   }

   p0[n-1] = 0.;
   if ( n==1 )
     p0[0] =  CHI( ddf/s2*(cu-(1.-l)*hs)/l, df);
   else
     for (j=0; j<Ntilde; j++)
       p0[n-1] += Pns[ (n-1)*NN + ihs*Ntilde+j ] * Tn( (2.*hs-b[ihs+1]-b[ihs])/(b[ihs+1]-b[ihs]), j);
     
   mn_minus = 1.; mn_plus = 0.;
   if ( n > 1) {
     for (i=0; i<M; i++)
       for (j=0; j<Ntilde; j++) {
         oben = 0.;
         unten = 0.;
         for (jj=0; jj<Ntilde; jj++) {
           oben += Pns[ (n-1)*NN + i*Ntilde+jj ] * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
           unten+= Pns[ (n-2)*NN + i*Ntilde+jj ] * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
         }
         if ( fabs(unten)<1e-16 )
           if ( fabs(oben)<1e-16 ) q = 0.;
           else q = 1.;
         else q = oben/unten;
         if ( q<mn_minus ) mn_minus = q;
         if ( q>mn_plus ) mn_plus = q;
       }
     *rho = (mn_minus + mn_plus)/2.;
     if ( fabs(mn_plus - mn_minus) < FINALeps ) {       
       *nstop = n;
       n = nmax + 1;
     }     
   } /* n > 1 */
 } /* n=1; n<=nmax; n++ */      

 Free(Pns);
 Free(zs);
 Free(ws);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);
 
 Free(S00);
 Free(p00);
 Free(VF0);

 return 0;
}


double seUR_sf_prerun_SIGMA_deluxe(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0)
{ double *ww, *zz, b1, b2, ddf2, *SF, rho, s2;
  int i, m, n, nstop, Nlocal;

 Nlocal = choose_N_for_se2(l, cl, cu);
  
 SF = vector(nmax);
 ww = vector(qm2);
 zz = vector(qm2);
 
 ddf2 = (double)(df2);  
 b1 = qCHI(     truncate/2., df2)/ddf2;
 b2 = qCHI(1. - truncate/2., df2)/ddf2; 
 gausslegendre(qm2, b1, b2, zz, ww); 
 for (i=0; i<qm2; i++) ww[i] *= ddf2 * chi( ddf2*zz[i], df2);
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm2; i++) {    
    s2 = zz[i];
    m = seUR_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop, &rho);
    if ( m != 0 ) warning("trouble with internal [package spc] function seUR_sf_deluxe");
    if ( nstop > 0 ) {
      for (n=0; n<nstop; n++)  p0[n] += ww[i] * SF[n];
      for (n=nstop; n<nmax; n++)  p0[n] += ww[i] * SF[nstop-1] * pow(rho, n-nstop+1);
    } else {
      for (n=0; n<nmax; n++)  p0[n] += ww[i] * SF[n];
    }
 }  
 
 Free(ww);
 Free(zz); 
 Free(SF);
 
 return 0;
}


double seUR_sf_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0)
{ double *ww, *zz, b1, b2, ddf2, *SF, s2;
  int i, m, n, Nlocal;

 Nlocal = choose_N_for_se2(l, cl, cu); 
  
 SF = vector(nmax);
 ww = vector(qm2);
 zz = vector(qm2);
 
 ddf2 = (double)(df2);  
 b1 = qCHI(     truncate/2., df2)/ddf2;
 b2 = qCHI(1. - truncate/2., df2)/ddf2; 
 gausslegendre(qm2, b1, b2, zz, ww); 
 for (i=0; i<qm2; i++) ww[i] *= ddf2 * chi( ddf2*zz[i], df2);
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm2; i++) {
   s2 = zz[i];
   m = seUR_sf(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF);
   if ( m != 0 ) warning("trouble with internal [package spc] function seUR_sf");
   for (n=0; n<nmax; n++)  p0[n] += ww[i] * SF[n];
 }  
 
 Free(ww);
 Free(zz); 
 Free(SF);
 
 return 0;
}


double seUR_Wq_prerun_SIGMA_deluxe(double l, double cl, double cu, double p, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate)
{ double *ww, *zz, b1, b2, ddf2, *SF, *p0, rho, *rhomany, *SFlast, Lp=-1., sf_level_adj, s2;
  int i, j, n, nstop, nstop_, nsm, nn, qnspecial=0, Nlocal;
  
 Nlocal = choose_N_for_se2(l, cl, cu); 
  
 p0 = vector(nmax);
 SF = vector(nmax); 
 rhomany = vector(qm2);
 SFlast = vector(qm2);
 ww = vector(qm2);
 zz = vector(qm2);
 
 ddf2 = (double)(df2);  
 b1 = qCHI(     truncate/2., df2)/ddf2;
 b2 = qCHI(1. - truncate/2., df2)/ddf2; 
 gausslegendre(qm2, b1, b2, zz, ww); 
 for (i=0; i<qm2; i++) ww[i] *= ddf2 * chi( ddf2*zz[i], df2);
 
 qnspecial = (qm2+1) / 2;
 
 s2 = zz[qnspecial];
 j = seUR_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nsm, &rho);
 if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
 n = nsm;
 
 if ( nsm < 1 ) { /* did not converge yet -- should be the rare case */
   nn = nmax;
   warning("The geometric tail approximation might not work.");
 } else {
   nstop = nsm;

   s2 = zz[qnspecial+1];
   j = seUR_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop_, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax ) {
       nstop = nstop_;
       i++;
       s2 = zz[qnspecial+i];
       j = seUR_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop_, &rho);
       if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }
 
   nstop = n;
   s2 = zz[qnspecial-1];
   j = seUR_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop_, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax ) {
       nstop = nstop_;
       i++;   
       s2 = zz[qnspecial-i];
       j = seUR_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop_, &rho);
       if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }  
   nn = nsm;
 }
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm2; i++) {
   s2 = zz[i];
   j = seUR_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nn, qm1, SF, &nstop, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
   if ( nstop < 1 ) {
     nstop = nn;
     warning("The geometric tail approximation might not work.");     
   }
   rhomany[i] = rho;        
   for (n=0; n<nstop; n++)  p0[n] += ww[i] * SF[n];
   if ( nstop < nn) {
     for (n=nstop; n<nn; n++) p0[n] += ww[i] * SF[nstop-1] * pow(rho, n-nstop+1);          
   }
   SFlast[i] = SF[nstop-1] * pow(rho, nn-nstop);
 }
 
 sf_level_adj = 1.-p;
 if ( p0[nn-1] <= sf_level_adj ) { 
   n = nn-1;
   while (  p0[n] <= sf_level_adj &&  n > 0  ) n--;
   if ( p0[n] > 1.-p ) Lp = (double)( n + 2 ); else Lp = 1.;
 } else {
   for (n=nn; n<nmax; n++) {
      p0[n] = 0.;
      for (i=0; i<qm2; i++) p0[n] += ww[i] * SFlast[i] * pow(rhomany[i], n-nn+1);
      if (  p0[n] <= sf_level_adj  ) {
        Lp = (double)( n + 1 );
        n = nmax+1;
      }
   }
 }
 
 Free(p0);
 Free(ww);
 Free(zz); 
 Free(SF);
 Free(SFlast);
 Free(rhomany);

 return Lp;
}


double seUR_iglarl_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate)
{ double *ww, *zz, b1, b2, result, ddf2, s2;
  int i;
  
  ww = vector(qm2);
  zz = vector(qm2);  
  ddf2 = (double)(df2);  
  b1 = qCHI(     truncate/2., df2)/ddf2;
  b2 = qCHI(1. - truncate/2., df2)/ddf2; 
  gausslegendre(qm2, b1, b2, zz, ww);  
  result = 0.;
  for (i=0; i<qm2; i++) {
    s2 = zz[i];
    result += ww[i] * ddf2 * chi( ddf2*s2, df2) * seUR_iglarl(l, s2*cl, s2*cu, s2*hs, sigma, df1, N, qm1);
  }
  Free(ww);
  Free(zz);
  
  return result;
}


double seUR_Wq(double l, double cl, double cu, double p, double hs, double sigma, int df, int N, int nmax, int qm)
{ double *S1s, *S2s, *Pns, *p0, *ws, *zs, *zch, *rside, *b, za=0., s2, ddf, xl, xu, q_minus=0., q_plus=0., dN, Hij, *S00, *p00, *VF0, mn_minus=1., mn_plus=0., oben, unten, q, enumerator=0., Wq=0.;
  int i, j, k, n, *ps, Ntilde, ihs, M, NN, ii, jj;

 M = ceil( (log(cl)-log(cu))/log(1.-l) );
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde;
 
 ihs = floor( (log(cl) - log(hs))/log(1.-l) );
 if ( ihs<0 ) ihs = 0;

 S1s = matrix(NN,NN);
 S2s = matrix(NN,NN);
 ps = ivector(NN);
 zch = matrix(M,Ntilde);
 rside = vector(NN+1);
 b   = vector(M+1);
 ws  = vector(qm);
 zs  = vector(qm);
 p0 = vector(nmax);
 Pns = matrix(nmax,NN);
 
 S00 = vector(NN);
 p00 = vector(nmax);
 VF0 = vector(NN+1);

/* interval borders b_i = cl/(1-l)^i */
 for (i=0; i<M; i++) b[i] = cl/pow(1.-l, (double)(i));
 b[M] = cu;

 /* Chebyshev nodes on [b_0,b_1],[b_1,b_2],...,[b_M-1,cu] */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     zch[ i*Ntilde+j ] = b[i] + (b[i+1]-b[i])/2.*(1.+cos(PI*(2.*j+1.)/2./dN));
   }

 /* P(L>1)(zch[i,j]) */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     rside[ i*Ntilde+j ] =  CHI( ddf/s2*(cu-(1.-l)*zch[ i*Ntilde+j ])/l, df);
   }
 rside[NN] = CHI( ddf/s2*(cu-(1.-l)*cl)/l, df); /* reflexion at cl */
 
 /* P(zch[i,j] -> zreflect) */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++)
     VF0[ i*Ntilde+j ] = CHI( ddf/s2*(cl-(1.-l)*zch[ i*Ntilde+j ])/l, df) ;
 VF0[NN] = CHI( ddf/s2*cl, df);
 
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     za = (1.-l)*zch[ i*Ntilde+j ];
     for (ii=0; ii<M; ii++)
       for (jj=0; jj<Ntilde; jj++) {
         if ( b[ii+1]<za ) S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = 0.;
         else {
           if ( za<b[ii] ) xl = b[ii]; else xl = za;
           xu = b[ii+1];
           if ( df!=2 ) {
             xl = sqrt(xl-za);
             xu = sqrt(xu-za);
           }
           gausslegendre(qm, xl, xu, zs, ws);
           Hij = 0.;
           for (k=0; k<qm; k++)
             if ( df==2 )
               Hij += ws[k]*Tn( (2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj) * exp(-zs[k]/s2/l);
             else
               Hij += ws[k]*Tn( (2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                      * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
           if ( df==2 ) Hij *= exp(za/s2/l)/s2/l;
           else         Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf, ddf/2.);
           S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = Hij;
         }
       }
   }

 za = (1.-l)*cl;
 for (ii=0; ii<M; ii++)
   for (jj=0; jj<Ntilde; jj++) {
     if ( b[ii+1]<za ) S00[ ii*Ntilde+jj ] = 0.;
     else {
       if ( za<b[ii] ) xl = b[ii]; else xl = za;
       xu = b[ii+1];
       if ( df!=2 ) {
         xl = sqrt(xl-za);
         xu = sqrt(xu-za);
       }
       gausslegendre(qm, xl, xu, zs, ws);
       Hij = 0.;
       for (k=0; k<qm; k++)
         if ( df==2 )
           Hij += ws[k]*Tn((2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj) * exp(-zs[k]/s2/l);
         else
           Hij += ws[k] * Tn((2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                  * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
       if ( df==2 ) Hij *= exp(za/s2/l)/s2/l;
       else         Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf, ddf/2.);
       S00[ ii*Ntilde+jj ] = Hij;
     }
   }  
 
 for (i=0; i<NN; i++)
   for (j=0; j<NN; j++) S2s[i*NN+j] = 0.;

 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++)
     for (jj=0; jj<Ntilde; jj++)
       S2s[ (i*Ntilde+j)*NN + i*Ntilde+jj ] = Tn( (2.*zch[ i*Ntilde+j ]-b[i+1]-b[i])/(b[i+1]-b[i]), jj); 

 LU_decompose(S2s,ps,NN);

 for (n=1; n<=nmax; n++) {
   if ( n==1) {
     for (i=0; i<M; i++)
       for (j=0; j<Ntilde; j++) {
         Pns[ i*Ntilde+j ] = 0.;
         for (jj=0; jj<Ntilde; jj++)
           Pns[ i*Ntilde+j ] += /*  usual Chebyshev approximation  */
             2./Ntilde * Tn( (2.*zch[i*Ntilde+jj]-b[i+1]-b[i])/(b[i+1]-b[i]), j) * rside[ i*Ntilde+jj ];
         if ( j==0 ) Pns[ i*Ntilde+j ] /= 2.;
       }
     p00[0] = rside[NN];
   }
   else {
     for (i=0; i<NN; i++) {
       rside[i] = VF0[i] * p00[n-2];
       for (j=0; j<NN; j++) rside[i] += S1s[ i*NN+j ] * Pns[ (n-2)*NN+j ];
     }
     LU_solve2(S2s, rside, ps, NN);     
     for (i=0; i<NN; i++) Pns[ (n-1)*NN+i ] = rside[i];    
     p00[n-1] = VF0[NN] * p00[n-2];     
     for (i=0 ;i<NN; i++) p00[n-1] += S00[i] * Pns[ (n-2)*NN+i ];
   }

   p0[n-1] = 0.;
   if ( n==1 )
     p0[0] =  CHI( ddf/s2*(cu-(1.-l)*hs)/l, df);
   else
     for (j=0; j<Ntilde; j++)
       p0[n-1] += Pns[ (n-1)*NN + ihs*Ntilde+j ] * Tn( (2.*hs-b[ihs+1]-b[ihs])/(b[ihs+1]-b[ihs]), j);

   if ( p0[n-1] < 1.-p ) {
     Wq = (double)n;
     n = nmax+1;
   } else {  
     mn_minus = 1.; mn_plus = 0.;
     if ( n > 1) {
       for (i=0; i<M; i++)
         for (j=0; j<Ntilde; j++) {
           oben = 0.;
           unten = 0.;
           for (jj=0; jj<Ntilde; jj++) {
             oben += Pns[ (n-1)*NN + i*Ntilde+jj ] * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
             unten+= Pns[ (n-2)*NN + i*Ntilde+jj ] * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
           }
           if ( fabs(unten)<1e-16 )
             if ( fabs(oben)<1e-16 ) q = 0.;
             else q = 1.;
           else q = oben/unten;
           if ( q<mn_minus ) mn_minus = q;
           if ( q>mn_plus ) mn_plus = q;
         }
       enumerator = log( (1.-p)/p0[n-1] );
       q_minus = (double)n + enumerator/log(mn_minus);
       q_plus  = (double)n + enumerator/log(mn_plus); 
       /*if ( fabs( (q_plus-q_minus)/q_minus )<FINALeps ) n = nmax+1;*/      
       if ( fabs( ceil(q_plus) - ceil(q_minus) ) < .5 ) {
	 Wq = ceil(q_plus);
	 n = nmax +1;
       }
     } /* n > 1 */
   } /* p0[n-1] >= 1.-p */
 } /* n=1; n<=nmax; n++ */      

 Free(Pns);
 Free(p0);
 Free(zs);
 Free(ws);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);
 
 Free(S00);
 Free(p00);
 Free(VF0);

 return Wq;
}


double seUR_crit(double l, double L0, double cl, double hs, double sigma, int df, int N, int qm)
{ double s1, s2, s3, ds, L1, L2, L3;

 s2 = hs;
 do {
   s1 = s2;
   L1 = L2;
   s2 += .2;
   L2 = seUR_iglarl(l, cl, s2, hs, sigma, df, N, qm);
 } while (L2<L0);
 do {
   s1 = s2;
   L1 = L2;
   s2 -= .02;
   L2 = seUR_iglarl(l, cl, s2, hs, sigma, df, N, qm);
 } while (L2>L0);

 do {
   s3 = s1 + (L0-L1)/(L2-L1) * (s2-s1);
   L3 = seUR_iglarl(l, cl, s3, hs, sigma, df, N, qm);
   ds = s3-s2; s1 = s2; L1 = L2; s2 = s3; L2 = L3;
 } while ( fabs(L0-L3)>1e-6 && fabs(ds)>1e-7 );

 return s3;
}


double stdeUR_crit(double l, double L0, double cl, double hs, double sigma, int df, int N, int qm)
{ double s1, s2, s3, ds, L1, L2, L3;

 s2 = hs;
 do {
   s1 = s2;
   L1 = L2;
   s2 += .2;
   L2 = stdeUR_iglarl(l, cl, s2, hs, sigma, df, N, qm);
 } while ( L2<L0 );
 do {
   s1 = s2;
   L1 = L2;
   s2 -= .02;
   L2 = stdeUR_iglarl(l, cl, s2, hs, sigma, df, N, qm);
 } while ( L2>L0 );

 do {
   s3 = s1 + (L0-L1)/(L2-L1) * (s2-s1);
   L3 = stdeUR_iglarl(l, cl, s3, hs, sigma, df, N, qm);
   ds = s3-s2; s1 = s2; L1 = L2; s2 = s3; L2 = L3;
 } while ( fabs(L0-L3)>1e-7 && fabs(ds)>1e-8 );

 return s3;
}


double seUR_crit_prerun_SIGMA(double l, double L0, double cl, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate)
{ double s1, s2, s3, ds, L1=0., L2=0., L3=0.;

 s2 = hs;
 do {
   s1 = s2;
   L1 = L2;
   s2 += .2;   
   L2 = seUR_iglarl_prerun_SIGMA(l, cl, s2, hs, sigma, df1, df2, N, qm1, qm2, truncate);
 } while ( L2 < L0 );
 do {
   s1 = s2;
   L1 = L2;
   s2 -= .02;   
   L2 = seUR_iglarl_prerun_SIGMA(l, cl, s2, hs, sigma, df1, df2, N, qm1, qm2, truncate);
 } while ( L2 > L0 );


 do {
   s3 = s1 + (L0-L1)/(L2-L1) * (s2-s1);
   L3 = seUR_iglarl_prerun_SIGMA(l, cl, s3, hs, sigma, df1, df2, N, qm1, qm2, truncate);
   ds = s3-s2; s1 = s2; L1 = L2; s2 = s3; L2 = L3;
 } while ( fabs(L0-L3)>1e-6 && fabs(ds)>1e-7 );

 return s3;
}


double seUR_q_crit(double l, int L0, double alpha, double cl, double hs, double sigma, int df, int N, int qm, double c_error, double a_error)
{ double s1, s2, s3, ds, p1, p2, p3, *SF;
  int result=1;

 SF  = vector(L0); 
  
 s2 = hs; p2 = 1.;
 do {
   p1 = p2;
   s2 += .2;
   result = seUR_sf(l, cl, s2, hs, sigma, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in seUR_q_crit [package spc]");
   p2 = 1. - SF[L0-1];
 } while ( p2 > alpha );

 s1 = s2 - .2;

 do {
   s3 = s1 + (alpha - p1)/( p2 - p1 ) * (s2-s1);
   result = seUR_sf(l, cl, s3, hs, sigma, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in seUR_q_crit [package spc]");
   p3 = 1. - SF[L0-1];
   ds = s3 - s2; s1 = s2; p1 = p2; s2 = s3; p2 = p3;
 } while ( fabs(alpha - p3)>a_error && fabs(ds)>c_error );

 Free(SF);
 
 return s3;
}


double seUR_q_crit_prerun_SIGMA(double l, int L0, double alpha, double cl, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate, int tail_approx, double c_error, double a_error)
{ double s1, s2, s3, ds, p1, p2, p3, *SF;
  int result=1;  
  
 SF  = vector(L0); 
 
 s2 = seUR_q_crit(l, L0, alpha, cl, hs, sigma, df1, N, qm1, c_error, a_error); 
 if ( tail_approx ) result = seUR_sf_prerun_SIGMA_deluxe(l, cl, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
 else result = seUR_sf_prerun_SIGMA(l, cl, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
 if ( result != 0 ) warning("trouble in seUR_q_crit_prerun_SIGMA [package spc]");
 p2 = 1. - SF[L0-1];
 
 if ( p2 > alpha ) {
   do {
     p1 = p2;
     s2 += .2;
     if ( tail_approx ) result = seUR_sf_prerun_SIGMA_deluxe(l, cl, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     else result = seUR_sf_prerun_SIGMA(l, cl, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     if ( result != 0 ) warning("trouble in seUR_q_crit_prerun_SIGMA [package spc]");
     p2 = 1. - SF[L0-1];
   } while ( p2 > alpha );
   s1 = s2 - .2;
 } else {
   do {
     p1 = p2;
     s2 -= .2;
     if ( tail_approx ) result = seUR_sf_prerun_SIGMA_deluxe(l, cl, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     else result = seUR_sf_prerun_SIGMA(l, cl, s2, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     if ( result != 0 ) warning("trouble in seUR_q_crit_prerun_SIGMA [package spc]");
     p2 = 1. - SF[L0-1];
   } while ( p2 <= alpha && s2 > hs );
   s1 = s2 + .2;
 }

 do {
   s3 = s1 + (alpha - p1)/( p2 - p1 ) * (s2-s1);
   if ( tail_approx ) result = seUR_sf_prerun_SIGMA_deluxe(l, cl, s3, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
   else result = seUR_sf_prerun_SIGMA(l, cl, s3, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
   if ( result != 0 ) warning("trouble in seUR_q_crit_prerun_SIGMA [package spc]");
   p3 = 1. - SF[L0-1];
   ds = s3 - s2; s1 = s2; p1 = p2; s2 = s3; p2 = p3;
 } while ( fabs(alpha - p3)>a_error && fabs(ds)>c_error );

 Free(SF);
 
 return s3;
}


double seLR_iglarl(double l, double cl, double cu, double hs, double sigma, int df, int N, int qm)
{ double *a, *g, *w, *z, *t, h, arl, Hij, xl, za, dN, ddf, s2, t0, t1, x0, x1, dummy;
  int i, j, k, qi, qj, M, Ntilde, NN, ii, it, jj;

 M = ceil( (log(cl)-log(cu))/log(1.-l) );
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde - 1.;

 a = matrix(NN, NN);
 g = vector(NN);
 t = vector(NN);
 w = vector(qm);
 z = vector(qm);

 for(i=0;i<M;i++) {
   t0 = cl/pow(1.-l,(double)(i));
   t1 = t0/(1.-l);
   if (t1>cu) t1 = cu;

   for (j=1;j<Ntilde;j++) { /* node_i,Ntilde-1 = node_i+1,0 */
     h = cos( PI/dN *(dN-j) );
     t[i*(Ntilde-1)+j] = t0 + (h+1.)/2.*(t1-t0);
     /* Chebyshev Gauss-Lobatto nodes on [t0,t1] */
   }
 }
 t[0] = cl;

 for (i=0;i<M;i++) {
   for (j=1;j<=Ntilde;j++) {
     ii = i*Ntilde + j-1;
     it = i*(Ntilde-1) + j-1;

     za = (1.-l)*t[it];
     if (za<cl) xl = cl; else xl = za;

     for (qi=0;qi<i-1;qi++)
       for (qj=1;qj<=Ntilde;qj++) {
         jj = qi*Ntilde + qj-1;
         a[ii*NN+jj] = 0.;
       }

     if (i>0) {
       qi = i-1;
       t0 = cl/pow(1.-l,(double)qi);
       t1 = t0/(1.-l);
       if (t1>cu) t1 = cu;
       if (t0<xl) x0 = xl; else x0 = t0;
       if (df==2)
         x1 = t1;
       else {
         if (x0-za>1e-10) x0 = sqrt(x0-za); else x0 = 0.;
         if (t1-za>1e-10) x1 = sqrt(t1-za); else x1 = 0.;
       }

       for (qj=1;qj<=Ntilde;qj++) {
         jj = qi*Ntilde + qj-1;

         if (j==1) a[ii*NN+jj] = - Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1);
         else {
           if (fabs(t1-x0)>1e-8) {
             gausslegendre(qm,x0,x1,z,w);
             Hij = 0.;
             for (k=0;k<qm;k++) {
               if (df==2)
                 Hij += w[k] * Tn( (2.*z[k]-t0-t1)/(t1-t0), qj-1) *
                        exp((za-z[k])/s2/l);
               if (df!=2)
                 Hij += w[k] * Tn( (2.*(z[k]*z[k]+za)-t0-t1)/(t1-t0) ,qj-1) *
                        2. * pow(z[k], ddf-1.) * exp(-ddf*z[k]*z[k]/2./s2/l);
             }
             if (df==2) Hij /= s2*l;
             else       Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf,ddf/2.);
             a[ii*NN+jj] = -Hij;
           }
           else a[ii*NN+jj] = 0.;
         }
       }
     }

     for (qi=i;qi<M;qi++) {
       t0 = cl/pow(1.-l,(double)qi);
       t1 = t0/(1.-l);
       if (t1>cu) t1 = cu;
       if (t0<xl) x0 = xl; else x0 = t0;
       if (df==2)
         x1 = t1;
       else {
        if (x0-za>1e-10) x0 = sqrt(x0-za); else x0 = 0.;
        if (t1-za>1e-10) x1 = sqrt(t1-za); else x1 = 0.;
       }

       if (i>0 && j==1 && qi==i) {
         for (qj=1;qj<=Ntilde;qj++) {
           jj = qi*Ntilde + qj-1;
           a[ii*NN+jj] = Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1);
         }
       }

       if (i>0 && j==1 && qi>i) {
         for (qj=1;qj<=Ntilde;qj++) {
           jj = qi*Ntilde + qj-1;
           a[ii*NN+jj] = 0.;
         }
       }

       if (i==0 || j>1) {
         for (qj=1;qj<=Ntilde;qj++) {
           jj = qi*Ntilde + qj-1;
           gausslegendre(qm,x0,x1,z,w);
           Hij = 0.;
           for (k=0;k<qm;k++) {
             if (df==2)
               Hij += w[k] * Tn( (2.*z[k]-t0-t1)/(t1-t0), qj-1) *
                      exp((za-z[k])/s2/l);
             if (df!=2)
               Hij += w[k] * Tn( (2.*(z[k]*z[k]+za)-t0-t1)/(t1-t0),qj-1) *
                      2. * pow(z[k], ddf-1.) * exp(-ddf*z[k]*z[k]/2./s2/l);
           }
           if (df==2) Hij /= s2*l;
           else       Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf,ddf/2.);
           if (qi==i) a[ii*NN+jj] = Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1) -
                                        Hij;
           else a[ii*NN+jj] = -Hij;
         }
       }
     }

/*   "reflection area" */
     if (i==0 || j>1) {
       t0 = cl/pow(1.-l, (double)(M-1.));
       t1 = cu;
       for (qj=1;qj<=Ntilde;qj++) {
         dummy = (cu-za)/l/s2;
         if (dummy>0.) {
           if (df==1) dummy = 2.*( 1. - PHI( sqrt(dummy), 0. ) );
           if (df==2) dummy = exp( -dummy );
           if (df>2)  dummy = 1. - CHI( df*dummy, df);
         }
         else dummy = 0.;
         jj = (M-1)*Ntilde + qj-1;
         a[ii*NN+jj] -= dummy;
       }
     }
   }
 }

 for (j=0;j<NN;j++) g[j] = 1.;
 for (j=1;j<M;j++) g[Ntilde*j] = 0.;

 LU_solve(a,g,NN);

 arl = 0.;
 for (i=0;i<M;i++) {
   t0 = cl/pow(1.-l,(double)i);
   t1 = t0/(1.-l);
   if (t1>cu) t1 = cu;

   if (t0<hs && hs<=t1)
     for (j=1;j<=Ntilde;j++) {
        ii = i*Ntilde + j-1;
        arl += g[ii] * Tn((2.*hs-t0-t1)/(t1-t0),j-1);
     }
 }

 Free(z);
 Free(w);
 Free(t);
 Free(g);
 Free(a);

 return arl;
}


double stdeLR_iglarl(double l, double cl, double cu, double hs, double sigma, int df, int N, int qm)
{ double *a, *g, *w, *z, *t, h, arl, Hij, xl, za, dN, ddf, s2, t0, t1, x0, x1, dummy, v;
  int i, j, k, qi, qj, M, Ntilde, NN, ii, it, jj;

 M = ceil( (log(cl)-log(cu))/log(1.-l) );
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde - 1.;

 a = matrix(NN, NN);
 g = vector(NN);
 t = vector(NN);
 w = vector(qm);
 z = vector(qm);

 for(i=0; i<M; i++) {
   t0 = cl/pow(1.-l,(double)(i));
   t1 = t0/(1.-l);
   if ( t1>cu ) t1 = cu;

   for (j=1; j<Ntilde; j++) { /* node_i,Ntilde-1 = node_i+1,0 */
     h = cos( PI/dN *(dN-j) );
     t[i*(Ntilde-1)+j] = t0 + (h+1.)/2.*(t1-t0); /* Chebyshev Gauss-Lobatto nodes on [t0,t1] */
   }
 }
 t[0] = cl;

 for (i=0; i<M; i++) {
   for (j=1; j<=Ntilde; j++) {
     ii = i*Ntilde + j-1;
     it = i*(Ntilde-1) + j-1;

     za = (1.-l)*t[it];
     if ( za<cl ) xl = cl; else xl = za;

     for (qi=0; qi<i-1; qi++)
       for (qj=1; qj<=Ntilde; qj++) {
         jj = qi*Ntilde + qj-1;
         a[ii*NN+jj] = 0.;
       }

     if ( i>0 ) {
       qi = i-1;
       t0 = cl/pow(1.-l,(double)qi);
       t1 = t0/(1.-l);
       if ( t1>cu ) t1 = cu;
       if ( t0<xl ) x0 = xl; else x0 = t0;
       x1 = t1;

       for (qj=1; qj<=Ntilde; qj++) {
         jj = qi*Ntilde + qj-1;

         if ( j==1 ) a[ii*NN+jj] = - Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1);
         else {
           if ( fabs(t1-x0)>1e-8 ) {
             gausslegendre(qm, x0, x1, z, w);
             Hij = 0.;
             for (k=0; k<qm; k++) {
	       v = (z[k] - za) / l;
               Hij += w[k] * Tn( (2.*z[k]-t0-t1)/(t1-t0), qj-1) * pow(v,ddf-1.)*exp(-ddf/2./s2*v*v);
             }
             Hij *= 2./l/gammafn(ddf/2.)/pow(2.*s2/ddf,ddf/2.);
             a[ii*NN+jj] = -Hij;
           }
           else a[ii*NN+jj] = 0.;
         }
       }
     }

     for (qi=i; qi<M; qi++) {
       t0 = cl/pow(1.-l,(double)qi);
       t1 = t0/(1.-l);
       if ( t1>cu ) t1 = cu;
       if ( t0<xl ) x0 = xl; else x0 = t0;
       x1 = t1;

       if ( i>0 && j==1 && qi==i ) {
         for (qj=1; qj<=Ntilde; qj++) {
           jj = qi*Ntilde + qj-1;
           a[ii*NN+jj] = Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1);
         }
       }

       if ( i>0 && j==1 && qi>i ) {
         for (qj=1; qj<=Ntilde; qj++) {
           jj = qi*Ntilde + qj-1;
           a[ii*NN+jj] = 0.;
         }
       }

       if ( i==0 || j>1 ) {
         for (qj=1; qj<=Ntilde; qj++) {
           jj = qi*Ntilde + qj-1;
           gausslegendre(qm, x0, x1, z, w);
           Hij = 0.;
           for (k=0; k<qm; k++) {
	     v = (z[k] - za) / l;
             Hij += w[k] * Tn( (2.*z[k]-t0-t1)/(t1-t0), qj-1) * pow(v,ddf-1.)*exp(-ddf/2./s2*v*v);
           }
           Hij *= 2./l/gammafn(ddf/2.)/pow(2.*s2/ddf,ddf/2.);
           if ( qi==i ) a[ii*NN+jj] = Tn((2.*t[it]-t0-t1)/(t1-t0),qj-1) - Hij;
           else a[ii*NN+jj] = -Hij;
         }
       }
     }

/*   "reflection area" */
     if ( i==0 || j>1 ) {
       t0 = cl/pow(1.-l, (double)(M-1.));
       t1 = cu;
       for (qj=1; qj<=Ntilde; qj++) {
	 dummy = 0.;
         v = (cu-za)/l;
         if ( v>0. ) dummy = 1. - CHI( ddf/s2*v*v, df);
         jj = (M-1)*Ntilde + qj-1;
         a[ii*NN+jj] -= dummy;
       }
     }
   }
 }

 for (j=0; j<NN; j++) g[j] = 1.;
 for (j=1; j<M; j++) g[Ntilde*j] = 0.;

 LU_solve(a,g,NN);

 arl = 0.;
 for (i=0; i<M; i++) {
   t0 = cl/pow(1.-l,(double)i);
   t1 = t0/(1.-l);
   if ( t1>cu ) t1 = cu;

   if ( t0<hs && hs<=t1 )
     for (j=1; j<=Ntilde; j++) {
        ii = i*Ntilde + j-1;
        arl += g[ii] * Tn((2.*hs-t0-t1)/(t1-t0),j-1);
     }
 }

 Free(z);
 Free(w);
 Free(t);
 Free(g);
 Free(a);

 return arl;
}


double lns2ewmaU_arl_igl(double l, double cl, double cu, double hs, double sigma, int df, int N)
{ double *a, *g, *w, *z, arl, lns, ddf, s2;
  int i, j, NN;

 NN = N + 1;
 s2 = sigma*sigma;
 ddf = (double)df;

 a = matrix(NN, NN);
 g = vector(NN);
 w = vector(N);
 z = vector(N);

 gausslegendre(N, cl, cu, z, w);

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) {
     lns = exp( (z[j]-(1.-l)*z[i])/l );
     a[i*NN+j] = -w[j]/l * chi( ddf/s2*lns, df)*ddf/s2*lns;
   }
   ++a[i*NN+i];
   lns = exp( (cl-(1.-l)*z[i])/l );
   a[i*NN+NN-1] = -CHI( ddf/s2*lns, df);
 }
 
 for (j=0; j<N; j++) {
   lns = exp( (z[j]-(1.-l)*cl)/l );
   a[N*NN+j] = -w[j]/l * chi( ddf/s2*lns, df)*ddf/s2*lns;
 }
 a[N*NN+N] = 1. - CHI( ddf/s2*exp(cl), df); 

 for (j=0; j<NN; j++) g[j] = 1.;
 LU_solve(a, g, NN);

 lns = exp(  (cl-(1.-l)*hs)/l );
 arl = 1. + CHI( ddf/s2*lns, df) * g[N];
 for (j=0; j<N; j++) {
   lns = exp( (z[j]-(1.-l)*hs)/l );
   arl += w[j]/l * chi( ddf/s2*lns, df)*ddf/s2*lns * g[j];
 }

 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}


double lns2ewmaU_crit(double l, double L0, double cl, double hs, double sigma, int df, int N)
{ double s1, s2, s3, ds, L1, L2, L3;

 L2 = 1.; 
 s2 = hs;
 do {
   s1 = s2;
   L1 = L2;
   s2 += .1;
   L2 = lns2ewmaU_arl_igl(l,cl,s2,hs,sigma,df,N);
 } while ( L2<L0 );
 
 if ( L2 > 10.*L0 ) {
   do {
     s1 = s2;
     L1 = L2;
     s2 -= .01;
     L2 = lns2ewmaU_arl_igl(l,cl,s2,hs,sigma,df,N);
   } while ( L2>L0 );
 }

 do {
   s3 = s1 + (L0-L1)/(L2-L1) * (s2-s1);
   L3 = lns2ewmaU_arl_igl(l,cl,s3,hs,sigma,df,N);
   ds = s3-s2; s1 = s2; L1 = L2; s2 = s3; L2 = L3;
 } while ( fabs(L0-L3)>1e-7 && fabs(ds)>1e-8 );

 return s3;
}


double lns2ewma2_arl_igl(double l, double cl, double cu, double hs, double sigma, int df, int N)
{ double *a, *g, *w, *z, arl, lns, ddf, s2;
  int i, j;

 s2 = sigma*sigma;
 ddf = (double)df;

 a = matrix(N,N);
 g = vector(N);
 w = vector(N);
 z = vector(N);

 gausslegendre(N, cl, cu, z, w);

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) {
     lns = exp( (z[j]-(1.-l)*z[i])/l );
     a[i*N+j] = -w[j]/l * chi( ddf/s2*lns, df)*ddf/s2*lns;
   }
   ++a[i*N+i];
 }

 for (j=0;j<N;j++) g[j] = 1.;
 LU_solve(a,g,N);

 arl = 1.;
 for (j=0; j<N; j++) {
   lns = exp( (z[j]-(1.-l)*hs)/l );
   arl += w[j]/l * chi( ddf/s2*lns, df)*ddf/s2*lns * g[j];
 }

 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}


double lns2ewma2_crit_cufix(double l, double cu, double L0, double hs, double sigma, int df, int N)
{ double s1, s2, s3, ds, L1, L2, L3;

 L2 = 1.;
 s2 = hs;
 do {
   s1 = s2;
   L1 = L2;
   s2 -= .1;
   L2 = lns2ewma2_arl_igl(l,s2,cu,hs,sigma,df,N);
 } while ( L2<L0 );
 
  if ( L2 > 10.*L0 ) {
   do {
     s1 = s2;
     L1 = L2;
     s2 += .01;
     L2 = lns2ewma2_arl_igl(l,s2,cu,hs,sigma,df,N);
   } while ( L2>L0 );
 }

 do {
   s3 = s1 + (L0-L1)/(L2-L1) * (s2-s1);
   L3 = lns2ewma2_arl_igl(l,s3,cu,hs,sigma,df,N);
   ds = s3-s2; s1 = s2; L1 = L2; s2 = s3; L2 = L3;
 } while ( fabs(L0-L3)>1e-7 && fabs(ds)>1e-8 );

 return s3;
}


int lns2ewma2_crit_unbiased(double l, double L0, double *cl, double *cu, double hs, double sigma, int df, int N)
{ double s1, s2, s3, ds, sl1, sl2, sl3, csl, Lm, Lp, mitte, ddf;

 ddf = (double)df;
 /*mitte = -1./ddf - 1./3./ddf/ddf + 2./15./ddf/ddf/ddf/ddf;*/
 mitte = E_log_gamma(ddf);
 
 csl = lns2ewma2_crit_sym(l, L0, hs, sigma, df, N);
 s1 = 2.*mitte - csl;
 Lm = lns2ewma2_arl_igl(l,csl,s1,hs,sigma-lmEPS,df,N);
 Lp = lns2ewma2_arl_igl(l,csl,s1,hs,sigma+lmEPS,df,N);
 sl1 = (Lp-Lm)/(2.*lmEPS);
 
 do {
   s2 = s1;
   sl2 = sl1;
   s1 -= .1;
   csl = lns2ewma2_crit_cufix(l,s1,L0,hs,sigma,df,N);
   Lm = lns2ewma2_arl_igl(l,csl,s1,hs,sigma-lmEPS,df,N);
   Lp = lns2ewma2_arl_igl(l,csl,s1,hs,sigma+lmEPS,df,N);
   sl1 = (Lp-Lm)/(2.*lmEPS);
 } while ( sl1>0. );

 do {
   s3 = s1 - sl1/(sl2-sl1) * (s2-s1);
   csl = lns2ewma2_crit_cufix(l,s3,L0,hs,sigma,df,N);
   Lm = lns2ewma2_arl_igl(l,csl,s3,hs,sigma-lmEPS,df,N);
   Lp = lns2ewma2_arl_igl(l,csl,s3,hs,sigma+lmEPS,df,N);
   sl3 = (Lp-Lm)/(2.*lmEPS);
   ds = s3-s2; s1 = s2; sl1 = sl2; s2 = s3; sl2 = sl3;
 } while ( fabs(sl3)>1e-7 && fabs(ds)>1e-8 );

 *cl = csl; *cu = s3;

 return 0;
}


double lns2ewma2_crit_sym(double l, double L0, double hs, double sigma, int df, int N)
{ double cu, cl1, cl2, cl3, L1, L2, L3, dl, mitte, ddf;

  ddf = (double)df;
  /*mitte = -1./ddf - 1./3./ddf/ddf + 2./15./ddf/ddf/ddf/ddf;*/
  mitte = E_log_gamma(ddf);
  
  L2 = 1.;
  cl2 = mitte;
  do {
    cl1 = cl2;
    L1 = L2;
    cl2 -= .1;
    cu = 2.*mitte - cl2;
    L2 = lns2ewma2_arl_igl(l, cl2, cu, hs, sigma, df, N);
  } while ( L2<L0 );

  do {
    cl3 = cl1 + (L0-L1)/(L2-L1) * (cl2-cl1);
    cu = 2.*mitte - cl3;
    L3 = lns2ewma2_arl_igl(l, cl3, cu, hs, sigma, df, N);
    dl = cl3-cl2; cl1 = cl2; L1 = L2; cl2 = cl3; L2 = L3;
    if ( L3 < 1. ) error("invalid ARL value");
  } while ( (fabs(L0-L3)>1e-7) && (fabs(dl)>1e-8) ); 

  return cl3;
}


double seLR_sf(double l, double cl, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0)
{ double *S1s, *S2s, *Pns, *ws, *zs, *zch, *rside, *b, za=0., s2, ddf, xl, xu, dN, Hij, *S00, *p00, *VF0;
  int i, j, k, n, *ps, Ntilde, ihs, M, NN, ii, jj;

 M = ceil( (log(cl)-log(cu))/log(1.-l) );
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde;
 
 ihs = floor( (log(cl) - log(hs))/log(1.-l) );
 if ( ihs<0 ) ihs = 0;

 S1s = matrix(NN,NN);
 S2s = matrix(NN,NN);
 ps = ivector(NN);
 zch = matrix(M,Ntilde);
 rside = vector(NN+1);
 b   = vector(M+1);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,NN);
 
 S00 = vector(NN);
 p00 = vector(nmax);
 VF0 = vector(NN+1);

/* interval borders b_i = cl/(1-l)^i */
 for (i=0; i<M; i++) b[i] = cl/pow(1.-l, (double)(i));
 b[M] = cu;

 /* Chebyshev nodes on [b_0,b_1],[b_1,b_2],...,[b_M-1,cu] */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     zch[ i*Ntilde+j ] = b[i] + (b[i+1]-b[i])/2.*(1.+cos(PI*(2.*j+1.)/2./dN));
   }

 /* P(L>1)(zch[i,j]) */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     rside[ i*Ntilde+j ] =  1. - CHI( ddf/s2*(cl-(1.-l)*zch[ i*Ntilde+j ])/l, df);
   }
 rside[NN] = 1. - CHI( ddf/s2*(cl-(1.-l)*cu)/l, df); /* reflexion at cu */
 
 /* P(zch[i,j] -> zreflect) */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++)
     VF0[ i*Ntilde+j ] = 1. - CHI( ddf/s2*(cu-(1.-l)*zch[ i*Ntilde+j ])/l, df) ;
 VF0[NN] = 1. - CHI( ddf/s2*cu, df);
 
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     za = (1.-l)*zch[ i*Ntilde+j ];
     for (ii=0; ii<M; ii++)
       for (jj=0; jj<Ntilde; jj++) {
         if ( b[ii+1]<za ) S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = 0.;
         else {
           if ( za<b[ii] ) xl = b[ii]; else xl = za;
           xu = b[ii+1];
           if ( df!=2 ) {
             xl = sqrt(xl-za);
             xu = sqrt(xu-za);
           }
           gausslegendre(qm, xl, xu, zs, ws);
           Hij = 0.;
           for (k=0; k<qm; k++)
             if ( df==2 )
               Hij += ws[k]*Tn( (2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj) * exp(-zs[k]/s2/l);
             else
               Hij += ws[k]*Tn( (2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                      * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
           if ( df==2 ) Hij *= exp(za/s2/l)/s2/l;
           else         Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf, ddf/2.);
           S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = Hij;
         }
       }
   }

 za = (1.-l)*cu;
 for (ii=0; ii<M; ii++)
   for (jj=0; jj<Ntilde; jj++) {
     if ( b[ii+1]<za ) S00[ ii*Ntilde+jj ] = 0.;
     else {
       if ( za<b[ii] ) xl = b[ii]; else xl = za;
       xu = b[ii+1];
       if ( df!=2 ) {
         xl = sqrt(xl-za);
         xu = sqrt(xu-za);
       }
       gausslegendre(qm, xl, xu, zs, ws);
       Hij = 0.;
       for (k=0; k<qm; k++)
         if ( df==2 )
           Hij += ws[k]*Tn((2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj) * exp(-zs[k]/s2/l);
         else
           Hij += ws[k] * Tn((2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                  * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
       if ( df==2 ) Hij *= exp(za/s2/l)/s2/l;
       else         Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf, ddf/2.);
       S00[ ii*Ntilde+jj ] = Hij;
     }
   }  
 
 for (i=0; i<NN; i++)
   for (j=0; j<NN; j++) S2s[i*NN+j] = 0.;

 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++)
     for (jj=0; jj<Ntilde; jj++)
       S2s[ (i*Ntilde+j)*NN + i*Ntilde+jj ] = Tn( (2.*zch[ i*Ntilde+j ]-b[i+1]-b[i])/(b[i+1]-b[i]), jj); 

 LU_decompose(S2s,ps,NN);

 for (n=1; n<=nmax; n++) {
   if ( n==1) {
     for (i=0; i<M; i++)
       for (j=0; j<Ntilde; j++) {
         Pns[ i*Ntilde+j ] = 0.;
         for (jj=0; jj<Ntilde; jj++)
           Pns[ i*Ntilde+j ] += /*  usual Chebyshev approximation  */
             2./Ntilde * Tn( (2.*zch[i*Ntilde+jj]-b[i+1]-b[i])/(b[i+1]-b[i]), j) * rside[ i*Ntilde+jj ];
         if ( j==0 ) Pns[ i*Ntilde+j ] /= 2.;
       }
     p00[0] = rside[NN];
   }
   else {
     for (i=0; i<NN; i++) {
       rside[i] = VF0[i] * p00[n-2];
       for (j=0; j<NN; j++) rside[i] += S1s[ i*NN+j ] * Pns[ (n-2)*NN+j ];
     }
     LU_solve2(S2s, rside, ps, NN);     
     for (i=0; i<NN; i++) Pns[ (n-1)*NN+i ] = rside[i];    
     p00[n-1] = VF0[NN] * p00[n-2];     
     for (i=0 ;i<NN; i++) p00[n-1] += S00[i] * Pns[ (n-2)*NN+i ];
   }

   p0[n-1] = 0.;
   if ( n==1 )
     p0[0] =  1. - CHI( ddf/s2*(cl-(1.-l)*hs)/l, df);
   else
     for (j=0; j<Ntilde; j++)
       p0[n-1] += Pns[ (n-1)*NN + ihs*Ntilde+j ] * Tn( (2.*hs-b[ihs+1]-b[ihs])/(b[ihs+1]-b[ihs]), j);
 }
 
 Free(Pns);
 Free(zs);
 Free(ws);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);
 
 Free(S00);
 Free(p00);
 Free(VF0);

 return 0;
}


double seLR_sf_deluxe(double l, double cl, double cu, double hs, double sigma, int df, int N, int nmax, int qm, double *p0, int *nstop, double *rho)
{ double *S1s, *S2s, *Pns, *ws, *zs, *zch, *rside, *b, za=0., s2, ddf, xl, xu, dN, Hij, *S00, *p00, *VF0, mn_minus=1., mn_plus=0., oben, unten, q;
  int i, j, k, n, *ps, Ntilde, ihs, M, NN, ii, jj;

 M = ceil( (log(cl)-log(cu))/log(1.-l) );
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde;
 
 ihs = floor( (log(cl) - log(hs))/log(1.-l) );
 if ( ihs<0 ) ihs = 0;

 S1s = matrix(NN,NN);
 S2s = matrix(NN,NN);
 ps = ivector(NN);
 zch = matrix(M,Ntilde);
 rside = vector(NN+1);
 b   = vector(M+1);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,NN);
 
 S00 = vector(NN);
 p00 = vector(nmax);
 VF0 = vector(NN+1);
 
 *nstop = 0;

/* interval borders b_i = cl/(1-l)^i */
 for (i=0; i<M; i++) b[i] = cl/pow(1.-l, (double)(i));
 b[M] = cu;

 /* Chebyshev nodes on [b_0,b_1],[b_1,b_2],...,[b_M-1,cu] */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     zch[ i*Ntilde+j ] = b[i] + (b[i+1]-b[i])/2.*(1.+cos(PI*(2.*j+1.)/2./dN));
   }

 /* P(L>1)(zch[i,j]) */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     rside[ i*Ntilde+j ] =  1. - CHI( ddf/s2*(cl-(1.-l)*zch[ i*Ntilde+j ])/l, df);
   }
 rside[NN] = 1. - CHI( ddf/s2*(cl-(1.-l)*cu)/l, df); /* reflexion at cu */
 
 /* P(zch[i,j] -> zreflect) */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++)
     VF0[ i*Ntilde+j ] = 1. - CHI( ddf/s2*(cu-(1.-l)*zch[ i*Ntilde+j ])/l, df) ;
 VF0[NN] = 1. - CHI( ddf/s2*cu, df);
 
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     za = (1.-l)*zch[ i*Ntilde+j ];
     for (ii=0; ii<M; ii++)
       for (jj=0; jj<Ntilde; jj++) {
         if ( b[ii+1]<za ) S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = 0.;
         else {
           if ( za<b[ii] ) xl = b[ii]; else xl = za;
           xu = b[ii+1];
           if ( df!=2 ) {
             xl = sqrt(xl-za);
             xu = sqrt(xu-za);
           }
           gausslegendre(qm, xl, xu, zs, ws);
           Hij = 0.;
           for (k=0; k<qm; k++)
             if ( df==2 )
               Hij += ws[k]*Tn( (2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj) * exp(-zs[k]/s2/l);
             else
               Hij += ws[k]*Tn( (2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                      * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
           if ( df==2 ) Hij *= exp(za/s2/l)/s2/l;
           else         Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf, ddf/2.);
           S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = Hij;
         }
       }
   }

 za = (1.-l)*cu;
 for (ii=0; ii<M; ii++)
   for (jj=0; jj<Ntilde; jj++) {
     if ( b[ii+1]<za ) S00[ ii*Ntilde+jj ] = 0.;
     else {
       if ( za<b[ii] ) xl = b[ii]; else xl = za;
       xu = b[ii+1];
       if ( df!=2 ) {
         xl = sqrt(xl-za);
         xu = sqrt(xu-za);
       }
       gausslegendre(qm, xl, xu, zs, ws);
       Hij = 0.;
       for (k=0; k<qm; k++)
         if ( df==2 )
           Hij += ws[k]*Tn((2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj) * exp(-zs[k]/s2/l);
         else
           Hij += ws[k] * Tn((2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                  * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
       if ( df==2 ) Hij *= exp(za/s2/l)/s2/l;
       else         Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf, ddf/2.);
       S00[ ii*Ntilde+jj ] = Hij;
     }
   }  
 
 for (i=0; i<NN; i++)
   for (j=0; j<NN; j++) S2s[i*NN+j] = 0.;

 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++)
     for (jj=0; jj<Ntilde; jj++)
       S2s[ (i*Ntilde+j)*NN + i*Ntilde+jj ] = Tn( (2.*zch[ i*Ntilde+j ]-b[i+1]-b[i])/(b[i+1]-b[i]), jj); 

 LU_decompose(S2s,ps,NN);

 for (n=1; n<=nmax; n++) {
   if ( n==1) {
     for (i=0; i<M; i++)
       for (j=0; j<Ntilde; j++) {
         Pns[ i*Ntilde+j ] = 0.;
         for (jj=0; jj<Ntilde; jj++)
           Pns[ i*Ntilde+j ] += /*  usual Chebyshev approximation  */
             2./Ntilde * Tn( (2.*zch[i*Ntilde+jj]-b[i+1]-b[i])/(b[i+1]-b[i]), j) * rside[ i*Ntilde+jj ];
         if ( j==0 ) Pns[ i*Ntilde+j ] /= 2.;
       }
     p00[0] = rside[NN];
   }
   else {
     for (i=0; i<NN; i++) {
       rside[i] = VF0[i] * p00[n-2];
       for (j=0; j<NN; j++) rside[i] += S1s[ i*NN+j ] * Pns[ (n-2)*NN+j ];
     }
     LU_solve2(S2s, rside, ps, NN);     
     for (i=0; i<NN; i++) Pns[ (n-1)*NN+i ] = rside[i];    
     p00[n-1] = VF0[NN] * p00[n-2];     
     for (i=0 ;i<NN; i++) p00[n-1] += S00[i] * Pns[ (n-2)*NN+i ];
   }

   p0[n-1] = 0.;
   if ( n==1 )
     p0[0] =  1. - CHI( ddf/s2*(cl-(1.-l)*hs)/l, df);
   else
     for (j=0; j<Ntilde; j++)
       p0[n-1] += Pns[ (n-1)*NN + ihs*Ntilde+j ] * Tn( (2.*hs-b[ihs+1]-b[ihs])/(b[ihs+1]-b[ihs]), j);
     
   mn_minus = 1.; mn_plus = 0.;
   if ( n > 1) {
     for (i=0; i<M; i++)
       for (j=0; j<Ntilde; j++) {
         oben = 0.;
         unten = 0.;
         for (jj=0; jj<Ntilde; jj++) {
           oben += Pns[ (n-1)*NN + i*Ntilde+jj ] * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
           unten+= Pns[ (n-2)*NN + i*Ntilde+jj ] * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
         }
         if ( fabs(unten)<1e-16 )
           if ( fabs(oben)<1e-16 ) q = 0.;
           else q = 1.;
         else q = oben/unten;
         if ( q<mn_minus ) mn_minus = q;
         if ( q>mn_plus ) mn_plus = q;
       }
     *rho = (mn_minus + mn_plus)/2.;
     if ( fabs(mn_plus - mn_minus) < FINALeps ) {       
       *nstop = n;
       n = nmax + 1;
     }     
   } /* n > 1 */
 } /* n=1; n<=nmax; n++ */ 

 Free(Pns);
 Free(zs);
 Free(ws);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);
 
 Free(S00);
 Free(p00);
 Free(VF0);

 return 0;
}


double seLR_sf_prerun_SIGMA_deluxe(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0)
{ double *ww, *zz, b1, b2, ddf2, *SF, rho, s2;
  int i, m, n, nstop, Nlocal;

 Nlocal = choose_N_for_se2(l, cl, cu);
  
 SF = vector(nmax);
 ww = vector(qm2);
 zz = vector(qm2);
 
 ddf2 = (double)(df2);  
 b1 = qCHI(     truncate/2., df2)/ddf2;
 b2 = qCHI(1. - truncate/2., df2)/ddf2; 
 gausslegendre(qm2, b1, b2, zz, ww); 
 for (i=0; i<qm2; i++) ww[i] *= ddf2 * chi( ddf2*zz[i], df2);
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm2; i++) {    
    s2 = zz[i];
    m = seLR_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop, &rho);
    if ( m != 0 ) warning("trouble with internal [package spc] function seLR_sf_deluxe");
    if ( nstop > 0 ) {
      for (n=0; n<nstop; n++)  p0[n] += ww[i] * SF[n];
      for (n=nstop; n<nmax; n++)  p0[n] += ww[i] * SF[nstop-1] * pow(rho, n-nstop+1);
    } else {
      for (n=0; n<nmax; n++)  p0[n] += ww[i] * SF[n];
    }
 }  
 
 Free(ww);
 Free(zz); 
 Free(SF);
 
 return 0;
}


double seLR_sf_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate, double *p0)
{ double *ww, *zz, b1, b2, ddf2, *SF, s2;
  int i, m, n, Nlocal;

 Nlocal = choose_N_for_se2(l, cl, cu);
  
 SF = vector(nmax);
 ww = vector(qm2);
 zz = vector(qm2);
 
 ddf2 = (double)(df2);  
 b1 = qCHI(     truncate/2., df2)/ddf2;
 b2 = qCHI(1. - truncate/2., df2)/ddf2; 
 gausslegendre(qm2, b1, b2, zz, ww); 
 for (i=0; i<qm2; i++) ww[i] *= ddf2 * chi( ddf2*zz[i], df2);
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm2; i++) {
   s2 = zz[i];
   m = seLR_sf(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF);
   if ( m != 0 ) warning("trouble with internal [package spc] function seLR_sf");
   for (n=0; n<nmax; n++)  p0[n] += ww[i] * SF[n];
 }  
 
 Free(ww);
 Free(zz); 
 Free(SF);
 
 return 0;
}


double seLR_Wq_prerun_SIGMA_deluxe(double l, double cl, double cu, double p, double hs, double sigma, int df1, int df2, int nmax, int qm1, int qm2, double truncate)
{ double *ww, *zz, b1, b2, ddf2, *SF, *p0, rho, *rhomany, *SFlast, Lp=-1., sf_level_adj, s2;
  int i, j, n, nstop, nstop_, nsm, nn, qnspecial=0, Nlocal;
  
 Nlocal = choose_N_for_se2(l, cl, cu); 
  
 p0 = vector(nmax);
 SF = vector(nmax); 
 rhomany = vector(qm2);
 SFlast = vector(qm2);
 ww = vector(qm2);
 zz = vector(qm2);
 
 ddf2 = (double)(df2);  
 b1 = qCHI(     truncate/2., df2)/ddf2;
 b2 = qCHI(1. - truncate/2., df2)/ddf2; 
 gausslegendre(qm2, b1, b2, zz, ww); 
 for (i=0; i<qm2; i++) ww[i] *= ddf2 * chi( ddf2*zz[i], df2);
 
 qnspecial = (qm2+1) / 2;
 
 s2 = zz[qnspecial];
 j = seLR_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nsm, &rho);
 if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
 n = nsm;
 
 if ( nsm < 1 ) { /* did not converge yet -- should be the rare case */
   nn = nmax;
   warning("The geometric tail approximation might not work.");
 } else {
   nstop = nsm;

   s2 = zz[qnspecial+1];
   j = seLR_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop_, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax ) {
       nstop = nstop_;
       i++;
       s2 = zz[qnspecial+i];
       j = seLR_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop_, &rho);
       if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }
 
   nstop = n;
   s2 = zz[qnspecial-1];
   j = seLR_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop_, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
   if ( nstop_ > nsm ) nsm = nstop_;
   if ( nstop_ < 1) nsm = nmax;
   if ( nstop_ >= nstop && nsm<nmax ) {
     i = 1;
     while ( nstop_ >= nstop && nsm<nmax ) {
       nstop = nstop_;
       i++;   
       s2 = zz[qnspecial-i];
       j = seLR_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nmax, qm1, SF, &nstop_, &rho);
       if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
       if ( nstop_ > nsm ) nsm = nstop_;
       if ( nstop_ < 1) nsm = nmax;
     }
   }  
   nn = nsm;
 }
 
 for (n=0; n<nmax; n++) p0[n] = 0.;
 
 for (i=0; i<qm2; i++) {
   s2 = zz[i];
   j = seLR_sf_deluxe(l, s2*cl, s2*cu, s2*hs, sigma, df1, Nlocal, nn, qm1, SF, &nstop, &rho);
   if ( j != 0 ) warning("trouble with internal [package spc] function seU_sf_deluxe");
   if ( nstop < 1 ) {
     nstop = nn;
     warning("The geometric tail approximation might not work.");     
   }
   rhomany[i] = rho;        
   for (n=0; n<nstop; n++)  p0[n] += ww[i] * SF[n];
   if ( nstop < nn) {
     for (n=nstop; n<nn; n++) p0[n] += ww[i] * SF[nstop-1] * pow(rho, n-nstop+1);          
   }
   SFlast[i] = SF[nstop-1] * pow(rho, nn-nstop);
 }
 
 sf_level_adj = 1.-p;
 if ( p0[nn-1] <= sf_level_adj ) { 
   n = nn-1;
   while (  p0[n] <= sf_level_adj &&  n > 0  ) n--;
   if ( p0[n] > 1.-p ) Lp = (double)( n + 2 ); else Lp = 1.;
 } else {
   for (n=nn; n<nmax; n++) {
      p0[n] = 0.;
      for (i=0; i<qm2; i++) p0[n] += ww[i] * SFlast[i] * pow(rhomany[i], n-nn+1);
      if (  p0[n] <= sf_level_adj  ) {
        Lp = (double)( n + 1 );
        n = nmax+1;
      }
   }
 }
 
 Free(p0);
 Free(ww);
 Free(zz); 
 Free(SF);
 Free(SFlast);
 Free(rhomany);

 return Lp;
}


double seLR_iglarl_prerun_SIGMA(double l, double cl, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate)
{ double *ww, *zz, b1, b2, result, ddf2, s2;
  int i;
  
  ww = vector(qm2);
  zz = vector(qm2);  
  ddf2 = (double)(df2);  
  b1 = qCHI(     truncate/2., df2)/ddf2;
  b2 = qCHI(1. - truncate/2., df2)/ddf2; 
  gausslegendre(qm2, b1, b2, zz, ww);  
  result = 0.;
  for (i=0; i<qm2; i++) {
    s2 = zz[i];
    result += ww[i] * ddf2 * chi( ddf2*s2, df2) * seLR_iglarl(l, s2*cl, s2*cu, s2*hs, sigma, df1, N, qm1);
  }
  Free(ww);
  Free(zz);
  
  return result;
}


double seLR_Wq(double l, double cl, double cu, double p, double hs, double sigma, int df, int N, int nmax, int qm)
{ double *S1s, *S2s, *Pns, *p0, *ws, *zs, *zch, *rside, *b, za=0., s2, ddf, xl, xu, q_minus=0., q_plus=0., dN, Hij, *S00, *p00, *VF0, mn_minus=1., mn_plus=0., oben, unten, q, enumerator=0., Wq=0.;
  int i, j, k, n, *ps, Ntilde, ihs, M, NN, ii, jj;

 M = ceil( (log(cl)-log(cu))/log(1.-l) );
 Ntilde = ceil( (double)N/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde;
 
 ihs = floor( (log(cl) - log(hs))/log(1.-l) );
 if ( ihs<0 ) ihs = 0;

 S1s = matrix(NN,NN);
 S2s = matrix(NN,NN);
 ps = ivector(NN);
 zch = matrix(M,Ntilde);
 rside = vector(NN+1);
 b   = vector(M+1);
 ws  = vector(qm);
 zs  = vector(qm);
 p0 = vector(nmax);
 Pns = matrix(nmax,NN);
 
 S00 = vector(NN);
 p00 = vector(nmax);
 VF0 = vector(NN+1);

/* interval borders b_i = cl/(1-l)^i */
 for (i=0; i<M; i++) b[i] = cl/pow(1.-l, (double)(i));
 b[M] = cu;

 /* Chebyshev nodes on [b_0,b_1],[b_1,b_2],...,[b_M-1,cu] */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     zch[ i*Ntilde+j ] = b[i] + (b[i+1]-b[i])/2.*(1.+cos(PI*(2.*j+1.)/2./dN));
   }

 /* P(L>1)(zch[i,j]) */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     rside[ i*Ntilde+j ] =  1. - CHI( ddf/s2*(cl-(1.-l)*zch[ i*Ntilde+j ])/l, df);
   }
 rside[NN] = 1. - CHI( ddf/s2*(cl-(1.-l)*cu)/l, df); /* reflexion at cu */
 
 /* P(zch[i,j] -> zreflect) */
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++)
     VF0[ i*Ntilde+j ] = 1. - CHI( ddf/s2*(cu-(1.-l)*zch[ i*Ntilde+j ])/l, df) ;
 VF0[NN] = 1. - CHI( ddf/s2*cu, df);
 
 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++) {
     za = (1.-l)*zch[ i*Ntilde+j ];
     for (ii=0; ii<M; ii++)
       for (jj=0; jj<Ntilde; jj++) {
         if ( b[ii+1]<za ) S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = 0.;
         else {
           if ( za<b[ii] ) xl = b[ii]; else xl = za;
           xu = b[ii+1];
           if ( df!=2 ) {
             xl = sqrt(xl-za);
             xu = sqrt(xu-za);
           }
           gausslegendre(qm, xl, xu, zs, ws);
           Hij = 0.;
           for (k=0; k<qm; k++)
             if ( df==2 )
               Hij += ws[k]*Tn( (2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj) * exp(-zs[k]/s2/l);
             else
               Hij += ws[k]*Tn( (2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                      * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
           if ( df==2 ) Hij *= exp(za/s2/l)/s2/l;
           else         Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf, ddf/2.);
           S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = Hij;
         }
       }
   }

 za = (1.-l)*cu;
 for (ii=0; ii<M; ii++)
   for (jj=0; jj<Ntilde; jj++) {
     if ( b[ii+1]<za ) S00[ ii*Ntilde+jj ] = 0.;
     else {
       if ( za<b[ii] ) xl = b[ii]; else xl = za;
       xu = b[ii+1];
       if ( df!=2 ) {
         xl = sqrt(xl-za);
         xu = sqrt(xu-za);
       }
       gausslegendre(qm, xl, xu, zs, ws);
       Hij = 0.;
       for (k=0; k<qm; k++)
         if ( df==2 )
           Hij += ws[k]*Tn((2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj) * exp(-zs[k]/s2/l);
         else
           Hij += ws[k] * Tn((2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                  * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/l);
       if ( df==2 ) Hij *= exp(za/s2/l)/s2/l;
       else         Hij /= gammafn(ddf/2.) * pow(2.*s2*l/ddf, ddf/2.);
       S00[ ii*Ntilde+jj ] = Hij;
     }
   }  
 
 for (i=0; i<NN; i++)
   for (j=0; j<NN; j++) S2s[i*NN+j] = 0.;

 for (i=0; i<M; i++)
   for (j=0; j<Ntilde; j++)
     for (jj=0; jj<Ntilde; jj++)
       S2s[ (i*Ntilde+j)*NN + i*Ntilde+jj ] = Tn( (2.*zch[ i*Ntilde+j ]-b[i+1]-b[i])/(b[i+1]-b[i]), jj); 

 LU_decompose(S2s,ps,NN);

 for (n=1; n<=nmax; n++) {
   if ( n==1) {
     for (i=0; i<M; i++)
       for (j=0; j<Ntilde; j++) {
         Pns[ i*Ntilde+j ] = 0.;
         for (jj=0; jj<Ntilde; jj++)
           Pns[ i*Ntilde+j ] += /*  usual Chebyshev approximation  */
             2./Ntilde * Tn( (2.*zch[i*Ntilde+jj]-b[i+1]-b[i])/(b[i+1]-b[i]), j) * rside[ i*Ntilde+jj ];
         if ( j==0 ) Pns[ i*Ntilde+j ] /= 2.;
       }
     p00[0] = rside[NN];
   }
   else {
     for (i=0; i<NN; i++) {
       rside[i] = VF0[i] * p00[n-2];
       for (j=0; j<NN; j++) rside[i] += S1s[ i*NN+j ] * Pns[ (n-2)*NN+j ];
     }
     LU_solve2(S2s, rside, ps, NN);     
     for (i=0; i<NN; i++) Pns[ (n-1)*NN+i ] = rside[i];    
     p00[n-1] = VF0[NN] * p00[n-2];     
     for (i=0 ;i<NN; i++) p00[n-1] += S00[i] * Pns[ (n-2)*NN+i ];
   }

   p0[n-1] = 0.;
   if ( n==1 )
     p0[0] =  1. - CHI( ddf/s2*(cl-(1.-l)*hs)/l, df);
   else
     for (j=0; j<Ntilde; j++)
       p0[n-1] += Pns[ (n-1)*NN + ihs*Ntilde+j ] * Tn( (2.*hs-b[ihs+1]-b[ihs])/(b[ihs+1]-b[ihs]), j);
     
   if ( p0[n-1] < 1.-p ) {
     Wq = (double)n;
     n = nmax+1;
   } else {  
     mn_minus = 1.; mn_plus = 0.;
     if ( n > 1) {
       for (i=0; i<M; i++)
         for (j=0; j<Ntilde; j++) {
           oben = 0.;
           unten = 0.;
           for (jj=0; jj<Ntilde; jj++) {
             oben += Pns[ (n-1)*NN + i*Ntilde+jj ] * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
             unten+= Pns[ (n-2)*NN + i*Ntilde+jj ] * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
           }
           if ( fabs(unten)<1e-16 )
             if ( fabs(oben)<1e-16 ) q = 0.;
             else q = 1.;
           else q = oben/unten;
           if ( q<mn_minus ) mn_minus = q;
           if ( q>mn_plus ) mn_plus = q;
         }
       enumerator = log( (1.-p)/p0[n-1] );
       q_minus = (double)n + enumerator/log(mn_minus);
       q_plus  = (double)n + enumerator/log(mn_plus); 
       /*if ( fabs( (q_plus-q_minus)/q_minus )<FINALeps ) n = nmax+1;*/     
       if ( fabs( ceil(q_plus) - ceil(q_minus) ) < .5 ) {
	 Wq = ceil(q_plus);
	 n = nmax +1;
       }
     } /* n > 1 */
   } /* p0[n-1] >= 1.-p */
 } /* n=1; n<=nmax; n++ */ 

 Free(Pns);
 Free(p0);
 Free(zs);
 Free(ws);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);
 
 Free(S00);
 Free(p00);
 Free(VF0);

 return Wq;
}


double seLR_crit(double l, double L0, double cu, double hs, double sigma, int df, int N, int qm)
{ double s1, s2, s3, ds, L1, L2, L3;

 s2 = hs;
 L2 = 0.;
 do {
   s1 = s2;
   L1 = L2;
   s2 *= 0.9;
   L2 = seLR_iglarl(l, s2, cu, hs, sigma, df, N, qm);
 } while ( L2<L0 );
 
 do {
   s3 = s1 + (L0-L1)/(L2-L1) * (s2-s1);
   L3 = seLR_iglarl(l, s3, cu, hs, sigma, df, N, qm);
   ds = s3-s2; s1 = s2; L1 = L2; s2 = s3; L2 = L3;   
 } while ( fabs(L0-L3)>1e-6 && fabs(ds)>1e-7 && s3>0.);

 return s3;
}


double stdeLR_crit(double l, double L0, double cu, double hs, double sigma, int df, int N, int qm)
{ double s1, s2, s3, ds, L1, L2, L3;

 s2 = hs;
 L2 = 0.;
 do {
   s1 = s2;
   L1 = L2;
   s2 -= .1;
   L2 = stdeLR_iglarl(l, s2, cu, hs, sigma, df, N, qm);
 } while ( L2<L0 && s2>0. );

 do {
   s3 = s1 + (L0-L1)/(L2-L1) * (s2-s1);
   L3 = stdeLR_iglarl(l, s3, cu, hs, sigma, df, N, qm);
   ds = s3-s2; s1 = s2; L1 = L2; s2 = s3; L2 = L3;
 } while ( fabs(L0-L3)>1e-7 && fabs(ds)>1e-8 && s3>0.);

 return s3;
}


double seLR_crit_prerun_SIGMA(double l, double L0, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate)
{ double s1, s2, s3, ds, L1=0., L2=0., L3=0.;

 s2 = hs;
 do {
   L1 = L2;
   /*s2 -= .1;*/
   s2 *= 0.9;
   L2 = seLR_iglarl_prerun_SIGMA(l, s2, cu, hs, sigma, df1, df2, N, qm1, qm2, truncate);
 } while ( L2 < L0 && s2 > 0. );

 s1 = s2 + .1;

 do {
   s3 = s1 + (L0-L1)/(L2-L1) * (s2-s1);
   L3 = seLR_iglarl_prerun_SIGMA(l, s3, cu, hs, sigma, df1, df2, N, qm1, qm2, truncate);
   ds = s3-s2; s1 = s2; L1 = L2; s2 = s3; L2 = L3;
 } while ( fabs(L0-L3)>1e-6 && fabs(ds)>1e-7 && s3>0.);

 return s3;
}


double seLR_q_crit(double l, int L0, double alpha, double cu, double hs, double sigma, int df, int N, int qm, double c_error, double a_error)
{ double s1, s2, s3, ds, p1, p2, p3, *SF;
  int result=1;

 SF  = vector(L0); 
  
 s2 = hs; p2 = 1.;
 do {
   p1 = p2;
   s2 -= .1;
   result = seLR_sf(l, s2, cu, hs, sigma, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in seLR_q_crit [package spc]");
   p2 = 1. - SF[L0-1];
 } while ( p2 > alpha && s2>0.);

 s1 = s2 + .1;

 do {
   s3 = s1 + (alpha - p1)/( p2 - p1 ) * (s2-s1);
   result = seLR_sf(l, s3, cu, hs, sigma, df, N, L0, qm, SF);
   if ( result != 0 ) warning("trouble in seLR_q_crit [package spc]");
   p3 = 1. - SF[L0-1];
   ds = s3 - s2; s1 = s2; p1 = p2; s2 = s3; p2 = p3;
 } while ( fabs(alpha - p3)>a_error && fabs(ds)>c_error );

 Free(SF);
 
 return s3;
}


double seLR_q_crit_prerun_SIGMA(double l, int L0, double alpha, double cu, double hs, double sigma, int df1, int df2, int N, int qm1, int qm2, double truncate, int tail_approx, double c_error, double a_error)
{ double s1, s2, s3, ds, p1, p2, p3, *SF;
  int result=1;  
  
 SF  = vector(L0); 
 
 s2 = seLR_q_crit(l, L0, alpha, cu, hs, sigma, df1, N, qm1, c_error, a_error); 
 if ( tail_approx ) result = seLR_sf_prerun_SIGMA_deluxe(l, s2, cu, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
 else result = seLR_sf_prerun_SIGMA(l, s2, cu, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
 if ( result != 0 ) warning("trouble in seLR_q_crit_prerun_SIGMA [package spc]");
 p2 = 1. - SF[L0-1];
 
 if ( p2 > alpha ) {
   do {
     p1 = p2;
     s2 -= .1;
     if ( tail_approx ) result = seLR_sf_prerun_SIGMA_deluxe(l, s2, cu, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     else result = seLR_sf_prerun_SIGMA(l, s2, cu, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     if ( result != 0 ) warning("trouble in seLR_q_crit_prerun_SIGMA [package spc]");
     p2 = 1. - SF[L0-1];
   } while ( p2 > alpha && s2 > 0. );
   s1 = s2 + .1;
 } else {
   do {
     p1 = p2;
     s2 += .1;
     if ( tail_approx ) result = seLR_sf_prerun_SIGMA_deluxe(l, s2, cu, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     else result = seLR_sf_prerun_SIGMA(l, s2, cu, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
     if ( result != 0 ) warning("trouble in seLR_q_crit_prerun_SIGMA [package spc]");
     p2 = 1. - SF[L0-1];
   } while ( p2 <= alpha && s2 < hs );
   s1 = s2 - .1;
 }

 do {
   s3 = s1 + (alpha - p1)/( p2 - p1 ) * (s2-s1);
   if ( tail_approx ) result = seLR_sf_prerun_SIGMA_deluxe(l, s3, cu, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
   else result = seLR_sf_prerun_SIGMA(l, s3, cu, hs, sigma, df1, df2, L0, qm1, qm2, truncate, SF);
   if ( result != 0 ) warning("trouble in seLR_q_crit_prerun_SIGMA [package spc]");
   p3 = 1. - SF[L0-1];
   ds = s3 - s2; s1 = s2; p1 = p2; s2 = s3; p2 = p3;
 } while ( fabs(alpha - p3)>a_error && fabs(ds)>c_error );

 Free(SF);
 
 return s3;
}


/* MEWMA: Rigdon (1995a,b) */


/* classical GL Nyström */
double mxewma_arl_0a(double lambda, double ce, int p, double hs, int N)
{ double *a, *g, *w, *z, arl, rr, r2;
  int i, j;

 a = matrix(N, N);
 g = vector(N);
 w = vector(N);
 z = vector(N);

 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 r2 = lambda*lambda;
 
 gausslegendre(N, 0., ce, z, w);

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*N+j] = -w[j] * nchi( z[j]/r2, p, rr*z[i] ) / r2;
   ++a[i*N+i];
 } 

 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);

 arl = 1.;
 for (j=0; j<N; j++) arl += w[j] * nchi(z[j]/r2, p, rr*hs)/r2 * g[j];

 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}


double mxewma_arl_f_0a(double lambda, double ce, int p, int N, double *g, double *w, double *z)
{ double *a, rr, r2;
  int i, j;

 a = matrix(N, N);

 ce *= lambda/(2.-lambda); 
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 r2 = lambda*lambda;
 
 gausslegendre(N, 0., ce, z, w);

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*N+j] = -w[j] * nchi( z[j]/r2, p, rr*z[i] ) / r2;
   ++a[i*N+i];
 } 

 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);

 Free(a);

 return 0.;
}


/* GL Nyström with changed arguments */
double mxewma_arl_0a2(double lambda, double ce, int p, double hs, int N)
{ double *a, *g, *w, *z, arl, rr, r2;
  int i, j;

 a = matrix(N, N);
 g = vector(N);
 w = vector(N);
 z = vector(N);

 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 r2 = lambda*lambda;
 
 gausslegendre(N, 0., sqrt(ce), z, w);

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*N+j] = -w[j] * nchi( z[j]*z[j]/r2, p, rr*z[i]*z[i] ) / r2 * 2.*z[j];
   ++a[i*N+i];
 } 

 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);

 arl = 1.;
 for (j=0; j<N; j++) arl += w[j] * nchi( z[j]*z[j]/r2, p, rr*hs)/r2 * g[j] * 2.*z[j];

 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}


double mxewma_arl_f_0a2(double lambda, double ce, int p, int N, double *g, double *w, double *z)
{ double *a, rr, r2;
  int i, j;

 a = matrix(N, N);

 ce *= lambda/(2.-lambda); 
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 r2 = lambda*lambda;
 
 gausslegendre(N, 0., sqrt(ce), z, w);

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*N+j] = -w[j] * nchi( z[j]*z[j]/r2, p, rr*z[i]*z[i] ) / r2 * 2.*z[j];
   ++a[i*N+i];
 } 

 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);

 Free(a);

 return 0.;
}

/* collocation */
double mxewma_arl_0b(double lambda, double ce, int p, double hs, int N, int qm)
{ double *a, *g, *w, *z, arl, rr, r2, xi, dN;
  int i, j, k;
  
 a = matrix(N, N);
 g = vector(N);
 w = vector(qm);
 z = vector(qm);
 
 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 r2 = lambda*lambda;
 dN = (double)N;
 
 gausslegendre(qm, 0, sqrt(ce), z, w);
 
  for (i=0; i<N; i++) {
   xi = ce/2. * ( 1. + cos(PI*(2.*(i+1.)-1.)/2./dN) );
   for (j=0; j<N; j++) {
     a[i*N+j] = Tn( (2.*xi-ce)/ce, j);
     for (k=0; k<qm; k++) a[i*N+j] -= w[k] * Tn( (2.*z[k]*z[k]-ce)/ce, j) * 2.*z[k] * nchi( z[k]*z[k]/r2, p, rr*xi ) / r2;
   }
 } 
 
 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);
 
 arl = 0.;
 for (j=0; j<N; j++) arl += g[j] * Tn( (2.*hs-ce)/ce ,j);
 
 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}


double mxewma_arl_f_0b(double lambda, double ce, int p, int N, int qm, double *g)
{ double *a, *w, *z, rr, r2, xi, dN;
  int i, j, k;
  
 a = matrix(N, N);
 w = vector(qm);
 z = vector(qm);
 
 ce *= lambda/(2.-lambda); 
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 r2 = lambda*lambda;
 dN = (double)N;
 
 gausslegendre(qm, 0, sqrt(ce), z, w);
 
  for (i=0; i<N; i++) {
   xi = ce/2. * ( 1. + cos(PI*(2.*(i+1.)-1.)/2./dN) );
   for (j=0; j<N; j++) {
     a[i*N+j] = Tn( (2.*xi-ce)/ce, j);
     for (k=0; k<qm; k++) a[i*N+j] -= w[k] * Tn( (2.*z[k]*z[k]-ce)/ce, j) * 2.*z[k] * nchi( z[k]*z[k]/r2, p, rr*xi ) / r2;
   }
 } 
 
 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);

 Free(a);
 Free(w);
 Free(z);

 return 0.;
}


/* Rigdon's approach -- Radau quadrature */
double mxewma_arl_0c(double lambda, double ce, int p, double hs, int N)
{ double *a, *g, *w, *z, arl, rr, r2;
  int i, j;

 a = matrix(N, N);
 g = vector(N);
 w = vector(N);
 z = vector(N);

 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 r2 = lambda*lambda;
 
 radau(N, 0., ce, z, w);

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*N+j] = -w[j] * nchi( z[j]/r2, p, rr*z[i] ) / r2;
   ++a[i*N+i];
 } 

 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);

 if ( hs > 1e-10 ) {
   arl = 1.;
   for (j=0; j<N; j++) arl += w[j] * nchi(z[j]/r2, p, rr*hs)/r2 * g[j];
 } else arl = g[0]; /* Rigdon's rationale behind the Radau quadrature */

 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}


double mxewma_arl_f_0c(double lambda, double ce, int p, int N, double *g, double *w, double *z)
{ double *a, rr, r2;
  int i, j;

 a = matrix(N, N);

 ce *= lambda/(2.-lambda); 
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 r2 = lambda*lambda;
 
 radau(N, 0., ce, z, w);

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*N+j] = -w[j] * nchi( z[j]/r2, p, rr*z[i] ) / r2;
   ++a[i*N+i];
 } 

 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);

 Free(a);

 return 0.;
}



/* Clenshaw–Curtis quadrature */ 
double mxewma_arl_0d(double lambda, double ce, int p, double hs, int N)
{ double *a, *g, *w, *z, arl, rr, r2, dN;
  int i, j;

 a = matrix(N, N);
 g = vector(N);
 w = vector(N);
 z = vector(N);
 
 dN = (double)N;
 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 r2 = lambda*lambda;
 
 /* nodes */
 for (i=0; i<N; i++) z[i] = ce * ( ( cos( i*PI/(dN-1.) ) + 1.)/2. );

 /* weights */
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*N+j] = cos( i*j*PI/(dN-1.) );
 } 
 for (j=0; j<N; j++) w[j] = iTn(1.,j) - iTn(-1,j);
 LU_solve(a, w, N);
 
 /* usual linear equation system */
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*N+j] = -w[j] * nchi( z[j]/r2, p, rr*z[i] ) / r2 * (ce/2.);
   ++a[i*N+i];
 } 
 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);
 
 arl = 1.;
 for (j=0; j<N; j++) arl += w[j] * nchi(z[j]/r2, p, rr*hs)/r2 * g[j] * (ce/2.);

 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}


double mxewma_arl_f_0d(double lambda, double ce, int p, int N, double *g, double *w, double *z)
{ double *a, rr, r2, dN;
  int i, j;

 a = matrix(N, N);
 
 dN = (double)N;
 ce *= lambda/(2.-lambda); 
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 r2 = lambda*lambda;
 
 /* nodes */
 for (i=0; i<N; i++) z[i] = ce * ( ( cos( i*PI/(dN-1.) ) + 1.)/2. );

 /* weights */
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*N+j] = cos( i*j*PI/(dN-1.) );
 } 
 for (j=0; j<N; j++) w[j] = iTn(1.,j) - iTn(-1,j);
 LU_solve(a, w, N);
 
 /* usual linear equation system */
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*N+j] = -w[j] * nchi( z[j]/r2, p, rr*z[i] ) / r2 * (ce/2.);
   ++a[i*N+i];
 } 
 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);
 
 Free(a);

 return 0.;
}



/* Markov chain (Runger and Prabhu 1996) */ 
double mxewma_arl_0e(double lambda, double ce, int p, double hs, int N)
{ double *a, *g, arl, rr, w, ncp, wl;
  int i, j;

 a = matrix(N, N);
 g = vector(N);

 ce = sqrt( ce * lambda/(2.-lambda) ); 
 hs = sqrt( hs * lambda/(2.-lambda) );
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 w = 2.*ce/(2.*N-1.);
 wl = w*w/(lambda*lambda);

 for (i=0; i<N; i++) {
   ncp = (w*i*i*w) * rr;
   a[i*N] = -nCHI( 0.25*wl, p, ncp );
   for (j=1; j<N; j++) a[i*N+j] = -( nCHI( (j+.5)*(j+.5)*wl, p, ncp ) - nCHI( (j-.5)*(j-.5)*wl, p, ncp ) );
   ++a[i*N+i];
 } 
 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);
 
 arl = g[ (int)floor(hs/w + .5) ];

 Free(a);
 Free(g);

 return arl;
}


double mxewma_arl_f_0e(double lambda, double ce, int p, int N, double *g, double *z)
{ double *a, rr, w, ncp, wl;
  int i, j;

 a = matrix(N, N);

 ce = sqrt( ce * lambda/(2.-lambda) ); 
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 w = 2.*ce/(2.*N-1.);
 wl = w*w/(lambda*lambda);

 for (i=0; i<N; i++) {
   ncp = (w*i*i*w)*rr;
   a[i*N] = -nCHI( 0.25*wl, p, ncp );
   for (j=1; j<N; j++) a[i*N+j] = -( nCHI( (j+.5)*(j+.5)*wl, p, ncp ) - nCHI( (j-.5)*(j-.5)*wl, p, ncp ) );
   ++a[i*N+i];
 } 
 for (j=0; j<N; j++) { g[j] = 1.; z[j] = w*(j+.5); }
 LU_solve(a, g, N);

 Free(a);

 return 0.;
}


double mxewma_psi0_e(double lambda, double ce, int p, int N, double *PSI)
{ double *a, rr, w, ncp, wl, rho=1., norm;
  int i, j, status, noofit;

 a = matrix(N, N);
    
 ce = sqrt( ce * lambda/(2.-lambda) ); 
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 w = 2.*ce/(2.*N-1.);
 wl = w*w/(lambda*lambda);

 for (i=0; i<N; i++) {
   ncp = (w*i*i*w) * rr;
   a[i] = nCHI( 0.25*wl, p, ncp );
   for (j=1; j<N; j++) a[j*N+i] = nCHI( (j+.5)*(j+.5)*wl, p, ncp ) - nCHI( (j-.5)*(j-.5)*wl, p, ncp );
 } 

 pmethod(N, a, &status, &rho, PSI, &noofit);
 
 norm = 0.;
 for (i=0; i<N; i++) norm += PSI[i];
 for (i=0; i<N; i++) PSI[i] /= norm;

 Free(a);

 return rho;
}


double mxewma_psiS0_e(double lambda, double ce, int p, int N, double *PSI)
{ double *a, *g, rr, w, ncp, wl, norm;
  int i, j;

 a = matrix(N, N);
 g = vector(N);

 ce = sqrt( ce * lambda/(2.-lambda) ); 
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 w = 2.*ce/(2.*N-1.);
 wl = w*w/(lambda*lambda);

 for (i=0; i<N; i++) {
   ncp = (w*i*i*w)*rr;
   a[i*N] = -nCHI( 0.25*wl, p, ncp );
   for (j=1; j<N; j++) a[i*N+j] = -( nCHI( (j+.5)*(j+.5)*wl, p, ncp ) - nCHI( (j-.5)*(j-.5)*wl, p, ncp ) );
   ++a[i*N+i];
 } 

 g[0] = 1.;
 for (j=1; j<N; j++) g[j] = 0.;
 solve(&N, a, g); 
 
 for (i=0; i<N; i++) PSI[i] = g[i];
 norm = 0.;
 for (i=0; i<N; i++) norm += PSI[i];
 for (i=0; i<N; i++) PSI[i] /= norm;
 
 Free(g);
 Free(a);

 return 1.;
}


/* Nyström with Simpson rule */
double mxewma_arl_0f(double lambda, double ce, int p, double hs, int N)
{ double *a, *g, *w, *z, arl, rr, r2, b;
  int i, j;

 a = matrix(N, N);
 g = vector(N);
 w = vector(N);
 z = vector(N);

 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 r2 = lambda*lambda;
 
 b = ce/((double)N-1.);
 for (i=0; i<N; i++) {
   z[i] = (double)i * b;
   if ( (i+1) % 2 == 0 )   w[i] = 4.;
   if ( (i+1) % 2 == 1 )   w[i] = 2.;
   if ( i==0 || i==(N-1) ) w[i] = 1.;
   w[i] *= b/3.;
 }

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*N+j] = -w[j] * nchi( z[j]/r2, p, rr*z[i] ) / r2;
   ++a[i*N+i];
 } 

 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);

 if ( hs > 1e-10 ) {
   arl = 1.;
   for (j=0; j<N; j++) arl += w[j] * nchi(z[j]/r2, p, rr*hs)/r2 * g[j];
 } else arl = g[0];

 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}


double mxewma_arl_f_0f(double lambda, double ce, int p, int N, double *g, double *w, double *z)
{ double *a, rr, r2, b;
  int i, j;

 a = matrix(N, N);

 ce *= lambda/(2.-lambda); 
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 r2 = lambda*lambda;
 
 b = ce/((double)N-1.);
 for (i=0; i<N; i++) {
   z[i] = (double)i * b;
   if ( (i+1) % 2 == 0 )   w[i] = 4.;
   if ( (i+1) % 2 == 1 )   w[i] = 2.;
   if ( i==0 || i==(N-1) ) w[i] = 1.;
   w[i] *= b/3.;
 }

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*N+j] = -w[j] * nchi( z[j]/r2, p, rr*z[i] ) / r2;
   ++a[i*N+i];
 } 

 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);

 Free(a);

 return 0.;
}


/* classical GL Nyström */
double mxewma_arl_1a(double lambda, double ce, int p, double delta, double hs, int N)
{ double rdc, r2, rr, a, b, *z0, *w0, *z1, *w1, *M, *g, term1, term2,arl, mean, sigma, eta, korr;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 g = vector(N2);  
 z0 = vector(N);
 w0 = vector(N);
 z1 = vector(N);
 w1 = vector(N);

 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);

 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 p1 = p - 1;
 
 gausslegendre(N,  0., 1., z0, w0);
 gausslegendre(N, -1., 1., z1, w1);     
 
 for (i=0; i<N; i++) {
   mean = rdc + (1.-lambda)*z1[i];
   for (j=0; j<N; j++) {
     eta = rr * ce * (1. - z1[i]*z1[i]) * z0[j];
     for (k=0; k<N; k++) {
       korr = ce * (1.-z1[k]*z1[k]) / r2;
       term1 = w1[k] * phi( ( z1[k]-mean)/sigma, 0.)/sigma * korr;
       for (l=0; l<N; l++) {
         term2 = w0[l] * nchi( korr*z0[l], p1, eta );
         /*M[i*N3 + j*N2 + k*N + l] = - term1 * term2;*/
         M[k*N3 + l*N2 + i*N + j] = - term1 * term2;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g); 
 
 b = 0.;
 a = 0.;
 mean = rdc + (1.-lambda)*b;
 eta = rr * ce*(1. - b*b)*a;
 arl = 1.;
 for (k=0; k<N; k++) {
    korr = ce*(1.-z1[k]*z1[k])/r2;
    term1 = w1[k] * phi( ( z1[k]-mean)/sigma, 0.)/sigma * korr;
    for (l=0; l<N; l++) {
      term2 = w0[l] * nchi( korr*z0[l], p1, eta );
      arl += term1 * term2 * g[k*N + l];
    } /* l = 0 .. N-1 */
  } /* k = 0 .. N-1 */
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(g);
 Free(M);
 
 return arl; 
}


double mxewma_arl_f_1a(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1)
{ double rdc, r2, rr, *M, term1, term2, mean, sigma, eta, korr;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);

 ce *= lambda/(2.-lambda); 
 
 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 p1 = p - 1;
 
 gausslegendre(N,  0., 1., z0, w0);
 gausslegendre(N, -1., 1., z1, w1);     
 
 for (i=0; i<N; i++) {
   mean = rdc + (1.-lambda)*z1[i];
   for (j=0; j<N; j++) {
     eta = rr * ce * (1. - z1[i]*z1[i]) * z0[j];
     for (k=0; k<N; k++) {
       korr = ce * (1.-z1[k]*z1[k]) / r2;
       term1 = w1[k] * phi( ( z1[k]-mean)/sigma, 0.)/sigma * korr;
       for (l=0; l<N; l++) {
         term2 = w0[l] * nchi( korr*z0[l], p1, eta );
         /*M[i*N3 + j*N2 + k*N + l] = - term1 * term2;*/
         M[k*N3 + l*N2 + i*N + j] = - term1 * term2;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g); 
 
 Free(M);
 
 return 0.; 
}


/* GL Nyström with inner argument changed */
double mxewma_arl_1a2(double lambda, double ce, int p, double delta, double hs, int N)
{ double rdc, r2, rr, a, b, *z0, *w0, *z1, *w1, *M, *g, term1, term2,arl, mean, sigma, eta, korr;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 g = vector(N2);  
 z0 = vector(N);
 w0 = vector(N);
 z1 = vector(N);
 w1 = vector(N);

 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);

 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 p1 = p - 1;
 
 gausslegendre(N,  0., 1., z0, w0);
 gausslegendre(N, -1., 1., z1, w1);     
 
 for (i=0; i<N; i++) {
   mean = rdc + (1.-lambda)*z1[i];
   for (j=0; j<N; j++) {
     eta = rr * ce * (1. - z1[i]*z1[i]) * z0[j]*z0[j];
     for (k=0; k<N; k++) {
       korr = ce * (1.-z1[k]*z1[k]) / r2;
       term1 = w1[k] * phi( ( z1[k]-mean)/sigma, 0.)/sigma * korr;
       for (l=0; l<N; l++) {
         term2 = w0[l] * nchi( korr*z0[l]*z0[l], p1, eta ) * 2.*z0[l];
         /*M[i*N3 + j*N2 + k*N + l] = - term1 * term2;*/
         M[k*N3 + l*N2 + i*N + j] = - term1 * term2;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 b = 0.;
 a = 0.;
 mean = rdc + (1.-lambda)*b;
 eta = rr * ce*(1. - b*b)*a*a;
 arl = 1.;
 for (k=0; k<N; k++) {
    korr = ce*(1.-z1[k]*z1[k])/r2;
    term1 = w1[k] * phi( ( z1[k]-mean)/sigma, 0.)/sigma * korr;
    for (l=0; l<N; l++) {
      term2 = w0[l] * nchi( korr*z0[l]*z0[l], p1, eta ) * 2.*z0[l];
      arl += term1 * term2 * g[k*N + l];
    } /* l = 0 .. N-1 */
  } /* k = 0 .. N-1 */
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(g);
 Free(M);
 
 return arl; 
}


double mxewma_arl_f_1a2(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1)
{ double rdc, r2, rr, *M, term1, term2, mean, sigma, eta, korr;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);

 ce *= lambda/(2.-lambda); 
 
 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 p1 = p - 1;
 
 gausslegendre(N,  0., 1., z0, w0);
 gausslegendre(N, -1., 1., z1, w1);     
 
 for (i=0; i<N; i++) {
   mean = rdc + (1.-lambda)*z1[i];
   for (j=0; j<N; j++) {
     eta = rr * ce * (1. - z1[i]*z1[i]) * z0[j]*z0[j];
     for (k=0; k<N; k++) {
       korr = ce * (1.-z1[k]*z1[k]) / r2;
       term1 = w1[k] * phi( ( z1[k]-mean)/sigma, 0.)/sigma * korr;
       for (l=0; l<N; l++) {
         term2 = w0[l] * nchi( korr*z0[l]*z0[l], p1, eta ) * 2.*z0[l];
         /*M[i*N3 + j*N2 + k*N + l] = - term1 * term2;*/
         M[k*N3 + l*N2 + i*N + j] = - term1 * term2;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 Free(M);
 
 return 0.; 
}


/* GL Nyström with both arguments changed -- sin */
double mxewma_arl_1a3(double lambda, double ce, int p, double delta, double hs, int N)
{ double rdc, r2, rr, a, b, *z0, *w0, *z1, *w1, *M, *g, term1, term2,arl, mean, sigma, eta, korr, vi, vk;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 g = vector(N2);  
 z0 = vector(N);
 w0 = vector(N);
 z1 = vector(N);
 w1 = vector(N);

 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);

 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 p1 = p - 1;
 
 gausslegendre(N,  0., 1., z0, w0);
 gausslegendre(N, -PI/2., PI/2., z1, w1);     
 
 for (i=0; i<N; i++) {
   vi = sin(z1[i]);
   mean = rdc + (1.-lambda)*vi;
   for (j=0; j<N; j++) {
     eta = rr * ce * (1. - vi*vi) * z0[j]*z0[j];
     for (k=0; k<N; k++) {
       vk = sin(z1[k]);
       korr = ce * (1.-vk*vk) / r2;
       term1 = w1[k] * phi( ( vk-mean)/sigma, 0.)/sigma * korr * cos(z1[k]);
       for (l=0; l<N; l++) {
         term2 = w0[l] * nchi( korr*z0[l]*z0[l], p1, eta ) * 2.*z0[l];
         /*M[i*N3 + j*N2 + k*N + l] = - term1 * term2;*/
         M[k*N3 + l*N2 + i*N + j] = - term1 * term2;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 b = 0.;
 a = 0.;
 mean = rdc + (1.-lambda)*sin(b);
 eta = rr * ce*(1. - b*b)*a*a;
 arl = 1.;
 for (k=0; k<N; k++) {
    vk = sin(z1[k]);   
    korr = ce*(1.-vk*vk)/r2;
    term1 = w1[k] * phi( ( vk-mean)/sigma, 0.)/sigma * korr * cos(z1[k]);
    for (l=0; l<N; l++) {
      term2 = w0[l] * nchi( korr*z0[l]*z0[l], p1, eta ) * 2.*z0[l];
      arl += term1 * term2 * g[k*N + l];
    } /* l = 0 .. N-1 */
  } /* k = 0 .. N-1 */
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(g);
 Free(M);
 
 return arl; 
}


double mxewma_arl_f_1a3(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1)
{ double rdc, r2, rr, *M, term1, term2, mean, sigma, eta, korr, vi, vk;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);

 ce *= lambda/(2.-lambda); 
 
 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 p1 = p - 1;
 
 gausslegendre(N,  0., 1., z0, w0);
 gausslegendre(N, -PI/2., PI/2., z1, w1);     
 
 for (i=0; i<N; i++) {
   vi = sin(z1[i]);
   mean = rdc + (1.-lambda)*vi;
   for (j=0; j<N; j++) {
     eta = rr * ce * (1. - vi*vi) * z0[j]*z0[j];
     for (k=0; k<N; k++) {
       vk = sin(z1[k]);
       korr = ce * (1.-vk*vk) / r2;
       term1 = w1[k] * phi( ( vk-mean)/sigma, 0.)/sigma * korr * cos(z1[k]);
       for (l=0; l<N; l++) {
         term2 = w0[l] * nchi( korr*z0[l]*z0[l], p1, eta ) * 2.*z0[l];
         /*M[i*N3 + j*N2 + k*N + l] = - term1 * term2;*/
         M[k*N3 + l*N2 + i*N + j] = - term1 * term2;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 Free(M);
 
 return 0.; 
}


/* GL Nyström with both arguments changed -- tan */
double mxewma_arl_1a4(double lambda, double ce, int p, double delta, double hs, int N)
{ double rdc, r2, rr, a, b, *z0, *w0, *z1, *w1, *M, *g, term1, term2,arl, mean, sigma, eta, korr, vi, vk;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 g = vector(N2);  
 z0 = vector(N);
 w0 = vector(N);
 z1 = vector(N);
 w1 = vector(N);

 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);

 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 p1 = p - 1;
 
 gausslegendre(N,  0., 1., z0, w0);
 gausslegendre(N, -PI/4., PI/4., z1, w1);     
 
 for (i=0; i<N; i++) {
   vi = tan(z1[i]);
   mean = rdc + (1.-lambda)*vi;
   for (j=0; j<N; j++) {
     eta = rr * ce * (1. - vi*vi) * z0[j]*z0[j];
     for (k=0; k<N; k++) {
       vk = tan(z1[k]);
       korr = ce * (1.-vk*vk) / r2;
       term1 = w1[k] * phi( ( vk-mean)/sigma, 0.)/sigma * korr / ( cos(z1[k])*cos(z1[k]) );
       for (l=0; l<N; l++) {
         term2 = w0[l] * nchi( korr*z0[l]*z0[l], p1, eta ) * 2.*z0[l];
         /*M[i*N3 + j*N2 + k*N + l] = - term1 * term2;*/
         M[k*N3 + l*N2 + i*N + j] = - term1 * term2;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 b = 0.;
 a = 0.;
 mean = rdc + (1.-lambda)*tan(b);
 eta = rr * ce*(1. - b*b)*a*a;
 arl = 1.;
 for (k=0; k<N; k++) {
    vk = tan(z1[k]);   
    korr = ce*(1.-vk*vk)/r2;
    term1 = w1[k] * phi( ( vk-mean)/sigma, 0.)/sigma * korr / ( cos(z1[k])*cos(z1[k]) );
    for (l=0; l<N; l++) {
      term2 = w0[l] * nchi( korr*z0[l]*z0[l], p1, eta ) * 2.*z0[l];
      arl += term1 * term2 * g[k*N + l];
    } /* l = 0 .. N-1 */
  } /* k = 0 .. N-1 */
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(g);
 Free(M);
 
 return arl; 
}


double mxewma_arl_f_1a4(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1)
{ double rdc, r2, rr, *M, term1, term2, mean, sigma, eta, korr, vi, vk;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);

 ce *= lambda/(2.-lambda); 
 
 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 p1 = p - 1;
 
 gausslegendre(N,  0., 1., z0, w0);
 gausslegendre(N, -PI/4., PI/4., z1, w1);     
 
 for (i=0; i<N; i++) {
   vi = tan(z1[i]);
   mean = rdc + (1.-lambda)*vi;
   for (j=0; j<N; j++) {
     eta = rr * ce * (1. - vi*vi) * z0[j]*z0[j];
     for (k=0; k<N; k++) {
       vk = tan(z1[k]);
       korr = ce * (1.-vk*vk) / r2;
       term1 = w1[k] * phi( ( vk-mean)/sigma, 0.)/sigma * korr / ( cos(z1[k])*cos(z1[k]) );
       for (l=0; l<N; l++) {
         term2 = w0[l] * nchi( korr*z0[l]*z0[l], p1, eta ) * 2.*z0[l];
         /*M[i*N3 + j*N2 + k*N + l] = - term1 * term2;*/
         M[k*N3 + l*N2 + i*N + j] = - term1 * term2;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 Free(M);
 
 return 0.; 
}


/* GL Nyström with both arguments changed -- sinh */
double mxewma_arl_1a5(double lambda, double ce, int p, double delta, double hs, int N)
{ double rdc, r2, rr, a, b, *z0, *w0, *z1, *w1, *M, *g, term1, term2,arl, mean, sigma, eta, korr, vi, vk, norm;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 g = vector(N2);  
 z0 = vector(N);
 w0 = vector(N);
 z1 = vector(N);
 w1 = vector(N);

 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);

 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 p1 = p - 1;
 
 gausslegendre(N,  0., 1., z0, w0);
 gausslegendre(N, -1., 1., z1, w1);
 norm = sinh(1.);
 
 for (i=0; i<N; i++) {
   vi = sinh(z1[i])/norm;
   mean = rdc + (1.-lambda)*vi;
   for (j=0; j<N; j++) {
     eta = rr * ce * (1. - vi*vi) * z0[j]*z0[j];
     for (k=0; k<N; k++) {
       vk = sinh(z1[k])/norm;
       korr = ce * (1.-vk*vk) / r2;
       term1 = w1[k] * phi( ( vk-mean)/sigma, 0.)/sigma * korr * cosh(z1[k])/norm;
       for (l=0; l<N; l++) {
         term2 = w0[l] * nchi( korr*z0[l]*z0[l], p1, eta ) * 2.*z0[l];
         /*M[i*N3 + j*N2 + k*N + l] = - term1 * term2;*/
         M[k*N3 + l*N2 + i*N + j] = - term1 * term2;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 b = 0.;
 a = 0.;
 mean = rdc + (1.-lambda)*sinh(b);
 eta = rr * ce*(1. - b*b)*a*a;
 arl = 1.;
 for (k=0; k<N; k++) {
    vk = sinh(z1[k])/norm;  
    korr = ce*(1.-vk*vk)/r2;
    term1 = w1[k] * phi( ( vk-mean)/sigma, 0.)/sigma * korr * cosh(z1[k])/norm;
    for (l=0; l<N; l++) {
      /*term2 = w0[l] * nchi( korr*z0[l]*z0[l], p1, eta ) * 2.*z0[l];*/
      term2 = w0[l] * chi( korr*z0[l]*z0[l], p1 ) * 2.*z0[l];
      arl += term1 * term2 * g[k*N + l];
    } /* l = 0 .. N-1 */
  } /* k = 0 .. N-1 */
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(g);
 Free(M);
 
 return arl; 
}


double mxewma_arl_f_1a5(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1)
{ double rdc, r2, rr, *M, term1, term2, mean, sigma, eta, korr, vi, vk, norm;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);

 ce *= lambda/(2.-lambda); 
 
 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 p1 = p - 1;
 
 gausslegendre(N,  0., 1., z0, w0);
 gausslegendre(N, -1., 1., z1, w1);
 norm = sinh(1.);
 
 for (i=0; i<N; i++) {
   vi = sinh(z1[i])/norm;
   mean = rdc + (1.-lambda)*vi;
   for (j=0; j<N; j++) {
     eta = rr * ce * (1. - vi*vi) * z0[j]*z0[j];
     for (k=0; k<N; k++) {
       vk = sinh(z1[k])/norm;
       korr = ce * (1.-vk*vk) / r2;
       term1 = w1[k] * phi( ( vk-mean)/sigma, 0.)/sigma * korr * cosh(z1[k])/norm;
       for (l=0; l<N; l++) {
         term2 = w0[l] * nchi( korr*z0[l]*z0[l], p1, eta ) * 2.*z0[l];
         /*M[i*N3 + j*N2 + k*N + l] = - term1 * term2;*/
         M[k*N3 + l*N2 + i*N + j] = - term1 * term2;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 Free(M);
 
 return 0.; 
}



/* GL Nyström with changed integration order */
/* arguments unchanged */
double mxewma_arl_1q(double lambda, double ce, int p, double delta, int N)
{ double *z0, *w0, *z1, *w1, *M, *g, term, arl, eta, mij, ncpij, lsd, korr, l2;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 g = vector(N2);  
 z0 = vector(N);
 w0 = vector(N);
 z1 = vector(N);
 w1 = vector(N);

 ce *= lambda/(2.-lambda); 

 lsd = lambda * sqrt(delta);
 eta = ((1.-lambda)/lambda) * ((1.-lambda)/lambda);
 l2 = lambda * lambda;
 p1 = p - 1;
 
 gausslegendre(N,  0., ce, z0, w0);
 gausslegendre(N, -1., 1., z1, w1);     
 
 for (i=0; i<N; i++) {   
   for (j=0; j<N; j++) {
     mij = lsd + (1.-lambda)*sqrt(z0[i])*z1[j];
     ncpij = eta * z0[i] * (1.-z1[j]*z1[j]);
     for (k=0; k<N; k++) {
       korr = w0[k] * sqrt(z0[k]) / l2;       
       for (l=0; l<N; l++) {
         term = korr * w1[l] * phi( (sqrt(z0[k])*z1[l] - mij)/lambda, 0.)/lambda * nchi( z0[k] * (1.-z1[l]*z1[l]) / l2, p1, ncpij );
         /*M[i*N3 + j*N2 + k*N + l] = - term;*/
         M[k*N3 + l*N2 + i*N + j] = - term;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 mij = lsd;
 arl = 1.;
 for (k=0; k<N; k++) {
    korr = w0[k] * sqrt(z0[k]) / l2;
    for (l=0; l<N; l++) {
      term = korr * w1[l] * phi( (sqrt(z0[k])*z1[l] - mij)/lambda, 0.)/lambda * chi( z0[k] * (1.-z1[l]*z1[l]) / l2, p1 );        
      arl += term * g[k*N + l];
    } /* l = 0 .. N-1 */
  } /* k = 0 .. N-1 */
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(g);
 Free(M);
 
 return arl; 
}


double mxewma_arl_f_1q(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1)
{ double *M, term, eta, mij, ncpij, lsd, korr, l2;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);

 ce *= lambda/(2.-lambda); 

 lsd = lambda * sqrt(delta);
 eta = ((1.-lambda)/lambda) * ((1.-lambda)/lambda);
 l2 = lambda * lambda;
 p1 = p - 1;
 
 gausslegendre(N,  0., ce, z0, w0);
 gausslegendre(N, -1., 1., z1, w1);     
 
 for (i=0; i<N; i++) {   
   for (j=0; j<N; j++) {
     mij = lsd + (1.-lambda)*sqrt(z0[i])*z1[j];
     ncpij = eta * z0[i] * (1.-z1[j]*z1[j]);
     for (k=0; k<N; k++) {
       korr = w0[k] * sqrt(z0[k]) / l2;       
       for (l=0; l<N; l++) {
         term = korr * w1[l] * phi( (sqrt(z0[k])*z1[l] - mij)/lambda, 0.)/lambda * nchi( z0[k] * (1.-z1[l]*z1[l]) / l2, p1, ncpij );
         /*M[i*N3 + j*N2 + k*N + l] = - term;*/
         M[k*N3 + l*N2 + i*N + j] = - term;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g); 

 Free(M);
 
 return 0.; 
}


/* ()^2 modification */
double mxewma_arl_1r(double lambda, double ce, int p, double delta, int N)
{ double *z0, *w0, *z1, *w1, *M, *g, term, arl, eta, mij, ncpij, lsd, korr, l2;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 g = vector(N2);  
 z0 = vector(N);
 w0 = vector(N);
 z1 = vector(N);
 w1 = vector(N);

 ce *= lambda/(2.-lambda); 
 ce = sqrt(ce);

 lsd = lambda * sqrt(delta);
 eta = ((1.-lambda)/lambda) * ((1.-lambda)/lambda);
 l2 = lambda * lambda;
 p1 = p - 1;
 
 gausslegendre(N,  0., ce, z0, w0);
 gausslegendre(N, -1., 1., z1, w1);     
 
 for (i=0; i<N; i++) {   
   for (j=0; j<N; j++) {
     mij = lsd + (1.-lambda)*z0[i]*z1[j];
     ncpij = eta * z0[i]*z0[i] * (1.-z1[j]*z1[j]);
     for (k=0; k<N; k++) {
       korr = 2. * w0[k] * z0[k]*z0[k] / l2;       
       for (l=0; l<N; l++) {
         term = korr * w1[l] * phi( (z0[k]*z1[l] - mij)/lambda, 0.)/lambda * nchi( z0[k]*z0[k] * (1.-z1[l]*z1[l]) / l2, p1, ncpij );
         /*M[i*N3 + j*N2 + k*N + l] = - term;*/
         M[k*N3 + l*N2 + i*N + j] = - term;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 mij = lsd;
 arl = 1.;
 for (k=0; k<N; k++) {
    korr = 2. * w0[k] * z0[k]*z0[k] / l2;
    for (l=0; l<N; l++) {
      term = korr * w1[l] * phi( (z0[k]*z1[l] - mij)/lambda, 0.)/lambda * chi( z0[k]*z0[k] * (1.-z1[l]*z1[l]) / l2, p1 );        
      arl += term * g[k*N + l];
    } /* l = 0 .. N-1 */
  } /* k = 0 .. N-1 */
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(g);
 Free(M);
 
 return arl; 
}


double mxewma_arl_f_1r(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1)
{ double *M, term, eta, mij, ncpij, lsd, korr, l2;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);

 ce *= lambda/(2.-lambda); 
 ce = sqrt(ce);

 lsd = lambda * sqrt(delta);
 eta = ((1.-lambda)/lambda) * ((1.-lambda)/lambda);
 l2 = lambda * lambda;
 p1 = p - 1;
 
 gausslegendre(N,  0., ce, z0, w0);
 gausslegendre(N, -1., 1., z1, w1);     
 
 for (i=0; i<N; i++) {   
   for (j=0; j<N; j++) {
     mij = lsd + (1.-lambda)*z0[i]*z1[j];
     ncpij = eta * z0[i]*z0[i] * (1.-z1[j]*z1[j]);
     for (k=0; k<N; k++) {
       korr = 2. * w0[k] * z0[k]*z0[k] / l2;       
       for (l=0; l<N; l++) {
         term = korr * w1[l] * phi( (z0[k]*z1[l] - mij)/lambda, 0.)/lambda * nchi( z0[k]*z0[k] * (1.-z1[l]*z1[l]) / l2, p1, ncpij );
         /*M[i*N3 + j*N2 + k*N + l] = - term;*/
         M[k*N3 + l*N2 + i*N + j] = - term;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 Free(M);
 
 return 0.; 
}


/* ... + sin() */
double mxewma_arl_1s(double lambda, double ce, int p, double delta, int N)
{ double *z0, *w0, *z1, *w1, *M, *g, term, arl, eta, mij, ncpij, lsd, korr, l2, wj, wl;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 g = vector(N2);  
 z0 = vector(N);
 w0 = vector(N);
 z1 = vector(N);
 w1 = vector(N);

 ce *= lambda/(2.-lambda); 
 ce = sqrt(ce);

 lsd = lambda * sqrt(delta);
 eta = ((1.-lambda)/lambda) * ((1.-lambda)/lambda);
 l2 = lambda * lambda;
 p1 = p - 1;
 
 gausslegendre(N,  0., ce, z0, w0);
 gausslegendre(N, -PI/2., PI/2., z1, w1);
 
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) {
     wj = sin( z1[j] );  
     mij = lsd + (1.-lambda)*z0[i]*wj;
     ncpij = eta * z0[i]*z0[i] * (1.-wj*wj);
     for (k=0; k<N; k++) {
       korr = 2. * w0[k] * z0[k]*z0[k] / l2;       
       for (l=0; l<N; l++) {
         wl = sin( z1[l] );
         term = korr * w1[l] * phi( (z0[k]*wl - mij)/lambda, 0.)/lambda * nchi( z0[k]*z0[k] * (1.-wl*wl) / l2, p1, ncpij ) *cos(z1[l]);
         /*M[i*N3 + j*N2 + k*N + l] = - term;*/
         M[k*N3 + l*N2 + i*N + j] = - term;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 mij = lsd;
 arl = 1.;
 for (k=0; k<N; k++) {
    korr = 2. * w0[k] * z0[k]*z0[k] / l2;
    for (l=0; l<N; l++) {
      wl = sin( z1[l] );
      term = korr * w1[l] * phi( (z0[k]*wl - mij)/lambda, 0.)/lambda * chi( z0[k]*z0[k] * (1.-wl*wl) / l2, p1 ) * cos(z1[l]);
      arl += term * g[k*N + l];
    } /* l = 0 .. N-1 */
  } /* k = 0 .. N-1 */
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(g);
 Free(M);
 
 return arl; 
}


double mxewma_arl_f_1s(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1)
{ double *M, term, eta, mij, ncpij, lsd, korr, l2, wj, wl;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);

 ce *= lambda/(2.-lambda); 
 ce = sqrt(ce);

 lsd = lambda * sqrt(delta);
 eta = ((1.-lambda)/lambda) * ((1.-lambda)/lambda);
 l2 = lambda * lambda;
 p1 = p - 1;
 
 gausslegendre(N,  0., ce, z0, w0);
 gausslegendre(N, -PI/2., PI/2., z1, w1);
 
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) {
     wj = sin( z1[j] );  
     mij = lsd + (1.-lambda)*z0[i]*wj;
     ncpij = eta * z0[i]*z0[i] * (1.-wj*wj);
     for (k=0; k<N; k++) {
       korr = 2. * w0[k] * z0[k]*z0[k] / l2;       
       for (l=0; l<N; l++) {
         wl = sin( z1[l] );
         term = korr * w1[l] * phi( (z0[k]*wl - mij)/lambda, 0.)/lambda * nchi( z0[k]*z0[k] * (1.-wl*wl) / l2, p1, ncpij ) *cos(z1[l]);
         /*M[i*N3 + j*N2 + k*N + l] = - term;*/
         M[k*N3 + l*N2 + i*N + j] = - term;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);

 Free(M);
 
 return 0.; 
}


/* ... + tan() */
double mxewma_arl_1t(double lambda, double ce, int p, double delta, int N)
{ double *z0, *w0, *z1, *w1, *M, *g, term, arl, eta, mij, ncpij, lsd, korr, l2, wj, wl;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 g = vector(N2);  
 z0 = vector(N);
 w0 = vector(N);
 z1 = vector(N);
 w1 = vector(N);

 ce *= lambda/(2.-lambda); 
 ce = sqrt(ce);

 lsd = lambda * sqrt(delta);
 eta = ((1.-lambda)/lambda) * ((1.-lambda)/lambda);
 l2 = lambda * lambda;
 p1 = p - 1;
 
 gausslegendre(N,  0., ce, z0, w0);
 gausslegendre(N, -PI/4., PI/4., z1, w1);
 
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) {
     wj = tan( z1[j] );  
     mij = lsd + (1.-lambda)*z0[i]*wj;
     ncpij = eta * z0[i]*z0[i] * (1.-wj*wj);
     for (k=0; k<N; k++) {
       korr = 2. * w0[k] * z0[k]*z0[k] / l2;       
       for (l=0; l<N; l++) {
         wl = tan( z1[l] );
         term = korr * w1[l] * phi( (z0[k]*wl - mij)/lambda, 0.)/lambda * nchi( z0[k]*z0[k] * (1.-wl*wl) / l2, p1, ncpij ) / ( cos(z1[l])*cos(z1[l]) );
         /*M[i*N3 + j*N2 + k*N + l] = - term;*/
         M[k*N3 + l*N2 + i*N + j] = - term;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
  /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 mij = lsd;
 arl = 1.;
 for (k=0; k<N; k++) {
    korr = 2. * w0[k] * z0[k]*z0[k] / l2;
    for (l=0; l<N; l++) {
      wl = tan( z1[l] );
      term = korr * w1[l] * phi( (z0[k]*wl - mij)/lambda, 0.)/lambda * chi( z0[k]*z0[k] * (1.-wl*wl) / l2, p1 ) / ( cos(z1[l])*cos(z1[l]) );
      arl += term * g[k*N + l];
    } /* l = 0 .. N-1 */
  } /* k = 0 .. N-1 */
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(g);
 Free(M);
 
 return arl; 
}


double mxewma_arl_f_1t(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1)
{ double *M, term, eta, mij, ncpij, lsd, korr, l2, wj, wl;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2); 

 ce *= lambda/(2.-lambda); 
 ce = sqrt(ce);

 lsd = lambda * sqrt(delta);
 eta = ((1.-lambda)/lambda) * ((1.-lambda)/lambda);
 l2 = lambda * lambda;
 p1 = p - 1;
 
 gausslegendre(N,  0., ce, z0, w0);
 gausslegendre(N, -PI/4., PI/4., z1, w1);
 
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) {
     wj = tan( z1[j] );  
     mij = lsd + (1.-lambda)*z0[i]*wj;
     ncpij = eta * z0[i]*z0[i] * (1.-wj*wj);
     for (k=0; k<N; k++) {
       korr = 2. * w0[k] * z0[k]*z0[k] / l2;       
       for (l=0; l<N; l++) {
         wl = tan( z1[l] );
         term = korr * w1[l] * phi( (z0[k]*wl - mij)/lambda, 0.)/lambda * nchi( z0[k]*z0[k] * (1.-wl*wl) / l2, p1, ncpij ) / ( cos(z1[l])*cos(z1[l]) );
         /*M[i*N3 + j*N2 + k*N + l] = - term;*/
         M[k*N3 + l*N2 + i*N + j] = - term;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g); 

 Free(M);
 
 return 0.; 
}


/* ... + sinh() */
double mxewma_arl_1u(double lambda, double ce, int p, double delta, int N)
{ double *z0, *w0, *z1, *w1, *M, *g, term, arl, eta, mij, ncpij, lsd, korr, l2, wj, wl, norm;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 g = vector(N2);  
 z0 = vector(N);
 w0 = vector(N);
 z1 = vector(N);
 w1 = vector(N);

 ce *= lambda/(2.-lambda); 
 ce = sqrt(ce);

 lsd = lambda * sqrt(delta);
 eta = ((1.-lambda)/lambda) * ((1.-lambda)/lambda);
 l2 = lambda * lambda;
 p1 = p - 1;
 
 gausslegendre(N,  0., ce, z0, w0); 
 gausslegendre(N, -1., 1., z1, w1);
 norm = sinh(1.);

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) {
     wj = sinh( z1[j] ) / norm;  
     mij = lsd + (1.-lambda)*z0[i]*wj;
     ncpij = eta * z0[i]*z0[i] * (1.-wj*wj);
     for (k=0; k<N; k++) {
       korr = 2. * w0[k] * z0[k]*z0[k] / l2;       
       for (l=0; l<N; l++) {
         wl = sinh( z1[l] ) / norm;
         term = korr * w1[l] * phi( (z0[k]*wl - mij)/lambda, 0.)/lambda * nchi( z0[k]*z0[k] * (1.-wl*wl) / l2, p1, ncpij ) * cosh(z1[l]) / norm;
         /*M[i*N3 + j*N2 + k*N + l] = - term;*/
         M[k*N3 + l*N2 + i*N + j] = - term;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 mij = lsd;
 arl = 1.;
 for (k=0; k<N; k++) {
    korr = 2. * w0[k] * z0[k]*z0[k] / l2;
    for (l=0; l<N; l++) {
      wl = sinh( z1[l] ) / norm;
      term = korr * w1[l] * phi( (z0[k]*wl - mij)/lambda, 0.)/lambda * chi( z0[k]*z0[k] * (1.-wl*wl) / l2, p1 ) * cosh(z1[l]) / norm;
      arl += term * g[k*N + l];
    } /* l = 0 .. N-1 */
  } /* k = 0 .. N-1 */
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(g);
 Free(M);
 
 return arl; 
}


double mxewma_arl_f_1u(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1)
{ double *M, term, eta, mij, ncpij, lsd, korr, l2, wj, wl, norm;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);

 ce *= lambda/(2.-lambda); 
 ce = sqrt(ce);

 lsd = lambda * sqrt(delta);
 eta = ((1.-lambda)/lambda) * ((1.-lambda)/lambda);
 l2 = lambda * lambda;
 p1 = p - 1;
 
 gausslegendre(N,  0., ce, z0, w0); 
 gausslegendre(N, -1., 1., z1, w1);
 norm = sinh(1.);

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) {
     wj = sinh( z1[j] ) / norm;  
     mij = lsd + (1.-lambda)*z0[i]*wj;
     ncpij = eta * z0[i]*z0[i] * (1.-wj*wj);
     for (k=0; k<N; k++) {
       korr = 2. * w0[k] * z0[k]*z0[k] / l2;       
       for (l=0; l<N; l++) {
         wl = sinh( z1[l] ) / norm;
         term = korr * w1[l] * phi( (z0[k]*wl - mij)/lambda, 0.)/lambda * nchi( z0[k]*z0[k] * (1.-wl*wl) / l2, p1, ncpij ) * cosh(z1[l]) / norm;
         /*M[i*N3 + j*N2 + k*N + l] = - term;*/
         M[k*N3 + l*N2 + i*N + j] = - term;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);

 Free(M);
 
 return 0.; 
}


double sign (double x) {
 double result;
 result = (double)(x > 1e-12) - (double)(x < -1e-12);
 return result;
}


/* collocation with two halfs in the same step + sin() */
double mxewma_arl_1b(double lambda, double ce, int p, double delta, double hs, int N, int qm0, int qm1)
{ double rdc, r2, rr, a, b, *z0, *w0, *z1, *w1, *M, *g, dN,
         term1, term2, term2a, term2b, innen, arl, mean, sigma, eta, u, u2, uu, v, v2, alpha;
  int r, s, i, j, k, l, N2, N3, p1;  
 
 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 g = vector(N2);  
 z0 = vector(qm0);
 w0 = vector(qm0);
 z1 = vector(qm1);
 w1 = vector(qm1);

 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);

 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 dN = (double)N;
 p1 = p - 1; 

 /* canonical Gauss-Legendre nodes and weights */
 gausslegendre(qm0,  0., 1., z0, w0);
 gausslegendre(qm1,  0., 1., z1, w1);
 
 for (s=0; s<N; s++) {   
   b = cos(PI*(2.*(s+1.)-1.)/2./dN); /* Chebyshev nodes */
   mean = rdc + (1.-lambda)*b;
   for (r=0; r<N; r++) {          
     a = 1/2. * ( 1. + cos(PI*(2.*(r+1.)-1.)/2./dN) ); /* Chebyshev nodes */
     eta = rr * ce*(1. - b*b)*a;
     for (i=0; i<N; i++) {
       for (j=0; j<N; j++) {
	 term1 = Tn( 2.*a-1., i) * Tn( b, j);
	 term2a = 0.;
	 term2b = 0.;
         for (k=0; k<qm1; k++) {
	   innen = 0.;
	   alpha = PI/2. * z1[k];
	   v = sin(alpha);
	   v2 = v*v;
	   if ( i==0 ) {
	     uu = ce * (1.-v2) / r2;
	     innen = nCHI(uu, p1, eta);
	   } else {
	     for (l=0; l<qm0; l++) {
	       u = z0[l];
	       u2 = u*u;
	       /*uu = ce * (1.-v2*v2) * u2 / r2;*/
	       uu = ce * (1. - v2) * u2/ r2;
	       innen += w0[l] * Tn( 2.*u2-1., i) * nchi( uu, p1, eta ) * 2.*u;
	     } /* l = 0 .. qm-1*/
	     innen *= ce * (1. - v2) / r2;
	   }
	   /*term2a += w1[k] * Tn(  v2, j) * phi( ( v2-mean)/sigma, 0.)/sigma * innen * ce * (1. - v2*v2) / r2 * 2.*v;
	   term2b += w1[k] * Tn( -v2, j) * phi( (-v2-mean)/sigma, 0.)/sigma * innen * ce * (1. - v2*v2) / r2 * 2.*v;*/
	   term2a += PI/2. * w1[k] * Tn(  v, j) * phi( ( v-mean)/sigma, 0.)/sigma * cos(alpha) * innen;
	   term2b += PI/2. * w1[k] * Tn( -v, j) * phi( (-v-mean)/sigma, 0.)/sigma * cos(alpha) * innen;
         } /* k = 0 .. qm-1*/
         term2 = term2a + term2b;         
         M[r*N3 + s*N2 + i*N + j] = term1 - term2;
       } /* j = 0 .. n-1 */
     } /* i = 0 .. n-1 */
   } /* r = 0 .. n-1 */
 } /* s = 0 .. n-1 */

 for (i=0; i<N2; i++) g[i] = 1.;
 LU_solve(M, g, N2);

 b = 0.;
 a = 0.;
 arl = 0.;
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) {
     arl += g[i*N + j] * Tn( 2.*a-1., i) * Tn( b, j);
   }
 }
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(g);
 Free(M);
 
 return arl; 
}


double mxewma_arl_f_1b(double lambda, double ce, int p, double delta, int N, int qm0, int qm1, double *g)
{ double rdc, r2, rr, a, b, *z0, *w0, *z1, *w1, *M, dN,
         term1, term2, term2a, term2b, innen, mean, sigma, eta, u, u2, uu, v, v2, alpha;
  int r, s, i, j, k, l, N2, N3, p1;  
 
 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 z0 = vector(qm0);
 w0 = vector(qm0);
 z1 = vector(qm1);
 w1 = vector(qm1);

 ce *= lambda/(2.-lambda); 

 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 dN = (double)N;
 p1 = p - 1; 

 /* canonical Gauss-Legendre nodes and weights */
 gausslegendre(qm0,  0., 1., z0, w0);
 gausslegendre(qm1,  0., 1., z1, w1);
 
 for (s=0; s<N; s++) {   
   b = cos(PI*(2.*(s+1.)-1.)/2./dN); /* Chebyshev nodes */
   mean = rdc + (1.-lambda)*b;
   for (r=0; r<N; r++) {          
     a = 1/2. * ( 1. + cos(PI*(2.*(r+1.)-1.)/2./dN) ); /* Chebyshev nodes */
     eta = rr * ce*(1. - b*b)*a;
     for (i=0; i<N; i++) {
       for (j=0; j<N; j++) {
	 term1 = Tn( 2.*a-1., i) * Tn( b, j);
	 term2a = 0.;
	 term2b = 0.;
         for (k=0; k<qm1; k++) {
	   innen = 0.;
	   alpha = PI/2. * z1[k];
	   v = sin(alpha);
	   v2 = v*v;
	   if ( i==0 ) {
	     uu = ce * (1.-v2) / r2;
	     innen = nCHI(uu, p1, eta);
	   } else {
	     for (l=0; l<qm0; l++) {
	       u = z0[l];
	       u2 = u*u;
	       /*uu = ce * (1.-v2*v2) * u2 / r2;*/
	       uu = ce * (1. - v2) * u2/ r2;
	       innen += w0[l] * Tn( 2.*u2-1., i) * nchi( uu, p1, eta ) * 2.*u;
	     } /* l = 0 .. qm-1*/
	     innen *= ce * (1. - v2) / r2;
	   }
	   /*term2a += w1[k] * Tn(  v2, j) * phi( ( v2-mean)/sigma, 0.)/sigma * innen * ce * (1. - v2*v2) / r2 * 2.*v;
	   term2b += w1[k] * Tn( -v2, j) * phi( (-v2-mean)/sigma, 0.)/sigma * innen * ce * (1. - v2*v2) / r2 * 2.*v;*/
	   term2a += PI/2. * w1[k] * Tn(  v, j) * phi( ( v-mean)/sigma, 0.)/sigma * cos(alpha) * innen;
	   term2b += PI/2. * w1[k] * Tn( -v, j) * phi( (-v-mean)/sigma, 0.)/sigma * cos(alpha) * innen;
         } /* k = 0 .. qm-1*/
         term2 = term2a + term2b;         
         M[r*N3 + s*N2 + i*N + j] = term1 - term2;
       } /* j = 0 .. n-1 */
     } /* i = 0 .. n-1 */
   } /* r = 0 .. n-1 */
 } /* s = 0 .. n-1 */

 for (i=0; i<N2; i++) g[i] = 1.;
 LU_solve(M, g, N2);

 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(M);
 
 return 0.; 
}


/* collocation with two halfs in the same step */
double mxewma_arl_1b3(double lambda, double ce, int p, double delta, double hs, int N, int qm0, int qm1)
{ double rdc, r2, rr, a, b, *z0, *w0, *z1, *w1, *M, *g, dN,
         term1, term2, term2a, term2b, innen, arl, mean, sigma, eta, u, u2, uu, v, v2, alpha;
  int r, s, i, j, k, l, N2, N3, p1;  
 
 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 g = vector(N2);  
 z0 = vector(qm0);
 w0 = vector(qm0);
 z1 = vector(qm1);
 w1 = vector(qm1);

 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);

 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 dN = (double)N;
 p1 = p - 1; 

 /* canonical Gauss-Legendre nodes and weights */
 gausslegendre(qm0,  0., 1., z0, w0);
 gausslegendre(qm1,  0., 1., z1, w1);
 
 for (s=0; s<N; s++) {   
   b = cos(PI*(2.*(s+1.)-1.)/2./dN); /* Chebyshev nodes */
   mean = rdc + (1.-lambda)*b;
   for (r=0; r<N; r++) {          
     a = 1/2. * ( 1. + cos(PI*(2.*(r+1.)-1.)/2./dN) ); /* Chebyshev nodes */
     eta = rr * ce*(1. - b*b)*a;
     for (i=0; i<N; i++) {
       for (j=0; j<N; j++) {
	 term1 = Tn( 2.*a-1., i) * Tn( b, j);
	 term2a = 0.;
	 term2b = 0.;
         for (k=0; k<qm1; k++) {
	   innen = 0.;
	   alpha = PI/4. * z1[k];
	   v = tan(alpha);
	   v2 = v*v;
	   if ( i==0 ) {
	     uu = ce * (1.-v2) / r2;
	     innen = nCHI(uu, p1, eta);
	   } else {
	     for (l=0; l<qm0; l++) {
	       u = z0[l];
	       u2 = u*u;
	       /*uu = ce * (1.-v2*v2) * u2 / r2;*/
	       uu = ce * (1. - v2) * u2/ r2;
	       innen += w0[l] * Tn( 2.*u2-1., i) * nchi( uu, p1, eta ) * 2.*u;
	     } /* l = 0 .. qm-1*/
	     innen *= ce * (1. - v2) / r2;
	   }
	   /*term2a += w1[k] * Tn(  v2, j) * phi( ( v2-mean)/sigma, 0.)/sigma * innen * ce * (1. - v2*v2) / r2 * 2.*v;
	   term2b += w1[k] * Tn( -v2, j) * phi( (-v2-mean)/sigma, 0.)/sigma * innen * ce * (1. - v2*v2) / r2 * 2.*v;*/
	   term2a += PI/4. * w1[k] * Tn(  v, j) * phi( ( v-mean)/sigma, 0.)/sigma / ( cos(alpha)*cos(alpha) ) * innen;
	   term2b += PI/4. * w1[k] * Tn( -v, j) * phi( (-v-mean)/sigma, 0.)/sigma / ( cos(alpha)*cos(alpha) ) * innen;
         } /* k = 0 .. qm-1*/
         term2 = term2a + term2b;         
         M[r*N3 + s*N2 + i*N + j] = term1 - term2;
       } /* j = 0 .. n-1 */
     } /* i = 0 .. n-1 */
   } /* r = 0 .. n-1 */
 } /* s = 0 .. n-1 */

 for (i=0; i<N2; i++) g[i] = 1.;
 LU_solve(M, g, N2);

 b = 0.;
 a = 0.;
 arl = 0.;
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) {
     arl += g[i*N + j] * Tn( 2.*a-1., i) * Tn( b, j);
   }
 }
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(g);
 Free(M);
 
 return arl; 
}


double mxewma_arl_f_1b3(double lambda, double ce, int p, double delta, int N, int qm0, int qm1, double *g)
{ double rdc, r2, rr, a, b, *z0, *w0, *z1, *w1, *M, dN,
         term1, term2, term2a, term2b, innen, mean, sigma, eta, u, u2, uu, v, v2, alpha;
  int r, s, i, j, k, l, N2, N3, p1;  
 
 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2); 
 z0 = vector(qm0);
 w0 = vector(qm0);
 z1 = vector(qm1);
 w1 = vector(qm1);

 ce *= lambda/(2.-lambda); 

 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 dN = (double)N;
 p1 = p - 1; 

 /* canonical Gauss-Legendre nodes and weights */
 gausslegendre(qm0,  0., 1., z0, w0);
 gausslegendre(qm1,  0., 1., z1, w1);
 
 for (s=0; s<N; s++) {   
   b = cos(PI*(2.*(s+1.)-1.)/2./dN); /* Chebyshev nodes */
   mean = rdc + (1.-lambda)*b;
   for (r=0; r<N; r++) {          
     a = 1/2. * ( 1. + cos(PI*(2.*(r+1.)-1.)/2./dN) ); /* Chebyshev nodes */
     eta = rr * ce*(1. - b*b)*a;
     for (i=0; i<N; i++) {
       for (j=0; j<N; j++) {
	 term1 = Tn( 2.*a-1., i) * Tn( b, j);
	 term2a = 0.;
	 term2b = 0.;
         for (k=0; k<qm1; k++) {
	   innen = 0.;
	   alpha = PI/4. * z1[k];
	   v = tan(alpha);
	   v2 = v*v;
	   if ( i==0 ) {
	     uu = ce * (1.-v2) / r2;
	     innen = nCHI(uu, p1, eta);
	   } else {
	     for (l=0; l<qm0; l++) {
	       u = z0[l];
	       u2 = u*u;
	       /*uu = ce * (1.-v2*v2) * u2 / r2;*/
	       uu = ce * (1. - v2) * u2/ r2;
	       innen += w0[l] * Tn( 2.*u2-1., i) * nchi( uu, p1, eta ) * 2.*u;
	     } /* l = 0 .. qm-1*/
	     innen *= ce * (1. - v2) / r2;
	   }
	   /*term2a += w1[k] * Tn(  v2, j) * phi( ( v2-mean)/sigma, 0.)/sigma * innen * ce * (1. - v2*v2) / r2 * 2.*v;
	   term2b += w1[k] * Tn( -v2, j) * phi( (-v2-mean)/sigma, 0.)/sigma * innen * ce * (1. - v2*v2) / r2 * 2.*v;*/
	   term2a += PI/4. * w1[k] * Tn(  v, j) * phi( ( v-mean)/sigma, 0.)/sigma / ( cos(alpha)*cos(alpha) ) * innen;
	   term2b += PI/4. * w1[k] * Tn( -v, j) * phi( (-v-mean)/sigma, 0.)/sigma / ( cos(alpha)*cos(alpha) ) * innen;
         } /* k = 0 .. qm-1*/
         term2 = term2a + term2b;         
         M[r*N3 + s*N2 + i*N + j] = term1 - term2;
       } /* j = 0 .. n-1 */
     } /* i = 0 .. n-1 */
   } /* r = 0 .. n-1 */
 } /* s = 0 .. n-1 */

 for (i=0; i<N2; i++) g[i] = 1.;
 LU_solve(M, g, N2);
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(M);
 
 return 0.; 
}


/* collocation with shrinked supports of the outer integral */
double mxewma_arl_1b2(double lambda, double ce, int p, double delta, double hs, int N, int qm0, int qm1)
{ double rdc, r2, rr, a, b, *z0, *w0, *z1, *w1, *M, *g, dN, lower, upper, xm, xw, alpha,
         term1, term2, innen, arl, mean, sigma, eta, u, u2, uu, v, v2;
  int r, s, i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 g = vector(N2);  
 z0 = vector(qm0);
 w0 = vector(qm0);
 z1 = vector(qm1);
 w1 = vector(qm1);

 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);

 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 dN = (double)N;
 p1 = p - 1;
 
 /* canonical Gauss-Legendre nodes and weights */
 gausslegendre(qm0,  0., 1., z0, w0);
 gausslegendre(qm1, -1., 1., z1, w1);
 
 for (s=0; s<N; s++) {   
   b = cos(PI*(2.*(s+1.)-1.)/2./dN); /* Chebyshev nodes */
   mean = rdc + (1.-lambda)*b;
   /* reasonable limits for the outer quadrature */
   lower = mean + sigma*qPHI(1e-9);
   if ( lower < -1. ) lower = -1.;
   upper = mean + sigma*qPHI(1.-1e-9);
   if ( upper >  1. ) upper = 1.;
   /* substitution sin(alpha) = v */
   lower = asin(lower);
   upper = asin(upper);
   /* constants for (-1,1) <-> (lower,upper) */
   xm = (lower+upper)/2.;
   xw = (upper-lower)/2.;
   for (r=0; r<N; r++) {          
     a = 1/2. * ( 1. + cos(PI*(2.*(r+1.)-1.)/2./dN) ); /* Chebyshev nodes */
     eta = rr * ce*(1. - b*b)*a;
     for (i=0; i<N; i++) {
       for (j=0; j<N; j++) {
	 term1 = Tn( 2.*a-1., i) * Tn( b, j);
	 term2 = 0.;
         for (k=0; k<qm1; k++) {
	   innen = 0.;
	   alpha = xm + xw*z1[k];
	   v = sin(alpha);
	   v2 = v*v;
	   if ( i==0 ) {
	     uu = ce * (1.-v2) / r2;
	     innen = nCHI(uu, p1, eta);
	   } else {
	     for (l=0; l<qm0; l++) {
	       u = z0[l];
	       u2 = u*u;
	       uu = ce * (1. - v2) * u2 / r2;
	       innen += w0[l] * Tn( 2.*u2-1., i) * nchi( uu, p1, eta ) * 2.*u;
	     }
	     innen *= ce * (1. - v2) / r2;
	   }
	   term2 += xw * w1[k] * Tn(  v, j) * phi( ( v-mean)/sigma, 0.)/sigma * cos(alpha) * innen;
         } /* k = 0 .. qm-1*/         
         M[r*N3 + s*N2 + i*N + j] = term1 - term2;
       } /* j = 0 .. n-1 */
     } /* i = 0 .. n-1 */
   } /* r = 0 .. n-1 */
 } /* s = 0 .. n-1 */

 for (i=0; i<N2; i++) g[i] = 1.;
 LU_solve(M, g, N2);

 b = 0.;
 a = 0.;
 arl = 0.;
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) {
     arl += g[i*N + j] * Tn( 2.*a-1., i) * Tn( b, j);
   }
 }
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(g);
 Free(M);
 
 return arl; 
}


double mxewma_arl_f_1b2(double lambda, double ce, int p, double delta, int N, int qm0, int qm1, double *g)
{ double rdc, r2, rr, a, b, *z0, *w0, *z1, *w1, *M, dN, lower, upper, xm, xw, alpha,
         term1, term2, innen, mean, sigma, eta, u, u2, uu, v, v2;
  int r, s, i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2); 
 z0 = vector(qm0);
 w0 = vector(qm0);
 z1 = vector(qm1);
 w1 = vector(qm1);

 ce *= lambda/(2.-lambda); 

 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 dN = (double)N;
 p1 = p - 1;
 
 /* canonical Gauss-Legendre nodes and weights */
 gausslegendre(qm0,  0., 1., z0, w0);
 gausslegendre(qm1, -1., 1., z1, w1);
 
 for (s=0; s<N; s++) {   
   b = cos(PI*(2.*(s+1.)-1.)/2./dN); /* Chebyshev nodes */
   mean = rdc + (1.-lambda)*b;
   /* reasonable limits for the outer quadrature */
   lower = mean + sigma*qPHI(1e-9);
   if ( lower < -1. ) lower = -1.;
   upper = mean + sigma*qPHI(1.-1e-9);
   if ( upper >  1. ) upper = 1.;
   /* substitution sin(alpha) = v */
   lower = asin(lower);
   upper = asin(upper);
   /* constants for (-1,1) <-> (lower,upper) */
   xm = (lower+upper)/2.;
   xw = (upper-lower)/2.;
   for (r=0; r<N; r++) {          
     a = 1/2. * ( 1. + cos(PI*(2.*(r+1.)-1.)/2./dN) ); /* Chebyshev nodes */
     eta = rr * ce*(1. - b*b)*a;
     for (i=0; i<N; i++) {
       for (j=0; j<N; j++) {
	 term1 = Tn( 2.*a-1., i) * Tn( b, j);
	 term2 = 0.;
         for (k=0; k<qm1; k++) {
	   innen = 0.;
	   alpha = xm + xw*z1[k];
	   v = sin(alpha);
	   v2 = v*v;
	   if ( i==0 ) {
	     uu = ce * (1.-v2) / r2;
	     innen = nCHI(uu, p1, eta);
	   } else {
	     for (l=0; l<qm0; l++) {
	       u = z0[l];
	       u2 = u*u;
	       uu = ce * (1. - v2) * u2 / r2;
	       innen += w0[l] * Tn( 2.*u2-1., i) * nchi( uu, p1, eta ) * 2.*u;
	     }
	     innen *= ce * (1. - v2) / r2;
	   }
	   term2 += xw * w1[k] * Tn(  v, j) * phi( ( v-mean)/sigma, 0.)/sigma * cos(alpha) * innen;
         } /* k = 0 .. qm-1*/         
         M[r*N3 + s*N2 + i*N + j] = term1 - term2;
       } /* j = 0 .. n-1 */
     } /* i = 0 .. n-1 */
   } /* r = 0 .. n-1 */
 } /* s = 0 .. n-1 */

 for (i=0; i<N2; i++) g[i] = 1.;
 LU_solve(M, g, N2);

 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(M);
 
 return 0.; 
}


/* collocation with two halfs in the same step + sinh() instead of sin() */
double mxewma_arl_1b4(double lambda, double ce, int p, double delta, double hs, int N, int qm0, int qm1)
{ double rdc, r2, rr, a, b, *z0, *w0, *z1, *w1, *M, *g, dN,
         term1, term2, term2a, term2b, innen, arl, mean, sigma, eta, u, u2, uu, v, v2, norm;
  int r, s, i, j, k, l, N2, N3, p1;  
 
 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 g = vector(N2);  
 z0 = vector(qm0);
 w0 = vector(qm0);
 z1 = vector(qm1);
 w1 = vector(qm1);

 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);

 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 dN = (double)N;
 p1 = p - 1; 

 /* canonical Gauss-Legendre nodes and weights */
 gausslegendre(qm0,  0., 1., z0, w0);
 gausslegendre(qm1,  0., 1., z1, w1);
 norm = sinh(1.);
 
 for (s=0; s<N; s++) {   
   b = cos(PI*(2.*(s+1.)-1.)/2./dN); /* Chebyshev nodes */
   mean = rdc + (1.-lambda)*b;
   for (r=0; r<N; r++) {          
     a = 1/2. * ( 1. + cos(PI*(2.*(r+1.)-1.)/2./dN) ); /* Chebyshev nodes */
     eta = rr * ce*(1. - b*b)*a;
     for (i=0; i<N; i++) {
       for (j=0; j<N; j++) {
	 term1 = Tn( 2.*a-1., i) * Tn( b, j);
	 term2a = 0.;
	 term2b = 0.;
         for (k=0; k<qm1; k++) {
	   innen = 0.;
	   v = sinh(z1[k])/norm;
	   v2 = v*v;
	   if ( i==0 ) {
	     uu = ce * (1.-v2) / r2;
	     innen = nCHI(uu, p1, eta);
	   } else {
	     for (l=0; l<qm0; l++) {
	       u = z0[l];
	       u2 = u*u;
	       /*uu = ce * (1.-v2*v2) * u2 / r2;*/
	       uu = ce * (1. - v2) * u2/ r2;
	       innen += w0[l] * Tn( 2.*u2-1., i) * nchi( uu, p1, eta ) * 2.*u;
	     } /* l = 0 .. qm-1*/
	     innen *= ce * (1. - v2) / r2;
	   }
	   /*term2a += w1[k] * Tn(  v2, j) * phi( ( v2-mean)/sigma, 0.)/sigma * innen * ce * (1. - v2*v2) / r2 * 2.*v;
	   term2b += w1[k] * Tn( -v2, j) * phi( (-v2-mean)/sigma, 0.)/sigma * innen * ce * (1. - v2*v2) / r2 * 2.*v;*/
	   term2a += w1[k] * Tn(  v, j) * phi( ( v-mean)/sigma, 0.)/sigma * cosh(z1[k])/norm * innen;
	   term2b += w1[k] * Tn( -v, j) * phi( (-v-mean)/sigma, 0.)/sigma * cosh(z1[k])/norm * innen;
         } /* k = 0 .. qm-1*/
         term2 = term2a + term2b;         
         M[r*N3 + s*N2 + i*N + j] = term1 - term2;
       } /* j = 0 .. n-1 */
     } /* i = 0 .. n-1 */
   } /* r = 0 .. n-1 */
 } /* s = 0 .. n-1 */

 for (i=0; i<N2; i++) g[i] = 1.;
 LU_solve(M, g, N2);

 b = 0.;
 a = 0.;
 arl = 0.;
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) {
     arl += g[i*N + j] * Tn( 2.*a-1., i) * Tn( b, j);
   }
 }
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(g);
 Free(M);
 
 return arl; 
}


double mxewma_arl_f_1b4(double lambda, double ce, int p, double delta, int N, int qm0, int qm1, double *g)
{ double rdc, r2, rr, a, b, *z0, *w0, *z1, *w1, *M, dN,
         term1, term2, term2a, term2b, innen, mean, sigma, eta, u, u2, uu, v, v2, norm;
  int r, s, i, j, k, l, N2, N3, p1;  
 
 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 z0 = vector(qm0);
 w0 = vector(qm0);
 z1 = vector(qm1);
 w1 = vector(qm1);

 ce *= lambda/(2.-lambda); 

 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 dN = (double)N;
 p1 = p - 1; 

 /* canonical Gauss-Legendre nodes and weights */
 gausslegendre(qm0,  0., 1., z0, w0);
 gausslegendre(qm1,  0., 1., z1, w1);
 norm = sinh(1.);
 
 for (s=0; s<N; s++) {   
   b = cos(PI*(2.*(s+1.)-1.)/2./dN); /* Chebyshev nodes */
   mean = rdc + (1.-lambda)*b;
   for (r=0; r<N; r++) {          
     a = 1/2. * ( 1. + cos(PI*(2.*(r+1.)-1.)/2./dN) ); /* Chebyshev nodes */
     eta = rr * ce*(1. - b*b)*a;
     for (i=0; i<N; i++) {
       for (j=0; j<N; j++) {
	 term1 = Tn( 2.*a-1., i) * Tn( b, j);
	 term2a = 0.;
	 term2b = 0.;
         for (k=0; k<qm1; k++) {
	   innen = 0.;
	   v = sinh(z1[k])/norm;
	   v2 = v*v;
	   if ( i==0 ) {
	     uu = ce * (1.-v2) / r2;
	     innen = nCHI(uu, p1, eta);
	   } else {
	     for (l=0; l<qm0; l++) {
	       u = z0[l];
	       u2 = u*u;
	       /*uu = ce * (1.-v2*v2) * u2 / r2;*/
	       uu = ce * (1. - v2) * u2/ r2;
	       innen += w0[l] * Tn( 2.*u2-1., i) * nchi( uu, p1, eta ) * 2.*u;
	     } /* l = 0 .. qm-1*/
	     innen *= ce * (1. - v2) / r2;
	   }
	   /*term2a += w1[k] * Tn(  v2, j) * phi( ( v2-mean)/sigma, 0.)/sigma * innen * ce * (1. - v2*v2) / r2 * 2.*v;
	   term2b += w1[k] * Tn( -v2, j) * phi( (-v2-mean)/sigma, 0.)/sigma * innen * ce * (1. - v2*v2) / r2 * 2.*v;*/
	   term2a += w1[k] * Tn(  v, j) * phi( ( v-mean)/sigma, 0.)/sigma * cosh(z1[k])/norm * innen;
	   term2b += w1[k] * Tn( -v, j) * phi( (-v-mean)/sigma, 0.)/sigma * cosh(z1[k])/norm * innen;
         } /* k = 0 .. qm-1*/
         term2 = term2a + term2b;         
         M[r*N3 + s*N2 + i*N + j] = term1 - term2;
       } /* j = 0 .. n-1 */
     } /* i = 0 .. n-1 */
   } /* r = 0 .. n-1 */
 } /* s = 0 .. n-1 */

 for (i=0; i<N2; i++) g[i] = 1.;
 LU_solve(M, g, N2);
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(M);
 
 return 0.; 
}


/* Radau/Gauß-Legendre Nyström -- Rigdon 1995b */
double mxewma_arl_1c(double lambda, double ce, int p, double delta, double hs, int N)
{ double rdc, r2, rr, a, b, *z0, *w0, *z1, *w1, *M, *g, term1, term2,arl, mean, sigma, eta, korr;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 g = vector(N2);  
 z0 = vector(N);
 w0 = vector(N);
 z1 = vector(N);
 w1 = vector(N);

 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);

 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 p1 = p - 1;
 
 radau(N, 0., 1., z0, w0);
 gausslegendre(N, -1., 1., z1, w1);
 
 for (i=0; i<N; i++) {
   mean = rdc + (1.-lambda)*z1[i];
   for (j=0; j<N; j++) {
     eta = rr * ce * (1. - z1[i]*z1[i]) * z0[j];
     for (k=0; k<N; k++) {
       korr = ce * (1.-z1[k]*z1[k]) / r2;
       term1 = w1[k] * phi( ( z1[k]-mean)/sigma, 0.)/sigma * korr;
       for (l=0; l<N; l++) {
         term2 = w0[l] * nchi( korr*z0[l], p1, eta );
         /*M[i*N3 + j*N2 + k*N + l] = - term1 * term2;*/
         M[k*N3 + l*N2 + i*N + j] = - term1 * term2;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 /* arl = g[ N*(N-1)/2 + 0 ]; */
 
 b = 0.;
 a = 0.;
 mean = rdc + (1.-lambda)*b;
 eta = rr * ce*(1. - b*b)*a;
 arl = 1.;
 for (k=0; k<N; k++) {
    korr = ce*(1.-z1[k]*z1[k])/r2;
    term1 = w1[k] * phi( ( z1[k]-mean)/sigma, 0.)/sigma * korr;
    for (l=0; l<N; l++) {
      term2 = w0[l] * nchi( korr*z0[l], p1, eta );
      arl += term1 * term2 * g[k*N + l];
    } /* l = 0 .. N-1 */
  } /* k = 0 .. N-1 */
  
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(g);
 Free(M);
 
 return arl; 
}   


double mxewma_arl_f_1c(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1)
{ double rdc, r2, rr, *M, term1, term2, mean, sigma, eta, korr;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);

 ce *= lambda/(2.-lambda); 

 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 p1 = p - 1;
 
 radau(N, 0., 1., z0, w0);
 gausslegendre(N, -1., 1., z1, w1);
 
 for (i=0; i<N; i++) {
   mean = rdc + (1.-lambda)*z1[i];
   for (j=0; j<N; j++) {
     eta = rr * ce * (1. - z1[i]*z1[i]) * z0[j];
     for (k=0; k<N; k++) {
       korr = ce * (1.-z1[k]*z1[k]) / r2;
       term1 = w1[k] * phi( ( z1[k]-mean)/sigma, 0.)/sigma * korr;
       for (l=0; l<N; l++) {
         term2 = w0[l] * nchi( korr*z0[l], p1, eta );
         /*M[i*N3 + j*N2 + k*N + l] = - term1 * term2;*/
         M[k*N3 + l*N2 + i*N + j] = - term1 * term2;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 Free(M);
 
 return 0.; 
}   


/* Clenshaw–Curtis quadrature */ 
double mxewma_arl_1d(double lambda, double ce, int p, double delta, double hs, int N)
{ double rdc, r2, rr, a, b, *z0, *w0, *z1, *w1, *M, *g, term1, term2,arl, mean, sigma, eta, korr, dN, *D;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 g = vector(N2);
 D = matrix(N, N);
 z0 = vector(N);
 w0 = vector(N);
 z1 = vector(N);
 w1 = vector(N);

 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);

 dN = (double)N;
 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 p1 = p - 1;
 
  /* nodes */
 for (i=0; i<N; i++) z0[i] = ( cos( i*PI/(dN-1.) ) + 1.)/2.;
 for (i=0; i<N; i++) z1[i] = cos( i*PI/(dN-1.) );
 /* weights */
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) D[i*N+j] = cos( i*j*PI/(dN-1.) );
 } 
 for (j=0; j<N; j++) w1[j] = iTn(1.,j) - iTn(-1,j);
 LU_solve(D, w1, N);
 for (j=0; j<N; j++) w0[j] = w1[j]/2.;
 
 for (i=0; i<N; i++) {
   mean = rdc + (1.-lambda)*z1[i];
   for (j=0; j<N; j++) {
     eta = rr * ce * (1. - z1[i]*z1[i]) * z0[j];
     for (k=0; k<N; k++) {
       korr = ce * (1.-z1[k]*z1[k]) / r2;
       term1 = w1[k] * phi( ( z1[k]-mean)/sigma, 0.)/sigma * korr;
       for (l=0; l<N; l++) {
         term2 = w0[l] * nchi( korr*z0[l], p1, eta );
         /*M[i*N3 + j*N2 + k*N + l] = - term1 * term2;*/
         M[k*N3 + l*N2 + i*N + j] = - term1 * term2;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 b = 0.;
 a = 0.;
 mean = rdc + (1.-lambda)*b;
 eta = rr * ce*(1. - b*b)*a;
 arl = 1.;
 for (k=0; k<N; k++) {
    korr = ce*(1.-z1[k]*z1[k])/r2;
    term1 = w1[k] * phi( ( z1[k]-mean)/sigma, 0.)/sigma * korr;
    for (l=0; l<N; l++) {
      term2 = w0[l] * nchi( korr*z0[l], p1, eta );
      arl += term1 * term2 * g[k*N + l];
    } /* l = 0 .. N-1 */
  } /* k = 0 .. N-1 */
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(D);
 Free(g);
 Free(M);
 
 return arl; 
}


double mxewma_arl_f_1d(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1)
{ double rdc, r2, rr, *M, term1, term2, mean, sigma, eta, korr, dN, *D;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 D = matrix(N, N);

 ce *= lambda/(2.-lambda); 

 dN = (double)N;
 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 p1 = p - 1;
 
  /* nodes */
 for (i=0; i<N; i++) z0[i] = ( cos( i*PI/(dN-1.) ) + 1.)/2.;
 for (i=0; i<N; i++) z1[i] = cos( i*PI/(dN-1.) );
 /* weights */
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) D[i*N+j] = cos( i*j*PI/(dN-1.) );
 } 
 for (j=0; j<N; j++) w1[j] = iTn(1.,j) - iTn(-1,j);
 LU_solve(D, w1, N);
 for (j=0; j<N; j++) w0[j] = w1[j]/2.;
 
 for (i=0; i<N; i++) {
   mean = rdc + (1.-lambda)*z1[i];
   for (j=0; j<N; j++) {
     eta = rr * ce * (1. - z1[i]*z1[i]) * z0[j];
     for (k=0; k<N; k++) {
       korr = ce * (1.-z1[k]*z1[k]) / r2;
       term1 = w1[k] * phi( ( z1[k]-mean)/sigma, 0.)/sigma * korr;
       for (l=0; l<N; l++) {
         term2 = w0[l] * nchi( korr*z0[l], p1, eta );
         /*M[i*N3 + j*N2 + k*N + l] = - term1 * term2;*/
         M[k*N3 + l*N2 + i*N + j] = - term1 * term2;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 Free(D);
 Free(M);
 
 return 0.; 
}


/* Markov chain (Runger and Prabhu 1996) and (Molnau, Runger, Montgomery, Skenner, Loredo, and Prabhu 2001) */ 
double mxewma_arl_1e(double lambda, double ce, int p, double delta, double hs, int N)
{ double *Q, *g, arl, rr, w, ncp, wl, dN, ce2, *V, *H, ci, z1, z2;
  int ix, jx, iy, jy, index, N2, X, Y, i, i0=0;

 dN = (double)N;
 ce = sqrt( ce * lambda/(2.-lambda) ); 
 w = 2.*ce/(2.*dN+1.);
  
 ce2 = ce*ce;
 N2 = 2*N+1;
 rr = (1.-lambda)/lambda * (1.-lambda)/lambda;
 wl = w*w/(lambda*lambda);
 
 index = 0;
 for (ix=0; ix<N2; ix++)
   for (iy=0; iy<N+1; iy++)
     index += ( (ix-dN)*(ix-dN) + iy*iy < ce2/(w*w) );
 
 V = matrix(N+1, N+1);
 for (iy=0; iy<N+1; iy++) {
   ncp = (w*iy*iy*w) * rr;
   V[iy*(N+1)] = nCHI( 0.25*wl, p-1, ncp );
   for (jy=1; jy<N+1; jy++) V[iy*(N+1)+jy] = nCHI( (jy+.5)*(jy+.5)*wl, p-1, ncp ) - nCHI( (jy-.5)*(jy-.5)*wl, p-1, ncp );
 }
 
 H = matrix(N2, N2);
 for (ix=0; ix<N2; ix++) {
   ci = -ce + (ix+.5)*w;
   for (jx=0; jx<N2; jx++) {
     z1 = ( -ce+(jx+1.)*w - (1.-lambda)*ci )/lambda - delta;
     z2 = ( -ce+     jx*w - (1.-lambda)*ci )/lambda - delta;
     H[ix*N2+jx] = PHI(z1, 0.) - PHI(z2, 0.);
   }
 }
 
 Q = matrix(index, index);
 g = vector(index);
 X = 0;
 for (ix=0; ix<N2; ix++) {
   for (iy=0; iy<N+1; iy++) {
     if ( (ix-dN)*(ix-dN) + iy*iy < ce2/(w*w) ) {
       X++;       
       if ( ix==N && iy==0 ) i0 = X-1;
       Y = 0;
       for (jx=0; jx<N2; jx++) {
         for (jy=0; jy<N+1; jy++) {
	   if ( (jx-dN)*(jx-dN) + jy*jy < ce2/(w*w) ) {
	     Y++;
	     Q[(X-1)*index + Y-1] = - H[ix*N2+jx] * V[iy*(N+1)+jy];
	     if ( X == Y ) ++Q[(X-1)*index + X-1];
	   }
         } /* l = 0 .. N-1 */
       } /* k = 0 .. N2-1 */
     }
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N2-1 */
 
 for (i=0; i<index; i++) g[i] = 1.;
 LU_solve(Q, g, index);

 arl = g[i0];

 Free(Q);
 Free(g);
 Free(V);
 Free(H);

 return arl;
}


double mxewma_arl_f_1e(double lambda, double ce, int p, double delta, int N, double *g, int *dQ)
{ double *Q, rr, w, ncp, wl, dN, ce2, *V, *H, ci, z1, z2;
  int ix, jx, iy, jy, index, N2, X, Y, i;

 dN = (double)N;
 ce = sqrt( ce * lambda/(2.-lambda) ); 
 w = 2.*ce/(2.*dN+1.);
  
 ce2 = ce*ce;
 N2 = 2*N+1;
 rr = (1.-lambda)/lambda * (1.-lambda)/lambda;
 wl = (w*w) / (lambda*lambda);
 
 index = 0;
 for (ix=0; ix<N2; ix++)
   for (iy=0; iy<N+1; iy++)
     index += ( (ix-dN)*(ix-dN) + iy*iy < ce2/(w*w) );
 *dQ = index;
   
 V = matrix(N+1, N+1);
 for (iy=0; iy<N+1; iy++) {
   ncp = (w*iy*iy*w) * rr;
   V[iy*(N+1)] = nCHI( 0.25*wl, p-1, ncp );
   for (jy=1; jy<N+1; jy++) V[iy*(N+1)+jy] = nCHI( (jy+.5)*(jy+.5)*wl, p-1, ncp ) - nCHI( (jy-.5)*(jy-.5)*wl, p-1, ncp );
 }
 
 H = matrix(N2, N2);
 for (ix=0; ix<N2; ix++) {
   ci = -ce + (ix+.5)*w;
   for (jx=0; jx<N2; jx++) {
     z1 = ( -ce+(jx+1.)*w - (1.-lambda)*ci )/lambda - delta;
     z2 = ( -ce+     jx*w - (1.-lambda)*ci )/lambda - delta;
     H[ix*N2+jx] = PHI(z1, 0.) - PHI(z2, 0.);
   }
 }
 
 Q = matrix(index, index);
 X = 0;
 for (ix=0; ix<N2; ix++) {
   for (iy=0; iy<N+1; iy++) {
     if ( (ix-dN)*(ix-dN) + iy*iy < ce2/(w*w) ) {
       X++;       
       /* if ( ix==N && iy==0 ) i0 = X; */
       Y = 0;
       for (jx=0; jx<N2; jx++) {
         for (jy=0; jy<N+1; jy++) {
	   if ( (jx-dN)*(jx-dN) + jy*jy < ce2/(w*w) ) {
	     Y++;
	     Q[(X-1)*index + Y-1] = - H[ix*N2+jx] * V[iy*(N+1)+jy];
	     if ( X == Y ) ++Q[(X-1)*index + X-1];
	   }
         } /* l = 0 .. N-1 */
       } /* k = 0 .. N2-1 */
     }
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N2-1 */
 
 for (i=0; i<index; i++) g[i] = 1.;
 LU_solve(Q, g, index);
 
 Free(Q);
 Free(V);
 Free(H);

 return 0.;
}


double mxewma_psi1_e(double lambda, double ce, int p, int N, double *PSI)
{ double *Q, rr, w, ncp, wl, dN, ce2, *V, *H, ci, z1, z2, rho, norm;
  int ix, jx, iy, jy, index, N2, X, Y, i, status, noofit;

 dN = (double)N;
 ce = sqrt( ce * lambda/(2.-lambda) ); 
 w = 2.*ce/(2.*dN+1.);
  
 ce2 = ce*ce;
 N2 = 2*N+1;
 rr = (1.-lambda)/lambda * (1.-lambda)/lambda;
 wl = w*w/(lambda*lambda);
 
 index = 0;
 for (ix=0; ix<N2; ix++)
   for (iy=0; iy<N+1; iy++)
     index += ( (ix-dN)*(ix-dN) + iy*iy < ce2/(w*w) );
 
 V = matrix(N+1, N+1);
 for (iy=0; iy<N+1; iy++) {
   ncp = (w*iy*iy*w) * rr;
   V[iy*(N+1)] = nCHI( 0.25*wl, p-1, ncp );
   for (jy=1; jy<N+1; jy++) V[iy*(N+1)+jy] = nCHI( (jy+.5)*(jy+.5)*wl, p-1, ncp ) - nCHI( (jy-.5)*(jy-.5)*wl, p-1, ncp );
 }
 
 H = matrix(N2, N2);
 for (ix=0; ix<N2; ix++) {
   ci = -ce + (ix+.5)*w;
   for (jx=0; jx<N2; jx++) {
     z1 = ( -ce+(jx+1.)*w - (1.-lambda)*ci )/lambda;
     z2 = ( -ce+     jx*w - (1.-lambda)*ci )/lambda;
     H[ix*N2+jx] = PHI(z1, 0.) - PHI(z2, 0.);
   }
 }
 
 Q = matrix(index, index);
 X = 0;
 for (ix=0; ix<N2; ix++) {
   for (iy=0; iy<N+1; iy++) {
     if ( (ix-dN)*(ix-dN) + iy*iy < ce2/(w*w) ) {
       X++;       
       Y = 0;
       for (jx=0; jx<N2; jx++) {
         for (jy=0; jy<N+1; jy++) {
	   if ( (jx-dN)*(jx-dN) + jy*jy < ce2/(w*w) ) {
	     Y++;
	     /*Q[(X-1)*index + Y-1] = - H[ix*N2+jx] * V[iy*(N+1)+jy];*/
         Q[(Y-1)*index + X-1] = H[ix*N2+jx] * V[iy*(N+1)+jy];
	     /* if ( X == Y ) ++Q[(X-1)*index + X-1]; */
	   }
         } /* l = 0 .. N-1 */
       } /* k = 0 .. N2-1 */
     }
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N2-1 */
 
 pmethod(index, Q, &status, &rho, PSI, &noofit);
 
 norm = 0.;
 for (i=0; i<index; i++) norm += PSI[i];
 for (i=0; i<index; i++) PSI[i] /= norm;  

 Free(Q);
 Free(V);
 Free(H);

 return rho;
}


double mxewma_psiS1_e(double lambda, double ce, int p, int N, double *PSI)
{ double *Q, *g, rr, w, ncp, wl, dN, ce2, *V, *H, ci, z1, z2, norm;
  int ix, jx, iy, jy, index, N2, X, Y, i, i0=0;

 dN = (double)N;
 ce = sqrt( ce * lambda/(2.-lambda) ); 
 w = 2.*ce/(2.*dN+1.);
  
 ce2 = ce*ce;
 N2 = 2*N+1;
 rr = (1.-lambda)/lambda * (1.-lambda)/lambda;
 wl = w*w/(lambda*lambda);
 
 index = 0;
 for (ix=0; ix<N2; ix++)
   for (iy=0; iy<N+1; iy++)
     index += ( (ix-dN)*(ix-dN) + iy*iy < ce2/(w*w) );
 
 V = matrix(N+1, N+1);
 for (iy=0; iy<N+1; iy++) {
   ncp = (w*iy*iy*w) * rr;
   V[iy*(N+1)] = nCHI( 0.25*wl, p-1, ncp );
   for (jy=1; jy<N+1; jy++) V[iy*(N+1)+jy] = nCHI( (jy+.5)*(jy+.5)*wl, p-1, ncp ) - nCHI( (jy-.5)*(jy-.5)*wl, p-1, ncp );
 }
 
 H = matrix(N2, N2);
 for (ix=0; ix<N2; ix++) {
   ci = -ce + (ix+.5)*w;
   for (jx=0; jx<N2; jx++) {
     z1 = ( -ce+(jx+1.)*w - (1.-lambda)*ci )/lambda;
     z2 = ( -ce+     jx*w - (1.-lambda)*ci )/lambda;
     H[ix*N2+jx] = PHI(z1, 0.) - PHI(z2, 0.);
   }
 }
 
 Q = matrix(index, index);
 g = vector(index);
 X = 0;
 for (ix=0; ix<N2; ix++) {
   for (iy=0; iy<N+1; iy++) {
     if ( (ix-dN)*(ix-dN) + iy*iy < ce2/(w*w) ) {
       X++;       
       if ( ix==N && iy==0 ) i0 = X-1;
       Y = 0;
       for (jx=0; jx<N2; jx++) {
         for (jy=0; jy<N+1; jy++) {
	   if ( (jx-dN)*(jx-dN) + jy*jy < ce2/(w*w) ) {
	     Y++;
	     Q[(X-1)*index + Y-1] = - H[ix*N2+jx] * V[iy*(N+1)+jy];
	     if ( X == Y ) ++Q[(X-1)*index + X-1];
	   }
         } /* l = 0 .. N-1 */
       } /* k = 0 .. N2-1 */
     }
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N2-1 */
 
 for (i=0; i<index; i++) g[i] = 0.;
 g[i0] = 1.;
 solve(&index, Q, g); 
 
 for (i=0; i<index; i++) PSI[i] = g[i];
 norm = 0.;
 for (i=0; i<index; i++) norm += PSI[i];
 for (i=0; i<index; i++) PSI[i] /= norm;
 
 Free(g);
 Free(Q);
 Free(V);
 Free(H);

 return 1.;
}


/* Simpson rule Nyström */
double mxewma_arl_1f(double lambda, double ce, int p, double delta, double hs, int N)
{ double rdc, r2, rr, *z0, *w0, *z1, *w1, *M, *g, term1, term2,arl, mean, sigma, eta, korr, b, a;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);
 g = vector(N2);  
 z0 = vector(N);
 w0 = vector(N);
 z1 = vector(N);
 w1 = vector(N);

 ce *= lambda/(2.-lambda); 
 hs *= lambda/(2.-lambda);

 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 p1 = p - 1;
 
 b = ce/((double)N-1.);
 for (i=0; i<N; i++) {
   z0[i] = (double)i * b;
   z1[i] = -1. + 2.*(double)i * b;
   if ( (i+1) % 2 == 0 )   w0[i] = 4.;
   if ( (i+1) % 2 == 1 )   w0[i] = 2.;
   if ( i==0 || i==(N-1) ) w0[i] = 1.;
   w0[i] *= b/3.;
   w1[i] = 2*w0[i];
 }
 
 for (i=0; i<N; i++) {
   mean = rdc + (1.-lambda)*z1[i];
   for (j=0; j<N; j++) {
     eta = rr * ce * (1. - z1[i]*z1[i]) * z0[j];
     for (k=0; k<N; k++) {
       korr = ce * (1.-z1[k]*z1[k]) / r2;
       term1 = w1[k] * phi( ( z1[k]-mean)/sigma, 0.)/sigma * korr;
       for (l=0; l<N; l++) {
         term2 = w0[l] * nchi( korr*z0[l], p1, eta );
         /*M[i*N3 + j*N2 + k*N + l] = - term1 * term2;*/
         M[k*N3 + l*N2 + i*N + j] = - term1 * term2;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 b = 0.;
 a = 0.;
 mean = rdc + (1.-lambda)*b;
 eta = rr * ce*(1. - b*b)*a;
 arl = 1.;
 for (k=0; k<N; k++) {
    korr = ce*(1.-z1[k]*z1[k])/r2;
    term1 = w1[k] * phi( ( z1[k]-mean)/sigma, 0.)/sigma * korr;
    for (l=0; l<N; l++) {
      term2 = w0[l] * nchi( korr*z0[l], p1, eta );
      arl += term1 * term2 * g[k*N + l];
    } /* l = 0 .. N-1 */
  } /* k = 0 .. N-1 */
 
 Free(w0);
 Free(z0);
 Free(w1);
 Free(z1);
 Free(g);
 Free(M);
 
 return arl; 
}  


double mxewma_arl_f_1f(double lambda, double ce, int p, double delta, int N, double *g, double *w0, double *z0, double *w1, double *z1)
{ double rdc, r2, rr, *M, term1, term2, mean, sigma, eta, korr, b;
  int i, j, k, l, N2, N3, p1;

 N2 = N*N;
 N3 = N2*N;
 
 M = matrix(N2, N2);

 ce *= lambda/(2.-lambda); 

 sigma = lambda/sqrt(ce);
 rdc = lambda*sqrt(delta/ce);
 r2 = lambda*lambda;
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda); 
 p1 = p - 1;
 
 b = ce/((double)N-1.);
 for (i=0; i<N; i++) {
   z0[i] = (double)i * b;
   z1[i] = -1. + 2.*(double)i * b;
   if ( (i+1) % 2 == 0 )   w0[i] = 4.;
   if ( (i+1) % 2 == 1 )   w0[i] = 2.;
   if ( i==0 || i==(N-1) ) w0[i] = 1.;
   w0[i] *= b/3.;
   w1[i] = 2*w0[i];
 }
 
 for (i=0; i<N; i++) {
   mean = rdc + (1.-lambda)*z1[i];
   for (j=0; j<N; j++) {
     eta = rr * ce * (1. - z1[i]*z1[i]) * z0[j];
     for (k=0; k<N; k++) {
       korr = ce * (1.-z1[k]*z1[k]) / r2;
       term1 = w1[k] * phi( ( z1[k]-mean)/sigma, 0.)/sigma * korr;
       for (l=0; l<N; l++) {
         term2 = w0[l] * nchi( korr*z0[l], p1, eta );
         /*M[i*N3 + j*N2 + k*N + l] = - term1 * term2;*/
         M[k*N3 + l*N2 + i*N + j] = - term1 * term2;
       } /* l = 0 .. N-1 */
     } /* k = 0 .. N-1 */
     ++M[i*N3 + j*N2 + i*N + j];
   } /* j = 0 .. N-1 */
 } /* i = 0 .. N-1 */

 for (j=0; j<N2; j++) g[j] = 1.;
 /*LU_solve(M, g, N2);*/
 solve(&N2, M, g);
 
 Free(M);
 
 return 0.; 
}  


double mxewma_crit(double lambda, double L0, int p, double hs, int N)
{ double c1, c2, c3, L1=0., L2=0., L3=0., dc;
  /*int numAlg;
  
  numAlg = remainder((double)p, 2.)==0;*/

  c2 = .5;
  L2 = 1.;
  do {
    c1 = c2;
    L1 = L2;
    c2 += 1.;
    /*if ( numAlg ) L2 = mxewma_arl_0a(lambda, c2, p, hs, N); else L2 = mxewma_arl_0b(lambda, c2, p, hs, N, qm);*/
    L2 = mxewma_arl_0a2(lambda, c2, p, hs, N);
  } while ( L2 < L0 );

  do {
    c3 = c1 + (L0-L1)/(L2-L1) * (c2-c1);
    /*if ( numAlg ) L3 = mxewma_arl_0a(lambda, c3, p, hs, N); else L3 = mxewma_arl_0b(lambda, c3, p, hs, N, qm);*/
    L3 = mxewma_arl_0a2(lambda, c3, p, hs, N);
    dc = c3 - c2; c1 = c2; L1 = L2; c2 = c3; L2 = L3;
  } while ( (fabs(L0-L3)>1e-8 ) && ( fabs(dc)>1e-10) );
  
  return c3;
} 


double mxewma_psi (double lambda, double ce, int p, int N, double *PSI, double *w, double *z)
{ double *a, rr, r2, rho, norm;
  int i, j, status, noofit;

 a = matrix(N, N); 

 ce *= lambda/(2.-lambda); 
 rr = ( (1.-lambda)/lambda ) * ( (1.-lambda)/lambda );
 r2 = lambda*lambda;
 
 gausslegendre(N, 0., sqrt(ce), z, w);

 for (i=0; i<N; i++)
   for (j=0; j<N; j++) a[i*N+j] = w[j] * nchi( z[i]*z[i]/r2, p, rr*z[j]*z[j] ) / r2 * 2.*z[j];
 
 pmethod(N, a, &status, &rho, PSI, &noofit);
 
 norm = 0.;
 for (i=0; i<N; i++) norm += w[i] * PSI[i] * 2.*z[i];
 for (i=0; i<N; i++) PSI[i] /= norm;  

 Free(a);

 return rho;
}


double mxewma_psiS(double lambda, double ce, int p, double hs, int N, double *PSI, double *w, double *z)
{ double *a, rr, r2, L0, *b;
  int i, j;

 if ( hs < 0. ) hs = 0.;  

 L0 = mxewma_arl_0a2(lambda, ce, p, hs, N);
 
 a = matrix(N, N);
 b = vector(N);
  
 ce *= lambda/(2.-lambda); 
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 r2 = lambda*lambda;
 
 gausslegendre(N, 0., sqrt(ce), z, w);
 
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[i*N+j] = -w[j] * nchi( z[i]*z[i]/r2, p, rr*z[j]*z[j] ) / r2 * 2.*z[j];
   ++a[i*N+i];
 }
 if ( hs < 1e-9 ) {
   for (i=0; i<N; i++) b[i] = chi( z[i]*z[i]/r2, p ) / r2 / L0;
 } else {
   for (i=0; i<N; i++) b[i] = nchi( z[i]*z[i]/r2, p, rr*hs*hs ) / r2 / L0;
 }
 LU_solve(a, b, N);

 for (i=0; i<N; i++) PSI[i] = b[i];
  
 Free(b);
 Free(a);
 
 return L0;
}


#define cond 0
#define cycl 1

#define GL 0
#define CO 1
#define RA 2
#define CC 3
#define MC 4
#define SR 5
#define CO2 6
#define GL2 7
#define GL3 8
#define GL4 9
#define GL5 10
#define CO3 11
#define CO4 12

#define nGL1 13
#define nGL2 14
#define nGL3 15
#define nGL4 16
#define nGL5 17


double mxewma_L_of_ab(double lambda, double ce, int p, double delta, int N, int qtype, double *g, double a, double b, double *w0, double *z0, double *w1, double *z1)
{ double LL=-1., ccee, rr, r2, rdc, sig, a_, b_, m, eta, korr, innen, norm;
  int i, j;
  
 ccee = ce * lambda/(2.-lambda); 
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 r2 = lambda*lambda;
  
  if ( fabs(delta)<1e-10 ) { /* in-control */
    if ( qtype==GL || qtype==RA || qtype==CC || qtype==SR ) {
      LL = 1.;
      for (j=0; j<N; j++) LL += w0[j] * g[j] * nchi( z0[j]/r2, p, rr*a ) / r2; 
      if ( qtype==CC ) LL *= ccee/2.;
    } 
    if ( qtype==GL2 ) {
      LL = 1.;
      for (j=0; j<N; j++) LL += w0[j] * g[j] * 2.*z0[j] * nchi( z0[j]*z0[j]/r2, p, rr*a ) / r2; 
    } 
    if ( qtype==CO ) {
      LL = 0.;
      for (j=0; j<N; j++) LL +=  Tn( (2.*a-ccee)/ccee, j) * g[j]; 
    }
    if ( qtype==MC ) {
      LL = 1. + g[0] * ( nCHI( z0[0]*z0[0]/r2, p, rr*a ) - 0. );
      for (j=1; j<N; j++) LL += g[j] * ( nCHI( z0[j]*z0[j]/r2, p, rr*a ) - nCHI( z0[j-1]*z0[j-1]/r2, p, rr*a ) ); 
    }
  } else { /* out-of-control */
    rdc = lambda * sqrt(delta/ccee);
    sig = lambda / sqrt(ccee);
    
    if ( fabs(ccee - a) < 1e-10 ) a_ = 1.; else a_ = (a - b*b/delta) / ( ccee - b*b/delta );
    b_ = b / sqrt( delta * ccee );
    m = rdc + (1.-lambda) * b_;
    eta = rr * ccee * ( 1. - b_*b_ ) * a_;
    if ( eta < 1e-10 ) eta = 0.;

    if ( qtype==GL || qtype==RA || qtype==CC || qtype==SR ) {
      LL = 1.;
      for (i=0; i<N; i++) {
        korr = ccee * ( 1. - z1[i]*z1[i] ) / r2;
        innen = 0.;
        for (j=0; j<N; j++) innen += w0[j] * nchi(korr*z0[j], p-1, eta) * g[i*N + j];
        LL += korr * w1[i] * phi( (z1[i]-m)/sig, 0. )/sig * innen;
      }
    } 
    if ( qtype==GL2 ) {
      LL = 1.;
      for (i=0; i<N; i++) {
        korr = ccee * ( 1. - z1[i]*z1[i] ) / r2;
        innen = 0.;
        for (j=0; j<N; j++) innen += w0[j] * nchi(korr*z0[j]*z0[j], p-1, eta) * g[i*N + j] * 2.*z0[j];
        LL += korr * w1[i] * phi( (z1[i]-m)/sig, 0. )/sig * innen;
      } 
    } 
    if ( qtype==GL3 ) {
      LL = 1.;
      for (i=0; i<N; i++) {
        korr = ccee * ( 1. - sin(z1[i])*sin(z1[i]) ) / r2;
        innen = 0.;
        for (j=0; j<N; j++) innen += w0[j] * nchi(korr*z0[j]*z0[j], p-1, eta) * g[i*N + j] * 2.*z0[j];
        LL += korr * w1[i] * phi( (sin(z1[i])-m)/sig, 0. )/sig * innen * cos(z1[i]);
      } 
    }
    if ( qtype==GL4 ) {
      LL = 1.;
      for (i=0; i<N; i++) {
        korr = ccee * ( 1. - tan(z1[i])*tan(z1[i]) ) / r2;
        innen = 0.;
        for (j=0; j<N; j++) innen += w0[j] * nchi(korr*z0[j]*z0[j], p-1, eta) * g[i*N + j] * 2.*z0[j];
        LL += korr * w1[i] * phi( (tan(z1[i])-m)/sig, 0. )/sig * innen / cos(z1[i])/cos(z1[i]);
      }
    } 
    if ( qtype==GL5 ) {
      norm = sinh(1.);
      LL = 1.;
      for (i=0; i<N; i++) {
        korr = ccee * ( 1. - sinh(z1[i])*sinh(z1[i])/norm/norm ) / r2;
        innen = 0.;
        for (j=0; j<N; j++) innen += w0[j] * nchi(korr*z0[j]*z0[j], p-1, eta) * g[i*N + j] * 2.*z0[j];
        LL += korr * w1[i] * phi( (sinh(z1[i])/norm-m)/sig, 0. )/sig * innen * cosh(z1[i])/norm;
      } 
    }    
    if ( qtype==CO || qtype==CO2 || qtype==CO3 || qtype==CO4 ) {
      LL = 0.;
      for (i=0; i<N; i++) {
        innen =0.;
        for (j=0; j<N; j++) innen += Tn(b_, j) * g[i*N + j]; 
        LL += Tn(2.*a_-1., i) * innen;
      }
    }
  }  
  
  return LL; 
}  


double mxewma_L_of_ag(double lambda, double ce, int p, double delta, int N, int qtype, double *g, double a, double gam, double *w0, double *z0, double *w1, double *z1)
{ double arl=-1., term, eta, mij, ncpij, lsd, korr, l2, wl, norm=1., korr2=1.;
  int k, l, p1;
  
  if ( qtype==nGL5 ) norm = sinh(1.);  
      
  lsd = lambda * sqrt(delta);
  eta = ((1.-lambda)/lambda) * ((1.-lambda)/lambda);
  l2 = lambda * lambda;
  p1 = p - 1;
  
  mij = lsd + (1.-lambda) * sqrt(a) * gam;
  ncpij = eta * a * (1.-gam*gam);
  
  if ( fabs(delta)<1e-10 ) { /* in-control */
    arl = -2.;
  } else {      
    arl = 1.;
    for (k=0; k<N; k++) {
      if ( qtype==nGL1 ) { korr = w0[k] * sqrt(z0[k]) / l2; }
      else { korr = 2. * w0[k] * z0[k]*z0[k] / l2; }
      for (l=0; l<N; l++) {
        wl = z1[l];
        korr2 = 1.;
        if ( qtype==nGL3 ) { wl = sin(z1[l]); korr2 = cos(z1[l]); }
        if ( qtype==nGL4 ) { wl = tan(z1[l]); korr2 = 1. / ( cos(z1[l])*cos(z1[l]) ); }
        if ( qtype==nGL5 ) { wl = sinh(z1[l]) / norm;  korr2 = cosh(z1[l]) / norm; }
        term = korr * w1[l] * phi( (sqrt(z0[k])*wl - mij)/lambda, 0.)/lambda * nchi( z0[k] * (1.-wl*wl) / l2, p1, ncpij ) * korr2;        
        arl += term * g[k*N + l];
      } /* l = 0 .. N-1 */
    } /* k = 0 .. N-1 */   
  }  
  
  return arl;
}    


double angle_unif_sphere(double x, int p)
{ double dp, result;
 dp = (double) p;
 if ( fabs(dp - 2.) < .001 ) result = 1./PI; else result = gammafn( dp/2. ) / gammafn( (dp-1.)/2. ) * pow(sin(x), dp - 2.) / sqrt(PI);
 return result;
}
  

double mxewma_ad(double lambda, double ce, int p, double delta, int N, int qm2, int psi_type, double hs, int qtype, int qm0, int qm1)
{ double *PSI, *ARL, *w1, *z1, *w2, *z2, *w3, *z3, *w4, *z4, *w5, *z5, zahl, ad, psi0, LL, ccee, rr, r2, xi, yj, sdelta;
  int i, j, N2;  
  
 PSI = vector(N);
 w1  = vector(N);
 z1  = vector(N);
 
 if ( hs < 0. ) hs = 0.;
 
 if ( psi_type == cond ) zahl = mxewma_psi (lambda, ce, p, N, PSI, w1, z1);
 if ( psi_type == cycl ) zahl = mxewma_psiS(lambda, ce, p, hs, N, PSI, w1, z1);
 
 ccee = ce * lambda/(2.-lambda); 
 rr = ((1.-lambda)/lambda)*((1.-lambda)/lambda);
 r2 = lambda*lambda;
 
 w3  = vector(qm2);
 z3  = vector(qm2);
 gausslegendre(qm2, 0., sqrt(ccee), z3, w3);

 w5  = vector(qm2);
 z5  = vector(qm2);
 gausslegendre(qm2, 0., PI, z5, w5);
 
 ad = 0.;
 
 if ( fabs(delta)<1e-10 ) { /* in-control */
   ARL = vector(N);
   w2  = vector(N);
   z2  = vector(N);
   
   if ( qtype == GL )  LL =  mxewma_arl_f_0a (lambda, ce, p, N, ARL, w2, z2);
   if ( qtype == GL2 ) LL =  mxewma_arl_f_0a2(lambda, ce, p, N, ARL, w2, z2);
   if ( qtype == CO )  LL =  mxewma_arl_f_0b (lambda, ce, p, N, qm0, ARL);   
   if ( qtype == RA )  LL =  mxewma_arl_f_0c (lambda, ce, p, N, ARL, w2, z2);
   if ( qtype == CC )  LL =  mxewma_arl_f_0d (lambda, ce, p, N, ARL, w2, z2);
   if ( qtype == MC )  LL =  mxewma_arl_f_0e (lambda, ce, p, N, ARL, z2);
   if ( qtype == SR )  LL =  mxewma_arl_f_0f (lambda, ce, p, N, ARL, w2, z2);   

   for (i=0; i<qm2; i++) {
     xi = z3[i]*z3[i];
     
     psi0 = 0.;
     if ( psi_type == cycl ) {
       if ( fabs(hs) <= 1e-10 ) psi0 = chi(xi/r2, p) / r2 / zahl;
       if ( fabs(hs) > 1e-10 )  psi0 = 2.*hs * nchi(xi/r2, p, rr*hs*hs) / r2 / zahl;
     }  
     for (j=0; j<N; j++) psi0 += w1[j] * PSI[j] * 2.*z1[j] * nchi( xi/r2, p, rr*z1[j]*z1[j] ) / r2;      
     if ( psi_type == cond ) psi0 /= zahl;
     
     LL = mxewma_L_of_ab(lambda, ce, p, 0., N, qtype, ARL, xi, 0., w2, z2, w2, z2);  
     
     /*printf("%2d\t%.4f\t%.4f\t%.4f\n", i, z3[i], psi0, LL);*/
     
     ad += w3[i] * 2.*z3[i] * psi0 * LL;
   }  
   
   if ( psi_type == cycl ) {
     psi0 = 1./zahl;     
     LL = mxewma_L_of_ab(lambda, ce, p, 0., N, qtype, ARL, hs, 0., w2, z2, w2, z2);
     ad += psi0 * LL;
   }
   
   Free(z2);
   Free(w2);
   Free(ARL);
   
 } else { /* out-of-control */
   sdelta = sqrt(delta);
   
   N2  = N*N;
   ARL = vector(N2);
   w2  = vector(N);   
   z2  = vector(N);
   w4  = vector(N);
   z4  = vector(N);
   
   if ( qtype == GL )  LL = mxewma_arl_f_1a (lambda, ce, p, delta, N, ARL, w2, z2, w4, z4);
   if ( qtype == GL2 ) LL = mxewma_arl_f_1a2(lambda, ce, p, delta, N, ARL, w2, z2, w4, z4);
   if ( qtype == GL3 ) LL = mxewma_arl_f_1a3(lambda, ce, p, delta, N, ARL, w2, z2, w4, z4);   
   if ( qtype == GL4 ) LL = mxewma_arl_f_1a4(lambda, ce, p, delta, N, ARL, w2, z2, w4, z4);
   if ( qtype == GL5 ) LL = mxewma_arl_f_1a5(lambda, ce, p, delta, N, ARL, w2, z2, w4, z4);
   
   if ( qtype == CO )  LL = mxewma_arl_f_1b (lambda, ce, p, delta, N, qm0, qm1, ARL);
   if ( qtype == CO2 ) LL = mxewma_arl_f_1b2(lambda, ce, p, delta, N, qm0, qm1, ARL);
   if ( qtype == CO3 ) LL = mxewma_arl_f_1b3(lambda, ce, p, delta, N, qm0, qm1, ARL);
   if ( qtype == CO4 ) LL = mxewma_arl_f_1b4(lambda, ce, p, delta, N, qm0, qm1, ARL);
   
   if ( qtype == RA )  LL = mxewma_arl_f_1c(lambda, ce, p, delta, N, ARL, w2, z2, w4, z4);
   if ( qtype == CC )  LL = mxewma_arl_f_1d(lambda, ce, p, delta, N, ARL, w2, z2, w2, z2);
   /*if ( qtype == MC )  LL = mxewma_arl_f_1e(lambda, ce, p, delta, N, ARL, &dQ);*/
   if ( qtype == SR )  LL = mxewma_arl_f_1f(lambda, ce, p, delta, N, ARL, w2, z2, w2, z2);
       
   for (i=0; i<qm2; i++) {
     xi = z3[i]*z3[i];
     
     psi0 = 0.;
     if ( psi_type == cycl ) {
       if ( fabs(hs) <= 1e-10 ) psi0 = chi(xi/r2, p) / r2 / zahl;
       if ( fabs(hs) > 1e-10 )  psi0 = 2.*hs * nchi(xi/r2, p, rr*hs*hs) / r2 / zahl;
     }  
     for (j=0; j<N; j++) psi0 += w1[j] * PSI[j] * 2.*z1[j] * nchi( xi/r2, p, rr*z1[j]*z1[j] ) / r2;      
     if ( psi_type == cond ) psi0 /= zahl;
     
     for (j=0; j<qm2; j++) {
       yj = z3[i] * sdelta * cos(z5[j]);
 
       LL = mxewma_L_of_ab(lambda, ce, p, delta, N, qtype, ARL, xi, yj, w2, z2, w4, z4); 
 
       ad += w3[i] * 2.*z3[i] * w5[j] * psi0 * angle_unif_sphere(z5[j], p) * LL;
     }
   }
   
   if ( psi_type == cycl ) {
     psi0 = 1./zahl;     
     LL = mxewma_L_of_ab(lambda, ce, p, delta, N, qtype, ARL, 0., 0., w2, z2, w4, z4);
     ad += psi0 * LL;
   }
   
   Free(z4);
   Free(w4);
   Free(z2);
   Free(z2);
   Free(ARL);
 } 
 
 Free(z3);
 Free(w3); 
 Free(z1);
 Free(w1);
 Free(PSI);
 
 return ad;
}  


double cos_unif_sphere(double g, int p)
{ double dp, result;
 dp = (double) p;
 if ( fabs(dp - 3.) < .001 ) result = 0.5; else result = gammafn( dp/2. ) / gammafn( (dp-1.)/2. ) * pow( 1. - g*g, dp/2. - 1.5) / sqrt(PI);
 return result;
}


double mxewma_ad_new(double lambda, double ce, int p, double delta, int N, int psi_type, double hs, int qtype)
{ double *PSI, *ARL, *w1, *z1, *w2, *z2, *w3, *z3, zahl, ad, psi0, LL, term, gj, korr, korr2, norm;
  int i, j, N2;  
  
 PSI = vector(N);
 w1  = vector(N);
 z1  = vector(N);
 
 if ( hs < 0. ) hs = 0.;
 
 if ( psi_type == cond ) zahl = mxewma_psi (lambda, ce, p, N, PSI, w1, z1);
 if ( psi_type == cycl ) zahl = mxewma_psiS(lambda, ce, p, hs, N, PSI, w1, z1);
 
 norm = sinh(1.);
 
 ad = 0.;
 
 if ( fabs(delta)<1e-10 ) { /* in-control */
   ad = -2.;   
 } else { /* out-of-control */
     
   N2  = N*N;
   ARL = vector(N2);
   w2  = vector(N);   
   z2  = vector(N);
   w3  = vector(N);
   z3  = vector(N);
   
   if ( qtype == nGL1 ) LL = mxewma_arl_f_1q(lambda, ce, p, delta, N, ARL, w2, z2, w3, z3);
   if ( qtype == nGL2 ) LL = mxewma_arl_f_1r(lambda, ce, p, delta, N, ARL, w2, z2, w3, z3);
   if ( qtype == nGL3 ) LL = mxewma_arl_f_1s(lambda, ce, p, delta, N, ARL, w2, z2, w3, z3);   
   if ( qtype == nGL4 ) LL = mxewma_arl_f_1t(lambda, ce, p, delta, N, ARL, w2, z2, w3, z3);
   if ( qtype == nGL5 ) LL = mxewma_arl_f_1u(lambda, ce, p, delta, N, ARL, w2, z2, w3, z3);
   
   for (i=0; i<N; i++) {
      term = 0.;
      if ( qtype == nGL1 ) { korr = 1.; }
      else { korr = 2. * z2[i]; }
      for (j=0; j<N; j++) {
        gj = z3[j];
        korr2 = 1.;
        if ( qtype == nGL3 ) { korr2 = cos(z3[j]); gj = sin(z3[j]); }
        if ( qtype == nGL4 ) { korr2 = 1. / ( cos(z3[j])*cos(z3[j]) ); gj = tan(z3[j]); }
        if ( qtype == nGL5 ) { korr2 = cosh(z3[j]) / norm; gj = sinh(z3[j]); }
        term += w3[j] * cos_unif_sphere(gj, p) * ARL[i*N +j] * korr2;
      }
      ad += term * w2[i] * PSI[i] * korr;
   }

   if ( psi_type == cycl ) {
     psi0 = 1./zahl;
     LL = mxewma_L_of_ag(lambda, ce, p, delta, N, qtype, ARL, 0., 0., w2, z2, w3, z3);
     ad += psi0 * LL;
   }
   
   Free(z3);
   Free(w3);
   Free(z2);
   Free(w2);
   Free(ARL);   
 }
 
 Free(z1);
 Free(w1);
 Free(PSI);
 
 return ad;
}


double mxewma_ad_e(double lambda, double ce, int p, double delta, int psi_type, int N)
{ double *PSI, *ARL, *z, ad, dN=1., ce_=0., w=0., ce2=0.;
  int i, ix, iy, N2, index, dQ;
 
 if ( fabs(delta)<1e-10 ) { /* in-control */
   PSI = vector(N);
   ARL = vector(N);
   z   = vector(N);
   
   if ( psi_type == cond ) ad = mxewma_psi0_e(lambda, ce, p, N, PSI);
   if ( psi_type == cycl ) ad = mxewma_psiS0_e(lambda, ce, p, N, PSI);
   
   ad = mxewma_arl_f_0e(lambda, ce, p, N, ARL, z);
   
   ad = 0.;
   for (i=0; i<N; i++) ad += PSI[i] * ARL[i];
   
   Free(z);
 } else { /* out-of-control */   
   dN = (double)N;
   ce_ = sqrt( ce * lambda/(2.-lambda) ); 
   w = 2.*ce_/(2.*dN+1.);
   ce2 = ce_*ce_;
   N2 = 2*N+1;
   index = 0;
   for (ix=0; ix<N2; ix++) for (iy=0; iy<N+1; iy++) index += ( (ix-dN)*(ix-dN) + iy*iy < ce2/(w*w) );   
   PSI = vector(index);
   ARL = vector(index);
   
   if ( psi_type == cond ) ad = mxewma_psi1_e(lambda, ce, p, N, PSI);
   if ( psi_type == cycl ) ad = mxewma_psiS1_e(lambda, ce, p, N, PSI);
   
   ad = mxewma_arl_f_1e(lambda, ce, p, delta, N, ARL, &dQ);
     
   ad = 0.;
   for (i=0; i<index; i++) ad += PSI[i] * ARL[i];
 }    

 Free(ARL);
 Free(PSI);
 
 return ad; 
}


double xseU_arl(double lx, double ls, double cx, double cs, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm)
{ double *Sx, *Pnx, *wx, *zx, *p0x, *p0, *S1s, *S2s, *Pns, *ws, *zs, *p0s, q, *zch, *rside, za=0., s2,
         arl_minus=0., arl, arl_plus=0., mn_minus=1., mn_plus=0.,
         mn_minusx, mn_minuss, mn_plusx, mn_pluss, ddf, xl, xu,
         oben, unten;
  int i, j, k, n, *ps;

 cx  *= sqrt( lx/(2.-lx) ); 
 hsx *= sqrt( lx/(2.-lx) );

 s2 = sigma*sigma;
 ddf = (double)df;

 Sx  = matrix(Nx,Nx);
 wx  = vector(Nx);
 zx  = vector(Nx);
 Pnx = matrix(nmax,Nx);
 p0x = vector(nmax);

 S1s = matrix(Ns,Ns);
 S2s = matrix(Ns,Ns);
 ps = ivector(Ns);
 zch = vector(Ns);
 rside = vector(Ns);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,Ns);
 p0s = vector(nmax);

 p0  = vector(nmax);

 gausslegendre(Nx,-cx,cx,zx,wx);

 for (i=0;i<Nx;i++) {
   za = (1.-lx)*zx[i];
   for (j=0;j<Nx;j++)
     Sx[i*Nx+j] = wx[j]/lx*phi( ((zx[j]-za)/lx-mu)/sigma, 0.)/sigma;
 }  

/* Chebyshev nodes on [0,cs] */
 for (i=0;i<Ns;i++) 
   zch[i] = cs/2.*(1.+cos(PI*(2.*(i+1.)-1.)/2./(double)Ns) );

/* P(L>1)(zch[i]) */
 for (i=0;i<Ns;i++)
   rside[i] = CHI( ddf/s2*(cs-(1.-ls)*zch[i])/ls, df);

 for (i=0;i<Ns;i++) {
   za = (1.-ls)*zch[i];
   if (df==2) { xl = za; xu = cs; }
   else       { xl = 0.; xu = sqrt(cs-za); }
   gausslegendre(qm,xl,xu,zs,ws);
   for (j=0;j<Ns;j++) {
     S1s[i*Ns+j] = 0.;
     for (k=0;k<qm;k++)
       if (df==2)
         S1s[i*Ns+j] += ws[k]*Tn((2.*zs[k]-cs)/cs, j) * exp((za-zs[k])/s2/ls); 
       else
         S1s[i*Ns+j] += ws[k]*Tn((2.*(zs[k]*zs[k]+za)-cs)/cs, j)
                      *2.*pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/ls);
     if (df==2) S1s[i*Ns+j] /= s2*ls;
     else       S1s[i*Ns+j] /= gammafn(ddf/2.) * pow(2.*s2*ls/ddf,ddf/2.);
   }
 }

 for (i=0;i<Ns;i++)
   for (j=0;j<Ns;j++) S2s[i*Ns+j] = Tn( (2.*zch[i]-cs)/cs, j);

 LU_decompose(S2s,ps,Ns);

 arl = 1.;
 for (n=1;n<=nmax;n++) {

   if (n==1)
     for (i=0;i<Nx;i++)
       Pnx[i] = PHI( (( cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.) - 
                PHI( ((-cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.);
   else
     for (i=0;i<Nx;i++) {
       Pnx[(n-1)*Nx+i] = 0.;
       for (j=0;j<Nx;j++)
         Pnx[(n-1)*Nx+i] += Sx[i*Nx+j] * Pnx[(n-2)*Nx+j];
     }

   p0x[n-1] = 0.;
   if (n==1)
     p0x[0] = PHI( (( cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.) - 
              PHI( ((-cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.);
   else
     for (j=0;j<Nx;j++)
       p0x[n-1] += wx[j]/lx * phi( ((zx[j]-(1.-lx)*hsx)/lx-mu)/sigma, 0.)/sigma 
                   * Pnx[(n-2)*Nx+j];


   if (n==1)
     for (i=0;i<Ns;i++) {
       Pns[i] = 0.;
       for (j=0;j<Ns;j++)
         Pns[i] += 2./Ns * Tn( (2.*zch[j]-cs)/cs, i) * rside[j];
       if (i==0) Pns[i] /= 2.;
     }
   else {
     for (i=0;i<Ns;i++) {
       rside[i] = 0.;
       for (j=0;j<Ns;j++) rside[i] += S1s[i*Ns+j] * Pns[(n-2)*Ns+j];
     }
     LU_solve2(S2s,rside,ps,Ns);
     for (i=0;i<Ns;i++) Pns[(n-1)*Ns+i] = rside[i];
   }

   p0s[n-1] = 0.;
   if (n==1)
     p0s[0] = CHI(ddf/s2*(cs-(1.-ls)*hss)/ls, df);
   else
     for (j=0;j<Ns;j++)
       p0s[n-1] += Pns[(n-1)*Ns+j] * Tn( (2.*hss-cs)/cs, j);


   p0[n-1] = p0x[n-1] * p0s[n-1];

   mn_minusx = 1.; mn_plusx = 0.;
   mn_minuss = 1.; mn_pluss = 0.;
   if (n>1) {
     for (i=0;i<Nx;i++) {
       if (Pnx[(n-1)*Nx+i]==0)
         if (Pnx[(n-1)*Nx+i]==0) q = 0.;
         else q = 1.;
       else q = Pnx[(n-1)*Nx+i]/Pnx[(n-2)*Nx+i];
      if ( q<mn_minusx ) mn_minusx = q;
      if ( q>mn_plusx ) mn_plusx = q;
     }

     for (i=0;i<Ns;i++) {
       oben = 0.; unten = 0.;
       for (j=0;j<Ns;j++) {
         oben += Pns[(n-1)*Ns+j] * Tn( (2.*zch[i]-cs)/cs, j);
         unten+= Pns[(n-2)*Ns+j] * Tn( (2.*zch[i]-cs)/cs, j);
       }
       if (fabs(unten)<1e-16)
         if (fabs(oben)<1e-16) q = 0.;
         else q = 1.;
       else q = oben/unten;
      if ( q<mn_minuss ) mn_minuss = q;
      if ( q>mn_pluss ) mn_pluss = q;
     }

     mn_minus = mn_minusx * mn_minuss;
     mn_plus  = mn_plusx * mn_pluss;

     arl_minus = arl + p0[n-1]/(1.-mn_minus);
     arl_plus = arl + p0[n-1]/(1.-mn_plus);
   }
   arl += p0[n-1];
   if ( fabs( (arl_plus-arl_minus)/arl_minus )<FINALeps ) n = nmax+1;
 }

 Free(p0);

 Free(p0s);
 Free(Pns);
 Free(zs);
 Free(ws);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);
 
 Free(p0x);
 Free(Pnx);
 Free(zx);
 Free(wx);
 Free(Sx);

 return (arl_plus+arl_minus)/2.;
}


double xseU_sf(double lx, double ls, double cx, double cs, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm, double *p0)
{ double *Sx, *Pnx, *wx, *zx, *p0x, *S1s, *S2s, *Pns, *ws, *zs, *p0s, *zch, *rside, za=0., s2, ddf, xl, xu;
  int i, j, k, n, *ps;

 cx  *= sqrt( lx/(2.-lx) ); 
 hsx *= sqrt( lx/(2.-lx) );

 s2 = sigma*sigma;
 ddf = (double)df;

 Sx  = matrix(Nx,Nx);
 wx  = vector(Nx);
 zx  = vector(Nx);
 Pnx = matrix(nmax,Nx);
 p0x = vector(nmax);

 S1s = matrix(Ns,Ns);
 S2s = matrix(Ns,Ns);
 ps = ivector(Ns);
 zch = vector(Ns);
 rside = vector(Ns);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,Ns);
 p0s = vector(nmax);

 gausslegendre(Nx,-cx,cx,zx,wx);

 for (i=0;i<Nx;i++) {
   za = (1.-lx)*zx[i];
   for (j=0;j<Nx;j++)
     Sx[i*Nx+j] = wx[j]/lx*phi( ((zx[j]-za)/lx-mu)/sigma, 0.)/sigma;
 }  

/* Chebyshev nodes on [0,cs] */
 for (i=0;i<Ns;i++) 
   zch[i] = cs/2.*(1.+cos(PI*(2.*(i+1.)-1.)/2./(double)Ns) );

/* P(L>1)(zch[i]) */
 for (i=0;i<Ns;i++)
   rside[i] = CHI( ddf/s2*(cs-(1.-ls)*zch[i])/ls, df);

 for (i=0;i<Ns;i++) {
   za = (1.-ls)*zch[i];
   if (df==2) { xl = za; xu = cs; }
   else       { xl = 0.; xu = sqrt(cs-za); }
   gausslegendre(qm,xl,xu,zs,ws);
   for (j=0;j<Ns;j++) {
     S1s[i*Ns+j] = 0.;
     for (k=0;k<qm;k++)
       if (df==2)
         S1s[i*Ns+j] += ws[k]*Tn((2.*zs[k]-cs)/cs, j) * exp((za-zs[k])/s2/ls); 
       else
         S1s[i*Ns+j] += ws[k]*Tn((2.*(zs[k]*zs[k]+za)-cs)/cs, j)
                      *2.*pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/ls);
     if (df==2) S1s[i*Ns+j] /= s2*ls;
     else       S1s[i*Ns+j] /= gammafn(ddf/2.) * pow(2.*s2*ls/ddf,ddf/2.);
   }
 }

 for (i=0;i<Ns;i++)
   for (j=0;j<Ns;j++) S2s[i*Ns+j] = Tn( (2.*zch[i]-cs)/cs, j);

 LU_decompose(S2s,ps,Ns);

 for (n=1;n<=nmax;n++) {

   if (n==1)
     for (i=0;i<Nx;i++)
       Pnx[i] = PHI( (( cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.) - 
                PHI( ((-cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.);
   else
     for (i=0;i<Nx;i++) {
       Pnx[(n-1)*Nx+i] = 0.;
       for (j=0;j<Nx;j++)
         Pnx[(n-1)*Nx+i] += Sx[i*Nx+j] * Pnx[(n-2)*Nx+j];
     }

   p0x[n-1] = 0.;
   if (n==1)
     p0x[0] = PHI( (( cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.) - 
              PHI( ((-cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.);
   else
     for (j=0;j<Nx;j++)
       p0x[n-1] += wx[j]/lx * phi( ((zx[j]-(1.-lx)*hsx)/lx-mu)/sigma, 0.)/sigma 
                   * Pnx[(n-2)*Nx+j];

   if (n==1)
     for (i=0;i<Ns;i++) {
       Pns[i] = 0.;
       for (j=0;j<Ns;j++)
         Pns[i] += 2./Ns * Tn( (2.*zch[j]-cs)/cs, i) * rside[j];
       if (i==0) Pns[i] /= 2.;
     }
   else {
     for (i=0;i<Ns;i++) {
       rside[i] = 0.;
       for (j=0;j<Ns;j++) rside[i] += S1s[i*Ns+j] * Pns[(n-2)*Ns+j];
     }
     LU_solve2(S2s,rside,ps,Ns);
     for (i=0;i<Ns;i++) Pns[(n-1)*Ns+i] = rside[i];
   }

   p0s[n-1] = 0.;
   if (n==1)
     p0s[0] = CHI(ddf/s2*(cs-(1.-ls)*hss)/ls, df);
   else
     for (j=0;j<Ns;j++)
       p0s[n-1] += Pns[(n-1)*Ns+j] * Tn( (2.*hss-cs)/cs, j);


   p0[n-1] = p0x[n-1] * p0s[n-1];
 }

 Free(p0s);
 Free(Pns);
 Free(zs);
 Free(ws);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);
 
 Free(p0x);
 Free(Pnx);
 Free(zx);
 Free(wx);
 Free(Sx);

 return 0;
}


double xseU_sf_deluxe(double lx, double ls, double cx, double cs, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm, double *p0, int *nstop, double *rho)
{ double *Sx, *Pnx, *wx, *zx, *p0x, *S1s, *S2s, *Pns, *ws, *zs, *p0s, q, *zch, *rside, za=0., s2,
         mn_minus=1., mn_plus=0., mn_minusx, mn_minuss, mn_plusx, mn_pluss, ddf, xl, xu, oben, unten;
  int i, j, k, n, *ps;

 cx  *= sqrt( lx/(2.-lx) ); 
 hsx *= sqrt( lx/(2.-lx) );

 s2 = sigma*sigma;
 ddf = (double)df;

 Sx  = matrix(Nx,Nx);
 wx  = vector(Nx);
 zx  = vector(Nx);
 Pnx = matrix(nmax,Nx);
 p0x = vector(nmax);

 S1s = matrix(Ns,Ns);
 S2s = matrix(Ns,Ns);
 ps = ivector(Ns);
 zch = vector(Ns);
 rside = vector(Ns);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,Ns);
 p0s = vector(nmax);

 gausslegendre(Nx,-cx,cx,zx,wx);
 
 *nstop = 0;

 for (i=0;i<Nx;i++) {
   za = (1.-lx)*zx[i];
   for (j=0;j<Nx;j++)
     Sx[i*Nx+j] = wx[j]/lx*phi( ((zx[j]-za)/lx-mu)/sigma, 0.)/sigma;
 }  

/* Chebyshev nodes on [0,cs] */
 for (i=0;i<Ns;i++) 
   zch[i] = cs/2.*(1.+cos(PI*(2.*(i+1.)-1.)/2./(double)Ns) );

/* P(L>1)(zch[i]) */
 for (i=0;i<Ns;i++)
   rside[i] = CHI( ddf/s2*(cs-(1.-ls)*zch[i])/ls, df);

 for (i=0;i<Ns;i++) {
   za = (1.-ls)*zch[i];
   if (df==2) { xl = za; xu = cs; }
   else       { xl = 0.; xu = sqrt(cs-za); }
   gausslegendre(qm,xl,xu,zs,ws);
   for (j=0;j<Ns;j++) {
     S1s[i*Ns+j] = 0.;
     for (k=0;k<qm;k++)
       if (df==2)
         S1s[i*Ns+j] += ws[k]*Tn((2.*zs[k]-cs)/cs, j) * exp((za-zs[k])/s2/ls); 
       else
         S1s[i*Ns+j] += ws[k]*Tn((2.*(zs[k]*zs[k]+za)-cs)/cs, j)
                      *2.*pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/ls);
     if (df==2) S1s[i*Ns+j] /= s2*ls;
     else       S1s[i*Ns+j] /= gammafn(ddf/2.) * pow(2.*s2*ls/ddf,ddf/2.);
   }
 }

 for (i=0;i<Ns;i++)
   for (j=0;j<Ns;j++) S2s[i*Ns+j] = Tn( (2.*zch[i]-cs)/cs, j);

 LU_decompose(S2s,ps,Ns);

 for (n=1;n<=nmax;n++) {

   if (n==1)
     for (i=0;i<Nx;i++)
       Pnx[i] = PHI( (( cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.) - 
                PHI( ((-cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.);
   else
     for (i=0;i<Nx;i++) {
       Pnx[(n-1)*Nx+i] = 0.;
       for (j=0;j<Nx;j++)
         Pnx[(n-1)*Nx+i] += Sx[i*Nx+j] * Pnx[(n-2)*Nx+j];
     }

   p0x[n-1] = 0.;
   if (n==1)
     p0x[0] = PHI( (( cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.) - 
              PHI( ((-cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.);
   else
     for (j=0;j<Nx;j++)
       p0x[n-1] += wx[j]/lx * phi( ((zx[j]-(1.-lx)*hsx)/lx-mu)/sigma, 0.)/sigma 
                   * Pnx[(n-2)*Nx+j];


   if (n==1)
     for (i=0;i<Ns;i++) {
       Pns[i] = 0.;
       for (j=0;j<Ns;j++)
         Pns[i] += 2./Ns * Tn( (2.*zch[j]-cs)/cs, i) * rside[j];
       if (i==0) Pns[i] /= 2.;
     }
   else {
     for (i=0;i<Ns;i++) {
       rside[i] = 0.;
       for (j=0;j<Ns;j++) rside[i] += S1s[i*Ns+j] * Pns[(n-2)*Ns+j];
     }
     LU_solve2(S2s,rside,ps,Ns);
     for (i=0;i<Ns;i++) Pns[(n-1)*Ns+i] = rside[i];
   }

   p0s[n-1] = 0.;
   if (n==1)
     p0s[0] = CHI(ddf/s2*(cs-(1.-ls)*hss)/ls, df);
   else
     for (j=0;j<Ns;j++)
       p0s[n-1] += Pns[(n-1)*Ns+j] * Tn( (2.*hss-cs)/cs, j);

   p0[n-1] = p0x[n-1] * p0s[n-1];

   mn_minusx = 1.; mn_plusx = 0.;
   mn_minuss = 1.; mn_pluss = 0.;
   if ( n > 1 ) {
     for (i=0;i<Nx;i++) {
       if (Pnx[(n-1)*Nx+i]==0)
         if (Pnx[(n-1)*Nx+i]==0) q = 0.;
         else q = 1.;
       else q = Pnx[(n-1)*Nx+i]/Pnx[(n-2)*Nx+i];
      if ( q<mn_minusx ) mn_minusx = q;
      if ( q>mn_plusx ) mn_plusx = q;
     }

     for (i=0;i<Ns;i++) {
       oben = 0.; unten = 0.;
       for (j=0;j<Ns;j++) {
         oben += Pns[(n-1)*Ns+j] * Tn( (2.*zch[i]-cs)/cs, j);
         unten+= Pns[(n-2)*Ns+j] * Tn( (2.*zch[i]-cs)/cs, j);
       }
       if (fabs(unten)<1e-16)
         if (fabs(oben)<1e-16) q = 0.;
         else q = 1.;
       else q = oben/unten;
      if ( q<mn_minuss ) mn_minuss = q;
      if ( q>mn_pluss ) mn_pluss = q;
     }

     mn_minus = mn_minusx * mn_minuss;
     mn_plus  = mn_plusx * mn_pluss;

     *rho = (mn_minus + mn_plus)/2.;
     if ( fabs(mn_plus - mn_minus) < FINALeps ) {       
       *nstop = n;
       n = nmax + 1;
     }     
   } /* n > 1 */
 } /* n=1; n<=nmax; n++ */ 

 Free(p0s);
 Free(Pns);
 Free(zs);
 Free(ws);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);
 
 Free(p0x);
 Free(Pnx);
 Free(zx);
 Free(wx);
 Free(Sx);

 return 0;
}


double xseU_Wq(double lx, double ls, double cx, double cs, double p, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm)
{ double *p0, *Sx, *Pnx, *wx, *zx, *p0x, *S1s, *S2s, *Pns, *ws, *zs, *p0s, q, *zch, *rside, za=0., s2,
         mn_minus=1., mn_plus=0., mn_minusx, mn_minuss, mn_plusx, mn_pluss, ddf, xl, xu, oben, unten, q_minus=0., q_plus=0., enumerator=0., Wq=0.;
  int i, j, k, n, *ps;

 cx  *= sqrt( lx/(2.-lx) ); 
 hsx *= sqrt( lx/(2.-lx) );

 s2 = sigma*sigma;
 ddf = (double)df;

 Sx  = matrix(Nx,Nx);
 wx  = vector(Nx);
 zx  = vector(Nx);
 Pnx = matrix(nmax,Nx);
 p0x = vector(nmax);

 S1s = matrix(Ns,Ns);
 S2s = matrix(Ns,Ns);
 ps = ivector(Ns);
 zch = vector(Ns);
 rside = vector(Ns);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,Ns);
 p0s = vector(nmax);
 
 p0 = vector(nmax);

 gausslegendre(Nx,-cx,cx,zx,wx);

 for (i=0;i<Nx;i++) {
   za = (1.-lx)*zx[i];
   for (j=0;j<Nx;j++)
     Sx[i*Nx+j] = wx[j]/lx*phi( ((zx[j]-za)/lx-mu)/sigma, 0.)/sigma;
 }  

/* Chebyshev nodes on [0,cs] */
 for (i=0;i<Ns;i++) 
   zch[i] = cs/2.*(1.+cos(PI*(2.*(i+1.)-1.)/2./(double)Ns) );

/* P(L>1)(zch[i]) */
 for (i=0;i<Ns;i++)
   rside[i] = CHI( ddf/s2*(cs-(1.-ls)*zch[i])/ls, df);

 for (i=0;i<Ns;i++) {
   za = (1.-ls)*zch[i];
   if (df==2) { xl = za; xu = cs; }
   else       { xl = 0.; xu = sqrt(cs-za); }
   gausslegendre(qm,xl,xu,zs,ws);
   for (j=0;j<Ns;j++) {
     S1s[i*Ns+j] = 0.;
     for (k=0;k<qm;k++)
       if (df==2)
         S1s[i*Ns+j] += ws[k]*Tn((2.*zs[k]-cs)/cs, j) * exp((za-zs[k])/s2/ls); 
       else
         S1s[i*Ns+j] += ws[k]*Tn((2.*(zs[k]*zs[k]+za)-cs)/cs, j)
                      *2.*pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/ls);
     if (df==2) S1s[i*Ns+j] /= s2*ls;
     else       S1s[i*Ns+j] /= gammafn(ddf/2.) * pow(2.*s2*ls/ddf,ddf/2.);
   }
 }

 for (i=0;i<Ns;i++)
   for (j=0;j<Ns;j++) S2s[i*Ns+j] = Tn( (2.*zch[i]-cs)/cs, j);

 LU_decompose(S2s,ps,Ns);

 for (n=1; n<=nmax; n++) {

   if (n==1)
     for (i=0;i<Nx;i++)
       Pnx[i] = PHI( (( cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.) - 
                PHI( ((-cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.);
   else
     for (i=0;i<Nx;i++) {
       Pnx[(n-1)*Nx+i] = 0.;
       for (j=0;j<Nx;j++)
         Pnx[(n-1)*Nx+i] += Sx[i*Nx+j] * Pnx[(n-2)*Nx+j];
     }

   p0x[n-1] = 0.;
   if (n==1)
     p0x[0] = PHI( (( cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.) - 
              PHI( ((-cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.);
   else
     for (j=0;j<Nx;j++)
       p0x[n-1] += wx[j]/lx * phi( ((zx[j]-(1.-lx)*hsx)/lx-mu)/sigma, 0.)/sigma 
                   * Pnx[(n-2)*Nx+j];


   if (n==1)
     for (i=0;i<Ns;i++) {
       Pns[i] = 0.;
       for (j=0;j<Ns;j++)
         Pns[i] += 2./Ns * Tn( (2.*zch[j]-cs)/cs, i) * rside[j];
       if (i==0) Pns[i] /= 2.;
     }
   else {
     for (i=0;i<Ns;i++) {
       rside[i] = 0.;
       for (j=0;j<Ns;j++) rside[i] += S1s[i*Ns+j] * Pns[(n-2)*Ns+j];
     }
     LU_solve2(S2s,rside,ps,Ns);
     for (i=0;i<Ns;i++) Pns[(n-1)*Ns+i] = rside[i];
   }

   p0s[n-1] = 0.;
   if (n==1)
     p0s[0] = CHI(ddf/s2*(cs-(1.-ls)*hss)/ls, df);
   else
     for (j=0;j<Ns;j++)
       p0s[n-1] += Pns[(n-1)*Ns+j] * Tn( (2.*hss-cs)/cs, j);

   p0[n-1] = p0x[n-1] * p0s[n-1];

   if ( p0[n-1] < 1.-p ) {
     Wq = (double)n;
     n = nmax+1;
   } else { 
     mn_minusx = 1.; mn_plusx = 0.;
     mn_minuss = 1.; mn_pluss = 0.;
     if ( n > 1 ) {
       for (i=0;i<Nx;i++) {
         if (Pnx[(n-1)*Nx+i]==0)
           if (Pnx[(n-1)*Nx+i]==0) q = 0.;
           else q = 1.;
         else q = Pnx[(n-1)*Nx+i]/Pnx[(n-2)*Nx+i];
        if ( q<mn_minusx ) mn_minusx = q;
        if ( q>mn_plusx ) mn_plusx = q;
       }

       for (i=0;i<Ns;i++) {
         oben = 0.; unten = 0.;
         for (j=0;j<Ns;j++) {
           oben += Pns[(n-1)*Ns+j] * Tn( (2.*zch[i]-cs)/cs, j);
           unten+= Pns[(n-2)*Ns+j] * Tn( (2.*zch[i]-cs)/cs, j);
         }
         if (fabs(unten)<1e-16)
           if (fabs(oben)<1e-16) q = 0.;
           else q = 1.;
         else q = oben/unten;
        if ( q<mn_minuss ) mn_minuss = q;
        if ( q>mn_pluss ) mn_pluss = q;
       }
  
       mn_minus = mn_minusx * mn_minuss;
       mn_plus  = mn_plusx * mn_pluss;
     
       enumerator = log( (1.-p)/p0[n-1] );
       q_minus = (double)n + enumerator/log(mn_minus);
       q_plus  = (double)n + enumerator/log(mn_plus);
       if ( fabs( ceil(q_plus) - ceil(q_minus) ) < .5 ) {
	 Wq = ceil(q_plus);
	 n = nmax +1;
       }
       
     } /* n > 1 */
   } /* p0[n-1] >= 1.-p */
 } /* n=1; n<=nmax; n++ */ 
 
 Free(p0);

 Free(p0s);
 Free(Pns);
 Free(zs);
 Free(ws);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);
 
 Free(p0x);
 Free(Pnx);
 Free(zx);
 Free(wx);
 Free(Sx);

 return Wq;
}


int xseU_crit(double lx, double ls, double L0, double *cx, double *cs, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm)
{ double x1, x2, dx, s1, s2, ds, xARL1, xARL2, sARL1, sARL2, xsARL22, xsARL12, xsARL21,
         f11, f22, f21, f12, d11, d22, d21, d12, nenner, zr=0., c0=-1.;

 x1 = xe_crit(ewma2,lx,2.*L0,zr,hsx,mu,fix,Nx,c0) - .1;
 x2 = x1 + .1;
 s1 = seU_crit(ls,2.*L0,hss,sigma,df,Ns,qm);
 s2 = s1 + .05;

 xARL2 = xe2_iglarl(lx,x2,hsx,mu,Nx);
 sARL2 = seU_iglarl(ls,s2,hss,sigma,df,Ns,qm);
 xsARL22 = xseU_arl(lx,ls,x2,s2,hsx,hss,mu,sigma,df,Nx,Ns,nmax,qm);
 do {
   xARL1 = xe2_iglarl(lx,x1,hsx,mu,Nx);
   sARL1 = seU_iglarl(ls,s1,hss,sigma,df,Ns,qm);
   xsARL21 = xseU_arl(lx,ls,x2,s1,hsx,hss,mu,sigma,df,Nx,Ns,nmax,qm);
   xsARL12 = xseU_arl(lx,ls,x1,s2,hsx,hss,mu,sigma,df,Nx,Ns,nmax,qm);

   /* difference quotient */
   f11 = (xsARL22 - xsARL12)/(x2-x1); f12 = (xsARL22 - xsARL21)/(s2-s1);
   f21 = (xARL2   -   xARL1)/(x2-x1); f22 = (sARL1   -   sARL2)/(s2-s1);

   /* inverse of the difference quotient */
   nenner = f11*f22 - f12*f21;
   d11 =  f22/nenner;  d12 = -f12/nenner;
   d21 = -f21/nenner;  d22 =  f11/nenner;

   dx = d11*(xsARL22-L0) + d12*(xARL2-sARL2);
   ds = d21*(xsARL22-L0) + d22*(xARL2-sARL2);

   x1 = x2;   s1 = s2;
   x2 -= dx;  s2 -= ds;

   xARL2 = xe2_iglarl(lx,x2,hsx,mu,Nx);
   sARL2 = seU_iglarl(ls,s2,hss,sigma,df,Ns,qm);
   xsARL22 = xseU_arl(lx,ls,x2,s2,hsx,hss,mu,sigma,df,Nx,Ns,nmax,qm);
 } while (  (fabs(L0-xsARL22)>1e-6 || fabs(xARL2-sARL2)>1e-6) && (fabs(x2-x1)>1e-8 || fabs(s2-s1)>1e-8)  );

 *cx = x2; *cs = s2;

 return 0;
}


int xseU_q_crit(double lx, double ls, int L0, double alpha, double *cx, double *cs, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int qm, double c_error, double a_error)

{ double x1, x2, dx, s1, s2, ds, xp1, xp2, sp1, sp2, xsp22, xsp12, xsp21,
         f11, f22, f21, f12, d11, d22, d21, d12, nenner, zr=0., *SF;
  int result=1;

 SF  = vector(L0); 
   
 x1 = xe_q_crit(ewma2, lx, L0, 1. - sqrt(1.-alpha), zr, hsx, mu, fix, Nx, c_error, a_error);
 x2 = x1 + .1;
 s1 = seU_q_crit(ls, L0, 1. - sqrt(1.-alpha), hss, sigma, df, Ns, qm, c_error, a_error);
 s2 = s1 + .05;
 
 result = xe2_sf(lx, x2, hsx, mu, Nx, L0, SF);
 if ( result != 0 ) warning("trouble with xseU_q_crit calling xe2_sf [package spc]");
 xp2 = 1. - SF[L0-1];
 result = seU_sf(ls, s2, hss, sigma, df, Ns, L0, qm, SF);
 if ( result != 0 ) warning("trouble with xseU_q_crit calling seU_sf [package spc]");
 sp2 = 1. - SF[L0-1];
 result = xseU_sf(lx, ls, x2, s2, hsx, hss, mu, sigma, df, Nx, Ns, L0, qm, SF);
 if ( result != 0 ) warning("trouble with xseU_q_crit calling xseU_sf [package spc]");
 xsp22 = 1. - SF[L0-1];

 do { 
   result = xe2_sf(lx, x1, hsx, mu, Nx, L0, SF);
   if ( result != 0 ) warning("trouble with xseU_q_crit calling xe2_sf [package spc]");
   xp1 = 1. - SF[L0-1];
   result = seU_sf(ls, s1, hss, sigma, df, Ns, L0, qm, SF);
   if ( result != 0 ) warning("trouble with xseU_q_crit calling seU_sf [package spc]");
   sp1 = 1. - SF[L0-1];
   result = xseU_sf(lx, ls, x2, s1, hsx, hss, mu, sigma, df, Nx, Ns, L0, qm, SF);
   if ( result != 0 ) warning("trouble with xseU_q_crit calling xseU_sf [package spc]");
   xsp21 = 1. - SF[L0-1];
   result = xseU_sf(lx, ls, x1, s2, hsx, hss, mu, sigma, df, Nx, Ns, L0, qm, SF);
   if ( result != 0 ) warning("trouble with xseU_q_crit calling xseU_sf [package spc]");
   xsp12 = 1. - SF[L0-1];   

   /* difference quotient */
   f11 = (xsp22 - xsp12)/(x2-x1); f12 = (xsp22 - xsp21)/(s2-s1);
   f21 = (xp2   -   xp1)/(x2-x1); f22 = (sp1   -   sp2)/(s2-s1);

   /* inverse of the difference quotient */
   nenner = f11*f22 - f12*f21;
   d11 =  f22/nenner;  d12 = -f12/nenner;
   d21 = -f21/nenner;  d22 =  f11/nenner;

   dx = d11*(xsp22-alpha) + d12*(xp2-sp2);
   ds = d21*(xsp22-alpha) + d22*(xp2-sp2);

   x1 = x2;   s1 = s2;
   x2 -= dx;  s2 -= ds;

   result = xe2_sf(lx, x2, hsx, mu, Nx, L0, SF);
   if ( result != 0 ) warning("trouble with xseU_q_crit calling xe2_sf [package spc]");
   xp2 = 1. - SF[L0-1];
   result = seU_sf(ls, s2, hss, sigma, df, Ns, L0, qm, SF);
   if ( result != 0 ) warning("trouble with xseU_q_crit calling seU_sf [package spc]");
   sp2 = 1. - SF[L0-1];
   result = xseU_sf(lx, ls, x2, s2, hsx, hss, mu, sigma, df, Nx, Ns, L0, qm, SF);
   if ( result != 0 ) warning("trouble with xseU_q_crit calling xseU_sf [package spc]");
   xsp22 = 1. - SF[L0-1];
   
 } while (  (fabs(alpha - xsp22)>a_error || fabs(xp2-sp2)>a_error) && (fabs(x2-x1)>c_error || fabs(s2-s1)>c_error)  );

 *cx = x2; *cs = s2;
 
 Free(SF);

 return 0;
}


int xse2fu_q_crit(double lx, double ls, int L0, double alpha, double *cx, double *csl, double csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int qm, double c_error, double a_error)
{ double x1, x2, dx, s1, s2, ds, xp1, xp2, sp1, sp2, xsp22, xsp12, xsp21,
         f11, f22, f21, f12, d11, d22, d21, d12, nenner, zr=0., *SF;
  int result=1;

 SF  = vector(L0); 

 x1 = xe_q_crit(ewma2, lx, L0, 1. - sqrt(1.-alpha), zr, hsx, mu, fix, Nx, c_error, a_error);
 x2 = x1 + .05;
 s1 = se2fu_q_crit(ls, L0, 1. - sqrt(1.-alpha), csu, hss, sigma, df, Ns, qm, c_error, a_error);
 s2 = s1 + .05;
 
 result = xe2_sf(lx, x2, hsx, mu, Nx, L0, SF);
 if ( result != 0 ) warning("trouble with xse2fu_q_crit calling xe2_sf [package spc]");
 xp2 = 1. - SF[L0-1];
 result = se2_sf(ls, s2, csu, hss, sigma, df, Ns, L0, qm, SF);
 if ( result != 0 ) warning("trouble with xse2fu_q_crit calling se2_sf [package spc]");
 sp2 = 1. - SF[L0-1];
 result = xse2_sf(lx, ls, x2, s2, csu, hsx, hss, mu, sigma, df, Nx, Ns, L0, qm, SF); 
 if ( result != 0 ) warning("trouble with xse2fu_q_crit calling xse2_sf [package spc]");
 xsp22 = 1. - SF[L0-1];

 do {
   result = xe2_sf(lx, x1, hsx, mu, Nx, L0, SF);
   if ( result != 0 ) warning("trouble with xse2fu_q_crit calling xe2_sf [package spc]");
   xp1 = 1. - SF[L0-1];
   result = se2_sf(ls, s1, csu, hss, sigma, df, Ns, L0, qm, SF);
   if ( result != 0 ) warning("trouble with xse2fu_q_crit calling se2_sf [package spc]");
   sp1 = 1. - SF[L0-1];
   result = xse2_sf(lx, ls, x2, s1, csu, hsx, hss, mu, sigma, df, Nx, Ns, L0, qm, SF);
   if ( result != 0 ) warning("trouble with xse2fu_q_crit calling xse2_sf [package spc]");
   xsp21 = 1. - SF[L0-1];
   result = xse2_sf(lx, ls, x1, s2, csu, hsx, hss, mu, sigma, df, Nx, Ns, L0, qm, SF);
   if ( result != 0 ) warning("trouble with xse2fu_q_crit calling xse2_sf [package spc]");
   xsp12 = 1. - SF[L0-1]; 

   /* difference quotient */
   f11 = (xsp22 - xsp12)/(x2-x1); f12 = (xsp22 - xsp21)/(s2-s1);
   f21 = (xp2   -   xp1)/(x2-x1); f22 = (sp1   -   sp2)/(s2-s1);

   /* inverse of the difference quotient */
   nenner = f11*f22 - f12*f21;
   d11 =  f22/nenner;  d12 = -f12/nenner;
   d21 = -f21/nenner;  d22 =  f11/nenner;

   dx = d11*(xsp22-alpha) + d12*(xp2-sp2);
   ds = d21*(xsp22-alpha) + d22*(xp2-sp2);

   x1 = x2;   s1 = s2;
   x2 -= dx;  s2 -= ds;

   result = xe2_sf(lx, x2, hsx, mu, Nx, L0, SF);
   if ( result != 0 ) warning("trouble with xse2fu_q_crit calling xe2_sf [package spc]");
   xp2 = 1. - SF[L0-1];
   result = se2_sf(ls, s2, csu, hss, sigma, df, Ns, L0, qm, SF);
   if ( result != 0 ) warning("trouble with xse2fu_q_crit calling se2_sf [package spc]");
   sp2 = 1. - SF[L0-1];
   result = xse2_sf(lx, ls, x2, s2, csu, hsx, hss, mu, sigma, df, Nx, Ns, L0, qm, SF);
   if ( result != 0 ) warning("trouble with xse2fu_q_crit calling xse2_sf [package spc]");
   xsp22 = 1. - SF[L0-1];
 } while (  (fabs(alpha - xsp22)>a_error || fabs(xp2-sp2)>a_error) && (fabs(x2-x1)>c_error || fabs(s2-s1)>c_error)  );

 *cx = x2; *csl = s2;
 
 Free(SF);

 return 0;
}


int xse2_q_crit(double lx, double ls, int L0, double alpha, double *cx, double *csl, double *csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int qm, double c_error, double a_error)
{ double s1, s2, s3, ds, sl1, sl2, sl3, Lm, Lp, x, cl, *SF;
  int result=1;

 SF  = vector(L0);

 cl = 0.; 
 result = xseU_q_crit(lx, ls, L0, alpha, &x, &s1, hsx, hss, mu, sigma, df, Nx, Ns, qm, c_error, a_error);
 if ( result != 0 ) warning("trouble with xse2_q_crit calling xseU_q_crit [package spc]");
 result = xseU_sf(lx, ls, x, s1, hsx, hss, mu, sigma-lmEPS, df, Nx, Ns, L0, qm, SF);
 if ( result != 0 ) warning("trouble with xse2_q_crit calling xseU_sf [package spc]");
 Lm = 1. - SF[L0-1];
 result = xseU_sf(lx, ls, x, s1, hsx, hss, mu, sigma+lmEPS, df, Nx, Ns, L0, qm, SF);
 if ( result != 0 ) warning("trouble with xse2_q_crit calling xseU_sf [package spc]");
 Lp = 1. - SF[L0-1];
 sl1 = (Lp-Lm)/(2.*lmEPS);

 s2 = s1 + .15; 
 result = xse2fu_q_crit(lx, ls, L0, alpha, &x, &cl, s2, hsx, hss, mu, sigma, df, Nx, Ns, qm, c_error, a_error);
 if ( result != 0 ) warning("trouble with xse2_q_crit calling xse2fu_q_crit [package spc]");
 result = xse2_sf(lx, ls, x, cl, s2, hsx, hss, mu, sigma-lmEPS, df, Nx, Ns, L0, qm, SF);
 if ( result != 0 ) warning("trouble with xse2_q_crit calling xse2_sf [package spc]");
 Lm = 1. - SF[L0-1];
 result = xse2_sf(lx, ls, x, cl, s2, hsx, hss, mu, sigma+lmEPS, df, Nx, Ns, L0, qm, SF);
 if ( result != 0 ) warning("trouble with xse2_q_crit calling xse2_sf [package spc]");
 Lp = 1. - SF[L0-1];
 sl2 = (Lp-Lm)/(2.*lmEPS);
 
 do {
   s3 = s1 - sl1/(sl2-sl1) * (s2-s1);   
   result = xse2fu_q_crit(lx, ls, L0, alpha, &x, &cl, s3, hsx, hss, mu, sigma, df, Nx, Ns, qm, c_error, a_error); 
   if ( result != 0 ) warning("trouble with xse2_q_crit calling xse2fu_q_crit [package spc]");
   result = xse2_sf(lx, ls, x, cl, s3, hsx, hss, mu, sigma-lmEPS, df, Nx, Ns, L0, qm, SF);
   if ( result != 0 ) warning("trouble with xse2_q_crit calling xse2_sf [package spc]");
   Lm = 1. - SF[L0-1];
   result = xse2_sf(lx, ls, x, cl, s3, hsx, hss, mu, sigma+lmEPS, df, Nx, Ns, L0, qm, SF);
   if ( result != 0 ) warning("trouble with xse2_q_crit calling xse2_sf [package spc]");
   Lp = 1. - SF[L0-1];
   sl3 = (Lp-Lm)/(2.*lmEPS);
   ds = s3-s2; s1 = s2; sl1 = sl2; s2 = s3; sl2 = sl3;
 } while ( fabs(sl3)>a_error && fabs(ds)>c_error );

 *cx = x; *csl = cl; *csu = s3;

 Free(SF);
 
 return 0;
}


double xse2_arl(double lx, double ls, double cx, double csl, double csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm)
{ double *Sx, *Pnx, *wx, *zx, *p0x, *p0, *S1s, *S2s, *Pns, *ws, *zs, *p0s, q, *zch, *rside, *b, za=0., s2, dN, Hij,
         arl_minus=0., arl, arl_plus=0., mn_minus=1., mn_plus=0., mn_minusx, mn_minuss, mn_plusx, mn_pluss, ddf, xl, xu, oben, unten;
  int i, j, k, n, *ps, Ntilde, ihs, M, NN, ii, jj;

 cx  *= sqrt( lx/(2.-lx) ); 
 hsx *= sqrt( lx/(2.-lx) );

 M = ceil( (log(csl)-log(csu))/log(1.-ls) );
 Ntilde = ceil( (double)Ns/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde;

 ihs = floor( (log(csl)-log(hss))/log(1.-ls) );
 if (ihs<0) ihs = 0;

 Sx  = matrix(Nx,Nx);
 wx  = vector(Nx);
 zx  = vector(Nx);
 Pnx = matrix(nmax,Nx);
 p0x = vector(nmax); 

 S1s = matrix(NN,NN);
 S2s = matrix(NN,NN);
 ps = ivector(NN);
 zch = matrix(M,Ntilde);
 rside = vector(NN);
 b   = vector(M+1);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,NN);
 p0s = vector(nmax);

 p0  = vector(nmax);

 gausslegendre(Nx,-cx,cx,zx,wx);

 for (i=0;i<Nx;i++) {
   za = (1.-lx)*zx[i];
   for (j=0;j<Nx;j++)
     Sx[i*Nx+j] = wx[j]/lx*phi( ((zx[j]-za)/lx-mu)/sigma, 0.)/sigma;
 }  

/* interval borders b_i = cl/(1-l)^i */
 for (i=0;i<M;i++) b[i] = csl/pow(1.-ls, (double)(i));
 b[M] = csu;

 /* Chebyshev nodes on [b_0,b_1],[b_1,b_2],...,[b_M-1,cu] */
 for (i=0;i<M;i++)
   for (j=0;j<Ntilde;j++) {
     zch[ i*Ntilde+j ] = b[i] + (b[i+1]-b[i])/2.*(1.+cos(PI*(2.*j+1.)/2./dN));
   }

 /* P(L>1)(zch[i,j]) */
 for (i=0;i<M;i++)
   for (j=0;j<Ntilde;j++) {
     rside[ i*Ntilde+j ] =  CHI( ddf/s2*(csu-(1.-ls)*zch[ i*Ntilde+j ])/ls, df)
                          - CHI( ddf/s2*(csl-(1.-ls)*zch[ i*Ntilde+j ])/ls, df);
   }

 for (i=0;i<M;i++)
   for (j=0;j<Ntilde;j++) {
     za = (1.-ls)*zch[ i*Ntilde+j ];
     for (ii=0;ii<M;ii++)
       for (jj=0;jj<Ntilde;jj++) {
         if (b[ii+1]<za) S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = 0.;
         else {
           if (za<b[ii]) xl = b[ii]; else xl = za;
           xu = b[ii+1];
           if (df!=2) {
             xl = sqrt(xl-za);
             xu = sqrt(xu-za);
           }
           gausslegendre(qm,xl,xu,zs,ws);
           Hij = 0.;
           for (k=0;k<qm;k++)
             if (df==2)
               Hij += ws[k]*Tn( (2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                      * exp(-zs[k]/s2/ls);
             else
               Hij +=
            ws[k]*Tn( (2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                 * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/ls);

           if (df==2) Hij *= exp(za/s2/ls)/s2/ls;
           else       Hij /= gammafn(ddf/2.) * pow(2.*s2*ls/ddf, ddf/2.);
           S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = Hij;
         }
       }
   }

 for (i=0;i<NN;i++)
   for (j=0;j<NN;j++) S2s[i*NN+j] = 0.;

 for (i=0;i<M;i++)
   for (j=0;j<Ntilde;j++)
     for (jj=0;jj<Ntilde;jj++)
       S2s[ (i*Ntilde+j)*NN + i*Ntilde+jj ] =
         Tn( (2.*zch[ i*Ntilde+j ]-b[i+1]-b[i])/(b[i+1]-b[i]), jj); 

 LU_decompose(S2s,ps,NN);

 arl = 1.;

 for (n=1;n<=nmax;n++) {

   if (n==1)
     for (i=0;i<Nx;i++)
       Pnx[i] = PHI( (( cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.) - 
                PHI( ((-cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.);
   else
     for (i=0;i<Nx;i++) {
       Pnx[(n-1)*Nx+i] = 0.;
       for (j=0;j<Nx;j++)
         Pnx[(n-1)*Nx+i] += Sx[i*Nx+j] * Pnx[(n-2)*Nx+j];
     }

   p0x[n-1] = 0.;
   if (n==1)
     p0x[0] = PHI( (( cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.) - 
              PHI( ((-cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.);
   else
     for (j=0;j<Nx;j++)
       p0x[n-1] += wx[j]/lx * phi( ((zx[j]-(1.-lx)*hsx)/lx-mu)/sigma, 0.)/sigma 
                   * Pnx[(n-2)*Nx+j];

   if (n==1)
     for (i=0;i<M;i++)
       for (j=0;j<Ntilde;j++) {
         Pns[ i*Ntilde+j ] = 0.;
         for (jj=0;jj<Ntilde;jj++)
           Pns[ i*Ntilde+j ] += /*  usual Chebyshev approximation  */
             2./Ntilde * Tn( (2.*zch[i*Ntilde+jj]-b[i+1]-b[i])/(b[i+1]-b[i]), j)
                  * rside[ i*Ntilde+jj ];
         if (j==0) Pns[ i*Ntilde+j ] /= 2.;
       }
   else {
     for (i=0;i<NN;i++) {
       rside[i] = 0.;
       for (j=0;j<NN;j++) rside[i] += S1s[ i*NN+j ] * Pns[ (n-2)*NN+j ];
     }
     LU_solve2(S2s,rside,ps,NN);
     for (i=0;i<NN;i++) Pns[ (n-1)*NN+i ] = rside[i];
   }

   p0s[n-1] = 0.;
   if (n==1)
     p0s[0] =  CHI( ddf/s2*(csu-(1.-ls)*hss)/ls, df)
             - CHI( ddf/s2*(csl-(1.-ls)*hss)/ls, df);
   else
     for (j=0;j<Ntilde;j++)
       p0s[n-1] += Pns[ (n-1)*NN + ihs*Ntilde+j ]
                * Tn( (2.*hss-b[ihs+1]-b[ihs])/(b[ihs+1]-b[ihs]), j);


   p0[n-1] = p0x[n-1] * p0s[n-1];

   mn_minusx = 1.; mn_plusx = 0.;
   mn_minuss = 1.; mn_pluss = 0.;
   if (n>1) {
     for (i=0;i<Nx;i++) {
       if (Pnx[(n-1)*Nx+i]==0)
         if (Pnx[(n-1)*Nx+i]==0) q = 0.;
         else q = 1.;
       else q = Pnx[(n-1)*Nx+i]/Pnx[(n-2)*Nx+i];
      if ( q<mn_minusx ) mn_minusx = q;
      if ( q>mn_plusx ) mn_plusx = q;
     }

     for (i=0;i<M;i++)
       for (j=0;j<Ntilde;j++) {
         oben = 0.;
         unten = 0.;
         for (jj=0;jj<Ntilde;jj++) {
           oben += Pns[ (n-1)*NN + i*Ntilde+jj ]
                 * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
           unten+= Pns[ (n-2)*NN + i*Ntilde+jj ]
                 * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
         }
         if (fabs(unten)<1e-16)
           if (fabs(oben)<1e-16) q = 0.;
           else q = 1.;
         else q = oben/unten;
         if ( q<mn_minuss ) mn_minuss = q;
         if ( q>mn_pluss ) mn_pluss = q;
       }

     mn_minus = mn_minusx * mn_minuss;
     mn_plus  = mn_plusx * mn_pluss;

     arl_minus = arl + p0[n-1]/(1.-mn_minus);
     arl_plus = arl + p0[n-1]/(1.-mn_plus);
   }
   arl += p0[n-1];

   if ( fabs( (arl_plus-arl_minus)/arl_minus )<FINALeps ) n = nmax+1;
 }

 Free(p0);

 Free(p0s);
 Free(Pns);
 Free(zs);
 Free(ws);
 Free(b);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);

 Free(p0x);
 Free(Pnx);
 Free(zx);
 Free(wx);
 Free(Sx);

 return (arl_plus+arl_minus)/2.;
}


double xse2_sf(double lx, double ls, double cx, double csl, double csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm, double *p0)
{ double *Sx, *Pnx, *wx, *zx, *p0x, *S1s, *S2s, *Pns, *ws, *zs, *p0s, *zch, *rside, *b, za=0., s2, dN, Hij, ddf, xl, xu;
  int i, j, k, n, *ps, Ntilde, ihs, M, NN, ii, jj;

 cx  *= sqrt( lx/(2.-lx) ); 
 hsx *= sqrt( lx/(2.-lx) );

 M = ceil( (log(csl)-log(csu))/log(1.-ls) );
 Ntilde = ceil( (double)Ns/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde;

 ihs = floor( (log(csl)-log(hss))/log(1.-ls) );
 if (ihs<0) ihs = 0;

 Sx  = matrix(Nx,Nx);
 wx  = vector(Nx);
 zx  = vector(Nx);
 Pnx = matrix(nmax,Nx);
 p0x = vector(nmax); 

 S1s = matrix(NN,NN);
 S2s = matrix(NN,NN);
 ps = ivector(NN);
 zch = matrix(M,Ntilde);
 rside = vector(NN);
 b   = vector(M+1);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,NN);
 p0s = vector(nmax);

 gausslegendre(Nx,-cx,cx,zx,wx);

 for (i=0;i<Nx;i++) {
   za = (1.-lx)*zx[i];
   for (j=0;j<Nx;j++)
     Sx[i*Nx+j] = wx[j]/lx*phi( ((zx[j]-za)/lx-mu)/sigma, 0.)/sigma;
 }  

/* interval borders b_i = cl/(1-l)^i */
 for (i=0;i<M;i++) b[i] = csl/pow(1.-ls, (double)(i));
 b[M] = csu;

 /* Chebyshev nodes on [b_0,b_1],[b_1,b_2],...,[b_M-1,cu] */
 for (i=0;i<M;i++)
   for (j=0;j<Ntilde;j++) {
     zch[ i*Ntilde+j ] = b[i] + (b[i+1]-b[i])/2.*(1.+cos(PI*(2.*j+1.)/2./dN));
   }

 /* P(L>1)(zch[i,j]) */
 for (i=0;i<M;i++)
   for (j=0;j<Ntilde;j++) {
     rside[ i*Ntilde+j ] =  CHI( ddf/s2*(csu-(1.-ls)*zch[ i*Ntilde+j ])/ls, df)
                          - CHI( ddf/s2*(csl-(1.-ls)*zch[ i*Ntilde+j ])/ls, df);
   }

 for (i=0;i<M;i++)
   for (j=0;j<Ntilde;j++) {
     za = (1.-ls)*zch[ i*Ntilde+j ];
     for (ii=0;ii<M;ii++)
       for (jj=0;jj<Ntilde;jj++) {
         if (b[ii+1]<za) S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = 0.;
         else {
           if (za<b[ii]) xl = b[ii]; else xl = za;
           xu = b[ii+1];
           if (df!=2) {
             xl = sqrt(xl-za);
             xu = sqrt(xu-za);
           }
           gausslegendre(qm,xl,xu,zs,ws);
           Hij = 0.;
           for (k=0;k<qm;k++)
             if (df==2)
               Hij += ws[k]*Tn( (2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                      * exp(-zs[k]/s2/ls);
             else
               Hij +=
            ws[k]*Tn( (2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                 * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/ls);

           if (df==2) Hij *= exp(za/s2/ls)/s2/ls;
           else       Hij /= gammafn(ddf/2.) * pow(2.*s2*ls/ddf, ddf/2.);
           S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = Hij;
         }
       }
   }

 for (i=0;i<NN;i++)
   for (j=0;j<NN;j++) S2s[i*NN+j] = 0.;

 for (i=0;i<M;i++)
   for (j=0;j<Ntilde;j++)
     for (jj=0;jj<Ntilde;jj++)
       S2s[ (i*Ntilde+j)*NN + i*Ntilde+jj ] =
         Tn( (2.*zch[ i*Ntilde+j ]-b[i+1]-b[i])/(b[i+1]-b[i]), jj); 

 LU_decompose(S2s,ps,NN);

 for (n=1;n<=nmax;n++) {

   if (n==1)
     for (i=0;i<Nx;i++)
       Pnx[i] = PHI( (( cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.) - 
                PHI( ((-cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.);
   else
     for (i=0;i<Nx;i++) {
       Pnx[(n-1)*Nx+i] = 0.;
       for (j=0;j<Nx;j++)
         Pnx[(n-1)*Nx+i] += Sx[i*Nx+j] * Pnx[(n-2)*Nx+j];
     }

   p0x[n-1] = 0.;
   if (n==1)
     p0x[0] = PHI( (( cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.) - 
              PHI( ((-cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.);
   else
     for (j=0;j<Nx;j++)
       p0x[n-1] += wx[j]/lx * phi( ((zx[j]-(1.-lx)*hsx)/lx-mu)/sigma, 0.)/sigma 
                   * Pnx[(n-2)*Nx+j];

   if (n==1)
     for (i=0;i<M;i++)
       for (j=0;j<Ntilde;j++) {
         Pns[ i*Ntilde+j ] = 0.;
         for (jj=0;jj<Ntilde;jj++)
           Pns[ i*Ntilde+j ] += /*  usual Chebyshev approximation  */
             2./Ntilde * Tn( (2.*zch[i*Ntilde+jj]-b[i+1]-b[i])/(b[i+1]-b[i]), j)
                  * rside[ i*Ntilde+jj ];
         if (j==0) Pns[ i*Ntilde+j ] /= 2.;
       }
   else {
     for (i=0;i<NN;i++) {
       rside[i] = 0.;
       for (j=0;j<NN;j++) rside[i] += S1s[ i*NN+j ] * Pns[ (n-2)*NN+j ];
     }
     LU_solve2(S2s,rside,ps,NN);
     for (i=0;i<NN;i++) Pns[ (n-1)*NN+i ] = rside[i];
   }

   p0s[n-1] = 0.;
   if (n==1)
     p0s[0] =  CHI( ddf/s2*(csu-(1.-ls)*hss)/ls, df)
             - CHI( ddf/s2*(csl-(1.-ls)*hss)/ls, df);
   else
     for (j=0;j<Ntilde;j++)
       p0s[n-1] += Pns[ (n-1)*NN + ihs*Ntilde+j ]
                * Tn( (2.*hss-b[ihs+1]-b[ihs])/(b[ihs+1]-b[ihs]), j);


   p0[n-1] = p0x[n-1] * p0s[n-1];   
 }
 
 Free(p0s);
 Free(Pns);
 Free(zs);
 Free(ws);
 Free(b);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);

 Free(p0x);
 Free(Pnx);
 Free(zx);
 Free(wx);
 Free(Sx);

 return 0;
}


double xse2_sf_deluxe(double lx, double ls, double cx, double csl, double csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm, double *p0, int *nstop, double *rho)
{ double *Sx, *Pnx, *wx, *zx, *p0x, *S1s, *S2s, *Pns, *ws, *zs, *p0s, q, *zch, *rside, *b, za=0., s2, dN, Hij,
         mn_minus=1., mn_plus=0., mn_minusx, mn_minuss, mn_plusx, mn_pluss, ddf, xl, xu, oben, unten;
  int i, j, k, n, *ps, Ntilde, ihs, M, NN, ii, jj;

 cx  *= sqrt( lx/(2.-lx) ); 
 hsx *= sqrt( lx/(2.-lx) );

 M = ceil( (log(csl)-log(csu))/log(1.-ls) );
 Ntilde = ceil( (double)Ns/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde;

 ihs = floor( (log(csl)-log(hss))/log(1.-ls) );
 if (ihs<0) ihs = 0;

 Sx  = matrix(Nx,Nx);
 wx  = vector(Nx);
 zx  = vector(Nx);
 Pnx = matrix(nmax,Nx);
 p0x = vector(nmax); 

 S1s = matrix(NN,NN);
 S2s = matrix(NN,NN);
 ps = ivector(NN);
 zch = matrix(M,Ntilde);
 rside = vector(NN);
 b   = vector(M+1);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,NN);
 p0s = vector(nmax);

 gausslegendre(Nx,-cx,cx,zx,wx);
 
 *nstop = 0;

 for (i=0;i<Nx;i++) {
   za = (1.-lx)*zx[i];
   for (j=0;j<Nx;j++)
     Sx[i*Nx+j] = wx[j]/lx*phi( ((zx[j]-za)/lx-mu)/sigma, 0.)/sigma;
 }  

/* interval borders b_i = cl/(1-l)^i */
 for (i=0;i<M;i++) b[i] = csl/pow(1.-ls, (double)(i));
 b[M] = csu;

 /* Chebyshev nodes on [b_0,b_1],[b_1,b_2],...,[b_M-1,cu] */
 for (i=0;i<M;i++)
   for (j=0;j<Ntilde;j++) {
     zch[ i*Ntilde+j ] = b[i] + (b[i+1]-b[i])/2.*(1.+cos(PI*(2.*j+1.)/2./dN));
   }

 /* P(L>1)(zch[i,j]) */
 for (i=0;i<M;i++)
   for (j=0;j<Ntilde;j++) {
     rside[ i*Ntilde+j ] =  CHI( ddf/s2*(csu-(1.-ls)*zch[ i*Ntilde+j ])/ls, df)
                          - CHI( ddf/s2*(csl-(1.-ls)*zch[ i*Ntilde+j ])/ls, df);
   }

 for (i=0;i<M;i++)
   for (j=0;j<Ntilde;j++) {
     za = (1.-ls)*zch[ i*Ntilde+j ];
     for (ii=0;ii<M;ii++)
       for (jj=0;jj<Ntilde;jj++) {
         if (b[ii+1]<za) S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = 0.;
         else {
           if (za<b[ii]) xl = b[ii]; else xl = za;
           xu = b[ii+1];
           if (df!=2) {
             xl = sqrt(xl-za);
             xu = sqrt(xu-za);
           }
           gausslegendre(qm,xl,xu,zs,ws);
           Hij = 0.;
           for (k=0;k<qm;k++)
             if (df==2)
               Hij += ws[k]*Tn( (2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                      * exp(-zs[k]/s2/ls);
             else
               Hij +=
            ws[k]*Tn( (2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                 * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/ls);

           if (df==2) Hij *= exp(za/s2/ls)/s2/ls;
           else       Hij /= gammafn(ddf/2.) * pow(2.*s2*ls/ddf, ddf/2.);
           S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = Hij;
         }
       }
   }

 for (i=0;i<NN;i++)
   for (j=0;j<NN;j++) S2s[i*NN+j] = 0.;

 for (i=0;i<M;i++)
   for (j=0;j<Ntilde;j++)
     for (jj=0;jj<Ntilde;jj++)
       S2s[ (i*Ntilde+j)*NN + i*Ntilde+jj ] =
         Tn( (2.*zch[ i*Ntilde+j ]-b[i+1]-b[i])/(b[i+1]-b[i]), jj); 

 LU_decompose(S2s,ps,NN);

 for (n=1;n<=nmax;n++) {

   if (n==1)
     for (i=0;i<Nx;i++)
       Pnx[i] = PHI( (( cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.) - 
                PHI( ((-cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.);
   else
     for (i=0;i<Nx;i++) {
       Pnx[(n-1)*Nx+i] = 0.;
       for (j=0;j<Nx;j++)
         Pnx[(n-1)*Nx+i] += Sx[i*Nx+j] * Pnx[(n-2)*Nx+j];
     }

   p0x[n-1] = 0.;
   if (n==1)
     p0x[0] = PHI( (( cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.) - 
              PHI( ((-cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.);
   else
     for (j=0;j<Nx;j++)
       p0x[n-1] += wx[j]/lx * phi( ((zx[j]-(1.-lx)*hsx)/lx-mu)/sigma, 0.)/sigma 
                   * Pnx[(n-2)*Nx+j];

   if (n==1)
     for (i=0;i<M;i++)
       for (j=0;j<Ntilde;j++) {
         Pns[ i*Ntilde+j ] = 0.;
         for (jj=0;jj<Ntilde;jj++)
           Pns[ i*Ntilde+j ] += /*  usual Chebyshev approximation  */
             2./Ntilde * Tn( (2.*zch[i*Ntilde+jj]-b[i+1]-b[i])/(b[i+1]-b[i]), j)
                  * rside[ i*Ntilde+jj ];
         if (j==0) Pns[ i*Ntilde+j ] /= 2.;
       }
   else {
     for (i=0;i<NN;i++) {
       rside[i] = 0.;
       for (j=0;j<NN;j++) rside[i] += S1s[ i*NN+j ] * Pns[ (n-2)*NN+j ];
     }
     LU_solve2(S2s,rside,ps,NN);
     for (i=0;i<NN;i++) Pns[ (n-1)*NN+i ] = rside[i];
   }

   p0s[n-1] = 0.;
   if (n==1)
     p0s[0] =  CHI( ddf/s2*(csu-(1.-ls)*hss)/ls, df)
             - CHI( ddf/s2*(csl-(1.-ls)*hss)/ls, df);
   else
     for (j=0;j<Ntilde;j++)
       p0s[n-1] += Pns[ (n-1)*NN + ihs*Ntilde+j ]
                * Tn( (2.*hss-b[ihs+1]-b[ihs])/(b[ihs+1]-b[ihs]), j);


   p0[n-1] = p0x[n-1] * p0s[n-1];

   mn_minusx = 1.; mn_plusx = 0.;
   mn_minuss = 1.; mn_pluss = 0.;
   if ( n>1 ) {
     for (i=0;i<Nx;i++) {
       if (Pnx[(n-1)*Nx+i]==0)
         if (Pnx[(n-1)*Nx+i]==0) q = 0.;
         else q = 1.;
       else q = Pnx[(n-1)*Nx+i]/Pnx[(n-2)*Nx+i];
      if ( q<mn_minusx ) mn_minusx = q;
      if ( q>mn_plusx ) mn_plusx = q;
     }

     for (i=0;i<M;i++)
       for (j=0;j<Ntilde;j++) {
         oben = 0.;
         unten = 0.;
         for (jj=0;jj<Ntilde;jj++) {
           oben += Pns[ (n-1)*NN + i*Ntilde+jj ]
                 * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
           unten+= Pns[ (n-2)*NN + i*Ntilde+jj ]
                 * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
         }
         if (fabs(unten)<1e-16)
           if (fabs(oben)<1e-16) q = 0.;
           else q = 1.;
         else q = oben/unten;
         if ( q<mn_minuss ) mn_minuss = q;
         if ( q>mn_pluss ) mn_pluss = q;
       }

     mn_minus = mn_minusx * mn_minuss;
     mn_plus  = mn_plusx * mn_pluss;

     *rho = (mn_minus + mn_plus)/2.;
     if ( fabs(mn_plus - mn_minus) < FINALeps ) {       
       *nstop = n;
       n = nmax + 1;
     }     
   } /* n > 1 */
 } /* n=1; n<=nmax; n++ */

 Free(p0s);
 Free(Pns);
 Free(zs);
 Free(ws);
 Free(b);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);

 Free(p0x);
 Free(Pnx);
 Free(zx);
 Free(wx);
 Free(Sx);

 return 0;
}


double xse2_Wq(double lx, double ls, double cx, double csl, double csu, double p, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm)
{ double *p0, *Sx, *Pnx, *wx, *zx, *p0x, *S1s, *S2s, *Pns, *ws, *zs, *p0s, q, *zch, *rside, *b, za=0., s2, dN, Hij,
         mn_minus=1., mn_plus=0., mn_minusx, mn_minuss, mn_plusx, mn_pluss, ddf, xl, xu, oben, unten, q_minus=0., q_plus=0., enumerator=0., Wq=0.;
  int i, j, k, n, *ps, Ntilde, ihs, M, NN, ii, jj;

 cx  *= sqrt( lx/(2.-lx) ); 
 hsx *= sqrt( lx/(2.-lx) );

 M = ceil( (log(csl)-log(csu))/log(1.-ls) );
 Ntilde = ceil( (double)Ns/(double)M );
 NN = M*Ntilde;
 s2 = sigma*sigma;
 ddf = (double)df;
 dN = (double)Ntilde;

 ihs = floor( (log(csl)-log(hss))/log(1.-ls) );
 if (ihs<0) ihs = 0;

 Sx  = matrix(Nx,Nx);
 wx  = vector(Nx);
 zx  = vector(Nx);
 Pnx = matrix(nmax,Nx);
 p0x = vector(nmax); 

 S1s = matrix(NN,NN);
 S2s = matrix(NN,NN);
 ps = ivector(NN);
 zch = matrix(M,Ntilde);
 rside = vector(NN);
 b   = vector(M+1);
 ws  = vector(qm);
 zs  = vector(qm);
 Pns = matrix(nmax,NN);
 p0s = vector(nmax);

 p0 = vector(nmax);
 
 gausslegendre(Nx,-cx,cx,zx,wx);

 for (i=0;i<Nx;i++) {
   za = (1.-lx)*zx[i];
   for (j=0;j<Nx;j++)
     Sx[i*Nx+j] = wx[j]/lx*phi( ((zx[j]-za)/lx-mu)/sigma, 0.)/sigma;
 }  

/* interval borders b_i = cl/(1-l)^i */
 for (i=0;i<M;i++) b[i] = csl/pow(1.-ls, (double)(i));
 b[M] = csu;

 /* Chebyshev nodes on [b_0,b_1],[b_1,b_2],...,[b_M-1,cu] */
 for (i=0;i<M;i++)
   for (j=0;j<Ntilde;j++) {
     zch[ i*Ntilde+j ] = b[i] + (b[i+1]-b[i])/2.*(1.+cos(PI*(2.*j+1.)/2./dN));
   }

 /* P(L>1)(zch[i,j]) */
 for (i=0;i<M;i++)
   for (j=0;j<Ntilde;j++) {
     rside[ i*Ntilde+j ] =  CHI( ddf/s2*(csu-(1.-ls)*zch[ i*Ntilde+j ])/ls, df)
                          - CHI( ddf/s2*(csl-(1.-ls)*zch[ i*Ntilde+j ])/ls, df);
   }

 for (i=0;i<M;i++)
   for (j=0;j<Ntilde;j++) {
     za = (1.-ls)*zch[ i*Ntilde+j ];
     for (ii=0;ii<M;ii++)
       for (jj=0;jj<Ntilde;jj++) {
         if (b[ii+1]<za) S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = 0.;
         else {
           if (za<b[ii]) xl = b[ii]; else xl = za;
           xu = b[ii+1];
           if (df!=2) {
             xl = sqrt(xl-za);
             xu = sqrt(xu-za);
           }
           gausslegendre(qm,xl,xu,zs,ws);
           Hij = 0.;
           for (k=0;k<qm;k++)
             if (df==2)
               Hij += ws[k]*Tn( (2.*zs[k]-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                      * exp(-zs[k]/s2/ls);
             else
               Hij +=
            ws[k]*Tn( (2.*(zs[k]*zs[k]+za)-b[ii+1]-b[ii])/(b[ii+1]-b[ii]), jj)
                 * 2. * pow(zs[k], ddf-1.) * exp(-ddf*zs[k]*zs[k]/2./s2/ls);

           if (df==2) Hij *= exp(za/s2/ls)/s2/ls;
           else       Hij /= gammafn(ddf/2.) * pow(2.*s2*ls/ddf, ddf/2.);
           S1s[ (i*Ntilde+j)*NN + ii*Ntilde+jj ] = Hij;
         }
       }
   }

 for (i=0;i<NN;i++)
   for (j=0;j<NN;j++) S2s[i*NN+j] = 0.;

 for (i=0;i<M;i++)
   for (j=0;j<Ntilde;j++)
     for (jj=0;jj<Ntilde;jj++)
       S2s[ (i*Ntilde+j)*NN + i*Ntilde+jj ] =
         Tn( (2.*zch[ i*Ntilde+j ]-b[i+1]-b[i])/(b[i+1]-b[i]), jj); 

 LU_decompose(S2s,ps,NN);

 for (n=1;n<=nmax;n++) {

   if (n==1)
     for (i=0;i<Nx;i++)
       Pnx[i] = PHI( (( cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.) - 
                PHI( ((-cx-(1.-lx)*zx[i])/lx-mu)/sigma, 0.);
   else
     for (i=0;i<Nx;i++) {
       Pnx[(n-1)*Nx+i] = 0.;
       for (j=0;j<Nx;j++)
         Pnx[(n-1)*Nx+i] += Sx[i*Nx+j] * Pnx[(n-2)*Nx+j];
     }

   p0x[n-1] = 0.;
   if (n==1)
     p0x[0] = PHI( (( cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.) - 
              PHI( ((-cx-(1.-lx)*hsx)/lx-mu)/sigma, 0.);
   else
     for (j=0;j<Nx;j++)
       p0x[n-1] += wx[j]/lx * phi( ((zx[j]-(1.-lx)*hsx)/lx-mu)/sigma, 0.)/sigma 
                   * Pnx[(n-2)*Nx+j];

   if (n==1)
     for (i=0;i<M;i++)
       for (j=0;j<Ntilde;j++) {
         Pns[ i*Ntilde+j ] = 0.;
         for (jj=0;jj<Ntilde;jj++)
           Pns[ i*Ntilde+j ] += /*  usual Chebyshev approximation  */
             2./Ntilde * Tn( (2.*zch[i*Ntilde+jj]-b[i+1]-b[i])/(b[i+1]-b[i]), j)
                  * rside[ i*Ntilde+jj ];
         if (j==0) Pns[ i*Ntilde+j ] /= 2.;
       }
   else {
     for (i=0;i<NN;i++) {
       rside[i] = 0.;
       for (j=0;j<NN;j++) rside[i] += S1s[ i*NN+j ] * Pns[ (n-2)*NN+j ];
     }
     LU_solve2(S2s,rside,ps,NN);
     for (i=0;i<NN;i++) Pns[ (n-1)*NN+i ] = rside[i];
   }

   p0s[n-1] = 0.;
   if (n==1)
     p0s[0] =  CHI( ddf/s2*(csu-(1.-ls)*hss)/ls, df)
             - CHI( ddf/s2*(csl-(1.-ls)*hss)/ls, df);
   else
     for (j=0;j<Ntilde;j++)
       p0s[n-1] += Pns[ (n-1)*NN + ihs*Ntilde+j ]
                * Tn( (2.*hss-b[ihs+1]-b[ihs])/(b[ihs+1]-b[ihs]), j);


   p0[n-1] = p0x[n-1] * p0s[n-1];

   if ( p0[n-1] < 1.-p ) {
     Wq = (double)n;
     n = nmax+1;
   } else { 
     mn_minusx = 1.; mn_plusx = 0.;
     mn_minuss = 1.; mn_pluss = 0.;
     if ( n>1 ) {
       for (i=0;i<Nx;i++) {
         if (Pnx[(n-1)*Nx+i]==0)
           if (Pnx[(n-1)*Nx+i]==0) q = 0.;
           else q = 1.;
         else q = Pnx[(n-1)*Nx+i]/Pnx[(n-2)*Nx+i];
        if ( q<mn_minusx ) mn_minusx = q;
        if ( q>mn_plusx ) mn_plusx = q;
       }

       for (i=0;i<M;i++)
         for (j=0;j<Ntilde;j++) {
           oben = 0.;
           unten = 0.;
           for (jj=0;jj<Ntilde;jj++) {
             oben += Pns[ (n-1)*NN + i*Ntilde+jj ]
                   * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
             unten+= Pns[ (n-2)*NN + i*Ntilde+jj ]
                   * Tn((2.*zch[i*Ntilde+j]-b[i+1]-b[i])/(b[i+1]-b[i]), jj);
           }
           if (fabs(unten)<1e-16)
             if (fabs(oben)<1e-16) q = 0.;
             else q = 1.;
           else q = oben/unten;
           if ( q<mn_minuss ) mn_minuss = q;
           if ( q>mn_pluss ) mn_pluss = q;
         }

       mn_minus = mn_minusx * mn_minuss;
       mn_plus  = mn_plusx * mn_pluss;

       enumerator = log( (1.-p)/p0[n-1] );
       q_minus = (double)n + enumerator/log(mn_minus);
       q_plus  = (double)n + enumerator/log(mn_plus);
       if ( fabs( ceil(q_plus) - ceil(q_minus) ) < .5 ) {
	 Wq = ceil(q_plus);
	 n = nmax +1;
       }   
     } /* n > 1 */
   } /* p0[n-1] >= 1.-p */
 } /* n=1; n<=nmax; n++ */

 Free(p0);
 
 Free(p0s);
 Free(Pns);
 Free(zs);
 Free(ws);
 Free(b);
 Free(rside);
 Free(zch);
 Free(ps);
 Free(S2s);
 Free(S1s);

 Free(p0x);
 Free(Pnx);
 Free(zx);
 Free(wx);
 Free(Sx);

 return Wq;
}


int xse2lu_crit(double lx, double ls, double L0, double *cx, double csl, double *csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm)
{ double x1, x2, dx, s1, s2, ds, xARL1, xARL2, sARL1, sARL2, xsARL22, xsARL12, xsARL21,
         f11, f22, f21, f12, d11, d22, d21, d12, nenner, zr=0, c0=-1.;

 x1 = xe_crit(ewma2,lx,2.*L0,zr,hsx,mu,fix,Nx,c0) - .1;
 x2 = x1 + .2;
 s1 = se2lu_crit(ls,2.*L0,csl,hss,sigma,df,Ns,qm) - .1;
 s2 = s1 + .2;

 xARL2 = xe2_iglarl(lx,x2,hsx,mu,Nx);
 sARL2 = se2_iglarl(ls,csl,s2,hss,sigma,df,Ns,qm);
 xsARL22 = xse2_arl(lx,ls,x2,csl,s2,hsx,hss,mu,sigma,df,Nx,Ns,nmax,qm);

 do {
   xARL1 = xe2_iglarl(lx,x1,hsx,mu,Nx);
   sARL1 = se2_iglarl(ls,csl,s1,hss,sigma,df,Ns,qm);
   xsARL21 = xse2_arl(lx,ls,x2,csl,s1,hsx,hss,mu,sigma,df,Nx,Ns,nmax,qm);
   xsARL12 = xse2_arl(lx,ls,x1,csl,s2,hsx,hss,mu,sigma,df,Nx,Ns,nmax,qm);

   /* difference quotient */
   f11 = (xsARL22 - xsARL12)/(x2-x1); f12 = (xsARL22 - xsARL21)/(s2-s1);
   f21 = (xARL2   -   xARL1)/(x2-x1); f22 = (sARL1   -   sARL2)/(s2-s1);

   /* inverse of the difference quotient */
   nenner = f11*f22 - f12*f21;
   d11 =  f22/nenner;  d12 = -f12/nenner;
   d21 = -f21/nenner;  d22 =  f11/nenner;

   dx = d11*(xsARL22-L0) + d12*(xARL2-sARL2);
   ds = d21*(xsARL22-L0) + d22*(xARL2-sARL2);

   x1 = x2;   s1 = s2;
   x2 -= dx;  s2 -= ds;

   xARL2 = xe2_iglarl(lx,x2,hsx,mu,Nx);
   sARL2 = se2_iglarl(ls,csl,s2,hss,sigma,df,Ns,qm);
   xsARL22 = xse2_arl(lx,ls,x2,csl,s2,hsx,hss,mu,sigma,df,Nx,Ns,nmax,qm);

 } while (  (fabs(L0-xsARL22)>1e-6 || fabs(xARL2-sARL2)>1e-6) && (fabs(x2-x1)>1e-7 || fabs(s2-s1)>1e-7)   );

 *cx = x2; *csu = s2;

 return 0;
}


int xse2fu_crit(double lx, double ls, double L0, double *cx, double *csl, double csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm)
{ double x1, x2, dx, s1, s2, ds, xARL1, xARL2, sARL1, sARL2, xsARL22, xsARL12, xsARL21, 
         f11, f22, f21, f12, d11, d22, d21, d12, nenner, zr=0, c0=-1.;

 x1 = xe_crit(ewma2,lx,2.*L0,zr,hsx,mu,fix,Nx,c0) - .1;
 x2 = x1 + .2;
 s1 = se2fu_crit(ls,2.*L0,csu,hss,sigma,df,Ns,qm) - .1;
 s2 = s1 + .2;

 xARL2 = xe2_iglarl(lx,x2,hsx,mu,Nx);
 sARL2 = se2_iglarl(ls,s2,csu,hss,sigma,df,Ns,qm);
 xsARL22 = xse2_arl(lx,ls,x2,s2,csu,hsx,hss,mu,sigma,df,Nx,Ns,nmax,qm);
 /*printf("cx = %.4f,\tcsk = %.4f,\tcsu = %.4f\t,\txARL = %.2f,\tsARL = %.2f,\txsARL = %.2f\n", x2, s2, csu, xARL2, sARL2, xsARL22);*/
 do {
   xARL1 = xe2_iglarl(lx,x1,hsx,mu,Nx);
   sARL1 = se2_iglarl(ls,s1,csu,hss,sigma,df,Ns,qm);
   xsARL21 = xse2_arl(lx,ls,x2,s1,csu,hsx,hss,mu,sigma,df,Nx,Ns,nmax,qm);
   xsARL12 = xse2_arl(lx,ls,x1,s2,csu,hsx,hss,mu,sigma,df,Nx,Ns,nmax,qm);

   /* difference quotient */
   f11 = (xsARL22 - xsARL12)/(x2-x1); f12 = (xsARL22 - xsARL21)/(s2-s1);
   f21 = (xARL2   -   xARL1)/(x2-x1); f22 = (sARL1   -   sARL2)/(s2-s1);

   /* inverse of the difference quotient */
   nenner = f11*f22 - f12*f21;
   d11 =  f22/nenner;  d12 = -f12/nenner;
   d21 = -f21/nenner;  d22 =  f11/nenner;

   dx = d11*(xsARL22-L0) + d12*(xARL2-sARL2);
   ds = d21*(xsARL22-L0) + d22*(xARL2-sARL2);

   x1 = x2;   s1 = s2;
   x2 -= dx;  s2 -= ds;

   xARL2 = xe2_iglarl(lx,x2,hsx,mu,Nx);
   sARL2 = se2_iglarl(ls,s2,csu,hss,sigma,df,Ns,qm);
   xsARL22 = xse2_arl(lx,ls,x2,s2,csu,hsx,hss,mu,sigma,df,Nx,Ns,nmax,qm);
   /*printf("cx = %.4f,\tcsk = %.4f,\tcsu = %.4f\t,\txARL = %.2f,\tsARL = %.2f,\txsARL = %.2f\n", x2, s2, csu, xARL2, sARL2, xsARL22);*/
 } while (  (fabs(L0-xsARL22)>1e-6 || fabs(xARL2-sARL2)>1e-6) && (fabs(x2-x1)>1e-8 || fabs(s2-s1)>1e-8)   );

 *cx = x2; *csl = s2;

 return 0;
}


int xse2_crit(double lx, double ls, double L0, double *cx, double *csl, double *csu, double hsx, double hss, double mu, double sigma, int df, int Nx, int Ns, int nmax, int qm)
{ double s1, s2, s3, ds, sl1, sl2, sl3, Lm, Lp, x, cl;
  int flag;

 cl = 0.;
 flag = xseU_crit(lx,ls,L0,&x,&s1,hsx,hss,mu,sigma,df,Nx,Ns,nmax,qm);
 /*printf("cx = %.4f,\tcsl = %.4f,\tcsu = %.4f\n", x, cl, s1);*/

 Lm = xseU_arl(lx,ls,x,s1,hsx,hss,mu,sigma-lmEPS,df,Nx,Ns,nmax,qm);
 Lp = xseU_arl(lx,ls,x,s1,hsx,hss,mu,sigma+lmEPS,df,Nx,Ns,nmax,qm);
 sl1 = (Lp-Lm)/(2.*lmEPS);
 s2 = s1 + .15;
 flag = xse2fu_crit(lx,ls,L0,&x,&cl,s2,hsx,hss,mu,sigma,df,Nx,Ns,nmax,qm);
 /*printf("cx = %.4f,\tcsl = %.4f,\tcsu = %.4f\n", x, cl, s2);*/
 Lm = xse2_arl(lx,ls,x,cl,s2,hsx,hss,mu,sigma-lmEPS,df,Nx,Ns,nmax,qm);
 Lp = xse2_arl(lx,ls,x,cl,s2,hsx,hss,mu,sigma+lmEPS,df,Nx,Ns,nmax,qm);
 sl2 = (Lp-Lm)/(2.*lmEPS);

 do {
   s3 = s1 - sl1/(sl2-sl1) * (s2-s1);
   flag = xse2fu_crit(lx,ls,L0,&x,&cl,s3,hsx,hss,mu,sigma,df,Nx,Ns,nmax,qm);
   Lm = xse2_arl(lx,ls,x,cl,s3,hsx,hss,mu,sigma-lmEPS,df,Nx,Ns,nmax,qm);
   Lp = xse2_arl(lx,ls,x,cl,s3,hsx,hss,mu,sigma+lmEPS,df,Nx,Ns,nmax,qm);
   sl3 = (Lp-Lm)/(2.*lmEPS);
   /*printf("cx = %.4f,\tcsl = %.4f,\tcsu = %.4f\t,\tslope = %.6f\n", x, cl, s3, sl3);*/
   ds = s3-s2; s1 = s2; sl1 = sl2; s2 = s3; sl2 = sl3;
 } while ( fabs(sl3)>1e-6 && fabs(ds)>1e-7 );

 *cx = x; *csl = cl; *csu = s3;

 return flag;
}


/* EWMA p under sampling by variables */

/* p = h(mu, sigma) */

double WK_h(double mu, double sigma, double LSL, double USL)
{ double result;
 result = PHI( (LSL-mu)/sigma, 0.) + PHI( (mu-USL)/sigma, 0.);
 return result;
}


/* d/dmu h(mu, sigma) */

double wk_h_mu(double mu, double sigma, double LSL, double USL)
{ double result;
 result = ( -phi( (LSL-mu)/sigma, 0.) + phi( (mu-USL)/sigma, 0.) )/sigma;
 return result;
}


/* d/dsigma h(mu, sigma) */

double wk_h_sigma(double mu, double sigma, double LSL, double USL)
{ double result;
 result = -( (LSL-mu)*phi( (LSL-mu)/sigma, 0.) + (mu-USL)*phi( (mu-USL)/sigma, 0.) )/sigma/sigma;
 return result;
}


/* mu = h^-1(p, sigma) */

double WK_h_invers_mu(double p, double sigma, double LSL, double USL)
{ double mu, old_mu, merror, perror;
 mu = sigma*qPHI(p) + USL;
 perror = WK_h(mu, sigma, LSL, USL) - p;
 do {
   old_mu = mu;
   mu = mu - perror / wk_h_mu(mu, sigma, LSL, USL);
   merror = mu - old_mu;
   perror = WK_h(mu, sigma, LSL, USL) - p;
 } while ( fabs(merror) > 1e-10 && fabs(perror) > 1e-12 );
 return mu;
}


/* sigma = h^-1(p, mu) */

double WK_h_invers_sigma(double p, double mu, double LSL, double USL)
{ double sigma, old_sigma, serror, perror;
 sigma = (mu-USL)/qPHI(p);
 perror = WK_h(mu, sigma, LSL, USL) - p;
 do {
   old_sigma = sigma;
   sigma = sigma - perror / wk_h_sigma(mu, sigma, LSL, USL);
   serror = sigma - old_sigma;
   perror = WK_h(mu, sigma, LSL, USL) - p;
 } while ( fabs(serror) > 1e-10 && fabs(perror) > 1e-12 );
 return sigma;
}


/* alpha, the upper limit of the cdf (and pdf) definite integral */
double wk_alpha(double p, double sigma, int n, double LSL, double USL)
{ double alpha, dn, zphalf;
 dn = (double)n;
 zphalf = qPHI(p/2.);
 alpha = (dn-1.)/sigma/sigma * (USL-LSL)*(USL-LSL)/4. / (zphalf*zphalf);
 return alpha;
}


/* cdf of h(xbar, sigma0=1) for X ~ N(mu, sigma) */

double cdf_phat(double p, double mu, double sigma, int n, double LSL, double USL)
{ double result, pstar, mu_of_p, dn, centre;
 dn = (double)n;
 result = 0.;
 if ( p >= 1. ) result = 1.;
 centre = (LSL+USL)/2.;
 /*pstar = WK_h(centre, sigma, LSL, USL);*/
 pstar = WK_h(centre, 1., LSL, USL);
 if ( pstar < p && p < 1. ) {
   mu_of_p = WK_h_invers_mu(p, 1., LSL, USL);
   result = PHI( (mu_of_p - mu)*sqrt(dn)/sigma, 0. ) - PHI( (-mu_of_p - mu)*sqrt(dn)/sigma, 0. );
 }
 return result;
}


/* pdf of h(xbar, sigma0=1) for X ~ N(mu, sigma) */

double pdf_phat(double p, double mu, double sigma, int n, double LSL, double USL)
{ double result, pstar, mu_of_p, dn, centre;
 dn = (double)n;
 result = 0.;
 centre = (LSL+USL)/2.;
 /*pstar = WK_h(centre, sigma, LSL, USL);*/
 pstar = WK_h(centre, 1., LSL, USL);
 if ( pstar < p && p < 1. ) {
   mu_of_p = WK_h_invers_mu(p, 1., LSL, USL);
   result = sqrt(dn)*( phi( (mu_of_p - mu)*sqrt(dn)/sigma, 0. ) + phi( (-mu_of_p - mu)*sqrt(dn)/sigma, 0. ) ) / wk_h_mu(mu_of_p, 1., LSL, USL)/sigma;
 }
 return result;
}


/* quantile function of h(xbar, sigma0=1) for X ~ N(mu, sigma) */

double qf_phat(double p0, double mu, double sigma, int n, double LSL, double USL)
{ double pstar, centre, c1, c2, c3, p1, p2, p3, dc, cstep;
 centre = (LSL+USL)/2.;
 pstar = WK_h(centre, sigma, LSL, USL);
 c2 = pstar;
 p2 = 0.;
 cstep = p0/1e3;
 do {
   c1 = c2;
   p1 = p2;
   c2 += cstep;
   p2 = cdf_phat(c2, mu, sigma, n, LSL, USL);
 } while ( p2 < p0 );
 if ( c2 <= pstar + cstep + 1e-9 ) {
   c1 = c2 - cstep/2.;
   p1 = cdf_phat(c1, mu, sigma, n, LSL, USL);
 }
 do {
   c3 = c1 + ( p0 - p1 )/( p2 - p1 ) * ( c2 - c1 );
   p3 = cdf_phat(c3, mu, sigma, n, LSL, USL);
   dc = c3 - c2; c1 = c2; p1 = p2; c2 = c3; p2 = p3;
 } while ( fabs( p0 - p3 )>1e-10 && fabs(dc)>1e-10 );
 return c3; 
}


/* integrand for cdf of h(xbar, s) for X ~ N(mu, sigma) */

double wk_cdf_i(double y, double p, double mu, double sigma, int n, double LSL, double USL)
{ double result, alpha, x, s, mu_p, dn, atrim;
  dn = (double)n;
  alpha = wk_alpha(p, sigma, n, LSL, USL);
  atrim = qCHI(0.9999999999, n-1);
  if ( atrim < alpha ) alpha = atrim;
  x = alpha - pow(y,2.);
  s = sigma * sqrt( x/(dn-1.) );
  mu_p = WK_h_invers_mu(p, s, LSL, USL);
  result = PHI( (mu_p-mu)*sqrt(dn)/sigma, 0.) - PHI( (-mu_p-mu)*sqrt(dn)/sigma, 0.);
  result *= chi(x, n-1) * 2*y;
 return result;
}

/* cdf of h(xbar, s) for X ~ N(mu, sigma) */

double cdf_phat2(double p, double mu, double sigma, int n, double LSL, double USL, int nodes)
{ double result, alpha, *w, *z, xl, xu, atrim;
  int i;  
 w = vector(nodes);
 z = vector(nodes); 
 result = 0.;
 if ( p >= 1. ) result = 1.;
 xl = 0.;
 if ( 0. < p && p < 1. ) {
   alpha = wk_alpha(p, sigma, n, LSL, USL);
   atrim = qCHI(0.9999999999, n-1);
   if ( atrim < alpha ) alpha = atrim;
   xu = pow(alpha,0.5);
   gausslegendre(nodes, xl, xu, z, w);
   for (i=0; i<nodes; i++) result += w[i] * wk_cdf_i(z[i], p, mu, sigma, n, LSL, USL);   
 } 
 Free(z);
 Free(w); 
 return result;
}

/* integrand for pdf of h(xbar, s) for X ~ N(mu, sigma) */

double wk_pdf_i(double y, double p, double mu, double sigma, int n, double LSL, double USL)
{ double result, alpha, x, s, mu_p, dn;
  dn = (double)n;
  alpha = wk_alpha(p, sigma, n, LSL, USL);
  x = alpha - y*y;
  s = sigma * sqrt( x/(dn-1.) );
  mu_p = WK_h_invers_mu(p, s, LSL, USL);  
  result = ( phi( (mu_p-mu)*sqrt(dn)/sigma, 0.) + phi( (-mu_p-mu)*sqrt(dn)/sigma, 0.) ) * sqrt(dn)/sigma;
  result /= wk_h_mu(mu_p, s, LSL, USL);
  result *= chi(x, n-1) * 2.*y;  
  return result;
}

/* pdf of h(xbar, s) for X ~ N(mu, sigma) */

double pdf_phat2(double p, double mu, double sigma, int n, double LSL, double USL, int nodes)
{ double result, alpha, *w, *z, xl, xu;
  int i;  
 w = vector(nodes);
 z = vector(nodes); 
 result = 0.;
 xl = 0.;
 if ( 0. < p && p < 1. ) {
   alpha = wk_alpha(p, sigma, n, LSL, USL);
   xu = sqrt(alpha);
   gausslegendre(nodes, xl, xu, z, w);
   for (i=0; i<nodes; i++) result += w[i] * wk_pdf_i(z[i], p, mu, sigma, n, LSL, USL);
 } 
 Free(z);
 Free(w); 
 return result;
}

/* quantile function of h(xbar, s) for X ~ N(mu, sigma) */

double qf_phat2(double p0, double mu, double sigma, int n, double LSL, double USL, int nodes)
{ double c1, c2, c3, p1, p2, p3, dc, cstep;
 c2 = 0.;
 p2 = 0.;
 cstep = p0/1e3;
 do {
   c1 = c2;
   p1 = p2;
   c2 += cstep;
   p2 = cdf_phat2(c2, mu, sigma, n, LSL, USL, nodes);
 } while ( p2 < p0 );
 if ( c2 <= cstep + 1e-9 ) {
   c1 = c2 - cstep/2.;
   p1 = cdf_phat2(c1, mu, sigma, n, LSL, USL, nodes);
 }
 do {
   c3 = c1 + ( p0 - p1 )/( p2 - p1 ) * ( c2 - c1 );
   p3 = cdf_phat2(c3, mu, sigma, n, LSL, USL, nodes);
   dc = c3 - c2; c1 = c2; p1 = p2; c2 = c3; p2 = p3;
 } while ( fabs( p0 - p3 )>1e-10 && fabs(dc)>1e-10 );
 return c3; 
}


/* collocation */

double ewma_phat_arl(double lambda, double ucl, double mu, double sigma, int n, double z0, double LSL, double USL, int N, int qm)
{ double *a, *g, *w, *z, arl, Hij, dN, xl, xu, za, ll, pstar, xi, centre;
  int i, j, k;

 dN = (double)N;
 a = matrix(N,N);
 g = vector(N);
 w = vector(qm);
 z = vector(qm);
 
 centre = (LSL+USL)/2.;
 /*pstar = WK_h(centre, sigma, LSL, USL);*/
 pstar = WK_h(centre, 1., LSL, USL);

 for (i=0; i<N; i++) {
   xi = pstar + (ucl - pstar)/2. * (1.+cos(PI*(2.*(i+1.)-1.)/2./dN));
   za = (1.-lambda)*xi;
   ll = za + lambda*pstar;
   xl = 0.; 
   xu = sqrt(ucl - ll);
   gausslegendre(qm, xl, xu, z, w);
   a[i*N] = 1. - cdf_phat( (ucl - za)/lambda, mu, sigma, n, LSL, USL);
   for (j=1; j<N; j++) {
     Hij = 0.;
     for (k=0; k<qm; k++) {
       Hij += w[k] * Tn( 2.*(z[k]*z[k] + ll - pstar)/(ucl - pstar) - 1. ,j) * 2.*z[k]*pdf_phat(z[k]*z[k]/lambda + pstar, mu, sigma, n, LSL, USL)/lambda;
     }
     a[i*N+j] = Tn( 2.*(xi - pstar)/(ucl - pstar) - 1., j) - Hij;
   }
 }

 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);
 arl = g[0];
 for (j=1;j<N;j++) arl += g[j] * Tn( 2.*(z0 - pstar)/(ucl - pstar)-1., j);

 Free(z);
 Free(w);
 Free(g);
 Free(a);

 return arl;
}


double ewma_phat_arl_be(double lambda, double ucl, double mu, double sigma, int n, double z0, double LSL, double USL, int N)
{ double *a, *g, w, arl, dN, pstar, centre;
  int i, j;
  
 dN = (double)N;
 
 a = matrix(N,N);
 g = vector(N);
 
 centre = (LSL+USL)/2.;
 /*pstar = WK_h(centre, sigma, LSL, USL);*/
 pstar = WK_h(centre, 1., LSL, USL);
 
 w = (ucl - pstar)/dN; 
 
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++)
     a[i*N+j] = - ( cdf_phat( pstar + ((j+1)*w-(1.-lambda)*(i+0.5)*w)/lambda, mu, sigma, n, LSL, USL) - 
                    cdf_phat( pstar + (    j*w-(1.-lambda)*(i+0.5)*w)/lambda, mu, sigma, n, LSL, USL) );
   ++a[i*N+i];
 }

 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);

 arl = 1.;
 for (j=0;j<N;j++)
   arl += ( cdf_phat( (pstar+(j+1)*w-(1.-lambda)*z0)/lambda, mu, sigma, n, LSL, USL) -
            cdf_phat( (pstar+    j*w-(1.-lambda)*z0)/lambda, mu, sigma, n, LSL, USL) ) * g[j];

 Free(g);
 Free(a);

 return arl;
}


/* collocation */

double ewma_phat_arl2(double lambda, double ucl, double mu, double sigma, int n, double z0, double LSL, double USL, int N, int qm, int M)
{ double *a, *g, *w, *z, arl, Hij, dN, xl, xu, za, xi, dM, x, FF;
  int i, j, k, nodes=30;
  
 dN = (double)N;
 dM = (double)M;
 
 a = matrix(N,N);
 g = vector(N);
 w = vector(qm);
 z = vector(qm);

 for (i=0; i<N; i++) {
   xi = ucl/2. * (1.+cos(PI*(2.*(i+1.)-1.)/2./dN));
   za = (1.-lambda)*xi;
   FF = cdf_phat2( (ucl-za)/lambda, mu, sigma, n, LSL, USL, nodes);
   a[i*N] = 1. - FF;   
   xl = 0.; 
   xu = pow(ucl - za, 1./dM);
   gausslegendre(qm, xl, xu, z, w);   
   for (j=1; j<N; j++) {
     Hij = 0.;
     for (k=0; k<qm; k++) {
       x = pow(z[k],dM) + za;
       Hij += w[k] * dTn( 2.*x/ucl-1. ,j)*2./ucl * cdf_phat2( (x-za)/lambda, mu, sigma, n, LSL, USL, nodes) * dM*pow(z[k],dM-1.);
     }
     a[i*N+j] = Tn( 2.*xi/ucl - 1., j) - (FF - Hij);
   }
 }

 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);
 arl = g[0];
 for (j=1;j<N;j++) arl += g[j] * Tn( 2.*z0/ucl-1., j);

 Free(z);
 Free(w);
 Free(g);
 Free(a);

 return arl;
}


double ewma_phat_arl2_be(double lambda, double ucl, double mu, double sigma, int n, double z0, double LSL, double USL, int N)
{ double *a, *g, w, arl, dN;
  int i, j, nodes=30;
  
 dN = (double)N;
 w = ucl/dN;
 
 a = matrix(N,N);
 g = vector(N);
 
 for (i=0; i<N; i++) {
   for (j=0; j<N; j++)
     a[i*N+j] = - ( cdf_phat2( ((j+1)*w-(1.-lambda)*(i+0.5)*w)/lambda, mu, sigma, n, LSL, USL, nodes) - 
                    cdf_phat2( (    j*w-(1.-lambda)*(i+0.5)*w)/lambda, mu, sigma, n, LSL, USL, nodes) );
   ++a[i*N+i];
 }

 for (j=0; j<N; j++) g[j] = 1.;
 LU_solve(a, g, N);

 arl = 1.;
 for (j=0;j<N;j++)
   arl += ( cdf_phat2( ((j+1)*w-(1.-lambda)*z0)/lambda, mu, sigma, n, LSL, USL, nodes) -
            cdf_phat2( (    j*w-(1.-lambda)*z0)/lambda, mu, sigma, n, LSL, USL, nodes) ) * g[j];

 Free(g);
 Free(a);

 return arl;
}


double ewma_phat_crit(double lambda, double L0, double mu, double sigma, int n, double z0, double LSL, double USL, int N, int qm)
{ double c1, c2, c3, L1, L2, L3, dc, pstar, cstep, centre;
 centre = (LSL+USL)/2.;
 pstar = WK_h(centre, sigma, LSL, USL);
 c2 = pstar;
 cstep = lambda/10.;
 do {
   c2 += cstep;
   L2 = ewma_phat_arl(lambda, c2, mu, sigma, n, z0, LSL, USL, N, qm);
 } while ( L2 < L0 );
 c1 = c2 - cstep;
 if ( c2 <= pstar + cstep + 1e-9 ) c1 = c2 - cstep/2.;
 L1 = ewma_phat_arl(lambda, c1, mu, sigma, n, z0, LSL, USL, N, qm);
 do {
   c3 = c1 + ( L0 - L1 )/( L2 - L1 ) * ( c2 - c1 );
   L3 = ewma_phat_arl(lambda, c3, mu, sigma, n, z0, LSL, USL, N, qm);
   dc = c3 - c2; c1 = c2; L1 = L2; c2 = c3; L2 = L3;
 } while ( fabs( L0 - L3 )>1e-6 && fabs(dc)>1e-12 );
 return c3;
}


double ewma_phat_crit2(double lambda, double L0, double mu, double sigma, int n, double z0, double LSL, double USL, int N, int qm, int M)
{ double c1, c2, c3, L1, L2, L3, dc, cstep;
 c2 = 0.;
 L2 = 0.;
 cstep = lambda/10.;
 do {
   c1 = c2;
   L1 = L2;
   c2 += cstep;
   L2 = ewma_phat_arl2(lambda, c2, mu, sigma, n, z0, LSL, USL, N, qm, M);
 } while ( L2 < L0 );
 if ( c2 <= cstep + 1e-9 ) {
   c1 = c2 - cstep/2.;
   L1 = ewma_phat_arl2(lambda, c1, mu, sigma, n, z0, LSL, USL, N, qm, M);
 }
 do {
   c3 = c1 + ( L0 - L1 )/( L2 - L1 ) * ( c2 - c1 );
   L3 = ewma_phat_arl2(lambda, c3, mu, sigma, n, z0, LSL, USL, N, qm, M);
   dc = c3 - c2; c1 = c2; L1 = L2; c2 = c3; L2 = L3;
 } while ( fabs( L0 - L3 )>1e-6 && fabs(dc)>1e-12 );
 return c3;
}


int N_of_l(double lambda)
{ int N;
 N = 20;
 if ( lambda < 1e-1 ) N = 40;
 if ( lambda < 1e-2 ) N = 60;
 if ( lambda < 1e-3 ) N = 120;
 if ( lambda < 1e-4 ) N = 200;
 return N;
}


double ewma_phat_lambda(double L0, double mu, double sigma, double max_l, double min_l, int n, double z0, double LSL, double USL, int qm)
{ double dn, cS, cE, ldelta, one, L1, L1_, lambda;
  int i, j, N;
 lambda = 1.;
 dn = (double)n;
 cS = qPHI( 1. - 1./(2.*L0) )/sqrt(dn)*sigma;
 cE = WK_h( cS, 1., LSL, USL );
 L1 = 1./( PHI( (-cS-mu)*sqrt(dn)/sigma, 0.) + 1. - PHI( (cS-mu)*sqrt(dn)/sigma, 0.) );
 ldelta = .1;
 one = 1; 
 for (j=0; j<4; j++) {
   for (i=0; i<20; i++) {
     lambda = lambda - ldelta*one;
     if ( lambda <= min_l ) { lambda = min_l; i = 23; }
     if ( lambda >= max_l ) { lambda = max_l; i = 23; }      
     N = N_of_l(lambda);      
     cE  = ewma_phat_crit(lambda, L0, 0., sigma, n, z0, LSL, USL, N, qm);      
     L1_ = ewma_phat_arl(lambda, cE, mu, sigma, n, z0, LSL, USL, N, qm);
     if ( L1_ > L1 && i < 23 ) i = 21;
     L1 = L1_;
   }
   ldelta /= 10.;
   one *= -1.;
 }
 if ( i < 23 ) lambda -= 10.*ldelta*one;
 return lambda;
}


double ewma_phat_lambda2(double L0, double mu, double sigma, double max_l, double min_l, int n, double z0, double LSL, double USL, int qm, int M)
{ double dn, cS, cE, ldelta, one, L1, L1_, lambda;
  int i, j, N;
 lambda = 1.;
 dn = (double)n;
 cS = qPHI( 1. - 1./(2.*L0) )/sqrt(dn)*sigma;
 cE = WK_h( cS, 1., LSL, USL );
 L1 = 1./( PHI( (-cS-mu)*sqrt(dn)/sigma, 0.) + 1. - PHI( (cS-mu)*sqrt(dn)/sigma, 0.) );
 ldelta = .1;
 one = 1; 
 for (j=0; j<4; j++) {
   for (i=0; i<20; i++) {
     lambda = lambda - ldelta*one;
     if ( lambda <= min_l ) { lambda = min_l; i = 23; }
     if ( lambda >= max_l ) { lambda = max_l; i = 23; }      
     N = N_of_l(lambda);      
     cE  = ewma_phat_crit2(lambda, L0, 0., sigma, n, z0, LSL, USL, N, qm, M);      
     L1_ = ewma_phat_arl2 (lambda, cE, mu, sigma, n, z0, LSL, USL, N, qm, M);
     if ( L1_ > L1 && i < 23 ) i = 21;
     L1 = L1_;
   }
   ldelta /= 10.;
   one *= -1.;
 }
 if ( i < 23 ) lambda -= 10.*ldelta*one;
 return lambda;
}


/* attributive EWMA */
double ewma_pU_arl(double lambda, double ucl, int n, double p, double z0, int d_res, int round_mode, int mid_mode)
{ double *a, *g, arl, zj=0, pju, pj;
  int i, j, k, N, NN/*, k_max*/;
  
 N = (int)ceil(ucl*d_res);
 /*N = (int)floor(ucl*d_res);*/
 NN = N + 1;
 a = matrix(NN, NN);
 g = vector(NN);

 for (i=0; i<=N; i++) for (j=0; j<=N; j++) a[i*NN+j] = 0.;
 
 for (i=0; i<=N; i++) {
   /*k_max = (int)ceil( (ucl+1. - (1.-lambda)*i)/lambda );*/
   for (k=0; k<=n; k++) {
     zj = (1.-lambda)*i/d_res + lambda*k;  
     pj = pdf_binom((double)k, n, p);
     switch (round_mode) {
	case -1:	/* round down as probably Gan did */
	    j = (int)floor(zj*d_res + 1e-9);
	    if ( j <= N ) a[j*NN+i] += -pj;
	break;
	case 0:		/* round down */
	    j = (int)floor(zj*d_res);
	    if ( j <= N ) a[j*NN+i] += -pj;
	break;
	case 1:		/* round up */
	    j = (int)ceil(zj*d_res);
	    if ( j <= N ) a[j*NN+i] += -pj;
	break;
	case 2:		/* round to nearest -- round half to even, IEEE 754 */
	    j = (int)round(zj*d_res);
	    if ( j <= N ) a[j*NN+i] += -pj;
	break;
	case 3:		/* round to nearest -- round half up  */
	    j = (int)floor(zj*d_res+.5);
	    if ( j <= N ) a[j*NN+i] += -pj;
	break;
	case 4:		/* distribute */
	    j = (int)floor(zj*d_res);
	    pju = zj - j/d_res;          
	    if ( j <= N ) a[j*NN+i]   += -(1.-pju)*pj;
	    if ( j <  N ) a[(j+1)*NN+i] += -pju*pj;
	break;
     }
   }
   ++a[i*NN+i];
 } 

 for (j=0; j<=N; j++) g[j] = 1.;
 /*LU_solve(a, g, NN);*/
 solve(&NN, a, g);
 
 arl = 1.;
 /*k_max = (int)ceil( (ucl+1. - (1.-lambda)*z0)/lambda );*/
 for (k=0; k<=n; k++) {
   zj = (1.-lambda)*z0 + lambda*k;
   pj = pdf_binom((double)k, n, p);
   switch (round_mode) {
      case -1:	/* round down as probably Gan did */
	  j = (int)floor(zj*d_res + 1e-9);
	  if ( j <= N ) arl += pj*g[j];
      break;
      case 0:	/* round down */
	  j = (int)floor(zj*d_res);
	  if ( j <= N ) arl += pj*g[j];
      break;
      case 1:	/* round up */
	  j = (int)ceil(zj*d_res);
	  if ( j <= N ) arl += pj*g[j];
      break;
      case 2:	/* round to nearest -- round half to even, IEEE 754 */
	  j = (int)round(zj*d_res);
	  if ( j <= N ) arl += pj*g[j];
      break;
      case 3:	/* round to nearest -- round half up  */
	  j = (int)floor(zj*d_res+.5);
	  if ( j <= N ) arl += pj*g[j];
      break;
      case 4:	/* distribute */
	  j = (int)floor(zj*d_res);
	  pju = zj - j/d_res;          
	  if ( j <= N ) arl += (1.-pju)*pj*g[j];
	  if ( j <  N ) arl += pju*pj*g[j+1];
      break;
   }
 }     

 Free(a);
 Free(g);

 return arl;
}


double ewma_cU_arl(double lambda, double ucl, double mu, double z0, int d_res, int round_mode, int mid_mode)
{ double *a, *g, arl, zj=0, pju, pj;
  int i, j, k, N, NN, k_max;
  
 N = (int)ceil(ucl*d_res);
 /*N = (int)floor(ucl*d_res);*/
 NN = N + 1;
 a = matrix(NN, NN);
 g = vector(NN);

 for (i=0; i<=N; i++) for (j=0; j<=N; j++) a[i*NN+j] = 0.;
 
 for (i=0; i<=N; i++) {
   k_max = (int)ceil( (ucl+1. - (1.-lambda)*i)/lambda );
   for (k=0; k<=k_max; k++) {
     zj = (1.-lambda)*i/d_res + lambda*k;     
     pj = pdf_pois((double)k, mu);
     switch (round_mode) {
	case -1:	/* round down as probably Gan did */
	    j = (int)floor(zj*d_res + 1e-9);
	    if ( j <= N ) a[j*NN+i] += -pj;
	break;
	case 0:		/* round down */
	    j = (int)floor(zj*d_res);
	    if ( j <= N ) a[j*NN+i] += -pj;
	break;
	case 1:		/* round up */
	    j = (int)ceil(zj*d_res);
	    if ( j <= N ) a[j*NN+i] += -pj;
	break;
	case 2:		/* round to nearest -- round half to even, IEEE 754 */
	    j = (int)round(zj*d_res);
	    if ( j <= N ) a[j*NN+i] += -pj;
	break;
	case 3:		/* round to nearest -- round half up  */
	    j = (int)floor(zj*d_res+.5);
	    if ( j <= N ) a[j*NN+i] += -pj;
	break;
	case 4:		/* distribute */
	    j = (int)floor(zj*d_res);
	    pju = zj - j/d_res;          
	    if ( j <= N ) a[j*NN+i]   += -(1.-pju)*pj;
	    if ( j <  N ) a[(j+1)*NN+i] += -pju*pj;
	break;
     }
   }
   ++a[i*NN+i];
 } 

 for (j=0; j<=N; j++) g[j] = 1.;
 /*LU_solve(a, g, NN);*/
 solve(&NN, a, g);
 
 arl = 1.;
 k_max = (int)ceil( (ucl+1. - (1.-lambda)*z0)/lambda );
 for (k=0; k<=k_max; k++) {
   zj = (1.-lambda)*z0 + lambda*k;
   pj = pdf_pois((double)k, mu);
   switch (round_mode) {
      case -1:	/* round down as probably Gan did */
	  j = (int)floor(zj*d_res + 1e-9);
	  if ( j <= N ) arl += pj*g[j];
      break;
      case 0:	/* round down */
	  j = (int)floor(zj*d_res);
	  if ( j <= N ) arl += pj*g[j];
      break;
      case 1:	/* round up */
	  j = (int)ceil(zj*d_res);
	  if ( j <= N ) arl += pj*g[j];
      break;
      case 2:	/* round to nearest -- round half to even, IEEE 754 */
	  j = (int)round(zj*d_res);
	  if ( j <= N ) arl += pj*g[j];
      break;
      case 3:	/* round to nearest -- round half up  */
	  j = (int)floor(zj*d_res+.5);
	  if ( j <= N ) arl += pj*g[j];
      break;
      case 4:	/* distribute */
	  j = (int)floor(zj*d_res);
	  pju = zj - j/d_res;          
	  if ( j <= N ) arl += (1.-pju)*pj*g[j];
	  if ( j <  N ) arl += pju*pj*g[j+1];
      break;
   }
 }     

 Free(a);
 Free(g);

 return arl;
}


double ewma_pL_arl(double lambda, double lcl, int n, double p, double z0, int d_res, int round_mode, int mid_mode)
{ double *a, *g, arl, zj=0., pju, pj, zold, dr;
  int i, j, k, l, u, N, NN/*, k_max*/;
  
 dr = (double)d_res; 
 l = (int)floor(lcl*d_res);
 u = (int)qf_binom(.999999, n, p);
 N = u - l;
 NN = N + 1;
 a = matrix(NN, NN);
 g = vector(NN);

 for (i=0; i<=N; i++) for (j=0; j<=N; j++) a[i*NN+j] = 0.;
 
 for (i=0; i<=N; i++) {
   for (k=0; k<=n; k++) {
     zold = ( (double)( l + i ) ) / dr;
     zj = (1.-lambda)*zold + lambda*k;  
     pj = pdf_binom((double)k, n, p);
     switch (round_mode) {
        case -1:	/* round down as probably Gan did */
	    j = (int)floor(zj*dr + 1e-9) - l;
	    if ( 0 <= j && j <= u ) a[j*NN+i] += -pj;
        break;
        case 0:		/* round down */
	    j = (int)floor(zj*dr) - l;
	    if ( 0 <= j && j <= u ) a[j*NN+i] += -pj;
        break;
        case 1:		/* round up */
	    j = (int)ceil(zj*dr) - l;
	    if ( 0 <= j && j <= u ) a[j*NN+i] += -pj;
        break;
        case 2:		/* round to nearest -- round half to even, IEEE 754 */
	    j = (int)round(zj*dr) - l;
	    if ( 0 <= j && j <= u ) a[j*NN+i] += -pj;
        break;
        case 3:		/* round to nearest -- round half up  */
	    j = (int)floor(zj*dr+.5) - l;
	    if ( 0 <= j && j <= u ) a[j*NN+i] += -pj;
        break;
        case 4:		/* distribute */
	    j = (int)floor(zj*dr) - l;
	    pju = zj - j/dr;
            if ( 0 <= j && j <= u ) a[j*NN+i]     += -(1.-pju)*pj;
	    if ( 0 <  j && j <= u ) a[(j+1)*NN+i] += -pju*pj;
        break;
     }
   }
   ++a[i*NN+i];
 } 

 for (j=0; j<=N; j++) g[j] = 1.;
 solve(&NN, a, g);
 
 arl = 1.;
 for (k=0; k<=n; k++) {
   zj = (1.-lambda)*z0 + lambda*k;
   pj = pdf_binom((double)k, n, p);
   switch (round_mode) {
      case -1:	/* round down as probably Gan did */
	  j = (int)floor(zj*dr + 1e-9) - l;
	  if ( 0 <= j && j <= u ) arl += pj*g[j];
      break;
      case 0:	/* round down */
	  j = (int)floor(zj*dr) - l;
	  if ( 0 <= j && j <= u ) arl += pj*g[j];
      break;
      case 1:	/* round up */
	  j = (int)ceil(zj*dr) - l;
	  if ( 0 <= j && j <= u ) arl += pj*g[j];
      break;
      case 2:	/* round to nearest -- round half to even, IEEE 754 */
	  j = (int)round(zj*dr) - l;
	  if ( 0 <= j && j <= u ) arl += pj*g[j];
      break;
      case 3:	/* round to nearest -- round half up  */
	  j = (int)floor(zj*dr+.5) - l;
	  if ( 0 <= j && j <= u ) arl += pj*g[j];
      break;
      case 4:	/* distribute */
	  j = (int)floor(zj*dr) - l;
	  pju = zj - j/dr;          
	  if ( 0 <= j && j <= u ) arl += (1.-pju)*pj*g[j];
	  if ( 0 <  j && j <= u ) arl += pju*pj*g[j+1];
      break;
   }
 }     

 Free(a);
 Free(g);

 return arl;
}


double ewma_cL_arl(double lambda, double lcl, double ucl, double mu, double z0, int d_res, int round_mode, int mid_mode)
{ double *a, *g, arl, zj=0., pju, pj, zold, dr;
  int i, j, k, l, u, N, NN, k_max;
  
 dr = (double)d_res; 
 l = (int)floor(lcl*d_res);
 u = (int)ceil(ucl*d_res);
 k_max = (int)qf_pois(.99999999, mu);
 N = u - l;
 NN = N + 1;
 a = matrix(NN, NN);
 g = vector(NN);

 for (i=0; i<=N; i++) for (j=0; j<=N; j++) a[i*NN+j] = 0.;
 
 for (i=0; i<=N; i++) {   
   for (k=0; k<=k_max; k++) {
     zold = ( (double)( l + i ) ) / dr;
     zj = (1.-lambda)*zold + lambda*k;  
     pj = pdf_pois((double)k, mu);
     switch (round_mode) {
        case -1:	/* round down as probably Gan did */
	    j = (int)floor(zj*dr + 1e-9) - l;
	    if ( 0 <= j && j < N ) a[j*NN+i] += -pj;
        if ( N <= j ) a[N*NN+i] += -pj;
        break;
        case 0:		/* round down */
	    j = (int)floor(zj*dr) - l;
	    if ( 0 <= j && j < N ) a[j*NN+i] += -pj;
        if ( N <= j ) a[N*NN+i] += -pj;
        break;
        case 1:		/* round up */
	    j = (int)ceil(zj*dr) - l;
	    if ( 0 <= j && j < N ) a[j*NN+i] += -pj;
        if ( N <= j ) a[N*NN+i] += -pj;
        break;
        case 2:		/* round to nearest -- round half to even, IEEE 754 */
	    j = (int)round(zj*dr) - l;
	    if ( 0 <= j && j < N ) a[j*NN+i] += -pj;
        if ( N <= j ) a[N*NN+i] += -pj;
        break;
        case 3:		/* round to nearest -- round half up  */
	    j = (int)floor(zj*dr+.5) - l;
	    if ( 0 <= j && j < N ) a[j*NN+i] += -pj;
        if ( N <= j ) a[N*NN+i] += -pj;
        break;
        case 4:		/* distribute */
	    j = (int)floor(zj*dr) - l;
	    pju = zj - j/dr;
        if ( 0 <= j && j < N ) a[j*NN+i]     += -(1.-pju)*pj;
	    if ( 0 <  j && j < N ) a[(j+1)*NN+i] += -pju*pj;
        if ( N <= j ) a[N*NN+i] += -pj;
        break;         
     }
   }
   ++a[i*NN+i];
 } 

 for (j=0; j<=N; j++) g[j] = 1.;
 solve(&NN, a, g);
 
 arl = 1.;
 for (k=0; k<=k_max; k++) {
   zj = (1.-lambda)*z0 + lambda*k;
   pj = pdf_pois((double)k, mu);
   switch (round_mode) {
      case -1:	/* round down as probably Gan did */
	  j = (int)floor(zj*dr + 1e-9) - l;
	  if ( 0 <= j && j < N ) arl += pj*g[j];
      if ( N <= j ) arl += pj*g[N];
      break;
      case 0:	/* round down */
	  j = (int)floor(zj*dr) - l;
	  if ( 0 <= j && j < N ) arl += pj*g[j];
      if ( N <= j ) arl += pj*g[N];
      break;
      case 1:	/* round up */
	  j = (int)ceil(zj*dr) - l;
	  if ( 0 <= j && j < N ) arl += pj*g[j];
      if ( N <= j ) arl += pj*g[N];
      break;
      case 2:	/* round to nearest -- round half to even, IEEE 754 */
	  j = (int)round(zj*dr) - l;
	  if ( 0 <= j && j < N ) arl += pj*g[j];
      if ( N <= j ) arl += pj*g[N];
      break;
      case 3:	/* round to nearest -- round half up  */
	  j = (int)floor(zj*dr+.5) - l;
	  if ( 0 <= j && j < N ) arl += pj*g[j];
      if ( N <= j ) arl += pj*g[N];
      break;
      case 4:	/* distribute */
	  j = (int)floor(zj*dr) - l;
	  pju = zj - j/dr;          
	  if ( 0 <= j && j < N ) arl += (1.-pju)*pj*g[j];
	  if ( 0 <  j && j < N ) arl += pju*pj*g[j+1];
      if ( N <= j ) arl += pj*g[N];
      break;
   }
 }     

 Free(a);
 Free(g);

 return arl;
}


double ewma_p2_arl(double lambda, double lcl, double ucl, int n, double p, double z0, int d_res, int round_mode, int mid_mode)
{ double *a, *g, arl, zj=0, pju, pj;
  int i, j, k, N, N1, N2, NN/*, k_max*/;
  
 N2 = (int)ceil(ucl*d_res);
 N1 = (int)floor(lcl*d_res);
 N = N2 - N1;
 NN = N + 1;
 a = matrix(NN, NN);
 g = vector(NN);

 for (i=0; i<=N; i++) for (j=0; j<=N; j++) a[i*NN+j] = 0.;
 
 for (i=0; i<=N; i++) {
   /*k_max = (int)ceil( (ucl+1. - (1.-lambda)*i)/lambda );*/
   for (k=0; k<=n; k++) {
     zj = (1.-lambda)*(N1+i)/d_res + lambda*k;  
     pj = pdf_binom((double)k, n, p);
     switch (round_mode) {
	case -1:	/* round down as probably Gan did */
	    j = (int)floor(zj*d_res + 1e-9) - N1;
	    if ( 0 <= j && j <= N ) a[i*NN+j] += -pj;
	break;
	case 0:		/* round down */
	    j = (int)floor(zj*d_res) - N1;
	    if ( 0 <= j && j <= N ) a[i*NN+j] += -pj;
	break;
	case 1:		/* round up */
	    j = (int)ceil(zj*d_res) - N1;
	    if ( 0 <= j && j <= N ) a[i*NN+j] += -pj;
	break;
	case 2:		/* round to nearest -- round half to even, IEEE 754 */
	    j = (int)round(zj*d_res) - N1;
	    if ( 0 <= j && j <= N ) a[i*NN+j] += -pj;
	break;
	case 3:		/* round to nearest -- round half up  */
	    j = (int)floor(zj*d_res+.5) - N1;
	    if ( 0 <= j && j <= N ) a[i*NN+j] += -pj;
	break;
	case 4:		/* distribute */
	    j = (int)floor(zj*d_res) - N1;
	    pju = zj - j/d_res;          
	    if ( 0 <= j && j <= N ) a[i*NN+j]   += -(1.-pju)*pj;
	    if ( 0 <= j && j <  N ) a[i*NN+j+1] += -pju*pj;
	break;
     }
   }
   ++a[i*NN+i];
 } 

 for (j=0; j<=N; j++) g[j] = 1.;
 LU_solve(a, g, NN);
 /*solve(&NN, a, g);*/
 
 arl = 1.;
 /*k_max = (int)ceil( (ucl+1. - (1.-lambda)*z0)/lambda );*/
 for (k=0; k<=n; k++) {
   zj = (1.-lambda)*z0 + lambda*k;
   pj = pdf_binom((double)k, n, p);
   switch (round_mode) {
      case -1:	/* round down as probably Gan did */
	  j = (int)floor(zj*d_res + 1e-9) - N1;
	  if ( 0 <= j && j <= N ) arl += pj*g[j];
      break;
      case 0:	/* round down */
	  j = (int)floor(zj*d_res) - N1;
	  if ( 0 <= j && j <= N ) arl += pj*g[j];
      break;
      case 1:	/* round up */
	  j = (int)ceil(zj*d_res) - N1;
	  if ( 0 <= j && j <= N ) arl += pj*g[j];
      break;
      case 2:	/* round to nearest -- round half to even, IEEE 754 */
	  j = (int)round(zj*d_res) - N1;
	  if ( 0 <= j && j <= N ) arl += pj*g[j];
      break;
      case 3:	/* round to nearest -- round half up  */
	  j = (int)floor(zj*d_res+.5) - N1;
	  if ( 0 <= j && j <= N ) arl += pj*g[j];
      break;
      case 4:	/* distribute */
	  j = (int)floor(zj*d_res) - N1;
	  pju = zj - j/d_res;          
	  if ( 0 <= j && j <= N ) arl += (1.-pju)*pj*g[j];
	  if ( 0 <= j && j <  N ) arl += pju*pj*g[j+1];
      break;
   }
 }     

 Free(a);
 Free(g);

 return arl;
}


double ewma_c2_arl(double lambda, double lcl, double ucl, double mu, double z0, int d_res, int round_mode, int mid_mode)
{ double *a, *g, arl, zj=0, lzi=0, pju, pj;
  int i, j, k, N, N1, N2, NN, k_max;
 
 N1 = (int)ceil(lcl*d_res);  
 N2 = (int)floor(ucl*d_res);
 if ( round_mode == 4 ) {
   N1 = (int)ceil(lcl*d_res);  
   N2 = (int)floor(ucl*d_res);
 }
 N = N2 - N1;
 NN = N + 1;
 a = matrix(NN, NN);
 g = vector(NN);

 for (i=0; i<=N; i++) for (j=0; j<=N; j++) a[i*NN+j] = 0.;
 
 for (i=0; i<=N; i++) {
   lzi = (1.-lambda)*(N1+i)/d_res;
   k_max = (int)ceil( (ucl+1. - lzi)/lambda );
   for (k=0; k<=k_max; k++) {
     zj = lzi + lambda*k;  
     pj = pdf_pois((double)k, mu);
     switch (round_mode) {
	case -1:	/* round down as probably Gan did */
	    j = (int)floor(zj*d_res + 1e-9) - N1;
	    if ( 0 <= j && j <= N ) a[i*NN+j] += -pj;
	break;
	case 0:		/* round down */
	    j = (int)floor(zj*d_res) - N1;
	    if ( 0 <= j && j <= N ) a[i*NN+j] += -pj;
	break;
	case 1:		/* round up */
	    j = (int)ceil(zj*d_res) - N1;
	    if ( 0 <= j && j <= N ) a[i*NN+j] += -pj;
	break;
	case 2:		/* round to nearest -- round half to even, IEEE 754 */
	    j = (int)round(zj*d_res) - N1;
	    if ( 0 <= j && j <= N ) a[i*NN+j] += -pj;
	break;
	case 3:		/* round to nearest -- round half up  */
	    j = (int)floor(zj*d_res+.5) - N1;
	    if ( 0 <= j && j <= N ) a[i*NN+j] += -pj;
	break;
	case 4:		/* distribute */
	    j = (int)floor(zj*d_res) - N1;
	    pju = zj - (j+N1)/d_res;          
	    if (  0 <= j && j <= N ) a[i*NN + j]   += -(1.-pju)*pj;
	    if ( -1 <= j && j <  N ) a[i*NN + j+1] += -pju*pj;
	break;
     }
   }
   ++a[i*NN+i];
 } 

 for (j=0; j<=N; j++) g[j] = 1.;
 LU_solve(a, g, NN);
 /*solve(&NN, a, g);*/
 
 arl = 1.;
 k_max = (int)ceil( (ucl+1. - (1.-lambda)*z0)/lambda );   
 for (k=0; k<=k_max; k++) {
   zj = (1.-lambda)*z0 + lambda*k;
   pj = pdf_pois((double)k, mu);
   switch (round_mode) {
      case -1:	/* round down as probably Gan did */
	  j = (int)floor(zj*d_res + 1e-9) - N1;
	  if ( 0 <= j && j <= N ) arl += pj*g[j];
      break;
      case 0:	/* round down */
	  j = (int)floor(zj*d_res) - N1;
	  if ( 0 <= j && j <= N ) arl += pj*g[j];
      break;
      case 1:	/* round up */
	  j = (int)ceil(zj*d_res) - N1;
	  if ( 0 <= j && j <= N ) arl += pj*g[j];
      break;
      case 2:	/* round to nearest -- round half to even, IEEE 754 */
	  j = (int)round(zj*d_res) - N1;
	  if ( 0 <= j && j <= N ) arl += pj*g[j];
      break;
      case 3:	/* round to nearest -- round half up  */
	  j = (int)floor(zj*d_res+.5) - N1;
	  if ( 0 <= j && j <= N ) arl += pj*g[j];
      break;
      case 4:	/* distribute */
	  j = (int)floor(zj*d_res) - N1;
	  pju = zj - (j+N1)/d_res;          
	  if ( 0 <= j && j <= N ) arl += (1.-pju)*pj*g[j];
	  if ( -1 <= j && j <  N ) arl += pju*pj*g[j+1];
      break;
   }
 }     

 Free(a);
 Free(g);

 return arl;
}

/* Markov chain model from Borror / Champ / Rigdon (1998) "Poisson EWMA control charts", JQT 30(4), 352-361 */

double cewma_2_arl(double lambda, double AL, double AU, double mu0, double z0, double mu, int N)
{ double hL, hU, w, *a, *g, arl;
  int i, j;

 a = matrix(N,N);
 g = vector(N);

 hL = mu0 - AL * sqrt( lambda*mu0 / (2.-lambda) );
 hU = mu0 + AU * sqrt( lambda*mu0 / (2.-lambda) );
 w  = (hU - hL)/N;

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[j*N+i] = - ( cdf_pois( hL+w/2./lambda*(2.*(j+1.)-(1.-lambda)*(2.*i+1.)), mu) - cdf_pois( hL+w/2./lambda*(2.*j-(1.-lambda)*(2.*i+1.)), mu) );
   ++a[i*N+i];
 }

 for (j=0; j<N; j++) g[j] = 1.;
 solve(&N, a, g);
 
 arl = 1.;
 for (j=0; j<N; j++) arl += ( cdf_pois( ( hL + (j+1.)*w - (1.-lambda)*z0 )/lambda, mu) - cdf_pois( ( hL + j*w - (1.-lambda)*z0 )/lambda, mu) ) * g[j];

 Free(a);
 Free(g);

 return arl;  
}


double cewma_2_arl_rando(double lambda, double AL, double AU, double gammaL, double gammaU, double mu0, double z0, double mu, int N)
{ double hL, hU, w, *a, *g, arl;
  int i, j;

 a = matrix(N,N);
 g = vector(N);

 hL = mu0 - AL * sqrt( lambda*mu0 / (2.-lambda) );
 hU = mu0 + AU * sqrt( lambda*mu0 / (2.-lambda) );
 w  = (hU - hL)/N;

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[j*N+i] = - ( cdf_pois( hL+w/2./lambda*(2.*(j+1.)-(1.-lambda)*(2.*i+1.)), mu) - cdf_pois( hL+w/2./lambda*(2.*j-(1.-lambda)*(2.*i+1.)), mu) );
   a[i] *= (1.-gammaL);
   a[(N-1)*N+i] *= (1.-gammaU);  
   ++a[i*N+i];
 }

 for (j=0; j<N; j++) g[j] = 1.;
 solve(&N, a, g);
 
 arl = 1. + (1.-gammaL) * ( cdf_pois( ( hL + w - (1.-lambda)*z0 )/lambda, mu) - cdf_pois( ( hL - (1.-lambda)*z0 )/lambda, mu) ) * g[0];
 for (j=1; j<N-1; j++) arl += ( cdf_pois( ( hL + (j+1.)*w - (1.-lambda)*z0 )/lambda, mu) - cdf_pois( ( hL + j*w - (1.-lambda)*z0 )/lambda, mu) ) * g[j];
 arl += (1.-gammaU) * ( cdf_pois( ( hL + N*w - (1.-lambda)*z0 )/lambda, mu) - cdf_pois( ( hL +(N-1.)*w - (1.-lambda)*z0 )/lambda, mu) ) * g[N-1];

 Free(a);
 Free(g);

 return arl;  
}


double cewma_2_arl_new(double lambda, double AL, double AU, double mu0, double z0, double mu, int N)
{ double hL, hU, w, *a, *g, arl, zi, zj, zi1, zi2, pj, px;
  int i, j, x, x0, x1;

 a = matrix(N,N);
 g = vector(N);

 hL = mu0 - AL * sqrt( lambda*mu0 / (2.-lambda) );
 hU = mu0 + AU * sqrt( lambda*mu0 / (2.-lambda) );
 w  = (hU - hL)/N;

 /*printf("\n\nLos geht es (hL = %.4f, hU = %.4f, w = %.4f)\n\n", hL, hU, w);*/
 
 for (i=0; i<N; i++) {
   zi = (1.-lambda) * ( hL + (2.*i+1.)/2.*w );
   x0 = (int)floor( (hL-zi)/lambda );
   if ( x0 < 0 ) x0 = 0;
   x1 =  (int)ceil( (hU-zi)/lambda );
   for (j=0; j<N; j++) a[j*N+i] = 0.;
   for (x=x0; x<=x1; x++) {
     px = pdf_pois((double)x, mu);  
     zj = zi + (double)x*lambda;
     j = (int)ceil( (zj-hL)/w ) - 1;     
     zi1 = ( hL + j*w - (double)x*lambda ) / (1.-lambda);
     zi2 = ( hL + (j+1.)*w - (double)x*lambda ) / (1.-lambda);     
     if ( zi1 <= hL + i*w ) {
       if ( hL + (i+1.)*w <= zi2 ) {
         pj = 1.;
         if ( j >= 0 && j < N ) a[j*N+i] += - pj * px;
         /*printf("(ii)\ti = %d,\tx = %d,\tj = %d\n", i, x, j);*/
       } else {
         pj = ( zi2 - ( hL + i*w ) ) / w;
         if ( j >= 0  && j < N  ) a[j*N+i] += - pj * px;
         if ( j >= -1 && j < N-1 ) a[(j+1)*N+i] += - (1.-pj) * px;
         /*printf("(iv)\ti = %d,\tx = %d,\tj = %d,\tpj = %.4f\n", i, x, j, pj);*/
       }
     } else {
       if ( hL + (i+1.)*w <= zi2 ) {
         pj = ( hL + (i+1.)*w - zi1 ) / w;
         if ( j >= 0 && j < N  ) a[j*N+i] += - pj * px;
         if ( j > 0  && j <= N ) a[(j-1)*N+i] += - (1.-pj) * px;           
         /*printf("(iii)\ti = %d,\tx = %d,\tj = %d,\tpj = %.4f\n", i, x, j, pj);*/
       } else {
         pj = 0.; /* should not be possible */
         /*printf("(i)\tshould not happen.\n");*/
       }
     }
     /*printf("i = %d,\tx = %d,\tj = %d,\tmj = %.4f,\tpj = %.4f,\tzj = %.4f,\tzj-mj = %.4f\n", i, x, j, mj, pj, zj, zj-mj);*/
   }
   /*printf("\n\n");*/
   ++a[i*N+i];
 }

 for (j=0; j<N; j++) g[j] = 1.;
 solve(&N, a, g);

 arl = 1.; 
 zi = (1.-lambda) * z0;     
 x0 = (int)floor( (hL-zi)/lambda );
 if ( x0 < 0 ) x0 = 0;
 x1 =  (int)ceil( (hU-zi)/lambda ); 
 i = (int)ceil( (z0-hL)/w ) - 1;
 for (x=x0; x<=x1; x++) {
   px = pdf_pois((double)x, mu);  
   zj = zi+(double)x*lambda;
   j = (int)ceil( (zj-hL)/w ) - 1;
   zi1 = ( hL + j*w - (double)x*lambda ) / (1.-lambda);
   zi2 = ( hL + (j+1.)*w - (double)x*lambda ) / (1.-lambda);     
   if ( zi1 <= hL + i*w ) {
     if ( hL + (i+1.)*w <= zi2 ) {
       pj = 1.;
       if ( j >= 0 && j < N ) arl += pj * px * g[j];
     } else {
       pj = ( zi2 - ( hL + i*w ) ) / w;
       if ( j >= 0  && j < N   ) arl += pj * px * g[j];
       if ( j >= -1 && j < N-1 ) arl += (1.-pj) * px * g[j+1];
     }
   } else {
     if ( hL + (i+1.)*w <= zi2 ) {
       pj = ( hL + (i+1.)*w - zi1 ) / w;
       if ( j >= 0 && j < N  ) arl += pj * px * g[j];
       if ( j > 0  && j <= N ) arl += (1.-pj) * px * g[j-1];
     } else {
       pj = 0.; /* should not be possible */  
     }
   }
 }

 Free(a);
 Free(g);

 return arl;  
}


double cewma_2_arl_rando_new(double lambda, double AL, double AU, double gammaL, double gammaU, double mu0, double z0, double mu, int N)
{ double hL, hU, w, *a, *g, arl, zi, zj, zi1, zi2, pj, qj, px;
  int i, j, x, x0, x1;

 a = matrix(N,N);
 g = vector(N);

 hL = mu0 - AL * sqrt( lambda*mu0 / (2.-lambda) );
 hU = mu0 + AU * sqrt( lambda*mu0 / (2.-lambda) );
 w  = (hU - hL)/N;
 
 for (i=0; i<N; i++) {
   zi = (1.-lambda) * ( hL + (2.*i+1.)/2.*w );
   x0 = (int)floor( (hL-zi)/lambda );
   if ( x0 < 0 ) x0 = 0;
   x1 =  (int)ceil( (hU-zi)/lambda );
   for (j=0; j<N; j++) a[j*N+i] = 0.;
   for (x=x0; x<=x1; x++) {
     px = pdf_pois((double)x, mu);  
     zj = zi + (double)x*lambda;
     j = (int)ceil( (zj-hL)/w ) - 1;     
     zi1 = ( hL + j*w - (double)x*lambda ) / (1.-lambda);
     zi2 = ( hL + (j+1.)*w - (double)x*lambda ) / (1.-lambda);     
     if ( zi1 <= hL + i*w ) {
       if ( hL + (i+1.)*w <= zi2 ) {
         pj = 1.;
         if ( j==0 )   pj *= (1.-gammaL);
         if ( j==N-1 ) pj *= (1.-gammaU);
         if ( j >= 0 && j < N ) a[j*N+i] += - pj * px;
       } else {
         pj = ( zi2 - (hL+i*w) ) / w;
         qj = 1. - pj;
         if ( j >= 0  && j < N ) {
           if ( j==0 )   pj *= (1.-gammaL);
           if ( j==N-1 ) pj *= (1.-gammaU);
           a[j*N+i] += - pj * px;
         }
         if ( j >= -1 && j < N-1 ) {           
           if ( j==-1 )  qj *= (1.-gammaL);
           if ( j==N-2 ) qj *= (1.-gammaU);
           a[(j+1)*N+i] += - qj * px;
         }
       }
     } else {
       if ( hL + (i+1.)*w <= zi2 ) {
         pj = ( hL + (i+1.)*w - zi1 ) / w;
         qj = 1. - pj;
         if ( j >= 0 && j < N ) {
           if ( j==0 )   pj *= (1.-gammaL);
           if ( j==N-1 ) pj *= (1.-gammaU);
           a[j*N+i] += - pj * px;
         }
         if ( j > 0 && j <= N ) {           
           if ( j==1 ) qj *= (1.-gammaL);
           if ( j==N ) qj *= (1.-gammaU);
           a[(j-1)*N+i] += - qj * px;
         }
       } else {
         pj = 0.; /* should not be possible */
       }
     }
   }  
   ++a[i*N+i];
 }

 for (j=0; j<N; j++) g[j] = 1.;
 solve(&N, a, g);

 arl = 1.; 
 zi = (1.-lambda) * z0;     
 x0 = (int)floor( (hL-zi)/lambda );
 if ( x0 < 0 ) x0 = 0;
 x1 =  (int)ceil( (hU-zi)/lambda ); 
 i = (int)ceil( (z0-hL)/w ) - 1;
 for (x=x0; x<=x1; x++) {
   px = pdf_pois((double)x, mu);  
   zj = zi+(double)x*lambda;
   j = (int)ceil( (zj-hL)/w ) - 1;
   zi1 = ( hL + j*w - (double)x*lambda ) / (1.-lambda);
   zi2 = ( hL + (j+1.)*w - (double)x*lambda ) / (1.-lambda);     
   if ( zi1 <= hL + i*w ) {
     if ( hL + (i+1.)*w <= zi2 ) {
       pj = 1.;
       if ( j==0 )   pj *= (1.-gammaL);
       if ( j==N-1 ) pj *= (1.-gammaU);
       if ( j >= 0 && j < N ) arl += pj * px * g[j];
     } else {
       pj = ( zi2 - (hL+i*w) ) / w;
       qj = 1. - pj;  
       if ( j >= 0  && j < N ) {
         if ( j==0 )   pj *= (1.-gammaL);
         if ( j==N-1 ) pj *= (1.-gammaU);
         arl += pj * px * g[j];
       } 
       if ( j >= -1 && j < N-1 ) {
         if ( j==-1 )  qj *= (1.-gammaL);
         if ( j==N-2 ) qj *= (1.-gammaU);
         arl += qj * px * g[j+1];
       }
     }
   } else {
     if ( hL + (i+1.)*w <= zi2 ) {
       pj = ( hL+(i+1.)*w - zi1 ) / w;
       qj = 1. - pj;
       if ( j >= 0 && j < N ) {
         if ( j==0 )   pj *= (1.-gammaL);
         if ( j==N-1 ) pj *= (1.-gammaU);
         arl += pj * px * g[j];
       }
       if ( j > 0  && j <= N ) {
         if ( j==1 ) qj *= (1.-gammaL);
         if ( j==N ) qj *= (1.-gammaU);
         arl += qj * px * g[j-1];
       }
     } else {
       pj = 0.; /* should not be possible */  
     }
   }
 }

 Free(a);
 Free(g);

 return arl;  
}


double cewma_U_arl(double lambda, double AU, double mu0, double z0, double mu, int N)
{ double hL, hU, w, *a, *g, arl;
  int i, j;

 a = matrix(N,N);
 g = vector(N);

 hL = 0.;
 hU = mu0 + AU * sqrt( lambda*mu0 / (2.-lambda) );
 w  = hU/N;

 for (i=0; i<N; i++) {
   for (j=0; j<N; j++) a[j*N+i] = - ( cdf_pois( hL+w/2./lambda*(2.*(j+1.)-(1.-lambda)*(2.*i+1.)) , mu) - cdf_pois( hL+w/2./lambda*(2.*j-(1.-lambda)*(2.*i+1.)) , mu) );
   ++a[i*N+i];
 }

 for (j=0; j<N; j++) g[j] = 1.;
 solve(&N, a, g);

 arl = 1.;
 for (j=0; j<N; j++) arl += ( cdf_pois( ( hL + (j+1.)*w - (1.-lambda)*z0 )/lambda, mu) - cdf_pois( ( hL + j*w - (1.-lambda)*z0 )/lambda, mu) ) * g[j];

 Free(a);
 Free(g);

 return arl;  
}


double cewma_U_arl_new(double lambda, double AU, double mu0, double z0, double mu, int N)   
{ double hL, hU, w, *a, *g, arl, zi, zj, zi1, zi2, pj, px;
  int i, j, x, x1;
  
 a = matrix(N,N);
 g = vector(N);

 hL = 0.;
 hU = mu0 + AU * sqrt( lambda*mu0 / (2.-lambda) );
 w  = hU/N;

 for (i=0; i<N; i++) {
   zi = (1.-lambda) * ( hL + (2.*i+1.)/2.*w );
   x1 =  (int)ceil( (hU-zi)/lambda );
   for (j=0; j<N; j++) a[j*N+i] = 0.;
   for (x=0; x<=x1; x++) {
     px = pdf_pois((double)x, mu);  
     zj = zi + (double)x*lambda;
     j = (int)ceil( (zj-hL)/w ) - 1;     
     zi1 = ( hL + j*w - (double)x*lambda ) / (1.-lambda);
     zi2 = ( hL + (j+1.)*w - (double)x*lambda ) / (1.-lambda);     
     if ( zi1 <= hL + i*w ) {
       if ( hL + (i+1.)*w <= zi2 ) {
         pj = 1.;
         if ( j >= 0 && j < N ) a[j*N+i] += - pj * px;
         /*printf("(ii)\ti = %d,\tx = %d,\tj = %d\n", i, x, j);*/
       } else {
         pj = ( zi2 - ( hL + i*w ) ) / w;
         if ( j >= 0  && j < N  ) a[j*N+i] += - pj * px;
         if ( j >= -1 && j < N-1 ) a[(j+1)*N+i] += - (1.-pj) * px;
         /*printf("(iv)\ti = %d,\tx = %d,\tj = %d,\tpj = %.4f\n", i, x, j, pj);*/
       }
     } else {
       if ( hL + (i+1.)*w <= zi2 ) {
         pj = ( hL + (i+1.)*w - zi1 ) / w;
         if ( j >= 0 && j < N  ) a[j*N+i] += - pj * px;
         if ( j > 0  && j <= N ) a[(j-1)*N+i] += - (1.-pj) * px;           
         /*printf("(iii)\ti = %d,\tx = %d,\tj = %d,\tpj = %.4f\n", i, x, j, pj);*/
       } else {
         pj = 0.; /* should not be possible */
         /*printf("(i)\tshould not happen.\n");*/
       }
     }
     /*printf("i = %d,\tx = %d,\tj = %d,\tmj = %.4f,\tpj = %.4f,\tzj = %.4f,\tzj-mj = %.4f\n", i, x, j, mj, pj, zj, zj-mj);*/
   }
   /*printf("\n\n");*/
   ++a[i*N+i];
 }

 for (j=0; j<N; j++) g[j] = 1.;
 solve(&N, a, g);

 arl = 1.; 
 zi = (1.-lambda) * z0;     
 x1 =  (int)ceil( (hU-zi)/lambda ); 
 i = (int)ceil( (z0-hL)/w ) - 1;
 for (x=0; x<=x1; x++) {
   px = pdf_pois((double)x, mu);  
   zj = zi+(double)x*lambda;
   j = (int)ceil( (zj-hL)/w ) - 1;
   zi1 = ( hL + j*w - (double)x*lambda ) / (1.-lambda);
   zi2 = ( hL + (j+1.)*w - (double)x*lambda ) / (1.-lambda);     
   if ( zi1 <= hL + i*w ) {
     if ( hL + (i+1.)*w <= zi2 ) {
       pj = 1.;
       if ( j >= 0 && j < N ) arl += pj * px * g[j];
     } else {
       pj = ( zi2 - ( hL + i*w ) ) / w;
       if ( j >= 0  && j < N   ) arl += pj * px * g[j];
       if ( j >= -1 && j < N-1 ) arl += (1.-pj) * px * g[j+1];
     }
   } else {
     if ( hL + (i+1.)*w <= zi2 ) {
       pj = ( hL + (i+1.)*w - zi1 ) / w;
       if ( j >= 0 && j < N  ) arl += pj * px * g[j];
       if ( j > 0  && j <= N ) arl += (1.-pj) * px * g[j-1];
     } else {
       pj = 0.; /* should not be possible */  
     }
   }
 }

 Free(a);
 Free(g);

 return arl;  
}


double cewma_L_arl(double lambda, double AL, double AU, double mu0, double z0, double mu, int N)
{ double hL, hU, w, *a, *g, arl;
  int i, j;

 a = matrix(N,N);
 g = vector(N);

 hL = mu0 - AL * sqrt( lambda*mu0 / (2.-lambda) );
 hU = mu0 + AU * sqrt( lambda*mu0 / (2.-lambda) );
 w  = (hU - hL)/N;

 for (i=0; i<N; i++) {
   for (j=0; j<(N-1); j++) a[j*N+i] = - ( cdf_pois( hL+w/2./lambda*(2.*(j+1.)-(1.-lambda)*(2.*i+1.)) , mu) - cdf_pois( hL+w/2./lambda*(2.*j-(1.-lambda)*(2.*i+1.)) , mu) );
   a[(N-1)*N+i] = - ( 1. - cdf_pois( hL+w/2./lambda*(2.*j-(1.-lambda)*(2.*i+1.)) , mu) );
   ++a[i*N+i];
 }

 for (j=0; j<N; j++) g[j] = 1.;
 solve(&N, a, g);

 arl = 1.;
 for (j=0; j<(N-1); j++) arl += ( cdf_pois( ( hL + (j+1.)*w - (1.-lambda)*z0 )/lambda, mu) - cdf_pois( ( hL + j*w - (1.-lambda)*z0 )/lambda, mu) ) * g[j];
 arl += ( 1. - cdf_pois( ( hL + j*w - (1.-lambda)*z0 )/lambda, mu) ) * g[N-1];
 
 Free(a);
 Free(g);

 return arl;  
}


double cewma_L_arl_new(double lambda, double AL, double AU, double mu0, double z0, double mu, int N)
{ double hL, hU, w, *a, *g, arl, zi, zj, zi1, zi2, pj, px;
  int i, j, x, x0, x1;
  
 a = matrix(N,N);
 g = vector(N);

 hL = mu0 - AL * sqrt( lambda*mu0 / (2.-lambda) );
 hU = mu0 + AU * sqrt( lambda*mu0 / (2.-lambda) );
 w  = (hU - hL)/N;

 for (i=0; i<N; i++) {
   zi = (1.-lambda) * ( hL + (2.*i+1.)/2.*w );
   x0 = (int)floor( (hL-zi)/lambda );
   if ( x0 < 0 ) x0 = 0;
   x1 =  (int)ceil( (hU-zi)/lambda );
   for (j=0; j<N; j++) a[j*N+i] = 0.;
   for (x=0; x<x1; x++) {
     px = pdf_pois((double)x, mu);  
     zj = zi + (double)x*lambda;
     j = (int)ceil( (zj-hL)/w ) - 1;     
     zi1 = ( hL + j*w - (double)x*lambda ) / (1.-lambda);
     zi2 = ( hL + (j+1.)*w - (double)x*lambda ) / (1.-lambda);     
     if ( zi1 <= hL + i*w ) {
       if ( hL + (i+1.)*w <= zi2 ) {
         pj = 1.;
         if ( j >= 0 && j < N ) a[j*N+i] += - pj * px;
         /*printf("(ii)\ti = %d,\tx = %d,\tj = %d\n", i, x, j);*/
       } else {
         pj = ( zi2 - ( hL + i*w ) ) / w;
         if ( j >= 0  && j < N  ) a[j*N+i] += - pj * px;
         if ( j >= -1 && j < N-1 ) a[(j+1)*N+i] += - (1.-pj) * px;
         if ( j >= -1 && j == N-1 ) a[(N-1)*N+i] += - (1.-pj) * px;
         /*printf("(iv)\ti = %d,\tx = %d,\tj = %d,\tpj = %.4f\n", i, x, j, pj);*/
       }
     } else {
       if ( hL + (i+1.)*w <= zi2 ) {
         pj = ( hL + (i+1.)*w - zi1 ) / w;
         if ( j >= 0 && j < N  ) a[j*N+i] += - pj * px;
         if ( j > 0  && j <= N ) a[(j-1)*N+i] += - (1.-pj) * px;           
         /*printf("(iii)\ti = %d,\tx = %d,\tj = %d,\tpj = %.4f\n", i, x, j, pj);*/
       } else {
         pj = 0.; /* should not be possible */
         /*printf("(i)\tshould not happen.\n");*/
       }
     }
   }
   a[(N-1)*N+i] += - ( 1. - cdf_pois( (double)x1-1., mu) );
   ++a[i*N+i];
 }
 
 
 for (j=0; j<N; j++) g[j] = 1.;
 solve(&N, a, g);

 arl = 1.; 
 zi = (1.-lambda) * z0;     
 x0 = (int)floor( (hL-zi)/lambda );
 if ( x0 < 0 ) x0 = 0;
 x1 =  (int)ceil( (hU-zi)/lambda ); 
 i = (int)ceil( (z0-hL)/w ) - 1;
 for (x=x0; x<x1; x++) {
   px = pdf_pois((double)x, mu);  
   zj = zi+(double)x*lambda;
   j = (int)ceil( (zj-hL)/w ) - 1;
   zi1 = ( hL + j*w - (double)x*lambda ) / (1.-lambda);
   zi2 = ( hL + (j+1.)*w - (double)x*lambda ) / (1.-lambda);     
   if ( zi1 <= hL + i*w ) {
     if ( hL + (i+1.)*w <= zi2 ) {
       pj = 1.;
       if ( j >= 0 && j < N ) arl += pj * px * g[j];
     } else {
       pj = ( zi2 - ( hL + i*w ) ) / w;
       if ( j >= 0  && j < N   ) arl += pj * px * g[j];
       if ( j >= -1 && j < N-1 ) arl += (1.-pj) * px * g[j+1];
       if ( j >= -1 && j == N-1 ) arl += (1.-pj) * px * g[N-1];
     }
   } else {
     if ( hL + (i+1.)*w <= zi2 ) {
       pj = ( hL + (i+1.)*w - zi1 ) / w;
       if ( j >= 0 && j < N  ) arl += pj * px * g[j];
       if ( j > 0  && j <= N ) arl += (1.-pj) * px * g[j-1];
     } else {
       pj = 0.; /* should not be possible */  
     }
   }
 }
 arl += ( 1. - cdf_pois( (double)x1-1., mu) ) * g[N-1];
 
 Free(a);
 Free(g);

 return arl;  
}


double cewma_U_crit(double lambda, double L0, double mu0, double z0, int N, int jmax)
{ double A, A1, L1;
  int j, dA;
   
 A1 = floor(mu0);
 if ( A1 < 1. ) A1 = 1.;
 L1 = 1.;
 A = 1.;
 for (j=1; j<=(int)A1; j++) {
   A = (double)j;
   L1 = cewma_U_arl(lambda, A, mu0, z0, mu0, N);
   if ( L1 > L0 ) j = (int)A1 + 1; 
 }
 A1 = A;
 
 if ( L1 > L0 ) {
   for (j=1; j<=jmax; j++) {
     for (dA=1; dA<20; dA++) {
       A = A1 + dA / pow(-10., (double)j);
       L1 = cewma_U_arl(lambda, A, mu0, z0, mu0, N);
       if ( ( fmod((double)j, 2.) > 0. && L1 < L0 ) || ( fmod((double)j, 2.) < 1. && L1 > L0 ) ) dA = 20;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 } else {
   for (j=1; j<=jmax; j++) {
     for (dA=1; dA<20; dA++) {
       A = A1 - dA / pow(-10., (double)j);
       L1 = cewma_U_arl(lambda, A, mu0, z0, mu0, N);
       if ( ( fmod((double)j, 2.) < 1. && L1 < L0 ) || ( fmod((double)j, 2.) > 0. && L1 > L0 ) ) dA = 20;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 }
 
 return A;
}


double cewma_U_crit_new(double lambda, double L0, double mu0, double z0, int N, int jmax)
{ double A, A1, L1;
  int j, dA;
    
 A1 = floor(mu0);
 if ( A1 < 1. ) A1 = 1.;
 L1 = 1.;
 A = 1.;
 for (j=1; j<=(int)A1; j++) {
   A = (double)j;
   L1 = cewma_U_arl_new(lambda, A, mu0, z0, mu0, N);
   if ( L1 > L0 ) j = (int)A1 + 1; 
 }
 A1 = A;
 
 if ( L1 > L0 ) {
   for (j=1; j<=jmax; j++) {
     for (dA=1; dA<20; dA++) {
       A = A1 + dA / pow(-10., (double)j);
       L1 = cewma_U_arl_new(lambda, A, mu0, z0, mu0, N);
       if ( ( fmod((double)j, 2.) > 0. && L1 < L0 ) || ( fmod((double)j, 2.) < 1. && L1 > L0 ) ) dA = 20;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 } else {
   for (j=1; j<=jmax; j++) {
     for (dA=1; dA<20; dA++) {
       A = A1 - dA / pow(-10., (double)j);
       L1 = cewma_U_arl_new(lambda, A, mu0, z0, mu0, N);
       if ( ( fmod((double)j, 2.) < 1. && L1 < L0 ) || ( fmod((double)j, 2.) > 0. && L1 > L0 ) ) dA = 20;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 }
 
 return A;
}


double cewma_L_crit(double lambda, double L0, double AU, double mu0, double z0, int N, int jmax)
{ double A, A1, L1;
  int j, dA;
   
 A1 = floor(mu0);
 if ( A1 < 1. ) A1 = 1.;
 L1 = 1.;
 A = 1.;
 for (j=1; j<=(int)A1; j++) {
   A = (double)j;
   L1 = cewma_L_arl(lambda, A, AU, mu0, z0, mu0, N);
   if ( L1 > L0 ) j = (int)A1 + 1; 
 }
 A1 = A;
 
 if ( L1 > L0 ) {
   for (j=1; j<=jmax; j++) {
     for (dA=1; dA<20; dA++) {
       A = A1 + dA / pow(-10., (double)j);
       L1 = cewma_L_arl(lambda, A, AU, mu0, z0, mu0, N);
       if ( ( fmod((double)j, 2.) > 0. && L1 < L0 ) || ( fmod((double)j, 2.) < 1. && L1 > L0 ) ) dA = 20;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 } else {
   for (j=1; j<=jmax; j++) {
     for (dA=1; dA<20; dA++) {
       A = A1 - dA / pow(-10., (double)j);
       L1 = cewma_L_arl(lambda, A, AU, mu0, z0, mu0, N);
       if ( ( fmod((double)j, 2.) < 1. && L1 < L0 ) || ( fmod((double)j, 2.) > 0. && L1 > L0 ) ) dA = 20;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 }
 
 return A;
}


double cewma_L_crit_new(double lambda, double L0, double AU, double mu0, double z0, int N, int jmax)
{ double A, A1, L1, ALmax;
  int j, dA; 

 ALmax = mu0 / sqrt( lambda*mu0 / (2.-lambda) ) - 0.00000000001;
       
 A1 = floor(mu0);
 if ( A1 < 1. ) A1 = 1.;
 if ( A1 > ALmax ) A1 = floor(ALmax);
 L1 = 1.;
 A = 1.;
 
 /*printf("\nALmax = %.4f,\tA1 = %.4f\n\n", ALmax, A1);*/
 
 for (j=1; j<=(int)A1; j++) {
   A = (double)j;   
   L1 = cewma_L_arl_new(lambda, A, AU, mu0, z0, mu0, N);
   /*printf("!!!A = %.4f,\tL1 = %.6f\n", A, L1);*/
   if ( L1 > L0 ) j = (int)A1 + 1; 
 }
 A1 = A;
 
 if ( L1 > L0 ) {
   for (j=1; j<=jmax; j++) {
     for (dA=1; dA<30; dA++) {
       A = A1 + dA / pow(-10., (double)j);
       if ( A > ALmax ) {
          A = ALmax - 1./pow(10., (double)j+1);
          dA = 30;
       }
       L1 = cewma_L_arl_new(lambda, A, AU, mu0, z0, mu0, N);
       /*printf("---A = %.4f,\tL1 = %.6f\n", A, L1);*/
       if ( ( fmod((double)j, 2.) > 0. && L1 < L0 ) || ( fmod((double)j, 2.) < 1. && L1 > L0 ) ) dA = 30;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 } else {
   for (j=1; j<=jmax; j++) {
     for (dA=1; dA<30; dA++) {
       A = A1 - dA / pow(-10., (double)j);
       if ( A > ALmax ) {
          A = ALmax - 1./pow(10., (double)j+1);
          dA = 30;
       }
       L1 = cewma_L_arl_new(lambda, A, AU, mu0, z0, mu0, N);
       /*printf("+++A = %.4f,\tL1 = %.6f\n", A, L1);*/
       if ( ( fmod((double)j, 2.) < 1. && L1 < L0 ) || ( fmod((double)j, 2.) > 0. && L1 > L0 ) ) dA = 30;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 }
 
 return A;
}


double cewma_2_crit_sym(double lambda, double L0, double mu0, double z0, int N, int jmax)
{ double A, A1, L1;
  int j, dA;
    
 A1 = floor(mu0);
 if ( A1 < 1. ) A1 = 1.;
 L1 = 1.;
 A = 1.;
 for (j=1; j<=(int)A1; j++) {
   A = (double)j;
   L1 = cewma_2_arl(lambda, A, A, mu0, z0, mu0, N);
   if ( L1 > L0 ) j = (int)A1 + 1; 
 }
 A1 = A;
 
 if ( L1 > L0 ) {
   for (j=1; j<=jmax; j++) {
     for (dA=1; dA<20; dA++) {
       A = A1 + dA / pow(-10., (double)j);
       L1 = cewma_2_arl(lambda, A, A, mu0, z0, mu0, N);
       if ( ( fmod((double)j, 2.) > 0. && L1 < L0 ) || ( fmod((double)j, 2.) < 1. && L1 > L0 ) ) dA = 20;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 } else {
   for (j=1; j<=jmax; j++) {
     for (dA=1; dA<20; dA++) {
       A = A1 - dA / pow(-10., (double)j);
       L1 = cewma_2_arl(lambda, A, A, mu0, z0, mu0, N);
       if ( ( fmod((double)j, 2.) < 1. && L1 < L0 ) || ( fmod((double)j, 2.) > 0. && L1 > L0 ) ) dA = 20;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 }
 
 return A;
}


double cewma_2_crit_sym_new(double lambda, double L0, double mu0, double z0, int N, int jmax)
{ double A, A1, L1;
  int j, dA;
  
 A1 = floor(mu0);
 if ( A1 < 1. ) A1 = 1.;
 L1 = 1.;
 A = 1.;
 for (j=1; j<=(int)A1; j++) {
   A = (double)j;
   L1 = cewma_2_arl_new(lambda, A, A, mu0, z0, mu0, N);
   /*printf("!!!A = %.4f,\tL1 = %.6f\n", A, L1);*/
   if ( L1 > L0 ) j = (int)A1 + 1; 
 }
 A1 = A;
 
 if ( L1 > L0 ) {
   for (j=0; j<=jmax; j++) {
     for (dA=1; dA<20; dA++) {
       A = A1 - dA / pow(-10., (double)j);
       L1 = cewma_2_arl_new(lambda, A, A, mu0, z0, mu0, N);
       /*printf("+++A = %.4f,\tL1 = %.6f\n", A, L1);*/
       if ( ( fmod((double)j, 2.) < 1. && L1 < L0 ) || ( fmod((double)j, 2.) > 0. && L1 > L0 ) ) dA = 20;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 } else {
   for (j=0; j<=jmax; j++) {
     for (dA=1; dA<20; dA++) {
       A = A1 + dA / pow(-10., (double)j);
       L1 = cewma_2_arl_new(lambda, A, A, mu0, z0, mu0, N);
       /*printf("---A = %.4f,\tL1 = %.6f\n", A, L1);*/
       if ( ( fmod((double)j, 2.) > 0. && L1 < L0 ) || ( fmod((double)j, 2.) < 1. && L1 > L0 ) ) dA = 20;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 }
 
 return A;
}


double cewma_2_crit_AL(double lambda, double L0, double AU, double mu0, double z0, int N, int jmax)
{ double A, A1, L1, ALmax;
  int j, dA;
  
 ALmax = mu0 / sqrt( lambda*mu0 / (2.-lambda) ) - 0.00000000001; 
    
 A1 = AU;
 L1 = cewma_2_arl(lambda, A1, AU, mu0, z0, mu0, N);

 if ( L1 > L0 ) {
   for (j=1; j<=jmax; j++) {
     for (dA=1; dA<30; dA++) {
       A = A1 + dA / pow(-10., (double)j);
       if ( A > ALmax ) {
         A = ALmax - 1./pow(10., (double)j+1);
         dA = 30;
       }
       L1 = cewma_2_arl(lambda, A, AU, mu0, z0, mu0, N);
       if ( ( fmod((double)j, 2.) > 0. && L1 < L0 ) || ( fmod((double)j, 2.) < 1. && L1 > L0 ) ) dA = 30;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 } else {
   for (j=1; j<=jmax; j++) {
     for (dA=1; dA<30; dA++) {
       A = A1 - dA / pow(-10., (double)j);
       if ( A > ALmax ) {
         A = ALmax - 1./pow(10., (double)j+1);
         dA = 30;
       }
       L1 = cewma_2_arl(lambda, A, AU, mu0, z0, mu0, N);
       if ( ( fmod((double)j, 2.) < 1. && L1 < L0 ) || ( fmod((double)j, 2.) > 0. && L1 > L0 ) ) dA = 30;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 }
 
 return A;
}


double cewma_2_crit_AL_new(double lambda, double L0, double AU, double mu0, double z0, int N, int jmax)
{ double A, A1, L1, ALmax;;
  int j, dA;

 ALmax = mu0 / sqrt( lambda*mu0 / (2.-lambda) ) - 0.00000000001;  
    
 A1 = AU;
 L1 = cewma_2_arl_new(lambda, A1, AU, mu0, z0, mu0, N);

 if ( L1 > L0 ) {
   for (j=1; j<=jmax; j++) {
     for (dA=1; dA<30; dA++) {
       A = A1 + dA / pow(-10., (double)j);
       if ( A > ALmax ) {
         A = ALmax - 1./pow(10., (double)j+1);
         dA = 30;
       }
       L1 = cewma_2_arl_new(lambda, A, AU, mu0, z0, mu0, N);
       if ( ( fmod((double)j, 2.) > 0. && L1 < L0 ) || ( fmod((double)j, 2.) < 1. && L1 > L0 ) ) dA = 30;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 } else {
   for (j=1; j<=jmax; j++) {
     for (dA=1; dA<30; dA++) {
       A = A1 - dA / pow(-10., (double)j);
       if ( A > ALmax ) {
         A = ALmax - 1./pow(10., (double)j+1);
         dA = 30;
       }
       L1 = cewma_2_arl_new(lambda, A, AU, mu0, z0, mu0, N);
       if ( ( fmod((double)j, 2.) < 1. && L1 < L0 ) || ( fmod((double)j, 2.) > 0. && L1 > L0 ) ) dA = 30;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 }
 
 return A;
}


double cewma_2_crit_AU(double lambda, double L0, double AL, double mu0, double z0, int N, int jmax)
{ double A, A1, L1;
  int j, dA;
    
 /*printf("\n\nGet AU for AL = %.4f and L0 = %.1f\n\n", AL, L0); */
  
 A1 = AL;
 L1 = cewma_2_arl(lambda, AL, A1, mu0, z0, mu0, N);
 /*printf("\n***A1 = %.4f,\tL1 = %.6f\n\n", A1, L1);*/

 if ( L1 > L0 ) {
   for (j=1; j<=jmax; j++) {
     for (dA=1; dA<20; dA++) {
       A = A1 + dA / pow(-10., (double)j);
       L1 = cewma_2_arl(lambda, AL, A, mu0, z0, mu0, N);
       /*printf("---A = %.4f,\tL1 = %.6f\n", A, L1);*/
       if ( ( fmod((double)j, 2.) > 0. && L1 > L0 ) || ( fmod((double)j, 2.) < 1. && L1 < L0 ) ) dA = 20;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 } else {
   for (j=1; j<=jmax; j++) {
     for (dA=1; dA<20; dA++) {
       A = A1 - dA / pow(-10., (double)j);
       L1 = cewma_2_arl(lambda, AL, A, mu0, z0, mu0, N);
       /*printf("+++A = %.4f,\tL1 = %.6f\n", A, L1);*/
       if ( ( fmod((double)j, 2.) < 1. && L1 < L0 ) || ( fmod((double)j, 2.) > 0. && L1 > L0 ) ) dA = 20;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 }
 
 return A;
}


double cewma_2_crit_AU_new(double lambda, double L0, double AL, double mu0, double z0, int N, int jmax)
{ double A, A1, L1;
  int j, dA;
  
 A1 = AL;
 L1 = cewma_2_arl_new(lambda, AL, A1, mu0, z0, mu0, N);
 /*printf("  A = %.4f,\tL1 = %.6f (L0=%.0f)\n\n", A1, L1, L0);*/

 if ( L1 < L0 ) {
   for (j=0; j<=jmax; j++) {
     for (dA=1; dA<20; dA++) {
       A = A1 + dA / pow(-10., (double)j);
       L1 = cewma_2_arl_new(lambda, AL, A, mu0, z0, mu0, N);
       /*printf("--A = %.4f,\tL1 = %.6f (L0=%.0f)\n", A, L1, L0);*/
       if ( ( fmod((double)j, 2.) < 1. && L1 > L0 ) || ( fmod((double)j, 2.) > 0. && L1 < L0 ) ) dA = 20;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 } else {
   for (j=0; j<=jmax; j++) {
     for (dA=1; dA<20; dA++) {
       A = A1 - dA / pow(-10., (double)j);
       L1 = cewma_2_arl_new(lambda, AL, A, mu0, z0, mu0, N);
       /*printf("++A = %.4f,\tL1 = %.6f (L0=%.0f)\n", A, L1, L0);*/
       if ( ( fmod((double)j, 2.) > 0. && L1 < L0 ) || ( fmod((double)j, 2.) < 1. && L1 > L0 ) ) dA = 20;
     }
     A1 = A;
   }
   if ( L1 < L0 ) A = A + pow(0.1, (double)jmax);
 }
 /*printf("\n");*/
 
 return A;
}


int cewma_2_crit_unb(double lambda, double L0, double mu0, double z0, int N, int jmax, double *AL, double *AU)
{ double A1, Lp, Lm, eps=.1, slope, lAL, lAU;
  int j, dA;
  
  A1 = cewma_2_crit_sym(lambda, L0, mu0, z0, N, jmax);
  Lp = cewma_2_arl(lambda, A1, A1, mu0, z0, mu0+eps, N);
  Lm = cewma_2_arl(lambda, A1, A1, mu0, z0, mu0-eps, N);
  slope = ( Lp - Lm ) / (2*eps);
  /*printf("\nA1 = %.4f,\tslope = %.6f\n\n", A1, slope);*/
  
  if ( slope > 0 ) {
    for (j=1; j<=jmax; j++) {
      for (dA=1; dA<20; dA++) {
        lAL = A1 - dA / pow(-10., (double)j);  
        lAU = cewma_2_crit_AU(lambda, L0, lAL, mu0, z0, N, jmax);
        Lp = cewma_2_arl(lambda, lAL, lAU, mu0, z0, mu0+eps, N);
        Lm = cewma_2_arl(lambda, lAL, lAU, mu0, z0, mu0-eps, N);
        slope = ( Lp - Lm ) / (2*eps);
        /*printf("--AL = %.4f,\tAU = %.4f,\tslope = %.6f\n", lAL, lAU, slope);*/
        if ( ( fmod((double)j, 2.) > 0. && slope < 0. ) || ( fmod((double)j, 2.) < 1. && slope > 0. ) ) dA = 20;
      }
      A1 = lAL;
    }
  } else {
    for (j=1; j<=jmax; j++) {
      for (dA=1; dA<20; dA++) {
        lAL = A1 + dA / pow(-10., (double)j);  
        lAU = cewma_2_crit_AU(lambda, L0, lAL, mu0, z0, N, jmax);
        Lp = cewma_2_arl(lambda, lAL, lAU, mu0, z0, mu0+eps, N);
        Lm = cewma_2_arl(lambda, lAL, lAU, mu0, z0, mu0-eps, N);
        slope = ( Lp - Lm ) / (2*eps);
        /*printf("++AL = %.4f,\tAU = %.4f,\tslope = %.6f\n", lAL, lAU, slope);*/
        if ( ( fmod((double)j, 2.) < 1. && slope < 0. ) || ( fmod((double)j, 2.) > 0. && slope > 0. ) ) dA = 20;
      }
      A1 = lAL;
    }  
  }
  
  *AL = lAL;
  *AU = lAU;
  
  return 0;
}


int cewma_2_crit_unb_new(double lambda, double L0, double mu0, double z0, int N, int jmax, double *AL, double *AU)
{ double A1, symA, Lp, Lm, eps=.01, slope, lAL, lAU, ALmin, AUlarge=10.;
  int j, dA;
  
  symA = cewma_2_crit_sym_new(lambda, L0, mu0, z0, N, jmax);
  Lp = cewma_2_arl_new(lambda, symA, symA, mu0, z0, mu0+eps, N);
  Lm = cewma_2_arl_new(lambda, symA, symA, mu0, z0, mu0-eps, N);
  slope = ( Lp - Lm ) / (2*eps);  
  ALmin = cewma_L_crit_new(lambda, L0, AUlarge, mu0, z0, N, jmax);
  /*printf("\nsymA = %.4f,\tslope = %.6f,\tALmin = %.4f\n\n", symA, slope, ALmin);*/
  
  lAL = symA;
  lAU = symA;
  A1 = symA;
  
  if ( slope > 0 ) {
    for (j=0; j<=jmax; j++) {
      for (dA=1; dA<30; dA++) {
        lAL = A1 + dA / pow(-10., (double)j);
        if ( lAL < ALmin ) {
          lAL = ALmin + 1./pow(10., (double)j+1);
          dA = 30;
        }
        lAU = cewma_2_crit_AU_new(lambda, L0, lAL, mu0, z0, N, jmax);
        Lp = cewma_2_arl_new(lambda, lAL, lAU, mu0, z0, mu0+eps, N);
        Lm = cewma_2_arl_new(lambda, lAL, lAU, mu0, z0, mu0-eps, N);
        slope = ( Lp - Lm ) / (2*eps);
        /*printf("--AL = %.4f,\tAU = %.4f,\tslope = %.6f\n", lAL, lAU, slope);*/
        if ( ( fmod((double)j, 2.) < 1. && slope < 0. ) || ( fmod((double)j, 2.) > 0. && slope > 0. ) ) dA = 30;
      }
      A1 = lAL;
    }
  } else {
    for (j=0; j<=jmax; j++) {
      /*printf("!!j = %d\n", j);*/
      for (dA=1; dA<30; dA++) {
        lAL = A1 - dA / pow(-10., (double)j);
        if ( lAL < ALmin ) {
          lAL = ALmin + 1./pow(10., (double)j+1);
          dA = 30;
        } else {
          if ( lAL > symA ) {
            /*printf("\nlAL (%.6f) too large\n", lAL);*/
            lAL = symA - 1./pow(10., (double)j+1);
            dA = 30;
          }  
        }
        lAU = cewma_2_crit_AU_new(lambda, L0, lAL, mu0, z0, N, jmax);
        Lp = cewma_2_arl_new(lambda, lAL, lAU, mu0, z0, mu0+eps, N);
        Lm = cewma_2_arl_new(lambda, lAL, lAU, mu0, z0, mu0-eps, N);
        slope = ( Lp - Lm ) / (2*eps);
        /*printf("++AL = %.4f,\tAU = %.4f,\tslope = %.6f\n", lAL, lAU, slope);*/
        if ( ( fmod((double)j, 2.) > 0. && slope < 0. ) || ( fmod((double)j, 2.) < 1. && slope > 0. ) ) dA = 30;
      }
      A1 = lAL;
    }  
  }
  
  *AL = lAL;
  *AU = lAU;
  
  return 0;
}


double cewma_2_get_gL(double lambda, double L0, double mu0, double z0, double AL, double AU, double gU, int N)
{ double gL1, gL2, gL3, L1, L2, L3, dL, dG;
 
 gL1 = 1.;   
 L1 = cewma_2_arl_rando_new(lambda, AL, AU, gL1, gU, mu0, z0, mu0, N);   
 while ( L1 < L0 ) {
   gL2 = gL1;
   L2 = L1;
   gL1 /= 2.;
   L1 = cewma_2_arl_rando_new(lambda, AL, AU, gL1, gU, mu0, z0, mu0, N);
 }
 
 do {
   gL3 = gL1 + (L0 - L1)/(L2 - L1) * (gL2-gL1);
   L3 = cewma_2_arl_rando_new(lambda, AL, AU, gL3, gU, mu0, z0, mu0, N);
   dG = gL3 - gL2; dL = L0 - L3;
   gL1 = gL2; L1 = L2; gL2 = gL3; L2 = L3;
 } while ( fabs(dL)>.00000000001 && fabs(dG)>.00000000001 );
 
 return gL3;
}


double cewma_2_get_gU(double lambda, double L0, double mu0, double z0, double AL, double AU, double gL, int N)
{ double gU1, gU2, gU3, L1, L2, L3, dL, dG;
 
 gU1 = 1.;   
 L1 = cewma_2_arl_rando_new(lambda, AL, AU, gL, gU1, mu0, z0, mu0, N);
 while ( L1 < L0 ) {
   gU2 = gU1;
   L2 = L1;
   gU1 /= 2.;
   L1 = cewma_2_arl_rando_new(lambda, AL, AU, gL, gU1, mu0, z0, mu0, N);
 }
 
 do {
   gU3 = gU1 + (L0 - L1)/(L2 - L1) * (gU2-gU1);
   L3 = cewma_2_arl_rando_new(lambda, AL, AU, gL, gU3, mu0, z0, mu0, N);
   dG = gU3 - gU2; dL = L0 - L3;
   gU1 = gU2; L1 = L2; gU2 = gU3; L2 = L3;
 } while ( fabs(dL)>.00000000001 && fabs(dG)>.00000000001 );
 
 return gU3;
}


int cewma_2_crit_unb_rando_new(double lambda, double L0, double mu0, double z0, int N, int jmax, double *AL, double *AU, double *gL, double *gU)
{ double A1, symA, Lp, Lm, eps=.01, slope, lAL1, lAL2, lAU1, lAU2, lAL, lAU, g1, g2, g3, u1, u2, u3, s1, s2, s3, miAL, maAL, miAU, maAU, L0act, dg, rdA, ALmin, AUlarge=10.;
  int j, dA, done=0;
  
  symA = cewma_2_crit_sym_new(lambda, L0, mu0, z0, N, jmax);
  Lp = cewma_2_arl_new(lambda, symA, symA, mu0, z0, mu0+eps, N);
  Lm = cewma_2_arl_new(lambda, symA, symA, mu0, z0, mu0-eps, N);
  slope = ( Lp - Lm ) / (2*eps);  
  ALmin = cewma_L_crit_new(lambda, L0, AUlarge, mu0, z0, N, jmax);
  /*printf("\nsymA = %.4f,\tslope = %.6f,\tALmin = %.4f\n\n", symA, slope, ALmin);*/
  
  lAL1 = symA;
  lAU1 = symA;
  A1 = symA;
  
  if ( slope > 0 ) {
    for (j=0; j<=jmax; j++) {
      for (dA=1; dA<30; dA++) {
        lAL2 = lAL1;
        lAU2 = lAU1;
        lAL1 = A1 + dA / pow(-10., (double)j);      
        if ( lAL1 < ALmin ) {
          lAL1 = ALmin + 1./pow(10., (double)j+1);
          dA = 30;
        }
        lAU1 = cewma_2_crit_AU_new(lambda, L0, lAL1, mu0, z0, N, jmax);
        Lp = cewma_2_arl_new(lambda, lAL1, lAU1, mu0, z0, mu0+eps, N);
        Lm = cewma_2_arl_new(lambda, lAL1, lAU1, mu0, z0, mu0-eps, N);
        slope = ( Lp - Lm ) / (2*eps);
        /*printf("--AL = %.4f,\tAU = %.4f,\tslope = %.6f\n", lAL1, lAU1, slope);*/
        if ( ( fmod((double)j, 2.) < 1. && slope < 0. ) || ( fmod((double)j, 2.) > 0. && slope > 0. ) ) dA = 30;
      }
      A1 = lAL1;
    }
  } else {
    for (j=0; j<=jmax; j++) {
      /*printf("j = %d\n", j);*/
      for (dA=1; dA<30; dA++) {
        lAL2 = lAL1;
        lAU2 = lAU1;
        lAL1 = A1 - dA / pow(-10., (double)j); 
        
        if ( lAL1 < ALmin ) {
          /*printf("\nlAL1 too small\n");*/
          lAL1 = ALmin + 1./pow(10., (double)j+1);
          dA = 30;
        } else {
          if ( lAL1 > symA ) {
            /*printf("\nlAL1 (%.6f) too large\n", lAL1);*/
            lAL1 = symA - 1./pow(10., (double)j+1);
            dA = 30;
          }                    
        }
        
        lAU1 = cewma_2_crit_AU_new(lambda, L0, lAL1, mu0, z0, N, jmax);
        Lp = cewma_2_arl_new(lambda, lAL1, lAU1, mu0, z0, mu0+eps, N);
        Lm = cewma_2_arl_new(lambda, lAL1, lAU1, mu0, z0, mu0-eps, N);
        slope = ( Lp - Lm ) / (2*eps);
        /*printf("++AL = %.4f,\tAU = %.4f,\tslope = %.6f\n", lAL1, lAU1, slope);*/
        if ( ( fmod((double)j, 2.) > 0. && slope < 0. ) || ( fmod((double)j, 2.) < 1. && slope > 0. ) ) dA = 30;
      }
      A1 = lAL1;
    }  
  }
  
  L0act = cewma_2_arl_new(lambda, lAL1, lAU1, mu0, z0, mu0, N);
  /*printf("\n\nlAL1 = %.8f,\tlAU1 = %.8f,\tslope1 = %.6f,\tL0act = %.3f\n", lAL1, lAU1, slope, L0act);
  printf("lAL2 = %.8f,\tlAU2 = %.8f\n", lAL2, lAU2);*/
  rdA = pow(10.,-(double)jmax);
  /*printf("rdA = %.8f\n\n", rdA);*/
  
  miAL = lAL1; maAL = lAL2;
  if ( lAL2 < miAL ) { miAL = lAL2; maAL = lAL1; }
  miAU = lAU1; maAU = lAU2;
  if ( lAU2 < miAU ) { miAU = lAU2; maAU = lAU1; }
  if ( (maAU - miAU)/rdA > 100. ) maAU += 20.*rdA;
  if ( (maAU - miAU)/rdA > 1000. ) maAU += 200.*rdA;
  
  for ( lAL=miAL; lAL<=maAL+rdA/10.; lAL+=rdA ) {
    /*printf("lAL = %.8f\n", lAL);*/
    miAU = cewma_2_crit_AU_new(lambda, L0, lAL, mu0, z0, N, jmax);
    /*for ( lAU=miAU; lAU<=maAU+rdA/10.; lAU+=rdA ) {*/
    for ( lAU=maAU; lAU>=miAU-rdA/10.; lAU-=rdA ) {
    /*for ( lAU=maAU; lAU>=miAU-rdA; lAU-=rdA ) {*/
      /*printf("lAU = %.8f\n", lAU);*/
      L0act = cewma_2_arl_rando_new(lambda, lAL, lAU, 0., 0., mu0, z0, mu0, N);
      if ( L0act < L0 ) {
        /*printf("L0act too small\n");*/
        done = 0;
      } else {
        L0act = cewma_2_arl_rando_new(lambda, lAL, lAU, 1., 1., mu0, z0, mu0, N);
        if ( L0act > L0 ) {
          /*printf("L0act too large\n");*/
          done = 0;
        } else {         
          g1 = 0.;
          L0act = cewma_2_arl_rando_new(lambda, lAL, lAU, g1, 1., mu0, z0, mu0, N);
          if ( L0act < L0 ) {
            /*printf("\nfull gL interval\n");*/
            u1 = cewma_2_get_gU(lambda, L0, mu0, z0, lAL, lAU, g1, N);
            Lp = cewma_2_arl_rando_new(lambda, lAL, lAU, g1, u1, mu0, z0, mu0+eps, N);
            Lm = cewma_2_arl_rando_new(lambda, lAL, lAU, g1, u1, mu0, z0, mu0-eps, N);
            s1 = ( Lp - Lm ) / (2*eps);
            L0act = cewma_2_arl_rando_new(lambda, lAL, lAU, g1, u1, mu0, z0, mu0, N);
            /*printf("g1 = %.6f,\tu1 = %.6f,\tL0act = %.3f,\ts1 = %.6f\n\n", g1, u1, L0act, s1);*/
          } else {
            /*printf("\nsearch for min gL\n");*/
            u1 = 1.;
            g1 = cewma_2_get_gL(lambda, L0, mu0, z0, lAL, lAU, u1, N);
            Lp = cewma_2_arl_rando_new(lambda, lAL, lAU, g1, u1, mu0, z0, mu0+eps, N);
            Lm = cewma_2_arl_rando_new(lambda, lAL, lAU, g1, u1, mu0, z0, mu0-eps, N);
            s1 = ( Lp - Lm ) / (2*eps);
            L0act = cewma_2_arl_rando_new(lambda, lAL, lAU, g1, u1, mu0, z0, mu0, N);
            /*printf("g1 = %.6f,\tu1 = %.6f,\tL0act = %.3f,\ts1 = %.6f\n\n", g1, u1, L0act, s1);*/
          }
        
          u2 = 0.;
          L0act = cewma_2_arl_rando_new(lambda, lAL, lAU, 1., u2, mu0, z0, mu0, N);
          if ( L0act < L0 ) {
            /*printf("\nfull gU interval\n");*/
            g2 = cewma_2_get_gL(lambda, L0, mu0, z0, lAL, lAU, u2, N);
            Lp = cewma_2_arl_rando_new(lambda, lAL, lAU, g2, u2, mu0, z0, mu0+eps, N);
            Lm = cewma_2_arl_rando_new(lambda, lAL, lAU, g2, u2, mu0, z0, mu0-eps, N);
            s2 = ( Lp - Lm ) / (2*eps);
            L0act = cewma_2_arl_rando_new(lambda, lAL, lAU, g2, u2, mu0, z0, mu0, N);
            /*printf("g2 = %.6f,\tu2 = %.6f,\tL0act = %.3f,\ts2 = %.6f\n\n", g2, u2, L0act, s2);*/
          } else {
            /*printf("\nsearch for min gU\n");*/
            g2 = 1.;
            u2 = cewma_2_get_gU(lambda, L0, mu0, z0, lAL, lAU, g2, N);
            Lp = cewma_2_arl_rando_new(lambda, lAL, lAU, g2, u2, mu0, z0, mu0+eps, N);
            Lm = cewma_2_arl_rando_new(lambda, lAL, lAU, g2, u2, mu0, z0, mu0-eps, N);
            s2 = ( Lp - Lm ) / (2*eps);
            L0act = cewma_2_arl_rando_new(lambda, lAL, lAU, g2, u2, mu0, z0, mu0, N);
            /*printf("g2 = %.6f,\tu2 = %.6f,\tL0act = %.3f,\ts2 = %.6f\n\n", g2, u2, L0act, s2);*/
          } 
      
          if ( s1*s2 > 0. ) {
            /*printf("slope does not change sign\n");*/
            done = 0;
          } else {
            do {
              g3 = g1 + (0. - s1)/(s2 - s1) * (g2-g1);
              u3 = cewma_2_get_gU(lambda, L0, mu0, z0, lAL, lAU, g3, N);
              L0act = cewma_2_arl_rando_new(lambda, lAL, lAU, g3, u3, mu0, z0, mu0, N);
              Lp = cewma_2_arl_rando_new(lambda, lAL, lAU, g3, u3, mu0, z0, mu0+eps, N);
              Lm = cewma_2_arl_rando_new(lambda, lAL, lAU, g3, u3, mu0, z0, mu0-eps, N);
              s3 = ( Lp - Lm ) / (2*eps);
              /*printf("g3 = %.6f,\tu3 = %.6f,\tL0act = %.3f,\ts3 = %.6f\n", g3, u3, L0act, s3);*/
              dg = g3 - g2;
              g1 = g2; s1 = s2; g2 = g3; s2 = s3;
            } while ( fabs(s3)>.00000000001 && fabs(dg)>.00000000001 );
            done = 1;
            break;
          } /* s1*s2 > 0 */
        } /* L0act too large */
      } /* L0act too small */
    } /* lAU=miAU; lAU<=maAU; lAU+=rdA */
    if ( done == 1 ) break;
  }
  
  *AL = lAL;
  *AU = lAU;
  *gL = g3;
  *gU = u3;
  
  return 0;
}


/* TEWMA stuff */

double tewma_arl(double lambda, int k, int lk, int uk, double z0, double mu)
{ double *a, *g, *DM, *DELL, *F2, arl, pij, term;
  int i, il, j, jl, m, km, l, l0, l1, M, N, kM;
  
 N = uk - lk + 1;
 a = matrix(N, N);
 g = vector(N);
 
 M = (int)qf_pois( 1.-1e-15, mu);
 DM = vector(M+1);
 kM = k*M;
 F2 = matrix(M+1, kM+1);
 /*for (i=0; i<=M; i++) for (j=0; j<=kM; j++) F2[i*kM+j] = 0.;*/
 for (i=0; i<=M; i++) {
   DM[i] = pdf_pois((double)i, mu);
   for (j=0; j<=k*i; j++) F2[i*kM+j] = pdf_binom((double)j, k*i, lambda);
 }
 
 DELL = vector(uk+1);
 
 for (i=0; i<N; i++) for (j=0; j<N; j++) a[i*N+j] = 0.;
 
 for (i=0; i<N; i++) {
   il = lk + i;
   for (l=0; l<=il; l++) DELL[l] = pdf_binom(l, il, 1.-lambda);
   for (j=0; j<N; j++) {
     jl = lk + j;
     l1 = jl;
     if ( l1 > il ) l1 = il;
     pij = 0.;
     for (m=0; m<=M; m++) {
       km = k*m;  
       l0 = jl - km;
       if ( l0 < 0 ) l0 = 0;
       if ( l0 <= l1 ) {
         term = 0.;  
         /*for (l=l0; l<=l1; l++) term += pdf_binom((double)(jl-l), km, lambda) * DELL[l];*/
         for (l=l0; l<=l1; l++) term += F2[m*kM+jl-l] * DELL[l];         
         term *= DM[m];
       }
       pij += term;
     } /* loop over m */
     a[j*N+i] = - pij;
   } /* loop over j */
   ++a[i*N+i];
 } /* loop over i */      
 
 for (j=0; j<N; j++) g[j] = 1.;
 /*LU_solve(a, g, N);*/
 solve(&N, a, g);

 arl = g[ (int)round(z0) - lk ];
  
 Free(F2);
 Free(DELL);
 Free(DM); 
 Free(a);
 Free(g);

 return arl;
}    


double tewma_arl_R(double lambda, int k, int lk, int uk, double gl, double gu, double z0, double mu)
{ double *a, *g, *DM, *DELL, *F2, arl, pij, term;
  int i, il, j, jl, m, km, l, l0, l1, M, N, kM;
  
 N = uk - lk + 1;
 a = matrix(N, N);
 g = vector(N);
 
 M = (int)qf_pois( 1.-1e-15, mu);
 DM = vector(M+1);
 kM = k*M;
 F2 = matrix(M+1, kM+1);
 /*for (i=0; i<=M; i++) for (j=0; j<=kM; j++) F2[i*kM+j] = 0.;*/
 for (i=0; i<=M; i++) {
   DM[i] = pdf_pois((double)i, mu);
   for (j=0; j<=k*i; j++) F2[i*kM+j] = pdf_binom((double)j, k*i, lambda);
 }
 
 DELL = vector(uk+1);
 
 for (i=0; i<N; i++) for (j=0; j<N; j++) a[i*N+j] = 0.;
 
 for (i=0; i<N; i++) {
   il = lk + i;
   for (l=0; l<=il; l++) DELL[l] = pdf_binom(l, il, 1.-lambda);
   for (j=0; j<N; j++) {
     jl = lk + j;
     l1 = jl;
     if ( l1 > il ) l1 = il;
     pij = 0.;
     for (m=0; m<=M; m++) {
       km = k*m;  
       l0 = jl - km;
       if ( l0 < 0 ) l0 = 0;
       if ( l0 <= l1 ) {
         term = 0.;  
         /*for (l=l0; l<=l1; l++) term += pdf_binom((double)(jl-l), km, lambda) * DELL[l];*/
         for (l=l0; l<=l1; l++) term += F2[m*kM+jl-l] * DELL[l];         
         term *= DM[m];
       }
       pij += term;
     } /* loop over m */
     if ( j == 0 )   pij *= 1. - gl;
     if ( j == N-1 ) pij *= 1. - gu;
     a[j*N+i] = - pij;
   } /* loop over j */
   ++a[i*N+i];
 } /* loop over i */

/* for (i=0; i<N; i++) {
   printf("Zeile i = %d\n", i);
   for (j=0; j<N; j++) printf("\tSpalte j = %d,\tElement = %.9f\n", j, a[j*N+i]);
   printf("\n");
 }*/
 
 for (j=0; j<N; j++) g[j] = 1.;
 /*LU_solve(a, g, N);*/
 solve(&N, a, g);

 arl = g[ (int)round(z0) - lk ];

/* printf("\n");
 for (i=0; i<N; i++) printf("Zeile i = %d,\tarl = %.9f\n", i, g[i]);
 printf("\n");*/
  
 Free(F2);
 Free(DELL);
 Free(DM); 
 Free(a);
 Free(g);

 return arl;
} 


double eewma_arl(int gX, int gY, int kL, int kU, double mu, double y0, int r0)
{ double *a, *g, arl, *DELL; 
  int amin, amax, bmin, bmax, hb, N, xmin, xmax, i, j, x, y, o, i0, M;
 
 bmin = gY * kL;
 bmax = gX + gY*(kU+1) - 1;
 amin = (gX+gY)*kL;
 amax = (gX+gY)*(kU+1) - 1;
 hb = bmax - bmin;
 N = hb+1;
 
 a = matrix(N, N);
 g = vector(N);
 
 for (i=0; i<N; i++) for (j=0; j<N; j++) a[i*N+j] = 0.;
 
 M = (int)ceil( (amax-bmin)/(double)gX );
 DELL = vector(M+1);
 for (i=0; i<=M; i++) DELL[i] = -pdf_pois((double)i, mu);
 
 for (i=0; i<=hb; i++) {
   xmin = (int)ceil( (amin-i-bmin)/(double)gX );
   xmax = (int)floor( (amax-i-bmin)/(double)gX );
   for (x=xmin; x<=xmax; x++) {
     o = gX*x + i + bmin;
     y = (int)floor( o/(double)(gX+gY) );
     j = o - gX*y - bmin;
     a[j*N+i] += -pdf_pois((double)x, mu);
     /*a[j*N+i] += DELL[x];*/
   }
 }
 
 for (i=0; i<N; i++) {
   g[i] = 1.;
   ++a[i*N+i];
 }
 
 /*LU_solve(a, g, N);*/
 solve(&N, a, g); 
 
 i0 = gY*(int)floor(y0) - bmin + r0;
 
 arl = g[i0];
 
 Free(DELL);
 Free(g);
 Free(a);
 
 return arl;
}

/* HIER */
double ccusum_U_arl(double mu, int km, int hm, int m, int i0)
{ double *a, *b1, *b2, *x, *y, *z, *phi, *psi, *g, px, al, ga, et, de, be, arl;
  int i, j, l, lmax, N, N1;
 
 N = hm + 1; 
 N1 = N - 1;
 
 a   = vector(2*N-1);
 b1  = vector(N);
 b2  = vector(N);
 x   = vector(N);
 y   = vector(N);
 z   = vector(N);
 phi = vector(N);
 psi = vector(N);
 g   = vector(N);

 lmax = (int)ceil( (hm + km) / m ) + 1;
 
 for (l=0; l<=lmax; l++) {
   px = pdf_pois((double)l, mu);
   i = km - l*m;   
   if ( 0 <= N+i-1 && N+i-1 < 2*N-1 ) a[N+i-1] = -px;
   if ( 0 <= i-1 && i-1 < N ) b2[i-1] = px;
 }
 a[N1] += 1.; 
 
 for (i=N-1; i>=0; i--) {
   b1[i] = 1.;
   if ( i > 0 ) b2[i-1] += b2[i];
 }

 x[0] = 1./a[N1];
 y[0] = 1./a[N1];
 phi[0] = b1[0]/a[N1];
 psi[0] = b2[0]/a[N1];
 
 for (i=1; i<N; i++) {
   al = 0.;
   for (j=0; j<i; j++) al += a[N1 + i - j] * x[j];
   ga = 0.;
   for (j=0; j<i; j++) ga += a[N1 - 1 - j] * y[j];
   et = -b1[i];
   for (j=0; j<i; j++) et += a[N1 + i - j] * phi[j];
   de = -b2[i];
   for (j=0; j<i; j++) de += a[N1 + i - j] * psi[j];

   be = 1. - al*ga;
   
   z[0] = -ga*x[0] / be;
   for (j=1; j<i; j++) z[j] = ( y[j-1] - ga*x[j] ) / be;
   z[i] = y[i-1] / be;   

   x[0] = x[0] / be;
   for (j=1; j<i; j++) x[j] = ( x[j] - al*y[j-1] ) / be;
   x[i] = -al*y[i-1] / be;
   
   for (j=0; j<=i; j++) y[j] = z[j];
       
   for (j=0; j<i; j++) {
     phi[j] = ( phi[j] - et*z[j] );
     psi[j] = ( psi[j] - de*z[j] );
   }
   phi[i] = -et*z[i];
   psi[i] = -de*z[i];
 }

 be = phi[0] / ( 1. - psi[0] );
 for (i=0; i<N; i++) g[i] = phi[i] + psi[i] * be;
 arl = g[i0];

 Free(g);
 Free(psi);
 Free(phi);
 Free(z);
 Free(y);
 Free(x);
 Free(b2);
 Free(b1);
 Free(a);

 return arl;  
}

int ccusum_U_crit(double A, double mu0, int km, int m, int i0)
{ int p10, hm, p;
  double L1;
    
  p10 = (int)ceil( log10((double)m) );
  hm = 2 * (int)qf_pois(1.-1./A, mu0) * m; /* twice the Shewhart critical value */
  L1 = ccusum_U_arl(mu0, km, hm, m, i0);
  /*printf("\n(*)hm = %d,\tL1 = %.3f\n", hm, L1);*/
    
  for (p=p10; p>=0; p--) {
    if ( L1 < A ) {
      while ( L1 < A ) {
        hm += pow(10., (double)p);
        L1 = ccusum_U_arl(mu0, km, hm, m, i0);
        /*printf("(+)\thm = %d,\tL1 = %.3f\n", hm, L1);*/
      }
    } else {
      while ( L1 >= A ) {
        hm -= pow(10., (double)p);
        if ( hm < km ) {
          p--;
          hm += pow(10., (double)p+1.) - pow(10., (double)p);
        }    
        L1 = ccusum_U_arl(mu0, km, hm, m, i0);
        /*printf("(-)\thm = %d,\tL1 = %.3f\n", hm, L1);*/
      }
    }
  }
  if ( L1 < A ) hm = hm + 1;      
  return hm;
}

double ccusum_U_arl_rando(double mu, int km, int hm, int m, double gamma, int i0)
{ double *a, *b1, *b2, *b3, *x, *y, *z, *phi, *psi, *xi, *rr, *g, *gx, px, al, ga, et, de, be, ka, lambda, arl, nen, zae;
  int i, j, l, lmax, N, N1;
  
 N  = hm; 
 N1 = N - 1;
 
 a   = vector(2*N-1);
 b1  = vector(N);
 b2  = vector(N);
 b3  = vector(N);
 x   = vector(N);
 y   = vector(N);
 z   = vector(N);
 phi = vector(N);
 psi = vector(N);
 xi  = vector(N);
 rr  = vector(N);
 g   = vector(N);
 gx  = vector(N);

 lmax = (int)ceil( (hm + km) / m ) + 1;
 
 for (l=0; l<=lmax; l++) {
   px = pdf_pois((double)l, mu);
   i = km - l*m;
   if ( 0 <= N+i-1 && N+i-1 < 2*N-1 ) a[N+i-1] = -px;
   if ( 0 <= i-1 && i-1 < N ) {
    b2[i-1] = px;
    rr[N-i] = px;
   }
   if (0 <= N+i && N+i < N ) b3[N+i] = (1.-gamma)*px;
 }
 a[N1] += 1.; 
 
 for (i=N-1; i>=0; i--) {
   b1[i] = 1.;
   if ( i > 0 ) b2[i-1] += b2[i];
 }

 x[0] = 1./a[N1];
 y[0] = 1./a[N1];
 phi[0] = b1[0]/a[N1];
 psi[0] = b2[0]/a[N1];
 xi[0] = b3[0]/a[N1];
 
 for (i=1; i<N; i++) {
   al = 0.;
   for (j=0; j<i; j++) al += a[N1 + i - j] * x[j];
   ga = 0.;
   for (j=0; j<i; j++) ga += a[N1 - 1 - j] * y[j];
   et = -b1[i];
   for (j=0; j<i; j++) et += a[N1 + i - j] * phi[j];
   de = -b2[i];
   for (j=0; j<i; j++) de += a[N1 + i - j] * psi[j];
   ka = -b3[i];
   for (j=0; j<i; j++) ka += a[N1 + i - j] * xi[j];

   be = 1. - al*ga;
   
   z[0] = -ga*x[0] / be;
   for (j=1; j<i; j++) z[j] = ( y[j-1] - ga*x[j] ) / be;
   z[i] = y[i-1] / be;   

   x[0] = x[0] / be;
   for (j=1; j<i; j++) x[j] = ( x[j] - al*y[j-1] ) / be;
   x[i] = -al*y[i-1] / be;
   
   for (j=0; j<=i; j++) y[j] = z[j];
       
   for (j=0; j<i; j++) {
     phi[j] = ( phi[j] - et*z[j] );
     psi[j] = ( psi[j] - de*z[j] );
     xi[j]  = (  xi[j] - ka*z[j] );
   }
   phi[i] = -et*z[i];
   psi[i] = -de*z[i];
   xi[i]  = -ka*z[i];
 }

 be = phi[0] / ( 1. - psi[0] );
 for (i=0; i<N; i++) g[i] = phi[i] + psi[i] * be;
 
 ka = xi[0] / ( 1. - psi[0] );
 for (i=0; i<N; i++) gx[i] = xi[i] + psi[i] * ka;
 
 nen = 0.;
 zae = 0.;
 for (i=0; i<N; i++) {
   zae += rr[i] * g[i];
   nen += rr[i] * gx[i];
 }
 lambda = ( 1. + zae ) / ( 1. - (1.-gamma)*(1.-a[N1]) - nen );

 for (i=0; i<N; i++) g[i] += lambda*gx[i];

 arl = g[i0];
 
 Free(gx);
 Free(g);
 Free(rr);
 Free(xi);
 Free(psi);
 Free(phi);
 Free(z);
 Free(y);
 Free(x);
 Free(b3);
 Free(b2);
 Free(b1);
 Free(a);

 return arl;  
}

int ccusum_U_rando_crit(double A, double mu0, int km, int m, int i0, int *hm, double *gamma)
{ double *a, *b1, *b2, *b3, *x, *y, *z, *phi, *psi, *xi, *rr, *g, *gx, px, al, ga, et, de, be, ka, lambda, L1, L2, L3, nen, zae, g1, g2, g3;
  int i, j, l, lmax, N, N1, ihm;
  
 ihm = ccusum_U_crit(A, mu0, km, m, i0);
  
 N  = ihm; 
 N1 = N - 1;
 
 a   = vector(2*N-1);
 b1  = vector(N);
 b2  = vector(N);
 b3  = vector(N);
 x   = vector(N);
 y   = vector(N);
 z   = vector(N);
 phi = vector(N);
 psi = vector(N);
 xi  = vector(N);
 rr  = vector(N);
 g   = vector(N);
 gx  = vector(N);

 lmax = (int)ceil( (ihm + km) / m ) + 1;
 
 for (l=0; l<=lmax; l++) {
   px = pdf_pois((double)l, mu0);
   i = km - l*m;
   if ( 0 <= N+i-1 && N+i-1 < 2*N-1 ) a[N+i-1] = -px;
   if ( 0 <= i-1 && i-1 < N ) {
    b2[i-1] = px;
    rr[N-i] = px;
   }
   if (0 <= N+i && N+i < N ) b3[N+i] = px;
 }
 a[N1] += 1.; 
 
 for (i=N-1; i>=0; i--) {
   b1[i] = 1.;
   if ( i > 0 ) b2[i-1] += b2[i];
 }

 x[0] = 1./a[N1];
 y[0] = 1./a[N1];
 phi[0] = b1[0]/a[N1];
 psi[0] = b2[0]/a[N1];
 xi[0] = b3[0]/a[N1];
 
 for (i=1; i<N; i++) {
   al = 0.;
   for (j=0; j<i; j++) al += a[N1 + i - j] * x[j];
   ga = 0.;
   for (j=0; j<i; j++) ga += a[N1 - 1 - j] * y[j];
   et = -b1[i];
   for (j=0; j<i; j++) et += a[N1 + i - j] * phi[j];
   de = -b2[i];
   for (j=0; j<i; j++) de += a[N1 + i - j] * psi[j];
   ka = -b3[i];
   for (j=0; j<i; j++) ka += a[N1 + i - j] * xi[j];

   be = 1. - al*ga;
   
   z[0] = -ga*x[0] / be;
   for (j=1; j<i; j++) z[j] = ( y[j-1] - ga*x[j] ) / be;
   z[i] = y[i-1] / be;   

   x[0] = x[0] / be;
   for (j=1; j<i; j++) x[j] = ( x[j] - al*y[j-1] ) / be;
   x[i] = -al*y[i-1] / be;
   
   for (j=0; j<=i; j++) y[j] = z[j];
       
   for (j=0; j<i; j++) {
     phi[j] = ( phi[j] - et*z[j] );
     psi[j] = ( psi[j] - de*z[j] );
     xi[j]  = (  xi[j] - ka*z[j] );
   }
   phi[i] = -et*z[i];
   psi[i] = -de*z[i];
   xi[i]  = -ka*z[i];
 }

 be = phi[0] / ( 1. - psi[0] );
 for (i=0; i<N; i++) g[i] = phi[i] + psi[i] * be;
 
 ka = xi[0] / ( 1. - psi[0] );
 for (i=0; i<N; i++) gx[i] = xi[i] + psi[i] * ka;
 
 nen = 0.;
 zae = 0.;
 for (i=0; i<N; i++) {
  zae += rr[i] * g[i];
  nen += rr[i] * gx[i];
 } 

 /* secant rule etc. */
 
 g1 = 1.;
 lambda = ( 1. + zae ) / ( 1. - (1.-a[N1]) - nen );
 L1 = g[i0] + lambda*gx[i0]; 
 /*printf("g1 = %.6f,\tL1 = %.6f,\tA = %.0f\n", g1, L1, A);*/
 
 while ( L1 >= A ) {
   g2 = g1;  
   L2 = L1;
   g1 /= 2.;
   lambda = ( 1. + zae ) / ( 1. - g1*(1.-a[N1]) - g1*nen );   
   L1 = g[i0] + g1 * lambda * gx[i0];   
   /*printf("g1 = %.6f,\tL1 = %.6f\n", g1, L1);*/
 }
 
 i = 0;
 do {
   i++;
   g3 = g1 + (A-L1)/(L2-L1) * (g2-g1);
   lambda = ( 1. + zae ) / ( 1. - g3*(1.-a[N1]) - g3*nen );
   L3 = g[i0] + g3 * lambda * gx[i0];
   /*printf("g3 = %.6f,\tL3 = %.6f\n", g3, L3);*/   
   g1 = g2; L1 = L2; g2 = g3; L2 = L3;
 } while ( fabs(g1-g2)>1e-9 && fabs(L3-A)>1e-9 && i < 100);
 
 *hm = ihm;
 *gamma = 1. - g3;
 
 Free(gx);
 Free(g);
 Free(rr);
 Free(xi);
 Free(psi);
 Free(phi);
 Free(z);
 Free(y);
 Free(x);
 Free(b3);
 Free(b2);
 Free(b1);
 Free(a);

 return 0;
}

/* HIER */
double ccusum_L_arl(double mu, int km, int hm, int m, int i0)
{ double *a, *b1, *b2, *x, *y, *z, *phi, *psi, *g, px, al, ga, et, de, be, arl;
  int i, j, l, lmax, N, N1;
 
 N = hm + 1; 
 N1 = N - 1;
 
 a   = vector(2*N-1);
 b1  = vector(N);
 b2  = vector(N);
 x   = vector(N);
 y   = vector(N);
 z   = vector(N);
 phi = vector(N);
 psi = vector(N);
 g   = vector(N);

 lmax = (int)ceil( (hm + km) / m ) + 1;
 
 for (l=0; l<=lmax; l++) {
   px = pdf_pois((double)l, mu);
   i = l*m - km;   
   if ( 0 <= N+i-1 && N+i-1 < 2*N-1 ) a[N+i-1] = -px;
   if ( 0 <= i-1 && i-1 < N ) b2[i-1] = px;
 }
 a[N1] += 1.;
 b2[N-1] = 1. - cdf_pois( (hm + km)/m, mu);
 
 for (i=N-1; i>=0; i--) {
   b1[i] = 1.;
   if ( i > 0 ) b2[i-1] += b2[i];
 }

 x[0] = 1./a[N1];
 y[0] = 1./a[N1];
 phi[0] = b1[0]/a[N1];
 psi[0] = b2[0]/a[N1];
 
 for (i=1; i<N; i++) {
   al = 0.;
   for (j=0; j<i; j++) al += a[N1 + i - j] * x[j];
   ga = 0.;
   for (j=0; j<i; j++) ga += a[N1 - 1 - j] * y[j];
   et = -b1[i];
   for (j=0; j<i; j++) et += a[N1 + i - j] * phi[j];
   de = -b2[i];
   for (j=0; j<i; j++) de += a[N1 + i - j] * psi[j];

   be = 1. - al*ga;
   
   z[0] = -ga*x[0] / be;
   for (j=1; j<i; j++) z[j] = ( y[j-1] - ga*x[j] ) / be;
   z[i] = y[i-1] / be;   

   x[0] = x[0] / be;
   for (j=1; j<i; j++) x[j] = ( x[j] - al*y[j-1] ) / be;
   x[i] = -al*y[i-1] / be;
   
   for (j=0; j<=i; j++) y[j] = z[j];
       
   for (j=0; j<i; j++) {
     phi[j] = ( phi[j] - et*z[j] );
     psi[j] = ( psi[j] - de*z[j] );
   }
   phi[i] = -et*z[i];
   psi[i] = -de*z[i];
 }

 be = phi[0] / ( 1. - psi[0] );
 for (i=0; i<N; i++) g[i] = phi[i] + psi[i] * be;
 arl = g[i0];

 Free(g);
 Free(psi);
 Free(phi);
 Free(z);
 Free(y);
 Free(x);
 Free(b2);
 Free(b1);
 Free(a);

 return arl;  
}

int ccusum_L_crit(double A, double mu0, int km, int m, int i0)
{ int p10, hm, p;
  double L1;
    
  p10 = (int)ceil( log10((double)m) );
  hm = 2 * (int)qf_pois(1.-1./A, mu0) * m;  /* borrow upper starting value (Shewhart) */
  L1 = ccusum_L_arl(mu0, km, hm, m, i0);  
    
  for (p=p10; p>=0; p--) {
    if ( L1 < A ) {
      while ( L1 < A ) {
        hm += pow(10., (double)p);
        L1 = ccusum_L_arl(mu0, km, hm, m, i0);        
      }
    } else {
      while ( L1 >= A ) {
        hm -= pow(10., (double)p);
        if ( hm < km ) {
          p--;
          hm += pow(10., (double)p+1.) - pow(10., (double)p);
        }    
        L1 = ccusum_L_arl(mu0, km, hm, m, i0);  
      }
    }
  }
  if ( L1 < A ) hm = hm + 1;      
  return hm;
}

double ccusum_L_arl_rando(double mu, int km, int hm, int m, double gamma, int i0)
{ double *a, *b1, *b2, *b3, *x, *y, *z, *phi, *psi, *xi, *rr, *g, *gx, px, al, ga, et, de, be, ka, lambda, arl, nen, zae;
  int i, j, l, lmax, N, N1;
  
 N  = hm; 
 N1 = N - 1;
 
 a   = vector(2*N-1);
 b1  = vector(N);
 b2  = vector(N);
 b3  = vector(N);
 x   = vector(N);
 y   = vector(N);
 z   = vector(N);
 phi = vector(N);
 psi = vector(N);
 xi  = vector(N);
 rr  = vector(N);
 g   = vector(N);
 gx  = vector(N);

 lmax = (int)ceil( (hm + km) / m ) + 1; 
 
 for (l=0; l<=lmax; l++) {
   px = pdf_pois((double)l, mu);
   i = l*m - km;
   if ( 0 <= N+i-1 && N+i-1 < 2*N-1 ) a[N+i-1] = -px;
   if ( 0 <= i-1 && i-1 < N ) {
    b2[i-1] = px;
    rr[N-i] = px;
   }
   if (0 <= N+i && N+i < N ) b3[N+i] = (1.-gamma)*px;
 }
 a[N1] += 1.;
 b2[N-1] = 1. - cdf_pois( (hm + km)/m, mu);
 
 for (i=N-1; i>=0; i--) {
   b1[i] = 1.;
   if ( i > 0 ) b2[i-1] += b2[i];
 }

 x[0] = 1./a[N1];
 y[0] = 1./a[N1];
 phi[0] = b1[0]/a[N1];
 psi[0] = b2[0]/a[N1];
 xi[0] = b3[0]/a[N1];
 
 for (i=1; i<N; i++) {
   al = 0.;
   for (j=0; j<i; j++) al += a[N1 + i - j] * x[j];
   ga = 0.;
   for (j=0; j<i; j++) ga += a[N1 - 1 - j] * y[j];
   et = -b1[i];
   for (j=0; j<i; j++) et += a[N1 + i - j] * phi[j];
   de = -b2[i];
   for (j=0; j<i; j++) de += a[N1 + i - j] * psi[j];
   ka = -b3[i];
   for (j=0; j<i; j++) ka += a[N1 + i - j] * xi[j];

   be = 1. - al*ga;
   
   z[0] = -ga*x[0] / be;
   for (j=1; j<i; j++) z[j] = ( y[j-1] - ga*x[j] ) / be;
   z[i] = y[i-1] / be;   

   x[0] = x[0] / be;
   for (j=1; j<i; j++) x[j] = ( x[j] - al*y[j-1] ) / be;
   x[i] = -al*y[i-1] / be;
   
   for (j=0; j<=i; j++) y[j] = z[j];
       
   for (j=0; j<i; j++) {
     phi[j] = ( phi[j] - et*z[j] );
     psi[j] = ( psi[j] - de*z[j] );
     xi[j]  = (  xi[j] - ka*z[j] );
   }
   phi[i] = -et*z[i];
   psi[i] = -de*z[i];
   xi[i]  = -ka*z[i];
 }

 be = phi[0] / ( 1. - psi[0] );
 for (i=0; i<N; i++) g[i] = phi[i] + psi[i] * be;
 
 ka = xi[0] / ( 1. - psi[0] );
 for (i=0; i<N; i++) gx[i] = xi[i] + psi[i] * ka;
 
 nen = 0.;
 zae = 0.;
 for (i=0; i<N; i++) {
   zae += rr[i] * g[i];
   nen += rr[i] * gx[i];
 }
 lambda = ( 1. + zae ) / ( 1. - (1.-gamma)*(1.-a[N1]) - nen );

 for (i=0; i<N; i++) g[i] += lambda*gx[i];

 arl = g[i0];
 
 Free(gx);
 Free(g);
 Free(rr);
 Free(xi);
 Free(psi);
 Free(phi);
 Free(z);
 Free(y);
 Free(x);
 Free(b3);
 Free(b2);
 Free(b1);
 Free(a);

 return arl;  
}

int ccusum_L_rando_crit(double A, double mu0, int km, int m, int i0, int *hm, double *gamma)
{ double *a, *b1, *b2, *b3, *x, *y, *z, *phi, *psi, *xi, *rr, *g, *gx, px, al, ga, et, de, be, ka, lambda, L1, L2, L3, nen, zae, g1, g2, g3;
  int i, j, l, lmax, N, N1, ihm;
  
 ihm = ccusum_L_crit(A, mu0, km, m, i0);
  
 N  = ihm; 
 N1 = N - 1;
 
 a   = vector(2*N-1);
 b1  = vector(N);
 b2  = vector(N);
 b3  = vector(N);
 x   = vector(N);
 y   = vector(N);
 z   = vector(N);
 phi = vector(N);
 psi = vector(N);
 xi  = vector(N);
 rr  = vector(N);
 g   = vector(N);
 gx  = vector(N);

 lmax = (int)ceil( (ihm + km) / m ) + 1;
 
 for (l=0; l<=lmax; l++) {
   px = pdf_pois((double)l, mu0);
   i = l*m - km;
   if ( 0 <= N+i-1 && N+i-1 < 2*N-1 ) a[N+i-1] = -px;
   if ( 0 <= i-1 && i-1 < N ) {
    b2[i-1] = px;
    rr[N-i] = px;
   }
   if (0 <= N+i && N+i < N ) b3[N+i] = px;
 }
 a[N1] += 1.;
 b2[N-1] = 1. - cdf_pois( (ihm + km)/m, mu0);
 
 for (i=N-1; i>=0; i--) {
   b1[i] = 1.;
   if ( i > 0 ) b2[i-1] += b2[i];
 }

 x[0] = 1./a[N1];
 y[0] = 1./a[N1];
 phi[0] = b1[0]/a[N1];
 psi[0] = b2[0]/a[N1];
 xi[0] = b3[0]/a[N1];
 
 for (i=1; i<N; i++) {
   al = 0.;
   for (j=0; j<i; j++) al += a[N1 + i - j] * x[j];
   ga = 0.;
   for (j=0; j<i; j++) ga += a[N1 - 1 - j] * y[j];
   et = -b1[i];
   for (j=0; j<i; j++) et += a[N1 + i - j] * phi[j];
   de = -b2[i];
   for (j=0; j<i; j++) de += a[N1 + i - j] * psi[j];
   ka = -b3[i];
   for (j=0; j<i; j++) ka += a[N1 + i - j] * xi[j];

   be = 1. - al*ga;
   
   z[0] = -ga*x[0] / be;
   for (j=1; j<i; j++) z[j] = ( y[j-1] - ga*x[j] ) / be;
   z[i] = y[i-1] / be;   

   x[0] = x[0] / be;
   for (j=1; j<i; j++) x[j] = ( x[j] - al*y[j-1] ) / be;
   x[i] = -al*y[i-1] / be;
   
   for (j=0; j<=i; j++) y[j] = z[j];
       
   for (j=0; j<i; j++) {
     phi[j] = ( phi[j] - et*z[j] );
     psi[j] = ( psi[j] - de*z[j] );
     xi[j]  = (  xi[j] - ka*z[j] );
   }
   phi[i] = -et*z[i];
   psi[i] = -de*z[i];
   xi[i]  = -ka*z[i];
 }

 be = phi[0] / ( 1. - psi[0] );
 for (i=0; i<N; i++) g[i] = phi[i] + psi[i] * be;
 
 ka = xi[0] / ( 1. - psi[0] );
 for (i=0; i<N; i++) gx[i] = xi[i] + psi[i] * ka;
 
 nen = 0.;
 zae = 0.;
 for (i=0; i<N; i++) {
  zae += rr[i] * g[i];
  nen += rr[i] * gx[i];
 } 

 /* secant rule etc. */
 
 g1 = 1.;
 lambda = ( 1. + zae ) / ( 1. - (1.-a[N1]) - nen );
 L1 = g[i0] + lambda*gx[i0]; 
 /*printf("g1 = %.6f,\tL1 = %.6f,\tA = %.0f\n", g1, L1, A);*/
 
 while ( L1 >= A ) {
   g2 = g1;  
   L2 = L1;
   g1 /= 2.;
   lambda = ( 1. + zae ) / ( 1. - g1*(1.-a[N1]) - g1*nen );   
   L1 = g[i0] + g1 * lambda * gx[i0];   
   /*printf("g1 = %.6f,\tL1 = %.6f\n", g1, L1);*/
 }
 
 i = 0;
 do {
   i++;
   g3 = g1 + (A-L1)/(L2-L1) * (g2-g1);
   lambda = ( 1. + zae ) / ( 1. - g3*(1.-a[N1]) - g3*nen );
   L3 = g[i0] + g3 * lambda * gx[i0];
   /*printf("g3 = %.6f,\tL3 = %.6f\n", g3, L3);*/   
   g1 = g2; L1 = L2; g2 = g3; L2 = L3;
 } while ( fabs(g1-g2)>1e-9 && fabs(L3-A)>1e-9 && i < 100);
 
 *hm = ihm;
 *gamma = 1. - g3;
 
 Free(gx);
 Free(g);
 Free(rr);
 Free(xi);
 Free(psi);
 Free(phi);
 Free(z);
 Free(y);
 Free(x);
 Free(b3);
 Free(b2);
 Free(b1);
 Free(a);

 return 0;
}

/* HIER */
double ccusum_2_arl(double mu, int km1, int hm1, int m1, int i01, int km2, int hm2, int m2, int i02)
{ double arl1, arl2, arl3, arl4, arl;

/* relation between 1- and 2-sided CUSUM schemes due to Lucas/Crosier 1982, Technometrics 24, 199-205;
   only for sufficiently small headstarts !!
*/
 
 arl1 = ccusum_U_arl(mu, km1, hm1, m1, 0);
 arl2 = ccusum_U_arl(mu, km1, hm1, m1, i01);
 arl3 = ccusum_L_arl(mu, km2, hm2, m2, 0);
 arl4 = ccusum_L_arl(mu, km2, hm2, m2, i02);

 arl = ( arl2*arl3 + arl1*arl4 - arl1*arl3 ) / ( arl1 + arl3 );
 return arl;
}

double ccusum_2_arl_rando(double mu, int km1, int hm1, int m1, double gamma1, int i01, int km2, int hm2, int m2, double gamma2, int i02)
{ double arl1, arl2, arl3, arl4, arl;

/* relation between 1- and 2-sided CUSUM schemes due to Lucas/Crosier 1982, Technometrics 24, 199-205;
   only for sufficiently small headstarts !!
*/

 arl1 = ccusum_U_arl_rando(mu, km1, hm1, m1, gamma1, 0);
 arl2 = ccusum_U_arl_rando(mu, km1, hm1, m1, gamma1, i01);
 arl3 = ccusum_L_arl_rando(mu, km2, hm2, m2, gamma2, 0);
 arl4 = ccusum_L_arl_rando(mu, km2, hm2, m2, gamma2, i02);
 
 arl = ( arl2*arl3 + arl1*arl4 - arl1*arl3 ) / ( arl1 + arl3 );
 
 return arl;
}



/* 2-sided tolerance limits factors */

/* Wald & Wolfowitz */

double r_Fww (int n, double r)
{ double x1, x2;
 x1 = 1./sqrt(n*1.) - r; x2 = x1 + 2.*r;
 return ( PHI(x2,0.) - PHI(x1,0.) );
}

double r_fww (int n, double r)
{ return(
   exp(-(1./n+r*r)/2.)*(exp(-r/sqrt(n*1.))+exp(r/sqrt(n*1.)))/sqrt(2.*PI)
  );
}

double rww(int n, double p)
{ double r;
 r = .5;
 do r = r - (r_Fww(n,r)-p)/r_fww(n,r);
 while ( fabs(r_Fww(n,r)-p) > 1e-8 );
 return r;
}

double kww(int n, double p, double a)
{ double k;
 k = rww(n,p);
 k *= sqrt( (n-1.) );
 k /= sqrt( qCHI(a,n-1) );
 return k;
}

/* exact by Gauss-Legendre quadrature */

double tl_rx_f(double x, double r)
{ return ( PHI(x+r,0.) - PHI(x-r,0.) );
}

double tl_rx(double x, double p)
{ double r1, r2, r3, f1, f2, f3;
 r1 = 1.; f1 = tl_rx_f(x,r1);
 r2 = .8; f2 = tl_rx_f(x,r2);
 do {
   r3 = r1 - (f1-p)*(r2-r1)/(f2-f1);
   f3 = tl_rx_f(x,r3);
   if (f3<p) { r1 = r3; f1 = f3; }
   else      { r2 = r3; f2 = f3; }
 } while ( (fabs(f3-p)>1e-8) && (fabs(r1-r2)>1e-8) ); 
 return r3;
}

double tl_niveau(int n, double p, double k, int m)
{ double ni, xmax, *w, *z, dn, rxi;
  int i;
 ni = 0.; 
 dn = (double) n;
 xmax = qPHI(1.-(1e-10)/2.)/sqrt(dn);
 w = vector(m);
 z = vector(m);
 gausslegendre(m,0.,xmax,z,w);
 for (i=0;i<m;i++) {
   rxi = tl_rx (z[i],p);
   ni += 2. * w[i] * (1-CHI((dn-1.)*rxi*rxi/k/k,n-1))
         * sqrt(dn)*phi(sqrt(dn)*z[i],0.);
 }
 Free(z);
 Free(w);
 return ni;
}

double tl_factor (int n, double p, double a, int m)
{ double k0, k1, k2, n0, n1, n2, dk;

 k1 = kww(n,p,a);
 k0 = k1 - .2; k1 += .2;
 n0 = tl_niveau(n,p,k0,m); 
 n1 = tl_niveau(n,p,k1,m);

 do {
   k2 = k0 + ( (1.-a) - n0 )/( n1 - n0 ) * ( k1 - k0);
   n2 = tl_niveau(n,p,k2,m);
/* Regula falsi */
   if ( n2 < (1.-a) ) { dk = k2-k0; k0 = k2; n0 = n2; }
   else               { dk = k1-k0; k1 = k2; n1 = n2; }
 } while ( ( fabs((1.-a)-n2) > 1e-8 ) && ( fabs(dk) > 1e-7 ) );
 return k2;
}


/* solution of Ax = b with nxn matrix A and and n-dim vectors x and b */
/* by means of LU decomposition etc.                                  */

int LU_decompose(double *a, int *ps, int n)
{ int i, j, k;
  int pii = 0;
  double pivot, biggest, mult, t, *lu, *scales;

  lu = matrix(n,n);
  scales = vector(n);

  for (i=0;i<n;i++) {
    biggest = 0.;
    for (j=0;j<n;j++)
      if (biggest < (t = fabs(lu[i*n+j] = a[i*n+j]))) biggest  = t;
    if (biggest != 0.) scales[i] = 1. / biggest;
    else {
      scales[i] = 0.;
      Free(lu); Free(scales);
      return(0);
    }
    ps[i] = i;
  }

  for (k=0;k<n-1;k++) {
    biggest = 0.;
    for (i=k;i<n;i++) {
      if (biggest < (t = fabs(lu[ps[i]*n+k]) * scales[ps[i]])) {
        biggest = t;
        pii = i;
      }
    }
    if (biggest == 0.) { Free(lu); Free(scales); return(0); }
    if (pii != k) {
      j = ps[k];
      ps[k] = ps[pii];
      ps[pii] = j;
    }
    pivot = lu[ps[k]*n+k];
    for (i=k+1;i<n;i++) {
      lu[ps[i]*n+k] = mult = lu[ps[i]*n+k] / pivot;
      if (mult != 0.) {
        for (j=k+1;j<n;j++)
          lu[ps[i]*n+j] -= mult * lu[ps[k]*n+j];
      }
    }
  }

  if (lu[ps[n-1]*n+n-1] == 0.) { Free(lu); Free(scales); return(0); }

  for (i=0;i<n;i++) for (j=0;j<n;j++) a[i*n+j] = lu[i*n+j];

  Free(lu); Free(scales);

  return(1);
}


void LU_solve(double *a, double *b, int n)
{ int i, j, *ps;
  double dot, *x;

  x = vector(n);
  ps = ivector(n);

  LU_decompose(a,ps,n);

  for (i=0;i<n;i++) {
    dot = 0.;
    for (j=0;j<i;j++)
      dot += a[ps[i]*n+j] * x[j];
    x[i] = b[ps[i]] - dot;
  }

  for (i=n-1;i>=0;i--) {
    dot = 0.;
    for (j=i+1;j<n;j++)
      dot += a[ps[i]*n+j] * x[j];
    x[i] = (x[i] - dot) / a[ps[i]*n+i];
  }

  for (i=0;i<n;i++) b[i] = x[i];

  Free(x); Free(ps);
}


void LU_solve2(double *a, double *b, int *ps, int n)
{ int i, j;
  double dot, *x;

  x = vector(n);

  for (i=0;i<n;i++) {
    dot = 0.;
    for (j=0;j<i;j++)
      dot += a[ps[i]*n+j] * x[j];
    x[i] = b[ps[i]] - dot;
  }

  for (i=n-1;i>=0;i--) {
    dot = 0.;
    for (j=i+1;j<n;j++)
      dot += a[ps[i]*n+j] * x[j];
    x[i] = (x[i] - dot) / a[ps[i]*n+i];
  }

  for (i=0;i<n;i++) b[i] = x[i];

  Free(x); 
}
