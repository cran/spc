#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>

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
#define sven 5
#define fink 6

/* export */

double xc_crit(int ctyp, double k, double L0, double hs, double m0, int N);

double xc1_iglarl(double k, double h, double hs, double mu, int N); 
double xc1_iglad (double k, double h, double mu0, double mu1, int N);

double xc2_iglarl(double k, double h, double hs, double mu, int N);
double xc2_iglad (double k, double h, double mu0, double mu1, int N);

double xcC_iglarl(double k, double h, double hs, double mu, int N);
double xcC_iglad (double k, double h, double mu0, double mu1, int N);

double xe_crit(int ctyp, double l, double L0, double zr, double hs,
               double m0, int ltyp, int N);

double xe1_iglarl(double l, double c, double zr, double hs, double mu, int N);
double xe1_iglad (double l, double c, double zr, double mu0, double mu1, int N);

double xe2_iglarl(double l, double c, double hs, double mu, int N);
double xe2_iglad (double l, double c, double mu0, double mu1, int N);
double xe2_arlm(double l, double c, double hs, int q, double mu0, double mu1, 
                int mode, int N, int nmax);

/* internal functions etc. */

static void gausslegendre(int n, double x1, double x2, double *x, double *w);

void LU_solve(double *a, double *b, int n);

void pmethod(int n, double *p, int *status, double *lambda,
             double x_[], int *noofit);

   int *ivector(long n);
double *vector (long n);
double *matrix(long m, long n);

double phi(double x, double mu);
double PHI(double x, double mu);

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

/* normal density */

double phi(double x, double mu)
{
 return dnorm(x,mu,1.,LOG);
}

/* normal cumulative distribution function */

double PHI(double x, double mu)
{
 return pnorm(x,mu,1.,TAIL,LOG);
}

/* roots and abscissae of Gauﬂ-Legendre quadrature */

#define GLeps 3e-11

void gausslegendre(int n, double x1, double x2, double *x, double *w)
/*
   The following algorithm is based on ideas of Knut Petras
   (see http://www-public.tu-bs.de:8080/~petras/).

   The nodes are derived by means of the Newton method.
   Afterwards, the weights are obtained by exploiting
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
{ double xw, xmid, z0, z1, diff, p0, p1, p2, a;
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
#define maxits          5000

void pmethod(int n, double *p, int *status, double *lambda, 
             double x_[], int *noofit)
{ int count, i, newi, oldi, D1;
  double newmu, oldmu, lastmu, *z, *y_;
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
 else { lastmu = newmu;   *noofit = maxits; }
}

/* ************************************************************************* */
/*      zero-state and steady-state ARl and critical value routines          */

double xc_crit(int ctyp, double k, double L0, double hs, double m0, int N)
{ double c1, c2, c3, L1, L2, L3, dc;
  int m = 1;

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

 do {
   c3 = c1 + (L0-L1)/(L2-L1) * (c2-c1);
   if (ctyp==cusum1) L3 = xc1_iglarl ( k,c3,hs,m0,N );
   if (ctyp==cusum2) L3 = xc2_iglarl ( k,c3,hs,m0,N );
   if (ctyp==cusumC) L3 = xcC_iglarl ( k,c3,hs,m0,N );
/* Regula falsi */
/*   if (L3<L0) { dc=c3-c1; c1 = c3; L1 = L3; }
   else       { dc=c2-c1; c2 = c3; L2 = L3; }*/
/* Sekantenverfahren */
   dc = c3-c2; c1 = c2; L1 = L2; c2 = c3; L2 = L3;
 } while ( (fabs(L0-L3)>1e-5) && (fabs(dc)>1e-6) );
 return c3;
}

double xe_crit(int ctyp, double l, double L0, double zr, double hs,
               double m0, int ltyp, int N)
{ double c1, c2, c3, L1, L2, L3, dc;
  int m = 1;

 c2 = 0.;
 do {
   c2 += .5;
   if (ctyp==ewma1) L2 = xe1_iglarl ( l,c2,zr,hs,m0,N );
   if (ctyp==ewma2 && ltyp==fix) L2 = xe2_iglarl ( l,c2,hs,m0,N );
   if (ctyp==ewma2 && ltyp>fix) {
     if (hs<0. && ltyp==fir) L2 = xe2_arlm ( l,c2,c2/2.,1,m0,m0,ltyp,N,5000 );
     if (hs<0. && ltyp==both) 
       L2 = xe2_arlm ( l,c2,c2/2.*sqrt(l*(2.-l)),1,m0,m0,ltyp,N,5000 ); 
     if (hs>=0.) L2 = xe2_arlm ( l,c2,hs,1,m0,m0,ltyp,N,5000 );
   }
 } while (L2<L0);

 c1 = c2 - .5;
 if (ctyp==ewma1) L1 = xe1_iglarl ( l,c1,zr,hs,m0,N );
 if (ctyp==ewma2 && ltyp==fix) L1 = xe2_iglarl ( l,c1,hs,m0,N );
 if (ctyp==ewma2 && ltyp>fix) {
   if (hs<0. && ltyp==fir) L1 = xe2_arlm ( l,c1,c1/2.,1,m0,m0,ltyp,N,5000 );
   if (hs<0. && ltyp==both)
     L1 = xe2_arlm ( l,c1,c1/2.*sqrt(l*(2.-l)),1,m0,m0,ltyp,N,5000 ); 
   if (hs>=0.) L1 = xe2_arlm ( l,c1,hs,1,m0,m0,ltyp,N,5000 );
 }

 do {
   c3 = c1 + (L0-L1)/(L2-L1) * (c2-c1);
   if (ctyp==ewma1) L3 = xe1_iglarl ( l,c3,zr,hs,m0,N );
   if (ctyp==ewma2 && ltyp==fix) L3 = xe2_iglarl ( l,c3,hs,m0,N );
   if (ctyp==ewma2 && ltyp>fix) {
     if (hs<0. && ltyp==fir) L3 = xe2_arlm ( l,c3,c3/2.,1,m0,m0,ltyp,N,5000 );
     if (hs<0. && ltyp==both)
       L3 = xe2_arlm ( l,c3,c3/2.*sqrt(l*(2.-l)),1,m0,m0,ltyp,N,5000 );
     if (hs>=0.) L3 = xe2_arlm ( l,c3,hs,1,m0,m0,ltyp,N,5000 );
   }
/* Regula falsi */
/*   if (L3<L0) { dc=c3-c1; c1 = c3; L1 = L3; }
   else       { dc=c2-c1; c2 = c3; L2 = L3; }*/
/* Sekantenverfahren */
   dc = c3-c2; c1 = c2; L1 = L2; c2 = c3; L2 = L3;
 } while ( (fabs(L0-L3)>1e-5) && (fabs(dc)>1e-6) );
 return c3;
}

double xc1_iglarl (double k, double h, double hs, double mu, int N)
{ double *a, d, *g, *w, *z, arl;
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

double xcC_iglarl (double k, double h, double hs, double mu, int N)
{ double *a, d, *g, *w, *z, arl;
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

double xe1_iglarl(double l, double c, double zr, double hs, double mu, int N)
{ double *a, d, *g, *w, *z, h, arl;
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
{ double *a, d, *g, *w, *z, h, arl;
  int i, j;

 a = matrix(N,N);
 g = vector(N);
 w = vector(N);
 z = vector(N);

 c  *= sqrt( l/(2.-l) ); 
 hs *= sqrt( l/(2.-l) );

 gausslegendre(N,-c,c,z,w);

 for (i=0;i<N;i++) {
   for (j=0;j<N;j++) a[i*N+j] = -w[j]/l * phi((z[j]-(1.-l)*z[i])/l,mu);
   ++a[i*N+i];
 }

 for (j=0;j<N;j++) g[j] = 1.;
 LU_solve(a,g,N);

 arl = 1.;
 for (j=0;j<N;j++)
   arl += w[j]/l * phi((z[j]-(1.-l)*hs)/l,mu) * g[j];

 Free(a);
 Free(g);
 Free(w);
 Free(z);

 return arl;
}

double xc1_iglad (double k, double h, double mu0, double mu1, int N)
{ double *a, d, *w, *z, *arl, *psi, rho, ad, norm;
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
{ double *a, d, *w, *z, *arl, *psi, rho, ad, norm;
  int *indx, i, j, status, noofit, NN, job=0;

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
{ double *a, d, *w, *z, h, *arl, *psi, rho, ad, norm;
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
{ double *a, d, *w, *z, h, *arl, *psi, rho, ad, norm;
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

double xe2_arlm(double l, double c, double hs, int q, double mu0, double mu1,
                int mode, int N, int nmax)
{ double *Smatrix, *p0, *fn, *w, *z,
         arl0, var0, rho, dn, rn, cn, rn0, cn0, delta=0.,
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

 /* in-control, i. e. n<=m-1 */

 for (n=1;n<=q-1;n++) {
  nn = (double) n;

  /* determine c_n and r_n, n=1,2,...,q-1 */
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
 
  /* determine f_n, n=1,2,...,q-1 */
  if (n==1) {
    for (i=0;i<N;i++)
      if (mode==sven)
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

 arl0 = 1.; var0 = 0.; rho = 0.;

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
      if (mode==sven)
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
      if (n==q && q>1) fn[(n-1)*N+i] /= p0[q-1];
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

double xc2_iglad (double k, double h, double mu0, double mu1, int N)
{ double *a, d, *arl, *psi, rho, ad, norm, 
         z1, z2, z11, z12, z21, z22, w;
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

/* Richardson extrapolation */
double xc2_igladR (double k, double h, double mu0, double mu1, int r)
{ double *a, *b, d, ad;
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
