#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define ewma1 0
#define ewma2 1
#define fix 0
#define vacl 1
#define fir 2
#define both 3
#define steiner 4
#define stat 5
#define fink 6

#define conditional 0
#define cyclical 1

extern double rho0;

double xe1_iglad(double l, double c, double zr, double mu0, double mu1, int N);
double xe1_arlm(double l, double c, double zr, double hs, int q, double mu0, double mu1, int mode, int N, int nmax);

double xe2_iglad (double l, double c, double mu0, double mu1, int N);
double xe2_igladc(double l, double c, double mu0, double mu1, double z0, int N);
double xe2_arlm  (double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax);
double xe2_arlmc(double l, double c, double hs, int q, double mu0, double mu1, int mode, int N, int nmax);

void xewma_ad(int *ctyp, double *l, double *c, double *zr, double *mu0, double *mu1, double *z0, int *ltyp, int *styp, int *r, double *ad)
{ int nmax=1000000, q2=200;
  double ad1=1., ad2=1., eps=.0001;
 if ( *styp==conditional ) {
   if ( *ctyp==ewma1 && *ltyp==fix ) *ad = xe1_iglad(*l,*c,*zr,*mu0,*mu1,*r);
   
   if ( *ctyp==ewma1 && *ltyp>fix )  { 
     ad1 = xe1_arlm(*l,*c,*zr,0.,200,*mu0,*mu1,*ltyp,*r,nmax);
     ad2 = xe1_arlm(*l,*c,*zr,0.,300,*mu0,*mu1,*ltyp,*r,nmax);
     if ( fabs(ad1-ad2) > eps ) {
       q2 = 300;
       while ( fabs(ad1-ad2) > eps ) {         
         ad1 = ad2;
         q2 += 100;
         ad2 = xe1_arlm(*l,*c,*zr,0.,q2,*mu0,*mu1,*ltyp,*r,nmax);
       }
     }
     *ad = ad2;
   }
   if ( *ctyp==ewma2 && *ltyp==fix ) *ad = xe2_iglad(*l,*c,*mu0,*mu1,*r);
   
   if ( *ctyp==ewma2 && *ltyp>fix )  {
     ad1 = xe2_arlm(*l,*c,0.,200,*mu0,*mu1,*ltyp,*r,nmax); 
     /*printf("q = 200,\tad = %.5f\n", ad1);*/
     ad2 = xe2_arlm(*l,*c,0.,300,*mu0,*mu1,*ltyp,*r,nmax);
     /*printf("q = 300,\tad = %.5f\n", ad2);*/
     if ( fabs(ad1-ad2) > eps ) {
       /*printf("\nTuning needed\n\n");*/
       q2 = 300;
       while ( fabs(ad1-ad2) > eps ) {         
         ad1 = ad2;
         q2 += 100;
         ad2 = xe2_arlm(*l,*c,0.,q2,*mu0,*mu1,*ltyp,*r,nmax);
         /*printf("q = %d,\tad = %.5f\n", q2, ad2);*/
       }
       /*printf("\n");*/
     }
     *ad = ad2;
   }
 } else {
   if ( *ctyp==ewma2 && *ltyp==fix ) *ad = xe2_igladc(*l, *c, *mu0, *mu1, *z0, *r);
   
   if ( *ctyp==ewma2 && *ltyp>fix )  {
     ad1 = xe2_arlmc(*l,*c,0.,200,*mu0,*mu1,*ltyp,*r,nmax);
     ad2 = xe2_arlmc(*l,*c,0.,300,*mu0,*mu1,*ltyp,*r,nmax);
     if ( fabs(ad1-ad2) > eps ) {
       q2 = 300;
       while ( fabs(ad1-ad2) > eps ) {         
         ad1 = ad2;
         q2 += 100;
         ad2 = xe2_arlmc(*l,*c,0.,q2,*mu0,*mu1,*ltyp,*r,nmax);
       }
     }
     *ad = ad2;
   }
 }
}
