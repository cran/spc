#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>

#define GL 0
#define Ra 1

double *vector (long n);
void gausslegendre(int n, double x1, double x2, double *x, double *w);
void radau(int n, double x1, double x2, double *x, double *w);

void quadrature_nodes_weights(int *n, double *x1, double *x2, int *type, double *nodes_weights)
{ double *knoten, *gewichte;
  int i; 
  
 knoten   = vector(*n); 
 gewichte = vector(*n);
 
 if ( *type==GL ) gausslegendre(*n, *x1, *x2, knoten, gewichte);
 if ( *type==Ra ) radau(*n, *x1, *x2, knoten, gewichte);
 
 for (i=0; i<*n; i++) {
   nodes_weights[i]    = knoten[i];
   nodes_weights[i+*n] = gewichte[i];
 }
 
 Free(gewichte);
 Free(knoten);
}
