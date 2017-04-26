#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"qspline.h"
#define RND (double)rand()/RAND_MAX

/*
struct qspline qspline_alloc(int n, double *x, double *y){
//allocate memory, calculate spline coefficients b and c	
}

double qspline_evaluate(struct qsline *s, double z){
//compute the spline interpolation function at a given point z
}

void qspline_free(struct qspline *s){
//free the memory occupied by structure qspline
}
*/

double linterp(int n, double *x, double *y, double z);


int main(){
	int n = 10;
	double *x = (double*)calloc(n,sizeof(double));
	double *y = (double*)calloc(n,sizeof(double));
	for(int i=0; i<n; i++) x[i] = i;
	for(int i=0; i<n; i++) y[i] = RND;
	printf("x\ty\n");
	for(int i=0; i<n; i++) printf("%g\t%g\n",x[i],y[i]);

	int N=40;
	double *x_interp = (double*)calloc(n,sizeof(double));
	for(int i=0; i<N; i++) x_interp[i] = (double)i/(double)N*(double)n;	

	double *y_linterp = (double*)calloc(n,sizeof(double));
	for(int i=0; i<N; i++) y_linterp[i] = linterp(n,x,y,x_interp[i]);
	printf("x_interp\ty_linterp\n");
	for(int i=0; i<N; i++) printf("%g\t%g\n",x_interp[i],y_linterp[i]);
	
	double *x2 = (double*)calloc(n,sizeof(double));
	for(int i=0; i<n; i++) x2[i] = i;

	double step = (x2[n-1]-x2[0])/100;
	qspline *S = qspline_alloc(n,x2,y);
	for(double z=x2[0]; z<=x2[n-1]; z+=step) 
		printf("%g %g\n", z, qspline_eval(S,z));	
	
	qspline_free(S);
	

	return 0;
}
