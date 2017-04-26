#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include"qspline.h"
#define RND (double)rand()/RAND_MAX

double linterp(int n, double *x, double *y, double z);

int main(){
	int n = 10;
	double *x = (double*)calloc(n,sizeof(double));
	double *y = (double*)calloc(n,sizeof(double));
	for(int i=0; i<n; i++) x[i] = i;
	for(int i=0; i<n; i++) y[i] = RND;
	for(int i=0; i<n; i++) printf("%g\t%g\n",x[i],y[i]);
	printf("\n");

	int N=100;
	double *x_interp = (double*)calloc(n,sizeof(double));
	for(int i=0; i<N; i++) x_interp[i] = (double)i/(double)N*(double)n;	

	/*double *y_linterp = (double*)calloc(n,sizeof(double));
	for(int i=0; i<N; i++) y_linterp[i] = linterp(n,x,y,x_interp[i]);
	for(int i=0; i<N; i++) printf("%g\t%g\n",x_interp[i],y_linterp[i]);
	printf("\n");*/
	
	double *x2 = (double*)calloc(n,sizeof(double));
	for(int i=0; i<n; i++) x2[i] = i;
	double step = (x2[n-1]-x2[0])/(double)N;

	for(double z=x2[0]; z<=x2[n-1]; z+=step)
		printf("%g %g\n", z, linterp(n,x,y,z));
	printf("\n");

	qspline *S = qspline_alloc(n,x2,y);
	for(double z=x2[0]; z<=x2[n-1]; z+=step) 
		printf("%g %g\n", z, qspline_eval(S,z));	
	printf("\n");
	
	return 0;
}
