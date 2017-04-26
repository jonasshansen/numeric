#include<stdlib.h>
#include<stdio.h>
#include<assert.h>
//#include"qspline.h"

typedef struct {int n; double *x, *y, *b, *c;} qspline;

qspline *qspline_alloc(int n, double *x, double *y){
	qspline *S = (qspline*)malloc(sizeof(qspline)); //spline
	S-> b = (double*)malloc((n-1)*sizeof(double)); //b_i
	S-> c = (double*)malloc((n-1)*sizeof(double)); //c_i
	S-> x = (double*)malloc(n*sizeof(double)); // copy of x_i
	S-> y = (double*)malloc(n*sizeof(double)); // copy of y_i
	for(int i=0; i<n; i++){
		S-> x[i] = x[i];	
		S-> y[i] = y[i];
	}
	S-> n = n;

	double p[n-1], h[n-1];

	for(int i=0; i<n-1; i++){
		h[i] = x[i+1] - x[i]; assert(h[i]>0);
		p[i] = (y[i+1] - y[i])/h[i];
	}
	
	S->c[0]=0; //recursion up
	for(int i=0;i<n-2;i++) S->c[i+1]=(p[i+1]-p[i]-S->c[i]*h[i])/h[i+1];

	S->c[n-2]/=2; //recursion down
	for(int i=n-3;i>=0;i--) S->c[i]=(p[i+1]-p[i]-S->c[i+1]*h[i+1])/h[i];

	for(int i=0;i<n-1;i++) S->b[i]=p[i]-S->c[i]*h[i];
	return S;
}

double qspline_eval(qspline *S, double z){
	double* x=S->x;
	assert(z>=x[0] && z<=x[S->n-1]);
	int i=0, j=S->n-1;
	while(j-i>1){int m=(i+j)/2; if(z>x[m]) i=m; else j=m;}
	double h=z-x[i];
	return S->y[i]+h*(S->b[i]+h*S->c[i]);
}

double qspline_deriv(qspline *S, double z){
	double *x=S->x;
	if(z<x[0] || z>x[S->n-1]) printf("qseval: out of range: z=%g",z);
	int i=0, j=S->n-1;
	while(j-i>1){int m=(i+j)/2; if(z>x[m]) i=m; else j=m;}
	double h=z-x[i];
	return S->b[i]+h*S->c[i];
}

void qspline_free(qspline *S){
free(S->x);
free(S->y);
free(S->b);
free(S->c);
free(S);
}

