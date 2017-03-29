#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#define RND (double)rand()/RAND_MAX

struct qspline qspline_alloc(int n, double *x, double *y){
	
}

int main(){
	int n = 5;
	gsl_vector * x = gsl_vector_alloc(n);
	for(int i=0; i<x->size; i++) gsl_vector_set(x,i,RND);
	for(int i=0; i<x->size; i++) printf("%g\n", gsl_vector_get(x,i));
	gsl_vector_free(x);

	struct qspline {int n; double *x, *y, *b, *c};

	return 0;
}
