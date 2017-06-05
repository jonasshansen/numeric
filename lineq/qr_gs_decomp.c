#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "qr_gs_decomp.h"

void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R){
/* QR-decomposition of matrix A, A is replaced with Q, R is filled */
	int m = A->size2;
	for(int i=0;i<m;i++){
		gsl_vector_view e = gsl_matrix_column(A,i);
		double r = gsl_blas_dnrm2(&e.vector);
		gsl_matrix_set(R,i,i,r);
		gsl_vector_scale(&e.vector,1/r); //normalization
		for(int j=i+1;j<m;j++){
			gsl_vector_view q = gsl_matrix_column(A,j);
			double s=0; gsl_blas_ddot(&e.vector,&q.vector,&s);
			gsl_blas_daxpy(-s,&e.vector,&q.vector); //orthogonalization
			gsl_matrix_set(R,i,j,s);
		}
	}
}

void qr_gs_decomp_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x){
	int m=R->size1;
	gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,x); // Q^T*b -> x
	for(int i=m-1;i>=0;i--){                  // backsubstitution -> x
		double s=0;
		for(int k=i+1;k<m;k++) s+=gsl_matrix_get(R,i,k)*gsl_vector_get(x,k);
		gsl_vector_set(x,i,(gsl_vector_get(x,i)-s)/gsl_matrix_get(R,i,i));
	}
}

void qr_gs_decomp_inverse(gsl_matrix *Q, gsl_matrix* R, gsl_matrix *B) {
	int n=Q->size1;
	gsl_vector *b = gsl_vector_calloc(n);
	gsl_vector *x = gsl_vector_calloc(n);
	for(int i=0;i<n;i++){
		gsl_vector_set(b,i,1.0);
		qr_gs_decomp_solve(Q,R,b,x);
		gsl_matrix_set_col(B,i,x);
		gsl_vector_set(b,i,0.0);
	}
	gsl_vector_free(b);
	gsl_vector_free(x);
}
