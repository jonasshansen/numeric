#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<gramschmidt_qr.h>
void lsfit(
	int m, double f(int i,double x),
	gsl_vector* x, gsl_vector* y, gsl_vector* dy,
	gsl_vector* c, gsl_matrix* S)
{
int n = x->size;

gsl_matrix *A    = gsl_matrix_alloc(n,m);
gsl_vector *b    = gsl_vector_alloc(n);
gsl_matrix *R    = gsl_matrix_alloc(m,m);
gsl_matrix* invR = gsl_matrix_alloc(m,m);
gsl_matrix *I    = gsl_matrix_alloc(m,m);

for(int i=0;i<n;i++){
	double xi  = gsl_vector_get(x ,i);
	double yi  = gsl_vector_get(y ,i);
	double dyi = gsl_vector_get(dy,i);
	gsl_vector_set(b,i,yi/dyi);
	for(int k=0;k<m;k++)gsl_matrix_set(A,i,k,f(k,xi)/dyi);
	}
gramschmidt_qr(A,R);
gramschmidt_qr_solve(A,R,b,c);

gsl_matrix_set_identity(I);
gramschmidt_qr_inverse(I,R,invR);
gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,invR,invR,0,S);

gsl_vector_free(b);
gsl_matrix_free(invR);
gsl_matrix_free(R);
gsl_matrix_free(A);
gsl_matrix_free(I);
}
