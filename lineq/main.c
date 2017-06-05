#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include"qr_gs_decomp.h"
#define RND ((double)rand()/RAND_MAX)
#define FMT "%7.3f"

void printm(gsl_matrix *A){
	for(int r=0;r<A->size1;r++){
		for(int c=0;c<A->size2;c++) printf(FMT,gsl_matrix_get(A,r,c));
		printf("\n");}
}

void printv(gsl_vector *b){
		for(int i=0;i<b->size;i++){
      printf(FMT,gsl_vector_get(b,i));
      printf("\n");
    }
}

int main(int argc, char** argv){

printf("\n---PART A1---\n");
size_t n=4,m=3; // size of random tall matrix A
gsl_matrix *A = gsl_matrix_calloc(n,m);
gsl_matrix *R = gsl_matrix_calloc(m,m);
for(int i=0;i<n;i++)for(int j=0;j<m;j++)gsl_matrix_set(A,i,j,RND);

printf("Random tall matrix A:\n");
printm(A);

qr_gs_decomp(A,R);

printf("Q:\n");
printm(A);
printf("R (upper triangular):\n");
printm(R);
printf("Q^T Q = 1:\n");
gsl_matrix *qtq = gsl_matrix_calloc(m,m);
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,A,A,0.0,qtq);
printm(qtq);
printf("QR (same as A):\n");
gsl_matrix *qr = gsl_matrix_calloc(n,m);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,R,0.0,qr);
printm(qr);


printf("\n---PART A2---\n");
gsl_matrix *A2 = gsl_matrix_calloc(m,m);
for(int i=0;i<m;i++)for(int j=0;j<m;j++)gsl_matrix_set(A2,i,j,RND);
gsl_matrix *R2 = gsl_matrix_calloc(m,m);
gsl_vector *b = gsl_vector_calloc(m);
for(int i=0;i<m;i++) gsl_vector_set(b,i,RND);
gsl_matrix *A2copy = gsl_matrix_calloc(m,m);
gsl_matrix_memcpy(A2copy,A2);

printf("Random square matrix A2:\n");
printm(A2);
printf("Random vector b:\n");
printv(b);

qr_gs_decomp(A2,R2);

gsl_vector *x = gsl_vector_calloc(m);
qr_gs_decomp_solve(A2,R2,b,x);

printf("Solution x:\n");
printv(x);
printf("Ax (equal to b):\n");
gsl_vector *Ax = gsl_vector_calloc(m);
gsl_blas_dgemv(CblasNoTrans,1.0,A2copy,x,0.0,Ax);
printv(Ax);


printf("\n---PART B---\n");
gsl_matrix *A3 = gsl_matrix_calloc(m,m);
for(int i=0;i<m;i++)for(int j=0;j<m;j++)gsl_matrix_set(A3,i,j,RND);
gsl_matrix *A3copy = gsl_matrix_calloc(m,m);
gsl_matrix_memcpy(A3copy,A3);
gsl_matrix *R3 = gsl_matrix_calloc(m,m);
gsl_matrix *B3 = gsl_matrix_calloc(m,m);

qr_gs_decomp(A3,R3);
qr_gs_decomp_inverse(A3,R3,B3);

printf("A:\n");
printm(A3copy);
printf("B (inverse of A):\n");
printm(A3);
printf("AB:\n");
gsl_matrix *AB = gsl_matrix_calloc(m,m);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A3copy,B3,0.0,AB);
printm(AB);



gsl_matrix_free(A);
gsl_matrix_free(R);
gsl_matrix_free(qtq);
gsl_matrix_free(qr);
gsl_matrix_free(A2);
gsl_matrix_free(R2);
gsl_vector_free(b);
gsl_matrix_free(A2copy);
gsl_vector_free(x);
gsl_vector_free(Ax);
gsl_matrix_free(A3);
gsl_matrix_free(R3);
gsl_matrix_free(B3);
gsl_matrix_free(A3copy);
gsl_matrix_free(AB);
}
