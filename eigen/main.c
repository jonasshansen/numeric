#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#define RND ((double)rand()/RAND_MAX)
#define FMT "%7.3f" // print format
#define max_print 11 // only print if matrix dim are less than this value

int jacobi(gsl_matrix*A, gsl_vector*e, gsl_matrix*V);

void matrix_print(gsl_matrix *A){
	for(int r=0;r<A->size1;r++){
		for(int c=0;c<A->size2;c++) fprintf(stderr,FMT,gsl_matrix_get(A,r,c));
		fprintf(stderr,"\n");
	}
}

void vector_print(gsl_vector *v){
	for(int i=0;i<v->size;i++) fprintf(stderr,FMT,gsl_vector_get(v,i));
	fprintf(stderr,"\n");
}


int main(int argc, char** argv){

int n=(argc>1? atoi(argv[1]):5); // convert input to string

gsl_matrix *A = gsl_matrix_alloc(n,n);
gsl_matrix *B = gsl_matrix_alloc(n,n);
for(int i=0;i<n;i++) for(int j=i;j<n;j++) {
	double x = RND;
	gsl_matrix_set(A,i,j,x);
	gsl_matrix_set(A,j,i,x);
}
gsl_matrix_memcpy(B,A);


gsl_matrix *V = gsl_matrix_alloc(n,n);
gsl_vector *e = gsl_vector_alloc(n);
int sweeps = jacobi(A,e,V);
printf("n=%i, sweeps=%i\n",n,sweeps);

if(n<max_print){
	fprintf(stderr,"\nRandom symmetric matrix A: \n"); matrix_print(B);
	fprintf(stderr,"\nResult of Jacobi diagonalization: \n");
	fprintf(stderr,"Number of sweeps: %d\n",sweeps);
	fprintf(stderr,"Eigenvalues:\n"); vector_print(e);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,B,V,0,A);
	gsl_blas_dgemm(CblasTrans  ,CblasNoTrans,1,V,A,0,B);
	fprintf(stderr, "V^T*A*V should be a diagonal matrix with the above eigenvalues on the diagonal:\n");
	matrix_print(B);
}

return 0;
}
