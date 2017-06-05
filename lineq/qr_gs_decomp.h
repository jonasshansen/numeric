#ifndef HAVE_GRAMSCHMIDT_H
void qr_gs_decomp(gsl_matrix*A, gsl_matrix*R);
void qr_gs_decomp_solve(gsl_matrix*Q,gsl_matrix*R,gsl_vector*b,gsl_vector*x);
void qr_gs_decomp_inverse(gsl_matrix *A, gsl_matrix *R, gsl_matrix* B);
#define HAVE_GRAMSCHMIDT_H
#endif
