typedef struct {int n; double *x, *y, *b, *c;} qspline;
qspline* qspline_alloc(int n,double* x,double* y);
void qspline_free(qspline *S);
double qspline_eval(qspline *S, double z);
