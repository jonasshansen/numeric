#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<stdio.h>

int newton(
	void f(gsl_vector* x,gsl_vector* fx),
	gsl_vector* x, double dx, double eps);

void vector_print(char* s,gsl_vector* v){
	printf("%s",s);
	for(int i=0;i<v->size;i++)printf("%10.4g ",gsl_vector_get(v,i));
	printf("\n");
	}


int main() {

// Rosenbrock
int ncalls=0;

void f(gsl_vector* p,gsl_vector* fx){
	ncalls++;
	double x0=gsl_vector_get(p,0), y0=gsl_vector_get(p,1);
	gsl_vector_set(fx,0, 2*(1-x0)*(-1)+100*2*(y0-x0*x0)*(-1)*2*x0);
	gsl_vector_set(fx,1, 100*2*(y0-x0*x0));
}

gsl_vector* x0=gsl_vector_alloc(2);
gsl_vector_set(x0,0,-2);
gsl_vector_set(x0,1,8);
gsl_vector* fx=gsl_vector_alloc(2);

printf("Root finding:\n");
printf("Extremum of the Rosenbrock's function:\n");
vector_print("Initial vector x: ",x0);
f(x0,fx);
vector_print("            f(x): ",fx);
newton(f,x0,1e-6,1e-3);

printf("ncalls = %i\n",ncalls);
gsl_vector_fprintf(stderr,x0,"%g");
vector_print("      solution x: ",x0);
f(x0,fx);
vector_print("            f(x): ",fx);


// Himmelblau
ncalls=0;

void g(gsl_vector* p,gsl_vector* gx){
	ncalls++;
	double x1=gsl_vector_get(p,0), y1=gsl_vector_get(p,1);
	gsl_vector_set(gx,1,2*(2*x1*(x1*x1+y1-11)+x1+y1*y1-7));
	gsl_vector_set(gx,1,2*(x1*x1+2*y1*(x1+y1*y1-7)+y1-11));
}

gsl_vector* x1=gsl_vector_alloc(2);
gsl_vector_set(x1,0,2.9);
gsl_vector_set(x1,1,2-1);
gsl_vector* gx=gsl_vector_alloc(2);

printf("Root finding:\n");
printf("Extremum of the Himmelblau's function:\n");
vector_print("Initial vector x: ",x1);
g(x1,gx);
vector_print("            f(x): ",gx);

/* // the following line takes too long to run
newton(g,x1,1e-4,1e-3);

printf("ncalls = %i\n",ncalls);
gsl_vector_fprintf(stderr,x1,"%g");
vector_print("      solution x: ",x1);
g(x1,gx);
vector_print("            g(x): ",gx);
*/

return 0;
}
