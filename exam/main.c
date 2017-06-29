#include<math.h>
#include<stdio.h>
#include<assert.h> // For various checks during compilation

double adapt24(double f(double),double a, double b,
	double acc, double eps, double f2, double f3, int nrec, double scaling){

	assert(nrec<1000000); // Assert that the number of recursions doesn't explode

	double f1=f(a+(b-a)/6), f4=f(a+5*(b-a)/6);
	double Q=(2*f1+f2+f3+2*f4)/6*(b-a), q=(f1+f4+f2+f3)/4*(b-a); // Higher and lower order rules
	double tolerance=acc+eps*fabs(Q), error=fabs(Q-q)/2; // Estimate error

	if(error < tolerance) return Q;
	else{ // If the error is larger than the tolerance, subdivide the interval
		double Q1=adapt24(f,a,(a+b)/2,acc/scaling,eps,f1,f2,nrec+1,scaling);
		double Q2=adapt24(f,(a+b)/2,b,acc/scaling,eps,f3,f4,nrec+1,scaling);
		return Q1+Q2; 
	}
}

double adapt(double f(double),double a,double b,
	double acc,double eps,double scaling){

	double f2=f(a+2*(b-a)/6),f3=f(a+4*(b-a)/6);
	int nrec=0;
	return adapt24(f,a,b,acc,eps,f2,f3,nrec,scaling);
}

void printresults(double Q,double exact,int calls,double err_est,double err_act){
	printf("              Q = %g\n",Q);
	printf("          exact = %g\n",exact);
	printf("          calls = %d\n",calls);
	printf("estimated error = %g\n",err_est);
	printf("   actual error = %g\n",err_act);
}


int main(){  // Uses gcc nested functions

	double scaling = sqrt(2.0); // Rescaling of the absolute accuracy goal

	int calls=0; // Depth of nesting (number of calls)

	// Integral of 1/sqrt(x) from 0 to 1
	double a=0,b=1,acc=0.001,eps=0.001;
	double f(double x){calls++; return 1/sqrt(x);}; // Nested function
	calls=0;
	scaling=sqrt(2.0);
	double Q=adapt(f,a,b,acc,eps,scaling);
	double exact=2;
	double err_est=acc+fabs(Q)*eps;
	double err_act=fabs(Q-exact);
	printf("Integrating 1/sqrt(x) from %g to %g with rescaling = %g\n",a,b,scaling);
	printresults(Q,exact,calls,err_est,err_act);

	calls=0;
	scaling=2;
	Q=adapt(f,a,b,acc,eps,scaling);
	err_est=acc+fabs(Q)*eps;
	err_act=fabs(Q-exact);
	printf("Integrating 1/sqrt(x) from %g to %g with rescaling = %g\n",a,b,scaling);	
	printresults(Q,exact,calls,err_est,err_act);


	// Integral of log(x)/sqrt(x) from 0 to 1
	a=0,b=1,acc=0.001,eps=0.001;
	double f2(double x){calls++; return log(x)/sqrt(x);}; // Nested function
	calls=0;
	scaling=sqrt(2.0);
	Q=adapt(f2,a,b,acc,eps,scaling);
	exact=-4;
	err_est=acc+fabs(Q)*eps;
	err_act=fabs(Q-exact);
	printf("\nIntegrating log(x)/sqrt(x) from %g to %g with rescaling = %g\n",a,b,scaling);
	printresults(Q,exact,calls,err_est,err_act);

	calls=0;
	scaling=2;
	Q=adapt(f2,a,b,acc,eps,scaling);
	err_est=acc+fabs(Q)*eps;
	err_act=fabs(Q-exact);
	printf("Integrating log(x)/sqrt(x) from %g to %g with rescaling = %g\n",a,b,scaling);	
	printresults(Q,exact,calls,err_est,err_act);


	// Integral of sqrt(x) from 0 to 1
	a=0,b=1,acc=0.001,eps=0.001;
	double f3(double x){calls++; return sqrt(x);}; // Nested function
	calls=0;
	scaling=sqrt(2.0);
	Q=adapt(f3,a,b,acc,eps,scaling);
	exact=0.66666666;
	err_est=acc+fabs(Q)*eps;
	err_act=fabs(Q-exact);
	printf("\nIntegrating sqrt(x) from %g to %g with rescaling = %g\n",a,b,scaling);
	printresults(Q,exact,calls,err_est,err_act);

	calls=0;
	scaling=2;
	Q=adapt(f3,a,b,acc,eps,scaling);
	err_est=acc+fabs(Q)*eps;
	err_act=fabs(Q-exact);
	printf("Integrating sqrt(x) from %g to %g with rescaling = %g\n",a,b,scaling);	
	printresults(Q,exact,calls,err_est,err_act);


	// Integral of 4*sqrt(1-(1-x)^2) from 0 to 1
	a=0,b=1,acc=0.001,eps=0.001;
	double f4(double x){calls++; return 4*sqrt(1-(1-x)*(1-x));}; // Nested function
	calls=0;
	scaling=sqrt(2.0);
	Q=adapt(f4,a,b,acc,eps,scaling);
	exact=3.14159265;
	err_est=acc+fabs(Q)*eps;
	err_act=fabs(Q-exact);
	printf("\nIntegrating 4*sqrt(1-(1-x)^2) from %g to %g with rescaling = %g\n",a,b,scaling);
	printresults(Q,exact,calls,err_est,err_act);

	calls=0;
	scaling=2;
	Q=adapt(f4,a,b,acc,eps,scaling);
	err_est=acc+fabs(Q)*eps;
	err_act=fabs(Q-exact);
	printf("Integrating 4*sqrt(1-(1-x)^2) from %g to %g with rescaling = %g\n",a,b,scaling);	
	printresults(Q,exact,calls,err_est,err_act);

return 0;
}
