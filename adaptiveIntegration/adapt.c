#include<math.h>
#include<assert.h> // For various checks during compilation
#include<stdio.h>

double adapt24(double f(double),double a, double b,
	double acc, double eps, double f2, double f3, int nrec){

	assert(nrec<1000000);
	double f1=f(a+(b-a)/6), f4=f(a+5*(b-a)/6);
	double Q=(2*f1+f2+f3+2*f4)/6*(b-a), q=(f1+f4+f2+f3)/4*(b-a);
	double tolerance=acc+eps*fabs(Q), error=fabs(Q-q)/2;
	if(error < tolerance) return Q;
	else{
		double Q1=adapt24(f,a,(a+b)/2,acc/sqrt(2.),eps,f1,f2,nrec+1);
		double Q2=adapt24(f,(a+b)/2,b,acc/sqrt(2.),eps,f3,f4,nrec+1);
		return Q1+Q2; 
	}
}

double adapt(double f(double),double a,double b,
	double acc,double eps){

	double f2=f(a+2*(b-a)/6),f3=f(a+4*(b-a)/6);
	int nrec=0;
	return adapt24(f,a,b,acc,eps,f2,f3,nrec);
}


int main(){  //uses gcc nested functions

	int calls=0;
	double a=0,b=1,acc=0.001,eps=0.001;
	double f(double x){calls++; return 1/sqrt(x);}; //nested function
	calls=0;
	double Q=adapt(f,a,b,acc,eps);
	double exact=2;
	printf("Integrating 1/sqrt(x) from %g to %g\n",a,b);
	printf("              Q = %g\n",Q);
	printf("          exact = %g\n",exact);
	printf("          calls = %d\n",calls);
	printf("estimated error = %g\n",acc+fabs(Q)*eps);
	printf("   actual error = %g\n",fabs(Q-exact));


	a=0,b=1,acc=0.001,eps=0.001;
	double f2(double x){calls++; return log(x)/sqrt(x);}; //nested function

	calls=0;
	Q=adapt(f2,a,b,acc,eps);
	exact=-4;
	printf("\nIntegrating log(x)/sqrt(x) from %g to %g\n",a,b);
	printf("              Q = %g\n",Q);
	printf("          exact = %g\n",exact);
	printf("          calls = %d\n",calls);
	printf("estimated error = %g\n",acc+fabs(Q)*eps);
	printf("   actual error = %g\n",fabs(Q-exact));


	a=0,b=1,acc=0.001,eps=0.001;
	double f3(double x){calls++; return sqrt(x);}; //nested function

	calls=0;
	Q=adapt(f3,a,b,acc,eps);
	exact=0.666666;
	printf("\nIntegrating sqrt(x) from %g to %g\n",a,b);
	printf("              Q = %g\n",Q);
	printf("          exact = %g\n",exact);
	printf("          calls = %d\n",calls);
	printf("estimated error = %g\n",acc+fabs(Q)*eps);
	printf("   actual error = %g\n",fabs(Q-exact));


	a=0,b=1,acc=0.001,eps=0.001;
	double f4(double x){calls++; return 4*sqrt(1-(1-x)*(1-x));}; //nested function

	calls=0;
	Q=adapt(f4,a,b,acc,eps);
	exact=3.14159265;
	printf("\nIntegrating 4*sqrt(1-(1-x)^2) from %g to %g\n",a,b);
	printf("              Q = %g\n",Q);
	printf("          exact = %g\n",exact);
	printf("          calls = %d\n",calls);
	printf("estimated error = %g\n",acc+fabs(Q)*eps);
	printf("   actual error = %g\n",fabs(Q-exact));

return 0;
}
