#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int ode_driver(
        void f(int n, double x, double*y, double*dydx),
        int n, double* xlist, double** ylist,
        double b, double h, double acc, double eps, int max);

void f(int n, double x, double* y, double* dydx){
	dydx[0]=y[1];
	dydx[1]=-y[0];
	return;
}


int main(){
	int n=2;
	int max=1000;
	double*xlist=(double*)calloc(max,sizeof(double));
	double**ylist=(double**)calloc(max,sizeof(double*));
	for(int i=0;i<max;i++) ylist[i]=(double*)calloc(n,sizeof(double));

	double pi=atan(1.0)*4;
	double a=0, b=8*pi, h=0.1, acc=0.01, eps=0.01;
	xlist[0]=a; ylist[0][0]=0; ylist[0][1]=1;
	int k = ode_driver(f,n,xlist,ylist,b,h,acc,eps,max);
	if(k<0)printf("max steps reached in ode_driver\n");
	printf("# m=0, S=4\n");
	for(int i=0;i<k;i++)printf("%g %g\n",xlist[i],ylist[i][0]);
	printf("# m=1, S=0\n");
	for(int i=0;i<k;i++)printf("%g %g\n",xlist[i],sin(xlist[i]));
  
	return 0;
}
