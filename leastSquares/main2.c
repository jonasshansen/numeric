#include<math.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

void lsfit(
        int m, double f(int i,double x),
        gsl_vector* x, gsl_vector* y, gsl_vector* dy,
        gsl_vector* c, gsl_matrix* S);

int main(){
	double x[] = {0.100,0.145,0.211,0.307,0.447,0.649,0.944,1.372,1.995,2.900};
	double y[] = {12.644,9.235,7.377,6.460,5.555,5.896,5.673,6.964,8.896,11.355};
	double dy[] = {0.858,0.359,0.505,0.403,0.683,0.605,0.856,0.351,1.083,1.002};
	int n=sizeof(x)/sizeof(x[0]);

	for(int i=0;i<n;i++)printf("%g %g %g\n",x[i],y[i],dy[i]);
	printf("\n\n");

	gsl_vector* vx = gsl_vector_alloc(n);
	gsl_vector* vy = gsl_vector_alloc(n);
	gsl_vector* vdy = gsl_vector_alloc(n);
	for(int i=0;i<n;i++){
		gsl_vector_set(vx,i,x[i]);
		gsl_vector_set(vy,i,y[i]);
		gsl_vector_set(vdy,i,dy[i]);
		}

	int m=3;
	double funs(int i, double x){
   		switch(i){
   			case 0: return 1/x; break;
   			case 1: return 1;   break;
   			case 2: return x;     break;
   			default: return NAN;
   			}
		}

	gsl_vector* c = gsl_vector_alloc(m);
	gsl_matrix* S = gsl_matrix_alloc(m,m);
	lsfit(m,funs,vx,vy,vdy,c,S);

	gsl_vector* dc = gsl_vector_alloc(m);
	for(int k=0;k<m;k++){
		double skk=gsl_matrix_get(S,k,k);
		gsl_vector_set(dc,k,sqrt(skk));
		}

	double fit(double x){
		double s=0;
		for(int k=0;k<m;k++)s+=gsl_vector_get(c,k)*funs(k,x);
		return s;
		}

	double fit_plus(double x){
		double s=0;
		for(int k=0;k<m;k++)s+=(gsl_vector_get(c,k)+gsl_vector_get(dc,k))*funs(k,x);
		return s;
		}

	double fit_minus(double x){
		double s=0;
		for(int k=0;k<m;k++)s+=(gsl_vector_get(c,k)-gsl_vector_get(dc,k))*funs(k,x);
		return s;
		}

	double dz=(x[n-1]-x[0])/90;
	double z=x[0]-dz/2;
	do{
		printf("%g %g %g %g\n",z,fit(z),fit_plus(z),fit_minus(z));
		z+=dz;
	}while(z<x[n-1]+dz);

return 0;
}
