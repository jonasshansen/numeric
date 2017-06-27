#include<math.h>
#include<stdlib.h>
#include<stdio.h>

void rkstep12(
        void f(int n,double x,double*y,double*dydx),
        int n, double x, double* y, double h, double* yh, double* dy);

int ode_driver(void f(int n,double x,double*y,double*dydx),
int n,double*xlist,double**ylist,
double b,double h,double acc,double eps,int max){

  int i,k=0; double x,*y,s,err,normy,tol,a=xlist[0],yh[n],dy[n];

  while(xlist[k]<b){

    x=xlist[k], y=ylist[k]; if(x+h>b) h=b-x;

    rkstep12(f,n,x,y,h,yh,dy);

    s=0; for(i=0;i<n;i++) s+=dy[i]*dy[i]; err  =sqrt(s);
    s=0; for(i=0;i<n;i++) s+=yh[i]*yh[i]; normy=sqrt(s);

    tol=(normy*eps+acc)*sqrt(h/(b-a));
    if(err<tol){ // accept step and continue
      k++; if(k>max-1) return -k;
      xlist[k]=x+h; for(i=0;i<n;i++)ylist[k][i]=yh[i];
      }
    if(err>0) h*=pow(tol/err,0.25)*0.95; else h*=2;
  }
  return k+1; // return the number of entries in xlist/ylist
}
