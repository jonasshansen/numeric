#include<stdlib.h>
#include<assert.h>

double linterp(int n, double *x, double *y, double z){
    assert(n>1 && z>=x[0] && z<=x[n-1]); //check that the argument values are correct

    int i = 0, j = n-1;
    while(j-i>1){int m=(i+j)/2; if(z>x[m]) i=m; else j=m;} // binary search algorithm

    double dx_i = x[i+1] - x[i];
    double dy_i = y[i+1] - y[i];
    double p_i = dy_i/dx_i;

    return y[i] + p_i*(z - x[i]);
}

