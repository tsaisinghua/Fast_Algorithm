#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

int main()
{
	int n, k, N = 8;
	double *X_re, *X_im, *x_re, *x_im;
	
	X_re = (double *) malloc(N*sizeof(double));
	X_im = (double *) malloc(N*sizeof(double));	
	x_re = (double *) malloc(N*sizeof(double));
	x_im = (double *) malloc(N*sizeof(double));
	
	for(n=0;n<N;++n)
	{
		x_re[n] = n+1;
		x_im[n] = 0.0;
	}
	
	for(k=0;k<N;++k)
	{
		X_re[k] = 0.0;
		X_im[k] = 0.0;
		// x_n = x_re + i* x_im
		// x_n *(cos(...)- i * sin(...))
		//real:x_re * cos(...) +  x_im * sin(...)
		//image:x_im * cos(...) -  x_re * sin(...)
		for(n=0;n<N;++n)
		{
			X_re[k] += x_re[n] * cos(2*M_PI*k*n/N) +  x_im[n] * sin(2*M_PI*k*n/N);
			X_im[k] += x_im[n] * cos(2*M_PI*k*n/N) -  x_re[n] * sin(2*M_PI*k*n/N);
		}
	}
	
	for(k=0; k<N; ++k)
	{
		if(X_im[k] >= 0)
			printf("%f + %f i\n", X_re[k], X_im[k]);
		else
			printf("%f - %f i\n", X_re[k], -X_im[k]);
	}
	return 0;
} 
