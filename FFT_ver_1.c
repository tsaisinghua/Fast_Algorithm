#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
int main()
{
	int i;
	double y_re[12], y_im[12], x_re[12], x_im[12];
	for(i=0;i<12;++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}
	Fast_Fourier_Transform(y_re, y_im, x_re, x_im, 12);
	for(i=0;i<12;++i)
	{
		printf("%f + %f i\n", y_re[i], y_im[i]);
	}
}
int Fast_Fourier_Transform(double *y_re, double *y_im, double *x_re, double *x_im, int N)
{
	double theta1, w_N_re, w_N_im;
	//double t_r1, t_r2;
	
	if(N==2) 
	{
		// y, y[0] = x[0]+x[1], y[1] = x[0] - x[1]
		y_re[0] = x_re[0] + x_re[1];
		y_im[0] = x_im[0] + x_im[1];
		y_re[1] = x_re[0] - x_re[1]; 
		y_im[1] = x_im[0] - x_im[1];
	}
	else if (N>2 && (N%2)==0)
	{
		int k;
		double *y_even_re, *y_even_im, *y_odd_re, *y_odd_im;
		double *x_even_re, *x_even_im, *x_odd_re, *x_odd_im;
		double w_re, w_im, w_N_re, w_N_im, a, b, temp;
		y_even_re = (double *) malloc( N/2 * sizeof(double));
		y_even_im = (double *) malloc( N/2 * sizeof(double));
		x_even_re = (double *) malloc( N/2 * sizeof(double));
		x_even_im = (double *) malloc( N/2 * sizeof(double));
		y_odd_re = (double *) malloc( N/2 * sizeof(double));
		y_odd_im = (double *) malloc( N/2 * sizeof(double));
		x_odd_re = (double *) malloc( N/2 * sizeof(double));
		x_odd_im = (double *) malloc( N/2 * sizeof(double));
		for(k=0;k<N/2;++k)
		{
			x_even_re[k] = x_re[2*k];
			x_even_im[k] = x_im[2*k];
			x_odd_re[k]  = x_re[2*k+1];
			x_odd_im[k]  = x_im[2*k+1];
		}
		Fast_Fourier_Transform(y_even_re, y_even_im, x_even_re, x_even_im, N/2);
		Fast_Fourier_Transform(y_odd_re, y_odd_im, x_odd_re, x_odd_im, N/2);
		// y_k = even_k + w_N^k odd_k = even_k + (a + bi)
		w_N_re =  cos(2.0*M_PI/N);
		w_N_im = -sin(2.0*M_PI/N);
		w_re   = 1.0;
		w_im   = 0.0;
		for(k=0;k<N/2;++k)
		{
			a = w_re*y_odd_re[k] - w_im*y_odd_im[k];
			b = w_re*y_odd_im[k] + w_im*y_odd_re[k];
			y_re[k]     = y_even_re[k] + a;
			y_im[k]     = y_even_im[k] + b;
			y_re[N/2+k] = y_even_re[k] - a;
			y_im[N/2+k] = y_even_im[k] - b;
			temp = w_re;
			w_re = w_re*w_N_re - w_im*w_N_im;
			w_im = temp*w_N_im + w_im*w_N_re;
		}
		free(y_even_re);
		free(x_even_re);
		free(y_even_im);
		free(x_even_im);
		free(y_odd_re);
		free(y_odd_im);
		free(x_odd_re);
		free(x_odd_im);
	}
	else if(N==3)
	{
		// y, y[0] = x[0]+x[1]+x[2], y[1] = x[0] +W_3 * x[1] + w_3^2 * x[2], y[2] = x[0] +W_3^2 * x[1] + w_3 * x[2]
		theta1 = 2.0*M_PI/3;
		w_N_re =  cos(theta1);
		w_N_im = -sin(theta1);
		y_re[0] = x_re[0] + x_re[1] + x_re[2];
		y_im[0] = x_im[0] + x_im[1] + x_im[2];
		//t_r1 = x_re[1];
		y_re[1] = x_re[0] + w_N_re * (x_re[1] + x_re[2]) - w_N_im * (x_im[1] - x_im[2]);
		y_im[1] = x_im[0] + w_N_re * (x_im[1] + x_im[2]) + w_N_im * (x_re[1] - x_re[2]);
		//t_r2 = x_re[2];
		y_re[2] = x_re[0] + w_N_re * (x_re[1] + x_re[2]) + w_N_im * (x_im[1] - x_im[2]);
		y_im[2] = x_im[0] + w_N_re * (x_im[1] + x_im[2]) - w_N_im * (x_re[1] - x_re[2]);
	}
}
