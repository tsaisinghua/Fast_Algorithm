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
	else if(N>2 && (N%2)==0)
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
	else(N>3 && (N%3)==0)
	{
		int k;
		double *y_30_re, *y_30_im, *y_31_re, *y_31_im, *y_32_re, *y_32_im;
		double *x_30_re, *x_30_im, *x_31_re, *x_31_im, *x_32_re, *x_32_im;
		double w_re, w_im, w_N_re, w_N_im, a, b, temp;
		y_30_re = (double *) malloc( N/3 * sizeof(double));
		y_30_im = (double *) malloc( N/3 * sizeof(double));
		x_30_re = (double *) malloc( N/3 * sizeof(double));
		x_30_im = (double *) malloc( N/3 * sizeof(double));
		y_31_re = (double *) malloc( N/3 * sizeof(double));
		y_31_im = (double *) malloc( N/3 * sizeof(double));
		x_31_re = (double *) malloc( N/3 * sizeof(double));
		x_31_im = (double *) malloc( N/3 * sizeof(double));
		for(k=0;k<N/3;++k)
		{
			x_30_re[k] = x_re[3*k];
			x_30_im[k] = x_im[3*k];
			x_31_re[k]  = x_re[3*k+1];
			x_31_im[k]  = x_im[3*k+1];
			x_32_re[k]  = x_re[3*k+2];
			x_32_im[k]  = x_im[3*k+2];
		}
		Fast_Fourier_Transform(y_30_re, y_30_im, x_30_re, x_30_im, N/3);
		Fast_Fourier_Transform(y_31_re, y_31_im, x_31_re, x_31_im, N/3);
		Fast_Fourier_Transform(y_32_re, y_32_im, x_32_re, x_32_im, N/3);
		// y_k = even_k + w_N^k odd_k = even_k + (a + bi)
		w_N_re =  cos(2.0*M_PI/N);
		w_N_im = -sin(2.0*M_PI/N);
		w_re   = 1.0;
		w_im   = 0.0;
		for(k=0;k<N/3;++k)
		{
			/*
			a = w_re*y_odd_re[k] - w_im*y_odd_im[k];
			b = w_re*y_odd_im[k] + w_im*y_odd_re[k];
			y_re[k]     = y_even_re[k] + a;
			y_im[k]     = y_even_im[k] + b;
			y_re[N/2+k] = y_even_re[k] - a;
			y_im[N/2+k] = y_even_im[k] - b;
			temp = w_re;
			w_re = w_re*w_N_re - w_im*w_N_im;
			w_im = temp*w_N_im + w_im*w_N_re;
			*/
			
			// multiply (w_re + w_im * i) on x[q]
			t_rq = x_re[q];
			x_re[q] = w_re * x_re[q] - w_im * x_im[q];
			x_im[q] = w_re * x_im[q] + w_im * t_rq;
			// multiply (w_re + w_im * i)^2 = w_re - w_im * i on x[r]
			t_rr = x_re[r];
			x_re[r] = (w_re * w_re - w_im * w_im) * x_re[r] - (w_re * w_im + w_re * w_im)*x_im[r];
			x_im[r] = (w_re * w_re - w_im * w_im) * x_im[r] + (w_re * w_im + w_re * w_im)*t_rr;
				
			y_re[k]     = y_30_re[k] + y_31_re[k] + y_32_re[k];
			y_im[k]     = y_30_im[k] + y_31_im[k] + y_32_im[k];
			y_re[N/3+k] = y_30_re[k] + w_re * y_31_re[k] - w_im * y_31_im[k] + (w_re * w_re - w_im * w_im) * y_32_re[k] - (w_re * w_im + w_re * w_im) * y_32_im[k];
			y_im[N/3+k] = y_30_im[k] + w_re * y_31_im[k] - w_im * y_31_re[k] + (w_re * w_re - w_im * w_im) * y_32_im[k] + (w_re * w_im + w_re * w_im) * y_32_re[k];
			y_re[2*N/3+k] = y_30_re[k] + ;
			y_im[2*N/3+k] = y_30_im[k] + ;
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
}
