//DST-I
// 因為 N = ( N + 1 ) * 2 而我們的 fft 只能做 N = 2, 3, 5 的公倍數，所以一開始的 N 要慎選！  

#include <stdio.h> 
#include <stdlib.h> 
#include <omp.h>
#include <math.h>
#include <string.h> /* memset */
#include <unistd.h> /* close */
#include <time.h>
#define DEBUG 0

int Fast_Fourier_Transform(double *y_re, double *y_im, double *x_re, double *x_im, int N);

int main()
{
	int i, p, q, r, N;
	double T1;
	clock_t t1, t2;
	
	N = 8;
	printf("N = %d\n",N);
	N = (N+1)*2;
	
	double *y_re, *y_im, *x_re, *x_im;
	y_re = (double *) malloc( N * sizeof(double));
	y_im = (double *) malloc( N * sizeof(double));
	x_re = (double *) malloc( N * sizeof(double));
	x_im = (double *) malloc( N * sizeof(double));
	memset( x_re, 0, N*sizeof(double) );
	memset( x_im, 0, N*sizeof(double) );
	
	#if DEBUG
	for(i=0;i<N;++i)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	}
	printf("1.================================\n");
	#endif
	
	//initial condition
	for(i=0;i<N/2-1;++i)
	{
		x_re[i+1] = i;
		x_im[i] = 0.0;
	}
	
	#if DEBUG
	for(i=0;i<N;++i)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	}
	printf("2.================================\n");
	#endif
	
	Fast_Fourier_Transform(y_re, y_im, x_re, x_im, N);
	for(i=1;i<N/2;++i)
	{
		printf("%f \n", -y_im[i]);
	}
	
	free(y_re);
	free(x_re);
	free(y_im);
	free(x_im);
	return 0;
}

int Fast_Fourier_Transform(double *y_re, double *y_im, double *x_re, double *x_im, int N)
{
	double theta, theta1;
	double w_re, w_im, w_N_re, w_N_im, w_2_re, w_2_im, w_N_2_re, w_N_2_im; 
	double w_3_re, w_3_im, w_4_re, w_4_im, t_rq, t_rr, t_rs, t_rt, temp, a, b;
	double u1_r, u1_i, u2_r, u2_i, u3_r, u3_i, u4_r, u4_i;
	int k;
	//double t_r1, t_r2;
	
	if(N == 2) 
	{
		// y, y[0] = x[0]+x[1], y[1] = x[0] - x[1]
		y_re[0] = x_re[0] + x_re[1];
		y_im[0] = x_im[0] + x_im[1];
		y_re[1] = x_re[0] - x_re[1];
		y_im[1] = x_im[0] - x_im[1];
	}
	else if(N>2 && (N%2)==0)
	{
		double *z_r, *z_i, *u_r, *u_i;			   // z °μ§1?a!A|s‥i u ??  
		z_r = (double *) malloc(N*sizeof(double)); // z_r:0~N/2-1: even, z_r:N/2~N-1: odd
		z_i = (double *) malloc(N*sizeof(double)); // z_i:0~N/2-1: even, z_i:N/2~N-1: odd 
		u_r = (double *) malloc(N*sizeof(double)); // z_r:0~N/2-1: even, z_r:N/2~N-1: odd
		u_i = (double *) malloc(N*sizeof(double)); // z_i:0~N/2-1: even, z_i:N/2~N-1: odd 
		
		for(k=0;k<N/2;++k)
		{
			z_r[k] = x_re[2*k];
			z_i[k] = x_im[2*k];
			z_r[N/2+k]  = x_re[2*k+1];
			z_i[N/2+k]  = x_im[2*k+1];
		}
		Fast_Fourier_Transform(u_r, u_i, z_r, z_i, N/2);
		Fast_Fourier_Transform(u_r+N/2, u_i+N/2, z_r+N/2, z_i+N/2, N/2);
		theta1 = 2.0*M_PI/N;
		w_N_re =  cos(theta1);
		w_N_im = -sin(theta1);
		w_re   = 1.0;
		w_im   = 0.0;
		for(k=0;k<N/2;++k)
		{
			a = w_re*u_r[N/2+k] - w_im*u_i[N/2+k];
			b = w_re*u_i[N/2+k] + w_im*u_r[N/2+k];
			y_re[k]     = u_r[k] + a;
			y_im[k]     = u_i[k] + b;
			y_re[N/2+k] = u_r[k] - a;
			y_im[N/2+k] = u_i[k] - b;
			temp = w_re;
			w_re = w_re*w_N_re - w_im*w_N_im;
			w_im = temp*w_N_im + w_im*w_N_re;
		}
		free(u_r); free(u_i); free(z_r); free(z_i);
		
	}	
	else if(N == 3)
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
	else if(N>3 && (N%3)==0)
	{
		double *z_r, *z_i, *u_r, *u_i;
		z_r = (double *) malloc(N*sizeof(double)); // z_r:0~N/2-1: even, z_r:N/2~N-1: odd
		z_i = (double *) malloc(N*sizeof(double)); // z_i:0~N/2-1: even, z_i:N/2~N-1: odd 
		u_r = (double *) malloc(N*sizeof(double)); // z_r:0~N/2-1: even, z_r:N/2~N-1: odd
		u_i = (double *) malloc(N*sizeof(double)); // z_i:0~N/2-1: even, z_i:N/2~N-1: odd 
		
		for(k=0;k<N/3;++k)
		{
			z_r[k] = x_re[3*k];
			z_i[k] = x_im[3*k];
			z_r[N/3+k]  = x_re[3*k+1];
			z_i[N/3+k]  = x_im[3*k+1];
			z_r[2*N/3+k]  = x_re[3*k+2];
			z_i[2*N/3+k]  = x_im[3*k+2];
		}
		Fast_Fourier_Transform(u_r, u_i, z_r, z_i, N/3);
		Fast_Fourier_Transform(u_r+N/3, u_i+N/3, z_r+N/3, z_i+N/3, N/3);
		Fast_Fourier_Transform(u_r+2*N/3, u_i+2*N/3, z_r+2*N/3, z_i+2*N/3, N/3);
		
		w_N_re =  cos(2.0*M_PI/3);
		w_N_im = -sin(2.0*M_PI/3);
		
		
		for(k=0;k<N/3;++k)
		{
			theta = 2.0*k*M_PI/N;
			w_re  =  cos(theta);
			w_im  = -sin(theta);
			w_2_re = w_re * w_re - w_im * w_im;
			w_2_im = w_re * w_im + w_re * w_im;
					
			u1_r = w_re*u_r[N/3+k] - w_im*u_i[N/3+k];
			u1_i = w_re*u_i[N/3+k] + w_im*u_r[N/3+k];
			// (w_re + i w_im)^2 = w_re*w_re - w_im*w_im + i (2*w_re*w_im)
			u2_r = w_2_re*u_r[2*N/3+k] - w_2_im*u_i[2*N/3+k];
			u2_i = w_2_re*u_i[2*N/3+k] + w_2_im*u_r[2*N/3+k];
			a = -0.5; b = -sqrt(3)/2;
			y_re[k]       = u_r[k] + u1_r + u2_r;
			y_im[k]       = u_i[k] + u1_i + u2_i;
			y_re[N/3+k]   = u_r[k] + (u1_r*a-u1_i*b) + (u2_r*a+u2_i*b);
			y_im[N/3+k]   = u_i[k] + (u1_i*a+u1_r*b) + (u2_i*a-u2_r*b);
			y_re[2*N/3+k] = u_r[k] + (u1_r*a+u1_i*b) + (u2_r*a-u2_i*b);
			y_im[2*N/3+k] = u_i[k] + (u1_i*a-u1_r*b) + (u2_i*a+u2_r*b);
			/* 
			temp = w_re;
			w_re = w_re*w_N_re - w_im*w_N_im;
			w_im = temp*w_N_im + w_im*w_N_re;
			*/ 
		}
		free(u_r); free(u_i); free(z_r); free(z_i);


	}
	else if(N == 5) 
	{
		
		double w_N_2_re, w_N_2_im;
		theta1 = 2.0*M_PI/5;
		w_N_re =  cos(theta1);
		w_N_im = -sin(theta1);
		// y, y[0] = x[0] +     x[1]     +     x[2]     +     x[3]     +     x[4], 
		//	  y[1] = x[0] + w_5   * x[1] + w_5^2 * x[2] + w_5^3 * x[3] + w_5^4 * x[4]
		//	  y[2] = x[0] + w_5^2 * x[1] + w_5^4 * x[2] + w_5   * x[3] + w_5^3 * x[4]
		//	  y[3] = x[0] + w_5^3 * x[1] + w_5   * x[2] + w_5^4 * x[3] + w_5^2 * x[4]
		//	  y[4] = x[0] + w_5^4 * x[1] + w_5^3 * x[2] + w_5^2 * x[3] + w_5   * x[4]
		
		w_N_2_re = w_N_re * w_N_re - w_N_im * w_N_im;
		w_N_2_im = 2 * w_N_re * w_N_im;
		
		y_re[0] = x_re[0] + x_re[1] + x_re[2] + x_re[3] + x_re[4];
		y_im[0] = x_im[0] + x_im[1] + x_im[2] + x_im[3] + x_im[4];
		
		y_re[1] = x_re[0] + w_N_re * (x_re[1] + x_re[4]) - w_N_im * (x_im[1] - x_im[4]) + w_N_2_re * (x_re[2] + x_re[3]) - w_N_2_im * (x_im[2] - x_im[3]);
		y_im[1] = x_im[0] + w_N_re * (x_im[1] + x_im[4]) + w_N_im * (x_re[1] - x_re[4]) + w_N_2_re * (x_im[2] + x_im[3]) + w_N_2_im * (x_re[2] - x_re[3]);
		
		y_re[2] = x_re[0] + w_N_re * (x_re[2] + x_re[3]) - w_N_im * (x_im[3] - x_im[2]) + w_N_2_re * (x_re[1] + x_re[4]) - w_N_2_im * (x_im[1] - x_im[4]);
		y_im[2] = x_im[0] + w_N_re * (x_im[2] + x_im[3]) + w_N_im * (x_re[3] - x_re[2]) + w_N_2_re * (x_im[1] + x_im[4]) + w_N_2_im * (x_re[1] - x_re[4]);	
		
		y_re[3] = x_re[0] + w_N_re * (x_re[2] + x_re[3]) - w_N_im * (x_im[2] - x_im[3]) + w_N_2_re * (x_re[4] + x_re[1]) - w_N_2_im * (x_im[4] - x_im[1]);
		y_im[3] = x_im[0] + w_N_re * (x_im[2] + x_im[3]) + w_N_im * (x_re[2] - x_re[3]) + w_N_2_re * (x_im[4] + x_im[1]) + w_N_2_im * (x_re[4] - x_re[1]);
		
		y_re[4] = x_re[0] + w_N_re * (x_re[1] + x_re[4]) - w_N_im * (x_im[4] - x_im[1]) + w_N_2_re * (x_re[3] + x_re[2]) - w_N_2_im * (x_im[3] - x_im[2]);
		y_im[4] = x_im[0] + w_N_re * (x_im[1] + x_im[4]) + w_N_im * (x_re[4] - x_re[1]) + w_N_2_re * (x_im[3] + x_im[2]) + w_N_2_im * (x_re[3] - x_re[2]);		
	}
	else if(N>5 && (N%5)==0)
	{
		double *z_r, *z_i, *u_r, *u_i;
		z_r = (double *) malloc(N*sizeof(double)); // z_r:0~N/5-1: 0th, z_r:N/5~2*N/5-1: 1 th
		z_i = (double *) malloc(N*sizeof(double)); // z_i:0~N/5-1: 0th, z_i:N/5~2*N/5-1: 1 th
		u_r = (double *) malloc(N*sizeof(double)); // z_r:0~N/5-1: 0th, z_r:N/5~2*N/5-1: 1 th
		u_i = (double *) malloc(N*sizeof(double)); // z_i:0~N/5-1: 0th, z_i:N/5~2*N/5-1: 1 th
		
		for(k=0;k<N/5;++k)
		{
			z_r[k] = x_re[5*k];
			z_i[k] = x_im[5*k];
			z_r[N/5+k] = x_re[5*k+1];
			z_i[N/5+k] = x_im[5*k+1];
			z_r[2*N/5+k] = x_re[5*k+2];
			z_i[2*N/5+k] = x_im[5*k+2];
			z_r[3*N/5+k] = x_re[5*k+3];
			z_i[3*N/5+k] = x_im[5*k+3];
			z_r[4*N/5+k] = x_re[5*k+4];
			z_i[4*N/5+k] = x_im[5*k+4];
		}
		Fast_Fourier_Transform(u_r, u_i, z_r, z_i, N/5);
		Fast_Fourier_Transform(u_r+N/5, u_i+N/5, z_r+N/5, z_i+N/5, N/5);
		Fast_Fourier_Transform(u_r+2*N/5, u_i+2*N/5, z_r+2*N/5, z_i+2*N/5, N/5);
		Fast_Fourier_Transform(u_r+3*N/5, u_i+3*N/5, z_r+3*N/5, z_i+3*N/5, N/5);
		Fast_Fourier_Transform(u_r+4*N/5, u_i+4*N/5, z_r+4*N/5, z_i+4*N/5, N/5);
		w_N_re =  cos(2.0*M_PI/5);
		w_N_im = -sin(2.0*M_PI/5);
		
		for(k=0;k<N/5;++k)
		{
			theta = 2.0*k*M_PI/N;
			w_re   =  cos(theta);
			w_im   = -sin(theta);
			
			w_2_re = w_re * w_re - w_im * w_im;
			w_2_im = w_re * w_im + w_re * w_im;
			w_3_re = w_re * w_re * w_re - 3 * w_re * w_im * w_im;
			w_3_im = 3 * w_re * w_re * w_im - w_im * w_im * w_im;
			w_4_re = w_2_re * w_2_re - w_2_im * w_2_im;
			w_4_im = 2 * w_2_re * w_2_im;
			
			u1_r = w_re * u_r[N/5+k] - w_im * u_i[N/5+k];
			u1_i = w_re * u_i[N/5+k] + w_im * u_r[N/5+k];
			
			u2_r = w_2_re * u_r[2*N/5+k] - w_2_im*u_i[2*N/5+k];
			u2_i = w_2_re * u_i[2*N/5+k] + w_2_im*u_r[2*N/5+k];
				
			u3_r = w_3_re * u_r[3*N/5+k] - w_3_im * u_i[3*N/5+k];
			u3_i = w_3_re * u_i[3*N/5+k] + w_3_im * u_r[3*N/5+k];
			
			u4_r = w_4_re * u_r[4*N/5+k] - w_4_im * u_i[4*N/5+k];
			u4_i = w_4_re * u_i[4*N/5+k] + w_4_im * u_r[4*N/5+k];
			
	
			w_N_2_re = w_N_re * w_N_re - w_N_im * w_N_im;
			w_N_2_im = 2 * w_N_re * w_N_im;	
			
			y_re[k] 	  = u_r[k] + u1_r + u2_r + u3_r + u4_r;
			y_im[k] 	  = u_i[k] + u1_i + u2_i + u3_i + u4_i;
			
			y_re[N/5+k]   = u_r[k] + w_N_re * (u1_r + u4_r) - w_N_im * (u1_i - u4_i) + w_N_2_re * (u2_r + u3_r) - w_N_2_im * (u2_i - u3_i);
			y_im[N/5+k]   = u_i[k] + w_N_re * (u1_i + u4_i) + w_N_im * (u1_r - u4_r) + w_N_2_re * (u2_i + u3_i) + w_N_2_im * (u2_r - u3_r);
			
			y_re[2*N/5+k] = u_r[k] + w_N_re * (u2_r + u3_r) - w_N_im * (u3_i - u2_i) + w_N_2_re * (u1_r + u4_r) - w_N_2_im * (u1_i - u4_i);
			y_im[2*N/5+k] = u_i[k] + w_N_re * (u2_i + u3_i) + w_N_im * (u3_r - u2_r) + w_N_2_re * (u1_i + u4_i) + w_N_2_im * (u1_r - u4_r);

			y_re[3*N/5+k] = u_r[k] + w_N_re * (u2_r + u3_r) - w_N_im * (u2_i - u3_i) + w_N_2_re * (u4_r + u1_r) - w_N_2_im * (u4_i - u1_i);
			y_im[3*N/5+k] = u_i[k] + w_N_re * (u2_i + u3_i) + w_N_im * (u2_r - u3_r) + w_N_2_re * (u4_i + u1_i) + w_N_2_im * (u4_r - u1_r);
			
			y_re[4*N/5+k] = u_r[k] + w_N_re * (u1_r + u4_r) - w_N_im * (u4_i - u1_i) + w_N_2_re * (u3_r + u2_r) - w_N_2_im * (u3_i - u2_i);
			y_im[4*N/5+k] = u_i[k] + w_N_re * (u1_i + u4_i) + w_N_im * (u4_r - u1_r) + w_N_2_re * (u3_i + u2_i) + w_N_2_im * (u3_r - u2_r);
			/* 
			temp = w_re;
			w_re = w_re*w_N_re - w_im*w_N_im;
			w_im = temp*w_N_im + w_im*w_N_re;
			*/ 
		}
		free(u_r); free(u_i); free(z_r); free(z_i);
	}
	else
	{
		printf("We didn't do this!\n");
	}
}
