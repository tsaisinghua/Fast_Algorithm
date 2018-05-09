#include <stdio.h>
#include <stdlib.h> 
#include <omp.h>
#include <math.h>
#include <time.h>
#define DEBUG 0

int main()
{
	int i, p, q, r, N;
	double T1;
	clock_t t1, t2;
	/*	
	printf("input 2^p 3^q 5^r : (p, q, r) = ");
	scanf("%d %d %d", &p, &q, &r);
	
	srand(time(NULL));
	p = rand() % 5;
	q = rand() % 5;
	r = rand() % 5;
	*/
	N = 2*2*3*3*3*5*5*5;
	
	//N = (int)(pow(2, p) * pow(3, q) * pow(5, r));
	printf("N=%d\n",N);
	
	double y_re[N], y_im[N], x_re[N], x_im[N];
	for(i=0;i<N;++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}
	
	t1 = clock();
	Fast_Fourier_Transform(y_re, y_im, x_re, x_im, N);
	t2 = clock();
	T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
	printf("FFT_ver1 of %d elements: %f secs\n",N, T1);
		
	#if DEBUG
	for(i=0;i<N;++i)
	{
		printf("%f + %f i\n", y_re[i], y_im[i]);
	}
	#endif
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
		
		//#pragma omp parallel for 
		for(k=0;k<N/2;++k)
		{
			//#pragma omp ordered
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
	else if(N>3 && (N%3)==0)
	{
		int k;
		double *y_30_re, *y_30_im, *y_31_re, *y_31_im, *y_32_re, *y_32_im;
		double *x_30_re, *x_30_im, *x_31_re, *x_31_im, *x_32_re, *x_32_im;
		double w_re, w_im, w_N_re, w_N_im, w_2_re, w_2_im, t_rq, t_rr, temp, theta;
		y_30_re = (double *) malloc( N/3 * sizeof(double));
		y_30_im = (double *) malloc( N/3 * sizeof(double));
		x_30_re = (double *) malloc( N/3 * sizeof(double));
		x_30_im = (double *) malloc( N/3 * sizeof(double));
		y_31_re = (double *) malloc( N/3 * sizeof(double));
		y_31_im = (double *) malloc( N/3 * sizeof(double));
		x_31_re = (double *) malloc( N/3 * sizeof(double));
		x_31_im = (double *) malloc( N/3 * sizeof(double));
		y_32_re = (double *) malloc( N/3 * sizeof(double));
		y_32_im = (double *) malloc( N/3 * sizeof(double));
		x_32_re = (double *) malloc( N/3 * sizeof(double));
		x_32_im = (double *) malloc( N/3 * sizeof(double));
		
		//#pragma omp parallel for ordered
		for(k=0;k<N/3;++k)
		{
			//#pragma omp ordered
			x_30_re[k] = x_re[3*k];
			x_30_im[k] = x_im[3*k];
			x_31_re[k] = x_re[3*k+1];
			x_31_im[k] = x_im[3*k+1];
			x_32_re[k] = x_re[3*k+2];
			x_32_im[k] = x_im[3*k+2];
		}
		Fast_Fourier_Transform(y_30_re, y_30_im, x_30_re, x_30_im, N/3);
		Fast_Fourier_Transform(y_31_re, y_31_im, x_31_re, x_31_im, N/3);
		Fast_Fourier_Transform(y_32_re, y_32_im, x_32_re, x_32_im, N/3);
		// y_k = even_k + w_N^k odd_k = even_k + (a + bi)
		w_N_re =  cos(2.0*M_PI/3);
		w_N_im = -sin(2.0*M_PI/3);
		
		//#pragma omp parallel for ordered
		for(k=0;k<N/3;++k)
		{
			//#pragma omp ordered
			theta = 2.0*k*M_PI/N;
			w_re  =  cos(theta);
			w_im  = -sin(theta);
			
			w_2_re = w_re * w_re - w_im * w_im;
			w_2_im = w_re * w_im + w_re * w_im;
			// multiply (w_re + w_im * i) on x[q]
			t_rq = y_31_re[k];
			y_31_re[k] = w_re * y_31_re[k] - w_im * y_31_im[k];
			y_31_im[k] = w_re * y_31_im[k] + w_im * t_rq;
			// multiply (w_re + w_im * i)^2 = w_re - w_im * i on x[r]
			t_rr = y_32_re[k];
			y_32_re[k] = w_2_re * y_32_re[k] - w_2_im * y_32_im[k];
			y_32_im[k] = w_2_re * y_32_im[k] + w_2_im * t_rr;
				
				
			y_re[k] 	  = y_30_re[k] + y_31_re[k] + y_32_re[k];
			y_im[k] 	  = y_30_im[k] + y_31_im[k] + y_32_im[k]; 
			y_re[N/3+k]   = y_30_re[k] + w_N_re * (y_31_re[k] +  y_32_re[k]) - w_N_im * (y_31_im[k] - y_32_im[k]);
			y_im[N/3+k]   = y_30_im[k] + w_N_re * (y_31_im[k] +  y_32_im[k]) + w_N_im * (y_31_re[k] - y_32_re[k]);
			y_re[2*N/3+k] = y_30_re[k] + w_N_re * (y_31_re[k] +  y_32_re[k]) + w_N_im * (y_31_im[k] - y_32_im[k]);
			y_im[2*N/3+k] = y_30_im[k] + w_N_re * (y_31_im[k] +  y_32_im[k]) - w_N_im * (y_31_re[k] - y_32_re[k]);
			/*
			temp = w_re;
			w_re = w_re*w_N_re - w_im*w_N_im;
			w_im = temp*w_N_im + w_im*w_N_re;
			*/
		}
		free(y_30_re);free(y_30_im);free(x_30_re);
		free(x_30_im);free(y_31_re);free(y_31_im);
		free(x_31_re);free(x_31_im);free(y_32_re);
		free(y_32_im);free(x_32_re);free(x_32_im);
	}
	else if(N==5) 
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
		int k;
		double *y_50_re, *y_50_im, *y_51_re, *y_51_im, *y_52_re, *y_52_im, *y_53_re, *y_53_im, *y_54_re, *y_54_im;
		double *x_50_re, *x_50_im, *x_51_re, *x_51_im, *x_52_re, *x_52_im, *x_53_re, *x_53_im, *x_54_re, *x_54_im;
		double w_re, w_im, w_N_re, w_N_im, w_2_re, w_2_im, w_N_2_re, w_N_2_im, w_3_re, w_3_im, w_4_re, w_4_im, t_rq, t_rr, t_rs, t_rt, temp, theta;
		y_50_re = (double *) malloc( N/5 * sizeof(double));
		y_50_im = (double *) malloc( N/5 * sizeof(double));
		x_50_re = (double *) malloc( N/5 * sizeof(double));
		x_50_im = (double *) malloc( N/5 * sizeof(double));
		y_51_re = (double *) malloc( N/5 * sizeof(double));
		y_51_im = (double *) malloc( N/5 * sizeof(double));
		x_51_re = (double *) malloc( N/5 * sizeof(double));
		x_51_im = (double *) malloc( N/5 * sizeof(double));
		y_52_re = (double *) malloc( N/5 * sizeof(double));
		y_52_im = (double *) malloc( N/5 * sizeof(double));
		x_52_re = (double *) malloc( N/5 * sizeof(double));
		x_52_im = (double *) malloc( N/5 * sizeof(double));
		y_53_re = (double *) malloc( N/5 * sizeof(double));
		y_53_im = (double *) malloc( N/5 * sizeof(double));
		x_53_re = (double *) malloc( N/5 * sizeof(double));
		x_53_im = (double *) malloc( N/5 * sizeof(double));
		y_54_re = (double *) malloc( N/5 * sizeof(double));
		y_54_im = (double *) malloc( N/5 * sizeof(double));
		x_54_re = (double *) malloc( N/5 * sizeof(double));
		x_54_im = (double *) malloc( N/5 * sizeof(double));
		
		//#pragma omp parallel for ordered
		for(k=0;k<N/5;++k)
		{
			//#pragma omp ordered
			x_50_re[k] = x_re[5*k];
			x_50_im[k] = x_im[5*k];
			x_51_re[k] = x_re[5*k+1];
			x_51_im[k] = x_im[5*k+1];
			x_52_re[k] = x_re[5*k+2];
			x_52_im[k] = x_im[5*k+2];
			x_53_re[k] = x_re[5*k+3];
			x_53_im[k] = x_im[5*k+3];
			x_54_re[k] = x_re[5*k+4];
			x_54_im[k] = x_im[5*k+4];
		}
		Fast_Fourier_Transform(y_50_re, y_50_im, x_50_re, x_50_im, N/5);
		Fast_Fourier_Transform(y_51_re, y_51_im, x_51_re, x_51_im, N/5);
		Fast_Fourier_Transform(y_52_re, y_52_im, x_52_re, x_52_im, N/5);
		Fast_Fourier_Transform(y_53_re, y_53_im, x_53_re, x_53_im, N/5);
		Fast_Fourier_Transform(y_54_re, y_54_im, x_54_re, x_54_im, N/5);
		// y_k = even_k + w_N^k odd_k = even_k + (a + bi)
		w_N_re =  cos(2.0*M_PI/5);
		w_N_im = -sin(2.0*M_PI/5);
		
		//#pragma omp parallel for ordered
		for(k=0;k<N/5;++k)
		{
			//#pragma omp ordered
			theta = 2.0*k*M_PI/N;
			w_re   =  cos(theta);
			w_im   = -sin(theta);
			
			w_2_re = w_re * w_re - w_im * w_im;
			w_2_im = w_re * w_im + w_re * w_im;
			w_3_re = w_re * w_re * w_re - 3 * w_re * w_im * w_im;
			w_3_im = 3 * w_re * w_re * w_im - w_im * w_im * w_im;
			w_4_re = w_2_re * w_2_re - w_2_im * w_2_im;
			w_4_im = 2 * w_2_re * w_2_im;
			
			// multiply (w_re + w_im * i) on x[q]
			t_rq = y_51_re[k];
			y_51_re[k] = w_re * y_51_re[k] - w_im * y_51_im[k];
			y_51_im[k] = w_re * y_51_im[k] + w_im * t_rq;
			// multiply (w_re + w_im * i)^2 = w_re - w_im * i on x[r]
			t_rr = y_52_re[k];
			y_52_re[k] = w_2_re * y_52_re[k] - w_2_im*y_52_im[k];
			y_52_im[k] = w_2_re * y_52_im[k] + w_2_im*t_rr;
				
			t_rs = y_53_re[k];
			y_53_re[k] = w_3_re * y_53_re[k] - w_3_im * y_53_im[k];
			y_53_im[k] = w_3_re * y_53_im[k] + w_3_im * t_rs;
			
			t_rt = y_54_re[k];
			y_54_re[k] = w_4_re * y_54_re[k] - w_4_im * y_54_im[k];
			y_54_im[k] = w_4_re * y_54_im[k] + w_4_im * t_rt;
	//=============================================================================================================================================================
			
			w_N_2_re = w_N_re * w_N_re - w_N_im * w_N_im;
			w_N_2_im = 2 * w_N_re * w_N_im;	
			
			y_re[k] 	  = y_50_re[k] + y_51_re[k] + y_52_re[k] + y_53_re[k] + y_54_re[k];
			y_im[k] 	  = y_50_im[k] + y_51_im[k] + y_52_im[k] + y_53_im[k] + y_54_im[k]; 
			
			y_re[N/5+k]   = y_50_re[k] + w_N_re * (y_51_re[k] + y_54_re[k]) - w_N_im * (y_51_im[k] - y_54_im[k]) + w_N_2_re * (y_52_re[k] + y_53_re[k]) - w_N_2_im * (y_52_im[k] - y_53_im[k]);
			y_im[N/5+k]   = y_50_im[k] + w_N_re * (y_51_im[k] + y_54_im[k]) + w_N_im * (y_51_re[k] - y_54_re[k]) + w_N_2_re * (y_52_im[k] + y_53_im[k]) + w_N_2_im * (y_52_re[k] - y_53_re[k]);
			
			y_re[2*N/5+k] = y_50_re[k] + w_N_re * (y_52_re[k] + y_53_re[k]) - w_N_im * (y_53_im[k] - y_52_im[k]) + w_N_2_re * (y_51_re[k] + y_54_re[k]) - w_N_2_im * (y_51_im[k] - y_54_im[k]);
			y_im[2*N/5+k] = y_50_im[k] + w_N_re * (y_52_im[k] + y_53_im[k]) + w_N_im * (y_53_re[k] - y_52_re[k]) + w_N_2_re * (y_51_im[k] + y_54_im[k]) + w_N_2_im * (y_51_re[k] - y_54_re[k]);

			y_re[3*N/5+k] = y_50_re[k] + w_N_re * (y_52_re[k] + y_53_re[k]) - w_N_im * (y_52_im[k] - y_53_im[k]) + w_N_2_re * (y_54_re[k] + y_51_re[k]) - w_N_2_im * (y_54_im[k] - y_51_im[k]);
			y_im[3*N/5+k] = y_50_im[k] + w_N_re * (y_52_im[k] + y_53_im[k]) + w_N_im * (y_52_re[k] - y_53_re[k]) + w_N_2_re * (y_54_im[k] + y_51_im[k]) + w_N_2_im * (y_54_re[k] - y_51_re[k]);
			
			y_re[4*N/5+k] = y_50_re[k] + w_N_re * (y_51_re[k] + y_54_re[k]) - w_N_im * (y_54_im[k] - y_51_im[k]) + w_N_2_re * (y_53_re[k] + y_52_re[k]) - w_N_2_im * (y_53_im[k] - y_52_im[k]);
			y_im[4*N/5+k] = y_50_im[k] + w_N_re * (y_51_im[k] + y_54_im[k]) + w_N_im * (y_54_re[k] - y_51_re[k]) + w_N_2_re * (y_53_im[k] + y_52_im[k]) + w_N_2_im * (y_53_re[k] - y_52_re[k]);
		}
		free(y_50_re);free(y_50_im);free(x_50_re);free(x_50_im);
		free(y_51_re);free(y_51_im);free(x_51_re);free(x_51_im);		
		free(y_52_re);free(y_52_im);free(x_52_re);free(x_52_im);
		free(y_53_re);free(y_53_im);free(x_53_re);free(x_53_im);		
		free(y_54_re);free(y_54_im);free(x_54_re);free(x_54_im);
	}
	else
	{
		printf("We didn't do this!\n");
	}
}
