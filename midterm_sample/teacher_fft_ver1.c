#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <time.h>

int main()
{
	int i, N, N2, K, *p;
	double *x_r, *x_i, *y_r, *y_i; // y = fft(x);
	clock_t t1, t2;
	
	printf("N=");
	scanf("%d", &N);
	printf("N=%d\n", N);
	
	x_r = (double *) malloc(N * sizeof(double));
	x_i = (double *) malloc(N * sizeof(double));
	y_r = (double *) malloc(N * sizeof(double));
	y_i = (double *) malloc(N * sizeof(double));
	
	for(i=0;i<N;++i)
	{
		x_r[i] = i;
		x_i[i] = 0;
	}
	t1 = clock();
	FFT_general(x_r, x_i, y_r, y_i, N);
	t2 = clock();
	printf("time=%f\n",1.0*(t2-t1)/CLOCKS_PER_SEC);
	system("pause");
	for(i=0;i<N;++i)
	{
		printf("y[%d] : %f + %f i\n", i, y_r[i], y_i[i]);
	}

	free(x_r);
	free(x_i);
	free(y_r);
	free(y_i);
	return 0;
} 
int FFT_general(double *x_r, double *x_i, double *y_r, double *y_i, int N)
{
	int k;
	double w_re, w_im, w_N_re, w_N_im, temp, u1_r, u1_i, u2_r, u2_i, a, b, c, d;
	/*
	p = (int *) malloc(N * sizeof(int));	
	N2 = N;
	K  = 0;
	while(N2>1)
	{
		if(N2 % 2 == 0) { N2 /= 2; p[K] = 2; }
		else if(N2 % 3 == 0) { N2 /= 3; p[K] = 3; }
		else if(N2 % 5 == 0) { N2 /= 5; p[K] = 5; }
		K++;
	}
	printf("K=%d\n", K); 
	for(i=0;i<K;++i)
	{
		printf("%d ", p[i]);
	}
	*/
	if(N == 2)
	{
		y_r[0] = x_r[0] + x_r[1];
		y_i[0] = x_i[0] + x_i[1];
		y_r[1] = x_r[0] - x_r[1]; 
		y_i[1] = x_i[0] - x_i[1];
	} else if (N == 3)
	{
		// w_3    = cos(2pi/3)-i sin(2pi/3)=(-1/2 - i*sqrt(3)/2) = a+bi
		// w_3^2  = (-1/2 + i*sqrt(3)/2) = a-bi
		// y[0] = x[0]+x[1]+x[2]
		// y[1] = x[0]+x[1]*w_3+x[2]*w_3^2
		// y[2] = x[0]+x[1]*w_3^2 + x[2]*w_3
		a  = -0.5; b = -sqrt(3)/2;
		y_r[0] = x_r[0] + x_r[1] + x_r[2];
		y_i[0] = x_i[0] + x_i[1] + x_i[2];
		y_r[1] = x_r[0] + (x_r[1]*a-x_i[1]*b) + (x_r[2]*a+x_i[2]*b); 
		y_i[1] = x_i[0] + (x_i[1]*a+x_r[1]*b) + (x_i[2]*a-x_r[2]*b);
		y_r[2] = x_r[0] + (x_r[1]*a+x_i[1]*b) + (x_r[2]*a-x_i[2]*b);
		y_i[2] = x_i[0] + (x_i[1]*a-x_r[1]*b) + (x_i[2]*a+x_r[2]*b);	
	} else if (N == 5)
	{
		y_r[0] = x_r[0] + x_r[1];
		y_i[0] = x_i[0] + x_i[1];
		y_r[1] = x_r[0] - x_r[1]; 
		y_i[1] = x_i[0] - x_i[1];
		y_r[2] = x_r[0] - x_r[1]; 
		y_i[2] = x_i[0] - x_i[1];	
		y_r[3] = x_r[0] - x_r[1]; 
		y_i[3] = x_i[0] - x_i[1];	
		y_r[4] = x_r[0] - x_r[1]; 
		y_i[4] = x_i[0] - x_i[1];	
	} else 
	{
		if(N % 2 == 0)
		{
			double *z_r, *z_i, *u_r, *u_i;
			z_r = (double *) malloc(N*sizeof(double)); // z_r:0~N/2-1: even, z_r:N/2~N-1: odd
			if(z_r==NULL) { printf("no memory!!\n"); return 0; }
			z_i = (double *) malloc(N*sizeof(double)); // z_i:0~N/2-1: even, z_i:N/2~N-1: odd 
			if(z_i==NULL) { printf("no memory!!\n"); return 0; }
			u_r = (double *) malloc(N*sizeof(double)); // z_r:0~N/2-1: even, z_r:N/2~N-1: odd
			if(u_r==NULL) { printf("no memory!!\n"); return 0; }
			u_i = (double *) malloc(N*sizeof(double)); // z_i:0~N/2-1: even, z_i:N/2~N-1: odd 
			if(u_i==NULL) { printf("no memory!!\n"); return 0; }
			for(k=0;k<N/2;++k)
			{
				z_r[k] = x_r[2*k];
				z_i[k] = x_i[2*k];
				z_r[N/2+k]  = x_r[2*k+1];
				z_i[N/2+k]  = x_i[2*k+1];
			}		
			FFT_general(z_r, z_i, u_r, u_i, N/2);
			FFT_general(z_r+N/2, z_i+N/2, u_r+N/2, u_i+N/2, N/2);
			w_N_re =  cos(2.0*M_PI/N);
			w_N_im = -sin(2.0*M_PI/N);
			w_re   = 1.0;
			w_im   = 0.0; 
			for(k=0;k<N/2;++k)
			{
				a = w_re*u_r[N/2+k] - w_im*u_i[N/2+k];
				b = w_re*u_i[N/2+k] + w_im*u_r[N/2+k];
				y_r[k]     = u_r[k] + a;
				y_i[k]     = u_i[k] + b;
				y_r[N/2+k] = u_r[k] - a;
				y_i[N/2+k] = u_i[k] - b;
				temp = w_re;
				w_re = w_re*w_N_re - w_im*w_N_im;
				w_im = temp*w_N_im + w_im*w_N_re;
			}
			free(u_r); free(u_i); free(z_r); free(z_i);
		} else if(N % 3 == 0)
		{
			double *z_r, *z_i, *u_r, *u_i;
			z_r = (double *) malloc(N*sizeof(double)); // z_r:0~N/2-1: even, z_r:N/2~N-1: odd
			z_i = (double *) malloc(N*sizeof(double)); // z_i:0~N/2-1: even, z_i:N/2~N-1: odd 
			u_r = (double *) malloc(N*sizeof(double)); // z_r:0~N/2-1: even, z_r:N/2~N-1: odd
			u_i = (double *) malloc(N*sizeof(double)); // z_i:0~N/2-1: even, z_i:N/2~N-1: odd 

			for(k=0;k<N/3;++k)
			{
				z_r[k] = x_r[3*k];
				z_i[k] = x_i[3*k];
				z_r[N/3+k]  = x_r[3*k+1];
				z_i[N/3+k]  = x_i[3*k+1];
				z_r[2*N/3+k]  = x_r[3*k+2];
				z_i[2*N/3+k]  = x_i[3*k+2];
			}
			FFT_general(z_r, z_i, u_r, u_i, N/3);
			FFT_general(z_r+N/3, z_i+N/3, u_r+N/3, u_i+N/3, N/3);
			FFT_general(z_r+2*N/3, z_i+2*N/3, u_r+2*N/3, u_i+2*N/3, N/3);
			w_N_re =  cos(2.0*M_PI/N);
			w_N_im = -sin(2.0*M_PI/N);
			w_re   = 1.0;
			w_im   = 0.0; 
			for(k=0;k<N/3;++k)
			{
				u1_r = w_re*u_r[N/3+k] - w_im*u_i[N/3+k];
				u1_i = w_re*u_i[N/3+k] + w_im*u_r[N/3+k];
				// (w_re + i w_im)^2 = w_re*w_re - w_im*w_im + i (2*w_re*w_im)
				u2_r = (w_re*w_re - w_im*w_im)*u_r[2*N/3+k] - (2*w_re*w_im)*u_i[2*N/3+k];
				u2_i = (w_re*w_re - w_im*w_im)*u_i[2*N/3+k] + (2*w_re*w_im)*u_r[2*N/3+k];
				a = -0.5; b = -sqrt(3)/2;
				y_r[k]       = u_r[k] + u1_r + u2_r;
				y_i[k]       = u_i[k] + u1_i + u2_i;
				y_r[N/3+k]   = u_r[k] + (u1_r*a-u1_i*b) + (u2_r*a+u2_i*b);
				y_i[N/3+k]   = u_i[k] + (u1_i*a+u1_r*b) + (u2_i*a-u2_r*b);
				y_r[2*N/3+k] = u_r[k] + (u1_r*a+u1_i*b) + (u2_r*a-u2_i*b);
				y_i[2*N/3+k] = u_i[k] + (u1_i*a-u1_r*b) + (u2_i*a+u2_r*b);
				temp = w_re;
				w_re = w_re*w_N_re - w_im*w_N_im;
				w_im = temp*w_N_im + w_im*w_N_re;
			}
			free(u_r); free(u_i); free(z_r); free(z_i);
		} else if(N % 5 == 0)
		{
			
		}
	}	
	//free(p);
	return 0;
} 

