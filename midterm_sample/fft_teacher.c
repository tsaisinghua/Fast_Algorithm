#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
int FFT(double *x_r, double *x_i, double *y_r, double *y_i, int N);
int Generate_N(int p, int q, int r);
int SFT(double *x_r, double *x_i, double *y_r, double *y_i, int N);
int Initial(double *x, double *y, int N);
int Print_Complex_Vector(double *x, double *y, int N);
int Bit_Increase(int *D, int b, int N);
int Bit_Reverse(int *D, int b, int N);
int Bit_Reverse_Integer(int N);
int FFTv3(double *x_r, double *x_i, double *y_r, double *y_i, int N);
int Group(int N);
int main()
{
	int D[4];
	//Bit_Reverse(D, 3, 4);
	//Bit_Reverse_Integer(27);
	Group(8);
	// y_k = sum(x_n * w^{-kn}, n=0..N-1)
	// w = cos(2*pi/N) + i sin(2*pi/N)
	int k, n, N, p, q, r;
	double *y_r, *y_i, *x_r, *x_i, w_r, w_i;
	clock_t t1, t2;
	// Input N
	printf("Please input p q r=");
	scanf("%d %d %d", &p, &q, &r);
	N = Generate_N(p, q, r);
	printf("N=2^%d 3^%d 5^%d = %d\n",p,q,r,N);

	// Create Memory for x and y
	x_r = (double *) malloc(N*sizeof(double));
	x_i = (double *) malloc(N*sizeof(double));
	y_r = (double *) malloc(N*sizeof(double));
	y_i = (double *) malloc(N*sizeof(double));
	
	// y_k = sum(x_n * w^{-kn}, n=0..N-1)
	// w = cos(2*pi/N) + i sin(2*pi/N)
	Initial(x_r, x_i, N);
	t1 = clock();
	SFT(x_r, x_i, y_r, y_i, N);
	t2 = clock();
	printf("%f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
	Print_Complex_Vector(y_r, y_i, N);
	system("pause");
	
	Initial(x_r, x_i, N);
	t1 = clock();
	FFTv3(x_r, x_i, y_r, y_i, N);
	t2 = clock();
	printf("%f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);	
	Print_Complex_Vector(y_r, y_i, N);	
	
	return 0;
	
}
int SFT(double *x_r, double *x_i, double *y_r, double *y_i, int N)
{
	int k, n; 
	double w_r, w_i;
	for(k=0;k<N;++k)
	{
		y_r[k] = 0.0;
		y_i[k] = 0.0;
		for(n=0;n<N;n++)
		{
			// w^{-kn} = cos(-kn*2*pi/N)+ i sin(-kn*2*pi/N)
			w_r = cos(-k*n*2*M_PI/N);
			w_i = sin(-k*n*2*M_PI/N);
			y_r[k] = y_r[k] + x_r[n]*w_r - x_i[n]*w_i; 
			y_i[k] = y_i[k] + x_r[n]*w_i + x_i[n]*w_r;
		}
	}
	return 0;
}
int FFT(double *x_r, double *x_i, double *y_r, double *y_i, int N)
{
	if(N==1)
	{
		y_r[0] = x_r[0];
		y_i[0] = x_i[0];
		return 0; 
	}
	// input : x = x_r + i * x_i
	// output: y = y_r + i * y_i
	int k, n;
	//double *u_r, *u_i, *v_r, *v_i, w_r, w_i;
	double *even_r, *even_i, *odd_r, *odd_i, w_r, w_i;
	double *even_FT_r, *even_FT_i, *odd_FT_r, *odd_FT_i;
	
	even_r = (double *) malloc(N*sizeof(double));
	even_i = (double *) malloc(N*sizeof(double));
	odd_r  = even_r + N/2;//(double *) malloc(N/2*sizeof(double));
	odd_i  = even_i + N/2;//(double *) malloc(N/2*sizeof(double));
	even_FT_r = (double *) malloc(N*sizeof(double));
	even_FT_i = (double *) malloc(N*sizeof(double));
	odd_FT_r = even_FT_r+N/2;//(double *) malloc(N/2*sizeof(double));
	odd_FT_i = even_FT_i+N/2;//(double *) malloc(N/2*sizeof(double));
	for(n=0;n<N/2;n++)
	{
		even_r[n] = x_r[2*n];
		even_i[n] = x_i[2*n];
		odd_r[n] = x_r[2*n+1];
		odd_i[n] = x_i[2*n+1];
	}
	FFT(even_r,even_i,even_FT_r,even_FT_i,N/2);
	FFT(odd_r,odd_i,odd_FT_r,odd_FT_i,N/2);
	for(k=0;k<N/2;++k)
	{
		// w^{-k} 
		w_r = cos(-k*2*M_PI/N);
		w_i = sin(-k*2*M_PI/N);
		printf("N=%d, w_r = %f, w_i = %f\n", N, w_r, w_i);
		y_r[k] = even_FT_r[k] + (w_r*odd_FT_r[k] - w_i*odd_FT_i[k]);
		y_i[k] = even_FT_i[k] + (w_r*odd_FT_i[k] + w_i*odd_FT_r[k]);
		y_r[k+N/2] = even_FT_r[k] - (w_r*odd_FT_r[k] - w_i*odd_FT_i[k]);
		y_i[k+N/2] = even_FT_i[k] - (w_r*odd_FT_i[k] + w_i*odd_FT_r[k]);
	}
	
	free(even_r);
	free(even_i);
	free(even_FT_r);
	free(even_FT_i);
	/*
	for(k=N/2;k<N-1;++k)
	{
		// w^{-k} 
		w_r = cos(-k*2*M_PI/N);
		w_i = sin(-k*2*M_PI/N);
		y_r[k] = even_r[k-N/2] + w_r*odd_r[k-N/2] - w_i*odd_i[k-N/2];
		y_i[k] = even_i[k-N/2] + w_r*odd_i[k-N/2] + w_i*odd_r[k-N/2];
	}
	*/
	/* 
	u_r = (double *) malloc(N*sizeof(double));
	u_i = (double *) malloc(N*sizeof(double));
	v_r = (double *) malloc(N*sizeof(double));
	v_i = (double *) malloc(N*sizeof(double));
	for(n=0;n<N/2;n++)
	{
		u_r[n] = x_r[2*n];
		u_i[n] = x_i[2*n];
		u_r[n+N/2] = x_r[2*n+1];
		u_i[n+N/2] = x_i[2*n+1];
	}
	FFT(u_r, u_i, v_r, v_i, N/2);
	FFT(u_r+N/2, u_i+N/2, v_r+N/2, v_i+N/2, N/2);
	*/
	return 0;
}
int Generate_N(int p, int q, int r)
{
	int i, N = 1;
	// N = 2^p 3^q 5^r
	for(i=0;i<p;++i) N = N * 2;
	for(i=0;i<q;++i) N = N * 3;
	for(i=0;i<r;++i) N = N * 5;	
	
	return N;
}
int Initial(double *x, double *y, int N)
{
	int n;
	for(n=0;n<N;++n)
	{
		x[n] = n; 
		y[n] = 0; 
	}

}
int Print_Complex_Vector(double *x, double *y, int N)
{
	int n;
	for(n=0;n<N;++n)
	{
		if(y[n]>=0) printf("%d : %f +%f i\n", n, x[n], y[n]);
		else printf("%d : %f %f i\n", n, x[n], y[n]);
	}

}
int Bit_Increase(int *D, int b, int N)
{
	//       D[3] D[2] D[1] D[0]
	// print 0000 0001 0010 0011 0100 0101 ...
	int i;
	for(i=0;i<N;++i) D[i] = 0;
	
	while(D[N-1] < b)
	{
		for(i=N-1;i>=0;i--) printf("%d", D[i]);
		printf("\n");
		D[0] = D[0] + 1;
		// check every bit , if D[i] = b, D[i] = 0 D[i+1] + 1
		i = 0;
		while(D[i]==b & i < N-1)
		{
			D[i] = 0;
			i = i + 1;
			D[i] = D[i] + 1;
		}
		system("pause");
	}
	return 0;
}
int Bit_Reverse(int *D, int b, int N)
{
	//       D[3] D[2] D[1] D[0]
	// print 0000 1000 0100 1100 0010 1010 ...
	int i;
	for(i=0;i<N;++i) D[i] = 0;
	
	while(D[0] < b)
	{
		for(i=N-1;i>=0;i--) printf("%d", D[i]);
		printf("\n");
		
		// check every bit , if D[i] = b, D[i] = 0 D[i+1] + 1
		i = N-1;
		while(D[i]==b-1 & i > 0)
		{
			D[i] = 0;
			i = i - 1;
			//D[i] = D[i] + 1;
		}
		D[i] = D[i] + 1;
		system("pause");
	}
	return 0;
}
int Bit_Reverse_Integer(int N)
{
	//		 N = 8
	//       D[2] D[1] D[0]
	// X print 000 100 010 110 001 101 ...
	// print    0   4   2   6   1   5
	// N = 27
	//          0 9 18 3 12 ...
	int i, j, M;
	i = 0; 
	j = 0;
	while(i < N)
	{
		printf("%d <-> %d", i, j);
		
		// check every bit , if D[i] = b, D[i] = 0 D[i+1] + 1
		// i = N-1;
		M = N/3;
		// while(D[i]==b-1 & i > 0)
		while(j >= 2*M & M > 0)
		{
			// D[i] = 0;
			j = j - 2*M;
			// i = i - 1;
			M = M / 3;
		}
		//D[i] = D[i] + 1;
		j = j + M;
		
		i = i + 1;
		system("pause");
	}
	return 0;
}
int FFTv3(double *x_r, double *x_i, double *y_r, double *y_i, int N)
{
	// input : x = x_r + i * x_i
	// output: y = y_r + i * y_i
	int k,n;
	for(n=0;n<N;++n)
	{
		y_r[n] = x_r[n];
		y_i[n] = x_i[n];
	}
	int i, j, M;
	double t_r, t_i;
	i = j = 0;
	while(i < N)
	{
		if(i < j)
		{
			// swap y[i], y[j]
			t_r = y_r[i];
			t_i = y_i[i];
			y_r[i] = y_r[j];
			y_i[i] = y_i[j];
			y_r[j] = t_r;
			y_i[j] = t_i;
		}
		M = N/2;
		while(j >= M & M > 0)
		{
			j = j - M;
			M = M / 2;
		}
		j = j + M;		
		i = i + 1;
	}
	n = 1;
	double u_r, u_i, w_r, w_i;
	while(n < N)
	{
		u_r = cos(-2.0*M_PI/(2*n));
		u_i = sin(-2.0*M_PI/(2*n));
		w_r = 1;
		w_i = 0;
		for(i=0;i<n;++i)
		{
			// 
			w_r = w_r * u_r - w_i * u_i;
			w_i = w_r * u_i + w_i * u_r;
			for(j=i;j<N;j=j+2*n)
			{
				k = j + n;
				t_r = w_r*y_r[k]-w_i*y_i[k];
				t_i = w_r*y_i[k]+w_i*y_r[k];
				// jth jth + w^(-i)*kth = jth + t
				// k   jth - w^(-i)*Kth = jth - t
				y_r[k] = y_r[k] - t_r;
				y_i[k] = y_i[k] - t_i;
				y_r[j] = y_r[j] + t_r;
				y_i[j] = y_i[j] + t_i;
			}
			
		}
		n = n * 2;
	}	
	return 0;
}
int Group(int N)
{
	int i, j, n;
	// N = 8,
	// ((0,1) (2,3)     (4,5) (6,7))  Big Group number = 1, small 4 
	// ((0,2) (4,6))   ((1,3) (5,7))  Big group number = 2, small 2
	// ((0,4)) ((1,5)) ((2,6)) ((3,7))      Big group number = 4, small 1
	n = 1;
	while(n < N)
	{
		
		printf("n=%d\n",n);
		for(i=0;i<n;++i)
		{
			printf("%d:",i);
			// 0
			// 0 1
			// 0 1 2 3
			for(j=i;j<N;j=j+2*n)
			{
				// (0,  ) (2, ) (4, ) (6, )
				// (0,    (4
				printf("(%d,%d)\n",j,j+n);
			}
			
		}
		n = n * 2;
	}
}
