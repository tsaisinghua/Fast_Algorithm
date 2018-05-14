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
	*/
	p = 3;
	q = 0;
	r = 0;
	N = pow(2, p) * pow(3, q) * pow(5, r);
	printf("N=%d\n",N);
	
	double y_re[N], y_im[N], x_re[N], x_im[N];
	for(i=0;i<N;++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}
	
	t1 = clock();
	bit_reverse(x_re, x_im, N);
	butterfly_2(x_re, x_im, N);
	t2 = clock();
	T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
	printf("FFT_ver2 of %d elements: %f\n",N, T1);
	
	//#if DEBUG
	for(i=0;i<N;++i)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	}
	//#endif
	return;
	 
}
int bit_reverse(double *x_re, double *x_im, int N)
{
    int m,p,q,k;
    double t;
    
    m = N/2;                      	// Bit-Reverse 每次要進位的數字 
    q = m;							// p = 1, q = m (第一個要交換的) 
    for(p=1;p<N-1;++p)
    {
		#if DEBUG
    	printf("p=%d <-> q=%d\n",p,q);
        printf("m=%d, q=%d, %d <-> %d\n",m,q,p,q);
        #endif
        
        if(p < q)
        {
            t = x_re[p];
            x_re[p] = x_re[q];
			x_re[q] = t;
            t = x_im[p];
            x_im[p] = x_im[q];
			x_im[q] = t;			 
        }
        k = m;						// k, 用來檢查第 log_2 k + 1 位是不是1 
        							// (1110) = 14 >= (8=k)
									// (0110) = 6  >= (4=k)
									// (0010) = 2  >= (2=k)
									// (0000) = 0  >= (1=k) X break --> (0001)
									// (0110) = 6  >= (8=k) X break --> (1110) 
        #if DEBUG
		printf("k=%d\n",k);
		#endif
		
		while(q >= k & k > 0)		// q >=k 第 (log_2 k + 1)位是1,  
        {
        	#if DEBUG
        	printf("1.q=%d, k=%d\n",q,k);
        	#endif
        	
            q = q-k;				// 1->0
            k = k/2;				// 檢查下一位 
        }
        q = q+k;
        
        #if DEBUG
        printf("2.q=%d, k=%d\n",q,k);
  		printf("========================\n");
		#endif      
    }
    return 0;
}

int butterfly(double *x_re, double *x_im, int N)
{
	int k, p, q, m;
	double w_re, w_im, w_N_re, w_N_im, t, theta; 
	//m = 1;

	for(m=1;m<N;m*=2)
	{
		w_re = 1.0;
		w_im = 0.0; 
		theta = M_PI/m;		//找下一個 W_N 要加的角度	
		for(k=0;k<m;++k) 
		{
			for(p=k;p<N;p+=2*m)
			{
				q = p + m;
				
        		#if DEBUG
				printf("(%d,%d) (%f,%f) FFT2 \n", p,q, w_re, w_im);
				#endif
				// multiply (w_re + w_im * i) on x[q]
				t = x_re[q]; 
				x_re[q] = w_re*x_re[q] - w_im*x_im[q];
				x_im[q] = w_re*x_im[q] + w_im*t; 
				
				t = x_re[p];
				x_re[p] = x_re[p] + x_re[q];
				x_re[q] = t       - x_re[q]; 
				t = x_im[p];
				x_im[p] = x_im[p] + x_im[q];
				x_im[q] = t       - x_im[q]; 
			}
			t    = w_re; 
			w_re = w_N_re*w_re - w_N_im*w_im;
			w_im = w_N_re*w_im + w_N_im*t;
		}
		//m = m * 2;
	}
	
	return;
}

//===========================================================================================

int butterfly_2(double *x_re, double *x_im, int N)
{
	int k, p, q, m;
	double w_re, w_im, w_N_re, w_N_im, t, theta; 
	//m = 1;
	for(m=1;m<N;m*=2)
	{
		for(k=0;k<m;++k) 
		{
			theta = k*M_PI/m;		//找下一個 W_N 要加的角度
			w_re = cos(theta);
			w_im = -sin(theta);
			
			for(p=k;p<N;p+=2*m)
			{
				q = p + m;
				
        		#if DEBUG
				printf("(%d,%d) (%f,%f) FFT2 \n", p,q, w_re, w_im);
				#endif
				// multiply (w_re + w_im * i) on x[q]
				t = x_re[q]; 
				x_re[q] = w_re*x_re[q] - w_im*x_im[q];
				x_im[q] = w_re*x_im[q] + w_im*t; 
				
				t = x_re[p];
				x_re[p] = x_re[p] + x_re[q];
				x_re[q] = t       - x_re[q]; 
				t = x_im[p];
				x_im[p] = x_im[p] + x_im[q];
				x_im[q] = t       - x_im[q]; 
			}
		/*
			t    = w_re; 
			w_re = w_N_re*w_re - w_N_im*w_im;
			w_im = w_N_re*w_im + w_N_im*t;
		*/
		}
		//m = m * 2;
	}	
	return;
}

//===========================================================================================

int Fast_Fourier_Transform(double *y_re, double *y_im, double *x_re, double *x_im, int N)
{
	if(N==2) 
	{
		// y, y[0] = x[0] + x[1], y[1] = x[0] - x[1]
		y_re[0] = x_re[0] + x_re[1];
		y_im[0] = x_im[0] + x_im[1];
		y_re[1] = x_re[0] - x_re[1]; 
		y_im[1] = x_im[0] - x_im[1];
	} 
	else 
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
}
