#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#define DEBUG 0

int main()
{
	int i, p, n = 0;
	int N = 12;	
	double x_re[N], x_im[N];
		
	for(i=0;i<N;++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}
	//bit_reverse(x_re, x_im, N);	
	//butterfly(x_re, x_im, N);		
	FFT(x_re, x_im, N);
	//bit_reverse(x_re, x_im, N, order[n]);
	for(i=0;i<N;++i)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	}
 
	return;
}

int FFT(double *x_re, double *x_im, int N)
{	
	int N0, i, p, n = 0;
	int order[N];
	N0 = N;
    p = 1;
    while (N0>1)
	{
        if ((N0%5)==0) 
		{
            p=5;            
        }
		else if ((N0%3)==0)
		{
            p=3;            
		}
		else if ((N0%2)==0)
		{
			p=2;			
		}
		else
		{
            p=1;
        }
		order[n] = p;
		N0 /= p;
		printf("order[%d] = %d, N0 = %d\n",n, p, N0);		
		n++;
	}
	n--;
	#if DEBUG 
	printf("n = %d\n",n);
	for(i=0;i<=n;i++)
	{
		printf("order[%d] = %d\n",i, order[i]);
	}
	#endif
	bit_reverse(x_re, x_im, N, order);
}

//========================================================================
int bit_reverse(double *x_re, double *x_im, int N, int *order)
{
    int m, t, p, q, k;
    int n = 0;
    m = N/order[n];
    q = m;
    
	printf("order[%d] = %d , m = %d\n",n, order[n], m);
	
	for(p=1;p<N-1;++p)
    {
    	
        printf("p=%d <-> q=%d\n",p,q);
    	#if DEBUG
        printf("m=%d, q=%d, %d <-> %d\n",m,q,p,q);
        #endif  
		/*  
        if(p < q)
        {
            t = x_re[p];
            x_re[p] = x_re[q];
			x_re[q] = t;
            t = x_im[p];
            x_im[p] = x_im[q];
			x_im[q] = t;
        }
        */
        k = m;
        n = 0;
        while(q >= (order[n]-1)*k & k > 0)
        {
            q = q-(order[n]-1)*k;
			n++;
            k /= order[n];
            //printf("%d,%d\n",q,k);
        }
        q = q+k;
		#if DEBUG
   		printf("2.q=%d, k=%d\n",q,k);
  		printf("========================\n");
		#endif  
	}
	
	
	return 0;
}
