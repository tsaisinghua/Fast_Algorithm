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
    int m, t, p, q, i, k;
    int n = 0;
    int *change_index;
    int *cyclic;
    m = N/order[n];
    
    change_index = (int*)malloc(2*N*sizeof(int));
    cyclic = (int*)malloc(N*sizeof(int));
    
	printf("order[%d] = %d , m = %d\n",n, order[n], m);
	
	q = m;
	i = 2;
	//第一組 p, q 不動 
		change_index[0] = 0;
		change_index[1] = 0;
		//最後一組 p, q 也不動 
		change_index[2*N-2] = N-1;
		change_index[2*N-1] = N-1;
	for(p=1;p<N-1;p++)
   	{
   		
		//中間的 p, q 需用 bit_reverse 查看	 
   		change_index[i] = p;
	   	change_index[i+1] = q;
	   	
		
	   	#if DEBUG 
	    printf("i_p=%d, p=%d, i_q=%d,  q=%d\n", i, change_index[i], i+1, change_index[i+1]);
	    printf("p=%d, q=%d\n", change_index[10], change_index[11]);
	    #endif
	    i += 2;	  	
	    		
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
	
	p = 1;
	n = 0;
	cyclic[0] = p;
	//printf("p=%d, q=%d\n", change_index[10], change_index[11]);
	while((change_index[2*p] != change_index[2*p+1]) && (change_index[2*p+1] != cyclic[0]))
	{
		cyclic[n] = change_index[2*p];
		cyclic[n+1] = change_index[2*p+1];
		p = change_index[2*p+1];
		n++;
	}
	n++;
	for(i=0;i<n;i++)
	{
		printf("cyclic[%d] = %d\n", i, cyclic[i]);	
	}

	
	free(change_index);
	free(cyclic);
	return 0;
}
