#include <stdio.h>
#include <stdlib.h> 
#include <omp.h>
#include <math.h>
#include <time.h>
#define DEBUG 0

int main()
{
	int i, p, n = 0;
	int N ;	
	double *x_re, *x_im;
	double T1;
	clock_t t1, t2;
	
	N = 27000000;	
	printf("N = %d\n",N);
	x_re = (double *) malloc( N * sizeof(double));
	x_im = (double *) malloc( N * sizeof(double));
	for(i=0;i<N;++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}
	
	t1 = clock();
	FFT(x_re, x_im, N);
	t2 = clock();
	T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
	printf("FFT_ver2 of %d elements: %f\n",N, T1);
	#if DEBUG 	
	for(i=0;i<N;++i)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	} 
	#endif
	free(x_re);
	free(x_im);
	return;
}

//========================================================================================

int FFT(double *x_re, double *x_im, int N)
{	
	int N0, i, p, n = 0;
	int *order;
	double T1;
	clock_t t1, t2;
	N0 = N;
    p = 1;
    order = (int *) malloc( N * sizeof(int));
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
		#if DEBUG 
		printf("order[%d] = %d, N0 = %d\n",n, p, N0);	
		#endif
		n++;
	}
	n--;
	#if DEBUG 
	printf("n = %d\n",n);
	for(i=0;i<=n;i++)
	{
		printf("order[%d] = %d\n", i, order[i]);
	}
	#endif
	
	t1 = clock();
	bit_reverse(x_re, x_im, N, order);
	t2 = clock();
	T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
	printf("FFT_ver2 of %d elements: %f\n",N, T1);
	
	
	i = 0;
	
	while(N0<N && i<=n)
	{
		#if DEBUG
		printf("N0 = %d, order[%d] = %d\n", N0, i, order[i]);
		#endif
		
		
		t1 = clock();
		butterfly(x_re, x_im, N, order[i], N0);
		t2 = clock();
		T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
		printf("FFT_ver2 of %d elements: %f\n",N, T1);
	
		N0 *= order[i];
		i++;
	}
	free(order);
	return 0;
}

//=============================================================================================================================

	//                order[0]=3(p=3)	 order[1]=2(p=2)   order[2]=2(p=2)
	// N = 12, N0 = 12       ->    N0 = 4       ->    N0 = 2       ->    N0 = 1
	// prime = 2 or 3 or 5

int butterfly(double *x_re, double *x_im, int N, int prime, int N0)
{
	int k, p, q, r, s, t, m;
	double w_re, w_im, w_N_re, w_N_im, w_2_re, w_2_im, w_3_re, w_3_im, w_4_re, w_4_im, theta, theta1; 
	double t_rp, t_rq, t_rr, t_rs, t_rt, t_ip, t_iq, t_ir, t_is, t_it;		
	N0 *= prime;
	
	if(prime == 2)
	{
		for(k=0;k<N0/2;k++) 
		{
			theta = 2.0*k*M_PI/N0;		//找下一個 W_N 要加的角度
			w_re = cos(theta);
			w_im = -sin(theta);
			
			for(p=k;p<N;p+=N0)
			{
				q = p + N0/2;
				
        		#if DEBUG
				printf("(%d,%d) (%f,%f) FFT2 \n", p,q, w_re, w_im);
				#endif
				// multiply (w_re + w_im * i) on x[q]
				t_rq = x_re[q]; 
				x_re[q] = w_re*x_re[q] - w_im*x_im[q];
				x_im[q] = w_re*x_im[q] + w_im*t_rq; 
				
				t_rp = x_re[p];
				x_re[p] = x_re[p] + x_re[q];
				x_re[q] = t_rp    - x_re[q]; 
				t_ip = x_im[p];
				x_im[p] = x_im[p] + x_im[q];
				x_im[q] = t_ip    - x_im[q]; 
			}
		}	
	}
	else if(prime == 3)
	{
		theta1 = 2.0*M_PI/3;				
		w_N_re =  cos(theta1);
		w_N_im = -sin(theta1);
		for(k=0;k<N0/3;k++)
		{	 			
			theta = 2.0*k*M_PI/N0;	 
            w_re =  cos(theta);				
            w_im = -sin(theta);				       	
											  
			for(p=k;p<N;p+=N0)				
			{
				q = p + N0/3;
				r = q + N0/3;
				
				w_2_re = w_re * w_re - w_im * w_im;
				w_2_im = w_re * w_im + w_re * w_im;	
				
				// multiply (w_re + w_im * i) on x[q]
				t_rq = x_re[q];
				x_re[q] = w_re * x_re[q] - w_im * x_im[q];
				x_im[q] = w_re * x_im[q] + w_im * t_rq;
				// multiply (w_re + w_im * i)^2 = w_re - w_im * i on x[r]
				t_rr = x_re[r];
				x_re[r] = w_2_re * x_re[r] - w_2_im*x_im[r];
				x_im[r] = w_2_re * x_im[r] + w_2_im*t_rr;
				
				t_rp = x_re[p];
				t_rq = x_re[q];
				t_rr = x_re[r];		
						
				t_ip = x_im[p];
				t_iq = x_im[q];	
				t_ir = x_im[r];			
				
				x_re[r] = t_rp + w_N_re * (t_rq + t_rr) + w_N_im * (t_iq - t_ir);
				x_re[q] = t_rp + w_N_re * (t_rq + t_rr) - w_N_im * (t_iq - t_ir); 
				x_re[p] = t_rp + t_rq + t_rr;				
				
				x_im[r] = t_ip + w_N_re * (t_iq + t_ir) - w_N_im * (t_rq - t_rr);
				x_im[q] = t_ip + w_N_re * (t_iq + t_ir) + w_N_im * (t_rq - t_rr);
				x_im[p] = t_ip + t_iq + t_ir;	
			}		
		}
	}
	else if(prime == 5)
	{
		theta1 = 2.0*M_PI/5;
		w_N_re =  cos(theta1);
		w_N_im = -sin(theta1);
		for(k=0;k<N0/5;k++)
		{	
			theta = 2.0*k*M_PI/N0;		
            w_re =  cos(theta);				
            w_im = -sin(theta);			
			 
			for(p=k;p<N;p+=N0)		
			{
				q = p + N0/5;
				r = p + 2*N0/5;
				s = p + 3*N0/5;
				t = p + 4*N0/5;
				
				w_2_re = w_re * w_re - w_im * w_im;
				w_2_im = w_re * w_im + w_re * w_im;
				
				w_3_re = w_re * w_re * w_re - 3 * w_re * w_im * w_im;
				w_3_im = 3 * w_re * w_re * w_im - w_im * w_im * w_im;
				
				w_4_re = w_2_re * w_2_re - w_2_im * w_2_im;
				w_4_im = 2 * w_2_re * w_2_im;			
				
				// multiply (w_re + w_im * i) on x[q]
				t_rq = x_re[q];
				x_re[q] = w_re * x_re[q] - w_im * x_im[q];
				x_im[q] = w_re * x_im[q] + w_im * t_rq;
				// multiply (w_re + w_im * i)^2 = (w_re * w_re - w_im * w_im) + (w_re * w_im + w_re * w_im) * i on x[r]
				t_rr = x_re[r];
				x_re[r] = w_2_re * x_re[r] - w_2_im * x_im[r];
				x_im[r] = w_2_re * x_im[r] + w_2_im * t_rr;
				
				t_rs = x_re[s];
				x_re[s] = w_3_re * x_re[s] - w_3_im * x_im[s];
				x_im[s] = w_3_re * x_im[s] + w_3_im * t_rs;
				
				t_rt = x_re[t];
				x_re[t] = w_4_re * x_re[t] - w_4_im * x_im[t];
				x_im[t] = w_4_re * x_im[t] + w_4_im * t_rt;
							
				t_rp = x_re[p];
				t_rq = x_re[q];
				t_rr = x_re[r];
				t_rs = x_re[s];
				t_rt = x_re[t];		
						
				t_ip = x_im[p];
				t_iq = x_im[q];	
				t_ir = x_im[r];	
				t_is = x_im[s];	
				t_it = x_im[t];			
									
				x_re[p] = t_rp + t_rq + t_rr + t_rs + t_rt;				
				x_re[q] = t_rp + w_N_re * (t_rq + t_rt) - w_N_im * (t_iq - t_it) + (w_N_re * w_N_re - w_N_im * w_N_im) * (t_rr + t_rs) - (2 * w_N_re * w_N_im) * (t_ir - t_is);
				x_re[r] = t_rp + w_N_re * (t_rr + t_rs) - w_N_im * (t_is - t_ir) + (w_N_re * w_N_re - w_N_im * w_N_im) * (t_rq + t_rt) - (2 * w_N_re * w_N_im) * (t_iq - t_it);
				x_re[s] = t_rp + w_N_re * (t_rr + t_rs) - w_N_im * (t_ir - t_is) + (w_N_re * w_N_re - w_N_im * w_N_im) * (t_rt + t_rq) - (2 * w_N_re * w_N_im) * (t_it - t_iq);
				x_re[t] = t_rp + w_N_re * (t_rq + t_rt) - w_N_im * (t_it - t_iq) + (w_N_re * w_N_re - w_N_im * w_N_im) * (t_rs + t_rr) - (2 * w_N_re * w_N_im) * (t_is - t_ir);
				
				x_im[p] = t_ip + t_iq + t_ir + t_is + t_it;
				x_im[q] = t_ip + w_N_re * (t_iq + t_it) + w_N_im * (t_rq - t_rt) + (w_N_re * w_N_re - w_N_im * w_N_im) * (t_ir + t_is) + (2 * w_N_re * w_N_im) * (t_rr - t_rs);
				x_im[r] = t_ip + w_N_re * (t_ir + t_is) + w_N_im * (t_rs - t_rr) + (w_N_re * w_N_re - w_N_im * w_N_im) * (t_iq + t_it) + (2 * w_N_re * w_N_im) * (t_rq - t_rt);		
				x_im[s] = t_ip + w_N_re * (t_ir + t_is) + w_N_im * (t_rr - t_rs) + (w_N_re * w_N_re - w_N_im * w_N_im) * (t_it + t_iq) + (2 * w_N_re * w_N_im) * (t_rt - t_rq);
				x_im[t] = t_ip + w_N_re * (t_iq + t_it) + w_N_im * (t_rt - t_rq) + (w_N_re * w_N_re - w_N_im * w_N_im) * (t_is + t_ir) + (2 * w_N_re * w_N_im) * (t_rs - t_rr);		
			}
		}
	}
	else
	{
		printf("we didn't do this! p = 1; N = %d\n", N);
	}
	return;
}

//========================================================================================

int bit_reverse(double *x_re, double *x_im, int N, int *order)
{
    int m, t, p, q, i, k, change_index_p;
    int n = 0;
    int *change_index;
    int *cyclic;
    m = N/order[n];
    
    change_index = (int*)malloc(2*N*sizeof(int));
    cyclic = (int*)malloc(N*sizeof(int));
    
    #if DEBUG
	printf("order[%d] = %d , m = %d\n",n, order[n], m);
	#endif
	
	//第一組 p, q 不動 
	change_index[0] = 0;
	change_index[1] = 0;
	//最後一組 p, q 也不動 
	change_index[2*N-2] = N-1;
	change_index[2*N-1] = N-1;
	
	q = m;
	i = 2;
	//中間的 p, q 需用 bit_reverse 查看	 
	for(p=1;p<N-1;p++)
   	{	
   		change_index[i] = p;
	   	change_index[i+1] = q;
	   	
		
	   	#if DEBUG 
	    printf("i_p=%d, p=%d, i_q=%d,  q=%d\n", i, change_index[i], i+1, change_index[i+1]);
	    printf("p=%d, q=%d\n", change_index[10], change_index[11]);
	    #endif
	    i += 2;	  	
	    
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
	
	change_index_p = 0;
	while(cyclic[0] != N-1)
	{
		change_index_p += 2;
		p = change_index[change_index_p];	// p 的起始值 (在 change_index[2] 的位置) 
		while(p == 0)
		{
			change_index_p += 2;
			p = change_index[change_index_p];
		}
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
		
		#if DEBUG
		for(i=0;i<n;i++)
		{
			printf("cyclic[%d] = %d\n", i, cyclic[i]);
		}
		#endif
		
		// x[p] 與 x[q] 互換，終止條件：下一個要換的，為一開始的 p, i.e. q = p (= cyclic[0])
		// Ex. 1 -> 4 -> 6 ->1   	
		// cyclic = [1 4 6]
		// i=0, 1 <-> 4, i=1, 4 <-> 6
		for(i=0;i<n-1;i++)
		{	
			t = x_re[cyclic[i]];
		    x_re[cyclic[i]] = x_re[cyclic[i+1]];
			x_re[cyclic[i+1]] = t;
		    t = x_im[cyclic[i]];
		    x_im[cyclic[i]] = x_im[cyclic[i+1]];
			x_im[cyclic[i+1]] = t;
			
			change_index[2*cyclic[i]] = 0;			//將交換過的 index 歸零 
			change_index[2*cyclic[i]+1] = 0;
			change_index[2*cyclic[i+1]] = 0;
			change_index[2*cyclic[i+1]+1] = 0;
		}
				
		#if DEBUG
		for(i=0;i<2*N;i++)
		{
			printf("change_index[%d]=%d\n", i, change_index[i]);
		}
		#endif
	}
	free(change_index);
	free(cyclic);
	return 0;
}






