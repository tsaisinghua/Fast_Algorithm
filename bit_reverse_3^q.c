#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#define DEBUG 0

int main()
{
	int i;
	double y_re[27], y_im[27], x_re[27], x_im[27];
	for(i=0;i<27;++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}
	bit_reverse(x_re, x_im, 27);	
	butterfly(x_re, x_im, 27);	
	for(i=0;i<27;++i)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	}
	return;
	 
}
int bit_reverse(double *x_re, double *x_im, int N)
{
    int m,p,q,k;
    double t;
    
    m = N/3;                      	// Bit-Reverse 每次要進位的數字 
    q = m;							// p = 1, q = m (第一個要交換的) 
    for(p=1;p<N-1;++p)
    {
    	printf("p=%d <-> q=%d\n",p,q);
    	#if DEBUG
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
        
        k = m;						// k, 用來檢查第 log_3 k + 1 位是不是1 
        #if DEBUG
		printf("k=%d\n",k);
		#endif
		
		while(q >= 2*k & k > 0)		// q >=k 第 (log_2 k + 1)位是1,  
        {
        	#if DEBUG
        	printf("1.q=%d, k=%d\n",q,k);
        	#endif
        	
            q = q-2*k;				// 1->0
            k = k/3;				// 檢查下一位 
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
	int k, p, q, r, m;
	double theta, theta1;
	double w_re, w_im, w_N_re, w_N_im;
	double t, t_rp, t_rq, t_rr;
	double t_ip, t_iq, t_ir;
	
	//m = 1;
	for(m=1;m<N;m*=3)
	{
		//w_re = 1.0;
		//w_im = 0.0;		
		theta1 = 2.0*M_PI/3;				// 找下一個 W_N 要加的角度 (W_3)^0 = 1 + i*0 
											// => (W_3)^1 = cos(2*M_PI/3.0) -i*sin(2*M_PI/3.0) 
											// => (W_3)^2 = cos(4*M_PI/3.0) -i*sin(4*M_PI/3.0) = cos(2*M_PI/3.0) + i*sin(2*M_PI/3.0)
											// => (W_3)^3 = 1 +i*0
		//printf("%f\n", theta);
		w_N_re =  cos(theta1);
		w_N_im = -sin(theta1);
		for(k=0;k<m;k++)
		{				
			theta = (2.0/3)*k*M_PI/m;		// 設定 w_re 及 w_im 的起始值 (m = 1, k = 0) => theta = 0; w_re = 1.0, w_im = 0.0 
            w_re =  cos(theta);				//							  (m = 3, k = 0) => theta = 0; w_re = 1.0, w_im = 0.0	
            w_im = -sin(theta);				//							  (m = 3, k = 1) => theta = (2/9)*pi; w_re = cos((2/9)*pi), w_im = -sin((2/9)*pi)          
			for(p=k;p<N;p+=3*m)				//							  (m = 3, k = 2) => theta = (4/9)*pi; w_re = cos((4/9)*pi), w_im = -sin((4/9)*pi)  
			{
				q = p + m;
				r = q + m;
				printf("(%d,%d,%d) (%f,%f) FFT2 \n", p, q, r, w_re, w_im);
				 
				// multiply (w_re + w_im * i) on x[q]
				t_rq = x_re[q];
				x_re[q] = w_re * x_re[q] - w_im * x_im[q];
				x_im[q] = w_re * x_im[q] + w_im * t_rq;
				// multiply (w_re + w_im * i)^2 = w_re - w_im * i on x[r]
				t_rr = x_re[r];
				x_re[r] = (w_re * w_re - w_im * w_im) * x_re[r] - (w_re * w_im + w_re * w_im)*x_im[r];
				x_im[r] = (w_re * w_re - w_im * w_im) * x_im[r] + (w_re * w_im + w_re * w_im)*t_rr;
				
				//FFT 3 points
				// (w_re + w_im * i)^2 = w_re - w_im * i
				// (w_re + w_im * i)^4 = w_re + w_im * i
				// x[p] = x[p] + x[q]+ x[r] => x_im[p] = (x_re[p] + x_re[q]+ x_re[r]) +(x_im[p] + x_im[q]+ x_im[r])*i
				// x[q] = x[p] + w*x[q]+ w*w*x[r] 
				// x_re[q] = x_re[p] + (w_re * x_re[q] - w_im * x_im[q]) + (w_re * x_re[r] + w_im * x_im[r])
				//		   = x_re[p] + w_re * (x_re[q] + x_re[r]) - w_im * (x_im[q] - x_im[r])
				// x_im[q] = x_im[p] + (w_re * x_im[q] + w_im * x_re[q]) + (w_re * x_im[r] - w_im * x_re[r])
				// 		   = x_im[p] + w_re * (x_im[q] +  x_im[r]) + w_im * (x_re[q]- x_re[r])
				
				// x[r] = x[p] + w*w*x[q]+ w*w*w*w*x[r]
				// x_re[r] = x_re[p] + (w_re * x_re[q] + w_im * x_im[q]) + (w_re * x_re[r] - w_im * x_im[r])
				//         = x_re[p] + w_re * (x_re[q] + x_re[r]) + w_im * (x_im[q] - x_im[r])
				// x_im[r] = x_im[p] + (w_re * x_im[q] - w_im * x_re[q]) + (w_re * x_im[r] + w_im * x_re[r])
				//		   = x_im[p] + w_re * (x_im[q] + x_im[r]) - w_im * (x_re[q] - x_re[r])
				
				t_rp = x_re[p];
				t_rq = x_re[q];
				t_rr = x_re[r];				
				t_ip = x_im[p];
				t_iq = x_im[q];	
				t_ir = x_im[r];			
				
				/*
				x_re[p] = x_re[p] + x_re[q] + x_re[r];
				x_re[q] = t_rp    + w_re * (x_re[q] + x_re[r]) - w_im * (x_im[q] - x_im[r]); 
				x_re[r] = t_rp    + w_re * (t_rq    + x_re[r]) + w_im * (x_im[q] - x_im[r]);
					
				x_im[p] = x_im[p] + x_im[q] + x_im[r];
				x_im[q] = t_ip    + w_re * (x_im[q] +  x_im[r]) + w_im * (t_rq - t_rr);
				x_im[r] = t_ip    + w_re * ( t_iq   +  x_im[r]) - w_im * (t_rq - t_rr);
				printf("(%f+%fi,%f+%fi,%f+%fi) \n", x_re[p], x_im[p], x_re[q], x_im[q], x_re[r], x_im[r]);				
				*/
				x_re[r] = t_rp + w_N_re * (t_rq + t_rr) + w_N_im * (t_iq - t_ir);
				x_re[q] = t_rp + w_N_re * (t_rq + t_rr) - w_N_im * (t_iq - t_ir); 
				x_re[p] = t_rp + t_rq + t_rr;				
				
				x_im[r] = t_ip + w_N_re * (t_iq + t_ir) - w_N_im * (t_rq - t_rr);
				x_im[q] = t_ip + w_N_re * (t_iq + t_ir) + w_N_im * (t_rq - t_rr);
				x_im[p] = t_ip + t_iq + t_ir;	
				
			}		
			/*
			t    = w_re;
			w_re = w_N_re * w_re - w_N_im * w_im;
			w_im = w_N_re * w_im + w_N_im *t;
			*/
		}	
		//m = m * 3;
	}
	return;
}



