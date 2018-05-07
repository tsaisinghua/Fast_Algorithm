#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#define DEBUG 0

int main()
{
	int i;
	double y_re[25], y_im[25], x_re[25], x_im[25];
	for(i=0;i<25;++i)
	{
		x_re[i] = i;
		x_im[i] = 0.0;
	}
	bit_reverse(x_re, x_im, 25);	
	butterfly(x_re, x_im, 25);	
	for(i=0;i<25;++i)
	{
		printf("%f + %f i\n", x_re[i], x_im[i]);
	}
	return;	 
}
int bit_reverse(double *x_re, double *x_im, int N)
{
    int m,p,q,k;
    double t;
    
    m = N/5;                      	// Bit-Reverse 每次要進位的數字 
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
        
        k = m;						// k, 用來檢查第 log_5 k + 1 位是不是1 
        #if DEBUG
		printf("k=%d\n",k);
		#endif
		
		while(q >= 4*k & k > 0)		// q >=k 第 (log_5 k + 1)位是1,  
        {
        	#if DEBUG
        	printf("1.q=%d, k=%d\n",q,k);
        	#endif
        	
            q = q-4*k;				// 1->0
            k = k/5;				// 檢查下一位 
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
	int k, p, q, r, s, t, m;
	double theta, theta1;
	double w_re, w_im, w_N_re, w_N_im;
	double temp, t_rp, t_rq, t_rr, t_rs, t_rt;
	double 	  t_ip, t_iq, t_ir, t_is, t_it;
	
	//m = 1;
	for(m=1;m<N;m*=5)
	{
		//w_re = 1.0;
		//w_im = 0.0;		
		theta1 = 2.0*M_PI/5;
		w_N_re =  cos(theta1);
		w_N_im = -sin(theta1);
		for(k=0;k<m;k++)
		{	
			theta = (2.0/5)*k*M_PI/m;		
            w_re =  cos(theta);				
            w_im = -sin(theta);			
			 
			for(p=k;p<N;p+=5*m)		
			{
				q = p + m;
				r = q + m;
				s = r + m;
				t = s + m;
				printf("(%d,%d,%d,%d,%d) (%f,%f) FFT2 \n", p, q, r, s, t, w_re, w_im);
				 
				// multiply (w_re + w_im * i) on x[q]
				t_rq = x_re[q];
				x_re[q] = w_re * x_re[q] - w_im * x_im[q];
				x_im[q] = w_re * x_im[q] + w_im * t_rq;
				// multiply (w_re + w_im * i)^2 = (w_re * w_re - w_im * w_im) + (w_re * w_im + w_re * w_im) * i on x[r]
				t_rr = x_re[r];
				x_re[r] = (w_re * w_re - w_im * w_im) * x_re[r] - (w_re * w_im + w_re * w_im)*x_im[r];
				x_im[r] = (w_re * w_re - w_im * w_im) * x_im[r] + (w_re * w_im + w_re * w_im)*t_rr;
				// multiply (w_re + w_im * i)^3 = [(w_re * w_re - w_im * w_im) + (w_re * w_im + w_re * w_im)*i]*(w_re + w_im*i) on x[s]
				//								= ((w_re * w_re - w_im * w_im) * w_re - (w_re * w_im + w_re * w_im) * w_im)
				//								  + ((w_re * w_re - w_im * w_im) * w_im + (w_re * w_im + w_re * w_im) * w_re)*i
				//								= (w_re * w_re * w_re - w_im * w_im * w_re - w_re * w_im * w_im - w_re * w_im * w_im)
				//								  + (w_re * w_re * w_im - w_im * w_im * w_im + w_re * w_im * w_re + w_re * w_im * w_re)*i
				//								= (w_re * w_re * w_re - 3 * w_re * w_im * w_im) + (3 * w_re * w_re * w_im - w_im * w_im * w_im)*i
				t_rs = x_re[s];
				x_re[s] = (w_re * w_re * w_re - 3 * w_re * w_im * w_im) * x_re[s] - (3 * w_re * w_re * w_im - w_im * w_im * w_im) * x_im[s];
				x_im[s] = (w_re * w_re * w_re - 3 * w_re * w_im * w_im) * x_im[s] + (3 * w_re * w_re * w_im - w_im * w_im * w_im) * t_rs;
				// multiply (w_re + w_im * i)^4 = (w_re + w_im * i)^2 * (w_re + w_im * i)^2 on x[t]
				//								= [(w_re * w_re - w_im * w_im) + (2 * w_re * w_im)*i]*[(w_re * w_re - w_im * w_im) + (2 * w_re * w_im)*i]
				//								= [(w_re * w_re - w_im * w_im)*(w_re * w_re - w_im * w_im) - (2 * w_re * w_im)*(2* w_re * w_im)]
				//								  + 2*(w_re * w_re - w_im * w_im)*(2 * w_re * w_im)*i
				t_rt = x_re[t];
				x_re[t] = ((w_re * w_re - w_im * w_im)*(w_re * w_re - w_im * w_im) - (2 * w_re * w_im)*(2* w_re * w_im))*x_re[t] - 2*(w_re * w_re - w_im * w_im)*(2 * w_re * w_im)*x_im[t];
				x_im[t] = ((w_re * w_re - w_im * w_im)*(w_re * w_re - w_im * w_im) - (2 * w_re * w_im)*(2* w_re * w_im))*x_im[t] + 2*(w_re * w_re - w_im * w_im)*(2 * w_re * w_im)*t_rt;
							
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
				
				// w_N = w_N_re + w_N_im*i
				// w_N^2 = (w_N_re*w_N_re-w_N_im*w_N_im) + (2*w_N_re*w_N_im)*i
				// w_N^3 = (w_N_re*w_N_re-w_N_im*w_N_im) - (2*w_N_re*w_N_im)*i
				// w_N^4 = w_N_re - w_N_im*i
				
				/*
					x_[q] = x_[p] + w_N * x_[q] + w_N^2 * x_[r] + w_N^3 * x_[s] + w_N^4 * x_[t] 
				=>  x_[q] = x_re[p] + (w_N_re * x_re[q] - w_N_im * x_im[q]) + ((w_N_re*w_N_re-w_N_im*w_N_im)*x_re[r] - (2*w_N_re*w_N_im)*x_im[r]) + ((w_N_re*w_N_re-w_N_im*w_N_im)*x_re[s] + (2*w_N_re*w_N_im)*x_im[s]) + (w_N_re * x_re[t] + w_N_im * x_im[t]))
					    +i*(x_im[p] + (w_N_re * x_im[q] + w_N_im * x_re[q]) + ((w_N_re*w_N_re-w_N_im*w_N_im)*x_im[r] + (2*w_N_re*w_N_im)*x_re[r]) + ((w_N_re*w_N_re-w_N_im*w_N_im)*x_im[s] + (2*w_N_re*w_N_im)*x_re[s]) + (w_N_re * x_im[t] - w_N_im * x_re[t]))
				=>  x_[q] = x_re[p] + w_N_re * (x_re[q] + x_re[t]) - w_N_im * (x_im[q] - x_im[t]) + (w_N_re*w_N_re-w_N_im*w_N_im)*(x_re[r] + x_re[s]) - (2*w_N_re*w_N_im) * (x_im[r] - x_im[s])
				        +i*(x_im[p] + w_N_re * (x_im[q] + x_im[t]) + w_N_im * (x_re[q] - x_re[t]) + (w_N_re*w_N_re-w_N_im*w_N_im)*(x_im[r] + x_im[s]) + (2*w_N_re*w_N_im) * (x_re[r] + x_re[s]))
				
					x_[r] = x_[p] + w_N^2 * x_[q] + w_N^4 * x_[r] + w_N * x_[s] + w_N^3 * x_[t] 
				=>	x_[r] = x_re[p] + (w_N_re * w_N_re-w_N_im*w_N_im)*x_re[q] - (2*w_N_re*w_N_im)*x_im[q]) + (w_N_re*x_re[r] + w_N_im * x_im[r]) + (w_N_re * x_re[s] - w_N_im * x_im[s]) + ((w_N_re*w_N_re-w_N_im*w_N_im)*x_re[t] + (2*w_N_re*w_N_im)*x_im[t])
						+i*(x_im[p] + (w_N_re * w_N_re-w_N_im*w_N_im)*x_im[q] + (2*w_N_re*w_N_im)*x_re[q]) + (w_N_re*x_im[r] - w_N_im * x_re[r]) + (w_N_re * x_im[s] + w_N_im * x_re[s]) + ((w_N_re*w_N_re-w_N_im*w_N_im)*x_im[t] - (2*w_N_re*w_N_im)*x_re[t]))
				=>	x_[r] = x_re[p] + w_N_re * (x_re[r] + x_re[s]) - w_N_im * (x_im[s] - x_im[r]) + (w_N_re * w_N_re-w_N_im*w_N_im)*(x_re[q] + x_re[t]) - (2*w_N_re*w_N_im)*(x_im[q] - x_im[t])
						+i*(x_im[p] + w_N_re * (x_im[r] + x_im[s]) + w_N_im * (x_re[s] - x_re[r]) + (w_N_re * w_N_re-w_N_im*w_N_im)*(x_im[q] + x_im[t]) + (2*w_N_re*w_N_im)*(x_re[q] - x_re[t]))
						
					x_[s] = x_[p] + w_N^3 * x_[q] + w_N * x_[r] + w_N^4 * x_[s] + w_N^2 * x_[t] 
				=>	x_[s] = x_re[p] + ((w_N_re*w_N_re-w_N_im*w_N_im)*x_re[q] + (2*w_N_re*w_N_im)*x_im[q]) + (w_N_re * x_re[r] - w_N_im * x_im[r]) + (w_N_re * x_re[s] + w_N_im * x_im[s]) + ((w_N_re*w_N_re-w_N_im*w_N_im)*x_re[t] - (2*w_N_re*w_N_im)*x_im[t])
						+i*(x_im[p] + ((w_N_re*w_N_re-w_N_im*w_N_im)*x_im[q] - (2*w_N_re*w_N_im)*x_re[q]) + (w_N_re * x_im[r] + w_N_im * x_re[r]) + (w_N_re * x_im[s] - w_N_im * x_re[s]) + ((w_N_re*w_N_re-w_N_im*w_N_im)*x_im[t] + (2*w_N_re*w_N_im)*x_re[t]))
				=>	x_[s] = x_re[p] + w_N_re * (x_re[r] + x_re[s]) - w_N_im * (x_im[r] - x_im[s]) + (w_N_re * w_N_re-w_N_im*w_N_im)*(x_re[t] + x_re[q]) - (2*w_N_re*w_N_im)*(x_im[t] - x_im[q])
						+i*(x_im[p] + w_N_re * (x_im[r] + x_im[s]) + w_N_im * (x_re[r] - x_re[s]) + (w_N_re * w_N_re-w_N_im*w_N_im)*(x_im[t] + x_im[q]) + (2*w_N_re*w_N_im)*(x_re[t] - x_re[q]))
					
					x_[t] = x_[p] + w_N^4 * x_[q] + w_N^3 * x_[r] + w_N^2 * x_[s] + w_N * x_[t] 
				=>	x_[t] = x_re[p] + (w_N_re * x_re[q] + w_N_im * x_im[q]) + ((w_N_re*w_N_re-w_N_im*w_N_im)*x_re[r] + (2*w_N_re*w_N_im)*x_im[r]) + ((w_N_re*w_N_re-w_N_im*w_N_im)*x_re[s] - (2*w_N_re*w_N_im)*x_im[s]) + (w_N_re * x_re[t] - w_N_im * x_im[t]) 
						+i*(x_im[p] + (w_N_re * x_im[q] - w_N_im * x_re[q]) + ((w_N_re*w_N_re-w_N_im*w_N_im)*x_im[r] + (2*w_N_re*w_N_im)*x_re[r]) + ((w_N_re*w_N_re-w_N_im*w_N_im)*x_im[s] + (2*w_N_re*w_N_im)*x_re[s]) + (w_N_re * x_im[t] + w_N_im * x_re[t]) )				
				=>	x_[t] = x_re[p] + w_N_re * (x_re[q] + x_re[t]) - w_N_im * (x_im[t] - x_im[q]) + (w_N_re * w_N_re-w_N_im*w_N_im)*(x_re[s] + x_re[r]) - (2*w_N_re*w_N_im)*(x_im[s] - x_im[r])
						+i*(x_im[p] + w_N_re * (x_im[q] + x_im[t]) + w_N_im * (x_re[t] - x_re[q]) + (w_N_re * w_N_re-w_N_im*w_N_im)*(x_im[s] + x_im[r]) + (2*w_N_re*w_N_im)*(x_re[s] - x_re[r]))
						
				*/
				
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
			/*
			temp = w_re;
			w_re = w_N_re * w_re - w_N_im * w_im;
			w_im = w_N_re * w_im + w_N_im *t;
			*/
		}	
	}
	return;
}



