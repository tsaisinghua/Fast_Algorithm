#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define DEBUG 0
int Modtest(int a, int b, int p);
int FFT(int *x, int *y, int w, int p, int n, int dir);
int iFFT(int *x, int *y, int w, int p, int n, int dir);
int main()
{
	int i, j, k, a, b, m, n, p, w, Iter_Max = pow(2,30);
	int *x, *X, *y, *Y, *z, *Z;
	srand(time(NULL));
	n = 150;
	x = (int *) malloc(n*sizeof(int));
	X = (int *) malloc(n*sizeof(int));
	y = (int *) malloc(n*sizeof(int));
	Y = (int *) malloc(n*sizeof(int));
	z = (int *) malloc(n*sizeof(int));
	Z = (int *) malloc(n*sizeof(int));
	
	// initial value
	for(i=0;i<n;++i)  
	{
		
		if(i<n/2) 
		{
	
			x[i] = rand() % 10;
			y[i] = rand() % 10;
		
		/*		
			x[0] = 5;
			x[1] = 4;
			x[2] = 3;
			x[3] = 2;
			x[4] = 1;
			y[0] = 0;
			y[1] = 9;
			y[2] = 8;
			y[3] = 7;
			y[4] = 6;
		*/
		} 
		else 
		{
			x[i] = 0;
			y[i] = 0;
			z[i] = 0;
		}
	}
	
	for(i=n/2-1;i>=0;i--)
	{
		printf("%d",x[i]);
	}
	printf("\n");
	for(i=n/2-1;i>=0;i--)
	{
		printf("%d",y[i]);
	}
	printf("\n");
	
	for(k=1;k<=Iter_Max;++k)
	{
		if(is_prime(n*k+1) && n*k+1>81*n/2) 
		{
			p = n*k+1;
			//printf("%d",p);
			for(w=2;w<p;++w)
			{
				a = w;
				for(i=2;i<p;++i)
				{	
					#if DEBUG
					printf("1.a= %d, w = %d\n",a, w);
					#endif
					a = Modtest(a, w, p);
					
					#if DEBUG					
					printf("2.a= %d, w = %d, p = %d\n",a, w, p);
					#endif
					if(a == 1) break;
				}
				#if DEBUG
				printf("i= %d, n = %d\n", i, n);
				#endif
				if(i==n) break;		 // i : W的次方, W=2,...,p-1 ; i = n =8 : w_{n=8}^{i=8}=1  
			}
			#if DEBUG
			printf("1.k= %d\n",k);
			#endif
			if(w<p)					 //一執行此塊，碰到後面的break就會直接跳出最外面的for迴圈，因為break就是要為了跳出迴圈用的！ 
			{
				printf("Find!\n");
				break;				
			}
		}
		#if DEBUG
		printf("2.k= %d\n",k);		// 當 k = 42 時，因為前面已break，所以這裡不會被執行！
		#endif
	}
	printf("p=%d for n=%d, w=%d, a=%d\n",p,n,w,a);
	a = 1;
	for(i=0;i<n-1;++i)
	{
		a = Modtest(a, w, p);					// w^n=1 (mod p) => 找 w^{n-1}=? (mod p) => ? = 反元素！ 
	}
	printf("p=%d for n=%d, w=%d, a=%d\n",p,n,w,a);
	system("pause");
	
	FFT(x,X,a,p,n,1);
	FFT(y,Y,a,p,n,1);
	
	#if DEBUG
	for(i=n-1;i>=0;i--)
	{
		printf("After:X[%d]=%d\n",i,X[i]);
	}
	for(i=n-1;i>=0;i--)
	{
		printf("After:Y[%d]=%d\n",i,Y[i]);
	}
	#endif
	for(i=0;i<n;++i)
	{
		Z[i] =Modtest(X[i], Y[i], p);
	}
	#if DEBUG
	for(i=n-1;i>=0;i--)
	{
		printf("After:Z[%d]=%d\n", i, Z[i]);
	}
	#endif
	
	FFT(Z,z,w,p,n,-1);
	
	for(i=n-1;i>=0;i--)
	{
		printf("%d ",z[i]);
	}
	printf("\n");
	for(i=0;i<n-1;i++)
	{
		z[i+1] += z[i]/10;
		z[i] = z[i] % 10;
	}
	for(i=n-1;i>=0;i--)
	{
		printf("%d",z[i]);
	}
	printf("\n");

	return 0;
}
int is_prime(int p)		// 若 P 為質數，則回傳 1，這樣 main 中的 for 才會執行！ 
{
	int a;
	for (a=2;a*a<=p;++a)
	{
		if(p % a == 0) return 0;
	}
	return 1;
}
int FFT(int *x, int *y, int w, int p, int n, int dir)
{
	int i, j, a, b, c, d, e, f;
	int m2, m3, m4, m5;
	int t2, t3, t4;
	int mtup1, mtup2, mtup3, mtup4;
	if(n==2)
	{	
		y[0] = (x[0] + x[1]) % p;
		y[1] = (x[0] % p + Modtest(w, x[1], p))%p;
	}	
	else if(n>2 && (n%2)==0)
	{
		int *z, *u;
		int k, t;
		z = (int *) malloc(n*sizeof(int));	// before
		u = (int *) malloc(n*sizeof(int));	// after
		for(k=0;k<n/2;++k)
		{
			z[k] = x[2*k];			// x_even
			z[n/2+k] = x[2*k+1];	// x_odd
		}
		m2 = Modtest(w, w, p);
		FFT(z, u, m2, p, n/2, 1);
		FFT(z+n/2, u+n/2, m2, p, n/2, 1);
		t = 1;		
		for (j=0; j<n/2; ++j)
		{
      		y[j] = u[j] % p + Modtest(t, u[n/2+j], p);
      		y[n/2+j] = u[j] % p - Modtest(t, u[n/2+j], p);
      		if(y[j]<0)
      		{
      			y[j] += p;
			}
			if(y[n/2+j]<0)
      		{
      			y[n/2+j] += p;
			}
      		t = Modtest(t, w, p);
		}
  		free(u);
		free(z);
		//return y;
	}
	else if(n==3)
	{
		m2 = Modtest(w, w, p);
		y[0] = (x[0] + x[1] + x[2]) % p;
		y[1] = x[0]%p + Modtest(w, x[1], p) + Modtest(m2, x[2], p);
		y[2] = x[0]%p + Modtest(m2, x[1], p) + Modtest(w, x[2], p);
	}
	else if(n>3 && (n%3)==0)
	{
		int *z, *u;
		int k, t;
		z = (int *) malloc(n*sizeof(int));	// before
		u = (int *) malloc(n*sizeof(int));	// after
		
		for(k=0;k<n/3;++k)
		{
			z[k] = x[3*k];
			z[n/3+k]  = x[3*k+1];
			z[2*n/3+k]  = x[3*k+2];
		}
		
		m3 = Modtest(Modtest(w, w, p), w, p);
		FFT(z, u, m3, p, n/3, 1);
		FFT(z+n/3, u+n/3, m3, p, n/3, 1);
		FFT(z+2*n/3, u+2*n/3, m3, p, n/3, 1);
		
		c=1;
		for(i=0;i<n/3;++i)
		{
			c=Modtest(c, w, p);
		}
		d=1;
		for(i=0;i<2*n/3;++i)
		{
			d=Modtest(d, w, p);
		}
		t = 1;
		for(j=0;j<n/3;++j)
		{
			t2 = Modtest(t, t, p);
      		y[j] = u[j] %p + Modtest(t, u[n/3+j], p) + Modtest(t2, u[2*n/3+j], p);
			y[n/3+j] = u[j] %p + Modtest(Modtest(t, u[n/3+j], p),c,p) + Modtest(Modtest(t2, u[2*n/3+j], p),d,p);
      		y[2*n/3+j] = u[j] %p + Modtest(Modtest(t, u[n/3+j], p),d,p) + Modtest(Modtest(t2, u[2*n/3+j], p),c,p);
      		
			if(y[j]<0)
      		{
      			y[j] += p;
			}
			if(y[n/3+j]<0)
      		{
      			y[n/3+j] += p;
			}
			if(y[2*n/3+j]<0)
      		{
      			y[2*n/3+j] += p;
			}
      		t = Modtest(t, w, p);
		}
		free(u); free(z);
	}
	else if(n==5)
	{
		m2 = Modtest(w, w, p);
		m3 = Modtest(m2, w, p);
		m4 = Modtest(m3, w, p);
		y[0] = (x[0]+x[1]+x[2]+x[3]+x[4]) % p;
		y[1] = (x[0]%p + Modtest(w, x[1], p) + Modtest(m2, x[2], p) + Modtest(m3, x[3], p) + Modtest(m4, x[4], p))%p;
		y[2] = (x[0]%p + Modtest(m2, x[1], p) + Modtest(m4, x[2], p) + Modtest(w, x[3], p) + Modtest(m3, x[4], p))%p;
		y[3] = (x[0]%p + Modtest(m3, x[1], p) + Modtest(w, x[2], p) + Modtest(m4, x[3], p) + Modtest(m2, x[4], p))%p;	
		y[4] = (x[0]%p + Modtest(m4, x[1], p) + Modtest(m3, x[2], p) + Modtest(m2, x[3], p) + Modtest(w, x[4], p))%p;	
	}
	else if(n>5 && (n%5)==0)
	{
		int *z, *u;
		int k, t;
		z = (int *) malloc(n*sizeof(int));	// before
		u = (int *) malloc(n*sizeof(int));	// after
		
		for(k=0;k<n/5;++k)
		{
			z[k] = x[5*k];
			z[n/5+k] = x[5*k+1];
			z[2*n/5+k] = x[5*k+2];
			z[3*n/5+k] = x[5*k+3];
			z[4*n/5+k] = x[5*k+4];
		}
		m5 = Modtest(Modtest(Modtest(Modtest(w, w, p), w, p), w, p), w, p);
		FFT(z, u, m5, p, n/5, 1);
		FFT(z+n/5, u+n/5, m5, p, n/5, 1);
		FFT(z+2*n/5, u+2*n/5, m5, p, n/5, 1);
		FFT(z+3*n/5, u+3*n/5, m5, p, n/5, 1);
		FFT(z+4*n/5, u+4*n/5, m5, p, n/5, 1);
		
		c = 1;
		for(i=0;i<n/5;++i)
		{
			c = Modtest(c, w, p);
		}
		d = 1;
		for(i=0;i<2*n/5;++i)
		{
			d = Modtest(d, w, p);
		}
		e = 1;
		for(i=0;i<3*n/5;++i)
		{
			e = Modtest(e, w, p);
		}
		f = 1;
		for(i=0;i<4*n/5;++i)
		{
			f=Modtest(f, w, p);
		}
		t = 1;
		for(j=0;j<n/5;++j)
		{
			t2 = Modtest(t, t, p);
			t3 = Modtest(t2, t, p);
			t4 = Modtest(t3, t, p);
			mtup1 = Modtest(t, u[n/5+j], p);
			mtup2 = Modtest(t2, u[2*n/5+j], p);
			mtup3 = Modtest(t3, u[3*n/5+j], p);
			mtup4 = Modtest(t4, u[4*n/5+j], p);
      		y[j] = (u[j] %p + mtup1 + mtup2 + Modtest(t3, u[3*n/5+j], p) +Modtest(t4, u[4*n/5+j], p))%p;
			y[n/5+j] = (u[j] %p + Modtest(mtup1,c,p) + Modtest(mtup2,d,p) + Modtest(mtup3,e,p) + Modtest(mtup4,f,p))%p;
      		y[2*n/5+j] = (u[j] %p + Modtest(mtup1,d,p) + Modtest(mtup2,f,p) + Modtest(mtup3,c,p) + Modtest(mtup4,e,p))%p;
      		y[3*n/5+j] = (u[j] %p + Modtest(mtup1,e,p) + Modtest(mtup2,c,p) + Modtest(mtup3,f,p) + Modtest(mtup4,d,p))%p;
      		y[4*n/5+j] = (u[j] %p + Modtest(mtup1,f,p) + Modtest(mtup2,e,p) + Modtest(mtup3,d,p) + Modtest(mtup4,c,p))%p;
			if(y[j]<0)
      		{
      			y[j] += p;
			}
			if(y[n/5+j]<0)
      		{
      			y[n/5+j] += p;
			}
			if(y[2*n/5+j]<0)
      		{
      			y[2*n/5+j] += p;
			}
			if(y[3*n/5+j]<0)
      		{
      			y[3*n/5+j] += p;
			}
			if(y[4*n/5+j]<0)
      		{
      			y[4*n/5+j] += p;
			}
      		t = Modtest(t, w, p);
		}
		free(u); free(z);
	}
	else printf("we didn't do this!\n");
	if(dir==-1)
	{
		a = n;		
		for(i=1;i<p;++i)
		{
			b = Modtest(a, n, p);
			if(b==1) break;
			a = b;
		}		
		printf("inverse of n: %d\n",a);
		for(i=0;i<n;++i)
		{
			y[i] = Modtest(y[i], a, p);
		}
	}
	/*
	else if(n==3)
	{
		y[0] = (x[0] + x[1] + x[2]) % p;
		y[1] = (x[0] + w*x[1] + w*w*x[2]) % p;
		y[2] = (x[0] + w*w*x[1] + w*x[2]) % p;
	}
	*/
	return 0;
}

int Modtest(int a, int b, int p)
{
	// result = (a*b) % p
	// a * b = a + a + a + ... + a (加 b 次 	 
	int result = 0;
	int ceil = 0;
	int i;
	if(a==0)
	{
		return 0;
	}
	if(a < p)
		ceil = p/a + 1;
	else
		ceil = 1;
	
	// 找餘數 
	int remainder = b % ceil;	
	int sum = 0;
		
	for(i = 0; i < b/ceil ; i++)
	{
		int r = (a * ceil) % p;
		sum = (sum + r) %p;
	}
	
	// 把剩餘的餘數求出來 
	int r = (a*remainder) % p;
	
	sum = (sum + r) %p;
	result = sum;
	return result;
}

