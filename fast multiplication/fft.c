#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define DEBUG 1
int Modtest(int a, int b, int p);
int FFT(int *x, int *y, int w, int p, int n, int dir);
int iFFT(int *x, int *y, int w, int p, int n, int dir);
int main()
{
	int i, j, k, a, b, m, n, p, w, Iter_Max = pow(2,30);
	int *x, *X, *y, *Y, *z, *Z;
	srand(time(NULL));
	n = 16;
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
		/*	
			x[i] = rand() % 10;
			y[i] = rand() % 10;
		*/	
			
			x[0] = 4;
			x[1] = 3;
			x[2] = 2;
			x[3] = 1;
			y[0] = 8;
			y[1] = 7;
			y[2] = 6;
			y[3] = 5;
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
					//printf("1.a= %d, w = %d\n",a, w);
					a = Modtest(a, w, p);
					//printf("2.a= %d, w = %d, p = %d\n",a, w, p);
					if(a == 1) break;
				}
				//printf("i= %d, n = %d\n", i, n);
				if(i==n) break;		 // i : W的次方, W=2,...,p-1 ; i = n =8 : w_{n=8}^{i=8}=1  
			}
			//printf("1.k= %d\n",k);
			if(w<p)					 //一執行此塊，碰到後面的break就會直接跳出最外面的for迴圈，因為break就是要為了跳出迴圈用的！ 
			{
				printf("Find!\n");
				break;				
			}
		}
		//printf("2.k= %d\n",k);		// 當 k = 42 時，因為前面已break，所以這裡不會被執行！
	}
	printf("p=%d for n=%d, w=%d, a=%d\n",p,n,w,a);
	a = 1;
	for(i=0;i<n-1;++i)
	{
		a = Modtest(a, w, p);					// w^n=1 (mod p) => 找 w^{n-1}=? (mod p) => ? = 反元素！ 
	}
	printf("p=%d for n=%d, w=%d, a=%d\n",p,n,w,a);
	system("pause");
	#if DEBUG
	for(i=n-1;i>=0;i--)
	{
		printf("Before:X[%d]=%d\n",i,X[i]);
	}
	for(i=n-1;i>=0;i--)
	{
		printf("Before:Y[%d]=%d\n",i,Y[i]);
	}
	#endif
	// DFT(x,X,a,p,n,1);
	FFT(x,X,a,p,n,1);
	// DFT(y,Y,a,p,n,1);
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
	for(i=n-1;i>=0;i--)
	{
		printf("After:Z[%d]=%d\n",i,Z[i]);
	}
	// DFT(Z,z,w,p,n,-1);
	iFFT(Z,z,w,p,n,-1);
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
/*
int DFT(int *x, int *y, int w, int p, int n, int dir)
{
	int i, j, s, a, b;
	a = 1;
	for(i=0;i<n;++i)
	{
		s = 0; 
		b = 1;
		for(j=0;j<n;++j)
		{
			s = s % p + Modtest(b, x[j], p);
			b = Modtest(b, a, p);
		}
		y[i] = s;
		a = Modtest(a, w, p);
	}
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
			//printf("2.y[i] = %d\n",y[i]);
		}
	}
	return 0;
}
*/
int FFT(int *x, int *y, int w, int p, int n, int dir)
{
	int i, j, a, b;
	
	if(n==2)
	{	
		y[0] = (x[0] + x[1]) %p;
		y[1] = x[0]%p + Modtest(w, x[1], p);
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
		FFT(z, u, Modtest(w, w, p), p, n/2, 1);
		FFT(z+n/2, u+n/2, Modtest(w, w, p), p, n/2, 1);
		t = 1;
		
		for (j=0; j<n/2; ++j)
		{
      		y[j] = (u[j] + t * u[n/2+j]) % p;
      		y[n/2+j] = (u[j] - t * u[n/2+j]) % p;
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
	else printf("we didn't do this!\n");
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

int iFFT(int *x, int *y, int w, int p, int n, int dir)
{
	int i, j, a, b;	
	if(n==2)
	{	
		y[0] = (x[0] + x[1]) %p;
		y[1] = (x[0] +w*x[1]) %p;
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
		iFFT(z, u, Modtest(w, w, p), p, n/2, 1);
		iFFT(z+n/2, u+n/2, Modtest(w, w, p), p, n/2, 1);
		t = 1;
		
		for(j=0; j<n/2; ++j)
		{
      		y[j] = u[j]% p + Modtest(t, u[n/2+j], p);
      		y[n/2+j] = u[j]% p - Modtest(t, u[n/2+j], p);
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

