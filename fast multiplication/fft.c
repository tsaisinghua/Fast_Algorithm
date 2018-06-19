#include <stdio.h>
#include <stdlib.h>
#include <time.h>
int main()
{
	int i, j, k, a, b, m, n, p, w, Iter_Max = 1000;
	int *x, *X, *y, *Y, *z, *Z;
	srand(time(NULL));
	n = 8;
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
		if(is_prime(n*k+1)&&n*k+1>81*n/2) 
		{
			p = n*k+1;
			//printf("%d",p);
			for(w=2;w<p;++w)
			{
				a = w;
				for(i=2;i<p;++i)
				{
					a = a*w % p;
					if(a == 1) break;
				}
				if(i==n) break;
			}
			if(w<p)
			{
				printf("Find!\n");
				break;
			}
		}
	}
	a = 1;
	for(i=0;i<n-1;++i)
	{
		a = a*w % p;
	}
	printf("p=%d for n=%d, w=%d, a=%d\n",p,n,w,a);
	system("pause");
	DFT(x,X,a,p,n,1);
	DFT(y,Y,a,p,n,1);
		
	for(i=0;i<n;++i)
	{
		Z[i] = X[i]*Y[i] % p;
	}
	DFT(Z,z,w,p,n,-1);
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
			s = (s + b*x[j]) % p;
			b = (b*a) % p;
		}
		y[i] = s;
		a = (a*w) % p;
	}
	if(dir==-1)
	{
		a = n;
		for(i=1;i<p;++i)
		{
			b = a*n % p;
			if(b==1) break;
			a = b;
		}
		printf("inverse of n: %d\n",a);
		for(i=0;i<n;++i)
		{
			y[i] = y[i]*a % p; 
		}
	}
	return 0;
}

