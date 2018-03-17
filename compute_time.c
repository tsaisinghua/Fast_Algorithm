#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <time.h>

int main()
{
	clock_t t1, t2;				// variables for computing clocks 
	double a=1.234, b=5.678;
	double T1, T2;
	int i, j, k, N=100000000;

	/* Addition and Subtraction x 3 + 1 loop */
	t1 = clock();
	for(i=0;i<N;++i)
	{
		a = a+b;
		b = a-b;
		a = a-b;
	}
	t2 = clock();
	T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
//	printf("(+,-) time x 3 + 1 loop:%f\n",T1);
//	printf("(a,b)=%f %f\n", a,b);
	
	/* Addition and Subtraction x 6 + 1 loop*/
	t1 = clock();
	for(i=0;i<N;++i)
	{
		a = a+b;
		b = a-b;
		a = a-b;
		a = a+b;
		b = a-b;
		a = a-b;
	}
	t2 = clock();
	T2 = (t2-t1)/(double) CLOCKS_PER_SEC;
//	printf("(+,-) time x 6 + 1 loop:%f\n",T2);
	printf("Real (+,-) time: %f\n",(T2-T1)/3.0);	
//	printf("(a,b)=%f %f\n", a,b);

	/* Multiplication & Division x 3 + 1 loop */
	t1 = clock();
	for(i=0;i<N;++i)
	{
		a = a*b;
		b = a/b;
		a = a/b;
	}
	t2 = clock();
	T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
//	printf("(*,/) time x 3 + 1 loop:%f\n",T1);
//	printf("(a,b)=%f %f\n", a,b);

	/* Multiplication & Division x 6 + 1 loop */
	t1 = clock();
	for(i=0;i<N;++i)
	{
		a = a*b;
		b = a/b;
		a = a/b;
		a = a*b;
		b = a/b;
		a = a/b;
	}
	t2 = clock();
	T2 = (t2-t1)/(double) CLOCKS_PER_SEC;
//	printf("(*,/) time x 6 + 1 loop:%f\n",T2);
	printf("Real (*,/) time: %f\n",(T2-T1)/3.0);	
//	printf("(a,b)=%f %f\n", a,b);

	/* sin x 1 + 1 loop */
	t1 = clock();
	for(i=0;i<N;++i)
	{
		a = sin(a);
	}
	t2 = clock();
	T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
//	printf("(sin) time + 1 loop:%f\n",T1);
//	printf("(a,b)=%f %f\n", a,b);

	/* sin x 2 + 1 loop */
	t1 = clock();
	for(i=0;i<N;++i)
	{
		a = sin(a);
		b = sin(b);
	}
	t2 = clock();
	T2 = (t2-t1)/(double) CLOCKS_PER_SEC;
//	printf("(sin) time x 2 + 1 loop:%f\n",T2);
	printf("(sin) time: %f\n",T2-T1);
//	printf("(a,b)=%f %f\n", a,b);

	return 0;
} 
