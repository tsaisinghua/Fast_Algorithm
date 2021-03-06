#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#define DEBUG 0
#define DEBUG2 0
int median1(int *x, int left, int right);
int median2(int *x, int left, int right);
int partition(int *x, int left, int right);
int median_even(int *x, int left, int right);
int quicksort2(int *x, int left, int right);

int main()
{
	clock_t t1, t2;				// variables for computing clocks 
	int *x, *y, s;
	double T1, T2;
	int i, j, N, m, odd_mid_pos;
	int mid;
	srand( time(NULL) );

	for(N=999999;N<=1000000;N++)
	{
		x = (int *) malloc( N * sizeof(int) );
		y = (int *) malloc( N * sizeof(int) );

		for(i=0;i<N;++i)
		{
			y[i] = x[i] = rand() % N;
		}

		#if DEBUG
		for(i=0;i<N;++i)
		{
			printf("BEFORE:y[%d]=%d\n",i,y[i]);
		}
		#endif
		
		if(N%2==1)
		{
			t1 = clock();
			mid = median_odd(y,0,N);
			t2 = clock();
			T2 = (t2-t1)/(double) CLOCKS_PER_SEC;
			printf("Quick Median of %d elements: %f\n",N, T2);
			
			#if DEBUG		
			for(i=0;i<N;++i)
			{
				printf("AFTER:y[%d]=%d\n",i,y[i]);
			}
			#endif
			if (median_odd(y,0,N) < 0)
			{			
				printf("ERROR!\n");
			}
			else
			{
				printf("N=%d is odd, the median:y[%d]=%d\n",N,((N-1)/2), mid);
			}
		}
		else
		{
			t1 = clock();
			median_even(y,0,N);
			t2 = clock();
			T2 = (t2-t1)/(double) CLOCKS_PER_SEC;
			printf("Quick Median of %d elements: %f\n",N, T2);
			
			#if DEBUG		
			for(i=0;i<N;++i)
			{
				printf("AFTER:y[%d]=%d\n",i,y[i]);
			}
			#endif
			printf("N=%d is even, the median:(y[%d]+y[%d])/2.0=%f\n",N,(N/2)-1,N/2,(y[(N/2)-1]+y[N/2])/2.0);
		}		
		free(x);
		free(y);
	} 	
	return 0;
}
//=======================================  函    式    區  ========================================

int partition(int *x, int left, int right)
{
	int i, j, k;
	int pivot, t;
	pivot = x[left];
    i = left+1;
    j = right-1;
	
	if(left < right-1)
	{
		while(1)
		{	
			#if DEBUG2
			printf("BEFORE:i=%d j=%d pivot=%d\n", i,j,pivot);
			#endif
			
			while(i < right && pivot >= x[i]) i++; // 往右邊找到第一個  pivot <  x[i]  
      		while(j >  left && pivot <  x[j]) j--; // 往左邊找到第一個  pivot >= x[j] 
      		
			#if DEBUG2
			printf("AFTER:i=%d j=%d pivot=%d\n", i,j,pivot);
			#endif
				      		
			if(i>=j) break;
			
      		t = x[i];
      		x[i] = x[j];
      		x[j] = t;
      		// x[i] 與 x[j] 互換 
      		#if DEBUG2
			for(k=left;k<right;++k)
			{
				printf("x[%d]=%d\n",k,x[k]);
			}
			system("pause");
			#endif
      	}
        x[left] = x[j];
        x[j] = pivot;
        #if DEBUG2
        printf("i=%d,j=%d\n",i,j);
		for(k=left;k<right;++k)
		{
			printf("x[%d]=%d\n",k,x[k]);
		}
		system("pause");
        #endif
        return j;
    }
    else if(left == right-1)
    {
    	#if DEBUG2
    	printf("left=%d, right=%d\n", left, right);
    	#endif
		return left;		
	}	
	else 
	{
		return -1;
	}
}
//=========================================================================================

int median_even(int *x, int left, int right)
{
	if(left == right-1) //return x[left];
		return 1;
	int midIndex1 = (right-left)/2-1;
	int midIndex2 = (right-left)/2;
	int pivIndex = -1;
	int pivIndex1;
	int pivIndex2;
	float median;
	while(pivIndex != midIndex1 && pivIndex != midIndex2) 
	{
		pivIndex = partition(x, left, right);
		#if DEBUG2
		printf("pivIndex=%d, midIndex1=%d, midIndex2=%d\n", pivIndex, midIndex1, midIndex2);
		#endif
		if(pivIndex < 0)
			return pivIndex;
		
		if(pivIndex < midIndex1)	
		{
			left = pivIndex + 1;
			#if DEBUG2
			printf("pivIndex < midIndex1, left=%d, right=%d\n", left, right);
			#endif
		}			
		else if(pivIndex > midIndex2)
		{ 
			right = pivIndex;
			#if DEBUG2
			printf("pivIndex > midIndex2, left=%d, right=%d\n", left, right);
			#endif
		}
		else if(pivIndex == midIndex1)
		{
			left = pivIndex + 1;
			break;
		}
		else if(pivIndex == midIndex2)
		{
			right = pivIndex;
			break;
		}
		else
		{
			printf("ERROR!\n");
			break;
		}
		
	}
	
	if(pivIndex == midIndex1)
	{
		pivIndex1 = pivIndex;
		while(pivIndex != midIndex2) 
		{
			pivIndex = partition(x, left, right);
			#if DEBUG2
			printf("pivIndex=%d, midIndex1=%d, midIndex2=%d\n", pivIndex, midIndex1, midIndex2);
			#endif
			if(pivIndex < 0)
				return pivIndex;
			
			if(pivIndex < midIndex2)	
			{
				left = pivIndex + 1;
				#if DEBUG2
				printf("pivIndex < midIndex2, left=%d, right=%d\n", left, right);
				#endif
			}				
			else if(pivIndex > midIndex2)
			{ 
				right = pivIndex;
				#if DEBUG2
				printf("pivIndex > midIndex2, left=%d, right=%d\n", left, right);
				#endif
			}
			else break;
		}
		pivIndex2 = pivIndex;
	}
	else
	{
		pivIndex2 = pivIndex;
		while(pivIndex != midIndex1) 
		{
			pivIndex = partition(x, left, right);
			#if DEBUG2
			printf("pivIndex=%d, midIndex1=%d, midIndex2=%d\n", pivIndex, midIndex1, midIndex2);
			#endif
			if(pivIndex < 0)
				return pivIndex;
			
			if(pivIndex < midIndex1)	
			{
				left = pivIndex + 1;
				#if DEBUG2
				printf("pivIndex < midIndex1, left=%d, right=%d\n", left, right);
				#endif
			}
				
			else if(pivIndex > midIndex1)
			{ 
				right = pivIndex;
				#if DEBUG2
				printf("pivIndex > midIndex1, left=%d, right=%d\n", left, right);
				#endif
			}
			else break;
		}
		pivIndex1 = pivIndex;
	}
}

//=========================================================================================

int median_odd(int *x, int left, int right)
{
	if(left == right-1) return x[left];
	int midIndex = (right-left-1)/2;
	int pivIndex = -1;
	while(pivIndex != midIndex)
	{
		pivIndex = partition(x, left, right);
		#if DEBUG2
		printf("pivIndex=%d, midIndex=%d\n", pivIndex, midIndex);
		#endif
		if(pivIndex < 0)
			return pivIndex;
		
		if(pivIndex < midIndex)	
		{
			left = pivIndex + 1;
			#if DEBUG2
			printf("pivIndex < midIndex, left=%d, right=%d\n", left, right);
			#endif
		}			
		else if(pivIndex > midIndex)
		{ 
			right = pivIndex;
			#if DEBUG2
			printf("pivIndex > midIndex, left=%d, right=%d\n", left, right);
			#endif
		}
		else break;
	}
	return x[pivIndex];
}
