#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#define DEBUG 1 
int partition(int *x, int left, int right);
int median2(int *x, int left, int right);


int main()
{
	clock_t t1, t2;				// variables for computing clocks 
	int *y, s;
	double T1, T2;
	int i, j, N, m, mid, odd_mid_pos;
	srand( time(NULL) );
	int x[11] = {3,0,2,5,7,9,10,8,6,10,3};
	for(N=11;N<=11;N*=2)
	{
		
		y = (int *) malloc( N * sizeof(int) );

		for(i=0;i<N;++i)
		{
			y[i] = x[i];
		}
		#if DEBUG					
		for(i=0;i<N;++i)			
		{
			printf("x[%d]=%d\n",i,x[i]);
		}
		#endif

		/*
		t1 = clock();
		mid = median2(y,0,N);
		t2 = clock();
		T2 = (t2-t1)/(double) CLOCKS_PER_SEC;
		printf("Quick Median2 of  %d elements: %f\n",N, T2);
		*/
		
	printf("=======================================\n");
	
		mid = median2(y,0,N);
		if (mid < 0)
			printf("ERROR!\n");
		
		else
		{
			#if DEBUG
			for(i=0;i<N;++i)
			{
				printf("y[%d]=%d\n",i,y[i]);
			}
			#endif
			printf("median2: y[%d]=%d\n",((N-1)/2), mid);
			
		}
		//free(x);
		free(y);
	} 	
	return 0;
}

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
			#if DEBUG
			printf("BEFORE:i=%d j=%d pivot=%d\n", i,j,pivot);
			#endif
			
			while(i < right && pivot >= x[i]) i++; // 往右邊找到第一個  pivot <  x[i]  
      		while(j >  left && pivot <  x[j]) j--; // 往左邊找到第一個  pivot >= x[j] 
      		
			#if DEBUG
			printf("AFTER:i=%d j=%d pivot=%d\n", i,j,pivot);
			#endif
				      		
			if(i>=j) break;
			
      		t = x[i];
      		x[i] = x[j];
      		x[j] = t;
      		// x[i] 與 x[j] 互換 
      		#if DEBUG
			for(k=left;k<right;++k)
			{
				printf("x[%d]=%d\n",k,x[k]);
			}
			system("pause");
			#endif
      	}
        x[left] = x[j];
        x[j] = pivot;
        #if DEBUG
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
    	printf("left=%d, right=%d\n", left, right);
		return left;		
	}	
	else 
	{
		return -1;
	}

}

//=========================================================================================

int median2(int *x, int left, int right)
{
	if(left == right-1) return x[left];
	int midIndex = (right-left-1)/2;
	int pivIndex = -1;
	while(pivIndex != midIndex)
	{
		pivIndex = partition(x, left, right);
		printf("pivIndex=%d, midIndex=%d\n", pivIndex, midIndex);
		if(pivIndex < 0)
			return pivIndex;
			
		if(pivIndex < midIndex)	
		{
			left = pivIndex + 1;
			printf("pivIndex < midIndex, left=%d, right=%d\n", left, right);
		}
			
		else if(pivIndex > midIndex)
		{ 
			right = pivIndex;
			printf("pivIndex > midIndex, left=%d, right=%d\n", left, right);
		}
		else break;
	}
	return x[pivIndex];
}

