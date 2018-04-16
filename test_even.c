#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#define DEBUG 1 
#define DEBUG2 0
int median2(int *x, int left, int right);
int median1(int *x, int left, int right);
int partition(int *x, int left, int right);
int median_even1(int *x, int left, int right);
int median_even2(int *x, int left, int right);
int median_even(int *x, int left, int right);
int quicksort2(int *x, int left, int right);

int main()
{
	clock_t t1, t2;				// variables for computing clocks 
	int *y, s;
	double T1, T2;
	int i, j, N=8, m, odd_mid_pos;
	double mid;
	int x[8]={0,7,5,9,2,10,6,3};
	srand( time(NULL) );

	//for(N=1000000;N<=1000000;N*=2)
	{
		//x = (int *) malloc( N * sizeof(int) );
		y = (int *) malloc( N * sizeof(int) );

		for(i=0;i<N;++i)
		{
			//y[i] = x[i] = rand() % N;
			y[i] = x[i];
		}
		#if DEBUG					
		for(i=0;i<N;++i)			
		{
			printf("x[%d]=%d\n",i,x[i]);
		}
		#endif
	
		//median 1
		#if DEBUG		
		for(i=0;i<N;++i)
		{
			printf("BEFORE:y[%d]=%d\n",i,y[i]);
		}
		#endif
		
		t1 = clock();
		median1(y,0,N);
		t2 = clock();
		T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
		printf("Quick Median1 of  %d elements: %f\n",N, T1);
		
		#if DEBUG		
		for(i=0;i<N;++i)
		{
			printf("AFTER:y[%d]=%d\n",i,y[i]);
		}
		#endif
		if(N%2==1)
		{
			odd_mid_pos = (N-1)/2;
			printf("N=%d is odd, the median1:y[%d]=%d\n",N, odd_mid_pos ,y[odd_mid_pos]);
		}
		else
		{
			m = N/2;
			printf("N=%d is even, the median1:(y[%d]+y[%d])/2.0=%f\n",N, m-1, m, (y[m-1]+y[m])/2.0);
		}
		
		printf("==============================================================\n");
		
		//median 2
		//重新指定 y = x,
		for(i=0;i<N;++i) y[i] = x[i];
		
		#if DEBUG
		for(i=0;i<N;++i)
		{
			printf("BEFORE:y[%d]=%d\n",i,y[i]);
		}
			#endif
		if (median_odd(y,0,N) < 0 || (median_even1(y,0,N) + median_even2(y,0,N)) < 0)
		{
			printf("ERROR!\n");
		}
		else
		{
			if(N%2==1)
			{
				t1 = clock();
				mid = median_odd(y,0,N);
				t2 = clock();
				T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
				printf("Odd Terms Quick Median2 of  %d elements: %f\n",N, T1);
				#if DEBUG		
				for(i=0;i<N;++i)
				{
					printf("AFTER:y[%d]=%d\n",i,y[i]);
				}
				#endif
				printf("N=%d is odd, the median2:y[%d]=%d\n",N,((N-1)/2), mid);
			}
			else
			{
				t1 = clock();
			//	mid = median_even(y,0,N);
				median_even(y,0,N);
				t2 = clock();
				T2 = (t2-t1)/(double) CLOCKS_PER_SEC;
				printf("Even Terms Quick Median2 of  %d elements: %f\n",N, T2);
				
				#if DEBUG		
				for(i=0;i<N;++i)
				{
					printf("AFTER:y[%d]=%d\n",i,y[i]);
				}
				#endif
				printf("N=%d is even, the median:(y[%d]+y[%d])/2.0=%f\n",N, (N/2)-1, N/2, (y[(N/2)-1]+y[N/2])/2.0);
			}		
		}
		
	//	free(x);
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
	if(left == right-1) return x[left];
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
		printf("pivIndex=%d, midIndex=%d\n", pivIndex, midIndex);
		#endif
		if(pivIndex < 0)
			return pivIndex;
		
		if(pivIndex < midIndex1)	
		{
			left = pivIndex + 1;
			#if DEBUG2
			printf("pivIndex < midIndex, left=%d, right=%d\n", left, right);
			#endif
		}
			
		else if(pivIndex > midIndex2)
		{ 
			right = pivIndex;
			#if DEBUG2
			printf("pivIndex > midIndex, left=%d, right=%d\n", left, right);
			#endif
		}
		else break;
	}
	
	if(pivIndex == midIndex1)
	{
		pivIndex1 = pivIndex;
		while(pivIndex != midIndex2) 
		{
			pivIndex = partition(x, left, right);
			#if DEBUG2
			printf("pivIndex=%d, midIndex=%d\n", pivIndex, midIndex);
			#endif
			if(pivIndex < 0)
				return pivIndex;
			
			if(pivIndex < midIndex2)	
			{
				left = pivIndex + 1;
				#if DEBUG2
				printf("pivIndex < midIndex, left=%d, right=%d\n", left, right);
				#endif
			}
				
			else if(pivIndex > midIndex2)
			{ 
				right = pivIndex;
				#if DEBUG2
				printf("pivIndex > midIndex, left=%d, right=%d\n", left, right);
				#endif
			}
			else break;
		}
		pivIndex2 = pivIndex;
	//	median = (x[pivIndex1]+x[pivIndex2])/2.0;
	}
	else
	{
		pivIndex2 = pivIndex;
		while(pivIndex != midIndex1) 
		{
			pivIndex = partition(x, left, right);
			#if DEBUG2
			printf("pivIndex=%d, midIndex=%d\n", pivIndex, midIndex);
			#endif
			if(pivIndex < 0)
				return pivIndex;
			
			if(pivIndex < midIndex1)	
			{
				left = pivIndex + 1;
				#if DEBUG2
				printf("pivIndex < midIndex, left=%d, right=%d\n", left, right);
				#endif
			}
				
			else if(pivIndex > midIndex1)
			{ 
				right = pivIndex;
				#if DEBUG2
				printf("pivIndex > midIndex, left=%d, right=%d\n", left, right);
				#endif
			}
			else break;
		}
		pivIndex1 = pivIndex;
	//	median = (x[pivIndex1]+x[pivIndex2])/2.0;
	}

	
//	return median;
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

//=========================================================================================

int median_even1(int *x, int left, int right)
{
	if(left == right-1) return x[left];
	int midIndex = (right-left)/2-1;
	//int midIndex2 = (right-left)/2;
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

//=========================================================================================

int median_even2(int *x, int left, int right)
{
	if(left == right-1) return x[left];
	int midIndex = (right-left)/2;
	//int midIndex2 = (right-left)/2;
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

//=========================================================================================

int quicksort2(int *x, int left, int right)
{
	int i, j, k;
	int pivot, t;
	
	if(left < right-1)
	{
		pivot = x[left];
    	i = left+1;
    	j = right-1;
       	while(1)
		{	
			while(i < right && pivot >= x[i]) i++; // 往右邊找到第一個  pivot <  x[i]  
      		while(j >  left && pivot <  x[j]) j--; // 往左邊找到第一個  pivot >= x[j] 
      		#if DEBUG2
			printf("i=%d j=%d pivot=%d\n", i,j,pivot);
			#endif
      		if(i>=j) break;
      		t = x[i];
      		x[i] = x[j];
      		x[j] = t;
      		#if DEBUG2
			for(k=left;k<right;++k)
			{
				printf("x[%d]=%d\n",k,x[k]);
			}
			system("pause");
			#endif
        }
        //t = x[left];
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
		quicksort2(x, left, j);
		quicksort2(x, j+1, right);
    }
    else 
    {
    	return 1;
	}	
}

//============================================================================== 
int median1(int *x, int left, int right)
{
	int i, j, k;
	int pivot, t;
	int m = (right-left-1)/2.0;
	if(left < right-1)
	{
		pivot = x[left];
    	i = left+1;
    	j = right-1;
       	while(1)
		{	
			while(i < right && pivot >= x[i]) i++; // 往右邊找到第一個  pivot <  x[i]  
      		while(j >  left && pivot <  x[j]) j--; // 往左邊找到第一個  pivot >= x[j] 
      		#if DEBUG2
			printf("left=%d, right=%d, i=%d j=%d pivot=%d\n", left, right, i, j, pivot);
			#endif
      		if(i>=j) break;
      		t = x[i];
      		x[i] = x[j];
      		x[j] = t;
      		#if DEBUG2
			for(k=left;k<right;++k)
			{
				printf("x[%d]=%d\n",k,x[k]);
			}
			system("pause");
			#endif
        }
        //t = x[left];
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
        if(j==m)		// j 為 pivot 的位置, m 為中位數的位置 
        	printf("median:x[%d]=%d\n",j, x[j]);
        else if(j>m)
			quicksort2(x, left, j);
		else	
			quicksort2(x, j+1, right);
    }
    else 
    {
    	return 1;
	}	
}

