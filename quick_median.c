#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#define DEBUG 1 
#define DEBUG2 0
int quicksort2(int *x, int left, int right);
int median1(int *x, int left, int right);
int median2(int *x, int left, int right);
int partition(int *x, int left, int right);

int main()
{
	clock_t t1, t2;				// variables for computing clocks 
	int *x, *y, s;
	double T1, T2;
	int i, j, N, m, odd_mid_pos;
	srand( time(NULL) );

	for(N=17;N<=17;N*=2)
	{
		x = (int *) malloc( N * sizeof(int) );
		y = (int *) malloc( N * sizeof(int) );

		//#pragma omp parallel
		{
			//#pragma omp parallel for
			for(i=0;i<N;++i)
			{
				y[i] = x[i] = rand() % N;
				//y[i] = x[i] = rand() % (N*N);
			}
		}
		#if DEBUG					
		for(i=0;i<N;++i)			
		{
			printf("x[%d]=%d\n",i,x[i]);
		}
		#endif
		
		
		//median 1
		t1 = clock();
		median1(y,0,N);
		t2 = clock();
		T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
		printf("Quick Sorting %d elements: %f\n",N, T1);
		
		#if DEBUG		
		for(i=0;i<N;++i)
		{
			printf("y[%d]=%d\n",i,y[i]);
		}
		#endif
		if(N%2==1)
		{
			odd_mid_pos = (N-1)/2;
			printf("N=%d is odd, the median:y[%d]=%d\n",N, odd_mid_pos ,y[odd_mid_pos]);
		}
		else
		{
			m = N/2;
			printf("N=%d is even, median:(y[%d]+y[%d])/2.0=%f\n",N, m-1, m, (y[m-1]+y[m])/2.0);
		}
		
		printf("==============================================================\n");
		
		//median 2
		//���s���w y = x,
		for(i=0;i<N;++i) y[i] = x[i];
		
		t1 = clock();
		median2(y,0,N);
		t2 = clock();
		T2 = (t2-t1)/(double) CLOCKS_PER_SEC;
		printf("Quick Sorting %d elements: %f\n",N, T2);
		
		#if DEBUG		
		for(i=0;i<N;++i)
		{
			printf("y[%d]=%d\n",i,y[i]);
		}
		#endif
		odd_mid_pos = 0;
		if(N%2==1)
		{
			odd_mid_pos = (N-1)/2;
			printf("N=%d is odd, the median:y[%d]=%d\n",N, odd_mid_pos ,y[odd_mid_pos]);
		}
		else
		{
			m = N/2;
			printf("N=%d is even, median:(y[%d]+y[%d])/2.0=%f\n",N, m-1, m, (y[m-1]+y[m])/2.0);
		}
		
		free(x);
		free(y);
	} 	
	return 0;
}

//=============================  ��    ��    ��  ===================================

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
			// x: 5(pivot) 4 8 3 2 7 3 2 10  -> i = 2, j = 7 (�洫 x[i], x[j]) 
			//    5(pivot) 4 2 3 2 7 3 8 10  -> i = 5, j = 6 (�洫 x[i], x[j])
			//    5(pivot) 4 2 3 2 3 7 8 10  -> i = 6, j = 5 (���洫�F!!!) 
			//    3                5 (location: j) 
			// x: 4(pivot) 4 8 3 2 7 3 4 10  -> i = 2, j = 7 (�洫 x[i], x[j]) 
			//    4(pivot) 4 4 3 2 7 3 8 10  -> i = 5, j = 6 (�洫 x[i], x[j])
			//    4(pivot) 4 4 3 2 3 7 8 10  -> i = 6, j = 5 (���洫�F!!!) 
			//    3                4 (location: j) 			
			// x: 8 2 1 8 7 8 9 4 5 8      i = 3, j = 9 (�洫 x[i], x[j])
			// x: 8 2 1 8 7 8 9 4 5 8
      		while(i < right && pivot >= x[i]) i++; // ���k����Ĥ@��  pivot <  x[i]  
      		while(j >  left && pivot <  x[j]) j--; // ��������Ĥ@��  pivot >= x[j] 
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
			while(i < right && pivot >= x[i]) i++; // ���k����Ĥ@��  pivot <  x[i]  
      		while(j >  left && pivot <  x[j]) j--; // ��������Ĥ@��  pivot >= x[j] 
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
        if(j==m)		// j �� pivot ����m, m ������ƪ���m 
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
//=========================================================================================
int partition(int *x, int left, int right)
{
	int i, j, k;
	int pivot, t; 
	pivot = x[left];
	i = left+1;
	j = right-1;
	
	if(left >= right-1)
	{
		return -1;
	}	
	
	while(left < right-1)
	{
			
	    while(i < right && pivot >= x[i]) i++; // ���k����Ĥ@��  pivot <  x[i]  
	    while(j >  left && pivot <  x[j]) j--; // ��������Ĥ@��  pivot >= x[j] 
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

//=========================================================================================

int median2(int *x, int left, int right)
{
	int i, j, k;
	int pivot, t;
	int midIndex = (right-left-1)/2;
	int pivIndex = -1;
	while(pivIndex != midIndex)
	{
		pivIndex = partition(x, left, right);
		
		if(pivIndex < midIndex)
			left = pivIndex + 1;
		else if(pivIndex > midIndex)
			right = pivIndex - 1;
		else break;
	}
	return x[pivIndex];	
}

