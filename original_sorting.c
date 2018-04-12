#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#include "myfun.h" 
//#define DEBUG 0
//input a vector x (pointer), 排序的左右邊界(Left, right) 【先給你看函數的樣子是什麼，實作在後面】 
int quicksort1(int *x, int left, int right);
int quicksort2(int *x, int left, int right);

int main()
{
	clock_t t1, t2;				// variables for computing clocks 
	int *x, *y, s, p;
	double T1;
	int i, j, N;

	srand( time(NULL) );

	for(N=10000;N<=10000;N*=2)
	{
		x = (int *) malloc( N * sizeof(int) );
		y = (int *) malloc( N * sizeof(int) );

		//#pragma omp parallel
		{
			//#pragma omp parallel for
			for(i=0;i<N;++i)
			{
				y[i] = x[i] = rand() % N;
				y[i] = x[i] = rand() % (N*N);
			}
		}
//		#if DEBUG
		for(i=0;i<N;++i)			// y[i]為要檢查的元素，i目前比j小，但若y[i]<y[j],則將兩者互換位置...i.e.大的擺前面，小的擺後面。 
		{
			//printf("x[%d]=%d\n",i,x[i]);
		}
//		#endif
		
		//氣泡排序法：倆倆檢查一下,  Operation counts: SWAP, COMPARE (N-1)+(N-2)+(N-3)+...+1
		t1 = clock();	//O(N^2) Operations.
		for(i=0;i<N;++i)
		{
			for(j=i+1;j<N;++j)
			{
				if(y[i]<y[j]) 	//由大到小排列；若想由小到大，就將 "<" 改成 ">" 即可。 
				{
					s = y[i];
					y[i] = y[j];
					y[j] = s;
				}
			}			
			//排列後，y[i]就會是第 i 大的數 
		
			/* un-comment this for checking the order in each iteration
			printf("interation: %d\n", i);
			for(j=0;j<N;++j)
			{
				printf("y[%d]=%d\n",j,y[j]);
			}
			system("pause");
			*/
			
		}
		
		t2 = clock();
		T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
		printf("Sorting %d elements: %f\n",N, T1);
		
		for(i=0;i<N;++i)
		{
			//printf("y[%d]=%d\n",i,y[i]);
		}
		
		
		//重新指定 y = x
		for(i=0;i<N;++i) y[i] = x[i];
		
		t1 = clock();
		quicksort1(y,0,N);
//		#if DEBUG
		for(i=0;i<N;++i)
		{
			//printf("y[%d]=%d\n",i,y[i]);
		}
//		#endif
		t2 = clock();
		T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
		printf("(1)Quick Sorting %d elements: %f\n",N, T1);

		//重新指定 y = x,
		for(i=0;i<N;++i) y[i] = x[i];
		
		t1 = clock();
		quicksort2(y,0,N);
//		#if DEBUG
		for(i=0;i<N;++i)
		{
			//printf("y[%d]=%d\n",i,y[i]);
		}
//		#endif
		t2 = clock();
		T1 = (t2-t1)/(double) CLOCKS_PER_SEC;
		printf("(2)Quick Sorting %d elements: %f\n",N, T1);
		
		free(x);
		free(y);
	} 

	return 0;
}

int quicksort1(int *x, int left, int right)
{
	int i, j, k, pivot, pivot_loc, N = right-left; 
	int *y;
	if(left < right-1)
	{
		// 先找一個 pivot, 然後把比 pivot 小的 放到前面, 比pivot 大的放到後面
		// 	  想辦法做 
		// 再去排 pivot 的左右兩邊 
		//   quicksort(x, left, pivot_location);
		//   quicksort(x, pivot_location+1, right); 
		y = (int *) malloc(N*sizeof(int));
		pivot_loc = left+(rand() % N);
		pivot = x[pivot_loc];
		x[pivot_loc] = x[left];
		x[left] = pivot;
		i = 0; j = N-1;
		// x: 5(pivot) 4 4 3 2 7 8 9 10 -> y: 4 4 3 2 (i=j=4) 10 9 8 7
		// x: 4(pivot) 4 4 3 2 7 8 9 10 -> "<"   y: 3 2 (i=j=2) 10 9 8 7 4 4
		// x: 4(pivot) 4 4 3 2 7 8 9 10 -> "<="  y: 4 4 3 2 (i=j=4) 10 9 8 7
		// x: 5 4 4 3 2 7 8 9 10 , pivot_loc = 5, pivot = 7 
		// x: 7 4 4 3 2 5 8 9 10 -> y: 
		 
		for(k=1;k<N;++k) 
		{
			if(x[left+k] <= pivot) 
			{
				y[i++] = x[left+k];          //  i = i + 1; 
				//i = i + 1;                 //  y[i] = x[left+k]; --> y[++i] = x[left+k];
			}
			else
			{
				y[j--] = x[left+k];
				// j = j - 1;
			}
		}
		y[i] = pivot;
		#if DEBUG2
		printf("%d %d %d %d %d\n",left,i,j,pivot,N);
	 	for(k=0;k<N;++k)
		{
			//printf("y[%d]=%d\n",k,y[k]);
		}
		system("pause");
		#endif
		for(k=0;k<N;++k)
		{
			x[left+k] = y[k];
		}
		free(y);
		quicksort1(x,left,left+i);
		quicksort1(x,left+i+1,right);	
	}
	else
	{
		return 1;
	}
}


