#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#include "myfun.h" 
//#define DEBUG 0
//input a vector x (pointer), �ƧǪ����k���(Left, right) �i�����A�ݨ�ƪ��ˤl�O����A��@�b�᭱�j 
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
		for(i=0;i<N;++i)			// y[i]���n�ˬd�������Ai�ثe��j�p�A���Yy[i]<y[j],�h�N��̤�����m...i.e.�j���\�e���A�p���\�᭱�C 
		{
			//printf("x[%d]=%d\n",i,x[i]);
		}
//		#endif
		
		//��w�ƧǪk�G�ǭ��ˬd�@�U,  Operation counts: SWAP, COMPARE (N-1)+(N-2)+(N-3)+...+1
		t1 = clock();	//O(N^2) Operations.
		for(i=0;i<N;++i)
		{
			for(j=i+1;j<N;++j)
			{
				if(y[i]<y[j]) 	//�Ѥj��p�ƦC�F�Y�Q�Ѥp��j�A�N�N "<" �令 ">" �Y�i�C 
				{
					s = y[i];
					y[i] = y[j];
					y[j] = s;
				}
			}			
			//�ƦC��Ay[i]�N�|�O�� i �j���� 
		
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
		
		
		//���s���w y = x
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

		//���s���w y = x,
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
		// ����@�� pivot, �M���� pivot �p�� ���e��, ��pivot �j�����᭱
		// 	  �Q��k�� 
		// �A�h�� pivot �����k���� 
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


