
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
			// x: 5(pivot) 4 8 3 2 7 3 2 10  -> i = 2, j = 7 (ユ传 x[i], x[j]) 
			//    5(pivot) 4 2 3 2 7 3 8 10  -> i = 5, j = 6 (ユ传 x[i], x[j])
			//    5(pivot) 4 2 3 2 3 7 8 10  -> i = 6, j = 5 (ぃユ传!!!) 
			//    3                5 (location: j) 
			// x: 4(pivot) 4 8 3 2 7 3 4 10  -> i = 2, j = 7 (ユ传 x[i], x[j]) 
			//    4(pivot) 4 4 3 2 7 3 8 10  -> i = 5, j = 6 (ユ传 x[i], x[j])
			//    4(pivot) 4 4 3 2 3 7 8 10  -> i = 6, j = 5 (ぃユ传!!!) 
			//    3                4 (location: j) 			
			// x: 8 2 1 8 7 8 9 4 5 8      i = 3, j = 9 (ユ传 x[i], x[j])
			// x: 8 2 1 8 7 8 9 4 5 8
      		while(i < right && pivot >= x[i]) i++; // ┕娩т材  pivot <  x[i]  
      		while(j >  left && pivot <  x[j]) j--; // ┕オ娩т材  pivot >= x[j] 
      		#if DEBUG2
			printf("%d %d %d\n", i,j,pivot);
			#endif
      		if(i>=j) break;
      		t = x[i];
      		x[i] = x[j];
      		x[j] = t;
      		#if DEBUG2
			for(k=left;k<right;++k)
			{
				//printf("x[%d]=%d\n",k,x[k]);
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
			//printf("x[%d]=%d\n",k,x[k]);
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

