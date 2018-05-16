//
//  main.c
//  midterm
//
//  Created by Yu Hsun Lee on 2015/4/15.
//  Copyright (c) 2015�~ Yu Hsun Lee. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "sort.h"
#include "fft.c"
int main() {
	int N;
	int n, p, q, r;
	double *y_r, *y_i, *x_r, *x_i;
	clock_t t1, t2;
	/*
    printf("hello midterm \n");
	
	printf("input 2^p 3^q 5^r : p q r =>");
	scanf("%d %d %d", &p,&q,&r);
	N = 1 << p;
	N=N*pow(3,q)*pow(5, r);
	*/
	N = 33554432;
	printf("N=%d\n",N);
	
	x_r = (double *) malloc(N*sizeof(double));
	x_i = (double *) malloc(N*sizeof(double));
	y_r = (double *) malloc(N*sizeof(double));
	y_i = (double *) malloc(N*sizeof(double));
	
	//initial data
	for(n=0;n<N;++n)
	{
		x_r[n] = n;
		x_i[n] = 0;
	}
	
	t1 = clock();
	fft(x_r, x_i, y_r, y_i, N);
	t2 = clock();
	
	printf("%f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC); //print times
	//print_complex(y_r, y_i, N);
	
	free(x_r);
	free(x_i);
	free(y_r);
	free(y_i);
	//sort(v,N);
    return 0;
}
