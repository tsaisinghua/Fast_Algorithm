#include <stdio.h>
#include <stdlib.h>

int main(void)
{
 	int a = 40000;
 	int b = 60000;
	int i;
 	int c;
 	int p=50000;
 	c=a*b %p;
 	
 	printf("c=%d\n", c);
 	
 	for(i=0;i<=b;++i)
 	{
 		a=(a+a) %p;	
	}
	printf("c=%d", c);
	return 0;	
}

