#include <stdio.h>

int main()
{
	bit_reverse(16);
	return;
}

int bit_reverse(int N)
{
    int m,p,q,k;
    m = N/2;                        // Bit-Reverse �C���n�i�쪺�Ʀr 
    q = m;							// p = 1, q = m (�Ĥ@�ӭn�洫��) 
    for(p=1;p<N-1;++p)
    {
        printf("%d <-> %d\n", p,q);
        if(p < q)
        {
            //swap p and q
        }
        k = m;						// k, �Ψ��ˬd�� log_2 k + 1 ��O���O1 
        while(q >= k & k > 0)		// q >=k �� (log_2 k + 1)��O1,  
        {
            q = q-k;				// 1->0
            k = k/2;				// �ˬd�U�@�� 
        }
        q = q+k;
    }
    return 0;
}
