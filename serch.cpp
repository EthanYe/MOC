#include<stdio.h>
void serch(double gui, int *k1, int *k2,double A[],int len)
{
	int start=1, end= len-1, mid;
	while (start <= end)
	{
		mid = (start + end) / 2;
		if (A[mid + 1] < gui)
			start = mid + 1;
		else if (A[mid - 1] > gui)
			end = mid - 1;
		else if (A[mid] <= gui&&A[mid + 1] >= gui)
		{
			*k1 = mid, *k2 = mid + 1;			
			break;
		}
		else if (A[mid - 1] <= gui&&A[mid] >= gui)
		{
			*k1 = mid - 1, *k2 = mid;
			break;
		}
	}
}