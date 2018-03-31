#include<stdio.h>
#include<math.h>
void GAUSS(double(*A)[300], double X[], double B[], int n)
{
	int i, j, k;                                 //For loop
	int linenumber;                           //Choose line
	double MAXA = 0;
	double exchange;                          //For exchange line in characteristic matrix
	double Mik;                               //For elimination
	double tx;                                 //For backward substitution
	for (k = 1; k < n; k++)
	{
		for (i = k; i <= n; i++)      //Choose line number
		{
			if (fabs(*(*(A + i) + k)) > MAXA)
			{
				MAXA = fabs(*(*(A + i) + k));
				linenumber = i;
			}
			else linenumber = k;
		} //Choose line number
		for (i = 1; i <= n; i++)       //Exchange line in characteristic matrix
		{

			exchange = *(*(A + k) + i); *(*(A + k) + i) = *(*(A + linenumber) + i); *(*(A + linenumber) + i) = exchange;
		}
		exchange = *(B + k); *(B + k) = *(B + linenumber); *(B + linenumber) = exchange;
		for (i = k + 1; i <= n; i++)
		{
			if (*(*(A + k) + k) == 0)
				continue;
			Mik = *(*(A + i) + k) / *(*(A + k) + k);
			for (j = k; j <= n; j++)
				*(*(A + i) + j) = *(*(A + i) + j) - Mik*(*(*(A + k) + j));
			*(B + i) = *(B + i) - Mik*(*(B + k));
		}
		//Backward substitution
		for (i = n; i > 0; i--)
		{
			if (*(*(A + i) + i) == 0)
				break;
			tx = 0;
			for (j = i + 1; j <= n; j++)
				tx = tx + *(*(A + i) + j) * *(X + j);
			*(X + i) = (*(B + i) - tx) / *(*(A + i) + i);
		}
		return;
	}
}