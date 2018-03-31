#include<stdio.h>
#include<math.h>
double linear(double x, double x1, double x2, double y1, double y2)
{
	double y;
	y = y1 + (x - x1) / (x2 - x1)*(y2 - y1);
	return y;
}