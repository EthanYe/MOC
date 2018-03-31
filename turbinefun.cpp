#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"parameters.h"
#include"linear_interpolation.h"
void gui(int i, double *guide)
{
	int k1,j1;
	double ser;
	*guide = -1;
	for (j1 = 1;*guide < 0;j1++)
	{
		if ((T >= TUR[i].TimeServo[j1][0]) && (T <= TUR[i].TimeServo[j1 + 1][0]))
		{
			if (TUR[i].TimeServo[j1 + 1][0] == TUR[i].TimeServo[j1][0])
				ser = TUR[i].TimeServo[j1][1];
			else ser = TUR[i].TimeServo[j1][1] + (T - TUR[i].TimeServo[j1][0]) / (TUR[i].TimeServo[j1 + 1][0] - TUR[i].TimeServo[j1][0])*(TUR[i].TimeServo[j1 + 1][1] - TUR[i].TimeServo[j1][1]);
			for (k1 = 1;*guide < 0;k1++)
				if ((ser >= TUR[i].ServoGui[k1][0]) && (ser <= TUR[i].ServoGui[k1 + 1][0]))
					*guide = TUR[i].ServoGui[k1][1] + (ser - TUR[i].ServoGui[k1][0]) / (TUR[i].ServoGui[k1 + 1][0] - TUR[i].ServoGui[k1][0])*(TUR[i].ServoGui[k1 + 1][1] - TUR[i].ServoGui[k1][1]);
		}
	}
}
double SPLH(double gui, double xi, int k0)
{
	int m;
	int k1, k2;
	int mid, start = 1, end = TUR[k0].LCC;
	double g1, g2;
	double y1, y2, y;
	while (start <= end)
	{
		mid = (start + end) / 2;
		if (TUR[k0].guide[mid + 1] < gui)
			start = mid + 1;
		else if (TUR[k0].guide[mid - 1] > gui)
			end = mid - 1;
		else if (TUR[k0].guide[mid] <= gui&&TUR[k0].guide[mid + 1] >= gui)
		{
			k1 = mid, k2 = mid + 1;
			g1 = TUR[k0].guide[mid], g2 = TUR[k0].guide[mid + 1];
			break;
		}
		else if (TUR[k0].guide[mid - 1] <= gui&&TUR[k0].guide[mid] >= gui)
		{
			k1 = mid - 1, k2 = mid;
			g1 = TUR[k0].guide[mid - 1], g2 = TUR[k0].guide[mid];
			break;
		}
	}
	for (m = 2; m <= TUR[k0].NCC; m++)
	{
		if ((xi >= TUR[k0].CC[k1].ANG[m] && xi <= TUR[k0].CC[k1].ANG[m - 1]) || (xi <= TUR[k0].CC[k1].ANG[m] && xi >= TUR[k0].CC[k1].ANG[m - 1]))
		{
			y1 = TUR[k0].CC[k1].MH[m - 1] / 6 / TUR[k0].CC[k1].h[m] * pow(TUR[k0].CC[k1].ANG[m] - xi, 3) + TUR[k0].CC[k1].MH[m] / 6 / TUR[k0].CC[k1].h[m] * pow(xi - TUR[k0].CC[k1].ANG[m - 1], 3) + (TUR[k0].CC[k1].WH[m - 1] / TUR[k0].CC[k1].h[m] - TUR[k0].CC[k1].MH[m - 1] / 6 * TUR[k0].CC[k1].h[m])*(TUR[k0].CC[k1].ANG[m] - xi) + (TUR[k0].CC[k1].WH[m] / TUR[k0].CC[k1].h[m] - TUR[k0].CC[k1].MH[m] / 6 * TUR[k0].CC[k1].h[m])*(xi - TUR[k0].CC[k1].ANG[m - 1]);
			break;
		}
	}
	if (m == TUR[k0].NCC + 1)
	{
		printf("insect curve failed\n");
  		exit(0);
	}
	for (m = 2; m <= TUR[k0].NCC; m++)
	{
		if ((xi >= TUR[k0].CC[k2].ANG[m] && xi <= TUR[k0].CC[k2].ANG[m - 1]) || (xi <= TUR[k0].CC[k2].ANG[m] && xi >= TUR[k0].CC[k2].ANG[m - 1]))
		{
			y2 = TUR[k0].CC[k2].MH[m - 1] / 6 / TUR[k0].CC[k2].h[m] * pow(TUR[k0].CC[k2].ANG[m] - xi, 3) + TUR[k0].CC[k2].MH[m] / 6 / TUR[k0].CC[k2].h[m] * pow(xi - TUR[k0].CC[k2].ANG[m - 1], 3) + (TUR[k0].CC[k2].WH[m - 1] / TUR[k0].CC[k2].h[m] - TUR[k0].CC[k2].MH[m - 1] / 6 * TUR[k0].CC[k2].h[m])*(TUR[k0].CC[k2].ANG[m] - xi) + (TUR[k0].CC[k2].WH[m] / TUR[k0].CC[k2].h[m] - TUR[k0].CC[k2].MH[m] / 6 * TUR[k0].CC[k2].h[m])*(xi - TUR[k0].CC[k2].ANG[m - 1]);
			break;
		}
	}
	if (m == TUR[k0].NCC + 1)
	{
		printf("insect curve failed\n");
		exit(0);
	}
	y = linear(gui, g1, g2, y1, y2);
	return y;
}

double SPLB(double gui, double xi, int k0)
{
	int m;
	int k1, k2;
	int mid, start = 1, end = TUR[k0].LCC;
	double g1, g2;
	double y1, y2, y;
	while (start <= end)
	{
		mid = (start + end) / 2;
		if (TUR[k0].guide[mid + 1] < gui)
			start = mid + 1;
		else if (TUR[k0].guide[mid - 1] > gui)
			end = mid - 1;
		else if (TUR[k0].guide[mid] <= gui&&TUR[k0].guide[mid + 1] >= gui)
		{
			k1 = mid, k2 = mid + 1;
			g1 = TUR[k0].guide[mid], g2 = TUR[k0].guide[mid + 1];
			break;
		}
		else if (TUR[k0].guide[mid - 1] <= gui&&TUR[k0].guide[mid] >= gui)
		{
			k1 = mid - 1, k2 = mid;
			g1 = TUR[k0].guide[mid - 1], g2 = TUR[k0].guide[mid];
			break;
		}
	}
	for (m = 2; m <= TUR[k0].NCC; m++)
	{
		if ((xi >= TUR[k0].CC[k1].ANG[m] && xi <= TUR[k0].CC[k1].ANG[m - 1]) || (xi <= TUR[k0].CC[k1].ANG[m] && xi >= TUR[k0].CC[k1].ANG[m - 1]))
		{
			y1 = TUR[k0].CC[k1].MB[m - 1] / 6 / TUR[k0].CC[k1].h[m] * pow(TUR[k0].CC[k1].ANG[m] - xi, 3) + TUR[k0].CC[k1].MB[m] / 6 / TUR[k0].CC[k1].h[m] * pow(xi - TUR[k0].CC[k1].ANG[m - 1], 3) + (TUR[k0].CC[k1].WB[m - 1] / TUR[k0].CC[k1].h[m] - TUR[k0].CC[k1].MB[m - 1] / 6 * TUR[k0].CC[k1].h[m])*(TUR[k0].CC[k1].ANG[m] - xi) + (TUR[k0].CC[k1].WB[m] / TUR[k0].CC[k1].h[m] - TUR[k0].CC[k1].MB[m] / 6 * TUR[k0].CC[k1].h[m])*(xi - TUR[k0].CC[k1].ANG[m - 1]);
			break;
		}
	}
	if (m == TUR[k0].NCC +1)
	{
		printf("insect curve failed\n");
		exit(0);
	}
	for (m = 2; m <= TUR[k0].NCC; m++)
	{
		if ((xi >= TUR[k0].CC[k2].ANG[m] && xi <= TUR[k0].CC[k2].ANG[m - 1]) || (xi <= TUR[k0].CC[k2].ANG[m] && xi >= TUR[k0].CC[k2].ANG[m - 1]))
		{
			y2 = TUR[k0].CC[k2].MB[m - 1] / 6 / TUR[k0].CC[k2].h[m] * pow(TUR[k0].CC[k2].ANG[m] - xi, 3) + TUR[k0].CC[k2].MB[m] / 6 / TUR[k0].CC[k2].h[m] * pow(xi - TUR[k0].CC[k2].ANG[m - 1], 3) + (TUR[k0].CC[k2].WB[m - 1] / TUR[k0].CC[k2].h[m] - TUR[k0].CC[k2].MB[m - 1] / 6 * TUR[k0].CC[k2].h[m])*(TUR[k0].CC[k2].ANG[m] - xi) + (TUR[k0].CC[k2].WB[m] / TUR[k0].CC[k2].h[m] - TUR[k0].CC[k2].MB[m] / 6 * TUR[k0].CC[k2].h[m])*(xi - TUR[k0].CC[k2].ANG[m - 1]);
			break;
		}
	}
	if (m == TUR[k0].NCC + 1)
	{
		printf("insect curve failed\n");
		exit(0);
	}
	y = linear(gui, g1, g2, y1, y2);
	return y;
}