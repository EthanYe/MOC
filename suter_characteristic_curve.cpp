#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"parameters.h"
#include"linear_interpolation.h"
#include"gauss.h"

int ncc;
double n_z[300][300], Q_z[300][300], M_z[300][300];      //Relative unit parameter    a,v,¦Â
double n01, Q01, M01;                      //Rated unit parameter
FILE *fpt,*fpc;
void sutertrans(int,int);


void suter(int m,char *pro)
{
	int m1, m2;
	char buf[256],name[20];
	//Calculate rate unit paras
	n01 = TUR[m].n0*TUR[m].D1 / sqrt(TUR[m].H0);
	Q01= TUR[m].Q0/ TUR[m].D1/ TUR[m].D1/ sqrt(TUR[m].H0);
	M01 = TUR[m].M0 / pow(TUR[m].D1, 3) / TUR[m].H0;
	//Open name_cc.txt file
	strcpy_s(name, pro);
	strcat_s(name,20, "_cc.txt");
	fopen_s(&fpt,name, "r");
	//Make sure the file exists
	if (fpt == NULL)                                  
	{
		printf("Openning file TDT.txt failed\n");
		exit(0);
	}
	printf("Open name_cc file successfully\n");
	//Define LCC and NCC
	fscanf_s(fpt, "%d%d", &TUR[m].LCC, &TUR[m].NCC);
	for (m1 = 1;m1 <= TUR[m].LCC;m1++)
		fscanf_s(fpt, "%lf", &TUR[m].guide[m1]);
	ncc = TUR[m].NCC;
	//input n
	do
		fgets(buf, 256, fpt);
	while (strcmp(buf, "[n]\n") != 0);
	for(m1=1;m1<= TUR[m].LCC;m1++)
		for (m2 = 1;m2 <=ncc;m2++)
			fscanf_s(fpt, "%lf", &n_z[m1][m2]);
	//input Q
	do
		fgets(buf, 256, fpt);
	while (strcmp(buf, "[Q]\n") != 0);
	for (m1 = 1;m1 <= TUR[m].LCC;m1++)
		for (m2 = 1;m2 <= ncc;m2++)
			fscanf_s(fpt, "%lf", &Q_z[m1][m2]);
	//input M
	do
		fgets(buf, 256, fpt);
	while (strcmp(buf, "[M]\n") != 0);
	for (m1 = 1;m1 <= TUR[m].LCC;m1++)
		for (m2 = 1;m2 <= ncc;m2++)
			fscanf_s(fpt, "%lf", &M_z[m1][m2]);
	//
	fopen_s(&fpc, "suterCC.xls", "w");
	fprintf(fpc, "ANG%t	WH%t	WB\n");
	//
	for (m1 = 1;m1 <= TUR[m].LCC;m1++)
		sutertrans(m,m1);
	//
	for (m2 = 0;m2 < TUR[m].LCC;m2++)
	{
		for (j = 0;j < TUR[m].NCC;j++)
			fprintf(fpc, "%lf	%lf	%lf\n", TUR[m].CC[m2].ANG[j], TUR[m].CC[m2].WH[j], TUR[m].CC[m2].WB[j]);
	}
	fclose(fpc);
}
void sutertrans(int m,int k)
{
	int m2;
	double ar[300] = { 0 }, ga[300] = { 0 }, beh[300] = { 0 }, beb[300] = { 0 };
	double A[300][300] = { 0 };
	double hd[300], hyh[300], hyb[300];

	
	rewind(fpt);
	for (j = 1; j <= ncc; j++)                        //Suter transformation
	{
		n_z[k][j] /= n01, Q_z[k][j] /= Q01, M_z[k][j] /= M01;  //Relative unit parameter    a,v,¦Â
		TUR[m].CC[k].ANG[j] = atan((Q_z[k][j]+ TUR[m].k1)/n_z[k][j]);	
		if (n_z[k][j] < 0)
			TUR[m].CC[k].ANG[j] += 3.1415926;
		TUR[m].CC[k].WH[j] = pow(TUR[m].guide[k]+ TUR[m].Cy,2) / (pow(Q_z[k][j], 2) + pow(n_z[k][j], 2)+ TUR[m].Ch);
		TUR[m].CC[k].WB[j] = (M_z[k][j]+ TUR[m].k2) * TUR[m].CC[k].WH[j];
	}
	//Sort curve
	for (m2 = 1;m2 < TUR[m].NCC;m2++)
		if (TUR[m].CC[k].ANG[m2] < TUR[m].CC[k].ANG[m2 + 1])
		{
			TUR[m].CC[k].ANG[0] = TUR[m].CC[k].ANG[m2+1];
			TUR[m].CC[k].WH[0] = TUR[m].CC[k].WH[m2 + 1];
			TUR[m].CC[k].WB[0] = TUR[m].CC[k].WB[m2 + 1];
			TUR[m].CC[k].ANG[m2+1] = TUR[m].CC[k].ANG[m2];
			TUR[m].CC[k].WH[m2 + 1] = TUR[m].CC[k].WH[m2];
			TUR[m].CC[k].WB[m2 + 1] = TUR[m].CC[k].WB[m2];
			for (j = m2 - 1;TUR[m].CC[k].ANG[0] > TUR[m].CC[k].ANG[j];j--)
			{
				TUR[m].CC[k].ANG[j + 1] = TUR[m].CC[k].ANG[j];
				TUR[m].CC[k].WH[j + 1] = TUR[m].CC[k].WH[j];
				TUR[m].CC[k].WB[j + 1] = TUR[m].CC[k].WB[j];
			}
			TUR[m].CC[k].ANG[j + 1] = TUR[m].CC[k].ANG[0];
			TUR[m].CC[k].WH[j + 1] = TUR[m].CC[k].WH[0];
			TUR[m].CC[k].WB[j + 1] = TUR[m].CC[k].WB[0];
		}
	TUR[m].CC[k].ANG[0] = 1.572797;
	TUR[m].CC[k].WH[0] = linear(TUR[m].CC[k].ANG[0], TUR[m].CC[k].ANG[1], TUR[m].CC[k].ANG[2], TUR[m].CC[k].WH[1], TUR[m].CC[k].WH[2]);
	TUR[m].CC[k].WB[0] = linear(TUR[m].CC[k].ANG[0], TUR[m].CC[k].ANG[1], TUR[m].CC[k].ANG[2], TUR[m].CC[k].WB[1], TUR[m].CC[k].WB[2]);
	TUR[m].CC[k].ANG[ncc+1] = 0.0;
	TUR[m].CC[k].WH[ncc + 1] = linear(TUR[m].CC[k].ANG[ncc + 1], TUR[m].CC[k].ANG[ncc-1], TUR[m].CC[k].ANG[ncc], TUR[m].CC[k].WH[ncc - 1], TUR[m].CC[k].WH[ncc]);
	TUR[m].CC[k].WB[ncc + 1] = linear(TUR[m].CC[k].ANG[ncc + 1], TUR[m].CC[k].ANG[ncc-1], TUR[m].CC[k].ANG[ncc], TUR[m].CC[k].WB[ncc - 1], TUR[m].CC[k].WB[ncc]);

//Spline
	for (j = 2; j <= ncc; j++)
	{
		TUR[m].CC[k].h[j] = TUR[m].CC[k].ANG[j] - TUR[m].CC[k].ANG[j - 1];
		hyh[j] = TUR[m].CC[k].WH[j] - TUR[m].CC[k].WH[j - 1];
		hyb[j] = TUR[m].CC[k].WB[j] - TUR[m].CC[k].WB[j - 1];
	}
	for (j = 2; j < ncc; j++)
	{
		hd[j] = TUR[m].CC[k].h[j] + TUR[m].CC[k].h[j + 1];
		ar[j] = TUR[m].CC[k].h[j + 1] / hd[j];
		ga[j] = 1 - ar[j];
		beh[j] = 6 / hd[j] * (hyh[j + 1] / TUR[m].CC[k].h[j + 1] - hyh[j] / TUR[m].CC[k].h[j]);
		beb[j] = 6 / hd[j] * (hyb[j + 1] / TUR[m].CC[k].h[j + 1] - hyb[j] / TUR[m].CC[k].h[j]);
	}
	for (j = 1; j <= ncc; j++)
	{
		if (j == 1)
		{
			A[1][1] = -TUR[m].CC[k].h[3];
			A[1][2] = TUR[m].CC[k].h[2] + TUR[m].CC[k].h[3];
			A[1][3] = -TUR[m].CC[k].h[2];
		}
		else if (j == ncc)
		{
			A[ncc][ncc] = -TUR[m].CC[k].h[ncc - 1];
			A[ncc][ncc - 1] = TUR[m].CC[k].h[ncc - 1] + TUR[m].CC[k].h[ncc];
			A[ncc][ncc - 2] = -TUR[m].CC[k].h[ncc];
		}
		else
		{
			A[j][j - 1] = ga[j];
			A[j][j] = 2;
			A[j][j + 1] = ar[j];
		}
	}
	GAUSS(A, TUR[m].CC[k].MH, beh, ncc);
	GAUSS(A, TUR[m].CC[k].MB, beb, ncc);
	return;
}