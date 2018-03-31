#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"input.h"
#include"parameters.h"
#include"turbinefun.h"
#include"linear_interpolation.h"
#include"serch.h"


double Rough(int I, int J,int K);             // λ
double fun(double ,double, double , double,int,double);
double ang[30][300], wh[30][300];
//FILE *fpt;
double Cond(double(*)[300], int);
double qiuni(double (*)[300], int );
void SteadyFlowCalculation()
{
	////////Element method
	double A[50][300] = { 0 };                    //Characteristic matrix
	double B[50] = { 0 };                          //Flow matrix
	double E[50] = { 0 };                        //Head increment
	double QK[30][30][2] = { 0 };                //Flow for iteration  [0] represent previous iteration flow or head
	double HK[50][2] = { 0 };                    //Head for iteration
	int IMAX = 50;                               //Maximum iterations
	double v = 2.5;                           //Initial velocity
	double S[50][50][2] = { 0 };                    //Coefficient of resistance of pipe
	double SY[50][50][2] = { 0 };                    
	double SJ[50][50][2] = { 0 };                   
	double p[50][50] = { 0 }, q[50][50] = { 0 }, r[50][50] = { 0 };
	int TotalN;                                 //Total node
	int k1, k2,k3;                                //For loop
	double Mik;                                //For elimination
	double tx;                                 //For backward substitution
	int linenumber;                            //Choose line
	double MAXA = 0;
	double exchange;                           //For exchange line in characteristic matrix
	double e;                                  //Accuracy
	int n;                                     //EN[n]
	double pa;
	int tn;                                    //Node number of turbine 
	double x0,x1, x2,C;
	double guide0;                                 //guide vane degree at 0s
	int m1, m2;
	

	double A0, A1, x, HT[50], VQ, F, FQ; 
	////////////////
	TotalN = TotalNode + TotalENode + 1;
	//fopen_s(&fpt, "intersuter.xls", "w");
	//fprintf(fpt, "ANG%t	WH%t	WB\n");

//	fopen_s(&fp1[40],"Result of steady flow calculation.xls","w");
	if (EN[1].K > 3 && EN[1].K < 7)
	{
		for (i = 1; i <= N; i++)
		{
			H[i][0] = HR0;
			Q[i][0] = 0;
			for (k2 = 1; k2 <= L[i].NN; k2++)
			{
				Q[i][k2] = Q[i][0];
				H[i][k2] = H[i][k2 - 1] - Rough(i, k2 - 1,0)*8.0 *L[i].m_Length / 3.1415926 / 3.1415926 / pow(L[i].D, 5.0) / g / L[i].NN * Q[i][k2] * Q[i][k2];
			}
		}
		fprintf(fp1[40], "%lf\t %lf \t%lf \n", T, Q[N][L[N].NN], H[N][L[N].NN]);
		return; 
	}

	//Calculate coefficient of resistance 
	for (i = 1; i <= N; i++)
	{
		for (j = 0;j < 2;j++)
		{
			SY[L[i].StartNode][L[i].EndNode][j] = 8.0 * Rough(i, 1, j)*L[i].m_Length / 3.1415926 / 3.1415926 / pow(L[i].D, 5.0) / g;
			SJ[L[i].StartNode][L[i].EndNode][j] = L[i].LLCoefficient[j] / 2 / g / L[i].AJD[0] / L[i].AJD[0];
			if (JN[L[i].EndNode - 1].type == 2)
				SJ[L[i].StartNode][L[i].EndNode][j] += 1.0 / 2.0 / g / CDA[0][L[i].EndNode - 1] / CDA[0][L[i].EndNode - 1];
	//		else if (JN[L[i].EndNode - 1].type == 4)
	//			SJ[L[i].StartNode][L[i].EndNode][j] += -1.0 / 2.0 / g / L[i+1].AJD[0] / L[i+1].AJD[0] + 1.0 / 2.0 / g / (3.14159*0.15*0.15 / 4) / (3.14159*0.15*0.15 / 4);
			S[L[i].StartNode][L[i].EndNode][j] = SJ[L[i].StartNode][L[i].EndNode][j] + SY[L[i].StartNode][L[i].EndNode][j];		
			S[L[i].EndNode][L[i].StartNode][j] = S[L[i].StartNode][L[i].EndNode][j];
		}
	}

	//Initial flow
	for (i = 1; i <= N; i++)
	{
		QK[L[i].StartNode][L[i].EndNode][0] =0.0491;
		QK[L[i].EndNode][L[i].StartNode][0] = -QK[L[i].StartNode][L[i].EndNode][0];
	}

	//Initial head
 	HK[1][0] = HR0 - QK[1][2][0] * QK[1][2][0] / 2 / g / L[1].AJD[0] / L[1].AJD[0];
	//HK[1][0] = HR0;
	for (j = 1; j <= N; j++)
	{
		if (HK[L[j].StartNode][0] != 0 && HK[L[j].EndNode][0] == 0)
		{
			if (L[j].EndNode > TotalNode + 1)
			{
				n = L[j].EndNode - TotalNode - 1;
				if (EN[n].K == 2)                   //Downstream reservior
					HK[L[j].EndNode][0] = HE[n];
				else
					HK[L[j].EndNode][0] = HK[L[j].StartNode][0] - S[L[j].StartNode][L[j].EndNode][0] * QK[L[j].StartNode][L[j].EndNode][0] * QK[L[j].StartNode][L[j].EndNode][0];
			}
			else HK[L[j].EndNode][0] = HK[L[j].StartNode][0] - S[L[j].StartNode][L[j].EndNode][0] * QK[L[j].StartNode][L[j].EndNode][0] * QK[L[j].StartNode][L[j].EndNode][0];
		}
		if (JN[L[j].EndNode - 1].type == 4)
		{
			tn = L[j].EndNode - 1;
			gui(tn,&guide0);
			serch(guide0,&m1,&m2, TUR[tn].guide, TUR[tn].LCC);
			for (i = 1;i <= TUR[tn].NCC;i++)
			{
				ang[tn][i]=linear(guide0, TUR[tn].guide[m1], TUR[tn].guide[m2], TUR[tn].CC[m1].ANG[i], TUR[tn].CC[m2].ANG[i]);
				wh[tn][i] = linear(guide0, TUR[tn].guide[m1], TUR[tn].guide[m2], TUR[tn].CC[m1].WH[i], TUR[tn].CC[m2].WH[i]);
			//	fprintf(fpt, "%lf%t	%lf%t	%lf\n", ang[tn][i], wh[tn][i] , wb[tn][i]);
			}
		}
	}


	/////////Iteration
	for (k1 = 1; k1 <= IMAX; k1++)
	{
		for (i = 1; i <= N; i++)
		{
			if (QK[L[i].StartNode][L[i].EndNode][0] > 0)
				k2 = 0;
			else k2 = 1;
			if (QK[L[i].StartNode][L[i].EndNode][0] == 0)
			{
				p[L[i].StartNode][L[i].EndNode] = 1.0;
				p[L[i].EndNode][L[i].StartNode] = p[L[i].StartNode][L[i].EndNode];
				q[L[i].StartNode][L[i].EndNode] = -p[L[i].StartNode][L[i].EndNode];
				q[L[i].EndNode][L[i].StartNode] = q[L[i].StartNode][L[i].EndNode];
				r[L[i].StartNode][L[i].EndNode] = p[L[i].StartNode][L[i].EndNode] * (HK[L[i].StartNode][0] - HK[L[i].EndNode][0] - S[L[i].StartNode][L[i].EndNode][0] * fabs(QK[L[i].StartNode][L[i].EndNode][0])*QK[L[i].StartNode][L[i].EndNode][0]);
				r[L[i].EndNode][L[i].StartNode] = -r[L[i].StartNode][L[i].EndNode];
			}
			else if (JN[L[i].EndNode - 1].type == 4)
			{
				
				tn = L[i].EndNode - 1;				
				VQ = QK[L[i].StartNode][L[i].EndNode][0] / TUR[tn].Q0;
				x =  atan((VQ + TUR[tn].k1) / 1);
				//x = atan2(VQ, 1);
				for (j = 2; j <= TUR[tn].NCC; j++)
				{
					if ((x >= ang[tn][j] && x <= ang[tn][j-1] )|| (x <= ang[tn][j] && x >= ang[tn][j - 1]))
					{
						A0 = wh[tn][j - 1], A1 = (wh[tn][j] - wh[tn][j - 1]) / (ang[tn][j] - ang[tn][j - 1]);
						double	WH = A0 + A1*fabs(x - ang[tn][j]);
						HT[tn] = TUR[tn].H0*(1 + pow(VQ, 2))*WH/(pow(guide0+ TUR[tn].Cy,2)-WH*TUR[tn].Ch);
						break;
					}
				}
				if (j == TUR[tn].NCC + 1)
					exit(0);
				F = HK[L[i].StartNode][0] - HK[L[i].EndNode][0] - fabs(QK[L[i].StartNode][L[i].EndNode][0]) / QK[L[i].StartNode][L[i].EndNode][0] * HT[tn] - S[L[i].StartNode][L[i].EndNode][k2] * fabs(QK[L[i].StartNode][L[i].EndNode][0])*QK[L[i].StartNode][L[i].EndNode][0];
				FQ = -TUR[tn].H0*(2 * (A0 + A1*x)*VQ / TUR[tn].Q0 + A1/TUR[tn].Q0) - 2 * S[L[i].StartNode][L[i].EndNode][k2] * fabs(QK[L[i].StartNode][L[i].EndNode][0]);
				p[L[i].StartNode][L[i].EndNode] = -1 / FQ;
				p[L[i].EndNode][L[i].StartNode] = p[L[i].StartNode][L[i].EndNode];
				q[L[i].StartNode][L[i].EndNode] = -p[L[i].StartNode][L[i].EndNode];
				q[L[i].EndNode][L[i].StartNode] = q[L[i].StartNode][L[i].EndNode];
				r[L[i].StartNode][L[i].EndNode] = -F / FQ;
				r[L[i].EndNode][L[i].StartNode] = -r[L[i].StartNode][L[i].EndNode];
			}
			else
			{
				p[L[i].StartNode][L[i].EndNode] = 1 / (2 * S[L[i].StartNode][L[i].EndNode][k2] * fabs(QK[L[i].StartNode][L[i].EndNode][0]));
				p[L[i].EndNode][L[i].StartNode] = p[L[i].StartNode][L[i].EndNode];
				q[L[i].StartNode][L[i].EndNode] = -p[L[i].StartNode][L[i].EndNode];
				q[L[i].EndNode][L[i].StartNode] = q[L[i].StartNode][L[i].EndNode];
				r[L[i].StartNode][L[i].EndNode] = p[L[i].StartNode][L[i].EndNode] * (HK[L[i].StartNode][0] - HK[L[i].EndNode][0] - S[L[i].StartNode][L[i].EndNode][k2] * fabs(QK[L[i].StartNode][L[i].EndNode][0])*QK[L[i].StartNode][L[i].EndNode][0]);
				r[L[i].EndNode][L[i].StartNode] = -r[L[i].StartNode][L[i].EndNode];
			}
		}
		////Calculate characteristic matrix A
		for (i = 1; i <= TotalN; i++)
		{
			for (j = 1; j <= TotalN; j++)
			{
				if (i == 1)                             //a[1][j] at upstream node
				{
					if (j == 1)
						for (k2 = 1; k2 <= TotalN; k2++)
							if (k2 != i)
								A[i][j] = A[i][j] + p[i][k2];
					else if (j > 1)
								A[i][j] = 0;
				}
				else if (i > 1 && i <= TotalNode + 1)
				{
					if (i == j)
					{
						for (k2 = 1; k2 <= TotalN; k2++)
							if (k2 != i)
								A[i][j] = A[i][j] + p[i][k2];
					}
					else A[i][j] = q[i][j];
				}
				else if (i > TotalNode + 1)      //End node
				{
					n = i - TotalNode - 1;
					if (EN[n].K == 1)            //Closing valve
					{
						if (QK[i][L[EN[n].StartPipe].StartNode] < 0)
							k2 = 0;
						else k2 = 1;
						pa = S[i][L[EN[n].StartPipe].StartNode][k2] * CDA[1][n] * CDA[1][n] * 2 * g;
						if (i == j)
						{
							A[i][j] = 1.0 + pa;
						}
						else A[i][L[EN[n].StartPipe].StartNode] = -1.0;
					}
					if (EN[n].K == 2)                      //Downstream reservior
					{
						if (i == j)
						{
							for (k2 = 1; k2 <= TotalN; k2++)
								if (k2 != i)
									A[i][j] = A[i][j] + p[i][k2];
						}
						else A[i][j] = 0;
					}
					else if (EN[n].K == 3)           //Non reflecting boundary
					{
						if (i == j)
						{
							for (k2 = 1; k2 <= TotalN; k2++)
								if (k2 != i)
									A[i][j] = A[i][j] + p[i][k2];
						}
						else A[i][j] = q[i][j];
					}
				}
			}
		}
		double kkk=Cond(A, TotalN);
		////Calculate extended matrix B
		for (i = 1; i <= TotalN; i++)
		{
			if (i == 1)
			{
				for (k2 = 1; k2 <= TotalN; k2++)
					//B[i] -= QK[1][k2][0]/2.0/g/ L[1].AJD[0] / L[1].AJD[0]+HK[1][0]-HR0+QK[1][k2][0]* QK[1][k2][0]/2/g/ L[1].AJD[0] / L[1].AJD[0];
					B[i] = 0;
			}
			else if (i > 1 && i <= TotalNode + 1)
			{
				for (k2 = 1; k2 <= TotalN; k2++)
					B[i] = B[i] - QK[i][k2][0] - r[i][k2];
			}
			else if (i > TotalNode + 1)
			{
				n = i - TotalNode - 1;
				if (EN[n].K == 1)                             //Closing valve
				{
					B[i] = HK[L[EN[n].StartPipe].StartNode][0] - (1.0 + pa)*HK[i][0];
				}
				if (EN[n].K == 2)                            //Downstream reservior
					B[i] = 0;
				else if (EN[n].K == 3)                        //Non reflecting boundary
				{
					for (k2 = 1; k2 <= TotalN; k2++)
						B[i] = B[i] - QK[i][k2][0] - r[i][k2];
				}
			}
		}
		//Elimination
		for (int k = 1; k < TotalN; k++)
		{
			for (i = k; i <= TotalN; i++)      //Choose lin  e number
			{
				if (fabs(A[i][k]) > MAXA)
				{
					MAXA = fabs(A[i][k]);
					linenumber = i;
				}
				else linenumber = k;
			} //Choose line number
			for (i = 1; i <= TotalN; i++)       //Exchange line in characteristic matrix
			{
				exchange = A[k][i]; A[k][i] = A[linenumber][i]; A[linenumber][i] = exchange;
			}
			exchange = B[k]; B[k] = B[linenumber]; B[linenumber] = exchange;
			for (i = k + 1; i <= TotalN; i++)
			{
				if (A[k][k] == 0)
					continue;
				Mik = A[i][k] / A[k][k];
				for (j = k; j <= TotalN; j++)
					A[i][j] = A[i][j] - Mik*A[k][j];
				B[i] = B[i] - Mik*B[k];
			}
		}
		//Backward substitution
		for (i = TotalN; i > 0; i--) 
		{
			if (A[i][i] == 0)
				break;
			tx = 0;
			for (j = i + 1; j <= TotalN; j++)
				tx = tx + A[i][j] * E[j];
			E[i] = (B[i] - tx) / A[i][i];
		}
		for (i = 1; i <= TotalN; i++)
		{
			HK[i][1] = HK[i][0] + E[i];
		}
		for (i = 1; i <= N; i++)
		{
			if (JN[L[i].EndNode - 1].type == 4)
			{
				tn = L[i].EndNode - 1;
				C = HK[L[i].StartNode][1] - HK[L[i].EndNode][1];
				x1 = QK[L[i].StartNode][L[i].EndNode][0];
				x0 = x1 / 1.01;
				for (k2 = 1; k2 <= IMAX; k2++)
				{
					if (x1 > 0)
						k3 = 0;
					else k3 = 1;
					x2 = fun(x1,x0, C, S[L[i].StartNode][L[i].EndNode][k3], i,guide0);
					e = fabs(x2-x1);
					if (e < 0.000001)
						break;
					else
						x0 = x1,x1 = x2;
				}
				QK[L[i].StartNode][L[i].EndNode][1] = x2;
			}
			else if (HK[L[i].StartNode][1] > HK[L[i].EndNode][1])
			{
				QK[L[i].StartNode][L[i].EndNode][1] = sqrt((HK[L[i].StartNode][1] - HK[L[i].EndNode][1]) / S[L[i].StartNode][L[i].EndNode][0]);
			}
			else
				QK[L[i].StartNode][L[i].EndNode][1] = -sqrt((HK[L[i].EndNode][1] - HK[L[i].StartNode][1]) / S[L[i].StartNode][L[i].EndNode][1]);
			QK[L[i].EndNode][L[i].StartNode][1] = -QK[L[i].StartNode][L[i].EndNode][1];
		}
		j = 0;
		for (i = 1; i <= TotalN; i++)
			{
				e = fabs(E[i] / HK[i][0]);
				if (e < 0.001)
					j++;
			}
		if (j == TotalN)
			break;
		else
			for (i = 1; i <= TotalN; i++)
				{
					B[i] = 0;
					E[i] = 0;
					HK[i][0] = HK[i][1];
					for (j = 1; j <= TotalN; j++)
					{
						A[i][j] = 0;
						p[i][j] = 0;
						q[i][j] = 0;
						r[i][j] = 0;
						QK[i][j][0] = QK[i][j][1];
					}
				}
		}
		/////////////

	if (k1 == IMAX + 1)
		printf("hengdingliu convergence failed\n");
	else 
		printf("hengdingliu converged successfully\n");

	for (i = 1; i <= N; i++)
	{
		//H[i][0] = HK[L[i].StartNode][1]- QK[L[i].StartNode][L[i].EndNode][1]* QK[L[i].StartNode][L[i].EndNode][1]/2/g/ L[i].AJD[0]/ L[i].AJD[0];
		H[i][0] = HK[L[i].StartNode][1];
		Q[i][0] = QK[L[i].StartNode][L[i].EndNode][1];
		for (k2 = 1; k2 <= L[i].NN; k2++)
		{
			Q[i][k2] = Q[i][0];
			if (Q[i][k2] > 0)
				k3 = 0;
			else k3 = 1;
			H[i][k2] = H[i][k2 - 1] - S[L[i].StartNode][L[i].EndNode][k3] * fabs(Q[i][k2]) * Q[i][k2]/ L[i].NN;
		}
	}
//		fprintf(fp1[40], "%lf\t %lf \t%lf \n", T, Q[N][L[N].NN], H[N][L[N].NN]);
	for (i = 1; i <= N; i++)
		fprintf(fp1[FileNumber + 1], "L%d\t%lf\t %lf \t%lf \n", i, Q[i][0],H[i][0]-HE[1],H[i][L[i].NN]-HE[1]);
}
 
double Rough(int I, int J,int K)
{
	double  R;
	if (K= 0)
		R = 8 * 9.81*L[I].PipeLoss[0] * L[I].AJD[J] * L[I].AJD[J] * (L[I].AJD[J] / L[I].SJD[J]) / L[I].DX[J] / L[I].NN;
	else
		R = 8 * 9.81*L[I].PipeLoss[1] * L[I].AJD[J] * L[I].AJD[J] * (L[I].AJD[J] / L[I].SJD[J]) / L[I].DX[J] / L[I].NN;
	return R;
}
double fun(double x,double x0,double c,double s,int m0,double guide0)
{
	double y;
	double v;
	double X;
	double A0, A1;
	double HT,HT0;
	int m1;
	m1 = L[m0].EndNode - 1;
	v = x/ TUR[m1].Q0;
	X = atan(v+TUR[m1].k1 / 1);
	//X = atan2(v, 1);
	for (j = 2; j <= TUR[m1].NCC; j++)
	{
		if ((X >= ang[m1][j] && X <= ang[m1][j - 1]) || (X <= ang[m1][j] && X >= ang[m1][j - 1]))
		{
			A0 = wh[m1][j - 1], A1 = (wh[m1][j] - wh[m1][j - 1]) / (ang[m1][j] - ang[m1][j - 1]);
			double	WH = A0 + A1*fabs(X - ang[m1][j]);
			HT = TUR[m1].H0*(1 + pow(v, 2))*WH / (pow(guide0 + TUR[m1].Cy, 2) - WH*TUR[m1].Ch);
 			break;
		}
	}
	v = x0 / TUR[m1].Q0;
	X = atan((v+ TUR[m1].k1)/ 1);
	//X = atan2(v, 1);
	for (j = 2; j <= TUR[m1].NCC; j++)
	{
		if ((X >= ang[m1][j] && X <= ang[m1][j - 1]) || (X <= ang[m1][j] && X >= ang[m1][j - 1]))
		{
			A0 = wh[m1][j - 1], A1 = (wh[m1][j] - wh[m1][j - 1]) / (ang[m1][j] - ang[m1][j - 1]);
			double	WH = A0 + A1*fabs(X - ang[m1][j]);
			HT0 = TUR[m1].H0*(1 + pow(v, 2))*WH / (pow(guide0 + TUR[m1].Cy, 2) - WH*TUR[m1].Ch);
			break;
		}
	}
 	y = x - (c - HT - s*x*fabs(x))*(x-x0)/((c - HT - s*x*fabs(x))-(c-HT0-s*x0*fabs(x0)));
	//-TUR[tn].H0*(2 * (A0 + A1*x)*VQ / TUR[tn].Q0 + A1*TUR[tn].Q0) - 2 * S[L[i].StartNode][L[i].EndNode][k2] * fabs(QK[L[i].StartNode][L[i].EndNode][0]);
	return y;
}
double Cond(double (*x)[300],int n)
{
	double c1=0, c2,cond;

	double flag;
	int i,j;
	/*Calculate c1*/
	for (i = 1;i <= n;i++)
	{
		flag = 0;
		for (j = 1;j <= n;j++)
		{
			flag += *(*(x + i) + j);
			if (c1 < flag)
				c1 = flag;
		}
	}
	/*Calculate c2*/
	c2=qiuni(x, n);
	cond = c1*c2;
	return cond;
}
double qiuni(double(*a)[300], int N)
{
	double b[30][30] = { 0 }, c[30][30] = {0}, t;
	int i, j, m;
	double flag,c2=0;
	//增广矩阵（A|E）存入二维数组b中
	for (i = 1;i<=N;i++)
		for (j = 1;j<=N;j++)
			b[i][j] = *(*(a+i)+j);
	//for (i = 1;i<=N;i++)
	//	for (j = N;j<2 * N;j++)
	//		b[i][j] = 0;
	for (i = 1;i<=N;i++)
		b[i][N + i] = 1;
	for (m = 1;m<=N;m++)          //对每行进行处理。
	{
		t = b[m][m];                  //预存b[m][m]。
		i = m;
		while (b[m][m] == 0)
		{
			b[m][m] = b[i + 1][m];
			i++;
		}
		if (i>m)
		{
			b[i][m] = t;                  //实现交换。
										  //交换其它各列相应位置的元素
			for (j = 1;j<m;j++)
			{
				t = b[m][j];
				b[m][j] = b[i][j];
				b[i][j] = t;
			}
			for (j = m + 1;j<=2 * N;j++)
			{
				t = b[m][j];
				b[m][j] = b[i][j];
				b[i][j] = t;
			}

		}
		for (i = m + 1;i<=N;i++)
			for (j = 2 * N ;j >= m;j--)
				b[i][j] -= b[i][m] * b[m][j] / b[m][m]; //m=0时，将第一行的-b[i][0]/b[0][0]倍加到以下各行。这样以下每行第一个元素b[i][0]就为0。
		for (j = 2 * N ;j >= m;j--)
			b[m][j] /= b[m][m];   //对第m行作行变换，同除以b[m][m]，使b[m][m]为1。
	}
	m  = N ;
	while (m>0)
	{
		for (i = 1;i<m;i++)
			for (j = 2 * N;j >= m;j--)           //千万注意，此处j必须递减，否则b[i][m]先变为0，后面的计算就无效！
				b[i][j] -= b[i][m] * b[m][j];
		m--;
	}

	for (i = 1;i<=N;i++)                         //将逆矩阵存入二维数组c中。
		for (j = 1;j<=N;j++)
			c[i][j] = b[i][N + j];
	for (i = 1;i <= N;i++)
	{
		flag = 0;
		for (j = 1;j <= N;j++)
		{
			flag += c[i][j];
			if (c2 < flag)
				c2 = flag;
		}
	}
	return c2;
}