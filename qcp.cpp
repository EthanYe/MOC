#include<stdio.h>
#include<math.h>
#include"parameters.h"

double QCP(int I, int J)
{
	double  ROUGH;
	if (Q[I][J - 1] >= 0)
		ROUGH = 8 * 9.81*L[I].PipeLoss[0] * L[I].AJD[J - 1] * L[I].AJD[J - 1] * (L[I].AJD[J - 1] / L[I].SJD[J - 1]) / L[I].DX[J] / L[I].NN;
	else
		ROUGH = 8 * 9.81*L[I].PipeLoss[1] * L[I].AJD[J - 1] * L[I].AJD[J - 1] * (L[I].AJD[J - 1] / L[I].SJD[J - 1]) / L[I].DX[J] / L[I].NN;

	double c = L[I].a / 9.81;
	double c1 = 0.0;
	double c2 = TimeStep*ROUGH*fabs(Q[I][J - 1])*L[I].SJD[J] / 8 / L[I].AJD[J] / L[I].AJD[J] / L[I].AJD[J];
	double c3 = 0.0;
	double qcp = 1.0 / ((c - c3) / L[I].AJD[J] + c*(c1 + c2))*(Q[I][J - 1] * (c + c3) / L[I].AJD[J] + H[I][J - 1]);
	return qcp;
}

double CQP(int I, int J)
{
	double ROUGH;
	if (Q[I][J - 1] >= 0)
		ROUGH = 8 * 9.81*L[I].PipeLoss[0] * L[I].AJD[J - 1] * L[I].AJD[J - 1] * (L[I].AJD[J - 1] / L[I].SJD[J - 1]) / L[I].DX[J] / L[I].NN;
	else
		ROUGH = 8 * 9.81*L[I].PipeLoss[1] * L[I].AJD[J - 1] * L[I].AJD[J - 1] * (L[I].AJD[J - 1] / L[I].SJD[J - 1]) / L[I].DX[J] / L[I].NN;

	double c = L[I].a / 9.81;
	double c1 = 0.0;
	double c2 = TimeStep*ROUGH*fabs(Q[I][J - 1])*L[I].SJD[J] / 8 / L[I].AJD[J] / L[I].AJD[J] / L[I].AJD[J];
	double c3 = 0.0;
	double cqp = 1.0 / ((c - c3) / L[I].AJD[J] + c*(c1 + c2));
	return cqp;
}

double QCM(int I, int J)
{
	double ROUGH;
	if (Q[I][J + 1] >= 0)
		ROUGH = 8 * 9.81*L[I].PipeLoss[0] * L[I].AJD[J] * L[I].AJD[J] * (L[I].AJD[J] / L[I].SJD[J]) / L[I].DX[J] / L[I].NN;
	else
		ROUGH = 8 * 9.81*L[I].PipeLoss[1] * L[I].AJD[J] * L[I].AJD[J] * (L[I].AJD[J] / L[I].SJD[J]) / L[I].DX[J] / L[I].NN;

	double c = L[I].a / 9.81;
	double c3 = 0.0;
	double c4 = 0.0;
	double c5 = TimeStep*ROUGH*fabs(Q[I][J + 1])*L[I].SJD[J] / 8 / L[I].AJD[J] / L[I].AJD[J] / L[I].AJD[J];
	double qcm = 1.0 / ((c + c3) / L[I].AJD[J] + c*(c4 + c5))*(Q[I][J + 1] * (c - c3) / L[I].AJD[J] - H[I][J + 1]);
	return qcm;
}

double CQM(int I, int J)
{
	double ROUGH;
	if (Q[I][J + 1] >= 0)
		ROUGH = 8 * 9.81*L[I].PipeLoss[0] * L[I].AJD[J] * L[I].AJD[J] * (L[I].AJD[J] / L[I].SJD[J]) / L[I].DX[J] / L[I].NN;
	else
		ROUGH = 8 * 9.81*L[I].PipeLoss[1] * L[I].AJD[J] * L[I].AJD[J] * (L[I].AJD[J] / L[I].SJD[J]) / L[I].DX[J] / L[I].NN;

	double c = L[I].a / 9.81;
	double c3 = 0.0;
	double c4 = 0.0;
	double c5 = TimeStep*ROUGH*fabs(Q[I][J + 1])*L[I].SJD[J] / 8 / L[I].AJD[J] / L[I].AJD[J] / L[I].AJD[J];
	double cqm = 1.0 / ((c + c3) / L[I].AJD[J] + c*(c4 + c5));
	return cqm;
}