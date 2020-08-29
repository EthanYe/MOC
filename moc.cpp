#include<stdio.h>
#include<math.h>
#include"parameters.h"
#include"boundary_condition.h"
#include"qcp.h"
void moc()
{
	double a1, A2;
	A2 = 1.0 / 2 / g / L[1].AJD[0] / L[1].AJD[0];

	//Calculation of internal surface  
	///////////////////////////////////////	
	for (i = 1; i <= N; i++)                           //Calculation of internal nodes of pipes
	{//12345678
		for (int J0 = 1; J0 <= L[i].NN - 1; J0++)
		{
			int I1 = i;
			HP[I1][J0] = (QCP(I1, J0) - QCM(I1, J0)) / (CQP(I1, J0) + CQM(I1, J0));
			QP[I1][J0] = QCP(I1, J0) - CQP(I1, J0)*HP[I1][J0];
		}
	}

	//Calculation of boundary of upstream reservoir
	a1 = CQM(1, 0)*A2;
	QP[1][0] = (sqrt(1 + 4 * a1*(QCM(1, 0) + CQM(1, 0)*HR0)) - 1) / 2 / a1;
	HP[1][0] = HR0 - QP[1][0] * QP[1][0] * A2;

	////////////////////////////////
	//Calculation of pipe interface
	for (i = 1; i < N; i++)
	{
		struct Node NODE = JN[i];
		NodeBoundary(&NODE);
	}

	//Calculation of end boundary 
	for (i = 1; i <= TotalENode; i++)
	{
		struct ENode EndNODE = EN[i];
		EndBoundary(&EndNODE);
	}

	//Calculation of end boundary end

	for (i = 1; i <= N; i++)
	{
		for (int m0 = 0; m0 <= L[i].NN; m0++)
		{
			Q[i][m0] = QP[i][m0];
			H[i][m0] = HP[i][m0];
		}
	}
}
