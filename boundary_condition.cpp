#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"parameters.h"
#include"gauss.h"
#include"linear_interpolation.h"
#include"turbinefun.h"
#include"qcp.h"





 // NodeBoundary
void NodeBoundary(struct Node *JP)
{
	///Junction
	static double a0[30][2] = {0}, ap[30],a1[30];
	double  b0, c0, ap1[30], ap2[30];
	/////////TURBINE
	static double n_T[50][4], Q_T[50][4], h_T[50][4];                //Real parameter for iterration ([50] represent the node number the turbine )
	static double  H_T= TUR[i].H0;                                    //Real parameter for iterration
	static int densi=1;
	double n_R, Q_R, M_R, H_R;                          //Real relative parameter
	double  aang, WHH,  WBB;                             //For suter transformation

	int IMAX = 500;                                       //Max iteration
	int k0;                                           //loop
	int I, J;
	double TAU;

	double Qi, ni;                                       //To compare with estimated q,n
	/////////////////////


	double guidet;
	double pw, ps;
	static double warea[30];
	static double sarea[30];
	static double dd1 = 0.0, dd2 = 0.0;
	////SURGE TANK
	double SC, SQ;


	I = JP->StartPipe[1];
	J = JP->EndPipe[1];
	if (JP->type == 1)               //Junction
	{
		if (T == TimeStep)
		{
			if (ap[I] == 0)
				ap[I] = 1.0 / 2.0 / g / L[I].AJD[0] / L[I].AJD[0];
			if(ap[J]==0)
				ap[J] = 1.0 / 2.0 / g / L[J].AJD[0] / L[J].AJD[0];
			for (k0 = 0;k0 < 2;k0++)
			{
				a0[I][k0] = ap[J] - ap[I] + (1.0 - 2.0*k0)*L[I].LLCoefficient[k0] * ap[I];
			}
		}
		b0 = 1 / CQP(I, L[I].NN) + 1 / CQM(J, 0);
		c0 = -QCP(I, L[I].NN) / CQP(I, L[I].NN) - QCM(J, 0) / CQM(J, 0);
		if (Q[I][L[I].NN] >= 0)
			k0 = 0;
		else k0 = 1;
		if (a0[I][k0] < 0.000001)
			QP[I][L[I].NN] = -c0 / b0;
		else
			QP[I][L[I].NN] = (-b0 + sqrt(b0*b0 - 4.0*a0[I][k0]*c0)) / (2.0*a0[I][k0]);
		QP[J][0] = QP[I][L[I].NN];
		HP[I][L[I].NN] = (QCP(I, L[I].NN)- QP[I][L[I].NN]) / CQP(I, L[I].NN);
		HP[J][0] = (QP[J][0]  - QCM(J, 0)) / CQM(J, 0);
	}	
	else if (JP->type == 2)          //in-line valve
	{
		if (T <= TC[0][i])
		{
			////QPA=QPB, QPA=TAU*CDA*sqrt(2*g(HA-HB)), QPA+, QPB-
			TAU = pow((1.0 - T / TC[0][i]), EM[0][i]);
			double A = 1.0;
			A = 1.0 / 2.0 / g / TAU / TAU / CDA[0][i] / CDA[0][i];
			double b = 1 / CQP(I, L[I].NN) + 1 / CQM(J, 0);
			double c = -QCP(I, L[I].NN) / CQP(I, L[I].NN) - QCM(J, 0) / CQM(J, 0);
			QP[I][L[I].NN] = (-b + sqrt(b*b - 4 * A*c)) / 2 / A;
			HP[I][L[I].NN] = (QCP(I, L[I].NN) - QP[I][L[I].NN]) / CQP(I, L[I].NN);
		}
		else
		{
			TAU = 0.0;
			QP[I][L[I].NN] = 0;
			HP[I][L[I].NN] = QCP(I, L[I].NN) / CQP(I, L[I].NN);
		}
		QP[J][0] = QP[I][L[I].NN];
		HP[J][0] = (QP[J][0] - QCM(J, 0)) / CQM(J, 0);
	}
	else if (JP->type == 3)                     //Bifurcated pipe
	{
		double l = 0.0,m=0.0,B=0.0;
		double hp;
//		int n = JP->StartPipeNumber + JP->EndPipeNumber;
		for (j = 1; j <= JP->StartPipeNumber; j++)
		{
			int sp = JP->StartPipe[j];
			l = l + QCP(sp, L[sp].NN);
			B = B +  CQP(sp, L[sp].NN);
		}
		for (j = 1; j <= JP->EndPipeNumber; j++)
		{
			int ep = JP->EndPipe[j];
			m = m - QCM(ep, 0);
			B = B + CQM(ep, 0);
		}
		hp = (l + m) / B;
		for (j = 1; j <= JP->StartPipeNumber; j++)
		{
			int sp = JP->StartPipe[j];
			HP[sp][L[sp].NN] = hp;
			QP[sp][L[sp].NN]= QCP(sp, L[sp].NN) - CQP(sp, L[sp].NN)*HP[sp][L[sp].NN];
		}
		for (j = 1; j <= JP->EndPipeNumber; j++)
		{
			int ep = JP->EndPipe[j];
			HP[ep][0] = hp;
			QP[ep][0] = QCM(ep, 0) + CQM(ep,0)*HP[ep][0];
		}
	}
	else if (JP->type == 4)                //Turbine boundary
	{

		

		///Initialize rated parameters
		if (T == TimeStep)
		{
			double area;
			for (j = 0; j < 4; j++)
				n_T[i][j] = TUR[i].n0, Q_T[i][j] =Q[i][1],h_T[i][j]=H[i][1];
			area = 3.14159*TUR[i].D1*TUR[i].D1 / 4;
			ap1[i]= 1.0 / 2.0 / g /area/area ;
			area = 3.14159*TUR[i].D2*TUR[i].D2 / 4;
			ap2[i] = 1.0 / 2.0 / g / area / area;
			a1[i] = ap2[i] - ap1[i];
			warea[i] = 3.14159*TUR[i].D1*0.02444 * 2 * g*3.14159*TUR[i].D1*0.02444;
			sarea[i] = 3.14159*0.15*0.15/4 * 2 * g* 3.14159*0.15*0.15 / 4;
			pw = H[I][L[I].NN] - TUR[i].Tz - Q[I][L[I].NN]*Q[I][L[I].NN] / warea[i];
			ps = H[J][0] - TUR[i].Tz - Q[I][L[I].NN]*Q[I][L[I].NN] / sarea[i];
		}

		//Caculate gui
		gui(i,&guidet);
		//////Estimate speed and flow
		n_T[i][3] = 3 * n_T[i][2] - 3 * n_T[i][1] + n_T[i][0];
		Q_T[i][3] = 3 * Q_T[i][2] - 3 * Q_T[i][1] + Q_T[i][0];
		h_T[i][3] = 3 * h_T[i][2] - 3 * h_T[i][1] + h_T[i][0];

		/////////Iteration
		Iteration:
		for ( k0 = 1; k0 <= IMAX; k0++)
		{
			n_R = n_T[i][3] / TUR[i].n0, Q_R = Q_T[i][3] / TUR[i].Q0;       //Relative n,Q
			aang = atan((Q_R + TUR[i].k1*sqrt (H_T/ TUR[i].H0)) / n_R);                     //Horizontal ordinate
			if (n_R < 0)
					aang += 3.1415926;
			if (Q_R < 0)
				int kkkk = 1;
			WHH=SPLH(guidet,aang,i );
			WBB=SPLB(guidet,aang,i);
			
			H_R = WHH*(n_R*n_R + Q_R*Q_R)/(pow(guidet + TUR[i].Cy, 2) - WHH*TUR[i].Ch);
			M_R = (WBB /WHH- TUR[i].k2)* H_R;
			H_T = H_R*TUR[i].H0;
			Qi = (QCP(I, L[I].NN)* CQM(J, 0) + QCM(J, 0)*CQP(I, L[I].NN) - CQM(J, 0) * CQP(I, L[I].NN)* H_T) / (CQM(J, 0) + CQP(I, L[I].NN));
			/////////
			if(T>TUR[i].cctime)
				ni = n_T[i][2] + TUR[i].n0 / TUR[i].Ta*TimeStep*M_R;
			else 
				ni = n_T[i][2];
			if ((fabs(Qi - Q_T[i][3]) < 0.00001)&&(fabs(ni - n_T[i][3]) < 0.00001))
				break;
			else
			{
				n_T[i][3] = ni;
				Q_T[i][3] = Qi;
				h_T[i][3] = H_T;
			}
		}
		/////////////////Iteration over

		if (k0 > IMAX)
		{
			printf("Iteration failed");
		}
		for (j = 0; j < 3; j++)
		{
			n_T[i][j] = n_T[i][j + 1];
			Q_T[i][j] = Q_T[i][j + 1];
			h_T[i][j] =h_T[i][j + 1];
		}
		QP[I][L[I].NN] = Qi;
		HP[I][L[I].NN] = QCP(I, L[I].NN) / CQP(I, L[I].NN) - 1.0 / CQP(I, L[I].NN)*QP[I][L[I].NN];

		QP[J][0] = Qi;
		HP[J][0] = QP[J][0] / CQM(J, 0) - QCM(J, 0) / CQM(J, 0); 
		double P = 1.0262536*TUR[i].M0*M_R*n_T[i][3];
		pw = HP[I][L[I].NN] - TUR[i].Tz;
		ps= HP[J][0] - TUR[i].Tz;
		if (densi == 1||densi%den==0)
			fprintf(fp1[FileNumber], "%lf\t %lf\t %lf\t%lf\t%lf \t%lf\t%lf \t%lf \n", T, guidet, n_T[i][3], Q_T[i][3], H_T, pw, ps,WHH);
		densi++;
		return;
	}
	else if (JP->type == 5)            //Surge tank
	{
		if (T == TimeStep)
		{
			ST[i].QT = 0.0;
			ST[i].ZR = H[I][L[I].NN];
			fprintf(SgTk[ST[i].id], "%lf\t %lf\t%lf \n", 0.0, ST[i].QT, ST[i].ZR);
			for (k0 = 1;k0 <= 10;k0++)
			{
				if (((ST[i].ZR > ST[i].height[k0]) && (ST[i].ZR < ST[i].height[k0 + 1])) || (ST[i].ZR < ST[i].height[k0]) && (ST[i].ZR > ST[i].height[k0 + 1]))
					ST[i].area = ST[i].ar[k0];
			}
		}
		for (k0 = 1;k0 <= 10;k0++)
		{
			if (((ST[i].ZR > ST[i].height[k0]) && (ST[i].ZR < ST[i].height[k0 + 1])) || (ST[i].ZR < ST[i].height[k0]) && (ST[i].ZR > ST[i].height[k0 + 1]))
				ST[i].areap = ST[i].ar[k0];
		}
		ST[i].W0 = TimeStep / (ST[i].area+ ST[i].areap);
		SC = CQP(I, L[I].NN) + CQM(J, 0);
		SQ = QCP(I, L[I].NN) - QCM(J, 0);
		if (ST[i].QT > 0)
			k0 = 0;
		else k0 = 1;
		ST[i].QTP = (SQ - SC*(ST[i].W0*ST[i].QT + ST[i].ZR)) / (SC*(ST[i].W0 + ST[i].X0[k0]*fabs(ST[i].QT)) + 1.0);  //Flow of surge tank
		ST[i].ZRP = ST[i].ZR + (ST[i].QTP + ST[i].QT)*ST[i].W0;  //Head of surge tank
		fprintf(SgTk[ST[i].id], "%lf\t %lf\t%lf \n", T,ST[i].QTP, ST[i].ZRP);

		ST[i].ZR = ST[i].ZRP;
		ST[i].QT = ST[i].QTP;
		ST[i].area = ST[i].areap;

		HP[I][L[I].NN] = (SQ - ST[i].QTP) / SC;
		QP[I][L[I].NN] = QCP(I, L[I].NN) - CQP(I, L[I].NN)*HP[I][L[I].NN];

		HP[J][0] = HP[I][L[I].NN];
		QP[J][0] = QCM(J, 0) + CQM(J, 0)*HP[J][0];
	}
}

 //End Boundary
void EndBoundary(struct ENode *EP)
{ 
	int I;
	double TAU;
	I = EP->StartPipe ;
	if (EP->K == 1)         //valve closing
	{
		if (T <= TC[1][i])
		{
			TAU = pow((1.0 - T / TC[1][i]), EM[1][i]);
		}
		else
		{
			TAU = 0.0;
		}
		double C = TAU*CDA[1][i];
		double CP = QCP(I, L[I].NN) / CQP(I, L[I].NN);
		double B1 = 1.0 / CQP(I, L[I].NN);
		QP[I][L[I].NN] = -C*C*g*B1 + C*sqrt(C*C*g*g*B1*B1 + 2 * g*CP);
		HP[I][L[I].NN] = CP - B1*QP[I][L[I].NN];
	}
	else if (EP->K == 3)               //Non reflecting boundary
	{
		QP[I][L[I].NN] = (QCP(I, L[I].NN) + QCM(I, L[I].NN - 2)) / 2.0;
		HP[I][L[I].NN] = QCP(I, L[I].NN) / CQP(I, L[I].NN) - 1.0 / CQP(I, L[I].NN)*QP[I][L[I].NN];

	}
	else if (EP->K == 2)                //Downstream reservior
	{
		if (Q[I][L[I].NN] >= 0)
		{
			HP[I][L[I].NN] = HE[i];
			QP[I][L[I].NN] = QCP(I, L[I].NN) - CQP(I, L[I].NN)*HP[I][L[I].NN];
		}
	}
	else if (EP->K == 4)             //Pressure pulse
	{
		if (T <= 0.01)
		{
			HP[I][L[I].NN] = HR0 + 20 * T;
			QP[I][L[I].NN] = QCP(I, L[I].NN) - CQP(I, L[I].NN)*HP[I][L[I].NN];
		}
		else if (T <= 0.02)
		{
			HP[I][L[I].NN] = HR0 + 0.2 - (T - 0.01) * 20;
			QP[I][L[I].NN] = QCP(I, L[I].NN) - CQP(I, L[I].NN)*HP[I][L[I].NN];
		}
		else
		{
			QP[I][L[I].NN] = 0;
			HP[I][L[I].NN] = QCP(I, L[I].NN) / CQP(I, L[I].NN);
		}
	}
	else if (EP->K == 5)          //Flow pulse
	{
		if (T <= 0.005)
		{
			QP[I][L[I].NN] = 1 * T;
			HP[I][L[I].NN] = QCP(I, L[I].NN) / CQP(I, L[I].NN) - 1.0 / CQP(I, L[I].NN)*QP[I][L[I].NN];
		}
		else if (T <= 0.005)
		{
			QP[I][L[I].NN] = 0.05 - (T - 0.05) * 1;
			HP[I][L[I].NN] = QCP(I, L[I].NN) / CQP(I, L[I].NN) - 1.0 / CQP(I, L[I].NN)*QP[I][L[I].NN];
		}
		else
		{
			QP[I][L[I].NN] = 0;
			HP[I][L[I].NN] = QCP(I, L[I].NN) / CQP(I, L[I].NN);
		}
	}
	else if (EP->K == 6)               //flow = 0
	{
		QP[I][L[I].NN] = 0;
		HP[I][L[I].NN] = QCP(I, L[I].NN) / CQP(I, L[I].NN) - 1.0 / CQP(I, L[I].NN)*QP[I][L[I].NN];
	}
	return;
}


