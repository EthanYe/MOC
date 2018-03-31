#include<stdio.h>
#include<math.h>
#include"input.h"
#include"steady_flow_calculation.h"
#include"moc.h"


int N;                        //Number of pipes
int i,j;                      //For loop

int FileNumber;               //Number of output file
int NodePosition[2][30];      //The node in output file

double g;                     //Gravity
double TimeStep;		      //Time step
double T;                      //Time 
double TotalTimes;            //Total Times
int    den;                   //Print density

double CDA[2][30];             //Initial opening of valve,1 represent JN,2 represent EN
double TC[2][30];              //closing valve time
double EM[2][30];

double HR0;                   //Upstream head
double HE[5];                 //Downstream head
double H[30][4000];            //The head at the previous moment                           
double Q[30][4000];            //The flow at the previous moment
double HP[30][4000];
double QP[30][4000];

struct Turbine
{
	int NCC;                         //Node number of characteristic curve
	int LCC;                         //Lines of characteristic curve
	double guide[300];               //guide vane degree
	struct chacurve
	{
		double h[300];                   //x[i]-x[i-1]
		double MH[300], MB[300];
		double WH[300], WB[300];
		double ANG[300];
	}CC[300];
	double k1 = 1.2, k2 = 4.2, Ch = 0.4, Cy = 0.2;
	double n0, Q0, M0, H0;           //Rated parameter
	double Ta;                       //Ta
	double D1;                      //Turbine diameter
	double D2;                       //draft tube diameter
	double Tz;                      //Erection elevation of unit
	double cctime;                  //ConditionConversionTime
	double ServoGui[20][2];
	double TimeServo[20][2];
}TUR[30];                         //The number in [] represent which node the turbine is at ;
int turbine=0;            

struct Node
{
	int StartPipe[10];          //StartPipe 
	int StartPipeNumber;        //StartPipeNumber
	int EndPipe[10];            //EndPipe
	int EndPipeNumber;          //EndPipeNumber
	int type;                   //Node type
} JN[30];
int TotalNode;

struct ENode
{
	int StartPipe;               //StartPipe 
	int K;                       //EndNode type
} EN[5];
int TotalENode;

struct Pipe
{
	double m_Length;             //Length of pipe
	double a;                 //wave speed
	double D;                 //Diameter of pipe
	double AJD[4000];          //Sectional area
	double SJD[4000];          //Sectional circumference
	double DX[4000];           //Space step
	double AverageR;          //Hydraulic radius
	int NN;                   //Segment number
	double PipeLoss[2];       //Loss coefficient in characteristic method
	double PRoughness[2];		  //Roughness
	double LLCoefficient[2];       //Local loss coefficient
	int StartNode;               //StartNode 
	int EndNode;                 //EndNNode
} L[30];

struct SurgeTank
{
	double X0[2];               //Resistance coefficient (in and out)
	double area;                //
	double areap;               //
	double height[10],ar[10];	//storage of relationship between height and area
	double QTP, QT;				//Flow of surge tank
	double ZRP, ZR;				//Head of surge tank
	double W0;
	int id;
} ST[30];

FILE *fp1[50];
FILE *SgTk[5];

//调压室问题
//无反射边界条件的恒定流计算问题
//恒定流流量初始化问题
//转动惯量和Ta值的问题
//甩负荷小开度与topsys比较问题

void main()
{
	input();                                //input the parameters
	T = TimeStep;                           //Current time
	SteadyFlowCalculation();

top1000:
	moc();                  //moc
	//Save result
	for (i = 0;i<FileNumber;i++)
		fprintf(fp1[i], "%lf\t %lf \t%lf \n", T, QP[NodePosition[0][i]][NodePosition[1][i]], HP[NodePosition[0][i]][NodePosition[1][i]]);

	T = T + TimeStep;
	if (T>9)
		T = T;
	if (T<TotalTimes) goto  top1000;
	_fcloseall();
}



