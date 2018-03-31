#pragma once

extern int N;                        //Number of pipes

extern int i, j;                     //For loop
extern int FileNumber;               //Number of file
extern int NodePosition[2][30];      //The node position in output file


extern double HR0;                   // Initial upstream head
extern double HE[5];                 //Downstream head
extern double g;                     //Gravity
extern double TimeStep;		         //Time step
extern double TotalTimes;            //Total Times
extern int    den;                   //Print density
extern double CDA[2][30];            //Initial opening of valve
extern double TC[2][30];             //closing valve time
extern double EM[2][30];
extern double T;


extern double H[30][4000];            //The head of the previous moment                           
extern double Q[30][4000];            //The flow at the previous moment
extern double HP[30][4000];
extern double QP[30][4000];

//////////////
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
		double ANG[300];                      //Node number of characteristic curve
	}CC[300];
	double k1, k2 , Ch , Cy ;
	double n0, Q0, M0, H0;           //Rated parameter
	double Ta;                       //Ta
	double D1;                      //Turbine diameter
	double D2;                       //draft tube diameter
	double Tz;                      //Erection elevation of unit
	double cctime;                  //ConditionConversionTime
	double ServoGui[20][2];
	double TimeServo[20][2];
};
extern struct Turbine TUR[30];         //The number in [] represent which node the turbine is at 
extern int turbine;
									 /////////////////

struct Node
{
	int StartPipe[10];                //StartPipe 
	int StartPipeNumber;              //StartPipeNumber
	int EndPipe[10];                  //EndPipe
	int EndPipeNumber;                //EndPipeNumber
	int type;                         //Node type
}; 
extern struct Node JN[30];
extern int TotalNode;

struct ENode
{
	int StartPipe;                    //StartPipe 
	int K;                            //EndNode type
};
extern struct ENode EN[5];
extern int TotalENode;


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
	double LLCoefficient[2];          //Local loss coefficient
	int StartNode;               //StartNode 
	int EndNode;                 //EndNNode
};
extern struct Pipe L[30];

struct SurgeTank
{
	double X0[2];               //Resistance coefficient (in and out)
	double area;                //EndNode type
	double areap;
	double height[10], ar[10];	//storage of relationship between height and area
	double QTP, QT;				//flow of surge tank
	double ZRP, ZR;				//Head of surge tank
	double W0;
	int id;
};
extern struct SurgeTank ST[30];
extern FILE *fp1[50];
extern FILE *SgTk[5];
	