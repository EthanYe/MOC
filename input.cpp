#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"parameters.h"
#include"suter_characteristic_curve.h"

FILE *fp,*fpst,*fpg;
char c[256] = {}, c1[256] = {}, buf[256],name[256];
char pro_name[20];            //project name(blh or tdt)
int m,m1,m2,k0;
int PipeNumber;
int NodeNumber;

void ControlParameter();
void Reservoir();
void Pipe();
void Valve();
void Turbine();
void SurgeTank();
void Guide();
void StartBoundary();
void EndBoundary();
void Output();

void input()
{	
////////Input option file ////////
	fopen_s(&fp, "option.txt", "r");
	//Make sure the file exists
	if (fp == NULL)  
	{
		printf("Openning file input.txt failed\n");
		exit(0);
	}

	fscanf_s(fp, "%d", &m);//Read option number
	//find [option]
	do
		fgets(buf, 256, fp);
	while (strcmp(buf, "[option]\n") != 0);
	//store pro_name
	do
	{
		fscanf_s(fp, "%s", buf, 256);
	} 
	while (buf[0] != '0' + m);
	fscanf_s(fp, "%s", pro_name, 256);
	fclose(fp);
////////Input option file over////////

////////Input file ////////
	strcpy_s(name, pro_name);
	strcat_s(name, "_input.txt");
	fopen_s(&fp,name, "r"); 
	//Make sure the file exists
	if (fp == NULL) 
	{
		printf("Openning file input.txt failed\n");
		exit(0);
	}
	printf("Open input file successfully\n");

	ControlParameter();
	Reservoir();
	Pipe();
	Valve();
	Turbine();
	SurgeTank();
	Guide();
	StartBoundary();
	EndBoundary();
	Output();
	fclose(fp);

	double AverageArea[30];                  //Average area of pipe
	for (i = 1; i <= N; i++)
	{
		AverageArea[i] = 3.14159*L[i].D * L[i].D / 4;
	}
	for (i = 1; i <= N; i++)
		L[i].AverageR = 0.25*L[i].D;

	for (i = 1; i <= N; i++)                   //Loss coefficient
	{
		for (j = 0; j<2; j++)
			L[i].PipeLoss[j] = L[i].m_Length * pow(L[i].PRoughness[j], 2) / (pow(AverageArea[i], 2)*pow(L[i].AverageR, 4 / 3.0));
	}

	for (i = 1; i <= N; i++)                        //Initialization of pipe section
	{
		for (int k0 = 0; k0 <= L[i].NN; k0++)
		{
			L[i].AJD[k0] = 3.14159*L[i].D * L[i].D / 4;//Area
			L[i].SJD[k0] = 3.14159*L[i].D;           //Circumference
			L[i].DX[k0] = L[i].m_Length / L[i].NN;     //Space step
		}
	}
}
void ControlParameter()//input [CONTROL PARAMETERS]
{
	rewind(fp);
	//Find [CONTROL PARAMETERS]
	do
		fgets(buf, 256, fp);
	while (strcmp(buf, "[CONTROL PARAMETERS]\n") != 0); 
	//input paras
	do
	{
		fscanf_s(fp, "%s", buf, 256);
		if (strcmp(buf, "Time") == 0)
		{
			fscanf_s(fp, "%*s%lf", &TimeStep);//Input Time step
		}
		else if (strcmp(buf, "TotalTimes") == 0)             
			fscanf_s(fp, "%lf", &TotalTimes); //Input TotalTimes
		else if (strcmp(buf, "Gravity") == 0)                 
			fscanf_s(fp, "%lf", &g);//Input Gravity
		else if (strcmp(buf, "Print") == 0)
			fscanf_s(fp, "%*s%d", &den);//Input Gravity
	} while (buf[0] != '[');
}
void Reservoir()	//input [RESERVOIRS]
{
	rewind(fp);
	//Find [RESERVOIRS]
	do
		fgets(buf, 256, fp);
	while (strcmp(buf, "[RESERVOIRS]\n") != 0);   
	//input paras
	do
	{
		fscanf_s(fp, "%s", buf, 256);
		if (strcmp(buf, "R1") == 0)
			fscanf_s(fp, "%lf", &HR0);                         //Pressure of upstream reservior
		if (buf[0] == 'E')
		{
			m = buf[1] - '0';
			fscanf_s(fp, "%lf", &HE[m]);                        //Pressure of downstream reservior
		}
	} while (buf[0] != '[');
}
void Pipe()	//input [PIPES]
{
	rewind(fp);
	//Find [PIPES]
	do
		fgets(buf, 256, fp);
	while (strcmp(buf, "[PIPES]\n") != 0);                   
	N = 0;//Initialize the number of pipe
	//input paras
	do
	{
		fscanf_s(fp, "%s", buf, 256);
		if (buf[0] == 'L')
		{
			N++;										//Calculate the number of pipe		   
			if (buf[2] > 47 && buf[2] < 58)
				PipeNumber = (buf[1] - '0') * 10 + (buf[2] - '0');
			else
				PipeNumber = buf[1] - '0';
				/////////////////////////Previous node
			fscanf_s(fp, "%s", buf, 256);
			if (buf[0] == 'J')
			{
				if (buf[2] > 47 && buf[2] < 58)
					NodeNumber = (buf[1] - '0') * 10 + buf[2] - '0';
				else
					NodeNumber = buf[1] - '0';               //NodeNumber
				JN[NodeNumber].EndPipeNumber++;           //EndPipeNumber
				int n = JN[NodeNumber].EndPipeNumber;
				JN[NodeNumber].EndPipe[n] = PipeNumber;   //EndPipe
				L[PipeNumber].StartNode = NodeNumber + 1;
			}
			else if (buf[0] == 'R')                           //Upsteam node                     
			{
				JN[0].EndPipeNumber++;                     //EndPipeNumber
				JN[0].EndPipe[1] = PipeNumber;             //Upsteam EndPipe
				L[PipeNumber].StartNode = 1;
			}
			////////////////////////End node
			fscanf_s(fp, "%s", buf, 256);
			if (buf[0] == 'J')
			{
				if (buf[2] > 47 && buf[2] < 58)
					NodeNumber = (buf[1] - '0') * 10 + buf[2] - '0';
				else
					NodeNumber = buf[1] - '0';                    //NodeNumber
				JN[NodeNumber].StartPipeNumber++;         //StartPipeNumber
				int n = JN[NodeNumber].StartPipeNumber;
				JN[NodeNumber].StartPipe[n] = PipeNumber; //StartPipe
				L[PipeNumber].EndNode = NodeNumber + 1;
			}
			else if (buf[0] == 'E')                           //Downstream                        
			{
				NodeNumber = buf[1] - '0';
				EN[NodeNumber].StartPipe = PipeNumber;   //Downstream StartPipe
			}
			fscanf_s(fp, "%lf%lf%lf%lf%lf%lf%lf", &L[PipeNumber].a, &L[PipeNumber].D, &L[PipeNumber].m_Length, &L[PipeNumber].PRoughness[0], &L[PipeNumber].PRoughness[1], &L[PipeNumber].LLCoefficient[0], &L[PipeNumber].LLCoefficient[1]);//Parameters of pipe
		}
	} while (buf[0] != '[');
	TotalNode = 0;
	for (i = 0; i <= N; i++)
		JN[i].type = 0;                                       //Default node type
	for (i = 1; i < N + 1; i++)
	{
		if (JN[i].StartPipeNumber != 0)
		{
			TotalNode++;                                      //TotalNode
			if (JN[i].StartPipeNumber == 1 && JN[i].EndPipeNumber == 1)
			{
				JN[i].type = 1;                          //one-to-one
			}
			else JN[i].type = 3;                             //Bifurcated pipe
		}
	}
	////////////////
	rewind(fp);
	do
		fgets(buf, 256, fp);
	while (strcmp(buf, "[PIPES]\n") != 0);                    //Find [PIPES]
	do
	{
		fscanf_s(fp, "%s", buf, 256);
		if (buf[0] == 'L')
		{
			if (buf[2] > 47 && buf[2] < 58)
				PipeNumber = (buf[1] - '0') * 10 + buf[2] - '0';
			else
				PipeNumber = buf[1] - '0';
			/////////////////////////Previous node
			fscanf_s(fp, "%s", buf, 256);
			fscanf_s(fp, "%s", buf, 256);
			if (buf[0] == 'E')                           //Downstream                        
			{
				NodeNumber = buf[1] - '0';
				L[PipeNumber].EndNode = NodeNumber + TotalNode + 1;   //Downstream StartPipe
			}
		}
	} while (buf[0] != '[');


}
void Valve()
{
	//input [VALVE]
	rewind(fp);
	do
		fgets(buf, 256, fp);
	while (strcmp(buf, "[VALVE]\n") != 0);                 //Find [VALVE]
	do
	{
		fscanf_s(fp, "%s", buf, 256);
		if (buf[0] == 'V')
		{
			fscanf_s(fp, "%s", buf, 256);
			if (buf[0] == 'J')                               //The node where the in-line valve is 
			{
				m = buf[1] - '0';
				fscanf_s(fp, "%lf%lf%lf", &CDA[0][m], &TC[0][m], &EM[0][m]);
				JN[m].type = 2;                              //Change the type of the node
			}
			if (buf[0] == 'E')                               //End-line valve
			{
				m = buf[1] - '0';
				fscanf_s(fp, "%lf%lf%lf", &CDA[1][m], &TC[1][m], &EM[1][m]);
			}
		}
	} while (buf[0] != '[');
	for (i = 1; i <= N; i++)                  	//Pipe segmentation
	{
		L[i].NN = (int)L[i].m_Length / TimeStep / L[i].a + 1;
		L[i].a = L[i].m_Length / TimeStep / L[i].NN;
	}
}
void Turbine()
{
	//input [TURBINWE]
	rewind(fp);
	do
		fgets(buf, 256, fp);
	while (strcmp(buf, "[TURBINWE]\n") != 0);                 //Find [TURBINWE]
	do
	{
		fscanf_s(fp, "%s", buf, 256);
		if (buf[0] == 'T'&& buf[1] == 'u')
		{
			turbine += 1;
			fscanf_s(fp, "%s", buf, 256);
			if (buf[0] == 'J')                               //The node where the in-line valve is set
			{
				m = buf[1] - '0';
				JN[m].type = 4;                              //Change the type of the node
			}
			fscanf_s(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf", &TUR[m].H0, &TUR[m].Q0, &TUR[m].n0, &TUR[m].M0, &TUR[m].D1, &TUR[m].D2, &TUR[m].Ta, &TUR[m].Tz,&TUR[m].cctime);
			suter(m,pro_name);									//suter transformation and spline
		}
	} while (buf[0] != '[');
}
void SurgeTank()
{
	//input [SURGE TANK]
	rewind(fp);
	do
		fgets(buf, 256, fp);
	while (strcmp(buf, "[SURGE TANK]\n") != 0);                 //Find [SURGE TANK]
	do
	{
		fscanf_s(fp, "%s", buf, 256);
		if (buf[0] == 'S'&&buf[1] == 'T')
		{
			k0 = buf[2] - '0';			//id
			////Write surge tank file////
			strcpy_s(c1, buf);
			strcat_s(c1, 256, ".xls");
			fopen_s(&SgTk[k0], c1, "w");
			fprintf(SgTk[k0], "Time  \tFlow  \tHead  \n");
			////Writing surge tank file over////

			
			strcpy_s(c, buf);
			strcat_s(c, 256, ".txt");           //File name of reading
			fscanf_s(fp, "%s", buf, 256);
			if (buf[0] == 'J')                               //The node where the surge tank is set 
			{
				m1 = buf[1] - '0';
				JN[m1].type = 5;                              //Change the type of the node
				ST[m1].id = k0;
			}
			fscanf_s(fp, "%lf%lf", &ST[m1].X0[0], &ST[m1].X0[1]);

			strcpy_s(name, pro_name);
			strcpy_s(buf,"_");
			strcat_s(buf, c);
			strcat_s(name, buf);
			fopen_s(&fpst, name, "r");           //Read
			do
			{
				fscanf_s(fpst, "%s", c, 256);
				if (c[0] == 'N')
				{
					m2 = c[4] - '0';
					fscanf_s(fpst, "%lf%lf", &ST[m1].height[m2], &ST[m1].ar[m2]);
				}
			} while (!feof(fpst));
			fclose(fpst);
		}
	} while (buf[0] != '[');
}
void Guide()
{
	strcpy_s(name, pro_name);
	strcat_s(name, "_guide.txt");
	fopen_s(&fpg, name,"r");
	//input [servomotor-guide]
	do
		fgets(buf, 256, fpg);
	while (strcmp(buf, "[servomotor-guide]\n") != 0);            //Find [servomotor-guide]
	do
	{
		fscanf_s(fpg,"%s",buf, 256);
		if(buf[0]=='N')
			if (buf[2] > 47 && buf[2] < 58)
				k0 = (buf[1] - '0') * 10 + buf[2] - '0';
			else
				k0 = buf[1] - '0';
		fscanf_s(fpg, "%lf%lf", &TUR[m].ServoGui[k0][0], &TUR[m].ServoGui[k0][1]);
	}
	while (buf[0]!='[');
	//input [time-servomotor]
	rewind(fpg);
	do
		fgets(buf, 256, fpg);
	while (strcmp(buf, "[time-servomotor]\n") != 0);            //Find [time-servomotor]
	do
	{
		fscanf_s(fpg, "%s", buf, 256);
		if (buf[0] == 'N')
			if (buf[2] > 47 && buf[2] < 58)
				k0 = (buf[1] - '0') * 10 + buf[2] - '0';
			else
				k0 = buf[1] - '0';
		fscanf_s(fpg, "%lf%lf", &TUR[m].TimeServo[k0][0], &TUR[m].TimeServo[k0][1]);
	} while (buf[0] != '[');
	fclose(fpg);
}

void StartBoundary()
{
	;
}
void EndBoundary()
{
	//input END BOUNDARY
	rewind(fp);
	do
		fgets(buf, 256, fp);
	while (strcmp(buf, "[END BOUNDARY]\n") != 0);            //Find [END BOUNDARY]
	TotalENode = 0;
	do
	{
		fscanf_s(fp, "%s", buf, 256);
		if (buf[0] == 'E')
		{
			m = buf[1] - '0';
			fscanf_s(fp, "%s", buf, 256);
			m1 = buf[0] - '0';
			if (m1 == 4)
				m1 = m1 + buf[2] - '1';
			EN[m].K = m1;
			TotalENode++;
		}
	} while (buf[0] != '[');
}
void Output()
{
	rewind(fp);
	FileNumber = 0;                                             //Initial the number of output_file
	do
		fgets(buf, 256, fp);
	while (strcmp(buf, "[OUTPUT FILE]\n") != 0);          //Find [OUTPUT FILE]
	do
	{
		fscanf_s(fp, "%s", buf, 256);
		if (buf[0] == 'L')
		{
			strcat_s(buf, 256, ".xls");                     //Name of output file
			fopen_s(&fp1[FileNumber], buf, "w");           //Open output file
			m = buf[1] - '0';
			NodePosition[0][FileNumber] = m;               //The node the file write
			if (buf[3] == 'S')
				NodePosition[1][FileNumber] = 0;
			else if (buf[3] == 'E')
				NodePosition[1][FileNumber] = L[m].NN;
			FileNumber++;
		}
		else if (buf[0] == 'E')
		{
			strcat_s(buf, 256, ".xls");
			fopen_s(&fp1[FileNumber], buf, "w");
			m = buf[1] - '0';
			NodePosition[0][FileNumber] = EN[m].StartPipe;
			NodePosition[1][FileNumber] = L[NodePosition[0][FileNumber]].NN;
			FileNumber++;
		}

	} while (buf[0] != '[');
	fopen_s(&fp1[FileNumber], "Unsteady Flow Message.xls", "w");
	fprintf(fp1[FileNumber], "Time \t Guide\t Speed \tFlow \tHead  \tSpiral case pressure \tDraft tube \n");
	fopen_s(&fp1[FileNumber + 1], "Pipe result of steady flow .xls", "w");
	fprintf(fp1[FileNumber + 1], "Pipe number  \tFlow  \tStarting pressure \tEnding pressure  \n");
}