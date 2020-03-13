#include <string>
#include <cctype>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "omp.h"
#include <chrono>

#include "main.h"
#include "UI.h"
#include "Reader.h"
#include "Connect.h"
#include "Metric.h"
#include "Solver.h"
//#include <metis.h>
//#include "./include/TECXXX.h"


using namespace std;

string fileName;

int main() {
	// =================================== EXECUTABLE 2 ====================================================
	double mach = 2.0;
	double AoA = 0.0;
	double cfl = 2.3;
	bool RK_M=1;
	int RK_step =5;
	int Order = 1;
	int iterMax = 2000;
	double convergeCrit = pow(10,-13);
	double convergeFixLimit = pow(10,-4);

	Reader_c FileContents;
	solver_c solve(mach,AoA,Order,RK_step,RK_M,cfl,iterMax,convergeCrit,convergeFixLimit);
	//solve.PrintStylz();
	fileName = "Zone"+to_string(solve.World.world_rank)+" (1).su2";
	FileContents.read_file_local(fileName);
	solve.ComputeLocalConnectivity(FileContents);
	solve.ComputeMetric(FileContents);
	solve.SumNorm(FileContents, 1);	

	solve.InitMPIBuffer(FileContents);
	solve.Compute();		//Necessary for no seg faults
	//solve.SetAnalyticalGradiant(0.69,0.0,0.0);
	//solve.PrintGradiant();
	//solve.HighlightZoneBorder();


	// To change
	switch (solve.World.world_rank){
		case 0 : FileContents.write_tecplot(FileContents,"./PlotOut/Tecio0", solve.elem2vol, solve.rho, solve.u, solve.v, solve.w);solve.PrintStylz();break;
		case 1 : FileContents.write_tecplot(FileContents,"./PlotOut/Tecio1", solve.elem2vol, solve.rho, solve.u, solve.v, solve.w);break;
		case 2 : FileContents.write_tecplot(FileContents,"./PlotOut/Tecio2", solve.elem2vol, solve.rho, solve.u, solve.v, solve.w);break;
		case 3 : FileContents.write_tecplot(FileContents,"./PlotOut/Tecio3", solve.elem2vol, solve.rho, solve.u, solve.v, solve.w);break;
	}

	string OFileName = "./PlotOut/AsciiFile"+to_string(solve.World.world_rank)+".dat";
	FileContents.write_tecplot_ASCII(OFileName,solve.p,solve.rho,solve.u,solve.v,solve.w);
	return 1;
}