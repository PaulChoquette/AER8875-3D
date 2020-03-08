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
	double AoA = 0;
	double cfl = 0.1;	//WTF
	bool RK_M=1;
	int RK_step = 1;
	int Order = 1;
	int iterMax = 300;
	double convergeCrit = pow(10,-13);

	Reader_c FileContents;
	solver_c solve(mach,AoA,Order,RK_step,RK_M,cfl,iterMax,convergeCrit);
	cout << "Zone "<<solve.World.world_rank<<" Starting ..." << endl;
	fileName = "Zone"+to_string(solve.World.world_rank)+".su2";
	cout<<fileName<<endl;
	FileContents.read_file_local(fileName);
	solve.ComputeLocalConnectivity(FileContents);
	solve.ComputeMetric(FileContents);
	solve.SumNorm(FileContents, 1);	
	//cout << "**************\nEND OF PREMILINARY OPERATIONS FOR "<<solve.World.world_rank<<"\n**************\n";

	solve.InitMPIBuffer(FileContents);
	solve.Compute();		//Necessary for no seg faults
	//solve.HighlightZoneBorder();
	
	// cout << "**************\nGG NO RE\n**************\n";


	// //PURE EVIIIIIIIIIIIIIIIIIIIIIIIIIIIIL
	switch (solve.World.world_rank){
		case 0 : FileContents.write_tecplot(FileContents,"./PlotOut/IhAtEmYsElF0", solve.p, solve.rho, solve.u, solve.v, solve.w);break;
		case 1 : FileContents.write_tecplot(FileContents,"./PlotOut/IhAtEmYsElF1", solve.p, solve.rho, solve.u, solve.v, solve.w);break;
		case 2 : FileContents.write_tecplot(FileContents,"./PlotOut/IhAtEmYsElF2", solve.p, solve.rho, solve.u, solve.v, solve.w);break;
		case 3 : FileContents.write_tecplot(FileContents,"./PlotOut/IhAtEmYsElF3", solve.p, solve.rho, solve.u, solve.v, solve.w);break;
	}

	return 1;
}
