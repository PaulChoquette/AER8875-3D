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
	Reader_c FileContents;
	string parametreFile = "CFDsimPI4.txt";
	double convergeFixLimit = pow(10,-6);		//Gèle les limiteurs après convergence jusqu'à cette résolution là. Évite é vventuelles oscillations
	//////////////////////////////////////////////////////////////////////////////
	solver_c solve(FileContents, convergeFixLimit);
	//solve.PrintStylz();
	FileContents.computePrmt(parametreFile);
	fileName = FileContents.su2pFilePath+to_string(solve.World.world_rank)+".su2";
	cout<<fileName<<endl;
	//fileName = "naca0012_euler_65x65x2_ONE0.su2";
	FileContents.read_file_local(fileName);
	solve.ComputeLocalConnectivity(FileContents);
	solve.ComputeMetric(FileContents);
	solve.SumNorm(FileContents, 1);

	solve.InitMPIBuffer(FileContents);
	solve.Compute();		//Necessary for no seg faults
	//solve.SetAnalyticalGradiant(0.69,0.0,0.0);
	//solve.PrintGradiant();
	//solve.LimitTecplot();
	//solve.HighlightZoneBorder();
	solve.ComputeCoefficient();

	FileContents.write_tecplot(FileContents,"./Tecio",solve.World.world_rank, solve.p, solve.rho, solve.u, solve.v, solve.w);
	// To change
	if (solve.World.world_rank==0){
		solve.PrintStylz();
	}

	string OFileName = "./AsciiFile"+to_string(solve.World.world_rank)+".dat";
	FileContents.write_tecplot_ASCII(OFileName,solve.p,solve.rho,solve.u,solve.v,solve.w);
	return 1;
}
