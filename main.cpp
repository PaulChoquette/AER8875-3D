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
//#include "Metric.h"
#include "Solver.h"

//#include "./include/TECXXX.h"


using namespace std;

int main() {
	cout << "Executeble 2 Starting ..." << endl;

	// =================================== EXECUTABLE 2 ====================================================
	Reader_c FileContents;
	Solver_c solve;
	solve.cfl = 1.0;
	FileContents.read_file("Zone0.su2");
	//solve.ComputeLocalConnectivity();

	// METRIC
	cout << "\nMetriques ..." << endl;
	//Metric_c metric;
	//metric.Compute(solve, FileContents);

	//metric.SumNorm(solve, FileContents, 1);


	

	
	cout << "**************\nEnd\n**************\n";
	double* p = new double[FileContents.nelem];
	double* Rho = new double[FileContents.nelem];
	double* u = new double[FileContents.nelem];
	double* v = new double[FileContents.nelem];
	double* w = new double[FileContents.nelem];

	for (int i = 0; i < FileContents.nelem; i++) {
		p[i] = 0.;
		Rho[i] = 1.;
		u[i] = 0.;
		v[i] = 0.;
		w[i] = 0.;
	}


//FileContents.write_tecplot(FileContents, "test2", p, Rho, u, v, w);

	return 0;
}
