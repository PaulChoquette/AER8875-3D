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
#include "METIS.h"
#include "Connect.h"
#include "Metric.h"
#include "Solver.h"


using namespace std;

int main() {
	cout << "Starting ..." << endl;
	string File_Name = "block.su2";
	//string File_Name = "naca0012_euler_65x65x2_O_1B.su2";
	Reader_c FileContents;

	FileContents.read_file(File_Name);

	// CONNECTIVITY
	Solver_c solve;
	solve.ComputeGlobalConnectivity(FileContents);
	solve.ComputeZoneConnectivity(FileContents);
	solve.ComputeLocalConnectivity();


    //solve.Display2DArray(solve.elem2vtk, solve.nzone, solve.zone2ncell[0], "elem2vtk");
    solve.Display3DArray(solve.elem2node, 0, solve.zone2ncell[0], 4, "elem2node");
	solve.Display3DArray(solve.zone2coord, 0, solve.zone2nnode[0], 2, "zone2coord");
	solve.Display1DArray(solve.elem2vtk[0], solve.zone2ncell[0], "elem2vtk");
    //solve.Display2DArray(FileContents.elem2node, solve.ncell_g, 5, "elem2node_g");
	//solve.Display1DArray(solve.esup1, solve.mesup1, "esup1");
	//solve.Display2DArray(solve.elem2elem_g, solve.ncell_g, 5, "elem2elem_g");
	//solve.Display2DArray(FileContents.elem2node, solve.nelem_g, 5, "elem2node_g");
	//solve.Display2DArray(FileContents.coord, solve.nnode_g, solve.ndime, "coord");
	//solve.Display1DArray(solve.psup1, solve.mpsup, "psup1");
	//solve.Display1DArray(solve.psup2, solve.nnode_g+1, "psup2");
	solve.Display2DArray(solve.elem2elem_g, solve.ncell_g, 5, "elem2elem");
	solve.Display1DArray(solve.zone2nnode, 4, "zone2nnode");
	solve.Display1DArray(solve.zone2nelem, 4, "zone2nelem");
	solve.Display2DArray(solve.zone2node, solve.nzone, solve.zone2nnode[0], "zone2node");
	solve.Display2DArray(solve.zone2elem, solve.nzone, solve.zone2nelem[0], "zone2elem");
	solve.Display2DArray(solve.nodeglobal2local, solve.nnode_g, solve.nzone, "nodeglobal2local");
	solve.Display2DArray(solve.elemglobal2local, solve.nelem_g, 2, "elemglobal2local");
	solve.Display3DArray(solve.belem2node, 0, solve.zone2nbelem[0],4, "belem2node");
	solve.Display3DArray(solve.elem2node, 0, solve.zone2ncell[0], 4, "elem2node");
	solve.Display3DArray(solve.elem2elem, 0, solve.zone2ncell[0], 4, "elem2elem");
	solve.MetricExemple(0);
	solve.SolverExemple();



	return 0;
}
