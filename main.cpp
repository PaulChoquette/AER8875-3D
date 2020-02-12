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




	Solver_c solve;
	solve.ConnectExemple(FileContents);
	solve.InitializeGlobal(FileContents);
	solve.Node2Elements(FileContents);
	solve.Node2Nodes(FileContents);
	solve.Element2Elements(FileContents);

	solve.InitializeLocal();
	solve.Zone2nnode();
	solve.Zone2nelem();
	solve.Zone2Nodes();
	solve.Zone2Elements();
	solve.NodeGlobal2Local();
	solve.ElementGlobal2Local();
	solve.Zone2Coord(FileContents);
	solve.Element2Nodes(FileContents);




	//solve.Display2DArray(FileContents.elem2node, solve.nelem_g, 5, "elem2node");
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
	solve.MetricExemple(0);
	solve.SolverExemple();



	return 0;
}
