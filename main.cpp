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
	//string File_Name = "block.su2";
	string File_Name = "test_justSquare.su2";
	//string File_Name = "2cube.su2";
	//string File_Name = "naca0012_euler_65x65x2_O_1B.su2";
	Reader_c FileContents;

	FileContents.read_file(File_Name);

	// CONNECTIVITY
	Solver_c solve;
	
	solve.ComputeGlobalConnectivity(FileContents);
	solve.Display2DArray(solve.elem2elem_g, solve.ncell_g, 6, "elem2elem_g");
	solve.ComputeZoneConnectivity(FileContents);
	solve.ComputeLocalConnectivity();

	 // DISPLAY OF GLOBAL
	cout << "\n========================================= DISPLAY OF GLOBAL ========================================= " << endl;
	//solve.Display2DArray(FileContents.elem2node, FileContents.ncell, 8, "elem2node_g");
	//solve.Display1DArray(solve.esup1, solve.mesup1, "esup1");
	//solve.Display1DArray(solve.esup2, solve.nnode_g+1, "esup2");
	//solve.Display1DArray(solve.psup1, solve.mpsup, "psup1");
	//solve.Display1DArray(solve.psup2, solve.nnode_g+1, "psup2");
	//solve.Display2DArray(solve.elem2elem_g, solve.ncell_g, 6, "elem2elem_g");

	//// DISPLAY OF ZONE
	cout << "\n========================================= DISPLAY OF ZONE ========================================= " << endl;
	int izone = 0; // iZone to Display
	//solve.Display1DArray(solve.zone2nnode, solve.nzone, "zone2nnode");
	//solve.Display1DArray(solve.zone2nelem, solve.nzone, "zone2nelem");
	//solve.Display2DArray(solve.zone2node, solve.nzone, solve.zone2nnode[izone], "zone2node");
	//solve.Display2DArray(solve.zone2elem, solve.nzone, solve.zone2nelem[izone], "zone2elem");
	//solve.Display2DArray(solve.nodeglobal2local, solve.nnode_g, solve.nzone, "nodeglobal2local");
	//solve.Display2DArray(solve.elemglobal2local, solve.nelem_g, 2, "elemglobal2local");
	//solve.Display3DArray(solve.belem2node, izone, solve.zone2nbelem[izone],4, "belem2node");
	//solve.Display3DArray(solve.elem2node, izone, solve.zone2ncell[izone], 4, "elem2node");
	//solve.Display3DArray(solve.elem2elem, 0, solve.zone2ncell[0], 4, "elem2elem");
	//solve.Display2DArray(solve.elem2vtk, solve.nzone, solve.zone2ncell[izone], "elem2vtk");
	//solve.Display3DArray(solve.zone2coord, 0, solve.zone2nnode[0], 2, "zone2coord");
	solve.Display3DArray(solve.face2node, 0, solve.zone2nface[0], 5, "face2node");
	solve.Display3DArray(solve.face2elem, 0, solve.zone2nface[0], 2, "face2elem");
	solve.Display3DArray(solve.face2fael, 0, solve.zone2nface[0], 2, "face2fael");
	solve.Display3DArray(solve.elem2face, 0, solve.zone2ncell[0], 6, "elem2face");
	
	solve.MetricExemple(0);
	solve.SolverExemple();



	return 0;
}
