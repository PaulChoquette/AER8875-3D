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


#include "./include/TECXXX.h"
using namespace std;

int main() {
	cout << "Starting ..." << endl;
	//string File_Name = "block.su2";
	//string File_Name = "test_justSquare.su2";
	//string File_Name = "test.su2";
	//string File_Name = "2cube.su2";
	string File_Name = "naca0012_euler_65x65x2_O_1B.su2";
	Reader_c FileContents;
	FileContents.read_file(File_Name);
	FileContents.check();
///////////////////////////////////////////////////////////////////////////////////////
	// CONNECTIVITY
	Solver_c solve;
	solve.cfl = 1.0;
	solve.ComputeGlobalConnectivity(FileContents);
//	solve.Display2DArray(solve.elem2elem_g, solve.ncell_g, 6, "elem2elem_g");
	solve.ComputeZoneConnectivity(FileContents);
	solve.ComputeLocalConnectivity();
//	solve.Display3DArray(solve.face2node, 0, solve.zone2nface[0], 5, "face2node");
//	solve.Display2DArray(solve.face2Nbr_of_node, 1, solve.zone2nface[0], "face2Nbr_of_node");
	// METRIC
	cout << "\nMetriques ..." << endl;
	Metric_c metric;
	metric.Compute(solve, FileContents);
//	solve.Display3DArray(metric.Face2Norm, 0, solve.zone2nface[0], 3, "Face2Norm");
//	solve.Display2DArray(metric.Face2Area, 1, solve.zone2nface[0], "Face2Area");
//	solve.Display2DArray(metric.Elem2Vol, 1, solve.zone2nelem[0], "Elem2Vol");
	metric.SumNorm(solve, FileContents, 1);
//	solve.Display3DArray(metric.Elem2DeltaS_xyz, 0, solve.zone2nelem[0], 3, "Elem2DeltaS_xyz");
//	solve.Display3DArray(metric.Elem2Center, 0, solve.zone2nelem[0], 3, "Elem2Center");
	cout << "Face2ElemCenter[i_zone][faceID][0] : "; cout << metric.Face2ElemCenter[0][2][0] <<endl;
	cout << "Face2ElemCenter[i_zone][faceID][1] : "; cout << metric.Face2ElemCenter[0][2][1] <<endl;


	cout << "Face2ElemCenter : " << endl;
	for(int face_i=0; face_i<solve.zone2nface[0]; face_i++)
	{
		cout << metric.Face2ElemCenter[0][face_i][0]; cout << " ; ";cout << metric.Face2ElemCenter[0][face_i][1] << endl;
	}




	/*
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
	solve.Display3DArray(solve.elem2node, izone, solve.zone2ncell[izone], 8, "elem2node");
	solve.Display1DArray(solve.zone2nnode, solve.nzone, "zone2nnode");
	solve.Display1DArray(solve.zone2nelem, solve.nzone, "zone2nelem");
	solve.Display2DArray(solve.zone2node, solve.nzone, solve.zone2nnode[izone], "zone2node");
	solve.Display2DArray(solve.zone2elem, solve.nzone, solve.zone2nelem[izone], "zone2elem");
	solve.Display2DArray(solve.nodeglobal2local, solve.nnode_g, solve.nzone, "nodeglobal2local");
	solve.Display2DArray(solve.elemglobal2local, solve.nelem_g, 2, "elemglobal2local");
	solve.Display3DArray(solve.belem2node, izone, solve.zone2nbelem[izone],4, "belem2node");
	solve.Display3DArray(solve.elem2elem, 0, solve.zone2ncell[0], 4, "elem2elem");
	solve.Display2DArray(solve.elem2vtk, solve.nzone, solve.zone2ncell[izone], "elem2vtk");
	solve.Display3DArray(solve.zone2coord, 0, solve.zone2nnode[0], 2, "zone2coord");
	\
	solve.Display3DArray(solve.face2node, 0, solve.zone2nface[0], 5, "face2node");
	solve.Display3DArray(solve.face2elem, 0, solve.zone2nface[0], 2, "face2elem");
	solve.Display3DArray(solve.face2fael, 0, solve.zone2nface[0], 2, "face2fael");
	solve.Display3DArray(solve.elem2face, 0, solve.zone2ncell[0], 6, "elem2face");
	////////////////////////////////////////////
	solve.MetricExemple(0);
	solve.SolverExemple();
	*/
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


	FileContents.write_tecplot(FileContents, "test2", p, Rho, u, v, w);

	return 0;
}
