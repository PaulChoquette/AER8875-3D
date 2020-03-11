#include <string>
#include <cctype>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "omp.h"
#include <chrono>
#include <time.h>

#include "main.h"
#include "UI.h"
#include "Reader.h"
#include "Connect.h"
#include "Metric.h"
#include "Solver.h"


//#include "./include/TECXXX.h"
#include <metis.h>

using namespace std;

void TEST_Exemple2D() {
	cout << "\n============================== TEST 2D CONNECTIVITY ============================== " << endl;
    string File_Name = "block.su2";
    Reader_c File;
	File.read_file(File_Name);
	Solver_c test2D;
	test2D.ComputeGlobalConnectivity(File);
	test2D.ComputeElem2Zone();
	test2D.ComputeZoneConnectivity(File);
	File.WriteAllZoneFile("2DZone", File, test2D);
	File.write_tecplot_METIS("METISZONE_Block.dat", File,test2D);
	int izone = 0; // iZone to Display
	// test2D.Display1DArray(test2D.zone2nnode, test2D.nzone, "zone2nnode");
	// test2D.Display1DArray(test2D.zone2nelem, test2D.nzone, "zone2nelem");
	// test2D.Display2DArray(test2D.zone2node, test2D.nzone, test2D.zone2nnode[izone], "zone2node");
	// test2D.Display2DArray(test2D.zone2elem, test2D.nzone, test2D.zone2nelem[izone], "zone2elem");
	test2D.Display2DArray(test2D.nodeglobal2local, test2D.nnode_g, test2D.nzone, "nodeglobal2local");
	test2D.Display2DArray(test2D.elemglobal2local, test2D.nelem_g, 2, "elemglobal2local");
	test2D.Display3DArray(test2D.zone2coord, izone, test2D.zone2nnode[0], 2, "zone2coord");
	test2D.Display2DArray(test2D.elem2vtk, test2D.nzone, test2D.zone2ncell[izone], "elem2vtk (all zone)");
	test2D.Display3DArray(test2D.elem2node, izone, test2D.zone2ncell[izone], 4, "elem2node");
	}

void TEST_NACA0012() {
	cout << "\n============================== TEST NACA0012 ============================== " << endl;
    string File_Name = "naca0012_euler_65x65x2_O_1B.su2";
    Reader_c File0012;
	File0012.read_file(File_Name);
	Solver_c NACA0012;
	NACA0012.ComputeGlobalConnectivity(File0012);
	NACA0012.ComputeMETIS(4, 4, File0012);
	NACA0012.ComputeZoneConnectivity(File0012);
	File0012.WriteAllZoneFile("NACA0012Zone", File0012, NACA0012);
	File0012.write_tecplot_METIS("METISZONE_NACA0012.dat", File0012,NACA0012);	
	}

int main() {
	cout << "Starting ..." << endl;

    TEST_Exemple2D();
	TEST_NACA0012();

	//string File_Name = "block.su2";
	//string File_Name = "test_justSquare.su2";
	//string File_Name = "test.su2";
	//string File_Name = "2cube.su2";
	string File_Name = "naca0012_euler_65x65x2_O_1B.su2";
//	string File_Name = "mesh_ONERAM6_inv_ffd.su2";
	Reader_c FileContents;
	FileContents.read_file(File_Name);
	Solver_c solve;
	solve.ComputeGlobalConnectivity(FileContents);
	solve.ComputeMETIS(2, 4, FileContents);
	solve.ComputeZoneConnectivity(FileContents);
	FileContents.WriteAllZoneFile("Zone", FileContents, solve);
	FileContents.write_tecplot_METIS("METISZONE_ONERA.dat", FileContents,solve);

	// =================================== EXECUTABLE 2 ====================================================

	solve.ComputeLocalConnectivity();


	// METRIC
	//cout << "\nMetriques ..." << endl;
	//Metric_c metric;
	//metric.Compute(solve, FileContents);
//	solve.Display3DArray(metric.Face2Norm, 0, solve.zone2nface[0], 3, "Face2Norm");
//	solve.Display2DArray(metric.Face2Area, 1, solve.zone2nface[0], "Face2Area");
//	solve.Display2DArray(metric.Elem2Vol, 1, solve.zone2nelem[0], "Elem2Vol");
	//metric.SumNorm(solve, FileContents, 1);
//	solve.Display3DArray(metric.Elem2DeltaS_xyz, 0, solve.zone2nelem[0], 3, "Elem2DeltaS_xyz");
//	solve.Display3DArray(metric.Elem2Center, 0, solve.zone2nelem[0], 3, "Elem2Center");
	// cout << "Face2ElemCenter[i_zone][faceID][0] : "; cout << metric.Face2ElemCenter[0][2][0] <<endl;
	// cout << "Face2ElemCenter[i_zone][faceID][1] : "; cout << metric.Face2ElemCenter[0][2][1] <<endl;

    //TEST_Exemple2D();
	
   
	

	//cout << "Face2ElemCenter : " << endl;
	// for(int face_i=0; face_i<solve.zone2nface[0]; face_i++)
	// {
	// 	cout << metric.Face2ElemCenter[0][face_i][0]; cout << " ; ";cout << metric.Face2ElemCenter[0][face_i][1] << endl;
	// }

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
	// solve.Display1DArray(solve.zone2nnode, solve.nzone, "zone2nnode");
	// solve.Display1DArray(solve.zone2nelem, solve.nzone, "zone2nelem");
	// solve.Display2DArray(solve.zone2node, solve.nzone, solve.zone2nnode[izone], "zone2node");
	// solve.Display2DArray(solve.zone2elem, solve.nzone, solve.zone2nelem[izone], "zone2elem");
	// solve.Display2DArray(solve.nodeglobal2local, solve.nnode_g, solve.nzone, "nodeglobal2local");
	// solve.Display2DArray(solve.elemglobal2local, solve.nelem_g, 2, "elemglobal2local");
	// solve.Display3DArray(solve.belem2node, izone, solve.zone2nbelem[izone],6, "belem2node");
	// solve.Display3DArray(solve.elem2node, izone, solve.zone2ncell[izone], 6, "elem2node");
	//solve.Display3DArray(solve.elem2elem, izone, solve.zone2ncell[izone], 4, "elem2elem");
	//solve.Display2DArray(solve.elem2vtk, solve.nzone, solve.zone2ncell[izone], "elem2vtk");
	//solve.Display3DArray(solve.zone2coord, 0, solve.zone2nnode[0], 2, "zone2coord");
	//solve.Display3DArray(solve.face2node, 0, solve.zone2nface[0], 5, "face2node");
	// solve.Display3DArray(solve.face2elem, 0, solve.zone2nface[0], 2, "face2elem");
	// solve.Display3DArray(solve.face2fael, 0, solve.zone2nface[0], 2, "face2fael");
	//solve.Display3DArray(solve.elem2face, 0, solve.zone2ncell[0], 6, "elem2face");
	cout << "**************\nEnd\n**************\n";
	


//FileContents.write_tecplot(FileContents, "test2", p, Rho, u, v, w);

	return 0;
}
