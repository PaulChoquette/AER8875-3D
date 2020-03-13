// Execute with CMD :  g++ *.cpp -lgomp -o out -I /home/charles/Documents/metis-5.1.0/include -L /home/charles/Documents/metis-5.1.0/lib -lmetis


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
#include "Reader.h"
#include "Connect.h"

//#include "./include/TECXXX.h"
#include <metis.h>

using namespace std;

void TEST_Exemple2D() {
	cout << "\n============================== TEST 2D CONNECTIVITY ============================== " << endl;
    string File_Name = "block.su2";
    Reader_c File;
	Connect_c test2D;
	File.read_file(File_Name);
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
	cout << "\n================================== NACA0012 ================================== " << endl;
    string File_Name = "naca0012_euler_65x65x2_O_1B.su2";
    Reader_c File0012;
	Connect_c NACA0012;
	File0012.read_file(File_Name);
	NACA0012.ComputeGlobalConnectivity(File0012);
	NACA0012.ComputeMETIS(4, 4, File0012);
	NACA0012.ComputeZoneConnectivity(File0012);
	File0012.WriteAllZoneFile("NACA0012Zone", File0012, NACA0012);
	File0012.write_tecplot_METIS("METISZONE_NACA0012.dat", File0012,NACA0012);	
	File0012.write_tecplot_OtherZone(0, 1 ,"Tecplot_Connectivity_NACA0012_0.dat", File0012,NACA0012);	
	File0012.write_tecplot_OtherZone(1, 0, "Tecplot_Connectivity_NACA0012_1.dat", File0012,NACA0012);	
	}

void TEST_ONERA() {
	cout << "\n==================================== ONERA ==================================== " << endl;
    string File_Name = "mesh_ONERAM6_inv_ffd.su2";
    Reader_c FileONERA;
	Connect_c ONERA;
	FileONERA.read_file(File_Name);
	ONERA.ComputeGlobalConnectivity(FileONERA);
	ONERA.ComputeMETIS(4, 3, FileONERA);
	ONERA.ComputeZoneConnectivity(FileONERA);
	FileONERA.WriteAllZoneFile("ONERAZone", FileONERA, ONERA);
	FileONERA.write_tecplot_METIS("METISZONE_ONERA.dat", FileONERA,ONERA);	
	FileONERA.write_tecplot_OtherZone(0, 1 ,"Tecplot_Connectivity_ONERA_0_1.dat", FileONERA,ONERA);	
	FileONERA.write_tecplot_OtherZone(1, 0, "Tecplot_Connectivity_ONERA_1_0.dat", FileONERA,ONERA);

	FileONERA.write_tecplot_OtherZone(0, 2 ,"Tecplot_Connectivity_ONERA_0_2.dat", FileONERA,ONERA);	
	FileONERA.write_tecplot_OtherZone(2, 0, "Tecplot_Connectivity_ONERA_2_0.dat", FileONERA,ONERA);
	
	FileONERA.write_tecplot_OtherZone(0, 3 ,"Tecplot_Connectivity_ONERA_1_2.dat", FileONERA,ONERA);	
	FileONERA.write_tecplot_OtherZone(3, 0, "Tecplot_Connectivity_ONERA_2_1.dat", FileONERA,ONERA);

	FileONERA.write_tecplot_OtherZone(1, 2 ,"Tecplot_Connectivity_ONERA_1_2.dat", FileONERA,ONERA);	
	FileONERA.write_tecplot_OtherZone(2, 1, "Tecplot_Connectivity_ONERA_2_1.dat", FileONERA,ONERA);

	FileONERA.write_tecplot_OtherZone(1, 3 ,"Tecplot_Connectivity_ONERA_1_2.dat", FileONERA,ONERA);	
	FileONERA.write_tecplot_OtherZone(3, 1, "Tecplot_Connectivity_ONERA_2_1.dat", FileONERA,ONERA);

	FileONERA.write_tecplot_OtherZone(2, 3 ,"Tecplot_Connectivity_ONERA_2_3.dat", FileONERA,ONERA);	
	FileONERA.write_tecplot_OtherZone(3, 2, "Tecplot_Connectivity_ONERA_3_2.dat", FileONERA,ONERA);
	}

int main() {
	cout << "==================================== NOICE STARTING ====================================" << endl;

    TEST_Exemple2D();
	TEST_NACA0012();
	TEST_ONERA();
	

	 // DISPLAY OF GLOBAL
	//mesh.Display2DArray(FileContents.elem2node, FileContents.ncell, 8, "elem2node_g");
	//mesh.Display1DArray(mesh.esup1, mesh.mesup1, "esup1");
	//mesh.Display1DArray(mesh.esup2, mesh.nnode_g+1, "esup2");
	//mesh.Display1DArray(mesh.psup1, mesh.mpsup, "psup1");
	//mesh.Display1DArray(mesh.psup2, mesh.nnode_g+1, "psup2");
	//mesh.Display2DArray(mesh.elem2elem_g, mesh.ncell_g, 6, "elem2elem_g");

	//// DISPLAY OF ZONE
	//int izone = 0; // iZone to Display
	// mesh.Display1DArray(mesh.zone2nnode, mesh.nzone, "zone2nnode");
	// mesh.Display1DArray(mesh.zone2nelem, mesh.nzone, "zone2nelem");
	// mesh.Display2DArray(mesh.zone2node, mesh.nzone, mesh.zone2nnode[izone], "zone2node");
	// mesh.Display2DArray(mesh.zone2elem, mesh.nzone, mesh.zone2nelem[izone], "zone2elem");
	// mesh.Display2DArray(mesh.nodeglobal2local, mesh.nnode_g, mesh.nzone, "nodeglobal2local");
	// mesh.Display2DArray(mesh.elemglobal2local, mesh.nelem_g, 2, "elemglobal2local");
	// mesh.Display3DArray(mesh.belem2node, izone, mesh.zone2nbelem[izone],6, "belem2node");
	// mesh.Display3DArray(mesh.elem2node, izone, mesh.zone2ncell[izone], 6, "elem2node");
	//mesh.Display3DArray(mesh.elem2elem, izone, mesh.zone2ncell[izone], 4, "elem2elem");
	//mesh.Display2DArray(mesh.elem2vtk, mesh.nzone, mesh.zone2ncell[izone], "elem2vtk");
	//mesh.Display3DArray(mesh.zone2coord, 0, mesh.zone2nnode[0], 2, "zone2coord");

	cout << "**************\nEnd\n**************\n";
	

	return 0;
}
