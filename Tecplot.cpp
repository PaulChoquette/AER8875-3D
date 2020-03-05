#include <string>
#include <cctype>
#include <iostream>
#include <sstream>
#include <fstream>
#include "main.h"
#include "Reader.h"
#include "/home/cd/TECIO/include/TECXXX.h"


void Reader_c::write_tecplot(Reader_c &FileContents, const char* out_filename, double* p, double* Rho, double* u, double* v, double* w) {

	//Declare and fill X Y Z data arrays
	double* x = new double[FileContents.npoint];
	double* y = new double[FileContents.npoint];
	double* z = new double[FileContents.npoint];

	for (int i = 0; i < FileContents.npoint; i++) {
		x[i] = FileContents.coord[i][0];
		y[i] = FileContents.coord[i][1];
		if (FileContents.ndime > 2) {
			z[i] = FileContents.coord[i][2];
		}

	}

	//Initialize tecini variables and use tecini142 to initialize binary file
	INTEGER4 FileFormat = 0;
	INTEGER4 FileType = 0;
	INTEGER4 Debug = 1;
	INTEGER4 IsDouble = 1;
	char fname_temp[100];
	strcpy(fname_temp, "./");
	strcat(fname_temp, out_filename);
	strcat(fname_temp, ".plt");

	const char* fname = fname_temp;
	
	tecini142("NACA_0012", "X,Y,Z,P,Rho,U,V,W", fname , ".", &FileFormat, &FileType, &Debug, &IsDouble);

	//Initialize teczne 142 variables and use teczne142 to make header for each zone
	int Varlocation[] = { 1,1,1,0,0,0,0,0 }; // 1 is for node-centered, 0 is for cell-centered
	INTEGER4 ZoneType = 5;
	INTEGER4 NumFaces = 0; //Notused
	INTEGER4 ICellMax = 0; //Notused
	INTEGER4 JCellMax = 0; //Notused
	INTEGER4 KCellMax = 0; //Notused
	const double SolTime = 0.0;
	INTEGER4 StrandID = 0;
	INTEGER4 ParentZone = 0;
	INTEGER4 IsBlock = 1; //always 1
	INTEGER4 NumFaceConnec = 0;
	INTEGER4 FaceNeighMode = 2;
	INTEGER4 ShareConnec = 0;

	teczne142("main", &ZoneType, &FileContents.npoint, &FileContents.nelem, &NumFaces, &ICellMax, &JCellMax, &KCellMax, &SolTime
		, &StrandID, &ParentZone, &IsBlock, &NumFaceConnec, &FaceNeighMode, NULL, NULL, NULL, NULL
		, Varlocation, NULL, &ShareConnec);

	//Variable writing. One call of tecdat142 for each variable including coordinates
	tecdat142(&FileContents.npoint, x, &IsDouble);
	tecdat142(&FileContents.npoint, y, &IsDouble);
	tecdat142(&FileContents.npoint, z, &IsDouble);
	tecdat142(&FileContents.nelem, p, &IsDouble);
	tecdat142(&FileContents.nelem, Rho, &IsDouble);
	tecdat142(&FileContents.nelem, u, &IsDouble);
	tecdat142(&FileContents.nelem, v, &IsDouble);
	tecdat142(&FileContents.nelem, w, &IsDouble);

	//connectivity list vector writing and call tecnode142 to write connectivity
	INTEGER4 CC = FileContents.nelem * 8;
	INTEGER4* NData = new INTEGER4[CC];

	for (int i = 0; i < FileContents.nelem; i++) {
		for (int j = 0; j < 8; j++) {
			int k = j - floor(j / vtklookup[1][FileContents.elem2vtk[i]][1]);
			NData[i * 8 + j] = FileContents.elem2node[i][k] + 1;
		}
	}

	tecnode142(&CC, NData);

	tecend142();

}