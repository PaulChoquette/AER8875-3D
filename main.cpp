#include <string>
#include <cctype>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "omp.h"
#include <chrono>
#include "main.h"
#include "Reader.h"
#include "./include/TECXXX.h"
using namespace std;

int main() {
	string File_Name = "naca0012_euler_65x65x2_O_1B.su2";

	Reader_c FileContents;

	FileContents.read_file(File_Name);

	double* x = new double[FileContents.npoint];
	double* y = new double[FileContents.npoint];
	double* z = new double[FileContents.npoint];
	double* p = new double [FileContents.nelem];

	for (int i = 0; i < FileContents.npoint; i++) {
		x[i] = FileContents.coord[i][0];
		y[i] = FileContents.coord[i][1];
		if (FileContents.ndime > 2) {
			z[i] = FileContents.coord[i][2];
		}
		
	}
	for (int i = 0; i < FileContents.nelem; i++) {
		p[i] = 0.;
	}
	



	TECIO OP;


	

	INTEGER4 FileFormat = 0;
	INTEGER4 FileType = 0;
	INTEGER4 Debug = 1;
	INTEGER4 IsDouble = 1;

	tecini142("lmao","X,Y,Z,P","./lmao.plt",".", &FileFormat, &FileType, &Debug, &IsDouble);

	int Varlocation[] = { 1,1,1,0 };
	INTEGER4 ZoneType = 5;
	OP.npoint = FileContents.npoint;
	OP.nelem = FileContents.nelem;
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

	teczne142("main", &ZoneType,&OP.npoint, &OP.nelem,&NumFaces,&ICellMax,&JCellMax,&KCellMax,&SolTime
		,&StrandID,&ParentZone, &IsBlock,&NumFaceConnec, &FaceNeighMode, NULL, NULL, NULL,NULL
		, Varlocation,NULL, &ShareConnec);

	//coordinates and data writing
	tecdat142(&OP.npoint, x, &IsDouble);
	tecdat142(&OP.npoint, y, &IsDouble);
	tecdat142(&OP.npoint, z, &IsDouble);
	tecdat142(&OP.nelem, p, &IsDouble);

	//connectivity list loop writing
	INTEGER4 CC = OP.nelem * 8;
	INTEGER4* NData = new INTEGER4[CC];

	for (int i = 0; i < FileContents.nelem; i++) {
		for (int j = 0; j < 8; j++) {
			int k = j - floor(j/vtklookup[FileContents.elem2vtk[i]][1]);
			NData[i * 8 + j] = FileContents.elem2node[i][k]+1;
		}
	}

	tecnode142(&CC,NData);

	tecend142();

	return 0;
}
