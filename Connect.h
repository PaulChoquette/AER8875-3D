//#pragma once
#ifndef CONNECT_H // include guard
#define CONNECT_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <fstream>
#include "Reader.h"
using namespace std;
class Reader_c;

class Connect_c
{

public:
	// ================================================= INITIALIZATION ====================================================
	// (Mettre les variables a initialiser ici)


	// Integrer :
	int nelem, nnode, nzelem, nbound, ncell, nface;
	int nzone, ndime, mesup1, mpsup, neind;
	int nelem_g, nnode_g, ncell_g;
	// 1D Array of Integrer :
	int* esup2; int* esup1; int* psup1; int* psup2; int* lpoin;

	// 2D Array of Integrer :
	int** face2Nbr_of_node;
	int** vtk2nnofa; int** vtk2facevtk;

	int* elem2vtk; int* zelem2jelem; int* boundIndex; int* zoneIndex;
	// 2D Array of Integrer :
	int** elem2node; int*** vtk2lpofa; int** elem2elem; int** face2elem; int** face2node; int** face2fael; int** elem2face;
	// Double :

	// 1D Array of Double :

	// 2D Array of Double :

	// 3D Array of Double :
	double** coord;

	// Vector
	vector<vector<vector<int>>> zone2idmark;

	// =========================================== CONNECTIVITY FUNCTION MEMBERS ============================================
	// ELEMENTS CONNECTIVITY:
	void InitializeGlobal();

	// ZONE ELEMENTS CONNECTIVITY:
	void Node2Elements();
	void Node2Nodes();
	void Element2Elements();

	// FACES CONNECTIVITY:
	void Findnface();
	void Face2ElementsNodes();
	void Element2Faces();
	void ComputeLocalConnectivity();



	// OTHER FUNCTION MEMBER
	void Display1DArray(int V[], int size, string name);
	void Display2DArray(int* V[], int nLine, int nColomn, string name);
	void Display2DArray(double* V[], int nLine, int nColomn, string name);
	void Display3DArray(int** V[], int ix, int nLine, int nColomn, string name);
	void Display3DArray(double** V[], int ix, int nLine, int nColomn, string name);
};

#endif
