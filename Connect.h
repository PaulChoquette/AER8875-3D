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
	// Exemple d'entrer pour tester le code :
	int elem2zone[16] = { 0,0,1,1,0,0,1,1,2,2,3,3,2,2,3,3 };

public:
	// ================================================= INITIALIZATION ====================================================
	// (Mettre les variables a initialiser ici)




	// Integrer :
	int nzone, ndime, mesup1, mpsup;
	int nelem_g, nnode_g, ncell_g;
	// 1D Array of Integrer :
	int* esup2; int* esup1; int* psup1; int* psup2; int* lpoin;
	int* zone2nelem; int* zone2nnode; int* checkzone;  int* zone2nbelem;
	// 2D Array of Integrer :
	int** vtk2nnofa; int** vtk2facevtk;
	int** elem2elem_g;
	int** zone2node; int** zone2elem; int** elemglobal2local; int** nodeglobal2local; int** zone2markelem; vector<vector<vector<int>>> zone2idmark;

	// 3D Array of Integrer :
	int*** elem2node; int*** vtk2lpofa; int*** belem2node;
	// Double :

	// 1D Array of Double :

	// 2D Array of Double :
	
	// 3D Array of Double :
	double*** zone2coord;

	// Vector
	

	// =========================================== CONNECTIVITY FUNCTION MEMBERS ============================================
	// ELEMENTS CONNECTIVITY:
	void ConnectExemple(Reader_c& read);
	void InitializeGlobal(Reader_c& read);
	void Node2Elements(Reader_c& read);
	void Node2Nodes(Reader_c& read);
	void Element2Elements(Reader_c& read);
	int GetnelemArNode(int inode);
	int GetielemArNode(int inode, int ielno);

	// FACES CONNECTIVITY:
	void Compute_nface();
	void Face2ElementsNodes();
	void Element2Faces();

	// ZONES CONNECTIVITY:
	void InitializeLocal();
	void Zone2nnode();
	void Zone2nelem();
	void Zone2Nodes();
	void Zone2Elements();
	void NodeGlobal2Local();
	void ElementGlobal2Local();
	void Zone2Coord(Reader_c& read);
	void Element2Nodes(Reader_c& read);
	void Zone2Zones();

	// OTHER FUNCTION MEMBER
	void Display1DArray(int V[], int size, string name);
	void Display2DArray(int* V[], int nLine, int nColomn, string name);
	void Display2DArray(double* V[], int nLine, int nColomn, string name);
	void Display3DArray(int** V[], int ix, int nLine, int nColomn, string name);
};

#endif 