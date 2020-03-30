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
	//int elem2zone[16] = { 0,0,1,1,0,0,1,1,2,2,3,3,2,2,3,3 };
	//int elem2zone[7] = { 0,0,0,0,0,0,0};
	//int elem2zone[2] = { 0, 1};
	int* elem2zone;

	//[1, 17, 32, 1, 2, 31, 32, 2, 3, 30, 31, 3, 4, 29, 30, 4, 28, 29, 1, 5, 17, 18, 1, 2, 5, 6, 2, 3, 6, 7, 3, 4, 7, 8, 4, 8, 27, 28, 5, 9, 18, 19, 5, 6, 9, 10, 6, 7, 10, 11, 7, 8, 11, 12, 8, 12, 26, 27, 9, 13, 19, 20, 9, 10, 13, 14, 10, 11, 14, 15, 11, 12, 15, 16, 12, 16, 25, 26, 13, 20, 21, 13, 14, 21, 22, 14, 15, 22, 23, 15, 16, 23, 24, 16, 24, 25, ]
	//[ 1, 17, 1, 2, 32, 2, 3, 31, 32, 3, 4, 30, 31, 4, 29, 30, 1, 5, 17, 18, 1, 2, 5, 6, 2, 3, 6, 7, 3, 4, 7, 8, 4, 8, 28, 29, 5, 9, 18, 19, 5, 6, 9, 10, 6, 7, 10, 11, 7, 8, 11, 12, 8, 12, 27, 28, 9, 13, 19, 20, 9, 10, 13, 14, 10, 11, 14, 15, 11, 12, 15, 16, 12, 16, 26, 27, 13, 20, 21, 13, 14, 21, 22, 14, 15, 22, 23, 15, 16, 23, 24, 25, 16, 24, 25, 26, ]
public:
	// ================================================= INITIALIZATION ====================================================
	// (Mettre les variables a initialiser ici)




	// Integrer :
	int nzone, ndime, mesup1, mpsup, neind; 
	int nelem_g, nnode_g, ncell_g;
	// 1D Array of Integrer :
	int* esup2; int* esup1; int* psup1; int* psup2; int* lpoin; int* eptr; int* eind;
	int* zone2nelem; int* zone2nnode; int* checkzone;  int* zone2nbelem; int* zone2nbound; int* zone2ncell; 
	// 2D Array of Integrer :
	int** vtk2nnofa; int** vtk2facevtk;
	int** elem2elem_g;
	int** zone2node; int** zone2elem; int** elemglobal2local; int** nodeglobal2local; int** zone2markelem;  int** zone2zone; int** zone2jzone; int** zone2ijzone;
	int** elem2vtk; int** zelem2jelem; int** zone2boundIndex; int** zone2zoneIndex; int** zone2esup1; int** zone2esup2; int** zone2psup1; int** zone2psup2; int** zone2lpoin; 
	// 3D Array of Integrer :
	int*** elem2node; int*** vtk2lpofa; int*** belem2node; 
	// Double :

	// 1D Array of Double :

	// 2D Array of Double :

	// 3D Array of Double :
	double*** zone2coord;

	// Vector
	vector<vector<vector<int>>> zone2idmark;

	// =========================================== CONNECTIVITY FUNCTION MEMBERS ============================================
	~Connect_c();
	// ELEMENTS CONNECTIVITY:
	void InitializeGlobal(Reader_c& read);
	void Node2Elements(Reader_c& read);
	void Node2Nodes(Reader_c& read);
	void Element2Elements(Reader_c& read);
	void METISElement2Nodes(Reader_c& read);
	void ComputeMETIS(int nzoneMETIS, int ncommon, Reader_c& read);
	int GetnelemArNode(int inode);
	int GetielemArNode(int inode, int ielno);
	void ComputeGlobalConnectivity(Reader_c& read);

	

	// ZONES CONNECTIVITY:
	void Findnzone();
	void InitializeLocal();
	void Zone2nnode();
	void Zone2nelem();
	void Zone2Nodes();
	void Zone2Elements();
	void NodeGlobal2Local();
	void ElementGlobal2Local();
	void Zone2Coord(Reader_c& read);
	void Zone2Bound(Reader_c& read);
	void Zone2jZone();
	void InitializeElem2Node(Reader_c& read);
	void Element2Nodes(Reader_c& read);
	void Zone2Zones();
	void ComputeZoneConnectivity(Reader_c& read);
    void ComputeElem2Zone();

	// OTHER FUNCTION MEMBER
	void Display1DArray(int V[], int size, string name);
	void Display2DArray(int* V[], int nLine, int nColomn, string name);
	void Display2DArray(double* V[], int nLine, int nColomn, string name);
	void Display3DArray(int** V[], int ix, int nLine, int nColomn, string name);
	void Display3DArray(double** V[], int ix, int nLine, int nColomn, string name);
};

#endif
