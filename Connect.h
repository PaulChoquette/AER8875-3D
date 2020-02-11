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
private:
	// Exemple d'entrer pour tester le code :
	int elem2zone[16] = { 0,0,1,1,0,0,1,1,2,2,3,3,2,2,3,3 };
	int nzone = 4;
	int ndime = 3;

public:
	// ================================================= INITIALIZATION ====================================================
	// (Mettre les variables a initialiser ici)




	// Integrer :
	int nelem_g, nnode_g;
	// 1D Array of Integrer :
	int* esup2; int* esup1; 
	int* zone2nelem; int* zone2nnode; int* checkzone; 
	// 2D Array of Integrer :
	int** elem2elem_g;
	int** zone2node; int** zone2elem; int** elemglobal2local; int** nodeglobal2local; int** zone2markelem; vector<vector<vector<int>>> zone2idmark;
	// 3D Array of Integrer :
	int*** elem2node; 
	// Double :

	// 1D Array of Double :

	// 2D Array of Double :
	
	// 3D Array of Double :
	double*** zone2coord;

	// Vector
	vector<vector<vector<int>>> belem2node;

	// =========================================== CONNECTIVITY FUNCTION MEMBERS ============================================
	// ELEMENTS CONNECTIVITY:
	void DataTransferExemple(Reader_c& read);
	void InitializeGlobal(Reader_c& read);
	void Node2Elements();
	void Node2Nodes();
	void Element2Elements();
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


};

#endif 