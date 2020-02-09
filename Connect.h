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

	// 1D Array of Integrer :

	// 2D Array of Integrer :

	// 3D Array of Integrer :

	// Double :

	// 1D Array of Double :

	// 2D Array of Double :

	// 3D Array of Double :


	// =========================================== CONNECTIVITY FUNCTION MEMBERS ============================================
	// ELEMENTS CONNECTIVITY:
	void DataTransferExemple(Reader_c& read);
	void Node2Elements();
	void Node2Nodes();
	void Element2Elements();

	// FACES CONNECTIVITY:
	void Compute_nface();
	void Face2ElementsNodes();
	void Element2Faces();

	// ZONES CONNECTIVITY:
	void Zone2nnode();
	void Zone2nelem();
	void Zone2Nodes();
	void Zone2Elements();
	void NodeGlobal2Local();
	void ElementGlobal2Local();
	void Zone2Coord();
	void Element2Nodes();
	void Zone2Zones();


};

#endif 