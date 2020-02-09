#ifndef READER_H // include guard
#define READER_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <functional> // std::divides
using namespace std;

class Reader_c {
	public:
		ifstream file;
		unsigned linen;
		string line;
		const char* cline;
		unsigned ndime, nelem, npoint, nbc, nhalo, ncell;
		unsigned elem, point, bc, step;
		unsigned nelemlinen, npointlinen, bclinen;

		unsigned ** elem2node_nh,** elem2node;
		unsigned* elem2vtk_nh, * elem2vtk;
		double** coord;
		unsigned*** bc_elem2node;
		unsigned** bc_elem2vtk;
		unsigned bcnl, bc_nelem, bc_e2n_counter;
		unsigned* BoundIndex; unsigned* bc_nelemv;

		bool OpenFile(string filename);
		void read_file(string filename);
		unsigned Readcnst(const string& line, const string& tofind);
		void Fill_E2N_VTK(const char* cline);
		void Fill_BC_E2N_VTK(const char* cline, unsigned bc);
		double** Fill_coord(const char* cline);
};

#endif 