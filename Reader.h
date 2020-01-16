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

class Reader {
	public:
		ifstream file;
		unsigned linen;
		string line;
		const char* cline;
		unsigned ndime, nelem, npoint;
		unsigned elem, point;
		unsigned nelemlinen, npointlinen;

		unsigned ** elem2node;
		unsigned* elem2vtk;
		double** coord;

		bool OpenFile(string filename);
		void read_file(string filename);
		unsigned Readcnst(const string& line);
		void Fill_E2N_VTK(const char* cline);
		double** Fill_coord(const char* cline)
};