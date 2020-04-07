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
class Connect_c;
class Reader_c {
	public:
		//Read Files
		string SimName;
		string su2FilePath;
		int Npartition;

		ifstream file;
		int linen;
		string line;
		const char* cline;
		int ndime, nelem, npoint, nbc, nhalo, ncell;
		int elem, point, bc, step;
		int nelemlinen, npointlinen, bclinen;

		int ** elem2node_nh,** elem2node;
		int* elem2vtk_nh, * elem2vtk;
		double** coord;
		int*** bc_elem2node;
		int** bc_elem2vtk;
		int bcnl, bc_nelem, bc_e2n_counter;
		int* BoundIndex; int* bc_nelemv;
		string* bound2tag;
		bool OpenFile(string filename);
		void read_file(string filename);
		int Readcnst(const string& line, const string& tofind);
		void Fill_E2N_VTK(const char* cline);
		void Fill_BC_E2N_VTK(const char* cline, int bc);
		double** Fill_coord(const char* cline);
		void check();
		//Write SU2++ file
		void write_file(string FileName, Reader_c& read, Connect_c& solve, int izone);
		void WriteAllZoneFile(string FileName, Reader_c& read,Connect_c& solve );
		//Write Tecplot output
		void write_tecplot(Reader_c &FileContents, const char* out_filename, double* p, double* Rho, double* u, double* v, double* w);
		void write_tecplot_METIS(string FileName, Reader_c & read, Connect_c& solve);
		void write_tecplot_Connectivity(int izone, string FileName, Reader_c & read, Connect_c& solve);
		void write_tecplot_OtherZone(int izone, int jzone, string FileName, Reader_c & read, Connect_c& solve);

		void computePrmt(string filename);
		bool OpenFile_2(string filename);
		ifstream file_2;
		int linen_2;
		string line_2;
		const char* cline_2;
		int inputInt(const string& line_2, const string& tofind);
		double inputDouble(const string& line_2, const string& tofind);
		string inputStr(const string& line_2, const string& tofind);
};

#endif
