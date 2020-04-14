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
class Solver_c;
class Reader_c {
	public:
		//Read Global Files
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

		//Local Only variables
		int zoneid, nzone, zlinen, znl, nzelem, z_nelem;
		int z_e2n_counter;
		int** pre_zelem2jelem;
		int* zelem2jelem, *zoneIndex;
		int*** z_elem2node;
		int** z_elem2vtk;
		int* z_nelemv;

		string* zone2tag;

		void read_file_local(string filename);
		void Fill_ZN_E2N_VTK(const char* cline, int zoneid);


		//Write SU2++ file
		void write_file(Reader_c& read, Solver_c& solve, int izone);
		void WriteAllZoneFile(Reader_c& read,Solver_c& solve );
		//Write Tecplot output
		void write_tecplot(Reader_c &FileContents, const char* out_filename,int world_id,void* p, void* Rho, void* u, void* v, void* w);
		void write_tecplot_ASCII(string FileName,double*p,double*rho,double*u,double*v,double*w);
		void write_tecplot_Dist(string what, string FileName,double* Cp, int wallFace, int wallNode, double** wallNode_coord, int** elem2face);







		///////////////////////////////////////////////////////////////////////////////
				string SimName;
				string su2FilePath,su2pFilePath;
				int Npartition;
		    double AoA;
		    double mach;
		    double gamma;
				string tempMethod;
				int Nstage;
				double cfl;
				string Smoothing;
				string spatMethod;
				int spatMethod_ordre;
				double grad;
				int iterMax;
				double convCrit;
				double AoA_i;
				double AoA_f;
				double Sref;
				double Cref;
				double xref;
				double yref;
				double zref;
				string coeffRef;

		///////////////////////////////////////////////////////////////////////////////

				ifstream file_2;
				int linen_2;
				string line_2;
				const char* cline_2;

				void computePrmt(string filename);
		    bool OpenFile_2(string filename);

				int inputInt(const string& line_2, const string& tofind);
				double inputDouble(const string& line_2, const string& tofind);
				string inputStr(const string& line_2, const string& tofind);
};

#endif
