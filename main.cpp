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

	double* p = new double[FileContents.nelem];
	double* Rho = new double[FileContents.nelem];
	double* u = new double[FileContents.nelem];
	double* v = new double[FileContents.nelem];
	double* w = new double[FileContents.nelem];

	for (int i = 0; i < FileContents.nelem; i++) {
		p[i] = 0.;
		Rho[i] = 1.;
		u[i] = 0.;
		v[i] = 0.;
		w[i] = 0.;
	}


	FileContents.write_tecplot(FileContents, "test2", p, Rho, u, v, w);

	return 0;
}
