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
#include "include/TECIO.h"
using namespace std;

int main() {
	string File_Name = "naca0012_euler_65x65x2_O_1B.su2";

	Reader_c FileContents;

	FileContents.read_file(File_Name);

	double* p[1];

	for (int i = 0; i < FileContents.nelem; i++) {
		p[0][i] = 0.;
	}

	int* e2n[8];
	for (int i = 0; i < 8; i++){
		e2n[i] = new int [FileContents.nelem];
	}


	//TECIO OP;
	//tecini142("lmao","p","lmao","/", (int32_t*) 0, (int32_t*) 2, (int32_t*) 0, (int32_t*) 1);
	//teczne142("main", (int32_t*) 5, OP.npoint, OP.nelem,0,0,0,0,0,0,0, (int32_t*) 1,0,0,0,0,0,NULL,0,NULL,NULL);
	//tecdat142(OP.nelem, p , (int32_t*) 1);
	//tecend142();

	return 0;
}
