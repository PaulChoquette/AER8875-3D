
#include <string>
#include <cctype>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "omp.h"
#include <chrono>

#include "main.h"
#include "UI.h"
#include "Reader.h"
#include "METIS.h"
#include "Connect.h"
#include "Metric.h"
#include "Solver.h"
using namespace std;

int main() {
	cout << "Starting ..." << endl;
	string File_Name = "mesh_ONERAM6_inv_ffd.su2";

	Reader_c FileContents;

	FileContents.read_file(File_Name);

	Connect_c exemple;
	exemple.DataTransferExemple(FileContents);

	return 0;
}

