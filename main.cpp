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
using namespace std;

int main() {
	string File_Name = "naca0012_euler_65x65_O_1B.su2";
	cout << "Starting ..." << endl;



	Reader FileContents;

	FileContents.read_file(File_Name);

	return 0;
}
