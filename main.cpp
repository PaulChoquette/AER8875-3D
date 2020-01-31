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
	string File_Name = "mesh_ONERAM6_inv_ffd.su2";
<<<<<<< Updated upstream
	//Hello World!
=======
<<<<<<< HEAD
	//Hello World2!
=======

>>>>>>> parent of 27be833... Test_1
>>>>>>> Stashed changes
	Reader FileContents;

	FileContents.read_file(File_Name);

	return 0;
}
