#include <string>
#include <cctype>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "omp.h"
#include <chrono>
#include "main.h"
using namespace std;
hhh

int main() {
	string File_Name = "mesh_ONERAM6_inv_ffd.su2";

	Reader FileContents;

	FileContents.read_file(File_Name);
	part(FileContents);
	return 0;
}
