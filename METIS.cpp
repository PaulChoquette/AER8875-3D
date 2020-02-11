#include <iostream>
#include <algorithm>
#include <cmath>
#include <array>
#include <string>
#include <fstream>
#include <cstdio>
#include <time.h>
#include <omp.h>
#include <iomanip>

#include "METIS.h"
#include "Connect.h"
#include "Reader.h"

using namespace std;

void METIS_c::exemple(Reader_c& read) {
	cout << "exemple = " << read.ncell << endl;

}

void METIS_c::Mesh_Partition(Connect_c& mesh) {

}