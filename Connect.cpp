#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <array>
#include <string>
#include <fstream>
#include <cstdio>
#include <time.h>
#include <omp.h>
#include <iomanip>

#include "Connect.h"
#include "Reader.h"
using namespace std;


// ================================================= ELEMENTS CONNECTIVITY ====================================================
void Connect_c::DataTransferExemple(Reader_c& read) {
	// Test pour tester le transfere de donner entre les class
	cout << "test de transfert de donner entre class: nelem = " << read.nelem;
}

void Connect_c::Node2Elements() {

}

void Connect_c::Node2Nodes() {

}
void Connect_c::Element2Elements() {

}

// ================================================= FACES CONNECTIVITY ====================================================
void Connect_c::Compute_nface() {
	
}
void Connect_c::Face2ElementsNodes() {

}
void Connect_c::Element2Faces() {

}

// ================================================= ZONES CONNECTIVITY ====================================================
void Connect_c::Zone2nnode() {

}
void Connect_c::Zone2nelem() {

}
void Connect_c::Zone2Nodes() {

}
void Connect_c::Zone2Elements() {

}
void Connect_c::NodeGlobal2Local() {

}
void Connect_c::ElementGlobal2Local() {

}
void Connect_c::Zone2Coord() {

}
void Connect_c::Element2Nodes() {

}
void Connect_c::Zone2Zones() {

}