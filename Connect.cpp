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

void Connect_c::InitializeGlobal(Reader_c& read) {
	nnode_g = read.npoint;
	nelem_g = read.nelem;
}

void Connect_c::Node2Elements() {

}

void Connect_c::Node2Nodes() {

}
void Connect_c::Element2Elements() {

}

int Connect_c::GetnelemArNode(int inode) {
	// Retourne le nombre d'element autour d'un noeud
	int nelno = esup2[inode] - esup2[inode - 1];
	return nelno;
}
int Connect_c::GetielemArNode(int inode, int ielno) {
	// Retourne l'index de l'element autour d'un noeud
	// ielno	: i element around that node you want to acces	(0 based integrer)
	int ielem = esup1[esup2[inode] + ielno];			
	return ielem;
}
// ================================================= FACES CONNECTIVITY ====================================================
void Connect_c::Compute_nface() {
	
}
void Connect_c::Face2ElementsNodes() {

}
void Connect_c::Element2Faces() {

}

// ================================================= ZONES CONNECTIVITY ====================================================
void Connect_c::InitializeLocal() {
	zone2nelem = new int[nzone];
	zone2nnode = new int[nzone];
	checkzone = new int[nzone];
	for (int izone = 0; izone < nzone; izone++) {
		zone2nelem[izone] = 0;
	}

	elemglobal2local = new int* [nelem_g];
	for (int ielem_g = 0; ielem_g < nelem_g; ielem_g++) {
		elemglobal2local[ielem_g] = new int[2];
	}

	nodeglobal2local = new int* [nnode_g];
	for (int inode_g = 0; inode_g < nnode_g; inode_g++) {
		nodeglobal2local[inode_g] = new int[nzone];
		for (int izone = 0; izone < nzone; izone++) {
			nodeglobal2local[inode_g][izone] = -1;
		}
	}

}
void Connect_c::Zone2nnode() {
	for (int izone = 0; izone < nzone; izone++) {
		checkzone[izone] = 0;
		zone2nnode[izone] = 0;
	}

	int ielem_g, izone, nelno;
	for (int inode = 0; inode < nnode_g; inode++) {
		nelno = GetnelemArNode(inode);

		for (int ielno = 0; ielno < nelno; ielno++) {
			ielem_g = GetielemArNode(inode, ielno);
			izone = elem2zone[ielem_g];
			//Verifie si on n’a pas deja compter le nœud dans la zone:
			if (checkzone[izone] == 0) {
				zone2nnode[izone] += 1;
				checkzone[izone] = 1;
			}
		}
		// Reset de checkzone
		for (int izone = 0; izone < nzone; izone++) {
			checkzone[izone] = 0;
		}
	}
}
void Connect_c::Zone2nelem() {
	for (int ielem = 0; ielem < nelem_g; ielem++) {
		// Jusqu'au ghostcell p-e
		int izone = elem2zone[ielem];
		zone2nelem[izone] += 1;

	}
}
void Connect_c::Zone2Nodes() {
	int* indexzone;
	indexzone = new int[nzone];
	for (int izone = 0; izone < nzone; izone++) {
		checkzone[izone] = 0;
		indexzone[izone] = 0;
	}

	int ielem_g, izone, nelno, idz;
	for (int inode = 0; inode < nnode_g; inode++) {
		nelno = GetnelemArNode(inode);

		for (int ielno = 0; ielno < nelno; ielno++) {
			ielem_g = GetielemArNode(inode, ielno);
			izone = elem2zone[ielem_g];
			//Verifie si on n’a pas deja compter le nœud dans la zone:
			if (checkzone[izone] == 0) {
				idz = indexzone[izone];
				zone2node[izone][idz] = inode;
				indexzone[izone] += 1;
				checkzone[izone] = 1;
			}
		}
		// Reset de checkzone
		for (int izone = 0; izone < nzone; izone++) {
			checkzone[izone] = 0;
		}
	}
	delete[] indexzone;
}
void Connect_c::Zone2Elements() {
	//Initialization:
	int* indexzone; indexzone = new int[nzone];
	zone2elem = new int*[nzone];
	for (int izone = 0; izone < nzone; izone++) {
		zone2elem[izone] = new int[zone2nelem[izone]];
		indexzone[izone] = 0;
	}

	for (int ielem = 0; ielem < nelem_g; ielem++) {
		int izone = elem2zone[ielem];
		int idz = indexzone[izone];
		zone2elem[izone][idz] = ielem;
		indexzone[izone] += 1;
	}
	delete[] indexzone;
}
void Connect_c::NodeGlobal2Local() {
	for (int izone = 0; izone < nzone; izone++) {
		for (int inode = 0; inode < zone2nnode[izone]; inode++) {
			int inode_g = zone2node[izone][inode];
			nodeglobal2local[inode_g][izone] = inode;
		}
	}
}
void Connect_c::ElementGlobal2Local() {
	for (int izone = 0; izone < nzone; izone++) {
		for (int ielem = 0; ielem < zone2nelem[izone]; ielem++) {
			int ielem_g = zone2elem[izone][ielem];
			elemglobal2local[ielem_g][0] = ielem;
			elemglobal2local[ielem_g][1] = izone;
		}
	}
}
void Connect_c::Zone2Coord(Reader_c& read) {
	zone2coord = new double** [nzone];
	for (int izone = 0; izone < nzone; izone++) {
		zone2coord[izone] = new double* [zone2nnode[izone]];
		for (int inode = 0; inode < zone2nnode[izone]; inode++) {
			zone2coord[izone][inode] = new double[ndime];
			int inode_g = zone2node[izone][inode];
			for (int iD = 0; iD < ndime; iD++) {
				zone2coord[izone][inode][iD] = read.coord[inode_g][iD];
			}
		}
	}

}
void Connect_c::Element2Nodes(Reader_c& read) {
	
	// Initialization:
	int* indexzone; 
	int** izone2jzone;
	indexzone = new int[nzone];
	zone2markelem = new int*[nzone];
	zone2idmark.resize(nzone);
	izone2jzone = new int* [nzone];
	elem2node = new int**[nzone];

	for (int izone = 0; izone < nzone; izone++) {
		indexzone[izone] = 0;
		zone2markelem[izone] = new int[nzone - 1];;
		zone2idmark[izone].resize(nzone - 1);
		izone2jzone[izone] = new int[nzone];

		int idz = 0;
		for (int jzone = 0; jzone < nzone; jzone++) {
			if (izone != jzone) {
				izone2jzone[izone][jzone] = idz;
				idz += 1;
			}
			else {
				izone2jzone[izone][jzone] = -1;
			}
		}

		elem2node[izone] = new int* [zone2nelem[izone]];
		for (int ielem = 0; ielem < zone2nelem[izone]; ielem++) {
			int nnoel = 999999;
			elem2node[izone][ielem] = new int[nnoel];
		}

		for (int ijzone = 0; ijzone < nzone - 1; ijzone++) {
			zone2markelem[izone][ijzone] = 0;
		}
	}

	// Debut du calcul:
	for (int ielem_g = 0; ielem_g < nelem_g; ielem_g++) {
		int ielem_z = elemglobal2local[ielem_g][0];
		int izone = elemglobal2local[ielem_g][1];
		int nnoel = 999999;
		int nfael = 999999;

		for (int inoel = 0; inoel < nnoel; inoel++) {
			int inode_g = read.elem2node[ielem_g][inoel];
			int inode_z = nodeglobal2local[inode_g][izone];
			elem2node[izone][ielem_z][inoel] = inode_z;
		}

		for (int ifael = 0; ifael < nfael; ifael++) {
			int jelem_g = elem2elem_g[ielem_g][ifael];

			// Si jelem nest pas un ghostcell:
			if (jelem_g <= nelem_g) {
				int jzone = elem2zone[jelem_g];
				int jelem_z = elemglobal2local[jelem_g][0];
				
				// Detection d'une frontiere auvec une autre zone:
				if (jzone != izone) {
					int idz = indexzone[izone];
					int nnofa = 999999;
					belem2node[izone][idz].resize(nnofa + 2);
					belem2node[izone][idz][0] = 999999; // vtk de la face
					belem2node[izone][idz][nnofa + 1] = jelem_z;
					int ijzone = izone2jzone[izone][jzone];
					zone2markelem[izone][ijzone] += 1;
					zone2idmark[izone][ijzone].push_back(idz);

					for (int inofa = 0; inofa < nnofa; inofa++) {
						int inoel = 999999; // Getlpofa()
						int inode_g = read.elem2node[ielem_g][inoel];
						int inode_z = nodeglobal2local[inode_g][izone];
						belem2node[izone][idz][inofa + 1] = inode_z;
					}
					indexzone[izone] += 1;
				}
			}
		}
	}	
}
void Connect_c::Zone2Zones() {

}