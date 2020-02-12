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
#include "main.h"
using namespace std;


// ================================================= ELEMENTS CONNECTIVITY ====================================================
void Connect_c::ConnectExemple(Reader_c& read) {
	// Test pour tester le transfere de donner entre les clase
	nzone = 4;
	
}

void Connect_c::InitializeGlobal(Reader_c& read) {
	nnode_g = read.npoint;
	nelem_g = read.nelem;
	ndime = read.ndime;
	ncell_g = read.ncell;

	vtk2nnofa = new int* [15];
	vtk2nnofa[3] = new int[2]{ 2,2};
	vtk2nnofa[5] = new int[3]{ 2,2,2 };
	vtk2nnofa[9] = new int[4]{ 2,2,2,2 };
	vtk2nnofa[12] = new int[6]{ 4,4,4,4,4,4 };
	
	vtk2lpofa = new int** [15];
	vtk2lpofa[3] = new int* [2];
	vtk2lpofa[3][0] = new int[2]{ 0,1 };
	vtk2lpofa[3][1] = new int[2]{ 1,0 };

	vtk2lpofa[5] = new int* [2];
	vtk2lpofa[5][0] = new int[3]{ 0,1,2 };
	vtk2lpofa[5][1] = new int[3]{ 1,2,0 };

	vtk2lpofa[9] = new int* [2];
	vtk2lpofa[9][0] = new int[4]{ 0,1,2,3 };
	vtk2lpofa[9][1] = new int[4]{ 1,2,3,0 };

	vtk2lpofa[12] = new int* [4];
	vtk2lpofa[12][0] = new int[6]{ 0,4,0,0,1,2 };
	vtk2lpofa[12][1] = new int[6]{ 1,5,3,4,5,6 };
	vtk2lpofa[12][2] = new int[6]{ 2,6,7,5,6,7 };
	vtk2lpofa[12][3] = new int[6]{ 3,7,4,1,2,3 };

	

	vtk2facevtk = new int* [15];
	vtk2facevtk[3] = new int[4]{ 3,3,3};
	vtk2facevtk[9] = new int[4]{ 3,3,3,3 };
	vtk2facevtk[12] = new int[6]{ 9,9,9,9,9,9 };

}
void Connect_c::Node2Elements(Reader_c& read) {
	mesup1 = 0;
	for (int ielem = 0; ielem < ncell_g; ielem++) {
		int vtk = read.elem2vtk[ielem];
		mesup1 += vtklookup[vtk][1];
	}

	// Initializing esup and other parameters :
	esup2 = new int[nnode_g + 1]();
	esup1 = new int[mesup1]();
	int ipoil;
	//
	for (int ielem = 0; ielem < ncell_g; ielem++) {
		int vtk = read.elem2vtk[ielem];
		int nnode = vtklookup[vtk][1];
		for (int inode = 0; inode < nnode; inode++) {
			ipoil = read.elem2node[ielem][inode]+1; // elem2node en 1 base 999999999999999999999
			esup2[ipoil - 1] = esup2[ipoil - 1] + 1;					
		}
	}
	for (int ipoin = 2; ipoin <= nnode_g+1; ipoin++) {
		esup2[ipoin - 1] = esup2[ipoin - 1] + esup2[ipoin - 1 - 1];
	}
	int ipoin, istor;
	for (int ielem = 0; ielem < ncell_g; ielem++) {
		int vtk = read.elem2vtk[ielem];
		int nnode = vtklookup[vtk][1];

		for (int inode = 0; inode <= nnode - 1; inode++) {
			ipoin = read.elem2node[ielem][inode]; // elem2node en 1 base 999999999999999999999
			istor = esup2[ipoin - 1] + 1; 
			esup2[ipoin - 1] = istor;
			esup1[istor - 1] = ielem + 1;
		}
	}
	for (int ipoin = nnode_g + 1; ipoin >= 2; --ipoin) {
		esup2[ipoin - 1] = esup2[ipoin - 1 - 1];
	}
	esup2[0] = 0;

	

}
void Connect_c::Node2Nodes(Reader_c& read) {
	lpoin = new int[nnode_g]();
	psup2 = new int[nnode_g + 1]();
	int ielem, istor{ 0 }, jpoin;


	for (int ipoin = 1; ipoin <= nnode_g; ipoin++) {
		for (int iesup = esup2[ipoin - 1] + 1; iesup <= esup2[ipoin + 1 - 1]; iesup++) {
			ielem = esup1[iesup - 1] - 1;															
			int vtk = read.elem2vtk[ielem];
			int nnode = vtklookup[vtk][1];
			for (int inode = 1; inode <= nnode; inode++) {
				jpoin = read.elem2node[ielem][inode - 1]; // elem2node en 1 base 999999999999999999999
				if (jpoin != ipoin && lpoin[jpoin - 1] != ipoin) {
					istor = istor + 1;
					lpoin[jpoin - 1] = ipoin;
				}
			}
		}
		psup2[ipoin + 1 - 1] = istor;
	}

	//delete[] lpoin;
	lpoin = new int[nnode_g]();
	istor = 0; jpoin = 0;
	mpsup = psup2[nnode_g]; 
	psup1 = new int[mpsup]();
	for (int ipoin = 1; ipoin <= nnode_g; ipoin++) {
		for (int iesup = esup2[ipoin - 1] + 1; iesup <= esup2[ipoin + 1 - 1]; iesup++) {
			ielem = esup1[iesup - 1] - 1;																		// 1 based
			int vtk = read.elem2vtk[ielem];
			int nnode = vtklookup[vtk][1];
			for (int inode = 1; inode <= nnode; inode++) {
				jpoin = read.elem2node[ielem][inode - 1]; // elem2node en 1 base 999999999999999999999
				if (jpoin != ipoin && lpoin[jpoin - 1] != ipoin) {
					istor = istor + 1;
					psup1[istor - 1] = jpoin;
					lpoin[jpoin - 1] = ipoin;
				}
			}
		}
	}
}
void Connect_c::Element2Elements(Reader_c& read) {
	elem2elem_g = new int* [ncell_g]();
	int nnofa, ilpofa, ipoin, jelem, nnofj, icoun, jpoin;
	vector<int> lhelp;
	for (int ielem = 0; ielem < ncell_g; ielem++) {
		int vtk = read.elem2vtk[ielem];
		int nfael = vtklookup[vtk][0];
		elem2elem_g[ielem] = new int[nfael]();
		for (int ifael = 0; ifael <= nfael - 1; ifael++) {
			nnofa = vtk2nnofa[vtk][ifael];
			lhelp.resize(nnofa);															// vecteur a enlever

			for (int inofa = 0; inofa <= nnofa - 1; inofa++) {
			
				ilpofa = vtk2lpofa[vtk][inofa][ifael]+1;
				lhelp[inofa] = read.elem2node[ielem][ilpofa - 1];	// elem2node en 1 base 999999999999999999999
				lpoin[lhelp[inofa] - 1] = 1;
			}
			ipoin = lhelp[0] - 1;
			for (int istor = esup2[ipoin]; istor <= esup2[ipoin + 1] - 1; istor++) {
				jelem = esup1[istor];
				int jvtk = read.elem2vtk[jelem-1];
				int nfael = vtklookup[jvtk][0];
				if (jelem != ielem + 1) {
					for (int jfael = 0; jfael <= nfael - 1; jfael++) {
						//nnofj = Getlnofa(jelem - 1, jfael);  // lnofa[jfael];
						nnofj = vtk2nnofa[jvtk][jfael];
						if (nnofj == nnofa) {
							icoun = 0;
							for (int jnofa = 0; jnofa <= nnofa - 1; jnofa++) {
								//jpoin = read.elem2node[jelem - 1][Getlpofa(jelem - 1, jnofa, jfael) + 1 - 1];
								jpoin = read.elem2node[jelem - 1][vtk2lpofa[jvtk][jnofa][jfael] + 1 - 1];
								icoun = icoun + lpoin[jpoin - 1];
							}
							if (icoun == nnofa) {
								elem2elem_g[ielem][ifael] = jelem - 1;
							}
						}
					}
				}
			}
			int vtk2 = read.elem2vtk[ielem];
			nnofa = vtk2nnofa[vtk2][ifael];
			for (int inova = 0; inova <= nnofa - 1; inova++) {
				lpoin[lhelp[inova] - 1] = 0;
			}
		}
	}

	//delete[] lpoin;
}
int Connect_c::GetnelemArNode(int inode) {
	// Retourne le nombre d'element autour d'un noeud
	int nelno = esup2[inode+1] - esup2[inode - 1+1];
	return nelno;
}
int Connect_c::GetielemArNode(int inode, int ielno) {
	// Retourne l'index de l'element autour d'un noeud
	// ielno	: i element around that node you want to acces	(0 based integrer)
	int ielem = esup1[esup2[inode] + ielno] -1;			
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
			if (ielem_g < nelem_g && checkzone[izone] == 0) {
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
	zone2node = new int*[nzone];
	for (int izone = 0; izone < nzone; izone++) {
		checkzone[izone] = 0;
		indexzone[izone] = 0;
		zone2node[izone] = new int[zone2nnode[izone]];
	}

	int ielem_g, izone, nelno, idz;
	for (int inode = 0; inode < nnode_g; inode++) {
		nelno = GetnelemArNode(inode);

		for (int ielno = 0; ielno < nelno; ielno++) {
			ielem_g = GetielemArNode(inode, ielno);
			izone = elem2zone[ielem_g];
			//Verifie si on n’a pas deja compter le nœud dans la zone:
			if (ielem_g < nelem_g && checkzone[izone] == 0 ) {
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
	
	zone2nbelem = new int[nzone];
	elem2node = new int** [nzone];
	for (int izone = 0; izone < nzone; izone++) {
		zone2nbelem[izone] = 0;

		elem2node[izone] = new int* [zone2nelem[izone]];
		for (int ielem = 0; ielem < zone2nelem[izone]; ielem++) {
			int vtk = read.elem2vtk[ielem];
			int nnoel = vtklookup[vtk][1];
			elem2node[izone][ielem] = new int[nnoel];
		}
	}


	// Debut du calcul de grandeur:
	for (int ielem_g = 0; ielem_g < nelem_g; ielem_g++) {
		int ielem_z = elemglobal2local[ielem_g][0];
		int izone = elemglobal2local[ielem_g][1];
		int vtk = read.elem2vtk[ielem_g];
		int nnoel = vtklookup[vtk][1];
		int nfael = vtklookup[vtk][0];

		for (int inoel = 0; inoel < nnoel; inoel++) {
			int inode_g = read.elem2node[ielem_g][inoel] - 1; //  elem2node est en 1 based 99999999999999
			int inode_z = nodeglobal2local[inode_g][izone];
			elem2node[izone][ielem_z][inoel] = inode_z;
		}
		for (int ifael = 0; ifael < nfael; ifael++) {
			int jelem_g = elem2elem_g[ielem_g][ifael]; 

			// Si jelem nest pas un ghostcell:
			if (jelem_g < nelem_g) {
				int jzone = elem2zone[jelem_g];
				int jelem_z = elemglobal2local[jelem_g][0];

				// Detection d'une frontiere auvec une autre zone:
				if (jzone != izone) {
					zone2nbelem[izone] += 1;
				}
			}
		}
	}



	// Initialization:
	int** izone2jzone;
	int* indexzone;
	indexzone = new int[nzone];
	belem2node = new int**[nzone];
	zone2markelem = new int*[nzone];
	zone2idmark.resize(nzone);
	izone2jzone = new int* [nzone];

	for (int izone = 0; izone < nzone; izone++) {
		zone2markelem[izone] = new int[nzone - 1];;
		zone2idmark[izone].resize(nzone - 1);
		izone2jzone[izone] = new int[nzone];
		belem2node[izone] = new int*[zone2nbelem[izone]];
		indexzone[izone] = 0;
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


		for (int ijzone = 0; ijzone < nzone - 1; ijzone++) {
			zone2markelem[izone][ijzone] = 0;
		}


	}






	// Debut du calcul:
	for (int ielem_g = 0; ielem_g < nelem_g; ielem_g++) {
		int ielem_z = elemglobal2local[ielem_g][0];
		int izone = elemglobal2local[ielem_g][1];
		int vtk = read.elem2vtk[ielem_g];
		int nnoel = vtklookup[vtk][1];
		int nfael = vtklookup[vtk][0];

		for (int inoel = 0; inoel < nnoel; inoel++) {
			int inode_g = read.elem2node[ielem_g][inoel] - 1; //  elem2node est en 1 based 99999999999999
			int inode_z = nodeglobal2local[inode_g][izone];
			elem2node[izone][ielem_z][inoel] = inode_z;
		}

		for (int ifael = 0; ifael < nfael; ifael++) {
			int jelem_g = elem2elem_g[ielem_g][ifael];

			// Si jelem nest pas un ghostcell:
			if (jelem_g < nelem_g) {
				int jzone = elem2zone[jelem_g];
				int jelem_z = elemglobal2local[jelem_g][0];
				
				// Detection d'une frontiere auvec une autre zone:
				if (jzone != izone) {
					int idz = indexzone[izone];

					int nnofa = vtk2nnofa[vtk][ifael];
					belem2node[izone][idz] = new int[nnofa + 2];
					belem2node[izone][idz][0] = vtk2facevtk[vtk][ifael]; // vtk de la face
					belem2node[izone][idz][nnofa + 1] = jelem_z;
					int ijzone = izone2jzone[izone][jzone];
					//zone2markelem[izone][ijzone] += 1;
					//zone2idmark[izone][ijzone].push_back(idz);

					for (int inofa = 0; inofa < nnofa; inofa++) {
						int inoel = vtk2lpofa[vtk][inofa][ifael]; // Getlpofa()
						int inode_g = read.elem2node[ielem_g][inoel] - 1; // elem2node en 1 based 999999999
						int inode_z = nodeglobal2local[inode_g][izone];
						belem2node[izone][idz][inofa + 1] = inode_z;
					}
					indexzone[izone] += 1;
					
				}
			}

		}
	}

	delete[] indexzone;
}
void Connect_c::Zone2Zones() {

}


void	Connect_c::Display1DArray(int V[], int size, string name) {
	cout << "\n" << name << " = [ ";
	for (int i = 0; i <= size - 1; i++) {
		cout << V[i] << ", ";
	}
	cout << "]\n";
};

void	Connect_c::Display2DArray(int* V[], int nLine, int nColomn, string name) {
	cout << "\n" << name << " = [ ";
	for (int i = 0; i <= nLine - 1; i++) {
		cout << "\n";
		for (int j = 0; j <= nColomn - 1; j++) {
			cout << V[i][j] << ", ";
		}
	}
	cout << "]\n";
}
void	Connect_c::Display2DArray(double* V[], int nLine, int nColomn, string name) {
	cout << "\n" << name << " = [ ";
	for (int i = 0; i <= nLine - 1; i++) {
		cout << "\n";
		for (int j = 0; j <= nColomn - 1; j++) {
			cout << "\t" << V[i][j] << ",";
		}
	}
	cout << "]\n";
}

void	Connect_c::Display3DArray(int** V[], int ix, int nLine, int nColomn, string name) {
	cout << "\n" << name << " = [ ";
	for (int i = 0; i <= nLine - 1; i++) {
		cout << "\n";
		for (int j = 0; j <= nColomn - 1; j++) {
			cout << V[ix][i][j] << ", ";
		}
	}
	cout << "]\n";
}