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
void Connect_c::InitializeGlobal(Reader_c& read) {
	//int nelem, nnode, nzelem, nbound, ncell, nface;
	nelem = read.nelem;
	nnode = read.npoint;
	nzelem = read.nzelem;
	nzone = read.nzone;
	ndime = read.ndime;
	ncell = read.ncell;
	nbound = read.nhalo;
	vtk2nnofa = new int* [16];
	vtk2lpofa = new int** [16];
	vtk2facevtk = new int* [16];

	// nnofa : number of node per face of shape
	// vtk2lpofa[vtk][ifael][inoel] = 
	// vtk2fael2noel 
	

	// =================== Bar (vtk = 3) ===================
	//							0 - 1
	vtk2nnofa[3] = new int[1]{ 2 };
	vtk2lpofa[3] = new int* [1];
	vtk2lpofa[3][0] = new int[2]{ 0,1 }; // face 0
	


	// ================= Triangle (vtk = 5) =================
	//						     2
	//						    / \	
	//					       0 _ 1
	if (ndime == 2) {
		vtk2nnofa[5] = new int[3]{ 2,2,2 };
		vtk2facevtk[5] = new int[3]{ 3,3,3 };

		vtk2lpofa[5] = new int* [3];
		vtk2lpofa[5][0] = new int[2]{ 0,1 }; // face 0
		vtk2lpofa[5][1] = new int[2]{ 1,2 }; // face 1
		vtk2lpofa[5][2] = new int[2]{ 2,0 }; // face 2
	}
	else { // 3D
		vtk2nnofa[5] = new int[1]{ 3 };

		vtk2lpofa[5] = new int* [1];
		vtk2lpofa[5][0] = new int[3]{ 0,1,2 }; // face 0
	}

	// =================== Quad (vtk = 9) ===================
	//						    3 - 2
	//						    |   | 
	//						    0 - 1
	if (ndime == 2) {
		vtk2nnofa[9] = new int[4]{ 2,2,2,2 };
		vtk2facevtk[9] = new int[4]{ 3,3,3,3 };

		vtk2lpofa[9] = new int* [4];
		vtk2lpofa[9][0] = new int[2]{ 0, 1 }; // face 0
		vtk2lpofa[9][1] = new int[2]{ 1, 2 }; // face 1
		vtk2lpofa[9][2] = new int[2]{ 2, 3 }; // face 2
		vtk2lpofa[9][3] = new int[2]{ 3, 0 }; // face 3
	}
	else { // 3D
		vtk2nnofa[9] = new int[1]{ 4 };  

		vtk2lpofa[9] = new int* [1];
		vtk2lpofa[9][0] = new int[4]{ 0, 1, 2, 3 }; // face 0
		
	}
	// =============== Tetrahedral (vtk = 10) ================
	//						  3 _ 2
	//						 / \ |
	//						0 __ 1

	vtk2nnofa[10] = new int[4]{ 3,3,3,3 };
	vtk2facevtk[10] = new int[4]{ 5,5,5,5 };

	vtk2lpofa[10] = new int* [4];
	vtk2lpofa[10][0] = new int[3]{ 0,1,3 }; // face 0
	vtk2lpofa[10][1] = new int[3]{ 1,2,3 }; // face 1
	vtk2lpofa[10][2] = new int[3]{ 0,3,2 }; // face 2
	vtk2lpofa[10][3] = new int[3]{ 0,2,1 }; // face 3

	// ================ Hexahedron (vtk = 12) ================
	//					      7 ___ 6
	//					     /|    /|
	//						4_3__ 5 2
	//						| /   | /
	//						0 ___ 1
	vtk2nnofa[12] = new int[6]{ 4,4,4,4,4,4 };
	vtk2facevtk[12] = new int[6]{ 9,9,9,9,9,9 };

	vtk2lpofa[12] = new int* [6];
	vtk2lpofa[12][0] = new int[4]{ 0,3,2,1}; // face 0
	vtk2lpofa[12][1] = new int[4]{ 4,5,6,7}; // face 1
	vtk2lpofa[12][2] = new int[4]{ 1,2,6,5}; // face 2
	vtk2lpofa[12][3] = new int[4]{ 2,3,7,6}; // face 3
	vtk2lpofa[12][4] = new int[4]{ 0,4,7,3}; // face 4
	vtk2lpofa[12][5] = new int[4]{ 0,1,5,4}; // face 5

	// ================ Pyramid (vtk = 14) ================		     
	//					   3 __ 4__ 2
	//				    	|  /\  |
	//	                    | /  \ |
	//	                    |/____\|
	//                     0        1

	vtk2nnofa[14] = new int[5]{ 3,3,3,3,4 };
	vtk2facevtk[14] = new int[5]{ 5,5,5,5,9 };

	vtk2lpofa[14] = new int* [5];
	vtk2lpofa[14][0] = new int[3]{ 0,1,4 }; // face 0
	vtk2lpofa[14][1] = new int[3]{ 1,2,4 }; // face 1
	vtk2lpofa[14][2] = new int[3]{ 2,3,4 }; // face 2
	vtk2lpofa[14][3] = new int[3]{ 3,0,4 }; // face 3
	vtk2lpofa[14][4] = new int[4]{ 0,3,2,1 }; // face 4

}


// ============================================= ZONE ELEMENTS CONNECTIVITY ================================================
// === ELEMENT === 
void Connect_c::Node2Elements(Reader_c& read) {
	int mesup1 = 0;
	for (int ielem = 0; ielem < ncell; ielem++) {
		int vtk = read.elem2vtk[ielem];
		mesup1 += vtklookup[ndime-2][vtk][1];
	}
	// Initializing esup and other parameters :
	esup2 = new int[nnode + 1]();
	esup1 = new int[mesup1]();
	int ipoil;
	
	for (int ielem = 0; ielem < ncell; ielem++) {
		int vtk = read.elem2vtk[ielem];
		int nnode = vtklookup[ndime-2][vtk][1];
		for (int inode = 0; inode < nnode; inode++) {
			ipoil = read.elem2node[ielem][inode] + 2; // elem2node en 0 base 
			esup2[ipoil - 1] = esup2[ipoil - 1] + 1;
		}
	}
	for (int ipoin = 2; ipoin <= nnode + 1; ipoin++) {
		esup2[ipoin - 1] = esup2[ipoin - 1] + esup2[ipoin - 1 - 1];
	}
	int ipoin, istor;
	for (int ielem = 0; ielem < ncell; ielem++) {
		int vtk = read.elem2vtk[ielem];
		int nnode = vtklookup[ndime-2][vtk][1];

		for (int inode = 0; inode <= nnode - 1; inode++) {
			ipoin = read.elem2node[ielem][inode] + 1; // elem2node en 0 base
			istor = esup2[ipoin - 1] + 1;
			esup2[ipoin - 1] = istor;
			esup1[istor - 1] = ielem + 1;
		}
	}
	for (int ipoin = nnode + 1; ipoin >= 2; --ipoin) {
		esup2[ipoin - 1] = esup2[ipoin - 1 - 1];
	}
	esup2[0] = 0;
	

}
void Connect_c::Node2Nodes(Reader_c& read) {
	lpoin = new int[nnode]();
	psup2 = new int[nnode + 1]();
	int istor{ 0 };
	for (int ipoin = 1; ipoin <= nnode; ipoin++) {
		for (int iesup = esup2[ipoin - 1] + 1; iesup <= esup2[ipoin + 1 - 1]; iesup++) {
			int ielem = esup1[iesup - 1] - 1;
			int vtk = read.elem2vtk[ielem];
			int nnode = vtklookup[ndime-2][vtk][1];
			for (int inode = 1; inode <= nnode; inode++) {
				int jpoin = read.elem2node[ielem][inode - 1] + 1;
				if (jpoin != ipoin && lpoin[jpoin - 1] != ipoin) {
					istor = istor + 1;
					lpoin[jpoin - 1] = ipoin;
				}
			}
		}
		psup2[ipoin + 1 - 1] = istor;
	}

	lpoin = new int[nnode]();
	istor = 0; int jpoin = 0;
	int mpsup = psup2[nnode];
	psup1 = new int[mpsup]();
	for (int ipoin = 1; ipoin <= nnode; ipoin++) {
		for (int iesup = esup2[ipoin - 1] + 1; iesup <= esup2[ipoin + 1 - 1]; iesup++) {
			int ielem = esup1[iesup - 1] - 1;																		// 1 based
			int vtk = read.elem2vtk[ielem];
			int nnode = vtklookup[ndime-2][vtk][1];
			for (int inode = 1; inode <= nnode; inode++) {
				jpoin = read.elem2node[ielem][inode - 1] + 1; 
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
	elem2elem = new int* [ncell];

	vector<int> lhelp;
	for (int ielem = 0; ielem < ncell; ielem++) {
		int vtk = read.elem2vtk[ielem];
		int nfael = vtklookup[ndime-2][vtk][0];
		elem2elem[ielem] = new int[nfael]();
		for (int ifael = 0; ifael <= nfael - 1; ifael++) {
			int nnofa = vtk2nnofa[vtk][ifael];
			lhelp.resize(nnofa);															// vecteur a enlever

			for (int inofa = 0; inofa <= nnofa - 1; inofa++) {
				int ilpofa = vtk2lpofa[vtk][ifael][inofa] + 1;
				lhelp[inofa] = read.elem2node[ielem][ilpofa - 1] + 1;	// elem2node en 0 base 
				lpoin[lhelp[inofa] - 1] = 1;
			}
			int ipoin = lhelp[0] - 1;
			for (int istor = esup2[ipoin]; istor <= esup2[ipoin + 1] - 1; istor++) {
				int jelem = esup1[istor];
				int jvtk = read.elem2vtk[jelem - 1];
				int nfael = vtklookup[ndime-2][jvtk][0];
				if (jelem != ielem + 1) {
					for (int jfael = 0; jfael <= nfael - 1; jfael++) {
						int nnofj = vtk2nnofa[jvtk][jfael];
						if (nnofj == nnofa) {
							int icoun = 0;
							for (int jnofa = 0; jnofa <= nnofa - 1; jnofa++) {
								int jpoin = read.elem2node[jelem - 1][vtk2lpofa[jvtk][jfael][jnofa] + 1 - 1] + 1; // elem2node en 0 base 
								icoun = icoun + lpoin[jpoin - 1];
							}
							if (icoun == nnofa) {
								elem2elem[ielem][ifael] = jelem - 1;

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

}

// ==== FACE ====
void Connect_c::Findnface(Reader_c& read) {
	int nfaeltot = 0;
	for (int ielem = 0; ielem < nelem; ielem++) {
		int vtk = read.elem2vtk[ielem];
		int nfael = vtklookup[ndime-2][vtk][0];
		nfaeltot += nfael;
	}
	int ngcell = nbound + nzelem;
	nface = 0.5 * (nfaeltot - ngcell) + ngcell;
	//cout << "nface = " << nface << endl;
}
void Connect_c::Face2ElementsNodes(Reader_c& read) {
	
	face2elem = new int* [nface];
	face2node = new int* [nface];
	face2fael = new int* [nface];
	//face2Nbr_of_node = new int [nface];

	// Creation d'une matrice d'aide identique a elem2elem
	int** elem2elemHelp;
	elem2elemHelp = new int* [ncell];
	for (int ielem = 0; ielem < ncell; ielem++) {
		int vtk = read.elem2vtk[ielem];
		int nfael = vtklookup[ndime-2][vtk][0];
		elem2elemHelp[ielem] = new int[nfael];
		for (int ifael = 0; ifael < nfael; ifael++) {
			elem2elemHelp[ielem][ifael] = elem2elem[ielem][ifael];
		}
	}
    int iface = 0; int LeftCheck = 1;
	for (int ielem = 0; ielem < nelem; ielem++) {
		int vtk = read.elem2vtk[ielem];
		int nfael = vtklookup[ndime-2][vtk][0];

		for (int ifael = 0; ifael <= nfael - 1; ifael++) {
			// jelem = elem2elemHelp[ielem][ifael] - 1;   // 1 based
			int jelem = elem2elemHelp[ielem][ifael];   // 0 based

			if (jelem != -2) {

				int nnofa = vtk2nnofa[vtk][ifael];
				face2node[iface] = new int[nnofa];
				//face2Nbr_of_node[iface] = nnofa;
				face2elem[iface] = new int[2];
				face2fael[iface] = new int[2];

				if (ielem < jelem) {         // element i est a gauche de j  // 1 based
					LeftCheck = 1;
					
					for (int inofa = 0; inofa < nnofa; inofa++) {	
						int inoel = vtk2lpofa[vtk][ifael][inofa];
						face2node[iface][inofa] = read.elem2node[ielem][inoel];
					}
					face2elem[iface][0] = ielem;
					face2elem[iface][1] = jelem;
					face2fael[iface][0] = ifael;
					iface++;
				}
				else if (ielem > jelem) {      // element i est a droite de j    // 1 based
					LeftCheck = 0;

					for (int inofa = 0; inofa < nnofa; inofa++) {
						int inoel = vtk2lpofa[vtk][ifael][inofa];
						face2node[iface][inofa] = read.elem2node[ielem][inoel];
					}
					face2elem[iface][0] = jelem;
					face2elem[iface][1] = ielem;
					face2fael[iface][1] = ifael;
					iface++;
				}
				int jvtk = read.elem2vtk[jelem];
				int njfael = vtklookup[ndime-2][jvtk][0];

				for (int jfael = 0; jfael < njfael; jfael++) {
					int kelem = elem2elemHelp[jelem][jfael];   // 0 based
					if (kelem == ielem) {
						elem2elemHelp[jelem][jfael] = -2; 
						face2fael[iface - 1][LeftCheck] = jfael;
						break;
					}
				}
			}
		}
	}
	delete[] elem2elemHelp;
}
void Connect_c::Element2Faces(Reader_c& read) 
{
	elem2face = new int* [ncell];
	for (int ielem = 0; ielem < ncell; ielem++) {
		int vtk = read.elem2vtk[ielem];
		int nfael = vtklookup[ndime-2][vtk][0];
		elem2face[ielem] = new int[nfael];
	}
	for (int iface = 0; iface < nface; iface++) 
	{	
		int ielemL = face2elem[iface][0];
		int ielemR = face2elem[iface][1];
		int ifaelL = face2fael[iface][0];
		int ifaelR = face2fael[iface][1];
		elem2face[ielemL][ifaelL] = iface;
		elem2face[ielemR][ifaelR] = iface;
	}
}
void Connect_c::ComputeLocalConnectivity(Reader_c& read) 
{
	cout << "Local Connectivity \tSTARTING...";
	Connect_c::InitializeGlobal(read);
	Connect_c::Node2Elements(read); 
	Connect_c::Node2Nodes(read);
	Connect_c::Element2Elements(read); 
	Connect_c::Findnface(read); 
	Connect_c::Face2ElementsNodes(read);
	Connect_c::Element2Faces(read);
	
	delete[] esup1;
	delete[] esup2;
	delete[] psup1;
	delete[] psup2;
	delete[] lpoin;

	/*Display3DArray(face2node, 0, nface[0], 5, "face2node");
	Display3DArray(face2elem, 0, nface[0], 2, "face2elem");
	Display3DArray(face2fael, 0, nface[0], 2, "face2fael");*/
	cout << "...............DONE" << endl;
}


// ================================================= ZONES CONNECTIVITY ====================================================




// ================================================= OTHER FUNCTION MEMBERS ====================================================
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
void	Connect_c::Display3DArray(double** V[], int ix, int nLine, int nColomn, string name) {
	cout << "\n" << name << " = [ ";
	for (int i = 0; i <= nLine - 1; i++) {
		cout << "\n";
		for (int j = 0; j <= nColomn - 1; j++) {
			cout << V[ix][i][j] << ", ";
		}
	}
	cout << "]\n";
}