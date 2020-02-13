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
			ipoil = read.elem2node[ielem][inode] + 1; // elem2node en 1 base 999999999999999999999
			esup2[ipoil - 1] = esup2[ipoil - 1] + 1;
		}
	}
	for (int ipoin = 2; ipoin <= nnode_g + 1; ipoin++) {
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
	delete[] psup1;
	delete[] psup2;
	delete[] lpoin;
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
void Connect_c::ComputeGlobalConnectivity(Reader_c& read) {
	Connect_c::InitializeGlobal(read);
	Connect_c::Node2Elements(read);
	Connect_c::Node2Nodes(read);
	Connect_c::Element2Elements(read);
}
// ============================================= ZONE ELEMENTS CONNECTIVITY ================================================
void Connect_c::Node2Elements(int izone) {
	int mesup1 = 0;
	for (int ielem = 0; ielem < zone2ncell[izone]; ielem++) {
		int vtk = elem2vtk[izone][ielem];
		mesup1 += vtklookup[vtk][1];
	}

	// Initializing esup and other parameters :
	zone2esup2[izone] = new int[zone2nnode[izone] + 1]();
	zone2esup1[izone] = new int[mesup1]();
	int ipoil;
	//
	for (int ielem = 0; ielem < zone2ncell[izone]; ielem++) {
		int vtk = elem2vtk[izone][ielem];
		int nnode = vtklookup[vtk][1];
		for (int inode = 0; inode < nnode; inode++) {
			ipoil = elem2node[izone][ielem][inode] + 2; // elem2node en 0 base 
			zone2esup2[izone][ipoil - 1] = zone2esup2[izone][ipoil - 1] + 1;
		}
	}
	for (int ipoin = 2; ipoin <= zone2nnode[izone] + 1; ipoin++) {
		zone2esup2[izone][ipoin - 1] = zone2esup2[izone][ipoin - 1] + zone2esup2[izone][ipoin - 1 - 1];
	}
	int ipoin, istor;
	for (int ielem = 0; ielem < zone2ncell[izone]; ielem++) {
		int vtk = elem2vtk[izone][ielem];
		int nnode = vtklookup[vtk][1];

		for (int inode = 0; inode <= nnode - 1; inode++) {
			ipoin = elem2node[izone][ielem][inode] + 1; // elem2node en 0 base
			istor = zone2esup2[izone][ipoin - 1] + 1;
			zone2esup2[izone][ipoin - 1] = istor;
			zone2esup1[izone][istor - 1] = ielem + 1;
		}
	}
	for (int ipoin = zone2nnode[izone] + 1; ipoin >= 2; --ipoin) {
		zone2esup2[izone][ipoin - 1] = zone2esup2[izone][ipoin - 1 - 1];
	}
	zone2esup2[izone][0] = 0;

	Display1DArray(zone2esup1[izone], mesup1, "esup1");
	Display1DArray(zone2esup2[izone], zone2nnode[izone] + 1, "esup2");
}
void Connect_c::Node2Nodes(int izone) {
	zone2lpoin[izone] = new int[zone2nnode[izone]]();
	zone2psup2[izone] = new int[zone2nnode[izone] + 1]();
	int istor{ 0 };
	for (int ipoin = 1; ipoin <= zone2nnode[izone]; ipoin++) {
		for (int iesup = zone2esup2[izone][ipoin - 1] + 1; iesup <= zone2esup2[izone][ipoin + 1 - 1]; iesup++) {
			int ielem = zone2esup1[izone][iesup - 1] - 1;
			int vtk = elem2vtk[izone][ielem];
			int nnode = vtklookup[vtk][1];
			for (int inode = 1; inode <= nnode; inode++) {
				int jpoin = elem2node[izone][ielem][inode - 1] + 1; // elem2node en 9 base
				if (jpoin != ipoin && zone2lpoin[izone][jpoin - 1] != ipoin) {
					istor = istor + 1;
					zone2lpoin[izone][jpoin - 1] = ipoin;
				}
			}
		}
		zone2psup2[izone][ipoin + 1 - 1] = istor;
	}

	zone2lpoin[izone] = new int[zone2nnode[izone]]();
	istor = 0; int jpoin = 0;
	int mpsup = zone2psup2[izone][zone2nnode[izone]];
	zone2psup1[izone] = new int[mpsup]();
	for (int ipoin = 1; ipoin <= zone2nnode[izone]; ipoin++) {
		for (int iesup = zone2esup2[izone][ipoin - 1] + 1; iesup <= zone2esup2[izone][ipoin + 1 - 1]; iesup++) {
			int ielem = zone2esup1[izone][iesup - 1] - 1;																		// 1 based
			int vtk = elem2vtk[izone][ielem];
			int nnode = vtklookup[vtk][1];
			for (int inode = 1; inode <= nnode; inode++) {
				jpoin = elem2node[izone][ielem][inode - 1] + 1; // elem2node en 9 base
				if (jpoin != ipoin && zone2lpoin[izone][jpoin - 1] != ipoin) {
					istor = istor + 1;
					zone2psup1[izone][istor - 1] = jpoin;
					zone2lpoin[izone][jpoin - 1] = ipoin;
				}
			}
		}
	}
}
void Connect_c::Element2Elements(int izone) {
	elem2elem[izone] = new int* [zone2ncell[izone]];

	vector<int> lhelp;
	for (int ielem = 0; ielem < zone2ncell[izone]; ielem++) {
		int vtk = elem2vtk[izone][ielem];
		int nfael = vtklookup[vtk][0];
		elem2elem[izone][ielem] = new int[nfael]();
		for (int ifael = 0; ifael <= nfael - 1; ifael++) {
			int nnofa = vtk2nnofa[vtk][ifael];
			lhelp.resize(nnofa);															// vecteur a enlever

			for (int inofa = 0; inofa <= nnofa - 1; inofa++) {
				int ilpofa = vtk2lpofa[vtk][inofa][ifael] + 1;
				lhelp[inofa] = elem2node[izone][ielem][ilpofa - 1] + 1;	// elem2node en 0 base 
				zone2lpoin[izone][lhelp[inofa] - 1] = 1;
			}
			int ipoin = lhelp[0] - 1;
			for (int istor = zone2esup2[izone][ipoin]; istor <= zone2esup2[izone][ipoin + 1] - 1; istor++) {
				int jelem = zone2esup1[izone][istor];
				int jvtk = elem2vtk[izone][jelem - 1];
				int nfael = vtklookup[jvtk][0];
				if (jelem != ielem + 1) {
					for (int jfael = 0; jfael <= nfael - 1; jfael++) {
						int nnofj = vtk2nnofa[jvtk][jfael];
						if (nnofj == nnofa) {
							int icoun = 0;
							for (int jnofa = 0; jnofa <= nnofa - 1; jnofa++) {
								int jpoin = elem2node[izone][jelem - 1][vtk2lpofa[jvtk][jnofa][jfael] + 1 - 1] + 1; // elem2node en 0 base 
								icoun = icoun + zone2lpoin[izone][jpoin - 1];
							}
							if (icoun == nnofa) {
								elem2elem[izone][ielem][ifael] = jelem - 1;

							}
						}
					}
				}
			}
			int vtk2 = elem2vtk[izone][ielem];
			nnofa = vtk2nnofa[vtk2][ifael];
			for (int inova = 0; inova <= nnofa - 1; inova++) {
				zone2lpoin[izone][lhelp[inova] - 1] = 0;
			}
		}
	}

}

// ================================================= FACES CONNECTIVITY ====================================================
void Connect_c::Findnface(int izone) {
	int nfaeltot = 0;
	for (int ielem = 0; ielem < zone2nelem[izone]; ielem++) {
		int vtk = elem2vtk[izone][ielem];
		int nfael = vtklookup[vtk][0];
		nfaeltot += nfael;
	}
	int ngcell = zone2nbound[izone] + zone2nbelem[izone];
	zone2nface[izone] = 0.5 * (nfaeltot - ngcell) + ngcell;
	//cout << "nface = " << zone2nface[izone] << endl;
}
void Connect_c::Face2ElementsNodes(int izone) {
	
	face2elem[izone] = new int* [zone2nface[izone]];
	face2node[izone] = new int* [zone2nface[izone]];
	face2fael[izone] = new int* [zone2nface[izone]];


	

	//for (int iface = 0; iface < zone2nface[izone]; iface++) {
	//	//face2elem[izone][iface] = new int[2];
	//	//face2node[izone][iface] = new int[ndime]();  ///// ne marchera pas tout suite 3D
	//	//face2fael[izone][iface] = new int[ndime];	///// ne marchera pas tout suite 3D
	//}

	// Creation d'une matrice d'aide
	int** elem2elemHelp;
	elem2elemHelp = new int* [zone2ncell[izone]];
	for (int ielem = 0; ielem < zone2ncell[izone]; ielem++) {
		int vtk = elem2vtk[izone][ielem];
		int nfael = vtklookup[vtk][0];
		elem2elemHelp[ielem] = new int[nfael];
		for (int ifael = 0; ifael < nfael; ifael++) {
			elem2elemHelp[ielem][ifael] = elem2elem[izone][ielem][ifael];
		}
	}
	
	//Display2DArray(e2elem, zone2ncell[izone], 4, "esuel");
	//Display2DArray(elem2elemHelp, zone2ncell[izone], 4, "esuel");
    int iface = 0; int LeftCheck = 1;
	for (int ielem = 0; ielem < zone2nelem[izone]; ielem++) {
		int vtk = elem2vtk[izone][ielem];
		int nfael = vtklookup[vtk][0];

		for (int ifael = 0; ifael <= nfael - 1; ifael++) {
			// jelem = elem2elemHelp[ielem][ifael] - 1;   // 1 based
			int jelem = elem2elemHelp[ielem][ifael];   // 0 based

			if (jelem != -2) {

				int nnofa = vtk2nnofa[vtk][ifael];
				face2node[izone][iface] = new int[nnofa];
				face2elem[izone][iface] = new int[2];
				face2fael[izone][iface] = new int[2];

				if (ielem < jelem) {         // element i est a gauche de j  // 1 based
					LeftCheck = 1;
					
					for (int inofa = 0; inofa < nnofa; inofa++) {	
						int inoel = vtk2lpofa[vtk][inofa][ifael];
						face2node[izone][iface][inofa] = elem2node[izone][ielem][inoel];
					}

					face2elem[izone][iface][0] = ielem;
					face2elem[izone][iface][1] = jelem;

					face2fael[izone][iface][0] = ifael;

					iface++;
				}
				else if (ielem > jelem) {      // element i est a droite de j    // 1 based
					LeftCheck = 0;

					for (int inofa = 0; inofa < nnofa; inofa++) {
						int inoel = vtk2lpofa[vtk][inofa][ifael];
						face2node[izone][iface][inofa] = elem2node[izone][ielem][inoel];
					}
					face2elem[izone][iface][0] = jelem;
					face2elem[izone][iface][1] = ielem;

					face2fael[izone][iface][1] = ifael;
					iface++;
				}

				int jvtk = elem2vtk[izone][jelem];
				int njfael = vtklookup[jvtk][0];

				for (int jfael = 0; jfael < njfael; jfael++) {
					//int kelem = elem2elemHelp[jelem][jfael] - 1;   // 1 based
					int kelem = elem2elemHelp[jelem][jfael];   // 0 based
					if (kelem == ielem) {
						elem2elemHelp[jelem][jfael] = -2;  // 999999999999 -1 

						if (LeftCheck == 1) { // element i est a gauche de k 
							face2fael[izone][iface - 1][1] = jfael;
						}
						else {
							face2fael[izone][iface - 1][0] = jfael;
						}

						break;
					}
				}
			}

		}
	}
	delete[] elem2elemHelp;

}
void Connect_c::Element2Faces(int izone) 
{
	elem2face[izone] = new int* [zone2ncell[izone]];
	for (int ielem = 0; ielem < zone2ncell[izone]; ielem++) {
		int vtk = elem2vtk[izone][ielem];
		int nfael = vtklookup[vtk][0];
		elem2face[izone][ielem] = new int[nfael];
	}
	
	for (int iface = 0; iface < zone2nface[izone]; iface++) 
	{
		int ielemL = face2elem[izone][iface][0];
		int ielemR = face2elem[izone][iface][1];
		int ifaelL = face2fael[izone][iface][0];
		int ifaelR = face2fael[izone][iface][1];
		elem2face[izone][ielemL][ifaelL] = iface;
		elem2face[izone][ielemR][ifaelR] = iface;
	}
}
void Connect_c::ComputeLocalConnectivity() 
{
	zone2esup1 = new int* [nzone];
	zone2esup2 = new int* [nzone];
	zone2psup1 = new int* [nzone];
	zone2psup2 = new int* [nzone];
	zone2lpoin = new int* [nzone];
	elem2elem = new int** [nzone];
	zone2nface = new int [nzone];
	face2elem = new int** [nzone];
	face2node = new int** [nzone];
	face2fael = new int** [nzone];
	elem2face = new int** [nzone];

	for (int izone = 0; izone < nzone; izone++) {
		Connect_c::Node2Elements(izone);
		Connect_c::Node2Nodes(izone);
		Connect_c::Element2Elements(izone);
		Connect_c::Findnface(izone);
		Connect_c::Face2ElementsNodes(izone);
		Connect_c::Element2Faces(izone);
	}

	//delete[] zone2esup1;
	//delete[] zone2esup2;
	delete[] zone2psup1;
	delete[] zone2psup2;
	delete[] zone2lpoin;

	/*Display3DArray(face2node, 0, zone2nface[0], 2, "face2node");
	Display3DArray(face2elem, 0, zone2nface[0], 2, "face2elem");
	Display3DArray(face2fael, 0, zone2nface[0], 2, "face2fael");*/
}


// ================================================= ZONES CONNECTIVITY ====================================================
void Connect_c::Findnzone() {
	nzone = 1;
	for (int ielem = 0; ielem < nelem_g; ielem++) {
		int izone = elem2zone[ielem];
		if (izone + 1 > nzone) {
			nzone = izone + 1;
		}
	}
}
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
void Connect_c::InitializeElem2Node(Reader_c& read) {
	zone2boundIndex = new int* [nzone];
	zone2nbelem = new int[nzone];
	zone2nbound = new int[nzone];
	zone2ncell = new int[nzone];
	
	
	for (int izone = 0; izone < nzone; izone++) {
		zone2boundIndex[izone] = new int[read.nbc + 1];
		zone2boundIndex[izone][0] = zone2nelem[izone];
		zone2nbelem[izone] = 0;
		zone2nbound[izone] = 0;

	}

	// Calcul de Taille de belem2node --> zone2nbelem:
	for (int ielem_g = 0; ielem_g < nelem_g; ielem_g++) {
		int izone = elemglobal2local[ielem_g][1];
		int vtk = read.elem2vtk[ielem_g];
		int nfael = vtklookup[vtk][0];

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
	
	
	// Calcul de Taille de boundary --> zone2nbound:
	for (int ibc = 1; ibc < read.nbc + 1; ibc++) {
		int ielem1 = read.BoundIndex[ibc - 1];
		int ielem2 = read.BoundIndex[ibc];
		for (int ighost = ielem1; ighost < ielem2; ighost++) {
			int jelem = elem2elem_g[ighost][0];
			int izone = elem2zone[jelem];			
			zone2nbound[izone] += 1;
		}
		for (int izone = 0; izone < nzone; izone++) {
			zone2boundIndex[izone][ibc] = zone2nelem[izone] + zone2nbound[izone];
		}	
	}

	for (int izone = 0; izone < nzone; izone++) {
		cout << "zone2nbound = " << zone2nbound[izone] << endl;
	}
	for (int izone = 0; izone < nzone; izone++) {
		cout << "zone2boundIndex = " << zone2boundIndex[izone][0] << " " << zone2boundIndex[izone][1] << " " <<  zone2boundIndex[izone][2] <<  endl;
	}

	// ============================= Initialization ============================================

	// Initialisation de elem2node, elem2vtk
	elem2node = new int** [nzone]; 
	elem2vtk = new int* [nzone];
	for (int izone = 0; izone < nzone; izone++) {
		zone2ncell[izone] = zone2nelem[izone] + zone2nbound[izone] + zone2nbelem[izone];
		elem2node[izone] = new int* [zone2ncell[izone]];
		elem2vtk[izone] = new int[zone2ncell[izone]];
		for (int ielem = 0; ielem < zone2nelem[izone]; ielem++) {
			int vtk = read.elem2vtk[ielem];
			int nnoel = vtklookup[vtk][1];
			elem2node[izone][ielem] = new int[nnoel];
			elem2vtk[izone][ielem] = vtk;
		}
	}


}
void Connect_c::Zone2Bound(Reader_c& read) {
	// VOir avec Paul, BoundIndex est en 1 based mais pk
	int* indexzone; indexzone = new int[nzone];
	for (int izone = 0; izone < nzone; izone++) {
		indexzone[izone] = zone2nelem[izone];
	}

	for (int ibc = 1; ibc < read.nbc+1; ibc++) {
		int ielem1 = read.BoundIndex[ibc-1];
		int ielem2 = read.BoundIndex[ibc];
		for (int ighost = ielem1; ighost < ielem2; ighost++) {
			int vtk = read.elem2vtk[ighost];
			int nnoel = vtklookup[vtk][1]; 
			int jelem = elem2elem_g[ighost][0];
			int izone = elem2zone[jelem];
			int idz = indexzone[izone];
			elem2node[izone][idz] = new int[nnoel];
			for (int inoel = 0; inoel < nnoel; inoel++) {
				int inode_g = read.elem2node[ighost][inoel] -1; // -1 elem2node est en 1 based 99999999999999
				int inode_z = nodeglobal2local[inode_g][izone];
				elem2node[izone][idz][inoel] = inode_z;
			}
			elem2vtk[izone][idz] = vtk;
			indexzone[izone] += 1;
		}
		cout << "ielem1 = " << ielem1 << endl;
		cout << "ielem2 = " << ielem2 << endl;
	}
	cout << read.BoundIndex[0] << " " << read.BoundIndex[1] << " " << read.BoundIndex[2] << " " << read.BoundIndex[3] << endl;
	delete[] indexzone;

}
void Connect_c::Element2Nodes(Reader_c& read) 
{	
	Connect_c::InitializeElem2Node(read);
	Connect_c::Zone2Bound(read);

	// Initialization:
	int** izone2jzone;
	int* indexzone;
	indexzone = new int[nzone];
	belem2node = new int**[nzone];
	izone2jzone = new int* [nzone];

	zone2markelem = new int*[nzone];
	//zone2idmark.resize(nzone);

	// Creation de izone2jzone
	for (int izone = 0; izone < nzone; izone++) {
		zone2markelem[izone] = new int[nzone - 1];;
		//zone2idmark[izone].resize(nzone - 1);
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
					int ighost = idz + zone2nelem[izone] + zone2nbound[izone];
					int nnofa = vtk2nnofa[vtk][ifael];
					int fvtk = vtk2facevtk[vtk][ifael];
					belem2node[izone][idz] = new int[nnofa + 2];
					elem2node[izone][ighost] = new int[nnofa];
					belem2node[izone][idz][0] = fvtk; // vtk de la face
					belem2node[izone][idz][nnofa + 1] = jelem_z;
					elem2vtk[izone][ighost] = fvtk;
					int ijzone = izone2jzone[izone][jzone];
					zone2markelem[izone][ijzone] += 1;
					//zone2idmark[izone][ijzone].push_back(idz);

					for (int inofa = 0; inofa < nnofa; inofa++) {
						int inoel = vtk2lpofa[vtk][inofa][ifael]; // Getlpofa()
						int inode_g = read.elem2node[ielem_g][inoel] - 1; // elem2node en 1 based 999999999
						int inode_z = nodeglobal2local[inode_g][izone];
						belem2node[izone][idz][inofa + 1] = inode_z;
						elem2node[izone][ighost][inofa] = inode_z;
					}
					indexzone[izone] += 1;	
				}
			}
		}
	}
	delete[] indexzone;
}
void Connect_c::Zone2Zones() 
{
	// Initialization:
	zone2zone = new int* [nzone];
	zone2ijzone = new int* [nzone];
	for (int izone = 0; izone < nzone; izone++) {
		zone2zone[izone] = new int[nzone];
		zone2ijzone[izone] = new int[nzone - 1];
		int ijzone = 0;
		for (int jzone = 0; jzone < nzone; jzone++) {
			zone2zone[izone][jzone] = -1;
			if (izone != jzone) {
				zone2ijzone[izone][ijzone] = jzone;
				ijzone += 1;
			}
		}
	}
	
	for (int izone = 0; izone < nzone; izone++) {
		int istep = 0;
		for (int ijzone = 0; ijzone < nzone - 1; ijzone++) {
			int jzone = zone2ijzone[izone][ijzone];
			int nmark = zone2markelem[izone][ijzone];
			if (nmark > 0) {
				zone2zone[izone][istep] = jzone;
				istep += 1;
			}
		}
	}

	/*Display2DArray(zone2markelem, nzone, nzone - 1, "zone2markelem");
	Display2DArray(zone2ijzone, nzone, nzone - 1, "zone2ijzone");
	Display2DArray(zone2zone, nzone, nzone, "zone2zone");*/
}


void Connect_c::ComputeZoneConnectivity(Reader_c& read) {
	Connect_c::Findnzone();
	Connect_c::InitializeLocal();
	Connect_c::Zone2nnode();
	Connect_c::Zone2nelem();
	Connect_c::Zone2Nodes();
	Connect_c::Zone2Elements();
	Connect_c::NodeGlobal2Local();
	Connect_c::ElementGlobal2Local();
	Connect_c::Zone2Coord(read);
	Connect_c::Element2Nodes(read);
	Connect_c::Zone2Zones();
}


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