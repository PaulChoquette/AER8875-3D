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
#include <metis.h>
using namespace std;

// ================================================= ELEMENTS CONNECTIVITY ====================================================
void Connect_c::InitializeGlobal(Reader_c& read) {
	nnode_g = read.npoint;
	nelem_g = read.nelem;
	ndime = read.ndime;
	ncell_g = read.ncell;

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
void Connect_c::Node2Elements(Reader_c& read) {
	mesup1 = 0;
	for (int ielem = 0; ielem < ncell_g; ielem++) {
		int vtk = read.elem2vtk[ielem];
		mesup1 += vtklookup[ndime-2][vtk][1];
	}

	// Initializing esup and other parameters :
	esup2 = new int[nnode_g + 1]();
	esup1 = new int[mesup1]();
	int ipoil;
	//
	for (int ielem = 0; ielem < ncell_g; ielem++) {
		int vtk = read.elem2vtk[ielem];
		int nnode = vtklookup[ndime-2][vtk][1];
		for (int inode = 0; inode < nnode; inode++) {
			ipoil = read.elem2node[ielem][inode] + 1 + 1; // elem2node en 0 based
			esup2[ipoil - 1] = esup2[ipoil - 1] + 1;
		}
	}
	for (int ipoin = 2; ipoin <= nnode_g + 1; ipoin++) {
		esup2[ipoin - 1] = esup2[ipoin - 1] + esup2[ipoin - 1 - 1];
	}
	int ipoin, istor;
	for (int ielem = 0; ielem < ncell_g; ielem++) {
		int vtk = read.elem2vtk[ielem];
		int nnode = vtklookup[ndime-2][vtk][1];

		for (int inode = 0; inode <= nnode - 1; inode++) {
			ipoin = read.elem2node[ielem][inode]+1; // elem2node en o base
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
			int nnode = vtklookup[ndime-2][vtk][1];
			for (int inode = 1; inode <= nnode; inode++) {
				jpoin = read.elem2node[ielem][inode - 1] + 1; // elem2node en 0 base
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
			int nnode = vtklookup[ndime-2][vtk][1];
			for (int inode = 1; inode <= nnode; inode++) {
				jpoin = read.elem2node[ielem][inode - 1] +1; // elem2node en 0 based
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
		int nfael = vtklookup[ndime-2][vtk][0];
		elem2elem_g[ielem] = new int[nfael];
		for (int ifael = 0; ifael <= nfael - 1; ifael++) {
			nnofa = vtk2nnofa[vtk][ifael];
			lhelp.resize(nnofa);															// vecteur a enlever

			for (int inofa = 0; inofa <= nnofa - 1; inofa++) {
			
				ilpofa = vtk2lpofa[vtk][ifael][inofa]+1;
				lhelp[inofa] = read.elem2node[ielem][ilpofa - 1] + 1;	// elem2node en 0 base
				lpoin[lhelp[inofa] - 1] = 1;
			}
			ipoin = lhelp[0] - 1;
			for (int istor = esup2[ipoin]; istor <= esup2[ipoin + 1] - 1; istor++) {
				jelem = esup1[istor];
				int jvtk = read.elem2vtk[jelem-1];
				int nfael = vtklookup[ndime-2][jvtk][0];
				if (jelem != ielem + 1) {
					for (int jfael = 0; jfael <= nfael - 1; jfael++) {
						//nnofj = Getlnofa(jelem - 1, jfael);  // lnofa[jfael];
						nnofj = vtk2nnofa[jvtk][jfael];
						if (nnofj == nnofa) {
							icoun = 0;
							for (int jnofa = 0; jnofa <= nnofa - 1; jnofa++) {
								//jpoin = read.elem2node[jelem - 1][Getlpofa(jelem - 1, jnofa, jfael) + 1 - 1];
								jpoin = read.elem2node[jelem - 1][vtk2lpofa[jvtk][jfael][jnofa] + 1 - 1] + 1;  // ICI
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
}
void Connect_c::METISElement2Nodes(Reader_c& read){
	neind = 0;
	for (int ielem = 0; ielem < nelem_g; ielem++) {
		int vtk = read.elem2vtk[ielem];
		neind += vtklookup[ndime-2][vtk][1];
	}


	eind = new int[neind];
	eptr = new int[nelem_g+1];
	eptr[0] = 0;
	int ipos = 0;
	for (int ielem = 0; ielem < nelem_g; ielem++){
		int vtk = read.elem2vtk[ielem];
		int nnoel = vtklookup[ndime-2][vtk][1];
		for (int inoel = 0; inoel < nnoel; inoel++){
			int inode = read.elem2node[ielem][inoel];
			eind[ipos] = inode + 1;
			ipos += 1;
		}
		eptr[ielem+1] = ipos;
	}
}
void Connect_c::ComputeMETIS(int nzoneMETIS, int ncommon, Reader_c& read) {
	double StartTime = omp_get_wtime();
	cout << "METIS Connectivity \tSTARTING...";
	Connect_c::METISElement2Nodes(read);

	//Connect_c::Display1DArray(eind,neind,"eind");
	//Connect_c::Display1DArray(eptr,nelem_g+1,"eptr");

	elem2zone = new int[nelem_g];

	int* npart;
	int objval;
	//int ncommon = 4; // Input qui depend du nombre de noeuds par face (nombre de noeud en commun)

	npart = new int[nnode_g];

	int success = METIS_PartMeshDual(&nelem_g, &nnode_g, &eptr[0], &eind[0], NULL, NULL,
		&ncommon, &nzoneMETIS, NULL, NULL, &objval,
		&elem2zone[0], &npart[0]);

	//Connect_c::Display1DArray(elem2zone,nelem_g,"elem2zone");

	delete[] eptr;
	delete[] eind;
	delete[] npart;
	double EndTime = omp_get_wtime();
	double WorkTime = EndTime - StartTime;
	cout << "...............DONE\t (took " << WorkTime << " sec)" << endl;
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
	double StartTime = omp_get_wtime();
	cout << "Global Connectivity \tSTARTING...";
	Connect_c::InitializeGlobal(read);
	Connect_c::Node2Elements(read);
	Connect_c::Node2Nodes(read);
	Connect_c::Element2Elements(read);
	delete[] psup1;
	delete[] psup2;
	delete[] lpoin;
	
	double EndTime = omp_get_wtime();
	double WorkTime = EndTime - StartTime;
	cout << "...............DONE\t (took " << WorkTime << " sec)" << endl;

	//Connect_c::Display2DArray(elem2elem_g,ncell_g,4,"elem2elem");
}

// ============================================= ZONE ELEMENTS CONNECTIVITY ================================================
// === ELEMENT === 


// ================================================= ZONES CONNECTIVITY ====================================================
void Connect_c::Findnzone() {
	// A SUPPRIMER
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
			//Verifie si on na pas deja compter le noeud dans la zone:
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
			//Verifie si on na pas deja compter le noeud dans la zone:
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
void Connect_c::Zone2jZone(){
	// Initialization:
	zone2jzone = new int* [nzone];
	zone2ijzone = new int* [nzone];
	for (int izone = 0; izone < nzone; izone++) {
		zone2jzone[izone] = new int[nzone - 1];
		zone2ijzone[izone] = new int[nzone];
		int ijzone = 0;
		for (int jzone = 0; jzone < nzone; jzone++) {
			if (izone != jzone) {
				zone2jzone[izone][ijzone] = jzone;
				zone2ijzone[izone][jzone] = ijzone;
				ijzone += 1;
			}
			else {
				zone2ijzone[izone][jzone] = -1;
			}
		}
	}

	//Display2DArray(zone2jzone, nzone, nzone - 1, "zone2jzone");
	//Display2DArray(zone2ijzone, nzone, nzone, "zone2ijzone");
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
		int nfael = vtklookup[ndime-2][vtk][0];

		for (int ifael = 0; ifael < nfael; ifael++) {
			int jelem_g = elem2elem_g[ielem_g][ifael];

			// Si jelem nest pas un ghostcell:
			if (jelem_g < nelem_g) {
				int jzone = elem2zone[jelem_g];

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
			int nnoel = vtklookup[ndime-2][vtk][1];
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
			int nnoel = vtklookup[ndime-2][vtk][1]; 
			int jelem = elem2elem_g[ighost][0];
			int izone = elem2zone[jelem];
			int idz = indexzone[izone];
			elem2node[izone][idz] = new int[nnoel];
		
			for (int inoel = 0; inoel < nnoel; inoel++) {
				int inode_g = read.elem2node[ighost][inoel]; // 
				int inode_z = nodeglobal2local[inode_g][izone];
				elem2node[izone][idz][inoel] = inode_z;
				
			}
			elem2vtk[izone][idz] = vtk;
			indexzone[izone] += 1;
		}
		// cout << "ielem1 = " << ielem1 << endl;
		// cout << "ielem2 = " << ielem2 << endl;
	}
	//cout << read.BoundIndex[0] << " " << read.BoundIndex[1] << " " << read.BoundIndex[2] << " " << read.BoundIndex[3] << endl;
	delete[] indexzone;

}
void Connect_c::Element2Nodes(Reader_c& read) 
{	
	Connect_c::InitializeElem2Node(read);
	Connect_c::Zone2Bound(read);
	Connect_c::Zone2jZone();
	// Initialization:
	int** izone2jzone;
	int* indexzone;
	indexzone = new int[nzone];
	belem2node = new int**[nzone];
	izone2jzone = new int* [nzone];
	zone2zoneIndex = new int*[nzone];
	zone2markelem = new int*[nzone];
	zelem2jelem = new int*[nzone];
	zone2idmark.resize(nzone);

	// Creation de izone2jzone
	for (int izone = 0; izone < nzone; izone++) {
		zone2markelem[izone] = new int[nzone - 1];
		zone2zoneIndex[izone] = new int[nzone];
		zone2zoneIndex[izone][0] = zone2nelem[izone] + zone2nbound[izone];
		zone2idmark[izone].resize(nzone - 1);
		izone2jzone[izone] = new int[nzone];
		belem2node[izone] = new int*[zone2nbelem[izone]];
		zelem2jelem[izone] = new int[zone2nbelem[izone]];
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
			//zone2zoneIndex[izone][jzone] = 0;
		}

		for (int ijzone = 0; ijzone < nzone - 1; ijzone++) {
			zone2markelem[izone][ijzone] = 0;
		}
	}


	// Debut du calcul:
	int*** elem2ghost;
	elem2ghost = new int**[nelem_g];
	for (int ielem_g = 0; ielem_g < nelem_g; ielem_g++){
		int vtk = read.elem2vtk[ielem_g];
		int nfael = vtklookup[ndime-2][vtk][0];
		elem2ghost[ielem_g] = new int*[nfael+1];
		elem2ghost[ielem_g][0] = new int[0]();

	}
	for (int ielem_g = 0; ielem_g < nelem_g; ielem_g++) {
		int ielem_z = elemglobal2local[ielem_g][0];
		int izone = elemglobal2local[ielem_g][1];
		int vtk = read.elem2vtk[ielem_g];
		int nnoel = vtklookup[ndime-2][vtk][1];
		int nfael = vtklookup[ndime-2][vtk][0];

		for (int inoel = 0; inoel < nnoel; inoel++) {
			int inode_g = read.elem2node[ielem_g][inoel]; 
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
					belem2node[izone][idz] = new int[nnofa + 2 + 1];
					belem2node[izone][idz][0] = fvtk; // vtk de la face
					belem2node[izone][idz][nnofa + 1] = jelem_z;
					belem2node[izone][idz][nnofa + 2] = ielem_z;
					int ijzone = izone2jzone[izone][jzone];
					zone2markelem[izone][ijzone] += 1;
					zone2idmark[izone][ijzone].push_back(idz);

					// int pos = elem2ghost[jelem_g][0][0];
					// if (pos > 0){
					//     for(int ipos=1; ipos <= pos; ipos++){
					//  		int kelem =  elem2ghost[jelem_g][ipos][0];
					//  		int kdz = elem2ghost[jelem_g][ipos][1];
					//  		if (ielem_g == kelem) {

					//  			belem2node[izone][idz][nnofa + 1] = kdz+ zone2nelem[izone] + zone2nbound[izone];
					//  			belem2node[jzone][kdz][nnofa + 1] = idz+ zone2nelem[izone] + zone2nbound[izone];
				 	// 	}

					//  	}
				    // }


					// int jpos = elem2ghost[jelem_g][0][0];
					// elem2ghost[jelem_g][jpos+1] = new int[2];
					// elem2ghost[jelem_g][jpos+1][0] = ielem_g;
					// elem2ghost[jelem_g][jpos+1][1] = idz;
					// elem2ghost[jelem_g][0][0] = jpos+1;

					

					for (int inofa = 0; inofa < nnofa; inofa++) {
						int inoel = vtk2lpofa[vtk][ifael][inofa]; // Getlpofa()
						int inode_g = read.elem2node[ielem_g][inoel]; // 
						int inode_z = nodeglobal2local[inode_g][izone];
						belem2node[izone][idz][inofa + 1] = inode_z;
						//elem2node[izone][ighost][inofa] = inode_z;
					}
					indexzone[izone] += 1;	
				}
			}
		}
	}
	//Display3DArray(belem2node, 0, zone2nbelem[0],6, "belem2node");
	//Display3DArray(belem2node, 1, zone2nbelem[1],6, "belem2node");
	// Construction des Zones Boundaries
	
	for (int izone = 0; izone < nzone; izone++){
		int eid = 0;
		for (int ijzone = 0; ijzone < nzone - 1; ijzone++){
			int njzone = zone2markelem[izone][ijzone];
			zone2zoneIndex[izone][ijzone+1] = zone2zoneIndex[izone][ijzone] + njzone;
			for (int jelem = 0; jelem < njzone; jelem++){	
				int idz = zone2idmark[izone][ijzone][jelem];
				int ielem  = eid + zone2nelem[izone] + zone2nbound[izone];
				//cout << "iZONE" << izone << " ZONE_TAG=" << ijzone << " IDZ=" << idz << endl;
				int vtk = belem2node[izone][idz][0];
				int nnofa = vtk2nnofa[vtk][0];
				int jelem_z = belem2node[izone][idz][nnofa+1];
				int ielem_z = belem2node[izone][idz][nnofa+2];
				elem2vtk[izone][ielem] = vtk;
				elem2node[izone][ielem] = new int[nnofa + 1];
				elem2node[izone][ielem][nnofa] = jelem_z; 

				for (int inofa = 0; inofa < nnofa; inofa++) {
					int inode_z = belem2node[izone][idz][inofa+1];
					elem2node[izone][ielem][inofa] = inode_z;
				}
				eid += 1;

			
				int jzone = zone2jzone[izone][ijzone];
				int jizone = zone2ijzone[jzone][izone];
				int jeid = 0;

				for (int kjzone = 0; kjzone < jizone; kjzone++){
					jeid += zone2markelem[jzone][kjzone];
				}
				for (int jghost = 0; jghost < njzone; jghost++){
					int jdz = zone2idmark[jzone][jizone][jghost];
					int jvtk = belem2node[jzone][jdz][0];
					int jnnofa = vtk2nnofa[jvtk][0];
					int jelem2_z = belem2node[jzone][jdz][jnnofa+1];
					int ielem2_z = belem2node[jzone][jdz][jnnofa+2];	
					if (jelem2_z == ielem_z && ielem2_z == jelem_z){
						int jelem_ghost = jeid + zone2nelem[jzone] + zone2nbound[jzone];
						elem2node[izone][ielem][nnofa] = jelem_ghost;	
						break;
					}	
						jeid += 1;
				}
			}
		}


	
	// for (int ijzone = 0; ijzone < nzone-1; ijzone++){
	// 	ijzone2ghostData[ijzone] = new int*[ndata];
	// 	for (ijkzone = 0; ijkzone < ndata; ijkzone++){
	// 		int njzone = 0;
	// 		ijzone2ghostData[ijzone][ijkzone] = new int[]
	// 	}
	// 	ndata -= 1;
	// }
			
	}
	delete[] indexzone;
	delete[] belem2node;
	delete[] elem2ghost;
}
void Connect_c::Zone2Zones() 
{
	// Initialization:
	zone2zone = new int* [nzone];
	//zone2jzone = new int* [nzone];
	for (int izone = 0; izone < nzone; izone++) {
		zone2zone[izone] = new int[nzone];
		//zone2jzone[izone] = new int[nzone - 1];
		int ijzone = 0;
		for (int jzone = 0; jzone < nzone; jzone++) {
			zone2zone[izone][jzone] = -1;
			// if (izone != jzone) {
			// 	zone2jzone[izone][ijzone] = jzone;
			// 	ijzone += 1;
			// }
		}
	}
	
	for (int izone = 0; izone < nzone; izone++) {
		int istep = 0;
		for (int ijzone = 0; ijzone < nzone - 1; ijzone++) {
			int jzone = zone2jzone[izone][ijzone];
			int nmark = zone2markelem[izone][ijzone];
			if (nmark > 0) {
				zone2zone[izone][istep] = jzone;
				istep += 1;
			}
		}
	}

	

	//Display2DArray(zone2markelem, nzone, nzone - 1, "zone2markelem");
	//Display2DArray(zone2jzone, nzone, nzone - 1, "zone2jzone");
	//Display2DArray(zone2zone, nzone, nzone, "zone2zone");*/
}


void Connect_c::ComputeZoneConnectivity(Reader_c& read) {
	double StartTime = omp_get_wtime();
	cout << "Zone Connectivity \tSTARTING...";
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
	// delete[] zone2node;
	// delete[] zone2elem;
	// delete[] nodeglobal2local;
	// delete[] elemglobal2local;
	double EndTime = omp_get_wtime();
	double WorkTime = EndTime - StartTime;
	cout << "...............DONE\t (took " << WorkTime << " sec)" << endl;	
}

void Connect_c::ComputeElem2Zone() {
	int elem2zoneBlock[16] = { 0,0,1,1,0,0,1,1,2,2,3,3,2,2,3,3 };
	elem2zone = new int[nelem_g];
	if (nelem_g != 16){
		cout << "ERROR elem2zone" << endl; 
	}

	for (int ielem=0; ielem < nelem_g; ielem++){

		elem2zone[ielem] = elem2zoneBlock[ielem];
	}

}

Connect_c::~Connect_c(void){

	
	
	// int* zone2nelem; int* zone2nnode;  int* zone2nbelem; int* zone2nbound; int* zone2ncell; 
	
	
	// int** zone2markelem;  int** zone2zone; int** zone2jzone; int** zone2ijzone;
	 //int** zelem2jelem; int** zone2boundIndex; int** zone2zoneIndex; int** zone2esup1; int** zone2esup2; int** zone2psup1; int** zone2psup2; int** zone2lpoin; 
	// // 3D Array of Integrer :
  // int*** belem2node; 

	
	delete[] esup2;
	delete[] esup1;

    for (int izone =0; izone < nzone; izone++){
		delete[] zone2elem[izone];
		delete[] zone2node[izone];
		delete[] elem2vtk[izone];
		for (int ielem=0; ielem < zone2ncell[izone]; ielem++){
			delete elem2node[izone][ielem];			
		}
		delete elem2node[izone];	
	}

	delete[] zone2node;
	delete[] zone2elem;
	delete[] elem2node;
	delete[] elem2vtk;
	delete[] zone2ncell;

    for (int ielem=0; ielem < ncell_g; ielem++){
		delete[] elem2elem_g[ielem];
	 	if (ielem < nelem_g){
	 		delete[] elemglobal2local[ielem];
	 	}
	}
	delete[] elem2elem_g;
	delete[] elemglobal2local;
 	
	for (int inode=0; inode < nnode_g; inode++){
	 	delete[] nodeglobal2local[inode];
	}
	delete[] nodeglobal2local;



	delete[] vtk2nnofa;
	delete[] vtk2lpofa;
	delete[] vtk2facevtk;


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