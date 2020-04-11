#include <string>
#include <cctype>
#include <iostream>
#include <sstream>
#include <fstream>
#include "main.h"
#include "Reader.h"
#include "/home/cd/TECIO/include/TECXXX.h"
#include <iomanip>
#include <sstream>


void Reader_c::write_tecplot(Reader_c &FileContents, const char* out_filename,int world_id, void* p, void* Rho, void* u, void* v, void* w) {

	//Declare and fill X Y Z data arrays
	double* x = new double[FileContents.npoint];
	double* y = new double[FileContents.npoint];
	double* z = new double[FileContents.npoint];

	for (int i = 0; i < FileContents.npoint; i++) {
		x[i] = FileContents.coord[i][0];
		y[i] = FileContents.coord[i][1];
		if (FileContents.ndime > 2) {
			z[i] = FileContents.coord[i][2];
		}

	}

	//Initialize tecini variables and use tecini142 to initialize binary file
	INTEGER4 FileFormat = 0;
	INTEGER4 FileType = 0;
	INTEGER4 Debug = 1;
	INTEGER4 IsDouble = 1;
	char fname_temp[100];
	strcpy(fname_temp, "./");
	strcat(fname_temp, out_filename);
	string Tem = to_string(world_id);
	strcat(fname_temp, Tem.c_str());
	strcat(fname_temp, ".plt");

	const char* fname = fname_temp;

	tecini142("NACA_0012", "X,Y,Z,P,Rho,U,V,W", fname , ".", &FileFormat, &FileType, &Debug, &IsDouble);

	//Initialize teczne 142 variables and use teczne142 to make header for each zone
	int Varlocation[] = { 1,1,1,0,0,0,0,0 }; // 1 is for node-centered, 0 is for cell-centered
	INTEGER4 ZoneType = 5;
	INTEGER4 NumFaces = 0; //Notused
	INTEGER4 ICellMax = 0; //Notused
	INTEGER4 JCellMax = 0; //Notused
	INTEGER4 KCellMax = 0; //Notused
	const double SolTime = 0.0;
	INTEGER4 StrandID = 0;
	INTEGER4 ParentZone = 0;
	INTEGER4 IsBlock = 1; //always 1
	INTEGER4 NumFaceConnec = 0;
	INTEGER4 FaceNeighMode = 2;
	INTEGER4 ShareConnec = 0;

	teczne142("main", &ZoneType, &FileContents.npoint, &FileContents.nelem, &NumFaces, &ICellMax, &JCellMax, &KCellMax, &SolTime
		, &StrandID, &ParentZone, &IsBlock, &NumFaceConnec, &FaceNeighMode, NULL, NULL, NULL, NULL
		, Varlocation, NULL, &ShareConnec);

	//Variable writing. One call of tecdat142 for each variable including coordinates
	tecdat142(&FileContents.npoint, x, &IsDouble);
	tecdat142(&FileContents.npoint, y, &IsDouble);
	tecdat142(&FileContents.npoint, z, &IsDouble);
	tecdat142(&FileContents.nelem, p, &IsDouble);
	tecdat142(&FileContents.nelem, Rho, &IsDouble);
	tecdat142(&FileContents.nelem, u, &IsDouble);
	tecdat142(&FileContents.nelem, v, &IsDouble);
	tecdat142(&FileContents.nelem, w, &IsDouble);

	//connectivity list vector writing and call tecnode142 to write connectivity
	INTEGER4 CC = FileContents.nelem * 8;
	INTEGER4* NData = new INTEGER4[CC];

	int nrepeat;

	for (int i = 0; i < FileContents.nelem; i++) {
		nrepeat = 8 - vtklookup[FileContents.ndime-2][FileContents.elem2vtk[i]][1];
		if (nrepeat == 4) {
			for (int j = 0; j < nrepeat; j++) {
				NData[i * 8 + j] = FileContents.elem2node[i][0] + 1;
			}

			NData[i * 8 + 4] = FileContents.elem2node[i][1] + 1;

			for (int k = 5; k < 8; k++) {
				NData[i * 8 + k] = FileContents.elem2node[i][k - nrepeat] + 1;
			}
		}
		else {
			for (int j = 0; j < nrepeat; j++) {
				NData[i * 8 + j] = FileContents.elem2node[i][0] + 1;
			}
			for (int k = nrepeat; k < 8; k++) {
				NData[i * 8 + k] = FileContents.elem2node[i][k - nrepeat] + 1;
			}
		}
	}

	tecnode142(&CC, NData);

	tecend142();

}

void Reader_c::write_tecplot_ASCII(string FileName,double*p,double*rho,double*u,double*v,double*w){
    fstream outFile;
    outFile.open(FileName, ios::out);
	outFile.precision(15);
	outFile<<scientific;


    outFile << "VARIABLES=\"X\",\"Y\",\"Z\",\"p\",\"rho\",\"u\",\"v\",\"w\"" << endl;
	//outFile << "VARIABLES=\"X\",\"Y\",\"Z\"" << endl;
    //outFile << "VARIABLES=\"X\",\"Y\",\"P\",\"U\",\"V\"" << endl;
    //Zone Carre pour le tecplot
    outFile << "ZONE T=\"Element0\"" << endl; //Changer le nbr elements
    //outFile << "Nodes=" << solve.nnode_g << ", Elements=" << solve.nelem_g << ", ZONETYPE=FEBRICK" << endl;
    //outFile << "Nodes=" << npoint << ", Elements=" << nelem << ", ZONETYPE=FETETRAHEDRON" << endl;
	outFile << "Nodes=" << npoint << ", Elements=" << nelem << ", ZONETYPE=FEBRICK" << endl;
    outFile << "DATAPACKING=BLOCK" << endl;
    outFile << "VARLOCATION = ([4,5,6,7,8] = CELLCENTERED)" << endl;
    string a;                      //ecrire les coordonnees de laxe x a la suite
    for (int j = 0; j < npoint; j++)
    {
        a = to_string(coord[j][0]);
        outFile << a << endl;
    }
    string b;                      // ecrire les coordonnees de laxe y a la suite
    for (int i = 0; i < npoint; i++)
    {
        b = to_string(coord[i][1]);
        outFile << b << endl;
    }
    string z;                      // ecrire les coordonnees de laxe y a la suite
    for (int i = 0; i < npoint; i++)
    {
        z = to_string(coord[i][2]);
        outFile << z << endl;
    }

		// Conservatives
	string C;
	for (int i = 0; i < nelem; i++)
    {
        outFile <<fixed<<setprecision(15)<< p[i] << endl;
    }
	for (int i = 0; i < nelem; i++)
    {
        outFile <<fixed<<setprecision(15)<< rho[i] << endl;
    }
	for (int i = 0; i < nelem; i++)
    {
        outFile <<fixed<<setprecision(15)<< u[i] << endl;
    }
	for (int i = 0; i < nelem; i++)
    {
        outFile <<fixed<<setprecision(15)<< v[i] << endl;
    }
	for (int i = 0; i < nelem; i++)
    {
        outFile <<fixed<<setprecision(15)<< w[i] << endl;
    }

        // Ecriture des noeuds de chaque elements pour les carres
    for (int ielem = 0; ielem <nelem; ielem++)
    {
        int vtk = elem2vtk[ielem];
        int nnoel = vtklookup[ndime-2][vtk][1];
        for (int icol = 0; icol < nnoel; icol++)
        {
            string icols = to_string(elem2node[ielem][icol]+1);
            outFile << icols << " ";
        }
        outFile << endl;
    }



    outFile.close();
}
