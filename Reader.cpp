#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include "Reader.h"
#include "Solver.h"
#include "main.h"
using namespace std;

void Reader_c::read_file(string filename) {
	cout << "Reading of Global SU2 Starting ..." << endl;
	//Initialize mesh constants
	ndime = 0; nelem = 0; npoint = 0; nhalo = 0;
	ncell = 0; bclinen = 0;
	//Initialize information reading counters;
	elem = 0; point = 0; line = ""; bc = 0; nbc = 0;

	//Open the file and continue if file is succesfully opened
	if (OpenFile(filename))
	{
		//Read file line by line
		while (getline(file, line))
		{
			//Tick line number and store line as character array
			linen++;
			cline = line.c_str();
			if (ndime == 0)
			{
				//Read number of dimensions of the mesh
				ndime = Readcnst(line, "NDIME= ");
				continue;
			}
			if (nelem == 0)
			{
				//Read number of elements in the mesh
				nelem = Readcnst(line, "NELEM= ");

				//Define elem2node matrix row size if nelem has been defined
				//Store line where elem2node definitions start
				if (nelem) {
					elem2node_nh = new int* [nelem];
					elem2vtk_nh = new int[nelem];
					nelemlinen = linen;
				}
				continue;
			}

			//If elem2node definition line has been defined
			if (nelemlinen)
			{
				//Check if line defines the element
				if (linen > nelemlinen&& linen <= nelemlinen + nelem)
				{
					Fill_E2N_VTK(cline);
					continue;
				}
			}


			if (npoint == 0)
			{
				//Read number of points in mesh
				npoint = Readcnst(line,"NPOIN= ");

				//Define coordinate matrix if number of points has been defined
				if (npoint)
				{
					npointlinen = linen;
					coord = new double* [npoint];
					for (int i = 0; i < npoint; i++) {
						coord[i] = new double[ndime];
					}
				}

				continue;
			}

			//If number of point definition line has been found
			if (npointlinen)
			{
				if (linen > npointlinen&& linen <= npointlinen + npoint) {
					Fill_coord(cline);
					continue;
				}
			}

			if (nbc == 0) {
				nbc = Readcnst(line, "NMARK= ");
				if (nbc != 0) {
					bound2tag = new string[nbc];
					bclinen = linen;
					BoundIndex = new int[nbc + 1];
					BoundIndex[0] = nelem;
					bc_elem2node = new int** [nbc];
					bc_elem2vtk = new int* [nbc];
					bc_nelemv = new int[nbc];
					bcnl = 2 * nbc; //2 lines per marker for marker name and elemn
					bc = 0;
					step = 0;
					continue;
				}
			}

			if (bclinen != 0) {
				if (linen > bclinen && linen <= bclinen + bcnl) {
					
					//3 steps : step 0 is reading marker tag, step 1 is reading marker elemn and final step is filling bc_elem2nnode
					if (step == 0) {
						bound2tag[bc] = line.substr(12,8);
						step++;
					}
					else if (step == 1) {
						bc_nelem = Readcnst(line, "MARKER_ELEMS= ");
						bc_nelemv[bc] = bc_nelem;
						BoundIndex[bc + 1] =BoundIndex[bc]+ bc_nelem;
						cout << "bc_nelem " << bc_nelem << endl;
						bcnl += bc_nelem; //add nelem since each elem has 1 line
						bc_elem2vtk[bc] = new int [bc_nelem];
						bc_elem2node[bc] = new int* [bc_nelem];

						bc_e2n_counter = 0; //counter for marker element number used in FillMarker
						step++;
					}
					else if (step == 2) {
						Fill_BC_E2N_VTK(cline, bc);
						if (bc_e2n_counter == bc_nelem) {
							bc++;
							step = 0;
						}
					}
					if (bc == nbc) {
						nhalo = bcnl - 2 * nbc;
						ncell = nelem + nhalo;
						cout << "TEST Read 6.1 " << endl;
						elem2node = new int* [ncell];
						cout << "TEST Read 6.2 " << endl;
						elem2vtk = new int [ncell];
						cout << "TEST Read 6.3 " << endl;
						for (int ielem = 0; ielem < nelem; ielem++)
						{
							elem2vtk[ielem] = elem2vtk_nh[ielem];
							elem2node[ielem] = new int[vtklookup[ndime-2][elem2vtk[ielem]][1]];
							for (int inode = 0; inode < vtklookup[ndime-2][elem2vtk[ielem]][1]; inode++)
							{
								elem2node[ielem][inode] = elem2node_nh[ielem][inode];
							}
						}

						int imelem = 0;
						for (int ibc = 0; ibc < nbc; ibc++)
						{
							for (int ibcen = 0; ibcen < bc_nelemv[ibc]; ibcen++)
							{
								elem2vtk[nelem + imelem] = bc_elem2vtk[ibc][ibcen];
								elem2node[nelem + imelem] = new int[vtklookup[ndime-2][elem2vtk[nelem + imelem]][1]];

								for (int j = 0; j < vtklookup[ndime-2][elem2vtk[nelem + imelem]][1]; j++)
								{
									elem2node[nelem + imelem][j] = bc_elem2node[ibc][ibcen][j];
								}
								imelem++;
							}
						}
					}
				}
			}

		}
		for (int i = 0; i < nelem;  i++) {
			delete[] elem2node_nh[i];
		}
		delete[] elem2node_nh;

		for (int i = 0; i < nbc; i++) {
			for ( int j = 0; j < bc_nelemv[i]; j++) {
				delete[] bc_elem2node[i][j];
			}
			delete[] bc_elem2node[i];
		}
		delete[] bc_elem2node;

		delete[] bc_elem2vtk;
		delete[] elem2vtk_nh;

		file.close();
	}
	else {
		//ERROR 1 : File could not be opened
	}
	cout << "... Reading ENDING" << endl;
}

bool Reader_c::OpenFile(string filename)
{
	file.open(filename);

	if (file.is_open()) {
		return true;
	}
	else {
		return false;
	}
}

int Reader_c::Readcnst(const string& line, const string& tofind)
{
	size_t cnstpos = line.find(tofind);
	int cnst = 0;
	if (cnstpos != string::npos) {
		cnstpos = line.find_last_of(" ") + 1;

		string cnstch;

		cnstch = line.substr(cnstpos, line.length() - cnstpos);
		cnst = stoul(cnstch);
	}
	else {
		return 0;
	}
	return cnst;
}

void Reader_c::Fill_E2N_VTK(const char* cline) {
	//This function reads a character line and extracts VTK index as well as allocating and storing elem2node rows
	char* end;
	int j = 0;
	for (int c = strtoul(cline, &end, 10); cline != end; c = strtoul(cline, &end, 10))
	{
		cline = end;

		//Where j is the integer counter. 0 is the vtk index and the rest is the
		if (j == 0) {
			elem2vtk_nh[elem] = c;
			elem2node_nh[elem] = new int [vtklookup[ndime-2][c][1]];
		}

		else if(j > 0 && j <= vtklookup[ndime-2][elem2vtk_nh[elem]][1])
		{
			elem2node_nh[elem][j - 1] = c;
		}

		else {
			break;
		}

		j++;

	}

	elem++;
}


void Reader_c::Fill_BC_E2N_VTK(const char* cline, int bc) {
	//This function reads a character line and extracts VTK index as well as allocating and storing elem2node rows
	char* end;
	int j = 0;
	for (int c = strtoul(cline, &end, 10); cline != end; c = strtoul(cline, &end, 10))
	{
		cline = end;

		//Where j is the integer counter. 0 is the vtk index and the rest is the
		if (j == 0) {
			bc_elem2vtk[bc][bc_e2n_counter] = c;
			bc_elem2node[bc][bc_e2n_counter] = new int[vtklookup[ndime-2][c][1]];
		}

		else if (j > 0 && j <= vtklookup[ndime-2][bc_elem2vtk[bc][bc_e2n_counter]][1])
		{
			bc_elem2node[bc][bc_e2n_counter][j - 1] = c;
		}

		else {
			break;
		}

		j++;
	}
	bc_e2n_counter++;
}

double** Reader_c::Fill_coord(const char* cline) {
	char* end;
	int j = 0;
	for (long double c = strtod(cline, &end); cline != end; c = strtod(cline, &end))
	{
		if (j < ndime) {
			cline = end;
			coord[point][j] = c;
		}
		else {
			break;
		}
		j++;
	}
	point++;

	return coord;
}

void Reader_c::check()
{
	cout << "-------Parametres de la simulation--------\n";
	cout << "ndime : ";cout << ndime;cout << "\n";
	cout << "nelem : ";cout << nelem;cout << "\n";
	cout << "npoint : ";cout << npoint;cout << "\n";
	cout << "nhalo : ";cout << nhalo;cout << "\n";
	cout << "elem : ";cout << elem;cout << "\n";
	cout << "point : ";cout << point;cout << "\n";
	cout << "line : ";cout << line;cout << "\n";
	cout << "bc : ";cout << bc;cout << "\n";
	cout << "\nElem2Node :\n";
	for(int i=0; i<elem;i++)
	{
		for(int j=0; j<vtklookup[1][elem2vtk[i]][1]; j++)
		{
			cout << elem2node[i][j]; cout << " ; ";
		}
		cout << "\n";
	}
	cout << "\ncoord :\n";
	for(int i=0; i<npoint; i++)
	{
		cout << coord[i][0]; cout << " ; "; cout << coord[i][1]; cout << " ; "; cout << coord[i][2]; cout << "\n";
	}
	cout << "\nelem2vtk :\n";
	for(int i=0; i<nelem; i++)
	{
		cout << elem2vtk[i];cout << "\n";
	}
}
void Reader_c::write_file(Reader_c& read, Solver_c& solve, int izone) {
	string filename = "Zone"+to_string(izone)+".su2";

	ofstream outfile(filename, std::ios_base::binary | std::ios_base::out);
	if (outfile.is_open()) {
		setprecision(16);
		outfile << "% \n";
		outfile << "ZONE= " << to_string(izone) << "\n";
		outfile << "% \n";
		outfile << "% Problem dimension \n";
		outfile << "% \n";
		outfile << "NDIME= " << read.ndime << "\n";
		outfile << "% \n";

		//=========================================== ELEM2NODES ===========================================
		outfile << "% Inner element connectivity \n";
		outfile << "% \n";
		outfile << "NELEM= "<< solve.zone2nelem[izone] << "\n";
		for (int ielem = 0; ielem < solve.zone2nelem[izone]; ielem++) {
			string temp_connec = "";
			for (int inoel = 0; inoel < vtklookup[1][solve.elem2vtk[izone][ielem]][1]; inoel++) {
				temp_connec += "    ";
				temp_connec += to_string(solve.elem2node[izone][ielem][inoel]);
				
			}
			temp_connec += "    ";
			outfile << solve.elem2vtk[izone][ielem]<< temp_connec <<"\n";
		}
		//============================================== COORD ==============================================
		outfile << "% \n";
		outfile << "% Node coordinates \n";
		outfile << "% \n";
		outfile << "NPOIN= " << solve.zone2nnode[izone] << "\n";
		for (int i = 0; i < solve.zone2nnode[izone]; i++) {
			outfile << fixed;
			outfile << setprecision(16) << solve.zone2coord[izone][i][0] <<"    "<< solve.zone2coord[izone][i][1] << "    " << solve.zone2coord[izone][i][2] << "\n";
		}
		// ============================================= BOUNDARY =============================================
		outfile << "% \n";
		outfile << "% Boundary elements \n";
		outfile << "% \n";
		
		outfile << "NMARK= " << read.nbc << "\n";
		for (int ibc = 0; ibc < read.nbc; ibc++) { 
			//int jzone = solve.zone2ijzone[izone][ijzone];
			outfile << "MARKER_TAG= " << read.bound2tag[ibc] << "\n";
			int nghost = solve.zone2boundIndex[izone][ibc + 1] - solve.zone2boundIndex[izone][ibc];
			outfile << "MARKER_ELEMS= " << nghost <<"\n";
			int ielem1 = solve.zone2boundIndex[izone][ibc];
			int ielem2 = solve.zone2boundIndex[izone][ibc+1];
			for (int ighost = ielem1; ighost < ielem2; ighost++) {
				string temp_connec = "";
				int vtk = solve.elem2vtk[izone][ighost];
				int nnoel = vtklookup[read.ndime-2][vtk][1]; 
				for (int inoel = 0; inoel < nnoel; inoel++) {
					temp_connec += "    ";
					temp_connec += to_string(solve.elem2node[izone][ighost][inoel]);

				}
				outfile << vtk << temp_connec << "\n";
			}
		}

		//=========================================== ZONE BOUNDARY ===========================================
		outfile << "% \n";
		outfile << "% Zone Boundary elements \n";
		outfile << "% \n";
		
		outfile << "NZONE= " << solve.nzone -1 << "\n";
		int ighost = 0;
		for (int ijzone = 0; ijzone < solve.nzone -1; ijzone++) { 
			int jzone = solve.zone2ijzone[izone][ijzone];
			outfile << "ZONE_TAG= " << to_string(jzone) << "\n";

			int njzone = solve.zone2markelem[izone][ijzone];
			outfile << "ZONE_ELEMS= " << njzone <<"\n";

			int ielem1 = solve.zone2zoneIndex[izone][ijzone];
			int ielem2 = solve.zone2zoneIndex[izone][ijzone+1];
			
			for (int ielem = ielem1; ielem < ielem2; ielem++) {
				string temp_connec = "";
				int vtk = solve.elem2vtk[izone][ielem];
				int nnofa = vtklookup[read.ndime-2][vtk][1]; 

				for (int inofa = 0; inofa < nnofa + 1; inofa++) {
					temp_connec += "    ";
					temp_connec += to_string(solve.elem2node[izone][ielem][inofa]);

				}
				outfile << to_string(vtk) << temp_connec  <<  "\n";
			}
		}
		outfile.close();


	}

}
void Reader_c::WriteAllZoneFile(Reader_c& read,Solver_c& solve ){
	for (int izone = 0; izone < solve.nzone; izone++) {
		write_file(read, solve, izone);
	}
}

void Reader_c::write_tecplot_METIS(Reader_c & read, Solver_c& solve){
	

	fstream outFile;
	outFile.open("METISZONE.dat", ios::out);
	outFile << "VARIABLES=\"X\",\"Y\",\"Z\",\"Zone\"" << endl;
	//outFile << "VARIABLES=\"X\",\"Y\",\"P\",\"U\",\"V\"" << endl;
	//Zone Carre pour le tecplot 

	outFile << "ZONE T=\"Element0\"" << endl; //Changer le nbr elements
	outFile << "Nodes=" << solve.nnode_g << ", Elements=" << solve.nelem_g << ", ZONETYPE=FEBRICK" << endl;
	outFile << "DATAPACKING=BLOCK" << endl;
	outFile << "VARLOCATION = ([4] = CELLCENTERED)" << endl;

	string a;                      //ecrire les coordonnees de laxe x a la suite
	for (int j = 0; j < solve.nnode_g; j++)
	{
		a = to_string(read.coord[j][0]);
		outFile << a << endl;
	}
	string b;                      // ecrire les coordonnees de laxe y a la suite
	for (int i = 0; i <= solve.nnode_g - 1; i++)
	{
		b = to_string(read.coord[i][1]);
		outFile << b << endl;
	}
	string z;                      // ecrire les coordonnees de laxe y a la suite
	for (int i = 0; i <= solve.nnode_g - 1; i++)
	{
		z = to_string(read.coord[i][2]);
		outFile << z << endl;
	}

	// Zone
	string m;
	for (int j = 0; j <= solve.nelem_g - 1; j++)
	{
		m = to_string(solve.elem2zone[j]);
		outFile << m << endl;
	}

		// Ecriture des noeuds de chaque elements pour les carres

	for (int ielem = 0; ielem <= solve.nelem_g - 1; ielem++)
	{
		int vtk = read.elem2vtk[ielem];
		int nnoel = vtklookup[ndime-2][vtk][1];
		
		for (int icol = 0; icol <= nnoel - 1; icol++)
		{
			string icols = to_string(read.elem2node[ielem][icol]+1);
			outFile << icols << " ";
		}
		outFile << endl;
		
	}


	outFile.close();



}