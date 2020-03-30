#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include "Reader.h"
#include "Solver.h"
#include "main.h"
using namespace std;

void Reader_c::read_file_local(string filename) {
	cout << "Reading of Local SU2 \tSTARTING...";
	//Initialize mesh constants
	ndime = 0; nelem = 0; npoint = 0; nhalo = 0;
	ncell = 0;
	//Initialize information reading counters;
	elem = 0; point = 0; line = ""; linen = 0; bc = 0; nbc = 0; bclinen = 0; nzone = 0; zlinen = 0;
	int bcempty = 0;
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
						//cout << "================ STEP 0 : " << line.substr(12,8) << endl;
						step++;
					}
					else if (step == 1) {
						bc_nelem = Readcnst(line, "MARKER_ELEMS= ");
						//cout << "================ Step 1 : bc_nelem =  " << bc_nelem << "  Line = " << line << endl;

						bc_nelemv[bc] = bc_nelem;
						BoundIndex[bc + 1] =BoundIndex[bc]+ bc_nelem;

						nhalo += bc_nelem;
						bc_elem2vtk[bc] = new int [bc_nelem];
						bc_elem2node[bc] = new int* [bc_nelem];

						bc_e2n_counter = 0; //counter for marker element number used in FillMarker
						step++;
						if (bc_nelem==0){
							step = 0;
							bc++;
							//bcnl += 1; //add nelem since each elem has 1 line
							continue;
						}
						else {

							bcnl += bc_nelem; //add nelem since each elem has 1 line
						}
					}
					else if (step == 2) {
						//cout << "================ Step 2 = " << endl;

						Fill_BC_E2N_VTK(cline, bc);
						if (bc_e2n_counter == bc_nelem) {
							bc++;
							step = 0;
						}
					}
				}
			}

			if (nzone == 0) {
				nzone = Readcnst(line, "NZONE= ");
				if (nzone != 0) {
					zone2tag = new string[nzone];
					zlinen = linen;
					zoneIndex = new int[nzone + 1];
					zoneIndex[0] = 0;                           //Indexation starts at 0
					z_elem2node = new int** [nzone];
					z_elem2vtk = new int* [nzone];
					pre_zelem2jelem = new int* [nzone];
					z_nelemv = new int[nzone];
					znl = 2 * nzone; //2 lines per marker for marker name and elemn
					zoneid = 0;
					step = 0;
					continue;
				}
			}
			if (zlinen != 0) {
				if (linen > zlinen && linen <= zlinen + znl) {

					//3 steps : step 0 is reading marker tag, step 1 is reading marker elemn and final step is filling bc_elem2nnode
					if (step == 0) {
						zone2tag[zoneid] = line.substr(10, 12);
						step++;
					}
					else if (step == 1) {
						z_nelem = Readcnst(line, "ZONE_ELEMS= ");
						z_nelemv[zoneid] = z_nelem;
						zoneIndex[zoneid + 1] = zoneIndex[zoneid] + z_nelem;
						znl += z_nelem; //add nelem since each elem has 1 line
						z_elem2vtk[zoneid] = new int[z_nelem];
						pre_zelem2jelem[zoneid] = new int[z_nelem];
						z_elem2node[zoneid] = new int* [z_nelem];

						z_e2n_counter = 0; //counter for marker element number used in FillMarker
						step++;
					}
					else if (step == 2) {
						Fill_ZN_E2N_VTK(cline, zoneid);
						if (z_e2n_counter == z_nelem) {
							zoneid++;
							step = 0;
						}
					}
					if (bc == nbc && zoneid == nzone) {
						nhalo = bcnl - 2 * nbc;
						nzelem = (znl - 2 * nzone);
						ncell = nelem + nhalo + nzelem;
						elem2node = new int* [ncell];
						elem2vtk = new int[ncell];
						zelem2jelem = new int[nzelem];
						for (int ielem = 0; ielem < nelem; ielem++)
						{
							elem2vtk[ielem] = elem2vtk_nh[ielem];
							elem2node[ielem] = new int[vtklookup[ndime - 2][elem2vtk[ielem]][1]];
							for (int inode = 0; inode < vtklookup[ndime - 2][elem2vtk[ielem]][1]; inode++)
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
								elem2node[nelem + imelem] = new int[vtklookup[ndime - 2][elem2vtk[nelem + imelem]][1]];

								for (int j = 0; j < vtklookup[ndime - 2][elem2vtk[nelem + imelem]][1]; j++)
								{
									elem2node[nelem + imelem][j] = bc_elem2node[ibc][ibcen][j];
								}
								imelem++;
							}
						}
						imelem = 0;
						for (int izone = 0; izone < nzone; izone++)
						{
							for (int izone_en = 0; izone_en < z_nelemv[izone]; izone_en++)
							{
								elem2vtk[nelem + nhalo + imelem] = z_elem2vtk[izone][izone_en];
								zelem2jelem[imelem] = pre_zelem2jelem[izone][izone_en];
								elem2node[nelem + nhalo + imelem] = new int[vtklookup[ndime - 2][elem2vtk[nelem + nhalo + imelem]][1]];

								for (int j = 0; j < vtklookup[ndime - 2][elem2vtk[nelem + nhalo + imelem]][1]; j++)
								{
									elem2node[nelem + nhalo + imelem][j] = z_elem2node[izone][izone_en][j];
								}
								imelem++;
							}
						}

					}

				}

			}


		}


		//Boundary condition array deletion

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

		//Zone 2 Zone array deletion

		for (int i = 0; i < nzone; i++) {
			for (int j = 0; j < z_nelemv[i]; j++) {
				delete[] z_elem2node[i][j];
			}
			delete[] z_elem2node[i];
		}
		delete[] z_elem2node;
		delete[] z_elem2vtk;
		delete[] elem2vtk_nh;
		file.close();
	}
	else {
		//ERROR 1 : File could not be opened
	}
cout << "...............DONE" << endl;
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

void Reader_c::Fill_ZN_E2N_VTK(const char* cline, int zoneid) {
	//This function reads a character line and extracts VTK index as well as allocating and storing elem2node rows
	char* end;
	int j = 0;
	for (int c = strtoul(cline, &end, 10); cline != end; c = strtoul(cline, &end, 10))
	{
		cline = end;
		//Where j is the integer counter. 0 is the vtk index and the rest is the
		if (j == 0) {
			z_elem2vtk[zoneid][z_e2n_counter] = c;
			z_elem2node[zoneid][z_e2n_counter] = new int[vtklookup[ndime - 2][c][1]];
		}

		else if (j > 0 && j <= vtklookup[ndime - 2][z_elem2vtk[zoneid][z_e2n_counter]][1])
		{
			z_elem2node[zoneid][z_e2n_counter][j - 1] = c;
		}

		else if (j == vtklookup[ndime - 2][z_elem2vtk[zoneid][z_e2n_counter]][1] + 1)
		{
			pre_zelem2jelem[zoneid][z_e2n_counter] = c;
		}

		else {
			break;
		}

		j++;
	}
	z_e2n_counter++;
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
	cout << "nzelem : "; cout << nzelem; cout << "\n";
	cout << "elem : ";cout << elem;cout << "\n";
	cout << "point : ";cout << point;cout << "\n";
	cout << "line : ";cout << line;cout << "\n";
	cout << "bc : ";cout << bc;cout << "\n";
	cout << "nzone : "; cout << nzone; cout << "\n";
	// cout << "\nElem2Node :\n";
	// for(int i=0; i<elem;i++)
	// {
	// 	for(int j=0; j<vtklookup[1][elem2vtk[i]][1]; j++)
	// 	{
	// 		cout << elem2node[i][j]; cout << " ; ";
	// 	}
	// 	cout << "\n";
	// }
	// cout << "\ncoord :\n";
	// for(int i=0; i<npoint; i++)
	// {
	// 	cout << coord[i][0]; cout << " ; "; cout << coord[i][1]; cout << " ; "; cout << coord[i][2]; cout << "\n";
	// }
	// cout << "\nelem2vtk :\n";
	// for(int i=0; i<ncell; i++)
	// {
	// 	cout << elem2vtk[i];cout << "\n";
	// }
	// cout << "\nzelem2jelem :\n";
	// for (int i = 0; i < nzelem; i++)
	// {
	// 	cout << zelem2jelem[i]; cout << "\n";
	// }
}


void Reader_c::computePrmt(string filename) {

	//Open the file and continue if file is succesfully opened
	if (OpenFile_2(filename))
	{
		SimName = "666";
		su2FilePath = "666";
		Npartition = 666;
    AoA = 666;
    mach = 666;
    gamma = 666;
		tempMethod = "666";
		Nstage = 666;
		cfl = 666;
		Smoothing = "666";
		spatMethod = "666";
		spatMethod_ordre = 666;
		iterMax = 666;
		convCrit = 666;
		AoA_i = 666;
		AoA_f = 666;
		//Read file line_2 by line_2
		while (getline(file_2, line_2))
		{
			//Tick line_2 number and store line_2 as character array
			linen_2++;
			cline_2 = line_2.c_str();
      //cout << cline_2 <<endl;
			if(SimName=="666"){
        SimName = inputStr(cline_2, "SimName= ");
				continue;
      }
			if(su2FilePath=="666"){
        su2FilePath = inputStr(cline_2, "su2FilePath= ");
				continue;
      }
			if(Npartition==666){
        Npartition = inputInt(cline_2, "NbrPartition= ");
				continue;
      }
			if(AoA==666){
        AoA = inputDouble(cline_2, "alpha= ");
				continue;
      }
      if(mach==666){
        mach = inputDouble(cline_2, "mach= ");
				continue;
      }
      if(gamma==666){
        gamma = inputDouble(cline_2, "gamma= ");
				continue;
      }
      if(tempMethod=="666"){
        tempMethod = inputStr(cline_2, "TemporelMethod= ");
				continue;
      }
			if(Nstage==666){
        Nstage = inputInt(cline_2, "stgaeNbr= ");
				continue;
      }
			if(cfl==666){
        cfl = inputDouble(cline_2, "cflVAR_str= ");
				continue;
      }
			if(Smoothing=="666"){
        Smoothing = inputStr(cline_2, "ResidualSmooth= ");
				continue;
      }
			if(spatMethod=="666"){
        spatMethod = inputStr(cline_2, "SpatialMethod= ");
				continue;
      }
			if(spatMethod_ordre==666){
        spatMethod_ordre = inputInt(cline_2, "SpatialMethodORDER= ");
				continue;
      }
			if(iterMax==666){
        iterMax = inputInt(cline_2, "iterMAX= ");
				continue;
      }
			if(convCrit==666){
        convCrit = inputDouble(cline_2, "convergenceCriterea= ");
				continue;
      }
			if(AoA_i==666){
        AoA_i = inputDouble(cline_2, "alpha_0= ");
				continue;
      }
			if(AoA_f==666){
        AoA_f = inputDouble(cline_2, "alpha_inf= ");
				continue;
      }
		}

		file_2.close();
	}
	else {
		//ERROR 1 : File could not be opened
	}
}
bool Reader_c::OpenFile_2(string filename)
{
	file_2.open(filename);

	if (file_2.is_open()) {
		return true;
	}
	else {
		return false;
	}
}
int Reader_c::inputInt(const string& line_2, const string& tofind)
{
	size_t cnstpos = line_2.find(tofind);
	int cnst = 0;
	if (cnstpos != string::npos) {
		//cout << "cnstpos : "; cout << cnstpos<<endl;
		cnstpos = line_2.find_last_of(" ") + 1;

		string cnstch;

		cnstch = line_2.substr(cnstpos, line_2.length() - cnstpos);
		cnst = stoul(cnstch);
		cout << tofind; cout << " "; cout << cnst<<endl;
	}
	else {
		return 666;
	}
	//cout << tofind; cout << " "; cout << cnst; cout << "\n";
	return cnst;
}
double Reader_c::inputDouble(const string& line_2, const string& tofind)
{
	size_t cnstpos = line_2.find(tofind);
	double cnst = 0;
	if (cnstpos != string::npos) {
		cnstpos = line_2.find_last_of(" ") + 1;

		string cnstch;

		cnstch = line_2.substr(cnstpos, line_2.length() - cnstpos);
		cnst = stod(cnstch);
		cout << tofind; cout << " "; cout << cnst<<endl;
	}
	else {
		return 666;
	}
	//cout << tofind; cout << " "; cout << cnst; cout << "\n";
	return cnst;
}

string Reader_c::inputStr(const string& line_2, const string& tofind)
{
	size_t cnstpos = line_2.find(tofind);
	string cnst = "0";
	if (cnstpos != string::npos) {
		cnstpos = line_2.find_last_of(" ") + 1;

		string cnstch;

		cnstch = line_2.substr(cnstpos, line_2.length() - cnstpos);
		cnst = cnstch;
		cout << tofind; cout << " "; cout << cnst<<endl;
	}
	else {
		return "666";
	}
	return cnst;
}
