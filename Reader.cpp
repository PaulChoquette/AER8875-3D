#include <iostream>
#include <sstream>
#include <fstream>
#include <string> 
#include "Reader.h" 
#include "main.h"
using namespace std;

void Reader::read_file(string filename) {

	//Initialize mesh constants
	ndime = 0; nelem = 0; npoint = 0; nhalo = 0;

	//Initialize information reading counters;
	elem = 0; point = 0; line = ""; bc = 0;

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
					elem2node_nh = new unsigned* [nelem];
					elem2vtk_nh = new unsigned[nelem];
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
					for (unsigned i = 0; i < npoint; i++) {
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
					bclinen = linen;
					BoundIndex = new unsigned[nbc + 1];
					BoundIndex[0] = nelem + 1;
					bc_elem2node = new unsigned** [nbc];
					bc_elem2vtk = new unsigned* [nbc];
					bc_nelemv = new unsigned[nbc];
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
						step++;
					}
					else if (step == 1) {
						bc_nelem = Readcnst(line, "MARKER_ELEMS= ");
						bc_nelemv[bc] = bc_nelem;
						BoundIndex[bc + 1] =BoundIndex[bc]+ bc_nelem;
						bcnl += bc_nelem; //add nelem since each elem has 1 line
						bc_elem2vtk[bc] = new unsigned [bc_nelem];
						bc_elem2node[bc] = new unsigned* [bc_nelem];

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

						elem2node = new unsigned* [ncell];
						elem2vtk = new unsigned [ncell];
						for (unsigned ielem = 0; ielem < nelem; ielem++)
						{
							elem2vtk[ielem] = elem2vtk_nh[ielem];
							elem2node[ielem] = new unsigned[vtklookup[elem2vtk[ielem]][1]];

							for (unsigned inode = 0; inode < vtklookup[elem2vtk[ielem]][1]; inode++)
							{		
								elem2node[ielem][inode] = elem2node_nh[ielem][inode];
							}
						}

						int imelem = 0;
						for (unsigned ibc = 0; ibc < nbc; ibc++)
						{
							for (unsigned ibcen = 0; ibcen < bc_nelemv[ibc]; ibcen++)
							{
								elem2vtk[nelem + imelem] = bc_elem2vtk[ibc][ibcen];
								elem2node[nelem + imelem] = new unsigned[vtklookup[elem2vtk[nelem + imelem]][1]];

								for (unsigned j = 0; j < vtklookup[elem2vtk[nelem + imelem]][1]; j++)
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
}

bool Reader::OpenFile(string filename)
{
	file.open(filename);

	if (file.is_open()) {
		return true;
	}
	else {
		return false;
	}
}

unsigned Reader::Readcnst(const string& line, const string& tofind)
{
	size_t cnstpos = line.find(tofind);
	unsigned cnst = 0;
	if (cnstpos != string::npos) {
		cnstpos = line.find_last_of(" ") + 1;

		string cnstch;

		cnstch = line.substr(cnstpos, line.length() - cnstpos);
		cnst = stoul(cnstch);
	}
	return cnst;
}

void Reader::Fill_E2N_VTK(const char* cline) {
	//This function reads a character line and extracts VTK index as well as allocating and storing elem2node rows
	char* end;
	unsigned j = 0;
	for (unsigned c = strtoul(cline, &end, 10); cline != end; c = strtoul(cline, &end, 10))
	{
		cline = end;

		//Where j is the integer counter. 0 is the vtk index and the rest is the 
		if (j == 0) {
			elem2vtk_nh[elem] = c;
			elem2node_nh[elem] = new unsigned [vtklookup[c][1]];
		}

		else if(j > 0 && j <= vtklookup[elem2vtk_nh[elem]][1])
		{
			elem2node_nh[elem][j - 1] = c + 1;
		}

		else {
			break;
		}

		j++;

	}

	elem++;
}


void Reader::Fill_BC_E2N_VTK(const char* cline, unsigned bc) {
	//This function reads a character line and extracts VTK index as well as allocating and storing elem2node rows
	char* end;
	unsigned j = 0;
	for (unsigned c = strtoul(cline, &end, 10); cline != end; c = strtoul(cline, &end, 10))
	{
		cline = end;

		//Where j is the integer counter. 0 is the vtk index and the rest is the 
		if (j == 0) {
			bc_elem2vtk[bc][bc_e2n_counter] = c;
			bc_elem2node[bc][bc_e2n_counter] = new unsigned[vtklookup[c][1]];
		}

		else if (j > 0 && j <= vtklookup[bc_elem2vtk[bc][bc_e2n_counter]][1])
		{
			bc_elem2node[bc][bc_e2n_counter][j - 1] = c + 1;
		}

		else {
			break;
		}

		j++;
	}
	bc_e2n_counter++;
}

double** Reader::Fill_coord(const char* cline) {
	char* end;
	unsigned j = 0;
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