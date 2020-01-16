#include <iostream>
#include <sstream>
#include <fstream>
#include <string> 
#include "Reader.h" 
#include "main.h"
using namespace std;

void Reader::read_file(string filename) {

	//Initialize mesh constants
	ndime = 0; nelem = 0; npoint = 0;

	//Initialize information reading counters;
	elem = 0; point = 0; line = "";

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
				ndime = Readcnst(line);
				continue;
			}

			if (nelem == 0)
			{
				//Read number of elements in the mesh
				nelem = Readcnst(line);

				//Define elem2node matrix row size if nelem has been defined
				//Store line where elem2node definitions start
				if (nelem) {
					elem2node = new unsigned* [nelem];
					vtk = new unsigned[nelem];
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
				npoint = Readcnst(line);

				//Define coordinate matrix if number of points has been defined
				if (npoint) 
				{
					npointlinen = linen;
					coord = new double* [nelem];
					for (int i = 0; i < nelem; i++) {
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
				}
			}
		}
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

unsigned Reader::Readcnst(const string& line)
{
	size_t cnstpos = line.find("= ");
	unsigned cnst = 0;
	if (cnstpos != string::npos) {
		char cnstch;
		const char* ccnstch = &cnstch;
		cnst = atoi(ccnstch);
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
			elem2vtk[elem] = c;
			elem2node[elem] = new unsigned [vtklookup[elem2vtk[elem]][1]];
		}

		else if(j > 0 && j <= vtklookup[elem2vtk[elem]][1])
		{
			elem2node[elem][j - 1] = c + 1;
		}

		else {
			break;
		}

		j++;

	}

	elem++;
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
		j++;
	}
	point++;

	return coord;
}