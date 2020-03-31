// Execute with CMD :  g++ *.cpp -lgomp -o out -I /home/charles/Documents/metis-5.1.0/include -L /home/charles/Documents/metis-5.1.0/lib -lmetis


#include <string>
#include <cctype>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "omp.h"
#include <chrono>
#include <time.h>

#include "main.h"
#include "Reader.h"
#include "Connect.h"

//#include "./include/TECXXX.h"
#include <metis.h>

using namespace std;

void compute(string simName, string su2, int partition)
{
	cout << "\n==================================== ONERA ==================================== " << endl;
  size_t lastindex = su2.find_last_of(".");
  string out = su2.substr(0, lastindex) + "_" + "zone";
	cout << " Simulation name : "; cout << simName << endl;
	cout << " Name of su2++ output files : "; cout << out << endl;
	cout << " Su2 file to partition : "; cout << su2 << endl;
	cout << " Nomber of partition to make : "; cout << partition << endl;
	cout << endl;
  Reader_c FileONERA;
	Connect_c ONERA;
	FileONERA.read_file(su2);
	ONERA.ComputeGlobalConnectivity(FileONERA);
	ONERA.ComputeMETIS(partition, 3, FileONERA);
	ONERA.ComputeZoneConnectivity(FileONERA);
	FileONERA.WriteAllZoneFile(out, FileONERA, ONERA);
}
// ======================================================================================================================================== //
// ================================================================= MAIN ================================================================= //
// ======================================================================================================================================== //
int main(int argc, char *argv[])
{
	if(argc<4)	{		cout << argc<<endl;		cout << "Input arguments are missing"<<endl;		exit( 666 );	}
	else	{		cout << "Number of input argument : "; cout << argc << endl;	}
	cout << "==================================== NOICE STARTING ====================================" << endl;
	string simName = argv[1];
	string su2 = argv[2];
	int partition = atoi(argv[3]);
	compute(simName, su2, partition);
	cout << "===================================== NOICE ENDING =====================================" << endl;
	return 0;
}
