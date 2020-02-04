#include <mpi.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>

using namespace std;

class Comm {
    public :
    int world_rank; //Zone ID
    int world_size; //Number of zones
    int nRnd;       //Number of communication rounds
    int nTgt;       //Number of target zones
    vector<int> TgtList;    //List of frontiere Zones
    string FileName;        //Comm Schedule file name

    // Public methods
    Comm(string);
    ~Comm();
    void Build_Rx2Loc(void);
    double* Assign2Local(double**,double**);                //Store received values in the 2cnd array as per Rx2Loc
    int** ExchangeInt(int**,int**,int*);                                                  //Echange 1 int array (column) to all neighbors
    double** ExchangeDouble(double**,double**,int*);                                      //Echange 1 double array (column) to all neighbors
    double*** ExchangeDoubleStack(double***TxArray,double***RxArray,int StackLen,int*Len);//Echange [StackLen] double array (column) to all neighbors

    private :
    vector<int> LocSched;   //Local communication schedule
    int** Rx2Loc;           //Map from received position to local position / target
    // Private methods
    void Init(void);        //Reads FileName and builds LocSched
};