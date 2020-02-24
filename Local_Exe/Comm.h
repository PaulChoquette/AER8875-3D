#include <mpi.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>

using namespace std;

// NOTE : tgtList HAS TO BE CLASSED ACCORDING TO THREAD'S ORDER


class Comm {
    public :
    int world_rank;             // Zone ID
    int world_size;             // Number of zones
    int ntgt;                   // Number of target zones
    vector<int> tgtList;        // List of frontiere Zones
    string FileName;            // Comm Schedule file name
    int* zone2nbelem;           // Number of cell per frontiere
    double** primitivesBuffer;  // Buffer for recieved primitives values
    double** gradientBuffer;    // Buffer for received gradient values
    int** rxOrder2localOrder;   // Local Index to store received values
    

    // Public methods
    Comm();
    ~Comm();
    void Init(string);              //Reads FileName and builds LocSched
    void InitBuffer(int*);          //Initialise buffers based on number of border cells with each zones
    void ExchangeCellOrder(int**);  //Exchanges local indexes with relevent zones
    void ExchangeMetrics();         //Exchances metrics (only cent2face) of borders
    void ExchangePrimitives(double**);//Exchanges primitives with relevent zones
    void ReclassPrimitives(double**rho,double**u,double**v,double**w,double**p);    //Updates primitives with new values
    void ExchangeGradients(double**);//Exchanges gradients

    // Debugging Methods
    void PrintCellOrder(void);
    void PrintP(void);

    private :
    // Private methods
    
};