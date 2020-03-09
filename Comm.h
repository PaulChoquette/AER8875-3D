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
    int* tgtList;               // List of frontiere Zones
    int* zone2nbelem;           // Number of cell per frontiere
    double** primitivesBuffer;  // Buffer for recieved primitives values
    double** oneDBuffer;        // Buffer for recieved primitives values
    double** gradientBuffer;    // Buffer for received gradient values
    int** rxOrder2localOrder;   // Local Index to store received values
    double SumResidu;
    

    // Public methods
    Comm();
    ~Comm();
    void Init();                    //Reads FileName and builds LocSched
    void InitBuffer(int,int*,int*); //Initialise buffers based on number of border cells with each zones
    void ExchangeCellOrder(int**);  //Exchanges local indexes with relevent zones
    void Exchange1DBuffer(double**);
    void ExchangePrimitives(double**);//Exchanges primitives with relevent zones
    void ReclassPrimitives(double**rho,double**u,double**v,double**w,double**p);    //Updates primitives with new values
    void ExchangeGradients(double**);//Exchanges gradients
    double UpdateConvergence(double);

    // Debugging Methods
    void PrintCellOrder(void);
    void PrintP(void);

    private :
    // Private methods
    
};