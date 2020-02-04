#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <Comm.h>
#include <mpi.h>

using namespace std;

string FileName = "Communication_Schedule.txt";

int main() 
{   
    Comm Cluster(FileName);
    Cluster.Build_Rx2Loc();

    //General Declaration Required : 
    int DatLen = 2; //Basically the number of element/zone (set to fixed value here)
    int ** TxArray;
    int ** RxArray;
    int * Len;  //Number of boundary element / zone
    TxArray = new int*[Cluster.nTgt];
    RxArray = new int*[Cluster.nTgt];
    Len = new int[Cluster.nTgt];


    // Would be updated to the elements values
    for (int i=0 ; i<Cluster.nTgt; ++i) {   //For Zones
        Len[i] = DatLen;
        TxArray[i] = new int[DatLen];
        RxArray[i] = new int[DatLen];
        for (int j=0 ; j<DatLen; ++j) { //For Elements
            TxArray[i][j] = rand() %100;
            RxArray[i][j] = -1;
        }
    }

    // Communication
    RxArray = Cluster.ExchangeInt(TxArray,RxArray,Len);  


    // Printing
    for (int i=0 ; i<Cluster.nTgt; ++i) {
        for (int j=0 ; j<DatLen; ++j) {
            cout<<"[Thread "<<Cluster.world_rank<<"] Received : "<<RxArray[i][j]<<" from "<<Cluster.TgtList[i]<<endl;
        }
        //cout<<"[Thread "<<Cluster.world_rank<<"] Received from : "<<Cluster.TgtList[i]<<endl;
        delete[] RxArray[i];
        delete[] TxArray[i];
    }
    delete[] RxArray;
    delete[] TxArray;
    delete[] Len;

    return 1;
};