#include <mpi.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <Comm.h>

using namespace std;

Comm::Comm(string FileName_in) {
    FileName = FileName_in;
    Init();
}

Comm::~Comm() {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    for (int i=0;i<nTgt;++i){
        delete[] Rx2Loc[i];
    }
    delete[] Rx2Loc;
}

void Comm::Build_Rx2Loc(void) {
    Rx2Loc = new int*[nTgt];
    for (int i=0;i<nTgt;++i){
        Rx2Loc[i] = new int[nTgt];
        //To build according to number of elements
    }
}

void Comm::Init() {
    // Start MPI 
    MPI_Init(NULL,NULL);    // Init MPI

    MPI_Comm_size(MPI_COMM_WORLD,&world_size); // Get world size

    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);  //Get ID
    //cout<<"Thread "<<world_rank<<" running!"<<endl;

    //Fetch Comm Schedule
    TgtList.resize(0);
    LocSched.resize(0);
    nTgt = 0;
    int SpaceIndx;
    int FirstIndx;
    int LastIndx;
    bool isDone = 0;
    string line;
    ifstream file (FileName);
    if (file.is_open()) {
        while (getline(file,line)) {
            if (line.substr(0,4)=="RND ") {
                nRnd = stoi(line.substr(3));
                LocSched.push_back(-1);
                isDone = 0;
            }
            else if (!isDone) {
                SpaceIndx = line.find_first_of(" ");
                FirstIndx = stoi(line.substr(0,SpaceIndx+1));
                LastIndx = stoi(line.substr(SpaceIndx+1));
                if (FirstIndx==world_rank) {
                    LocSched[nRnd] = LastIndx;
                    TgtList.push_back(LastIndx);
                    nTgt += 1;
                    isDone = 1;

                }
                else if (LastIndx==world_rank) {
                    LocSched[nRnd] = FirstIndx;
                    TgtList.push_back(FirstIndx);
                    nTgt += 1;
                    isDone = 1;
                }
            }
        }
        file.close();
    }
    else {
        cout<<"ERROR : UNABLE TO OPEN "<<FileName<<endl;
    }
    nRnd += 1;
}

// Ugly, but works...
int** Comm::ExchangeInt(int**TxArray,int**RxArray,int*Len) {
    int Rnd = 0;
    int TgtNum = 0;
    bool isTx;
    while (Rnd<nRnd) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (LocSched[Rnd]!=-1) {    //If talks to someone 
        isTx = (world_rank==max(world_rank,LocSched[Rnd])); // High ID Tx first
            if (isTx) {
                //      Data_for_target Pointer_size Type    Target
                MPI_Send(TxArray[TgtNum],Len[TgtNum],MPI_INT,LocSched[Rnd],0,MPI_COMM_WORLD);
                MPI_Recv(RxArray[TgtNum],Len[TgtNum],MPI_INT,LocSched[Rnd],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            else {
                MPI_Recv(RxArray[TgtNum],Len[TgtNum],MPI_INT,LocSched[Rnd],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Send(TxArray[TgtNum],Len[TgtNum],MPI_INT,LocSched[Rnd],0,MPI_COMM_WORLD);
            }
            TgtNum += 1;
        }
        Rnd += 1; 
    }
    return RxArray;
}

double** Comm::ExchangeDouble(double**TxArray,double**RxArray,int*Len) {
    int Rnd = 0;
    int TgtNum = 0;
    bool isTx;
    while (Rnd<nRnd) {
        if (LocSched[Rnd]!=-1) {
        isTx = (world_rank==max(world_rank,LocSched[Rnd])); // High Tx first
            if (isTx) {
                //cout<<world_rank<<" Sending to "<<LocSched[Rnd]<<endl;
                MPI_Send(TxArray[TgtNum],Len[TgtNum],MPI_DOUBLE,LocSched[Rnd],0,MPI_COMM_WORLD);
                MPI_Recv(RxArray[TgtNum],Len[TgtNum],MPI_DOUBLE,LocSched[Rnd],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            else {
                //cout<<world_rank<<" Receiving from "<<LocSched[Rnd]<<endl;
                MPI_Recv(RxArray[TgtNum],Len[TgtNum],MPI_DOUBLE,LocSched[Rnd],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Send(TxArray[TgtNum],Len[TgtNum],MPI_DOUBLE,LocSched[Rnd],0,MPI_COMM_WORLD);
            }
            TgtNum += 1;
        }
        //cout<<" BARRIER HIT [ "<<world_rank<<" ]"<<endl;
        MPI_Barrier(MPI_COMM_WORLD);
        Rnd += 1; 
    }
    return RxArray;
}

double*** Comm::ExchangeDoubleStack(double***TxArray,double***RxArray,int StackLen,int*Len) {
    // Array[Zone][Dimention/Variable][element]
    // For rho/U/V/W... :
    //Len = number of elements transmitted per zone
    //StackLen = number of variables sent
    int Rnd = 0;
    int TgtNum = 0;
    bool isTx;
    while (Rnd<nRnd) {
        if (LocSched[Rnd]!=-1) {
        isTx = (world_rank==max(world_rank,LocSched[Rnd])); // High Tx first
            if (isTx) {
                for (int i=0;i<StackLen;++i) {
                    MPI_Send(TxArray[TgtNum][i],Len[TgtNum],MPI_DOUBLE,LocSched[Rnd],0,MPI_COMM_WORLD);
                    MPI_Recv(RxArray[TgtNum][i],Len[TgtNum],MPI_DOUBLE,LocSched[Rnd],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                }
            }
            else {
                for (int i=0;i<StackLen;++i) {
                    MPI_Recv(RxArray[TgtNum][i],Len[TgtNum],MPI_DOUBLE,LocSched[Rnd],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    MPI_Send(TxArray[TgtNum][i],Len[TgtNum],MPI_DOUBLE,LocSched[Rnd],0,MPI_COMM_WORLD);
                }
            }
            TgtNum += 1;
        }
        //cout<<" BARRIER HIT [ "<<world_rank<<" ]"<<endl;
        MPI_Barrier(MPI_COMM_WORLD);
        Rnd += 1; 
    }
    return RxArray;
}


