#include <mpi.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <Comm.h>

using namespace std;

Comm::Comm() {
}

Comm::~Comm() {
    MPI_Barrier(MPI_COMM_WORLD);
    //Delete Dynamic Variables
    for (int izone;izone<ntgt;++izone) {
        delete[] gradientBuffer[izone];
        delete[] primitivesBuffer[izone];
        delete[] rxOrder2localOrder[izone];
    }
    delete[] gradientBuffer;
    delete[] primitivesBuffer;
    delete[] zone2nbelem;
    delete[] rxOrder2localOrder;
    //Kill MPI process
    MPI_Finalize();
}

void Comm::Init(string FileName_in) {
    FileName = FileName_in;
    // Start MPI 
    MPI_Init(NULL,NULL);    // Init MPI
    MPI_Comm_size(MPI_COMM_WORLD,&world_size); // Get world size
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);  //Get world ID

    //Fetch Comm Schedule
    tgtList.resize(0);
    ntgt = 0;
    int SpaceIndx;
    int FirstIndx;
    int LastIndx;
    bool isDone = 0;
    string line;
    ifstream file (FileName);
    if (file.is_open()) {
        while (getline(file,line)) {
            if (line.substr(0,4)=="RND ") {
                isDone = 0;
            }
            else if (!isDone) {
                SpaceIndx = line.find_first_of(" ");
                FirstIndx = stoi(line.substr(0,SpaceIndx+1));
                LastIndx = stoi(line.substr(SpaceIndx+1));
                if (FirstIndx==world_rank) {
                    tgtList.push_back(LastIndx);
                    ntgt += 1;
                    isDone = 1;

                }
                else if (LastIndx==world_rank) {
                    tgtList.push_back(FirstIndx);
                    ntgt += 1;
                    isDone = 1;
                }
            }
        }
        file.close();
    }
    else {
        cout<<"ERROR : UNABLE TO OPEN "<<FileName<<endl;
    }
}

void Comm::InitBuffer(int* zone2nbelem_in) {
    //Initialise primitive buffers
    zone2nbelem = new int[ntgt];
    for (int i=0;i<ntgt;++i) {
        zone2nbelem[i] = zone2nbelem_in[i];
    }
    gradientBuffer =  new double*[ntgt];
    primitivesBuffer =  new double*[ntgt];
    rxOrder2localOrder = new int*[ntgt];
    for (int izone;izone<ntgt;++izone) {
        gradientBuffer[izone] = new double[zone2nbelem[izone]*15];
        primitivesBuffer[izone] = new double[zone2nbelem[izone]*5];
        rxOrder2localOrder[izone] = new int[zone2nbelem[izone]];
    }
}



void Comm::ExchangeCellOrder(int** localBorderID) {
    // Init requests
    MPI_Request send_request, recv_request;
    MPI_Status status;

    // Send Buffer
    for (int itgt=0;itgt<ntgt;++itgt) {
        MPI_Isend(localBorderID[itgt],zone2nbelem[itgt],MPI_INT,tgtList[itgt],1,MPI_COMM_WORLD,&send_request);
        MPI_Wait(&send_request,&status);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int itgt=0;itgt<ntgt;++itgt) {
        MPI_Irecv(rxOrder2localOrder[itgt],zone2nbelem[itgt],MPI_INT,tgtList[itgt],1,MPI_COMM_WORLD,&recv_request);
        MPI_Wait(&recv_request,&status);
    }
}

void Comm::ExchangePrimitives(double** primitivesTx) {
    // Init requests
    MPI_Request send_request, recv_request;
    MPI_Status status;
    int tempSize;
    
    // Send Buffer
    for (int itgt=0;itgt<ntgt;++itgt) {
        tempSize = zone2nbelem[itgt]*5;
        MPI_Isend(primitivesTx[itgt],tempSize,MPI_DOUBLE,tgtList[itgt],0,MPI_COMM_WORLD,&send_request);
        MPI_Wait(&send_request,&status);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int itgt=0;itgt<ntgt;++itgt) {
        tempSize = zone2nbelem[itgt]*5;
        MPI_Irecv(primitivesBuffer[itgt],tempSize,MPI_DOUBLE,tgtList[itgt],0,MPI_COMM_WORLD,&recv_request);
        MPI_Wait(&recv_request,&status);
    }
}

void Comm::ReclassPrimitives(double**rho,double**u,double**v,double**w,double**p) {
    for (int izone;izone<ntgt;++izone) {
        for (int ibelem=0;ibelem<zone2nbelem[izone];++ibelem) {       
            //*rho[rxOrder2localOrder[izone][ibelem]] = primitivesBuffer[izone][ibelem]; DOESN'T WORK
            // Yes, I like living on the edge 
            *(rho[0]+rxOrder2localOrder[izone][ibelem]) = primitivesBuffer[izone][ibelem];  //WORKS
            *(u[0]+rxOrder2localOrder[izone][ibelem]) = primitivesBuffer[izone][ibelem+zone2nbelem[izone]];
            *(v[0]+rxOrder2localOrder[izone][ibelem]) = primitivesBuffer[izone][ibelem+zone2nbelem[izone]*2];
            *(w[0]+rxOrder2localOrder[izone][ibelem]) = primitivesBuffer[izone][ibelem+zone2nbelem[izone]*3];
            *(p[0]+rxOrder2localOrder[izone][ibelem]) = primitivesBuffer[izone][ibelem+zone2nbelem[izone]*4];
        }
    }
}

void Comm::ExchangeGradients(double**gradientsTx) {
    // Init requests
    MPI_Request send_request, recv_request;
    MPI_Status status;
    int tempSize;

    // Send Buffer
    for (int itgt=0;itgt<ntgt;++itgt) {
        tempSize = zone2nbelem[itgt]*15;
        MPI_Isend(gradientsTx[itgt],tempSize,MPI_DOUBLE,tgtList[itgt],0,MPI_COMM_WORLD,&send_request);
        MPI_Wait(&send_request,&status);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int itgt=0;itgt<ntgt;++itgt) {
        tempSize = zone2nbelem[itgt]*15;
        MPI_Irecv(gradientBuffer[itgt],tempSize,MPI_DOUBLE,tgtList[itgt],0,MPI_COMM_WORLD,&recv_request);
        MPI_Wait(&recv_request,&status);
    }
}

void Comm::ExchangeMetrics(){
    
}

void Comm::PrintCellOrder(void){
    for (int izone=0;izone<ntgt;++izone) {
        for (int ibelem=0;ibelem<zone2nbelem[izone];++ibelem) {
            cout<<"T: "<<world_rank<<"| z:"<<tgtList[izone]<<" Lcell : "<<rxOrder2localOrder[izone][ibelem]<<endl;
        }
    }
}

void Comm::PrintP(void){
    for (int izone=0;izone<ntgt;++izone) {
        for (int ibelem=0;ibelem<zone2nbelem[izone];++ibelem) {
            cout<<"T: "<<world_rank<<"| z:"<<tgtList[izone]<<" Lcell : "<<primitivesBuffer[izone][ibelem+zone2nbelem[izone]*4]<<endl;
        }
    }
}