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
    for (int izone=0;izone<ntgt;++izone) {
        delete[] gradientBuffer[izone];
        delete[] primitivesBuffer[izone];
        delete[] rxOrder2localOrder[izone];
        delete[] oneDBuffer[izone];
    }
    delete[] tgtList;
    delete[] gradientBuffer;
    delete[] primitivesBuffer;
    delete[] zone2nbelem;
    delete[] rxOrder2localOrder;
    delete[] oneDBuffer;
    //Kill MPI process
    MPI_Finalize();
}

void Comm::Init() {
    // Start MPI
    MPI_Init(NULL,NULL);    // Init MPI
    MPI_Comm_size(MPI_COMM_WORLD,&world_size); // Get world size
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);  //Get world ID
}

void Comm::InitBuffer(int ntgt_in,int* tgtList_in,int* zone2nbelem_in) {
    ntgt = ntgt_in;
    tgtList = new int[ntgt];
    for (int itgt=0;itgt<ntgt;++itgt) {
        tgtList[itgt] = tgtList_in[itgt];
    }

    //Initialise primitive buffers
    zone2nbelem = new int[ntgt];
    for (int i=0;i<ntgt;++i) {
        zone2nbelem[i] = zone2nbelem_in[i];
    }
    gradientBuffer =  new double*[ntgt];
    primitivesBuffer =  new double*[ntgt];
    oneDBuffer =  new double*[ntgt];
    rxOrder2localOrder = new int*[ntgt];
    for (int izone=0;izone<ntgt;++izone) {
        gradientBuffer[izone] = new double[zone2nbelem[izone]*15];
        primitivesBuffer[izone] = new double[zone2nbelem[izone]*5];
        rxOrder2localOrder[izone] = new int[zone2nbelem[izone]];
        oneDBuffer[izone] = new double[zone2nbelem[izone]*3];
    }
}

void Comm::ExchangeCellOrder(int** localBorderID) {
    // Init requests
    MPI_Request Request[int(ntgt*2)];
    MPI_Barrier(MPI_COMM_WORLD);
    // Send Buffer
    for (int itgt=0;itgt<ntgt;++itgt) {
        //cout<<zone2nbelem[itgt]<<endl;
        MPI_Isend(localBorderID[itgt],zone2nbelem[itgt],MPI_INT,tgtList[itgt],1,MPI_COMM_WORLD,&Request[itgt]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int itgt=0;itgt<ntgt;++itgt) {
        MPI_Irecv(rxOrder2localOrder[itgt],zone2nbelem[itgt],MPI_INT,tgtList[itgt],1,MPI_COMM_WORLD,&Request[itgt+ntgt]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Waitall(ntgt*2,Request,MPI_STATUSES_IGNORE);
}

void Comm::ExchangePrimitives(double** primitivesTx) {
    // Init requests
    MPI_Request Request[int(ntgt*2)];
    MPI_Barrier(MPI_COMM_WORLD);
    int tempSize;
    // Send Buffer
    for (int itgt=0;itgt<ntgt;++itgt) {
        tempSize = zone2nbelem[itgt]*5;
        MPI_Isend(primitivesTx[itgt],tempSize,MPI_DOUBLE,tgtList[itgt],2,MPI_COMM_WORLD,&Request[itgt]);
    }
    for (int itgt=0;itgt<ntgt;++itgt) {
        tempSize = zone2nbelem[itgt]*5;
        MPI_Irecv(primitivesBuffer[itgt],tempSize,MPI_DOUBLE,tgtList[itgt],2,MPI_COMM_WORLD,&Request[itgt+ntgt]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Waitall(ntgt*2,Request,MPI_STATUSES_IGNORE);

}

void Comm::ReclassPrimitives(double**rho,double**u,double**v,double**w,double**p) {
    for (int izone=0;izone<ntgt;++izone) {
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
    MPI_Request Request[int(ntgt*2)];
    MPI_Barrier(MPI_COMM_WORLD);
    int tempSize;

    // Send Buffer
    for (int itgt=0;itgt<ntgt;++itgt) {
        tempSize = zone2nbelem[itgt]*15;
        MPI_Isend(gradientsTx[itgt],tempSize,MPI_DOUBLE,tgtList[itgt],3,MPI_COMM_WORLD,&Request[itgt]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int itgt=0;itgt<ntgt;++itgt) {
        tempSize = zone2nbelem[itgt]*15;
        MPI_Irecv(gradientBuffer[itgt],tempSize,MPI_DOUBLE,tgtList[itgt],3,MPI_COMM_WORLD,&Request[itgt+ntgt]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Waitall(ntgt*2,Request,MPI_STATUSES_IGNORE);
}

void Comm::Exchange1DBuffer(double ** Send1DBuff){
        // Init requests
    MPI_Request Request[int(ntgt*2)];
    MPI_Barrier(MPI_COMM_WORLD);
    int tempSize;
    // Send Buffer
    for (int itgt=0;itgt<ntgt;++itgt) {
        tempSize = zone2nbelem[itgt]*3;
        MPI_Isend(Send1DBuff[itgt],tempSize,MPI_DOUBLE,tgtList[itgt],42,MPI_COMM_WORLD,&Request[itgt]);
    }
    for (int itgt=0;itgt<ntgt;++itgt) {
        tempSize = zone2nbelem[itgt]*3;
        MPI_Irecv(oneDBuffer[itgt],tempSize,MPI_DOUBLE,tgtList[itgt],42,MPI_COMM_WORLD,&Request[itgt+ntgt]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Waitall(ntgt*2,Request,MPI_STATUSES_IGNORE);
}

double Comm::UpdateConvergence(double LocalRes){
    double GlobalRes,TempRes;
    if (world_rank!=0) {
        MPI_Request Request;
        MPI_Isend(&LocalRes,1,MPI_DOUBLE,0,4,MPI_COMM_WORLD,&Request);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Wait(&Request,MPI_STATUS_IGNORE);
        MPI_Irecv(&GlobalRes,1,MPI_DOUBLE,0,5,MPI_COMM_WORLD,&Request);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Wait(&Request,MPI_STATUS_IGNORE);
        SumResidu = GlobalRes;
    }
    else {
        MPI_Request Request[world_size-1];double resArr[world_size-1];

        GlobalRes = LocalRes;
        for (int itgt=1;itgt<world_size;++itgt) {
            MPI_Irecv(&resArr[itgt-1],1,MPI_DOUBLE,itgt,4,MPI_COMM_WORLD,&Request[itgt-1]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Wait(Request,MPI_STATUS_IGNORE);
        for (int itgt=1;itgt<world_size;++itgt) {
            GlobalRes += resArr[itgt-1];
        }
        SumResidu = GlobalRes;
        for (int itgt=1;itgt<world_size;++itgt) {
            MPI_Isend(&SumResidu,1,MPI_DOUBLE,itgt,5,MPI_COMM_WORLD,&Request[itgt-1]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Wait(Request,MPI_STATUS_IGNORE);
    }
    return SumResidu;
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
double Comm::UpdateCoefficient(double LocalForce){
    double GlobalForce,TempForce,SumForce;
    if (world_rank!=0) {
        MPI_Request Request;
        MPI_Isend(&LocalForce,1,MPI_DOUBLE,0,6,MPI_COMM_WORLD,&Request);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Wait(&Request,MPI_STATUS_IGNORE);
        MPI_Irecv(&GlobalForce,1,MPI_DOUBLE,0,7,MPI_COMM_WORLD,&Request);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Wait(&Request,MPI_STATUS_IGNORE);
        SumForce = GlobalForce;
    }
    else {
        MPI_Request Request[world_size-1];double resArr[world_size-1];

        GlobalForce = LocalForce;
        for (int itgt=1;itgt<world_size;++itgt) {
            MPI_Irecv(&resArr[itgt-1],1,MPI_DOUBLE,itgt,6,MPI_COMM_WORLD,&Request[itgt-1]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Wait(Request,MPI_STATUS_IGNORE);
        for (int itgt=1;itgt<world_size;++itgt) {
            GlobalForce += resArr[itgt-1];
        }
        SumForce = GlobalForce;
        for (int itgt=1;itgt<world_size;++itgt) {
            MPI_Isend(&SumForce,1,MPI_DOUBLE,itgt,7,MPI_COMM_WORLD,&Request[itgt-1]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Wait(Request,MPI_STATUS_IGNORE);
    }

    return SumForce;
}
