/* MPI DEMO 
The following example demonstrates a non-blocking implementation of
inter-zone communication. 
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <Comm.h>
#include <mpi.h>

using namespace std;

string FileName = "Communication_Schedule.txt";


//Simulation Settings
int nelem = 69;
int nbelem = 2;
int nzone;


//Dynamic Variables
double*rho,*u,*v,*w,*p;
int*zone2nbelem,**localOrder;
double**primitivesSendBuffer;



int main() 
{   
    Comm Cluster;
    Cluster.Init(FileName);
    nzone = Cluster.ntgt;

    //Init dynamic variables
    rho = new double[nelem];u = new double[nelem];v = new double[nelem];w = new double[nelem];p = new double[nelem];
    zone2nbelem = new int[nzone];
    
    primitivesSendBuffer = new double*[nzone];
    localOrder = new int*[nzone];
    for (int izone=0;izone<nzone;++izone) {
        primitivesSendBuffer[izone] = new double[nbelem*5];
        localOrder[izone] = new int[nbelem];
    }


    // Fill data arrays
    for (int ielem=0;ielem<nelem;++ielem) {
        rho[ielem] = double(Cluster.world_rank);
        u[ielem] = double(Cluster.world_rank);
        v[ielem] = double(Cluster.world_rank);
        w[ielem] = double(Cluster.world_rank);
        p[ielem] = double(Cluster.world_rank);
    }

    for (int izone=0;izone<nzone;++izone) {
        zone2nbelem[izone] = nbelem;
        for (int ibelem=0;ibelem<nbelem;++ibelem) {
            localOrder[izone][ibelem] = Cluster.world_rank + Cluster.world_rank + ibelem;
        }
        for (int ibelem=0;ibelem<nbelem;++ibelem) {
            primitivesSendBuffer[izone][ibelem*5] = double(Cluster.world_rank);
            primitivesSendBuffer[izone][ibelem*5+1] = double(Cluster.world_rank);
            primitivesSendBuffer[izone][ibelem*5+2] = double(Cluster.world_rank);
            primitivesSendBuffer[izone][ibelem*5+3] = double(Cluster.world_rank);
            primitivesSendBuffer[izone][ibelem*5+4] = double(Cluster.world_rank);
        }
    }
    Cluster.InitBuffer(zone2nbelem);
    cout<<"Thread : "<<Cluster.world_rank<<" | Initialisation successful"<<endl;
    
    Cluster.ExchangeCellOrder(localOrder);
    cout<<"Thread : "<<Cluster.world_rank<<" | Initial Handshake successful"<<endl;

    //Cluster.PrintCellOrder();

    Cluster.ExchangePrimitives(primitivesSendBuffer);
    //Cluster.PrintP();

    // Reclass primitives
    Cluster.ReclassPrimitives(&rho,&u,&v,&w,&p);

    //Print rho
    for (int ielem=0;ielem<nelem;++ielem) {
        cout<<"T:"<<Cluster.world_rank<<"| ie:"<<ielem<<" Val:"<<p[ielem]<<endl;
    }


    //Delete dynamic variables
    for (int izone=0;izone<nzone;++izone) {
        delete[] localOrder[izone];
        delete[] primitivesSendBuffer[izone];
    }
    delete[] rho;delete[] u;delete[] v;delete[] w;delete[] p;delete[] localOrder;delete[] primitivesSendBuffer;delete[] zone2nbelem;

    return 1;
};