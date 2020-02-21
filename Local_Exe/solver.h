#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <Comm.h>
#include <supposed_metric_n_connec.h>

using namespace std;

class solver : public supposed_metric_n_connec {    // TBD wheter public or private
    public:
    solver();                   // Constructeur
    ~solver();                  // Destructeur
    void Compute();             // solve problem

    private:
    //Private variables
    Comm World;                 // MPI Cluster
    double* rho;
    double* u;
    double* v;
    double* w;
    double* p;
    double** flux_c;            //convective flux  flux[iface][variable (rho=0,u,v,w,p=4)]
    double** flux_d;            //dissipative flux flux[iface][variable (rho=0,u,v,w,p=4)]
    double** residu_c;          //convective residuals residu[ielem][variable (rho=0,u,v,w,p=4)]
    double** residu_d;          //dissipative residuals residu[ielem][variable (rho=0,u,v,w,p=4)]
    double** residu_d_hyb;      //dissipative residuals residu[ielem][variable (rho=0,u,v,w,p=4)]
    double** W_0;               //copy of initial state for RK
    double** primitivesSendBuffer;  //Buffer for Tx
    double inf_speed,inf_speed_x,inf_speed_y,inf_speed_z;   // speed norm at infinity
    double RKM_coef[6][2][5] =     // RKM 5 coeffs [Stages][order-1][step]
    {   {{0,0,0,0,0                 },{0,0,0,0,0}},                     //Order 0 ; Filler
        {{1.0,0,0,0,0               },{1.0,0,0,0,0}},  //Order 1 ; Euler
        {{0.0,1.0,0,0,0                 },{0.0,1.0,0,0,0}},  //Order 2 ; Not in blasek
        {{0.1481,0.4000,1.0,0,0         },{0.1918,0.4929,1.0,0,0}},
        {{0.0833,0.2069,0.4265,1.0,0    },{0.1084,0.2602,0.5052,1.0,0}},
        {{0.0533,0.1263,0.2375,0.4414,1.0},{0.0695,0.1602,0.2898,0.5060,1.0}}};
    double RKH_coef[6][2][5] =     // RKH [Stages][0=alpha,1=beta][step]
    {   {{0,0,0,0,0                     },{0,0,0,0,0}},                     //Order 0 ; Filler
        {{0,0,0,0,0                     },{0,0,0,0,0}},       
        {{0,0,0,0,0                     },{0,0,0,0,0}},  
        {{0,0,0,0,0                     },{0,0,0,0,0}},
        {{0,0,0,0,0                     },{0,0,0,0,0}},
        {{0.2742,0.2069,0.5020,0.5142,1.0},{1.0,0.0,0.56,0.0,0.44}}};   //Only order 5 in blasek

    //Private methods
    void Initialisation();      // Initialise field to infinity
    void UpdateBound();         // Update boundary conditions
    void ExchangePrimitive();   // MPI exchange primitive values
    void ExchangeGradiants();   // MPI exchange gradiant values
    void TimeStepEul();         // Euler explicit time integration
    void TimeStepRkM();         // Runge-Kutta Multistage time integration
    void TimeStepRkH();         // Runge-Kutta Hybride time integration
    void ComputeFluxO1();       // Calcul des flux (Roe) ordre 1
	void ComputeFluxO1Conv() ;  // Calcul flux convectifs Ordre 1
    void ComputeFluxO2();       // Calcul des flux (Roe) ordre 2
    void InitMPI();             // Declare necessary structures for MPI
    //void computeGradiants();  // [May not be needed, depending on choosen implem]
    void ComputeResidu();       // Calcul des résidu
    double CheckConvergence();  // Somme des résidus
	double P2E(double);         // Calcul pression vers Energie
	double E2P(double);         // Calcul  Energie vers pression
	
};