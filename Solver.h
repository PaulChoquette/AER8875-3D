#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <Comm.h>
#include <Metric.h>

using namespace std;

class solver_c : public Metric_c {    // TBD wheter public or private
    public:
    solver_c(Reader_c&, double);                   // Constructeur
    ~solver_c();                                // Destructeur
    void InitMPIBuffer(Reader_c&);              // Declare necessary structures for MPI
    void Compute();                             // solve problem
    void PrintStylz();
    Comm World;                                 // MPI Cluster [World.world_rank to get thread ID]
    double invrho;
    double eTempo;
    // Simulation parametrisation
    string smoothOrNah;
    double mach,AoA,cfl,convergeCrit;
    int Order,RK_step,iterMax;
    double residuRel,convergeFixLimit;          // Residu, limit at which limitors will be fixed
    bool RK_M;
    double* rho;
    double* u;
    double* v;
    double* w;
    double* p;
    int iteration;
    //Debugging Methods
    void PrintPress();
    void HighlightZoneBorder(); //Set density in MPI-received zone borders to 0
    void SetAnalyticalGradiant(double,double,double);
    void PrintGradiant();
    void LimitTecplot();



    private:
    //Private variables
    int nbc;
    string* bound2tag;
    int* BoundIndex;            //Limits of boundary cells id -> [izone] = start, indexes up to n+1
    int* ZBoundIndex;           //Limits of zone boundary cells id -> [izone] = start, indexes up to n+1
    int ntgt;                   //Number of MPI targets
    int* elem2vtk;              //Overwrites the partially-deleted one in connect_c
    double** flux_c;            //convective F_lux  F_lux[iface][variable (rho=0,u,v,w,p=4)]
    double** flux_d;            //dissipative F_lux F_lux[iface][variable (rho=0,u,v,w,p=4)]
    double** residu_c;          //convective residuals residu[ielem][variable (rho=0,u,v,w,p=4)]
    double** F_lux;
    double** SmootyRezi;
    double** residu_d;          //dissipative residuals residu[ielem][variable (rho=0,u,v,w,p=4)]
    double** residu_d_hyb;      //dissipative residuals residu[ielem][variable (rho=0,u,v,w,p=4)]
    double** cons_;
    double** cons;               //copy of initial state for RK
    double** limit;             //Limitors per element and primary value
    double** primitivesSendBuffer;  //Buffer for Tx
    double** gradientSendBuffer;    //Gradient buffer for Tx
    double*** gradient;             // Array containing gradiant [ielem][dimention][conservative]
    double inf_speed,inf_speed_x,inf_speed_y,inf_speed_z;   // speed norm at infinity
    double RKM_coef[6][2][5] =     // RKM 5 coeffs [Stages][order-1][step]
    {   {{0,0,0,0,0                     },{0,0,0,0,0}},                     //Order 0 ; Filler
        {{1.0,0,0,0,0                   },{1.0,0,0,0,0}},                   //Order 1 ; Euler
        {{0.0,1.0,0,0,0                 },{0.0,1.0,0,0,0}},                 //Order 2 ; Not in blasek
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
    void ExchangeMetrics();     // MPI exchange needed metrics for order 2
    void ExchangePrimitive();   // MPI exchange primitive values
    void ExchangeGradiants();   // MPI exchange gradiant values
    void UpdateBound();         // Update boundary conditions
    void TimeStepEul();         // Euler explicit time integration
    void TimeStepRkM();         // Runge-Kutta Multistage time integration
    void TimeStepRkH();         // Runge-Kutta Hybride time integration
    void ComputeGrandientsNLimit();// Computes gradiants WITH limitors INCLUDED. Bool true
    void UpwindFlux(int,double,double,double,double,double,double,double,double,double,double);// Compute dissipative F_lux [local]
    void RoeDissipation(int,double,double,double,double,double,double,double,double,double,double);// Compute dissipative F_lux [local]
    void ComputeFluxO1();       // Calcul des F_lux (Roe) ordre 1
	void ComputeFluxO1Conv();   // Calcul F_lux convectifs Ordre 1
    void ComputeFluxO2();       // Calcul des F_lux (Roe) ordre 2
    void ComputeFluxO2Conv();   // Calcul F_lux convectifs Ordre 2
    void ResidualSmoothing();   //
    void ComputeResidu();       // Calcul des résidu
    void ComputeResiduConv();   // Calcul des résidus convectifs seulement
    double CheckConvergence();  // Somme des résidus
	double P2E(double,double,double,double,double);// Calcul pression vers Energie
	double E2P(double,double,double,double,double);// Calcul  Energie vers pression
    void WriteResidu();
    void SaveConservative();
    void SavePrimitive(int ielem);
    void SavePrimitiveRK(int ielem);
    void SaveFlux();
    void SaveFlux_Hyb();
};
