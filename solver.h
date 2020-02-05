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
    void compute();             // solve problem

    private:
    //Private variables
    double* rho;
    double* u;
    double* v;
    double* w;
    double* p;
    double** flux_c;    //convective flux  flux[iface][variable (rho=0,u,v,w,p=4)]
    double** flux_d;    //dissipative flux flux[iface][variable (rho=0,u,v,w,p=4)]
    double** residu_c;  //convective residuals residu[ielem][variable (rho=0,u,v,w,p=4)]
    double** residu_d;  //dissipative residuals residu[ielem][variable (rho=0,u,v,w,p=4)]
    double** residu_d_hyb;  //dissipative residuals residu[ielem][variable (rho=0,u,v,w,p=4)]
    double** W_0;       //copy of initial state for RK
    double inf_speed;   // speed norm at infinity
    double RKM_coef[2][5] = // RKM 5 coeffs [order-1][step]
    {{0.0533,0.1263,0.2375,0.4414,1.0},{0.0695,0.1602,0.2898,0.5060,1.0}};
    double RKH_coef[2][5] = // RKH [0=alpha,1=beta][step]
    {{0.2742,0.2069,0.5020,0.5142,1.0},{1.0,0.0,0.56,0.0,0.44}};

    //Private methods
    void initialisation();      // Initialise field to infinity
    void updateBound();         // Update boundary conditions
    void exchangePrimitive();   // MPI exchange primitive values
    void exchangeGradiants();   // MPI exchange gradiant values
    void timeStepEul();         // Euler explicit time integration
    void timeStepRkM();         // Runge-Kutta Multistage time integration
    void timeStepRkH();         // Runge-Kutta Hybride time integration
    void computeFluxO1();       // Calcul des flux (Roe) ordre 1
    void computeFluxO2();       // Calcul des flux (Roe) ordre 2
    //void computeGradiants();  // [May not be needed, depending on choosen implem]
    void computeResidu();       // Calcul des résidu
    double checkConvergence();  // Somme des résidus
};