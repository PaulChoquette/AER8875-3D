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
    double** flux_c;  //convective flux  flux[iface][variable (rho=0,u,v,w,p=4)]
    double** flux_d;  //dissipative flux flux[iface][variable (rho=0,u,v,w,p=4)]
    double** residu_c;//convective residuals residu[ielem][variable (rho=0,u,v,w,p=4)]
    double** residu_d;//dissipative residuals residu[ielem][variable (rho=0,u,v,w,p=4)]
    double inf_speed;   // speed norm at infinity

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