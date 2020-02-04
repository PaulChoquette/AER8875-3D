#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <solver.h>
#include <Comm.h>
#include <math.h>
#include <supposed_metric_n_connec.h>

// constructor
solver::solver() {
    inf_speed = mach*sqrt(1.4);
    initialisation();
}

// destructor
solver::~solver() {
// DELETE ALL CREATED ARRAYS
    delete[] rho,u,v,w,p; 
    for (int i=0;i<nface;++i) {
        delete[] flux_c[i];
        delete[] flux_d[i];
    }
    for (int i=0;i<ncell;++i) {
        delete[] residu_c[i];
        delete[] residu_d[i];
    }
    delete[] flux_c,flux_d,residu_c,residu_d;
}

// call all other methods in order while not converged
// currently implemented for euler explicit & order 1 scheme
void solver::compute() {
    double Residu = 2346874653674854; 
    double Old_Residu = 590456789087654678;
    double Residu_initial;
    int iteration = 0;
    Comm Cluster(Registre_file);    // Initialise MPI

    updateBound();
    while (iteration<iterMax) {
        ++iteration;
        computeFluxO1();
        computeResidu();
        if (iteration==1) {
            Residu_initial = checkConvergence();
        }
        else {
            Residu = checkConvergence();
            printf("Iteration : %d, \tREL R :  %e, \tABS R :  %e\n",iteration,Residu/Residu_initial,Residu);
        }
        
        if ((abs(Residu/Residu_initial)<convergeCrit)&&(abs(Old_Residu/Residu_initial)<convergeCrit)) {
            break;
        }
        Old_Residu = Residu;
        timeStepEul();
        exchangePrimitive();
        updateBound();
    }
}

// set every element to infinity
void solver::initialisation() {
// initialise arrays
    rho = new double [ncell]; 
    u = new double [ncell]; 
    v = new double [ncell]; 
    w = new double [ncell]; 
    p = new double [ncell]; 
    flux_c = new double*[nface];
    flux_d = new double*[nface];
    for (int i=0;i<nface;++i) {
        flux_c[i] = new double [5];
        flux_d[i] = new double [5];
    }
    residu_c = new double*[ncell];
    residu_d = new double*[ncell];
    for (int i=0;i<ncell;++i) {
        residu_c[i] = new double [5];
        residu_d[i] = new double [5];
        // initialise conservatives
        rho[i] = 1.0;
        u[i] = inf_speed*cos(AoA);
        v[i] = inf_speed*sin(AoA);
        w[i] = 0;   // TBD?
        p[i] = 1.0;
    }
}

// update boundary conditions
void solver::updateBound() {

}

// exchange primitive values between zones
void solver::exchangePrimitive() {

}

// exchange gradiant values between zones [Order 2 only]
void solver::exchangeGradiants() {

}

// euler explicit time integration
void solver::timeStepEul() {
    double c;           //local sound speed
    double sumLambda;   // see blasek, sum of spectral radiuses
    double dTi,dTi_sans_V;         //dTi local time step
    for (int ielem=0;ielem<nelem;++ielem) {
        c = sqrt(1.4*p[ielem]/rho[ielem]);
        sumLambda = (fabs(u[ielem])+c)*elem2spectral[ielem][0]+(fabs(v[ielem])+c)*elem2spectral[ielem][1]+(fabs(w[ielem])+c)*elem2spectral[ielem][2];
        dTi_sans_V = cfl/sumLambda;
        //dTi = dTi_sans_V*elem2vol[ielem];
        rho[ielem] -= (residu_c[ielem][0]+residu_d[ielem][0])*dTi_sans_V;        // CAREFULL WITH SIGN OF DISSIPATIVE FLUX
        u[ielem] -= (residu_c[ielem][1]+residu_d[ielem][1])*dTi_sans_V;
        v[ielem] -= (residu_c[ielem][2]+residu_d[ielem][2])*dTi_sans_V;
        w[ielem] -= (residu_c[ielem][3]+residu_d[ielem][3])*dTi_sans_V;
        p[ielem] -= (residu_c[ielem][4]+residu_d[ielem][4])*dTi_sans_V;
    }
}

// Runge-Kutta Multistage time integration
void solver::timeStepRkM() {

}

// Runge-Kutta Hybrid time integration
void solver::timeStepRkH() {

}

// Roe fluxes, order 1 [REMEMBER TO SPLIT CONVECTIVE AND DIFFUSIVE FLUXES]
void solver::computeFluxO1() {

}

// Roe fluxes, order 2 [REMEMBER TO SPLIT CONVECTIVE AND DIFFUSIVE FLUXES]
void solver::computeFluxO2() {

}

// Roe fluxes, order 2
void solver::computeResidu() {

}

// Roe fluxes, order 2
double solver::checkConvergence() {
    double residuSum = 0;
    for (int ielem=0;ielem<nelem;++ielem) { //Sum of residu in rho
        residuSum += pow((residu_c[ielem][0]+residu_d[ielem][0])/elem2vol[ielem],2);        // CAREFULL WITH SIGN OF DISSIPATIVE FLUX
    }
    return sqrt(residuSum);
}