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
    for (int i=0;i<nelem;++i) {
        delete[] residu_c[i];
        delete[] residu_d[i];
        delete[] residu_d_hyb[i];
        delete[] W_0[i];
    }
    delete[] flux_c,flux_d,residu_c,residu_d,W_0,residu_d_hyb;
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

    for (int i=0;i<ncell;++i) {
        // initialise conservatives
        rho[i] = 1.0;
        u[i] = inf_speed*cos(AoA);
        v[i] = inf_speed*sin(AoA);
        w[i] = 0;   // TBD?
        p[i] = 1.0;
    }

    W_0 = new double*[nelem];
    residu_c = new double*[nelem];
    residu_d = new double*[nelem];
    residu_d_hyb = new double*[nelem];
     for (int i=0;i<nelem;++i) {
        W_0[i] = new double [5];
        residu_c[i] = new double [5];
        residu_d[i] = new double [5];
        residu_d_hyb[i] = new double [5];
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
    double dTi_sans_V;         //dTi local time step
    for (int ielem=0;ielem<nelem;++ielem) {
        c = sqrt(1.4*p[ielem]/rho[ielem]);
        sumLambda = (fabs(u[ielem])+c)*elem2spectral[ielem][0]+(fabs(v[ielem])+c)*elem2spectral[ielem][1]+(fabs(w[ielem])+c)*elem2spectral[ielem][2];
        dTi_sans_V = cfl/sumLambda;
        //dTi = dTi_sans_V*elem2vol[ielem];
        rho[ielem] -= (residu_c[ielem][0]+residu_d[ielem][0])*dTi_sans_V;        // CAREFULL WITH SIGN OF DISSIPATIVE FLUX
        dTi_sans_V = dTi_sans_V/rho[ielem];
        u[ielem] -= (residu_c[ielem][1]-residu_d[ielem][1])*dTi_sans_V;
        v[ielem] -= (residu_c[ielem][2]-residu_d[ielem][2])*dTi_sans_V;
        w[ielem] -= (residu_c[ielem][3]-residu_d[ielem][3])*dTi_sans_V;
        p[ielem] -= (residu_c[ielem][4]-residu_d[ielem][4])*dTi_sans_V;
    }
}

// Runge-Kutta Multistage time integration
void solver::timeStepRkM() {
    int RK_step = 5;
    double c;           //local sound speed
    double sumLambda;   // see blasek, sum of spectral radiuses
    double dTi_sans_V;         //dTi local time step
    // update W_0 with initial results
    for (int ielem=0;ielem<nelem;++ielem) {
        W_0[ielem][0] = rho[ielem];
        W_0[ielem][1] = u[ielem];
        W_0[ielem][2] = v[ielem];
        W_0[ielem][3] = w[ielem];
        W_0[ielem][4] = p[ielem]; 
    }

    for (int k=0;k<RK_step;++k) {
        for (int ielem=0;ielem<nelem;++ielem) {
            c = sqrt(1.4*p[ielem]/rho[ielem]);
            sumLambda = (fabs(u[ielem])+c)*elem2spectral[ielem][0]+(fabs(v[ielem])+c)*elem2spectral[ielem][1]+(fabs(w[ielem])+c)*elem2spectral[ielem][2];
            dTi_sans_V = cfl/sumLambda;
            rho[ielem] = W_0[ielem][0] - RKM_coef[Order-1][k]*dTi_sans_V*(residu_c[ielem][0]-residu_d[ielem][0]);
            dTi_sans_V = dTi_sans_V/rho[ielem];
            u[ielem] = W_0[ielem][1] - RKM_coef[Order-1][k]*dTi_sans_V*(residu_c[ielem][1]-residu_d[ielem][1]);
            v[ielem] = W_0[ielem][2] - RKM_coef[Order-1][k]*dTi_sans_V*(residu_c[ielem][2]-residu_d[ielem][2]);
            w[ielem] = W_0[ielem][3] - RKM_coef[Order-1][k]*dTi_sans_V*(residu_c[ielem][3]-residu_d[ielem][3]);
            p[ielem] = W_0[ielem][4] - RKM_coef[Order-1][k]*dTi_sans_V*(residu_c[ielem][4]-residu_d[ielem][4]);
        }
        // update residu
        if (k!=RK_step) {
            timeStepEul();
            exchangePrimitive();
            updateBound();
            if (Order==1){computeFluxO1();}else {computeFluxO2();}
            computeResidu();
        }
    }
}

// Runge-Kutta Hybrid time integration
void solver::timeStepRkH() {
    int RK_step = 5;
    double c;           //local sound speed
    double sumLambda;   // see blasek, sum of spectral radiuses
    double dTi_sans_V;         //dTi local time step
    // update W_0 with initial results
    for (int ielem=0;ielem<nelem;++ielem) {
        W_0[ielem][0] = rho[ielem];
        W_0[ielem][1] = u[ielem];
        W_0[ielem][2] = v[ielem];
        W_0[ielem][3] = w[ielem];
        W_0[ielem][4] = p[ielem]; 
    }

    for (int k=0;k<RK_step;++k) {
        // update residu_d_hyb
        if (k==0) {
            for (int ielem=0;ielem<nelem;++ielem) {
                for (int jres=0;jres<5;++jres) {
                    residu_d_hyb[ielem][jres] = residu_d[ielem][jres];
                }
            }
        }
        else if((k==2)||(k==4)) {
            for (int ielem=0;ielem<nelem;++ielem) {
                for (int jres=0;jres<5;++jres) {
                    residu_d_hyb[ielem][jres] = RKH_coef[1][k]*residu_d[ielem][jres]+(1-RKH_coef[1][k])*residu_d_hyb[ielem][jres];
                }
            }
        }
        // do step
        for (int ielem=0;ielem<nelem;++ielem) {
            c = sqrt(1.4*p[ielem]/rho[ielem]);
            sumLambda = (fabs(u[ielem])+c)*elem2spectral[ielem][0]+(fabs(v[ielem])+c)*elem2spectral[ielem][1]+(fabs(w[ielem])+c)*elem2spectral[ielem][2];
            dTi_sans_V = cfl/sumLambda;
            rho[ielem] = W_0[ielem][0] - RKH_coef[0][k]*dTi_sans_V*(residu_c[ielem][0]-residu_d[ielem][0]);
            dTi_sans_V = dTi_sans_V/rho[ielem];
            u[ielem] = W_0[ielem][1] - RKH_coef[0][k]*dTi_sans_V*(residu_c[ielem][1]-residu_d_hyb[ielem][1]);
            v[ielem] = W_0[ielem][2] - RKH_coef[0][k]*dTi_sans_V*(residu_c[ielem][2]-residu_d_hyb[ielem][2]);
            w[ielem] = W_0[ielem][3] - RKH_coef[0][k]*dTi_sans_V*(residu_c[ielem][3]-residu_d_hyb[ielem][3]);
            p[ielem] = W_0[ielem][4] - RKH_coef[0][k]*dTi_sans_V*(residu_c[ielem][4]-residu_d_hyb[ielem][4]);
        }
        // update residu
        if (k!=RK_step) {
            timeStepEul();
            exchangePrimitive();
            updateBound();
            if (Order==1){computeFluxO1();}else {computeFluxO2();}  // Diffusive fluxes do not need to be computed at every step!
            computeResidu();
        }
    }
}

// Roe fluxes, order 1 [REMEMBER TO SPLIT CONVECTIVE AND DIFFUSIVE FLUXES]
void solver::computeFluxO1() {
	for (int iface = 0; iface < nface; iface++) {    
	

			int ielemL, ielemR;
			double dp, du, dv, dV, drho,dw;
			double rhoL, rhoR, VL, VR, UL, UR,WR,WL, uL, uR, vL, vR,wR,wL, pL, pR, cL, cR, HL, HR, nx, ny, rhobar, ubar, vbar,wbar, hbar, cbar, Vbar, qbar, SR1, SR2, SR3, delta;
			double rhoAvg, uAvg, vAvg,wAvg, pAvg, Vavg, Havg, Uavg;
			double imassFlux, imomentumFlux[2], ienergyFlux, F1mass, F1mom1, F1mom2, F1mom3, F1energy, F234mass, F234mom1, F234mom2,F234mom3, F234energy, F5mass, F5mom1, F5mom2, F5mom3, F5energy;
			double AWW1, AWW2, AWW3, AWW4,AWW5, Fcmass, Fcmom1, Fcmom2,Fcmom3, Fcenergy;
			double c1, c2, c3, c4, c5;

			ielemL = face2elem[iface][0] - 1;		// numero de l'element Gauche 
			ielemR = face2elem[iface][1] - 1;		// numero de l'element droit

			rhoL = rho[ielemL];
			rhoR = rho[ielemR];
			drho = rhoR - rhoL;

			//Calcul des vitesses à droite et à gauche de chaque élément
			uL = Velocity[ielemL][0];    //A changer pour la bonne variable
			vL = Velocity[ielemL][1];
			wL = Velocity[ielemL][2];
			uR = Velocity[ielemR][0];
			vR = Velocity[ielemR][1];
			wR = Velocity[ielemR][1];
			du = uR - uL;
			dv = vR - vL;
			dw = wR - wL; 

			//Calcul des normales
			nx = face2norm[iface][0];  
			ny = face2norm[iface][1];
			nz = face2norm[iface][2];
				

			UL = sqrt(uL * uL + vL * vL + wL*wL); //ajouter w??
			UR = sqrt(uR * uR + vR * vR + wR*wR);


			VL = nx * uL + ny * vL + nz*wL;     //ajouter w?
			VR = nx * uR + ny * vR + nz*wR;

			dV = VR - VL;
				


			pL = p[ielemL];
			pR = p[ielemR];
			dp = pR - pL;


			//Calcul des moyennes des variables 
			rhoAvg = 0.5 * (rhoL + rhoR);
			uAvg = 0.5 * (uL + uR);
			vAvg = 0.5 * (vL + vR);
			wAvg = 0.5 * (wL + wR);
			pAvg = 0.5 * (pL + pR);
			Vavg = uAvg * nx + vAvg * ny + wAvg*nz;  
			Uavg = sqrt(uAvg * uAvg + vAvg * vAvg + wAvg * wAvg);  //ajouter w?

			HL = 0.5 * UL * UL + pL / rhoL / (gamma - 1) + pL / rhoL;
			HR = 0.5 * UR * UR + pR / rhoR / (gamma - 1) + pR / rhoR;
			Havg = 0.5 * (HL + HR);

			//Calcul constantes pour simplifier la compilation 
			c1=sqrt(rhoL); 
			c2=sqrt(rhoR);
			c3= 1/(c1+c2);
			c4= 1/(2 * cbar * cbar); 
			c5= 1/(cbar * cbar);
			
			
			//Calcul des variables bar du schéma ROE
			rhobar = sqrt(rhoL * rhoR);
			ubar = (uL * c1 + uR * c2) * c3;
			vbar = (vL * c1 + vR * c2) * c3; 
			wbar = (wL * c1 + wR * c2) * c3;
			hbar = (HL * c1 + HR * c2) * c3;
			qbar = (ubar * ubar) + (vbar * vbar)+(wbar*wbar); // c'est le qbar au carré
			cbar = sqrt((gamma - 1) * (hbar - qbar * 0.5)); 
			Vbar = ubar * nx + vbar * ny+ wbar * nz;

			cL = sqrt(gamma * pL / rhoL);
			cR = sqrt(gamma * pR / rhoR);


			// Harten’s entropy correction
			SR1 = abs(Vbar - cbar);
			SR2 = abs(Vbar);
			SR3 = abs(Vbar + cbar);

			delta = 0.1 * (cL + cR) / 2;    // Potentielle erreur pour vitesse dus on locale

			if (SR1 < delta && SR1>-delta) {         
				SR1 = (SR1 + delta * delta) / (2 * delta);
			}
			if (SR2 < delta && SR2>-delta) {
				SR2 = (SR2 + delta * delta) / (2 * delta);
			}
			if (SR3 < delta && SR3>-delta) {
				SR3 = (SR3 + delta * delta) / (2 * delta);
			}
				

			//Calcul des différents termes du flux dissipatif
			F1mass = SR1 * ((dp - rhobar * cbar * dV) *c4) * 1;
			F1mom1 = (SR1 * ((dp - rhobar * cbar * dV)*c4)) * (ubar - (cbar * nx));
			F1mom2 = (SR1 * ((dp - rhobar * cbar * dV) *c4)) * (vbar - (cbar * ny));
			F1mom3 = (SR1 * ((dp - rhobar * cbar * dV) *c4)) * (wbar - (cbar * nz));
			F1energy = (SR1 * ((dp - rhobar * cbar * dV) *c4)) * (hbar - (cbar * Vbar));


			F234mass = SR2 * ((drho - (dp *c5)) * 1 + rhobar * 0); 
			F234mom1 = SR2 * ((drho - (dp *c5)) * ubar + rhobar * (du - dV * nx));
			F234mom2 = SR2 * ((drho - (dp *c5)) * vbar + rhobar * (dv - dV * ny));
			F234mom3 = SR2 * ((drho - (dp *c5)) * wbar + rhobar * (dw - dV * nz));
			F234energy = SR2 * ((drho - (dp*c5)) * (qbar / 2) + rhobar * (ubar * du + vbar * dv +wbar * dw - Vbar * dV)); 


			F5mass = SR3 * ((dp + rhobar * cbar * dV) *c4) * 1;
			F5mom1 = SR3 * ((dp + rhobar * cbar * dV) *c4) * (ubar + cbar * nx);
			F5mom2 = SR3 * ((dp + rhobar * cbar * dV) *c4) * (vbar + cbar * ny);
			F5mom3 = SR3 * ((dp + rhobar * cbar * dV) *c4) * (wbar + cbar * nz);
			F5energy = SR3 * ((dp + rhobar * cbar * dV)*c4) * (hbar + cbar * Vbar);


			//Sommes des termes pour avoir le flux dissipatif
			AWW1 = 0.5 * (F1mass + F234mass + F5mass);
			AWW2 = 0.5 * (F1mom1 + F234mom1 + F5mom1);
			AWW3 = 0.5 * (F1mom2 + F234mom2 + F5mom2);
			AWW4 = 0.5 * (F1mom3 + F234mom3 + F5mom3);
			AWW5 = 0.5 * (F1energy + F234energy + F5energy);

			//Calcul du flux conservatif
			Fcmass = rhoAvg * Vavg;
			Fcmom1 = rhoAvg * Vavg * uAvg + pAvg * nx;
			Fcmom2 = rhoAvg * Vavg * vAvg + pAvg * ny;
			Fcmom3 = rhoAvg * Vavg * wAvg + pAvg * nz;
			Fcenergy = rhoAvg * Vavg * Havg;
			
			// Flux dans les bonnes variables
			flux_c[iface][0] = Fcmass;   //rho
			flux_c[iface][1] = Fcmom1;   //u
			flux_c[iface][2] = Fcmom2;   //v
			flux_c[iface][3] = Fcmom3;   //w
			flux_c[iface][4] = Fcenergy;   //p
			

			
			flux_d[iface][0] = AWW1;   //rho
			flux_d[iface][1] = AWW2;   //u
			flux_d[iface][2] = AWW3;   //v
			flux_d[iface][3] = AWW4;   //w
			flux_d[iface][4] = AWW5;   //p
			
		}
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
        residuSum += pow((residu_c[ielem][0]-residu_d[ielem][0])/elem2vol[ielem],2);        // CAREFULL WITH SIGN OF DISSIPATIVE FLUX
    }
    return sqrt(residuSum);
}