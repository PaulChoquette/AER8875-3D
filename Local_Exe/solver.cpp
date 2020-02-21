#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <solver.h>
#include <math.h>

// constructor
solver::solver() {
    inf_speed = mach*sqrt(1.4);
    inf_speed_x = inf_speed*cos(AoA);
    inf_speed_y = inf_speed*sin(AoA);
    inf_speed_z = 0;
    Initialisation();
}

// destructor
solver::~solver() {
// DELETE ALL CREATED ARRAYS
    delete[] rho;
    delete[] u;
    delete[] v;
    delete[] w;
    delete[] p; 
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
    delete[] flux_c;
    delete[] flux_d;
    delete[] residu_c;
    delete[] residu_d;
    delete[] W_0;
    delete[] residu_d_hyb;


    // Delete MPI-related variables
    for (int izone=0;izone<World.ntgt;++izone) {
        delete[] primitivesSendBuffer[izone];
    }
    delete[] primitivesSendBuffer;
}

// call all other methods in order while not converged
// currently implemented for euler explicit & order 1 scheme
void solver::Compute() {
    double Residu = 2346874653674854; 
    double Old_Residu = 590456789087654678;
    double Residu_initial;
    int iteration = 0;
    InitMPI();
    UpdateBound();
    while (iteration<iterMax) {
        ++iteration;
        ComputeFluxO1();
        ComputeResidu();
        if (iteration==1) {
            Residu_initial = CheckConvergence();
        }
        else {
            Residu = CheckConvergence();
            printf("Iteration : %d, \tREL R :  %e, \tABS R :  %e\n",iteration,Residu/Residu_initial,Residu);
        }
        
        if ((abs(Residu/Residu_initial)<convergeCrit)&&(abs(Old_Residu/Residu_initial)<convergeCrit)) {
            break;
        }
        Old_Residu = Residu;
        TimeStepEul();
        ExchangePrimitive();
        UpdateBound();
    }
}

// set every element to infinity
void solver::Initialisation() {
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
void solver::UpdateBound() {
 int ig,ir;// ig = ghost index, ir = real index
    double c0,l_r0c0,r0c0,udotn; //local sound speed
    
    // See if possibility to reduce ifs...
    for (int iface=-1;iface<-1;++iface) {   // TO EDIT
        ir = face2elem[iface][0];
        ig = face2elem[iface][1];           // ASK CONFIRMATION

        if (1){ //If slipwall
            // Slipwall
            udotn = 2.0*(face2norm[iface][0]*u[ir]+face2norm[iface][1]*v[ir]+face2norm[iface][2]*w[ir]);
            p[ig] = p[ir];  //There are higher accruacy methods
            rho[ig] = rho[ir];
            u[ig] = u[ir]-udotn*face2norm[0][iface];    // p.273, eq 8.10
            v[ig] = v[ir]-udotn*face2norm[1][iface];
            w[ig] = w[ir]-udotn*face2norm[2][iface];
        }
        else if (1) {//If symmetry
            // Symmetry
            p[ig] = p[ir];
            rho[ig] = rho[ir];
            u[ig] = u[ir];  //p.285 8.38
            v[ig] = v[ir];
            w[ig] = w[ir];
        }
        else {
            c0 = sqrt(1.4*(p[ir]/rho[ir]));
            udotn = face2norm[iface][0]*u[ir]+face2norm[iface][1]*v[ir]+face2norm[iface][2]*w[ir];
            if ((pow(u[ir],2)+pow(v[ir],2)+pow(w[ir],2))>=pow(c0,2)) {  //if supersonic
                if (udotn>0) {    //inlet
                // Supersonic Inflow
                    p[ig] = 1.0;
                    rho[ig] = 1.0;
                    u[ig] = inf_speed_x;
                    v[ig] = inf_speed_y;
                    w[ig] = inf_speed_z;
                }
                else {
                    p[ig] = p[ir];
                    rho[ig] = rho[ir];
                    u[ig] = u[ir];
                    v[ig] = v[ir];
                    w[ig] = w[ir];
                }
            }
            else {      // If subsonic
                l_r0c0 = 1/(c0*rho[ir]);
                if (udotn>0) {    //inlet
                    // Subsonic Inflow
                    p[ig] = 0.5*(1.0+p[ir]-c0*rho[ir]*(face2norm[iface][0]*(inf_speed_x-u[ir])+face2norm[iface][1]*(inf_speed_y-v[ir])+face2norm[iface][2]*(inf_speed_z-w[ir])));
                    rho[ig] = rho[ir]+(p[ig]-1.0)/pow(c0,2);
                    u[ig] = inf_speed_x-face2norm[iface][0]*(1.0-p[ig])*l_r0c0;
                    v[ig] = inf_speed_y-face2norm[iface][1]*(1.0-p[ig])*l_r0c0;
                    w[ig] = inf_speed_z-face2norm[iface][2]*(1.0-p[ig])*l_r0c0;
                }
                else {
                    // Subsonic Outflow
                    p[ig] = 1.0;
                    rho[ig] = rho[ir]+(1.0-p[ir])/pow(c0,2);
                    u[ig] = u[ir]+face2norm[iface][0]*(p[ir]-1.0)*l_r0c0;
                    v[ig] = v[ir]+face2norm[iface][1]*(p[ir]-1.0)*l_r0c0;
                    w[ig] = w[ir]+face2norm[iface][2]*(p[ir]-1.0)*l_r0c0;
                }
            }

        }
    }
}

// Initialise MPI
void solver::InitMPI() {
    int* zone2nbelem;   //number of element / zone boundary
    int ntgt;           //Number of zones to talk to


    World.Init(Registre_file);    // Initialise MPI

    //Build zone2nbelem; number of element / zone boundary
    zone2nbelem = new int[ntgt];
    for (int izone=0;izone<ntgt;++izone) {
        zone2nbelem[izone] = 0;// TO EDIT; = TO NUMBER OF BOUNDARY ELEMENTS BETWEEN THE 2 ZONES
    }
    World.InitBuffer(zone2nbelem);

    // Create transmition buffer
    primitivesSendBuffer = new double*[World.ntgt];
    for (int izone=0;izone<World.ntgt;++izone) {
        primitivesSendBuffer[izone] = new double[zone2nbelem[izone]];
    }

    // Handshake order in which boundary elements will be sent between zones
    //World.ExchangeCellOrder();    // TO EDIT

    delete[] zone2nbelem;
}

// exchange primitive values between zones
void solver::ExchangePrimitive() {
    int ibelemIndx; //Local boundary element index
    // Populate Tx Buffer
    for (int izone=0;izone<World.ntgt;++izone) {
        for (int ibelem=0;ibelem<World.zone2nbelem[izone];++ibelem) {
            ibelemIndx = 0;// TO EDIT 
            primitivesSendBuffer[izone][ibelem] = rho[ibelemIndx];
            primitivesSendBuffer[izone][ibelem+World.zone2nbelem[izone]] = u[ibelemIndx];
            primitivesSendBuffer[izone][ibelem+World.zone2nbelem[izone]*2] = v[ibelemIndx];
            primitivesSendBuffer[izone][ibelem+World.zone2nbelem[izone]*3] = w[ibelemIndx];
            primitivesSendBuffer[izone][ibelem+World.zone2nbelem[izone]*4] = p[ibelemIndx];
        }
    }
    //Send & Store
    World.ExchangePrimitives(primitivesSendBuffer);
    World.ReclassPrimitives(&rho,&u,&v,&w,&p);
}

// exchange gradiant values between zones [Order 2 only]
void solver::ExchangeGradiants() {

}

// euler explicit time integration
void solver::TimeStepEul() {
    double c, eTempo,invrho;           //local sound speed
    double sumLambda;   // see blasek, sum of spectral radiuses
    double dTi_sans_V;         //dTi local time step
    for (int ielem=0;ielem<nelem;++ielem) {
        c = sqrt(1.4*p[ielem]/rho[ielem]);
        sumLambda = (fabs(u[ielem])+c)*elem2spectral[ielem][0]+(fabs(v[ielem])+c)*elem2spectral[ielem][1]+(fabs(w[ielem])+c)*elem2spectral[ielem][2];
        dTi_sans_V = cfl/sumLambda;
		invrho=1/rho[ielem];
        rho[ielem] -= ((residu_c[ielem][0]+residu_d[ielem][0])*dTi_sans_V);        // CAREFULL WITH SIGN OF DISSIPATIVE FLUX
        u[ielem] -= ((residu_c[ielem][1]-residu_d[ielem][1])*dTi_sans_V)*invrho;
        v[ielem] -= ((residu_c[ielem][2]-residu_d[ielem][2])*dTi_sans_V)*invrho;
        w[ielem] -= ((residu_c[ielem][3]-residu_d[ielem][3])*dTi_sans_V)*invrho;
		
		eTempo = P2E(p[ielem])-(residu_c[ielem][4]-residu_d[ielem][4])*dTi_sans_V;
        p[ielem] = E2P(eTempo);
    }
}

// Runge-Kutta Multistage time integration
void solver::TimeStepRkM() {
    int OIndx = Order-1;
    double c,invrho,eTempo;           //local sound speed
    double sumLambda;   // see blasek, sum of spectral radiuses
    double dTi_sans_V;         //dTi local time step
    
    // update W_0 with initial results
    for (int ielem=0;ielem<nelem;++ielem) {
        W_0[ielem][0] = rho[ielem];
        W_0[ielem][1] = rho[ielem]*u[ielem];
        W_0[ielem][2] = rho[ielem]*v[ielem];
        W_0[ielem][3] = rho[ielem]*w[ielem];
        W_0[ielem][4] = P2E(p[ielem]);         //To edit
    }

    for (int k=0;k<RK_step;++k) {
        for (int ielem=0;ielem<nelem;++ielem) {
            c = sqrt(1.4*p[ielem]/rho[ielem]);
            sumLambda = (fabs(u[ielem])+c)*elem2spectral[ielem][0]+(fabs(v[ielem])+c)*elem2spectral[ielem][1]+(fabs(w[ielem])+c)*elem2spectral[ielem][2];
            dTi_sans_V = cfl/sumLambda;
            rho[ielem] = W_0[ielem][0] - RKM_coef[RK_step][OIndx][k]*dTi_sans_V*(residu_c[ielem][0]-residu_d[ielem][0]);
			invrho=1/rho[ielem];
            u[ielem] = (W_0[ielem][1] - RKM_coef[RK_step][OIndx][k]*dTi_sans_V*(residu_c[ielem][1]-residu_d[ielem][1]))*invrho;
            v[ielem] = (W_0[ielem][2] - RKM_coef[RK_step][OIndx][k]*dTi_sans_V*(residu_c[ielem][2]-residu_d[ielem][2]))*invrho;
            w[ielem] = (W_0[ielem][3] - RKM_coef[RK_step][OIndx][k]*dTi_sans_V*(residu_c[ielem][3]-residu_d[ielem][3]))*invrho;
			
			eTempo= (W_0[ielem][4] - RKM_coef[RK_step][OIndx][k]*dTi_sans_V*(residu_c[ielem][4]-residu_d[ielem][4]))*invrho;
            p[ielem] = E2P(eTempo);    // to edit
        }
        // update residu
        if (k!=RK_step) {
            ExchangePrimitive();
            UpdateBound();
            if (Order==1){ComputeFluxO1();}else {ComputeFluxO2();}
            ComputeResidu();
        }
    }
}



// Runge-Kutta Hybrid time integration
void solver::TimeStepRkH() {
    double c, eTempo, invrho;           //local sound speed
    double sumLambda;   // see blasek, sum of spectral radiuses
    double dTi_sans_V;  //dTi local time step
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
                    residu_d_hyb[ielem][jres] = RKH_coef[RK_step][1][k]*residu_d[ielem][jres]+(1-RKH_coef[RK_step][1][k])*residu_d_hyb[ielem][jres];
                }
            }
        }
        // do step
        for (int ielem=0;ielem<nelem;++ielem) {
            c = sqrt(1.4*p[ielem]/rho[ielem]);
            sumLambda = (fabs(u[ielem])+c)*elem2spectral[ielem][0]+(fabs(v[ielem])+c)*elem2spectral[ielem][1]+(fabs(w[ielem])+c)*elem2spectral[ielem][2];
            dTi_sans_V = cfl/sumLambda;
            rho[ielem] = W_0[ielem][0] - RKM_coef[RK_step][OIndx][k]*dTi_sans_V*(residu_c[ielem][0]-residu_d_hyb[ielem][0]);
			invrho=1/rho[ielem];
            u[ielem] = (W_0[ielem][1] - RKM_coef[RK_step][OIndx][k]*dTi_sans_V*(residu_c[ielem][1]-residu_d_hyb[ielem][1]))*invrho;
            v[ielem] = (W_0[ielem][2] - RKM_coef[RK_step][OIndx][k]*dTi_sans_V*(residu_c[ielem][2]-residu_d_hyb[ielem][2]))*invrho;
            w[ielem] = (W_0[ielem][3] - RKM_coef[RK_step][OIndx][k]*dTi_sans_V*(residu_c[ielem][3]-residu_d_hyb[ielem][3]))*invrho;
			
			eTempo= (W_0[ielem][4] - RKM_coef[RK_step][OIndx][k]*dTi_sans_V*(residu_c[ielem][4]-residu_d_hyb[ielem][4]))*invrho;
            p[ielem] = E2P(eTempo);    // to edit
        }
        // update residu
        if (k!=RK_step) {
            ExchangePrimitive();
            UpdateBound();
			if (k==1 || k==3) {
				if (Order==1){ComputeFluxO1();}else {ComputeFluxO2();}  // Diffusive fluxes do not need to be computed at every step!
			}
			else {
				if (Order==1){ComputeFluxO1Conv();}else {ComputeFluxO2Conv();}  // Diffusive fluxes do not need to be computed at every step!
			}
            ComputeResidu();
        }
    }
}
void solver::ComputeFluxO1Conv() {
	for (int iface = 0; iface < nface; iface++) {    
			int ielemL, ielemR;
			double dp, du, dv, dV, drho,dw;
			double rhoL, rhoR, VL, VR, UL, UR,WR,WL, uL, uR, vL, vR,wR,wL, pL, pR, cL, cR, HL, HR, nx, ny, nz;
			double rhoAvg, uAvg, vAvg,wAvg, pAvg, Vavg, Havg, Uavg;
			double Fcmass, Fcmom1, Fcmom2,Fcmom3, Fcenergy;

			ielemL = face2elem[iface][0] - 1;		// numero de l'element Gauche 
			ielemR = face2elem[iface][1] - 1;		// numero de l'element droit

			rhoL = rho[ielemL];
			rhoR = rho[ielemR];
			drho = rhoR - rhoL;

			//Calcul des vitesses à droite et à gauche de chaque élément
			uL = u[ielemL];    //A changer pour la bonne variable
			vL = v[ielemL];
			wL = w[ielemL];
			uR = u[ielemR];
			vR = v[ielemR];
			wR = w[ielemR];
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
			


			
		}
	}





// Roe fluxes, order 1 [REMEMBER TO SPLIT CONVECTIVE AND DIFFUSIVE FLUXES]
void solver::ComputeFluxO1() {
	for (int iface = 0; iface < nface; iface++) {    
			int ielemL, ielemR;
			double dp, du, dv, dV, drho,dw;
			double rhoL, rhoR, VL, VR, UL, UR,WR,WL, uL, uR, vL, vR,wR,wL, pL, pR, cL, cR, HL, HR, nx, ny, nz, rhobar, ubar, vbar,wbar, hbar, cbar, Vbar, qbar, SR1, SR2, SR3, delta;
			double rhoAvg, uAvg, vAvg,wAvg, pAvg, Vavg, Havg, Uavg;
			double imassFlux, imomentumFlux[2], ienergyFlux, F1mass, F1mom1, F1mom2, F1mom3, F1energy, F234mass, F234mom1, F234mom2,F234mom3, F234energy, F5mass, F5mom1, F5mom2, F5mom3, F5energy;
			double AWW1, AWW2, AWW3, AWW4,AWW5, Fcmass, Fcmom1, Fcmom2,Fcmom3, Fcenergy;
			double c1, c2, c3, c4, c5;
            double gamma=1.4;

			ielemL = face2elem[iface][0] - 1;		// numero de l'element Gauche 
			ielemR = face2elem[iface][1] - 1;		// numero de l'element droit

			rhoL = rho[ielemL];
			rhoR = rho[ielemR];
			drho = rhoR - rhoL;

			//Calcul des vitesses à droite et à gauche de chaque élément
			uL = u[ielemL];    //A changer pour la bonne variable
			vL = v[ielemL];
			wL = w[ielemL];
			uR = u[ielemR];
			vR = v[ielemR];
			wR = w[ielemR];
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
			F234energy = SR2 * ((drho - (dp*c5)) * (qbar * 0.5) + rhobar * (ubar * du + vbar * dv +wbar * dw - Vbar * dV)); 


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
void solver::ComputeFluxO2() {

}

void solver::ComputeResidu() {

}

double solver::CheckConvergence() {
    double residuSum = 0;
    for (int ielem=0;ielem<nelem;++ielem) {                                                 //Sum of residu in rho
        residuSum += pow((residu_c[ielem][0]-residu_d[ielem][0])/elem2vol[ielem],2);        // CAREFULL WITH SIGN OF DISSIPATIVE FLUX
    }
    return sqrt(residuSum);
}

double solver::P2E(double p ) {
	
	
	
	
	
	
	}
	
double solver::E2P(double e) {
	
	
	
	
	
	
	}
	