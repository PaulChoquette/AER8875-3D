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
    
    for (int ielem=0;ielem<nelem;++ielem) {
        for (int idim=0;idim<3;++idim) {
            delete[] gradient[ielem][idim];
        }
        delete[] gradient[ielem];
    }
    delete[] gradient;

    // Delete MPI-related variables
    for (int izone=0;izone<World.ntgt;++izone) {
        if (Order==2){delete[] gradientSendBuffer[izone];}
        delete[] primitivesSendBuffer[izone];
    }
    if (Order==2){delete[] gradientSendBuffer;}
    delete[] primitivesSendBuffer;
}

// call all other methods in order while not converged
void solver::Compute() {
    Initialisation();
    double Residu = 2346874653674854; 
    double Old_Residu = 590456789087654678;
    double Residu_initial;
    int iteration = 0;
    InitMPI();
    UpdateBound();
    while (iteration<iterMax) {
        ++iteration;
        if (Order==1){ComputeFluxO1();}else {ComputeFluxO2();}
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
        if (RK_step==1){TimeStepEul();}
        else if(RK_M==0) {TimeStepRkH();}else {TimeStepRkM();}
        ExchangePrimitive();
        if (Order==2){ExchangeGradiants();}
        UpdateBound();
    }
}

// set every element to infinity and initialise other stuff
void solver::Initialisation() {
// initialise arrays
    rho = new double [ncell]; 
    u = new double [ncell]; 
    v = new double [ncell]; 
    w = new double [ncell]; 
    p = new double [ncell]; 
    if (Order==2) {
        gradient = new double**[ncell];     //Needs to be ncell as boundaries [between zones] need gradients
        for (int ielem=0;ielem<ncell;++ielem) {
            gradient[ielem] = new double*[3];
            for (int idim=0;idim<3;++idim) {
                gradient[ielem][idim] = new double[5];
            }
        }
    }

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

    // Create transmition buffers
    if (Order==2){gradientSendBuffer = new double*[World.ntgt];}
    primitivesSendBuffer = new double*[World.ntgt];
    for (int izone=0;izone<World.ntgt;++izone) {
        if (Order==2){primitivesSendBuffer[izone] = new double[zone2nbelem[izone]*5];}
        gradientSendBuffer[izone] = new double[zone2nbelem[izone]*15];
    }
    // Handshake order in which boundary elements will be sent between zones
    //World.ExchangeCellOrder();    // TO EDIT
    //To use ExchangeMetrics, a buffer needs to be created first...
    if (Order==2){World.ExchangeMetrics();}
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
    int ibelemIndx; //Local boundary element index
    // Populate Tx Buffer
    for (int izone=0;izone<World.ntgt;++izone) {
        for (int ibelem=0;ibelem<World.zone2nbelem[izone];++ibelem) {
            ibelemIndx = 0;// TO EDIT 
            for (int idim;idim<3;++idim) {
                gradientSendBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5] = gradient[ibelemIndx][idim][0];
                gradientSendBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5+World.zone2nbelem[izone]] = gradient[ibelemIndx][idim][1];
                gradientSendBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5+World.zone2nbelem[izone]*2] = gradient[ibelemIndx][idim][2];
                gradientSendBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5+World.zone2nbelem[izone]*3] = gradient[ibelemIndx][idim][3];
                gradientSendBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5+World.zone2nbelem[izone]*4] = gradient[ibelemIndx][idim][4];
            }
        }
    }
    // Exchange Values
    World.ExchangeGradients(gradientSendBuffer);
    //Store received values in local arrays
    for (int izone=0;izone<World.ntgt;++izone) {
        for (int ibelem=0;ibelem<World.zone2nbelem[izone];++ibelem) {
            ibelemIndx = World.rxOrder2localOrder[izone][ibelem];
            for (int idim;idim<3;++idim) {
                gradient[ibelemIndx][idim][0] = World.gradientBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5];
                gradient[ibelemIndx][idim][1] = World.gradientBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5+World.zone2nbelem[izone]];
                gradient[ibelemIndx][idim][2] = World.gradientBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5+World.zone2nbelem[izone]*2];
                gradient[ibelemIndx][idim][3] = World.gradientBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5+World.zone2nbelem[izone]*3];
                gradient[ibelemIndx][idim][4] = World.gradientBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5+World.zone2nbelem[izone]*4];
            }
        }
    }
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
		eTempo = P2E(p[ielem],rho[ielem],u[ielem],v[ielem],w[ielem])-(residu_c[ielem][4]-residu_d[ielem][4])*dTi_sans_V;
        p[ielem] = E2P(eTempo,rho[ielem],u[ielem],v[ielem],w[ielem]);
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
        W_0[ielem][4] = P2E(p[ielem],rho[ielem],u[ielem],v[ielem],w[ielem]);
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
            p[ielem] = E2P(eTempo,rho[ielem],u[ielem],v[ielem],w[ielem]);
        }
        // update residu
        if (k!=RK_step) {
            ExchangePrimitive();
            UpdateBound();
            if (Order==1){ComputeFluxO1();}else {ComputeGrandientsNLimit();ExchangeGradiants();ComputeFluxO2();}
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
        W_0[ielem][1] = rho[ielem]*u[ielem];
        W_0[ielem][2] = rho[ielem]*v[ielem];
        W_0[ielem][3] = rho[ielem]*w[ielem];
        W_0[ielem][4] = P2E(p[ielem],rho[ielem],u[ielem],v[ielem],w[ielem]);
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
            rho[ielem] = W_0[ielem][0] - RKH_coef[RK_step][0][k]*dTi_sans_V*(residu_c[ielem][0]-residu_d_hyb[ielem][0]);
			invrho=1/rho[ielem];
            u[ielem] = (W_0[ielem][1] - RKH_coef[RK_step][0][k]*dTi_sans_V*(residu_c[ielem][1]-residu_d_hyb[ielem][1]))*invrho;
            v[ielem] = (W_0[ielem][2] - RKH_coef[RK_step][0][k]*dTi_sans_V*(residu_c[ielem][2]-residu_d_hyb[ielem][2]))*invrho;
            w[ielem] = (W_0[ielem][3] - RKH_coef[RK_step][0][k]*dTi_sans_V*(residu_c[ielem][3]-residu_d_hyb[ielem][3]))*invrho;
			
			eTempo= (W_0[ielem][4] - RKM_coef[RK_step][0][k]*dTi_sans_V*(residu_c[ielem][4]-residu_d_hyb[ielem][4]))*invrho;
            p[ielem] = E2P(eTempo,rho[ielem],u[ielem],v[ielem],w[ielem]);
        }
        // update residu
        if (k!=RK_step) {
            ExchangePrimitive();
            UpdateBound();
			if (k==1 || k==3) {
				if (Order==1){ComputeFluxO1();}else {ComputeGrandientsNLimit();ExchangeGradiants();ComputeFluxO2();}  // Diffusive fluxes do not need to be computed at every step!
                ComputeResidu();
            }
			else {
				if (Order==1){ComputeFluxO1Conv();}else {ComputeGrandientsNLimit();ExchangeGradiants();ComputeFluxO2Conv();}  // Diffusive fluxes do not need to be computed at every step!
			    ComputeResiduConv();
            }
        }
    }
}


void solver::ComputeFluxO1Conv() {
    for (int iface = 0; iface < nface; iface++) {    
        int ielemL, ielemR;
        double rhoL, rhoR, uL, uR, vL, vR, wR, wL, pL, pR;

        ielemL = face2elem[iface][0] - 1;		// numero de l'element Gauche  [-1 NEEDED!?!?!]
        ielemR = face2elem[iface][1] - 1;		// numero de l'element droit   [-1 NEEDED!?!?!]

        //Update L/R values [order 2]
        rhoL = rho[ielemL];
        rhoR = rho[ielemR];
        uL = u[ielemL];   
        vL = v[ielemL];
        wL = w[ielemL];
        uR = u[ielemR];
        vR = v[ielemR];
        wR = w[ielemR];
        pL = p[ielemL];
        pR = p[ielemR];

        UpwindFlux(iface,rhoL,uL,vL,wL,pL,rhoR,uR,vR,wR,pR);
    }  
}

// Roe fluxes, order 1 [REMEMBER TO SPLIT CONVECTIVE AND DIFFUSIVE FLUXES]
void solver::ComputeFluxO1() {
    for (int iface = 0; iface < nface; iface++) {    
        int ielemL, ielemR;
        double rhoL, rhoR, uL, uR, vL, vR, wR, wL, pL, pR;

        ielemL = face2elem[iface][0] - 1;		// numero de l'element Gauche  [-1 NEEDED!?!?!]
        ielemR = face2elem[iface][1] - 1;		// numero de l'element droit   [-1 NEEDED!?!?!]

        //Update L/R values [order 2]
        rhoL = rho[ielemL];
        rhoR = rho[ielemR];
        uL = u[ielemL];   
        vL = v[ielemL];
        wL = w[ielemL];
        uR = u[ielemR];
        vR = v[ielemR];
        wR = w[ielemR];
        pL = p[ielemL];
        pR = p[ielemR];

        UpwindFlux(iface,rhoL,uL,vL,wL,pL,rhoR,uR,vR,wR,pR);
        RoeDissipation(iface,rhoL,uL,vL,wL,pL,rhoR,uR,vR,wR,pR);   
    }   
}

// Computes grandients for order 2
void solver::ComputeGrandientsNLimit() {
    int ielem0,ielem1,locFaceIndx;
    double TempDemiSurVol,sign;       
    double delta_2,UminRho,UmaxRho,UminU,UmaxU,UminV,UmaxV,UminP,UmaxP,UminW,UmaxW; //Barth & Jespersen p:166 Blasek
    int locElem2;   //Index [0/1]
    double MachinePrecision = pow(10,-15);
    double GradU[3][5] = {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
    double psi[5] = {1,1,1,1,1};
 
    for (int ielem0=0;ielem0<nelem;++ielem0) {
        UminRho=rho[ielem0];
        UmaxRho=rho[ielem0];
        UminU=u[ielem0];
        UmaxU=u[ielem0];
        UminV=v[ielem0];
        UmaxV=v[ielem0];
        UminW=w[ielem0];
        UmaxW=w[ielem0];
        UminP=p[ielem0];
        UmaxP=p[ielem0];

        // Compute gradient
        int MaxElemPoints=0;//Edit to good termination index
        for (int ineighbor=0;ineighbor<MaxElemPoints;++ineighbor) {       //Edit to good termination index
            ielem1 = elem2elem[ielem0][ineighbor];
            if (ielem1>-1) {
                locFaceIndx = elem2face[ielem0][ineighbor];
                locElem2 = int(face2elem[locFaceIndx][0]!=ielem0);      //TBD if still is valid
                sign = ((1-locElem2)*2.0-1);
                // x
                GradU[0][0] += sign*(rho[ielem0]+rho[ielem1])*face2norm[locFaceIndx][0]*face2area[locFaceIndx];
                GradU[0][1] += sign*(u[ielem0]+u[ielem1])*face2norm[locFaceIndx][0]*face2area[locFaceIndx];
                GradU[0][2] += sign*(v[ielem0]+v[ielem1])*face2norm[locFaceIndx][0]*face2area[locFaceIndx];
                GradU[0][3] += sign*(w[ielem0]+w[ielem1])*face2norm[locFaceIndx][0]*face2area[locFaceIndx];
                GradU[0][4] += sign*(p[ielem0]+p[ielem1])*face2norm[locFaceIndx][0]*face2area[locFaceIndx];

                // y
                GradU[1][0] += sign*(rho[ielem0]+rho[ielem1])*face2norm[locFaceIndx][1]*face2area[locFaceIndx];
                GradU[1][1] += sign*(u[ielem0]+u[ielem1])*face2norm[locFaceIndx][1]*face2area[locFaceIndx];
                GradU[1][2] += sign*(v[ielem0]+v[ielem1])*face2norm[locFaceIndx][1]*face2area[locFaceIndx];
                GradU[1][3] += sign*(w[ielem0]+w[ielem1])*face2norm[locFaceIndx][1]*face2area[locFaceIndx];
                GradU[1][4] += sign*(p[ielem0]+p[ielem1])*face2norm[locFaceIndx][1]*face2area[locFaceIndx];

                // z
                GradU[2][0] += sign*(rho[ielem0]+rho[ielem1])*face2norm[locFaceIndx][2]*face2area[locFaceIndx];
                GradU[2][1] += sign*(u[ielem0]+u[ielem1])*face2norm[locFaceIndx][2]*face2area[locFaceIndx];
                GradU[2][2] += sign*(v[ielem0]+v[ielem1])*face2norm[locFaceIndx][2]*face2area[locFaceIndx];
                GradU[2][3] += sign*(w[ielem0]+w[ielem1])*face2norm[locFaceIndx][2]*face2area[locFaceIndx];
                GradU[2][4] += sign*(p[ielem0]+p[ielem1])*face2norm[locFaceIndx][2]*face2area[locFaceIndx];

                //Update Umax/Umin
                UminRho=min(UminRho,rho[ielem1]);
                UmaxRho=max(UmaxRho,rho[ielem1]);
                UminU=min(UminU,u[ielem1]);
                UmaxU=max(UmaxU,u[ielem1]);
                UminV=min(UminV,v[ielem1]);
                UmaxV=max(UmaxV,v[ielem1]);
                UminW=min(UminW,w[ielem1]);
                UmaxW=max(UmaxW,w[ielem1]);
                UminP=min(UminP,p[ielem1]);
                UmaxP=max(UmaxP,p[ielem1]);
            }
        }
        TempDemiSurVol = 0.5/elem2vol[ielem0];
        for (int i=0;i<3;++i){
            for (int j=0;j<5;++j) {
                GradU[i][j] = GradU[i][j]*TempDemiSurVol;
            }
        }
        
        // Compute Limiter  [NOTE THAT FOR EFFICIENCY AN ALTERNATIVE FORMULATION TO IFS SHOULD BE USED ONCE VALIDATED]
        for (int ineighbor=0;ineighbor<MaxElemPoints;++ineighbor) {   //For surrounding elements
            ielem1 = elem2elem[ielem0][ineighbor];
            locFaceIndx = elem2face[ielem0][ineighbor];
            locElem2 = int(face2elem[locFaceIndx][0]!=ielem0);        //TBD if still is valid, see index of cent2face below
            // rho
            delta_2 = 0.5*(GradU[0][0]*cent2face[locFaceIndx][0][locElem2]+GradU[1][0]*cent2face[locFaceIndx][1][locElem2]+GradU[2][0]*cent2face[locFaceIndx][2][locElem2]);
            if ((delta_2>0)) {
                psi[0]  = min(psi[0],fabs((UmaxRho-rho[ielem0])/(delta_2+MachinePrecision)));
            }
            else if ((delta_2<0)) {
                psi[0]  = min(psi[0],fabs((UminRho-rho[ielem0])/(delta_2-MachinePrecision)));
            }

            // U
            delta_2 = 0.5*(GradU[0][1]*cent2face[locFaceIndx][0][locElem2]+GradU[1][1]*cent2face[locFaceIndx][1][locElem2]+GradU[2][1]*cent2face[locFaceIndx][2][locElem2]);
            if ((delta_2>0)) {
                psi[1]  = min(psi[1],fabs((UmaxU-u[ielem0])/(delta_2+MachinePrecision)));
            }
            else if ((delta_2<0)) {
                psi[1]  = min(psi[1],fabs((UminU-u[ielem0])/(delta_2-MachinePrecision)));
            }

            // V
            delta_2 = 0.5*(GradU[0][2]*cent2face[locFaceIndx][0][locElem2]+GradU[1][2]*cent2face[locFaceIndx][1][locElem2]+GradU[2][2]*cent2face[locFaceIndx][2][locElem2]);
            if ((delta_2>0)) {
                psi[2]  = min(psi[2],fabs((UmaxV-v[ielem0])/(delta_2+MachinePrecision)));
            }
            else if ((delta_2<0)) {
                psi[2]  = min(psi[2],fabs((UminV-v[ielem0])/(delta_2-MachinePrecision)));
            }

            // w
            delta_2 = 0.5*(GradU[0][3]*cent2face[locFaceIndx][0][locElem2]+GradU[1][3]*cent2face[locFaceIndx][1][locElem2]+GradU[2][3]*cent2face[locFaceIndx][2][locElem2]);
            if ((delta_2>0)) {
                psi[3]  = min(psi[3],fabs((UmaxW-w[ielem0])/(delta_2+MachinePrecision)));
            }
            else if ((delta_2<0)) {
                psi[3]  = min(psi[3],fabs((UminW-w[ielem0])/(delta_2-MachinePrecision)));
            }

            // E
            delta_2 = 0.5*(GradU[0][4]*cent2face[locFaceIndx][0][locElem2]+GradU[1][4]*cent2face[locFaceIndx][1][locElem2]+GradU[2][4]*cent2face[locFaceIndx][2][locElem2]);
            if((delta_2>0)) {
                psi[4]  = min(psi[4],fabs((UmaxP-p[ielem0])/(delta_2+MachinePrecision)));
            }
            else if ((delta_2<0)) {
                psi[4]  = min(psi[4],fabs((UminP-p[ielem0])/(delta_2-MachinePrecision)));
            }
        }

        //Store & reset Gradients
        for (int idim=0;idim<3;++idim) {
            for (int ivar;ivar<5;++ivar) {
                gradient[ielem0][idim][ivar] = GradU[idim][ivar]*psi[ivar];
                GradU[idim][ivar]=0;
            }
        }
        //Reset limitors
        for (int ivar;ivar<5;++ivar) {
            psi[ivar] = 1;
        }
    }
}



// Roe fluxes, order 2 [REMEMBER TO SPLIT CONVECTIVE AND DIFFUSIVE FLUXES]
void solver::ComputeFluxO2Conv() {  
for (int iface = 0; iface < nface; iface++) {    
        int ielemL, ielemR;
        double rhoL, rhoR, uL, uR, vL, vR, wR, wL, pL, pR;

        ielemL = face2elem[iface][0] - 1;		// numero de l'element Gauche  [-1 NEEDED!?!?!]
        ielemR = face2elem[iface][1] - 1;		// numero de l'element droit   [-1 NEEDED!?!?!]

        //Update L/R values [order 2]
        rhoL = rho[ielemL]+gradient[ielemL][0][0]*cent2face[iface][0][0]+gradient[ielemL][1][0]*cent2face[iface][1][0]+gradient[ielemL][2][0]*cent2face[iface][2][0];
        uL = u[ielemL]+gradient[ielemL][0][1]*cent2face[iface][0][0]+gradient[ielemL][1][1]*cent2face[iface][1][0]+gradient[ielemL][2][1]*cent2face[iface][2][0];
        vL = v[ielemL]+gradient[ielemL][0][2]*cent2face[iface][0][0]+gradient[ielemL][1][2]*cent2face[iface][1][0]+gradient[ielemL][2][2]*cent2face[iface][2][0];
        wL = w[ielemL]+gradient[ielemL][0][3]*cent2face[iface][0][0]+gradient[ielemL][1][3]*cent2face[iface][1][0]+gradient[ielemL][2][3]*cent2face[iface][2][0];
        pL = p[ielemL]+gradient[ielemL][0][4]*cent2face[iface][0][0]+gradient[ielemL][1][4]*cent2face[iface][1][0]+gradient[ielemL][2][4]*cent2face[iface][2][0];

        rhoR = rho[ielemR]+gradient[ielemR][0][0]*cent2face[iface][0][1]+gradient[ielemR][1][0]*cent2face[iface][1][1]+gradient[ielemR][2][0]*cent2face[iface][2][1];
        uR = u[ielemR]+gradient[ielemR][0][1]*cent2face[iface][0][1]+gradient[ielemR][1][1]*cent2face[iface][1][1]+gradient[ielemR][2][1]*cent2face[iface][2][1];
        vR = v[ielemR]+gradient[ielemR][0][2]*cent2face[iface][0][1]+gradient[ielemR][1][2]*cent2face[iface][1][1]+gradient[ielemR][2][2]*cent2face[iface][2][1];
        wR = w[ielemR]+gradient[ielemR][0][3]*cent2face[iface][0][1]+gradient[ielemR][1][3]*cent2face[iface][1][1]+gradient[ielemR][2][3]*cent2face[iface][2][1];
        pR = p[ielemR]+gradient[ielemR][0][4]*cent2face[iface][0][1]+gradient[ielemR][1][4]*cent2face[iface][1][1]+gradient[ielemR][2][4]*cent2face[iface][2][1];
        UpwindFlux(iface,rhoL,uL,vL,wL,pL,rhoR,uR,vR,wR,pR);
    }   
}

void solver::ComputeFluxO2() {
    for (int iface = 0; iface < nface; iface++) {    
        int ielemL, ielemR;
        double rhoL, rhoR, uL, uR, vL, vR, wR, wL, pL, pR;

        ielemL = face2elem[iface][0] - 1;		// numero de l'element Gauche  [-1 NEEDED!?!?!]
        ielemR = face2elem[iface][1] - 1;		// numero de l'element droit   [-1 NEEDED!?!?!]

        //Update L/R values [order 2]
        rhoL = rho[ielemL]+gradient[ielemL][0][0]*cent2face[iface][0][0]+gradient[ielemL][1][0]*cent2face[iface][1][0]+gradient[ielemL][2][0]*cent2face[iface][2][0];
        uL = u[ielemL]+gradient[ielemL][0][1]*cent2face[iface][0][0]+gradient[ielemL][1][1]*cent2face[iface][1][0]+gradient[ielemL][2][1]*cent2face[iface][2][0];
        vL = v[ielemL]+gradient[ielemL][0][2]*cent2face[iface][0][0]+gradient[ielemL][1][2]*cent2face[iface][1][0]+gradient[ielemL][2][2]*cent2face[iface][2][0];
        wL = w[ielemL]+gradient[ielemL][0][3]*cent2face[iface][0][0]+gradient[ielemL][1][3]*cent2face[iface][1][0]+gradient[ielemL][2][3]*cent2face[iface][2][0];
        pL = p[ielemL]+gradient[ielemL][0][4]*cent2face[iface][0][0]+gradient[ielemL][1][4]*cent2face[iface][1][0]+gradient[ielemL][2][4]*cent2face[iface][2][0];

        rhoR = rho[ielemR]+gradient[ielemR][0][0]*cent2face[iface][0][1]+gradient[ielemR][1][0]*cent2face[iface][1][1]+gradient[ielemR][2][0]*cent2face[iface][2][1];
        uR = u[ielemR]+gradient[ielemR][0][1]*cent2face[iface][0][1]+gradient[ielemR][1][1]*cent2face[iface][1][1]+gradient[ielemR][2][1]*cent2face[iface][2][1];
        vR = v[ielemR]+gradient[ielemR][0][2]*cent2face[iface][0][1]+gradient[ielemR][1][2]*cent2face[iface][1][1]+gradient[ielemR][2][2]*cent2face[iface][2][1];
        wR = w[ielemR]+gradient[ielemR][0][3]*cent2face[iface][0][1]+gradient[ielemR][1][3]*cent2face[iface][1][1]+gradient[ielemR][2][3]*cent2face[iface][2][1];
        pR = p[ielemR]+gradient[ielemR][0][4]*cent2face[iface][0][1]+gradient[ielemR][1][4]*cent2face[iface][1][1]+gradient[ielemR][2][4]*cent2face[iface][2][1];
        UpwindFlux(iface,rhoL,uL,vL,wL,pL,rhoR,uR,vR,wR,pR);
        RoeDissipation(iface,rhoL,uL,vL,wL,pL,rhoR,uR,vR,wR,pR);  
    }   
}

void solver::UpwindFlux(int iface, double rhoL,double uL,double vL,double wL,double pL, double rhoR,double uR,double vR,double wR,double pR) {
        double  UL, UR, WR, WL, HL, HR, nx, ny, nz;
        double rhoAvg, uAvg, vAvg, wAvg, pAvg, Vavg, Havg;
        double Fcmass, Fcmom1, Fcmom2,Fcmom3, Fcenergy;
        double gamma = 1.4;
        //Calcul des normales
        nx = face2norm[iface][0];  
        ny = face2norm[iface][1];
        nz = face2norm[iface][2];

        UL = sqrt(uL * uL + vL * vL + wL*wL); //ajouter w??
        UR = sqrt(uR * uR + vR * vR + wR*wR);

        //Calcul des moyennes des variables 
        rhoAvg = 0.5 * (rhoL + rhoR);
        uAvg = 0.5 * (uL + uR);
        vAvg = 0.5 * (vL + vR);
        wAvg = 0.5 * (wL + wR);
        pAvg = 0.5 * (pL + pR);
        Vavg = uAvg * nx + vAvg * ny + wAvg*nz;  

        HL = 0.5 * UL * UL + pL /rhoL /(gamma - 1) + pL / rhoL;
        HR = 0.5 * UR * UR + pR /rhoR /(gamma - 1) + pR / rhoR;
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
        flux_c[iface][4] = Fcenergy; //e
}

void solver::RoeDissipation(int iface, double rhoL,double uL,double vL,double wL,double pL, double rhoR,double uR,double vR,double wR,double pR) {
    double dp, du, dv, dV, drho,dw;
    double VL, VR, UL, UR, cL, cR, HL, HR, nx, ny, nz, rhobar, ubar, vbar, wbar, hbar, cbar, Vbar, qbar, SR1, SR2, SR3, delta;
    double F1mass, F1mom1, F1mom2, F1mom3, F1energy, F234mass, F234mom1, F234mom2,F234mom3, F234energy, F5mass, F5mom1, F5mom2, F5mom3, F5energy;
    double AWW1, AWW2, AWW3, AWW4,AWW5;
    double c1, c2, c3, c4, c5;
    double gamma=1.4;

    drho = rhoR - rhoL;
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
    dp = pR - pL;

    HL = 0.5 * UL * UL + pL / rhoL / (gamma - 1) + pL / rhoL;
    HR = 0.5 * UR * UR + pR / rhoR / (gamma - 1) + pR / rhoR;

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
    
    flux_d[iface][0] = AWW1;   //rho
    flux_d[iface][1] = AWW2;   //u
    flux_d[iface][2] = AWW3;   //v
    flux_d[iface][3] = AWW4;   //w
    flux_d[iface][4] = AWW5;   //e
}


void solver::ComputeResidu() {

}

void solver::ComputeResiduConv() {

}

// Check if converged
double solver::CheckConvergence() {
    double residuSum = 0;
    for (int ielem=0;ielem<nelem;++ielem) {                                                 //Sum of residu in rho
        residuSum += pow((residu_c[ielem][0]-residu_d[ielem][0])/elem2vol[ielem],2);        // CAREFULL WITH SIGN OF DISSIPATIVE FLUX
    }
    return sqrt(residuSum);
}

// Conversion from pressure to energy
double solver::P2E(double p_loc,double rho_loc,double u_loc,double v_loc,double w_loc) {
	double e_loc;double lSurgammaM1 = 2.5;
    e_loc = p_loc*lSurgammaM1/rho_loc+0.5*sqrt(pow(u_loc,2)+pow(v_loc,2)+pow(w_loc,2));
	return e_loc;
}
	
// Conversion from energy to pressure
double solver::E2P(double e_loc,double rho_loc,double u_loc,double v_loc,double w_loc) {
	double p_loc;double gammaM1 = 0.4;
    p_loc = gammaM1*rho_loc*(e_loc-0.5*sqrt(pow(u_loc,2)+pow(v_loc,2)+pow(w_loc,2)));
	return p_loc;
}
	


void solver::ResidualSmoothing() {
    
}