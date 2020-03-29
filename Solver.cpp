#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <Solver.h>
#include <math.h>
#include <main.h>

// constructor
solver_c::solver_c(Reader_c& FileContents, double convergeFixLimit_in)
{
  convergeFixLimit = convergeFixLimit_in;
  string parametreFile = "CFDsimPI4.txt";
  cout << "---------- Parametres... ----------\n";
  FileContents.computePrmt(parametreFile);
  cout << "**************\nEnd\n**************\n";
  smoothOrNah = FileContents.Smoothing;
  mach = FileContents.mach;
  AoA = FileContents.AoA;
  cfl = FileContents.cfl;
  if(FileContents.tempMethod=="Runge-Kutta")
  {
    RK_M=1;
    RK_step = FileContents.Nstage;
  }
  else if(FileContents.tempMethod=="Euler_explicite")
  {
    RK_step = 1;
  }
  else if(FileContents.tempMethod=="Runge-Kutta_Hybride")
  {
    RK_M = 0;
    RK_step = FileContents.Nstage;
  }
  else
  {
    cout << "ERROR METHOD NOT SUPPORTED";
    exit( 3 );
  }
  Order = FileContents.spatMethod_ordre;
  iterMax  = FileContents.iterMax;
  convergeCrit = FileContents.convCrit;
  inf_speed = mach*sqrt(1.4);
  inf_speed_x = inf_speed*cos(AoA);
  inf_speed_y = inf_speed*sin(AoA);
  inf_speed_z = 0;
  iteration = 0;
  World.Init();    // Initialise MPI
  ofstream myfile ("ResiduLog.txt");  //Flash log
  if (myfile.is_open())
  {
      myfile.close();
  }
}
// destructor
solver_c::~solver_c() {
// DELETE ALL CREATED ARRAYS
    delete[] bound2tag;
    delete[] BoundIndex;
    delete[] ZBoundIndex;
    delete[] elem2vtk;


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
        delete[] SmootyRezi[i];
        delete[] residu_d_hyb[i];
        delete[] F_lux[i];
        delete[] SmootyRezi[i];
        delete[] W_0[i];
        delete[] cons[i];
        delete[] cons_[i];
    }
    delete[] flux_c;
    delete[] flux_d;
    delete[] residu_c;
    delete[] residu_d;
    delete[] SmootyRezi;
    delete[] F_lux;
    delete[] W_0;
    delete[] cons;
    delete[] cons_;
    delete[] residu_d_hyb;

    if (Order==2) {
        for (int ielem=0;ielem<ncell;++ielem) {
            for (int idim=0;idim<ndime;++idim) {
                delete[] gradient[ielem][idim];
            }
            delete[] gradient[ielem];
        }
        delete[] gradient;
        for (int ielem=0;ielem<nelem;++ielem) {
            delete[] limit[ielem];
        }
        delete[] limit;
    }

    // Delete MPI-related variables
    for (int izone=0;izone<ntgt;++izone) {
        if (Order==2){delete[] gradientSendBuffer[izone];}
        delete[] primitivesSendBuffer[izone];
    }
    if (Order==2){delete[] gradientSendBuffer;}
    delete[] primitivesSendBuffer;
}


// call all other methods in order while not converged
void solver_c::Compute() {
    double Residu = pow(69,42);
    double Old_Residu = pow(69,42);
    double Residu_initial,ResiduLocal;
    Initialisation();
    UpdateBound();

    while (iteration<iterMax) {
        ++iteration;
        if (Order==1){ComputeFluxO1();}
        else {ComputeGrandientsNLimit();ExchangeGradiants();ComputeFluxO2();}
        ComputeResidu();
        ResiduLocal = CheckConvergence();
        Residu = World.UpdateConvergence(ResiduLocal);
        if (iteration==1) {
            Residu_initial = Residu;
        }
        residuRel = Residu/Residu_initial;
        if (World.world_rank==0) {
        WriteResidu();
        printf("Iteration : %d, \tREL R :  %e, \tABS R :  %e\n",iteration,residuRel,Residu);
        }
        if ((abs(Residu/Residu_initial)<convergeCrit)&&(abs(Old_Residu/Residu_initial)<convergeCrit)) {
            break;
        }
        Old_Residu = Residu;
        if (RK_step==1){TimeStepEul();}
        else if(RK_M==0) {TimeStepRkH();}
        else {TimeStepRkM();}
        ExchangePrimitive();
        UpdateBound();
    }
}

// set every element to infinity and initialise other stuff
void solver_c::Initialisation()
{
// initialise arrays
  rho = new double [ncell];
  u = new double [ncell];
  v = new double [ncell];
  w = new double [ncell];
  p = new double [ncell];
  if (Order==2) {
      gradient = new double**[ncell];     //Needs to be ncell as boundaries [between zones] need gradients
      limit = new double*[nelem];
      for (int ielem=0;ielem<ncell;++ielem) {
          gradient[ielem] = new double*[ndime];
          for (int idim=0;idim<ndime;++idim) {
              gradient[ielem][idim] = new double[5];
          }
      }
      for (int ielem=0;ielem<nelem;++ielem) {
          limit[ielem] = new double[5];
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
      u[i] = inf_speed_x;
      v[i] = inf_speed_y;
      w[i] = inf_speed_z;
      p[i] = 1.0;
  }
  W_0 = new double*[nelem];
  cons = new double*[nelem];
  cons_ = new double*[nelem];
  SmootyRezi = new double*[nelem];
  F_lux = new double*[nelem];
  residu_c = new double*[nelem];
  residu_d = new double*[nelem];
  residu_d_hyb = new double*[nelem];
  for (int i=0;i<nelem;++i)
  {
    W_0[i] = new double [5];
    cons[i] = new double [5];
    cons_[i] = new double [5];
    SmootyRezi[i] = new double [5];
    F_lux[i] = new double [5];
    residu_c[i] = new double [5];
    residu_d[i] = new double [5];
    residu_d_hyb[i] = new double [5];
  }
  SaveConservative();
}

// update boundary conditions
void solver_c::UpdateBound() {
 int ig,ir,indxMin,indxMax,iface;// ig = ghost index, ir = real index
    double c0,l_r0c0,r0c0,udotn; //local sound speed
    string BoundType;

    for (int ibc=0;ibc<nbc;++ibc) {
        BoundType = bound2tag[ibc];
        indxMin = BoundIndex[ibc];
        indxMax = BoundIndex[ibc+1];
        if (BoundType.substr(0,4)=="wall"){ //If slipwall
            // Slipwall
            for (ig=indxMin;ig<indxMax;++ig) {
                iface = elem2face[ig][0];
                ir = elem2elem[ig][0];
                udotn = 2.0*(face2norm[iface][0]*u[ir]+face2norm[iface][1]*v[ir]+face2norm[iface][2]*w[ir]);
                p[ig] = p[ir];  //There are higher accruacy methods
                rho[ig] = rho[ir];
                u[ig] = u[ir]-udotn*face2norm[iface][0];    // p.273, eq 8.10
                v[ig] = v[ir]-udotn*face2norm[iface][1];
                w[ig] = w[ir]-udotn*face2norm[iface][2];
            }
        }
        else if (BoundType.substr(0,3)=="sym") {//If symmetry
            // Symmetry aka slip wall
            for (ig=indxMin;ig<indxMax;++ig) {
                iface = elem2face[ig][0];
                ir = elem2elem[ig][0];
                udotn = 2.0*(face2norm[iface][0]*u[ir]+face2norm[iface][1]*v[ir]+face2norm[iface][2]*w[ir]);
                p[ig] = p[ir];
                rho[ig] = rho[ir];
                u[ig] = u[ir]-udotn*face2norm[iface][0];    // p.273, eq 8.10
                v[ig] = v[ir]-udotn*face2norm[iface][1];
                w[ig] = w[ir]-udotn*face2norm[iface][2];
            }
        }
        else {
            for (ig=indxMin;ig<indxMax;++ig) {
                iface = elem2face[ig][0];
                ir = elem2elem[ig][0];

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
                        rho[ig] = 1.0+(p[ig]-1.0)/pow(c0,2);
                        u[ig] = inf_speed_x-face2norm[iface][0]*(1.0-p[ig])*l_r0c0;
                        v[ig] = inf_speed_y-face2norm[iface][1]*(1.0-p[ig])*l_r0c0;
                        w[ig] = inf_speed_z-face2norm[iface][2]*(1.0-p[ig])*l_r0c0;

                    }
                    else {
                        // Subsonic Outflow
                        p[ig] = 1.0;
                        rho[ig] = rho[ir]+(p[ig]-p[ir])/pow(c0,2);
                        u[ig] = u[ir]+face2norm[iface][0]*(p[ir]-p[ig])*l_r0c0;
                        v[ig] = v[ir]+face2norm[iface][1]*(p[ir]-p[ig])*l_r0c0;
                        w[ig] = w[ir]+face2norm[iface][2]*(p[ir]-p[ig])*l_r0c0;
                    }
                }
            }
        }
    }
}

// Initialise MPI
void solver_c::InitMPIBuffer(Reader_c& Read) {
    int* zone2nbelem,*tgtList,**localBorderID;   //number of element / zone boundary

    ntgt = Read.nzone;
    //Build zone2nbelem; number of element / zone boundary
    zone2nbelem = new int[ntgt];
    tgtList = new int[ntgt];    //Build communication target list
    for (int izone=0;izone<ntgt;++izone) {
        zone2nbelem[izone] = Read.z_nelemv[izone];// = TO NUMBER OF BOUNDARY ELEMENTS BETWEEN THE 2 ZONES
        tgtList[izone] = stoi(Read.zone2tag[izone]);// = ZONE ID
    }
    World.InitBuffer(ntgt,tgtList,zone2nbelem);

    // Create transmition buffers
    if (Order==2){gradientSendBuffer = new double*[World.ntgt];}
    primitivesSendBuffer = new double*[World.ntgt];
    for (int izone=0;izone<World.ntgt;++izone) {
        if (Order==2){gradientSendBuffer[izone] = new double[zone2nbelem[izone]*15];}
        primitivesSendBuffer[izone] = new double[zone2nbelem[izone]*5];
    }

    // Handshake order in which boundary elements will be sent between zones
    localBorderID = new int*[ntgt];
    int IndxMax,IndxMin;
	for (int izone=0;izone<World.ntgt;++izone) {
		IndxMin = Read.zoneIndex[izone];
		IndxMax = Read.zoneIndex[izone+1];
        localBorderID[izone] = new int[zone2nbelem[izone]];
		for (int j=IndxMin;j<IndxMax;++j){
			localBorderID[izone][j-IndxMin] = Read.zelem2jelem[j];
		}
	}
    World.ExchangeCellOrder(localBorderID);

    for (int izone=0;izone<World.ntgt;++izone) {
        delete[] localBorderID[izone];
	}
    delete[] localBorderID;
    delete[] zone2nbelem;delete[] tgtList;


    // Copy needed informations from Reader
    nbc = Read.nbc;
    bound2tag = new string[nbc];
    BoundIndex = new int[nbc+1];
    for (int ibc=0;ibc<nbc;++ibc) {
        BoundIndex[ibc] = Read.BoundIndex[ibc];
        bound2tag[ibc] = Read.bound2tag[ibc];
    }
    BoundIndex[nbc] = Read.BoundIndex[nbc];

    int indx1,indx2=BoundIndex[nbc];
    ZBoundIndex = new int[nzone+1];
	for (int izone=0;izone<nzone;++izone) {
		indx1 = indx2;
		indx2 += Read.z_nelemv[izone];
        ZBoundIndex[izone] = indx1;
	}
    ZBoundIndex[nzone] = indx2;
    elem2vtk = new int[ncell];
    for(int ielem=0;ielem<ncell;++ielem) {
        elem2vtk[ielem] = Read.elem2vtk[ielem];
    }
    //To use ExchangeMetrics, a buffer needs to be created first...
    if (Order==2){ExchangeMetrics();}
}

// exchange needed metrics values between zones
void solver_c::ExchangeMetrics() {
    int ielem,iface,ibelemIndx,jbelemIndx,Nbr_of_faceIndx,BaseZoneIndxMin,BaseZoneIndxMax; //Local boundary element index
    int Pre_ind=ndime-2;
    // Create needed buffer
    double **metricBuffer;
    metricBuffer = new double*[World.ntgt];
    for (int izone=0;izone<World.ntgt;++izone) {
        metricBuffer[izone] = new double[World.zone2nbelem[izone]*3];
    }
    // face2elemCenter : Only need to send [:][0][:]     || face2elemCenter[iface][0:side 0,1:side 1][idim]
        // Populate Tx Buffer
    for (int izone=0;izone<World.ntgt;++izone) {
        for (int ibelem=0;ibelem<World.zone2nbelem[izone];++ibelem) {
            ibelemIndx = ibelem + ZBoundIndex[izone];
            iface = elem2face[ibelemIndx][0];
            metricBuffer[izone][ibelem] = face2elemCenter[iface][0][0];
            metricBuffer[izone][ibelem+World.zone2nbelem[izone]] = face2elemCenter[iface][0][1];
            metricBuffer[izone][ibelem+World.zone2nbelem[izone]*2] = face2elemCenter[iface][0][2];
        }
    }
        //Send & Store
    World.Exchange1DBuffer(metricBuffer);
    for (int izone=0;izone<World.ntgt;++izone) {
        BaseZoneIndxMin = ZBoundIndex[izone];
        BaseZoneIndxMax = ZBoundIndex[izone+1];
        for (int ibelem=0;ibelem<World.zone2nbelem[izone];++ibelem) {
        // HIGHLY INNEFICIENT, SENDING DIRECTLY THE GHOST ID WOULD BE MUCHHHHHHHHHHH FASTER [Would allow use of ReclassPrimitives]
            jbelemIndx = World.rxOrder2localOrder[izone][ibelem];   //Real Element in current zone

            iface = elem2face[jbelemIndx][0];
            face2elemCenter[iface][1][0] = World.oneDBuffer[izone][ibelem];
            face2elemCenter[iface][1][1] = World.oneDBuffer[izone][ibelem+World.zone2nbelem[izone]];
            face2elemCenter[iface][1][2] = World.oneDBuffer[izone][ibelem+World.zone2nbelem[izone]*2];
        }
    }
    // Delete needed buffer
    for (int izone=0;izone<World.ntgt;++izone) {
        delete[] metricBuffer[izone];
    }
    delete[] metricBuffer;

}

// exchange primitive values between zones
void solver_c::ExchangePrimitive() {
    int ielem,ibelemIndx,jbelemIndx,Nbr_of_faceIndx,BaseZoneIndxMin,BaseZoneIndxMax; //Local boundary element index
    int Pre_ind=ndime-2;
    // Populate Tx Buffer
    for (int izone=0;izone<World.ntgt;++izone) {
        for (int ibelem=0;ibelem<World.zone2nbelem[izone];++ibelem) {
            ibelemIndx = ibelem + ZBoundIndex[izone];
            ielem = elem2elem[ibelemIndx][0];
            primitivesSendBuffer[izone][ibelem] = rho[ielem];
            primitivesSendBuffer[izone][ibelem+World.zone2nbelem[izone]] = u[ielem];
            primitivesSendBuffer[izone][ibelem+World.zone2nbelem[izone]*2] = v[ielem];
            primitivesSendBuffer[izone][ibelem+World.zone2nbelem[izone]*3] = w[ielem];
            primitivesSendBuffer[izone][ibelem+World.zone2nbelem[izone]*4] = p[ielem];
        }
    }
    //Send & Store
    World.ExchangePrimitives(primitivesSendBuffer);
//    World.ReclassPrimitives(&rho,&u,&v,&w,&p);
// PROBLEMATIC IF AN ELEMENTS SHARE TWO NEIGBHOORS IN THE SAME ZONE---------------------------------------------
    for (int izone=0;izone<World.ntgt;++izone) {
        BaseZoneIndxMin = ZBoundIndex[izone];
        BaseZoneIndxMax = ZBoundIndex[izone+1];
        for (int ibelem=0;ibelem<World.zone2nbelem[izone];++ibelem) {
        // HIGHLY INNEFICIENT, SENDING DIRECTLY THE GHOST ID WOULD BE MUCHHHHHHHHHHH FASTER [Would allow use of ReclassPrimitives]
            jbelemIndx = World.rxOrder2localOrder[izone][ibelem];   //Real Element in current zone

            rho[jbelemIndx] = World.primitivesBuffer[izone][ibelem];
            u[jbelemIndx] = World.primitivesBuffer[izone][ibelem+World.zone2nbelem[izone]];
            v[jbelemIndx] = World.primitivesBuffer[izone][ibelem+World.zone2nbelem[izone]*2];
            w[jbelemIndx] = World.primitivesBuffer[izone][ibelem+World.zone2nbelem[izone]*3];
            p[jbelemIndx] = World.primitivesBuffer[izone][ibelem+World.zone2nbelem[izone]*4];
        }
    }
}

// exchange gradiant values between zones [Order 2 only]
void solver_c::ExchangeGradiants() {    // TO EDIT-----------------------------
    int ielem,ibelemIndx,jbelemIndx,Nbr_of_faceIndx,BaseZoneIndxMin,BaseZoneIndxMax; //Local boundary element index
    int Pre_ind=ndime-2;
    // Populate Tx Buffer
    for (int izone=0;izone<World.ntgt;++izone) {
        for (int ibelem=0;ibelem<World.zone2nbelem[izone];++ibelem) {
            for (int idim;idim<3;++idim) {
                ibelemIndx = ibelem + ZBoundIndex[izone];
                ielem = elem2elem[ibelemIndx][0];
                gradientSendBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5] = gradient[ielem][idim][0];
                gradientSendBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5+World.zone2nbelem[izone]*1] = gradient[ielem][idim][1];
                gradientSendBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5+World.zone2nbelem[izone]*2] = gradient[ielem][idim][2];
                gradientSendBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5+World.zone2nbelem[izone]*3] = gradient[ielem][idim][3];
                gradientSendBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5+World.zone2nbelem[izone]*4] = gradient[ielem][idim][4];
            }
        }
    }
    //Send & Store
    World.ExchangeGradients(gradientSendBuffer);
// PROBLEMATIC IF AN ELEMENTS SHARE TWO NEIGBHOORS IN THE SAME ZONE---------------------------------------------
    for (int izone=0;izone<World.ntgt;++izone) {
        BaseZoneIndxMin = ZBoundIndex[izone];
        BaseZoneIndxMax = ZBoundIndex[izone+1];
        for (int ibelem=0;ibelem<World.zone2nbelem[izone];++ibelem) {
        // HIGHLY INNEFICIENT, SENDING DIRECTLY THE GHOST ID WOULD BE MUCHHHHHHHHHHH FASTER [Would allow use of ReclassPrimitives]
            jbelemIndx = World.rxOrder2localOrder[izone][ibelem];   //Real Element in current zone

            for (int idim;idim<3;++idim) {
                gradient[jbelemIndx][idim][0] = World.gradientBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5];
                gradient[jbelemIndx][idim][1] = World.gradientBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5+World.zone2nbelem[izone]*1];
                gradient[jbelemIndx][idim][2] = World.gradientBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5+World.zone2nbelem[izone]*2];
                gradient[jbelemIndx][idim][3] = World.gradientBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5+World.zone2nbelem[izone]*3];
                gradient[jbelemIndx][idim][4] = World.gradientBuffer[izone][ibelem+World.zone2nbelem[izone]*idim*5+World.zone2nbelem[izone]*4];
            }
        }
    }
}
void solver_c::SaveFlux()
{
  for(int ielem=0; ielem<nelem; ++ielem)
  {
    F_lux[ielem][0] = (residu_c[ielem][0]-residu_d[ielem][0]);
    F_lux[ielem][1] = (residu_c[ielem][1]-residu_d[ielem][1]);
    F_lux[ielem][2] = (residu_c[ielem][2]-residu_d[ielem][2]);
    F_lux[ielem][3] = (residu_c[ielem][3]-residu_d[ielem][3]);
    F_lux[ielem][4] = (residu_c[ielem][4]-residu_d[ielem][4]);
  }
}

void solver_c::SaveConservative()
{
  for(int ielem=0; ielem<nelem; ielem++)
  {
    cons[ielem][0] = rho[ielem];
    cons[ielem][1] = rho[ielem]*u[ielem];
    cons[ielem][2] = rho[ielem]*v[ielem];
    cons[ielem][3] = rho[ielem]*w[ielem];
    cons[ielem][4] = rho[ielem]*P2E(p[ielem],rho[ielem],u[ielem],v[ielem],w[ielem]);
    cons_[ielem][0] = cons[ielem][0];
    cons_[ielem][1] = cons[ielem][1];
    cons_[ielem][2] = cons[ielem][2];
    cons_[ielem][3] = cons[ielem][3];
    cons_[ielem][4] = cons[ielem][4];
  }
}
void solver_c::SavePrimitive(int ielem)
{
  rho[ielem] = cons[ielem][0];
  invrho = 1.0/rho[ielem];
  u[ielem] = cons[ielem][1]*invrho;
  v[ielem] = cons[ielem][2]*invrho;
  w[ielem] = cons[ielem][3]*invrho;
  eTempo = cons[ielem][4]*invrho;
  p[ielem] = E2P(eTempo,rho[ielem],u[ielem],v[ielem],w[ielem]);
}
void solver_c::SavePrimitiveRK(int ielem)
{
  rho[ielem] = cons_[ielem][0];
  invrho = 1.0/rho[ielem];
  u[ielem] = cons_[ielem][1]*invrho;
  v[ielem] = cons_[ielem][2]*invrho;
  w[ielem] = cons_[ielem][3]*invrho;
  eTempo = cons_[ielem][4]*invrho;
  p[ielem] = E2P(eTempo,rho[ielem],u[ielem],v[ielem],w[ielem]);
}

// euler explicit time integration
void solver_c::TimeStepEul()
{
  double c;           //local sound speed
  double sumLambda;   // see blasek, sum of spectral radiuses
  double dTi_sans_V;         //dTi local time step
  SaveFlux();
  ResidualSmoothing();
  for (int ielem=0;ielem<nelem;++ielem)
  {
    c = sqrt(1.4*p[ielem]/rho[ielem]);
    eTempo = P2E(p[ielem],rho[ielem],u[ielem],v[ielem],w[ielem]);
    sumLambda = (fabs(u[ielem])+c)*elem2deltaSxyz[ielem][0]+(fabs(v[ielem])+c)*elem2deltaSxyz[ielem][1]+(fabs(w[ielem])+c)*elem2deltaSxyz[ielem][2];
    dTi_sans_V = cfl/sumLambda;
    cons[ielem][0] -= (F_lux[ielem][0])*dTi_sans_V;
    cons[ielem][1] -= (F_lux[ielem][1])*dTi_sans_V;
    cons[ielem][2] -= (F_lux[ielem][2])*dTi_sans_V;
    cons[ielem][3] -= (F_lux[ielem][3])*dTi_sans_V;
    cons[ielem][4] -= (F_lux[ielem][4])*dTi_sans_V;
    SavePrimitive(ielem);
  }
}

// Runge-Kutta Multistage time integration
void solver_c::TimeStepRkM()
{
  int OIndx = Order-1;
  double c,invrho;           //local sound speed
  double sumLambda;   // see blasek, sum of spectral radiuses
  double dTi_sans_V;         //dTi local time step
  for (int k=0;k<RK_step-1;++k)
  {
    SaveFlux();
    ResidualSmoothing();
    for (int ielem=0;ielem<nelem;++ielem)
    {
      c = sqrt(1.4*p[ielem]/rho[ielem]);
      sumLambda = (fabs(u[ielem])+c)*elem2deltaSxyz[ielem][0]+(fabs(v[ielem])+c)*elem2deltaSxyz[ielem][1]+(fabs(w[ielem])+c)*elem2deltaSxyz[ielem][2];
      dTi_sans_V = cfl/sumLambda;

      cons_[ielem][0] = cons[ielem][0] - dTi_sans_V*RKM_coef[RK_step][OIndx][k]*(F_lux[ielem][0]);
      cons_[ielem][1] = cons[ielem][1] - dTi_sans_V*RKM_coef[RK_step][OIndx][k]*(F_lux[ielem][1]);
      cons_[ielem][2] = cons[ielem][2] - dTi_sans_V*RKM_coef[RK_step][OIndx][k]*(F_lux[ielem][2]);
      cons_[ielem][3] = cons[ielem][3] - dTi_sans_V*RKM_coef[RK_step][OIndx][k]*(F_lux[ielem][3]);
      cons_[ielem][4] = cons[ielem][4] - dTi_sans_V*RKM_coef[RK_step][OIndx][k]*(F_lux[ielem][4]);
      SavePrimitiveRK(ielem);
    }
    // update residu
    if ((k+1)!=RK_step)
    {
      ExchangePrimitive();
      UpdateBound();
      if (Order==1){ComputeFluxO1();}
      else {ComputeGrandientsNLimit();ExchangeGradiants();ComputeFluxO2();}
      ComputeResidu();
    }
  }
  SaveFlux();
  ResidualSmoothing();
  for (int ielem=0;ielem<nelem;++ielem)
  {
    c = sqrt(1.4*p[ielem]/rho[ielem]);
    sumLambda = (fabs(u[ielem])+c)*elem2deltaSxyz[ielem][0]+(fabs(v[ielem])+c)*elem2deltaSxyz[ielem][1]+(fabs(w[ielem])+c)*elem2deltaSxyz[ielem][2];
    dTi_sans_V = cfl/sumLambda;

    cons[ielem][0] = cons[ielem][0] - dTi_sans_V*RKM_coef[RK_step][OIndx][RK_step - 1]*(F_lux[ielem][0]);
    cons[ielem][1] = cons[ielem][1] - dTi_sans_V*RKM_coef[RK_step][OIndx][RK_step - 1]*(F_lux[ielem][1]);
    cons[ielem][2] = cons[ielem][2] - dTi_sans_V*RKM_coef[RK_step][OIndx][RK_step - 1]*(F_lux[ielem][2]);
    cons[ielem][3] = cons[ielem][3] - dTi_sans_V*RKM_coef[RK_step][OIndx][RK_step - 1]*(F_lux[ielem][3]);
    cons[ielem][4] = cons[ielem][4] - dTi_sans_V*RKM_coef[RK_step][OIndx][RK_step - 1]*(F_lux[ielem][4]);
    SavePrimitive(ielem);
  }
}

// Runge-Kutta Hybrid time integration
void solver_c::TimeStepRkH() {
    double c, eTempo, invrho;           //local sound speed
    double sumLambda;   // see blasek, sum of spectral radiuses
    double dTi_sans_V;  //dTi local time step
    int kk;
    // update W_0 with initial results
    SaveConservative();

    for (int k=0;k<RK_step-1;++k)
    {
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
        SaveFlux_Hyb();
        ResidualSmoothing();
        // do step
        for (int ielem=0;ielem<nelem;++ielem)
        {
          c = sqrt(1.4*p[ielem]/rho[ielem]);
          sumLambda = (fabs(u[ielem])+c)*elem2deltaSxyz[ielem][0]+(fabs(v[ielem])+c)*elem2deltaSxyz[ielem][1]+(fabs(w[ielem])+c)*elem2deltaSxyz[ielem][2];
          dTi_sans_V = cfl/sumLambda;

          cons_[ielem][0] = cons[ielem][0] - RKH_coef[RK_step][0][k]*dTi_sans_V*(F_lux[ielem][0]);
          cons_[ielem][1] = cons[ielem][1] - RKH_coef[RK_step][0][k]*dTi_sans_V*(F_lux[ielem][1]);
          cons_[ielem][2] = cons[ielem][2] - RKH_coef[RK_step][0][k]*dTi_sans_V*(F_lux[ielem][2]);
          cons_[ielem][3] = cons[ielem][3] - RKH_coef[RK_step][0][k]*dTi_sans_V*(F_lux[ielem][3]);
          cons_[ielem][4] = cons[ielem][4] - RKH_coef[RK_step][0][k]*dTi_sans_V*(F_lux[ielem][4]);
          SavePrimitiveRK(ielem);
        }
        // update residu
        if ((k+1)!=RK_step)
        {
          ExchangePrimitive();
          UpdateBound();
    			if ((k==1) || (k==3))
          {
    				if (Order==1){ComputeFluxO1();}
            else {ComputeGrandientsNLimit();ExchangeGradiants();ComputeFluxO2();}  // Diffusive fluxes do not need to be computed at every step!
            ComputeResidu();
          }
    			else
          {
    				if (Order==1){ComputeFluxO1Conv();}
            else {ComputeGrandientsNLimit();ExchangeGradiants();ComputeFluxO2Conv();}  // Diffusive fluxes do not need to be computed at every step!
    			  ComputeResiduConv();
            }
        }
        kk = k;
    }
    kk += 1;

    // update residu_d_hyb
    if (kk==0)
    {
        for (int ielem=0;ielem<nelem;++ielem)
        {
            for (int jres=0;jres<5;++jres)
            {
                residu_d_hyb[ielem][jres] = residu_d[ielem][jres];
            }
        }
    }
    else if((kk==2)||(kk==4))
    {
        for (int ielem=0;ielem<nelem;++ielem)
        {
            for (int jres=0;jres<5;++jres)
            {
                residu_d_hyb[ielem][jres] = RKH_coef[RK_step][1][kk]*residu_d[ielem][jres]+(1-RKH_coef[RK_step][1][kk])*residu_d_hyb[ielem][jres];
            }
        }
    }
    SaveFlux_Hyb();
    ResidualSmoothing();
    for (int ielem=0;ielem<nelem;++ielem)
    {
      c = sqrt(1.4*p[ielem]/rho[ielem]);
      sumLambda = (fabs(u[ielem])+c)*elem2deltaSxyz[ielem][0]+(fabs(v[ielem])+c)*elem2deltaSxyz[ielem][1]+(fabs(w[ielem])+c)*elem2deltaSxyz[ielem][2];
      dTi_sans_V = cfl/sumLambda;

      cons[ielem][0] = cons[ielem][0] - RKH_coef[RK_step][0][kk]*dTi_sans_V*(F_lux[ielem][0]);
      cons[ielem][1] = cons[ielem][1] - RKH_coef[RK_step][0][kk]*dTi_sans_V*(F_lux[ielem][1]);
      cons[ielem][2] = cons[ielem][2] - RKH_coef[RK_step][0][kk]*dTi_sans_V*(F_lux[ielem][2]);
      cons[ielem][3] = cons[ielem][3] - RKH_coef[RK_step][0][kk]*dTi_sans_V*(F_lux[ielem][3]);
      cons[ielem][4] = cons[ielem][4] - RKH_coef[RK_step][0][kk]*dTi_sans_V*(F_lux[ielem][4]);
      SavePrimitive(ielem);
    }
    if ((kk+1)!=RK_step)
    {
      ExchangePrimitive();
      UpdateBound();
      if ((kk==1) || (kk==3))
      {
        if (Order==1){ComputeFluxO1();}
        else {ComputeGrandientsNLimit();ExchangeGradiants();ComputeFluxO2();}  // Diffusive fluxes do not need to be computed at every step!
        ComputeResidu();
      }
      else
      {
        if (Order==1){ComputeFluxO1Conv();}
        else {ComputeGrandientsNLimit();ExchangeGradiants();ComputeFluxO2Conv();}  // Diffusive fluxes do not need to be computed at every step!
        ComputeResiduConv();
      }
    }
}

 /*
void solver_c::TimeStepRkH() {
    double c, eTempo, invrho;           //local sound speed
    double sumLambda;   // see blasek, sum of spectral radiuses
    double dTi_sans_V;  //dTi local time step
    // update W_0 with initial results
    for (int ielem=0;ielem<nelem;++ielem) {
        W_0[ielem][0] = rho[ielem];
        W_0[ielem][1] = rho[ielem]*u[ielem];
        W_0[ielem][2] = rho[ielem]*v[ielem];
        W_0[ielem][3] = rho[ielem]*w[ielem];
        W_0[ielem][4] = rho[ielem]*P2E(p[ielem],rho[ielem],u[ielem],v[ielem],w[ielem]);
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
            sumLambda = (fabs(u[ielem])+c)*elem2deltaSxyz[ielem][0]+(fabs(v[ielem])+c)*elem2deltaSxyz[ielem][1]+(fabs(w[ielem])+c)*elem2deltaSxyz[ielem][2];
            dTi_sans_V = cfl/sumLambda;
            rho[ielem] = W_0[ielem][0] - RKH_coef[RK_step][0][k]*dTi_sans_V*(residu_c[ielem][0]-residu_d_hyb[ielem][0]);
			invrho=1/rho[ielem];
            u[ielem] = (W_0[ielem][1] - RKH_coef[RK_step][0][k]*dTi_sans_V*(residu_c[ielem][1]-residu_d_hyb[ielem][1]))*invrho;
            v[ielem] = (W_0[ielem][2] - RKH_coef[RK_step][0][k]*dTi_sans_V*(residu_c[ielem][2]-residu_d_hyb[ielem][2]))*invrho;
            w[ielem] = (W_0[ielem][3] - RKH_coef[RK_step][0][k]*dTi_sans_V*(residu_c[ielem][3]-residu_d_hyb[ielem][3]))*invrho;

			eTempo = (W_0[ielem][4] - RKH_coef[RK_step][0][k]*dTi_sans_V*(residu_c[ielem][4]-residu_d_hyb[ielem][4]))*invrho;
            p[ielem] = E2P(eTempo,rho[ielem],u[ielem],v[ielem],w[ielem]);
        }
        // update residu
        if ((k+1)!=RK_step) {
            ExchangePrimitive();
            UpdateBound();
			if ((k==1) || (k==3)) {
				if (Order==1){ComputeFluxO1();}
                else {ComputeGrandientsNLimit();ExchangeGradiants();ComputeFluxO2();}  // Diffusive fluxes do not need to be computed at every step!
                ComputeResidu();
            }
			else {
				if (Order==1){ComputeFluxO1Conv();}
                else {ComputeGrandientsNLimit();ExchangeGradiants();ComputeFluxO2Conv();}  // Diffusive fluxes do not need to be computed at every step!
			    ComputeResiduConv();
            }
        }
    }
}
 */


void solver_c::SaveFlux_Hyb()
{
  for(int ielem=0; ielem<nelem; ++ielem)
  {
    F_lux[ielem][0] = (residu_c[ielem][0]-residu_d_hyb[ielem][0]);
    F_lux[ielem][1] = (residu_c[ielem][1]-residu_d_hyb[ielem][1]);
    F_lux[ielem][2] = (residu_c[ielem][2]-residu_d_hyb[ielem][2]);
    F_lux[ielem][3] = (residu_c[ielem][3]-residu_d_hyb[ielem][3]);
    F_lux[ielem][4] = (residu_c[ielem][4]-residu_d_hyb[ielem][4]);
  }
}

void solver_c::ComputeFluxO1Conv() {
    for (int iface = 0; iface < nface; iface++) {
        int ielemL, ielemR;
        double rhoL, rhoR, uL, uR, vL, vR, wR, wL, pL, pR;

        ielemL = face2elem[iface][0];
        ielemR = face2elem[iface][1];

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
void solver_c::ComputeFluxO1() {
    for (int iface = 0; iface < nface; iface++) {
        int ielemL, ielemR;
        double rhoL, rhoR, uL, uR, vL, vR, wR, wL, pL, pR;

        ielemL = face2elem[iface][0];
        ielemR = face2elem[iface][1];

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
void solver_c::ComputeGrandientsNLimit() {
    int jelem,jface,Nbr_of_face,Pre_ind = ndime-2;
    double TempDemiSurVol,sign;
    double delta_2,UminRho,UmaxRho,UminU,UmaxU,UminV,UmaxV,UminP,UmaxP,UminW,UmaxW; //Barth & Jespersen p:166 Blasek
    int faceSide;   //Index [0/1]
    double MachinePrecision = pow(10,-15);
    double GradU[3][5] = {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
    double U[5];
    double UMax[5];
    double UMin[5];
    double psi[5] = {1,1,1,1,1};

    bool UpdateLim = ((residuRel>convergeFixLimit)||(iteration<10));

    for (int ielem=0;ielem<nelem;++ielem) {
        U[0] = rho[ielem];
        U[1] = u[ielem];
        U[2] = v[ielem];
        U[3] = w[ielem];
        U[4] = p[ielem];
        for (int ivar=0;ivar<5;++ivar) {
            UMax[ivar] = U[ivar];
            UMin[ivar] = U[ivar];
        }

        // Compute gradient
	    Nbr_of_face = vtklookup[Pre_ind][elem2vtk[ielem]][0];
        for (int ineighbor=0;ineighbor<Nbr_of_face;++ineighbor) {
            jelem = elem2elem[ielem][ineighbor];
            jface = elem2face[ielem][ineighbor];
            sign = ((double(face2elem[jface][0]==ielem))*2.0-1.0);  //jelem??

            for (int idime=0;idime<ndime;++idime) {
                GradU[idime][0] += sign*(rho[ielem]+rho[jelem])*face2norm[jface][idime]*face2area[jface];
                GradU[idime][1] += sign*(u[ielem]+u[jelem])*face2norm[jface][idime]*face2area[jface];
                GradU[idime][2] += sign*(v[ielem]+v[jelem])*face2norm[jface][idime]*face2area[jface];
                GradU[idime][3] += sign*(w[ielem]+w[jelem])*face2norm[jface][idime]*face2area[jface];
                GradU[idime][4] += sign*(p[ielem]+p[jelem])*face2norm[jface][idime]*face2area[jface];
            }

            //Update Umax/Umin
            UMin[0]=min(UMin[0],rho[jelem]);
            UMax[0]=max(UMax[0],rho[jelem]);
            UMin[1]=min(UMin[1],u[jelem]);
            UMax[1]=max(UMax[1],u[jelem]);
            UMin[2]=min(UMin[2],v[jelem]);
            UMax[2]=max(UMax[2],v[jelem]);
            UMin[3]=min(UMin[3],w[jelem]);
            UMax[3]=max(UMax[3],w[jelem]);
            UMin[4]=min(UMin[4],p[jelem]);
            UMax[4]=max(UMax[4],p[jelem]);
        }
        TempDemiSurVol = 0.5/elem2vol[ielem];
        for (int idime=0;idime<ndime;++idime){
            for (int ivar=0;ivar<5;++ivar) {
                GradU[idime][ivar] = GradU[idime][ivar]*TempDemiSurVol;
            }
        }

        if (UpdateLim) {
            // Compute Limiter
            for (int ineighbor=0;ineighbor<Nbr_of_face;++ineighbor) {   //For surrounding elements
                jelem = elem2elem[ielem][ineighbor];
                jface = elem2face[ielem][ineighbor];
                faceSide = int(face2elem[jface][0]==jelem);        //jelem??

                for (int ivar=0;ivar<5;++ivar) {
                    delta_2 = 0.5*(GradU[0][ivar]*face2elemCenter[jface][faceSide][0]+GradU[1][ivar]*face2elemCenter[jface][faceSide][1]+GradU[2][ivar]*face2elemCenter[jface][faceSide][2]);
                    if ((delta_2>0)) {
                        psi[ivar]  = min(psi[ivar],(UMax[ivar]-U[ivar])/(fabs(delta_2)+MachinePrecision));
                    }
                    else if ((delta_2<0)) {
                        psi[ivar]  = min(psi[ivar],-(UMin[ivar]-U[ivar])/(fabs(delta_2)+MachinePrecision));
                    }
                }
            }
            //store
            for (int iflux=0;iflux<5;++iflux) {
                limit[ielem][iflux] = psi[iflux];
            }
        }
        else {  //Load last-computed limitor
            for (int iflux=0;iflux<5;++iflux) {
                psi[iflux] = limit[ielem][iflux];
            }
        }

        //Store & reset Gradients
        for (int idim=0;idim<3;++idim) {
            for (int ivar=0;ivar<5;++ivar) {
                gradient[ielem][idim][ivar] = GradU[idim][ivar]*psi[ivar];
                GradU[idim][ivar]=0;
            }
        }
        //Reset limitors
        for (int ivar=0;ivar<5;++ivar) {
            psi[ivar] = 1;
        }
    }
}



// Roe fluxes, order 2 [REMEMBER TO SPLIT CONVECTIVE AND DIFFUSIVE FLUXES]
void solver_c::ComputeFluxO2Conv() {
for (int iface = 0; iface < nface; iface++) {
        int ielemL, ielemR;
        double rhoL, rhoR, uL, uR, vL, vR, wR, wL, pL, pR;

        ielemL = face2elem[iface][0];
        ielemR = face2elem[iface][1];

        //Update L/R values [order 2]
        rhoL = rho[ielemL]+gradient[ielemL][0][0]*face2elemCenter[iface][0][0]+gradient[ielemL][1][0]*face2elemCenter[iface][0][1]+gradient[ielemL][2][0]*face2elemCenter[iface][0][2];
        uL = u[ielemL]+gradient[ielemL][0][1]*face2elemCenter[iface][0][0]+gradient[ielemL][1][1]*face2elemCenter[iface][0][1]+gradient[ielemL][2][1]*face2elemCenter[iface][0][2];
        vL = v[ielemL]+gradient[ielemL][0][2]*face2elemCenter[iface][0][0]+gradient[ielemL][1][2]*face2elemCenter[iface][0][1]+gradient[ielemL][2][2]*face2elemCenter[iface][0][2];
        wL = w[ielemL]+gradient[ielemL][0][3]*face2elemCenter[iface][0][0]+gradient[ielemL][1][3]*face2elemCenter[iface][0][1]+gradient[ielemL][2][3]*face2elemCenter[iface][0][2];
        pL = p[ielemL]+gradient[ielemL][0][4]*face2elemCenter[iface][0][0]+gradient[ielemL][1][4]*face2elemCenter[iface][0][1]+gradient[ielemL][2][4]*face2elemCenter[iface][0][2];

        rhoR = rho[ielemR]+gradient[ielemR][0][0]*face2elemCenter[iface][1][0]+gradient[ielemR][1][0]*face2elemCenter[iface][1][1]+gradient[ielemR][2][0]*face2elemCenter[iface][1][2];
        uR = u[ielemR]+gradient[ielemR][0][1]*face2elemCenter[iface][1][0]+gradient[ielemR][1][1]*face2elemCenter[iface][1][1]+gradient[ielemR][2][1]*face2elemCenter[iface][1][2];
        vR = v[ielemR]+gradient[ielemR][0][2]*face2elemCenter[iface][1][0]+gradient[ielemR][1][2]*face2elemCenter[iface][1][1]+gradient[ielemR][2][2]*face2elemCenter[iface][1][2];
        wR = w[ielemR]+gradient[ielemR][0][3]*face2elemCenter[iface][1][0]+gradient[ielemR][1][3]*face2elemCenter[iface][1][1]+gradient[ielemR][2][3]*face2elemCenter[iface][1][2];
        pR = p[ielemR]+gradient[ielemR][0][4]*face2elemCenter[iface][1][0]+gradient[ielemR][1][4]*face2elemCenter[iface][1][1]+gradient[ielemR][2][4]*face2elemCenter[iface][1][2];
        UpwindFlux(iface,rhoL,uL,vL,wL,pL,rhoR,uR,vR,wR,pR);
    }
}

void solver_c::ComputeFluxO2() {
    for (int iface = 0; iface < nface; iface++) {
        int ielemL, ielemR;
        double rhoL, rhoR, uL, uR, vL, vR, wR, wL, pL, pR;

        ielemL = face2elem[iface][0];
        ielemR = face2elem[iface][1];

        //Update L/R values [order 2]
        rhoL = rho[ielemL]+gradient[ielemL][0][0]*face2elemCenter[iface][0][0]+gradient[ielemL][1][0]*face2elemCenter[iface][0][1]+gradient[ielemL][2][0]*face2elemCenter[iface][0][2];
        uL = u[ielemL]+gradient[ielemL][0][1]*face2elemCenter[iface][0][0]+gradient[ielemL][1][1]*face2elemCenter[iface][0][1]+gradient[ielemL][2][1]*face2elemCenter[iface][0][2];
        vL = v[ielemL]+gradient[ielemL][0][2]*face2elemCenter[iface][0][0]+gradient[ielemL][1][2]*face2elemCenter[iface][0][1]+gradient[ielemL][2][2]*face2elemCenter[iface][0][2];
        wL = w[ielemL]+gradient[ielemL][0][3]*face2elemCenter[iface][0][0]+gradient[ielemL][1][3]*face2elemCenter[iface][0][1]+gradient[ielemL][2][3]*face2elemCenter[iface][0][2];
        pL = p[ielemL]+gradient[ielemL][0][4]*face2elemCenter[iface][0][0]+gradient[ielemL][1][4]*face2elemCenter[iface][0][1]+gradient[ielemL][2][4]*face2elemCenter[iface][0][2];

        rhoR = rho[ielemR]+gradient[ielemR][0][0]*face2elemCenter[iface][1][0]+gradient[ielemR][1][0]*face2elemCenter[iface][1][1]+gradient[ielemR][2][0]*face2elemCenter[iface][1][2];
        uR = u[ielemR]+gradient[ielemR][0][1]*face2elemCenter[iface][1][0]+gradient[ielemR][1][1]*face2elemCenter[iface][1][1]+gradient[ielemR][2][1]*face2elemCenter[iface][1][2];
        vR = v[ielemR]+gradient[ielemR][0][2]*face2elemCenter[iface][1][0]+gradient[ielemR][1][2]*face2elemCenter[iface][1][1]+gradient[ielemR][2][2]*face2elemCenter[iface][1][2];
        wR = w[ielemR]+gradient[ielemR][0][3]*face2elemCenter[iface][1][0]+gradient[ielemR][1][3]*face2elemCenter[iface][1][1]+gradient[ielemR][2][3]*face2elemCenter[iface][1][2];
        pR = p[ielemR]+gradient[ielemR][0][4]*face2elemCenter[iface][1][0]+gradient[ielemR][1][4]*face2elemCenter[iface][1][1]+gradient[ielemR][2][4]*face2elemCenter[iface][1][2];
        UpwindFlux(iface,rhoL,uL,vL,wL,pL,rhoR,uR,vR,wR,pR);
        RoeDissipation(iface,rhoL,uL,vL,wL,pL,rhoR,uR,vR,wR,pR);
    }
}

void solver_c::UpwindFlux(int iface, double rhoL,double uL,double vL,double wL,double pL, double rhoR,double uR,double vR,double wR,double pR) {
    double EL, ER, nx, ny, nz;
    double rhoAvg, uAvg, vAvg, wAvg, pAvg, Vavg, Eavg;
    double Fcmass, Fcmom1, Fcmom2,Fcmom3, Fcenergy;
    //Calcul des normales
    nx = face2norm[iface][0];
    ny = face2norm[iface][1];
    nz = face2norm[iface][2];

    //Calcul des moyennes des variables
    rhoAvg = 0.5 * (rhoL + rhoR);
    uAvg = 0.5 * (uL + uR);
    vAvg = 0.5 * (vL + vR);
    wAvg = 0.5 * (wL + wR);
    pAvg = 0.5 * (pL + pR);
    Vavg = uAvg * nx + vAvg * ny + wAvg*nz;

    EL = P2E(pL,rhoL,uL,vL,wL);
    ER = P2E(pR,rhoR,uR,vR,wR);
    Eavg = 0.5 * (ER + EL);

    //Calcul du flux conservatif
    Fcmass = rhoAvg * Vavg;
    Fcmom1 = rhoAvg * Vavg * uAvg + pAvg * nx;
    Fcmom2 = rhoAvg * Vavg * vAvg + pAvg * ny;
    Fcmom3 = rhoAvg * Vavg * wAvg + pAvg * nz;
    Fcenergy = rhoAvg * Vavg * Eavg + pAvg * Vavg;

    // Flux dans les bonnes variables
    flux_c[iface][0] = Fcmass;   //rho
    flux_c[iface][1] = Fcmom1;   //u
    flux_c[iface][2] = Fcmom2;   //v
    flux_c[iface][3] = Fcmom3;   //w
    flux_c[iface][4] = Fcenergy; //e
}

void solver_c::RoeDissipation(int iface, double rhoL,double uL,double vL,double wL,double pL, double rhoR,double uR,double vR,double wR,double pR) {
    double dp, du, dv, dV, drho,dw,q2_t;
    double VL, VR, cL, cR, HL, HR, nx, ny, nz, rhobar, ubar, vbar, wbar, hbar, cbar, Vbar, qbar, SR1, SR2, SR3, delta;
    double F1mass, F1mom1, F1mom2, F1mom3, F1energy, F234mass, F234mom1, F234mom2,F234mom3, F234energy, F5mass, F5mom1, F5mom2, F5mom3, F5energy;
    double AWW1, AWW2, AWW3, AWW4,AWW5;
    double c1, c2, c3, c4, c5;
    double gamma=1.4,cmoy;

    drho = rhoR - rhoL;
    du = uR - uL;
    dv = vR - vL;
    dw = wR - wL;
    dp = pR - pL;

    //Calcul des normales
    nx = face2norm[iface][0];
    ny = face2norm[iface][1];
    nz = face2norm[iface][2];

    VL = nx * uL + ny * vL + nz*wL;     //ajouter w?
    VR = nx * uR + ny * vR + nz*wR;

    dV = VR - VL;


    HL = P2E(pL,rhoL,uL,vL,wL) + pL / rhoL;
    HR = P2E(pR,rhoR,uR,vR,wR) + pR / rhoR;


    //Calcul constantes pour simplifier la compilation
    c1=sqrt(rhoL);
    c2=sqrt(rhoR);
    c3= 1/(c1+c2);

    //Calcul des variables bar du schéma ROE
    rhobar = sqrt(rhoL * rhoR);
    ubar = (uL * c1 + uR * c2) * c3;
    vbar = (vL * c1 + vR * c2) * c3;
    wbar = (wL * c1 + wR * c2) * c3;
    hbar = (HL * c1 + HR * c2) * c3;
    qbar = (ubar * ubar) + (vbar * vbar)+(wbar*wbar); // c'est le qbar au carré
    cbar = sqrt((gamma - 1) * (hbar - qbar * 0.5));
    Vbar = ubar * nx + vbar * ny + wbar * nz;

    c4= 1/(2 * cbar * cbar);
    c5= 1/(cbar * cbar);

    cmoy = sqrt(gamma*(pL+pR)/(rhoL+rhoL));

    // Harten’s entropy correction
    SR1 = fabs(Vbar - cbar);
    SR2 = fabs(Vbar);
    SR3 = fabs(Vbar + cbar);

    delta = 0.1 * (cmoy);

    if (fabs(SR1) <= delta) {
        SR1 = (SR1 * SR1 + delta * delta) / (2 * delta);
    }
    if (fabs(SR2) <= delta) {
        SR2 = (SR2 * SR2 + delta * delta) / (2 * delta);
    }
    if (fabs(SR3) <= delta) {
        SR3 = (SR3 * SR3 + delta * delta) / (2 * delta);
    }

    //Calcul des différents termes du flux dissipatif
    F1mass = SR1 * ((dp - rhobar * cbar * dV) *c4) * 1;
    F1mom1 = (SR1 * ((dp - rhobar * cbar * dV)*c4)) * (ubar - (cbar * nx));
    F1mom2 = (SR1 * ((dp - rhobar * cbar * dV) *c4)) * (vbar - (cbar * ny));
    F1mom3 = (SR1 * ((dp - rhobar * cbar * dV) *c4)) * (wbar - (cbar * nz));
    F1energy = (SR1 * ((dp - rhobar * cbar * dV) *c4)) * (hbar - (cbar * Vbar));

    F234mass = SR2 * ((drho - (dp *c5)) );
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


void solver_c::ComputeResidu() {
    int iface,vtk,Nbr_of_face;
    double fluxSign;
    int Pre_ind = ndime-2;
    for(int ielem=0;ielem<nelem;ielem++) {
        //Reset Residu
        for (int iflux=0;iflux<5;++iflux) {
            residu_c[ielem][iflux] = 0;
            residu_d[ielem][iflux] = 0;
        }
        //Update it
        vtk = elem2vtk[ielem];
	    Nbr_of_face = vtklookup[Pre_ind][vtk][0];
        for (int jelemRel=0;jelemRel<Nbr_of_face;++jelemRel) {
            iface = elem2face[ielem][jelemRel];
            fluxSign = double(face2elem[iface][0]==ielem)*2.0-1.0;
            for (int iflux=0;iflux<5;++iflux) {
                residu_c[ielem][iflux] += fluxSign*flux_c[iface][iflux]*face2area[iface];
                residu_d[ielem][iflux] += fluxSign*flux_d[iface][iflux]*face2area[iface];
            }
        }
    }
}

void solver_c::ComputeResiduConv() {
    int iface,vtk,Nbr_of_face;
    double fluxSign;
    int Pre_ind = ndime-2;
    for(int ielem=0;ielem<nelem;ielem++) {
        //Reset Residu
        for (int iflux=0;iflux<5;++iflux) {
            residu_c[ielem][iflux] = 0;
        }
        //Update it
        vtk = elem2vtk[ielem];
	    Nbr_of_face = vtklookup[Pre_ind][vtk][0];
        for (int jelemRel=0;jelemRel<Nbr_of_face;++jelemRel) {
            iface = elem2face[ielem][jelemRel];
            fluxSign = double(face2elem[iface][0]==ielem)*2.0-1.0;
            for (int iflux=0;iflux<5;++iflux) {
                residu_c[ielem][iflux] += fluxSign*flux_c[iface][iflux]*face2area[iface];
            }
        }
    }
}

// Check if converged
double solver_c::CheckConvergence() {
    double residuSum = 0;
    for (int ielem=0;ielem<nelem;++ielem) {                                                 //Sum of residu in rho
        residuSum += pow((residu_c[ielem][0]-residu_d[ielem][0])/elem2vol[ielem],2);        // CAREFULL WITH SIGN OF DISSIPATIVE FLUX
    }
    return sqrt(residuSum);
}

// Conversion from pressure to energy
double solver_c::P2E(double p_loc,double rho_loc,double u_loc,double v_loc,double w_loc) {
	double e_loc;double lSurgammaM1 = 2.5;
    e_loc = p_loc*lSurgammaM1/rho_loc+0.5*(pow(u_loc,2)+pow(v_loc,2)+pow(w_loc,2));
	return e_loc;
}

// Conversion from energy to pressure
double solver_c::E2P(double e_loc,double rho_loc,double u_loc,double v_loc,double w_loc) {
	double p_loc;double gammaM1 = 0.4;
    p_loc = gammaM1*rho_loc*(e_loc-0.5*(pow(u_loc,2)+pow(v_loc,2)+pow(w_loc,2)));
	return p_loc;
}



void solver_c::ResidualSmoothing()
{
  if(smoothOrNah=="Yes")
  {
    int smooteration = 2;
    int NeiBoar;
    double S0, S1, S2, S3, S4;
    double epsilon = 0.6;
    for (int iter = 0; iter < smooteration; iter++)
		{
			for (int ielem = 0; ielem < nelem; ++ielem)
			{
        S0 = 0.0;        S1 = 0.0;        S2 = 0.0;
        S3 = 0.0;        S4 = 0.0;
        SmootyRezi[ielem][0] = 0.0;        SmootyRezi[ielem][1] = 0.0;
        SmootyRezi[ielem][2] = 0.0;        SmootyRezi[ielem][3] = 0.0;
        SmootyRezi[ielem][4] = 0.0;
        int vtk = elem2vtk[ielem];
	    	int Nbr_of_face = vtklookup[ndime-2][vtk][0];
        int FennTom = 0;
        for(int iface=0; iface<Nbr_of_face; ++iface)
        {
          NeiBoar = elem2elem[ielem][iface];
          if(NeiBoar<nelem)
          {
            S0 += epsilon * SmootyRezi[NeiBoar][0];
            S1 += epsilon * SmootyRezi[NeiBoar][1];
						S2 += epsilon * SmootyRezi[NeiBoar][2];
						S3 += epsilon * SmootyRezi[NeiBoar][3];
						S4 += epsilon * SmootyRezi[NeiBoar][4];
          }
          else
          {
            FennTom +=1;
          }
        }

				Nbr_of_face -= FennTom;
				SmootyRezi[ielem][0] = (F_lux[ielem][0] + S0) / (1.0 + epsilon * Nbr_of_face);
				SmootyRezi[ielem][1] = (F_lux[ielem][1] + S1) / (1.0 + epsilon * Nbr_of_face);
				SmootyRezi[ielem][2] = (F_lux[ielem][2] + S2) / (1.0 + epsilon * Nbr_of_face);
				SmootyRezi[ielem][3] = (F_lux[ielem][3] + S3) / (1.0 + epsilon * Nbr_of_face);
				SmootyRezi[ielem][4] = (F_lux[ielem][4] + S4) / (1.0 + epsilon * Nbr_of_face);
      }
    }
    for (int i = 0; i < nelem; ++i)
		{
			F_lux[i][0] = SmootyRezi[i][0];
			F_lux[i][1] = SmootyRezi[i][1];
			F_lux[i][2] = SmootyRezi[i][2];
			F_lux[i][3] = SmootyRezi[i][3];
			F_lux[i][4] = SmootyRezi[i][4];
		}
  }
}


void solver_c::WriteResidu(){
    ofstream myfile ("ResiduLog.txt",ios::app);
    if (myfile.is_open())
    {
        myfile << residuRel<<endl;
        myfile.close();
    }
}


void solver_c::PrintPress() {
    for (int ielem=0;ielem<ncell;++ielem) {
        cout<<p[ielem]<<endl;;
    }
}

void solver_c::HighlightZoneBorder() {
    int ielem,jbelemIndx,Nbr_of_faceIndx,BaseZoneIndxMin,BaseZoneIndxMax; //Local boundary element index
    int Pre_ind=ndime-2;
    // Populate Tx Buffer
    for (int izone=0;izone<World.ntgt;++izone) {
        for (int ibelem=0;ibelem<World.zone2nbelem[izone];++ibelem) {
            jbelemIndx = ibelem + ZBoundIndex[izone];
            ielem = elem2elem[jbelemIndx][0];
            primitivesSendBuffer[izone][ibelem] = 0;
            primitivesSendBuffer[izone][ibelem+World.zone2nbelem[izone]] = u[ielem];
            primitivesSendBuffer[izone][ibelem+World.zone2nbelem[izone]*2] = v[ielem];
            primitivesSendBuffer[izone][ibelem+World.zone2nbelem[izone]*3] = w[ielem];
            primitivesSendBuffer[izone][ibelem+World.zone2nbelem[izone]*4] = p[ielem];
            //rho[ielem] = 0;  //Sets THIS ZONE's borders to 0
        }
    }
    //Send & Store
    World.ExchangePrimitives(primitivesSendBuffer);
// // //    World.ReclassPrimitives(&rho,&u,&v,&w,&p);
// // // PROBLEMATIC IF AN ELEMENTS SHARE TWO NEIGBHOORS IN THE SAME ZONE---------------------------------------------
    for (int izone=0;izone<World.ntgt;++izone) {
        // BaseZoneIndxMin = ZBoundIndex[izone];
        // BaseZoneIndxMax = ZBoundIndex[izone+1];
        for (int ibelem=0;ibelem<World.zone2nbelem[izone];++ibelem) {
        // HIGHLY INNEFICIENT, SENDING DIRECTLY THE GHOST ID WOULD BE MUCHHHHHHHHHHH FASTER [Would allow use of ReclassPrimitives]
            ielem = World.rxOrder2localOrder[izone][ibelem];   //Real Element in current zone
            //cout<<World.world_rank<<"|"<<"From:"<<World.tgtList[izone]<<":=:"<<ielem<<endl;       //Good reception of SU2+ verified
	        // Nbr_of_faceIndx = vtklookup[0][elem2vtk[ielem]][0]-1;
            // for (int jelem=Nbr_of_faceIndx;jelem>-1;--jelem) {
            //     jbelemIndx = elem2elem[ielem][jelem];
            //     if((BaseZoneIndxMin<=jbelemIndx)&&(jbelemIndx<BaseZoneIndxMax)) {
            //         break;
            //     }
            // }
            rho[ielem] = World.primitivesBuffer[izone][ibelem];
        }
    }
}

void solver_c::SetAnalyticalGradiant(double dx,double dy,double dz) {
    int realElem;
    for (int ielem=0;ielem<nelem;++ielem) { //Real Elements Only!
        rho[ielem] = elem2center[ielem][0]*dx+elem2center[ielem][1]*dy+elem2center[ielem][2]*dz;
        u[ielem] = elem2center[ielem][0]*dx+elem2center[ielem][1]*dy+elem2center[ielem][2]*dz;
        v[ielem] = elem2center[ielem][0]*dx+elem2center[ielem][1]*dy+elem2center[ielem][2]*dz;
        w[ielem] = elem2center[ielem][0]*dx+elem2center[ielem][1]*dy+elem2center[ielem][2]*dz;
        p[ielem] = elem2center[ielem][0]*dx+elem2center[ielem][1]*dy+elem2center[ielem][2]*dz;
    }
    for (int ielem=nelem;ielem<ncell;++ielem) { //Set null gradiant at borders
        realElem = elem2elem[ielem][0];
        rho[ielem] = rho[realElem];
        u[ielem] = u[realElem];
        v[ielem] = v[realElem];
        w[ielem] = w[realElem];
        p[ielem] = p[realElem];
    }
    //Cheat so that limitors aren't computed and are used as 1
    residuRel = 0;
    iteration = 100;
    for (int ielem=0;ielem<nelem;ielem++) {
        for (int ivar=0;ivar<5;++ivar) {
            limit[ielem][ivar] = 1.0;
        }
    }
    ComputeGrandientsNLimit();
}

void solver_c::PrintGradiant() {
    int Nbr_of_face;bool hasGhost;
    for (int ielem=0;ielem<nelem;++ielem) { //Real Elements Only!
        hasGhost = 0;
	    Nbr_of_face = vtklookup[Order-2][elem2vtk[ielem]][0];
        for (int jface=Nbr_of_face-1;jface>0;jface--) {
            if (elem2elem[ielem][jface]>=nelem) {
                hasGhost = 1;
                //cout<<elem2elem[ielem][jface]<<endl;
                break;
            }
        }
        if (hasGhost){
            //cout<<World.world_rank<<"| dx:"<<gradient[ielem][0][0]<<" dy:"<<gradient[ielem][1][0]<<" dz:"<<gradient[ielem][2][0]<<" GHOSTED"<<endl;
        }
        else {
            cout<<World.world_rank<<"| dx:"<<gradient[ielem][0][0]<<" dy:"<<gradient[ielem][1][0]<<" dz:"<<gradient[ielem][2][0]<<endl;
        }
    }
}







void solver_c::PrintStylz() {
    // Mandatory ascii art
    system("Color 2B");
	cout<<"        "<< "\x1B[31m" << "|      ███▄    █ ▒█████   ██▓ ▄████▄ ▓█████   " << "\x1B[0m"<<"     "<< "\x1B[33m" << "______"<< "\x1B[0m"<<endl;
	cout<<"       "<< "\x1B[31m" << "/ \\     ██ ▀█   █▒██▒  ██▒▓██▒▒██▀ ▀█ ▓█   ▀  " << "\x1B[0m"<<" "<< "\x1B[33m" << ".---<__. \\ \\"<< "\x1B[0m"<<endl;
	cout<<"      "<< "\x1B[31m" << "/ "<< "\x1B[34m" << "_ " << "\x1B[0m"<<""<< "\x1B[31m" << "\\   ▓██  ▀█ ██▒██░  ██▒▒██▒▒▓█    ▄▒███    " << "\x1B[0m"<<" "<< "\x1B[33m" << "`---._  \\ \\ \\"<< "\x1B[0m"<<endl;
	cout<<"     "<< "\x1B[31m" << "|" << "\x1B[0m"<<"\x1B[34m"<<".o '."<< "\x1B[31m" << "|" << "\x1B[0m"<<"  "<< "\x1B[31m" << "▓██▒  ▐▌██▒██   ██░░██░▒▓▓▄ ▄██▒▓█  ▄   " << "\x1B[0m"<<" "<< "\x1B[33m" << ",----`- `.))"<< "\x1B[0m"<<endl;
	cout<<"     "<< "\x1B[31m" << "|" << "\x1B[0m"<<"\x1B[34m" <<"'._.'"<< "\x1B[31m" << "|" << "\x1B[0m"<<"  "<< "\x1B[31m" << "▒██░   ▓██░ ████▓▒░░██░▒ ▓███▀ ░▒████▒  " << "\x1B[0m"<<""<< "\x1B[33m" << "/ ,--.   )  |"<< "\x1B[0m"<<endl;
	cout<<"     "<< "\x1B[31m" << "|     |  ░ ▒░   ▒ ▒░ ▒░▒░▒░ ░▓  ░ ░▒ ▒  ░░ ▒░ ░  " << "\x1B[0m"<<""<< "\x1B[33m" << "/_/    >     |"<< "\x1B[0m"<<endl;
	cout<<"   "<< "\x1B[31m" << ",'|  |  |`.  ░ ░░   ░ ▒░ ░ ▒ ▒░  ▒ ░  ░  ▒   ░ ░  ░" << "\x1B[0m"<<" "<< "\x1B[33m" << "|,\\__-'      |"<< "\x1B[0m"<<endl;
	cout<<"  "<< "\x1B[31m" << "/  |  |  |  \\     ░   ░ ░░ ░ ░ ▒   ▒ ░░          ░ " << "\x1B[0m"<<"    "<< "\x1B[33m" << "\\_           \\"<< "\x1B[0m"<<endl;
	cout<<"  "<< "\x1B[31m" << "|,-'--|--'-.|           ░    ░ ░   ░  ░ ░        ░  ░" << "\x1B[0m"<<"    "<< "\x1B[33m" << "~~-___      )"<<"\x1B[0m"<<endl;
	cout<<"     "<< "\x1B[33m" << " /|||\\ "<< "\x1B[33m" << "            "<< "\x1B[31m" << "                    ░          " << "\x1B[0m"<<"          "<< "\x1B[33m" << "\\      \\"<< "\x1B[0m"<<endl;
    //cout<<" "<< "\x1B[43m" << "Hello World!\n" << "\x1B[0m"<<""<<endl;
}
