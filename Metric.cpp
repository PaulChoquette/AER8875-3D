#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <array>
#include <string>
#include <fstream>
#include <cstdio>
#include <time.h>
#include <omp.h>
#include <iomanip>
#include <typeinfo>
#include <functional> // std::divides
#include "Metric.h"
#include "Reader.h"
#include "Solver.h"
#include "main.h"
using namespace std;
////////////////////////////////////////////////////////////////////////////////
void Metric_c::Compute(Solver_c& solve, Reader_c& read)
{
	cout << "-------------------------------------- COMPUTE Normales & Surfaces --------------------------------------"<<endl;
  Norm_Area(solve, read);
	cout << "-------------------------------------- COMPUTE Volumes --------------------------------------"<<endl;
	Volume(solve, read);
	cout << "-------------------------------------- COMPUTE Elem2DeltaS_xyz --------------------------------------"<<endl;
	Elem2DeltaS(solve, read);
	cout << "-------------------------------------- COMPUTE Face2Center --------------------------------------"<<endl;
	Face2Center(solve, read);
  /*
  Metric_c.Face2Center(Nelem_real, Nface);
	*/
}
////////////////////////////////////////////////////////////////////////////////
void Metric_c::Norm_Area(Solver_c& solve, Reader_c& read)
{
	//cout << "1" << endl;
	Face2Norm	= new double** [solve.nzone];
	Face2Area	= new double* [solve.nzone];
	//cout << "2" << endl;
	for(int i_zone=0; i_zone<solve.nzone; i_zone++)
	{
		//cout << "3" << endl;cout << "i_zone : ";cout << i_zone << endl;
		//cout << "zone2nface : ";cout << solve.zone2nface[i_zone] << endl;
		Face2Norm[i_zone] = new double* [solve.zone2nface[i_zone]];
		Face2Area[i_zone] = new double [solve.zone2nface[i_zone]];
		//cout << "4" << endl;
		for(int i_face=0; i_face<solve.zone2nface[i_zone]; i_face++)
		{
			Face2Norm[i_zone][i_face] = new double[read.ndime];
			//cout << "5" << endl;
			int nnode = solve.face2Nbr_of_node[i_zone][i_face];
			//cout << "nnode : ";cout << nnode << endl;
      if(nnode == 3)
      {
        int pt1 = solve.face2node[i_zone][i_face][0];
        int pt2 = solve.face2node[i_zone][i_face][1];
        int pt3 = solve.face2node[i_zone][i_face][2];
				/*
				cout << "---- FACE triangle # : "; cout << i_face << endl;
				cout << "pt1 : "; cout << pt1 << endl;
				cout << "pt1_x : "; cout << solve.zone2coord[i_zone][pt1][0];
				cout << " ; pt1_y : "; cout << solve.zone2coord[i_zone][pt1][1];
				cout << " ; pt1_z : "; cout << solve.zone2coord[i_zone][pt1][2] << endl;
				cout << "pt2 : "; cout << pt2 << endl;
				cout << "pt2_x : "; cout << solve.zone2coord[i_zone][pt2][0];
				cout << " ; pt2_y : "; cout << solve.zone2coord[i_zone][pt2][1];
				cout << " ; pt2_z : "; cout << solve.zone2coord[i_zone][pt2][2] << endl;
				cout << "pt3 : "; cout << pt3 << endl;
				cout << "pt3_x : "; cout << solve.zone2coord[i_zone][pt3][0];
				cout << " ; pt3_y : "; cout << solve.zone2coord[i_zone][pt3][1];
				cout << " ; pt3_z : "; cout << solve.zone2coord[i_zone][pt3][2] << endl;
				*/
        double delta_xy_a = (solve.zone2coord[i_zone][pt1][0] - solve.zone2coord[i_zone][pt2][0])*(solve.zone2coord[i_zone][pt1][1] + solve.zone2coord[i_zone][pt2][1]);
        double delta_xy_b = (solve.zone2coord[i_zone][pt2][0] - solve.zone2coord[i_zone][pt3][0])*(solve.zone2coord[i_zone][pt2][1] + solve.zone2coord[i_zone][pt3][1]);
        double delta_xy_c = (solve.zone2coord[i_zone][pt3][0] - solve.zone2coord[i_zone][pt1][0])*(solve.zone2coord[i_zone][pt3][1] + solve.zone2coord[i_zone][pt1][1]);
        double delta_yz_a = (solve.zone2coord[i_zone][pt1][1] - solve.zone2coord[i_zone][pt2][1])*(solve.zone2coord[i_zone][pt1][2] + solve.zone2coord[i_zone][pt2][2]);
        double delta_yz_b = (solve.zone2coord[i_zone][pt2][1] - solve.zone2coord[i_zone][pt3][1])*(solve.zone2coord[i_zone][pt2][2] + solve.zone2coord[i_zone][pt3][2]);
        double delta_yz_c = (solve.zone2coord[i_zone][pt3][1] - solve.zone2coord[i_zone][pt1][1])*(solve.zone2coord[i_zone][pt3][2] + solve.zone2coord[i_zone][pt1][2]);
        double delta_zx_a = (solve.zone2coord[i_zone][pt1][2] - solve.zone2coord[i_zone][pt2][2])*(solve.zone2coord[i_zone][pt1][0] + solve.zone2coord[i_zone][pt2][0]);
        double delta_zx_b = (solve.zone2coord[i_zone][pt2][2] - solve.zone2coord[i_zone][pt3][2])*(solve.zone2coord[i_zone][pt2][0] + solve.zone2coord[i_zone][pt3][0]);
        double delta_zx_c = (solve.zone2coord[i_zone][pt3][2] - solve.zone2coord[i_zone][pt1][2])*(solve.zone2coord[i_zone][pt3][0] + solve.zone2coord[i_zone][pt1][0]);
        double S_x = 0.5*(delta_yz_a + delta_yz_b + delta_yz_c);
        double S_y = 0.5*(delta_zx_a + delta_zx_b + delta_zx_c);
        double S_z = 0.5*(delta_xy_a + delta_xy_b + delta_xy_c);
				/*
				cout << "S_x";cout << S_x  << endl;
				cout << "S_y";cout << S_y << endl;
				cout << "S_z";cout << S_z << endl;
				*/
        Face2Area[i_zone][i_face] = sqrt(S_x*S_x + S_y*S_y + S_z*S_z);
        //cout << "Area : ";cout << Face2Area[i_zone][i_face] << endl;
				Face2Norm[i_zone][i_face][0] = S_x/Face2Area[i_zone][i_face];
        Face2Norm[i_zone][i_face][1] = S_y/Face2Area[i_zone][i_face];
        Face2Norm[i_zone][i_face][2] = S_z/Face2Area[i_zone][i_face];
      }
      else if(nnode == 4)
      {
        int pt1 = solve.face2node[i_zone][i_face][0];
        int pt2 = solve.face2node[i_zone][i_face][1];
        int pt3 = solve.face2node[i_zone][i_face][2];
        int pt4 = solve.face2node[i_zone][i_face][3];
        double delta_x_a = (solve.zone2coord[i_zone][pt4][0] - solve.zone2coord[i_zone][pt2][0]);
        double delta_x_b = (solve.zone2coord[i_zone][pt3][0] - solve.zone2coord[i_zone][pt1][0]);
        double delta_y_a = (solve.zone2coord[i_zone][pt4][1] - solve.zone2coord[i_zone][pt2][1]);
        double delta_y_b = (solve.zone2coord[i_zone][pt3][1] - solve.zone2coord[i_zone][pt1][1]);
        double delta_z_a = (solve.zone2coord[i_zone][pt4][2] - solve.zone2coord[i_zone][pt2][2]);
        double delta_z_b = (solve.zone2coord[i_zone][pt3][2] - solve.zone2coord[i_zone][pt1][2]);
/* On met des (-) ci dessous mais on sait pas trop pourquoi.... ca regle le probleme */
        double S_x = -0.5*(delta_y_a*delta_z_b - delta_z_a*delta_y_b);
        double S_y = -0.5*(delta_z_a*delta_x_b - delta_x_a*delta_z_b);
        double S_z = -0.5*(delta_x_a*delta_y_b - delta_y_a*delta_x_b);
        Face2Area[i_zone][i_face] = sqrt(S_x*S_x + S_y*S_y + S_z*S_z);
        Face2Norm[i_zone][i_face][0] = S_x/Face2Area[i_zone][i_face];
        Face2Norm[i_zone][i_face][1] = S_y/Face2Area[i_zone][i_face];
        Face2Norm[i_zone][i_face][2] = S_z/Face2Area[i_zone][i_face];
				/*
				if(i_face==22)
				{
					cout << "FACE ID" << endl;
					cout << "Face2Norm[i_zone][i_face][0] : "; cout << Face2Norm[i_zone][i_face][0] << endl;
					cout << "Face2Norm[i_zone][i_face][1] : "; cout << Face2Norm[i_zone][i_face][1] << endl;
					cout << "Face2Norm[i_zone][i_face][2] : "; cout << Face2Norm[i_zone][i_face][2] << endl;
				}
				*/
      }
		}
	}
}
void Metric_c::Volume(Solver_c& solve, Reader_c& read)
{
	Elem2Vol	= new double* [solve.nzone];
	for(int i_zone=0; i_zone<solve.nzone; i_zone++)
	{
		Elem2Vol[i_zone]	= new double [solve.zone2nface[i_zone]];
		Elem2Vol[i_zone] = new double [solve.zone2nelem[i_zone]];
		//cout << "solve.zone2nelem : " ;cout << solve.zone2nelem[i_zone] << endl;
		for(int i_elem=0; i_elem<solve.zone2nelem[i_zone]; i_elem++)
		{
			int ElemIdentifier = solve.elem2vtk[i_zone][i_elem];
			//cout << "----------------- ElemIdentifier : "; cout << ElemIdentifier << endl;
	    int Nbr_of_face = vtklookup[0][ElemIdentifier][0];
			//cout << "Nbr_of_face : "; cout << Nbr_of_face << endl;
	    double TroisVolume = 0.0;
	    for(int i_face=0; i_face<Nbr_of_face; i_face++)
	    {
				int faceID = solve.elem2face[i_zone][i_elem][i_face];
				//cout << "Face ID : "; cout << faceID <<endl;
				int nnode = solve.face2Nbr_of_node[i_zone][faceID];
				//cout << "nnode : "; cout << nnode << endl;
	      double r_mid_x = 0.0;
	      double r_mid_y = 0.0;
	      double r_mid_z = 0.0;
	      for(int Node_i=0; Node_i<nnode; Node_i++)
	      {
					int nodeID = solve.face2node[i_zone][faceID][Node_i];
					//cout << "nodeID : "; cout << nodeID << endl;
	        r_mid_x += solve.zone2coord[i_zone][nodeID][0];
	        r_mid_y += solve.zone2coord[i_zone][nodeID][1];
	        r_mid_z += solve.zone2coord[i_zone][nodeID][2];
	      }
	      r_mid_x = r_mid_x/nnode;
	      r_mid_y = r_mid_y/nnode;
	      r_mid_z = r_mid_z/nnode;
				/*
				cout << "r_mid_x : " ;cout << r_mid_x << endl;
				cout << "r_mid_y: " ;cout << r_mid_y << endl;
				cout << "r_mid_z : " ;cout << r_mid_z << endl;
				cout << "Face2Norm : "; cout << Face2Norm[i_zone][faceID][0];cout << " ; "; cout << Face2Norm[i_zone][faceID][1];cout << " ; "; cout << Face2Norm[i_zone][faceID][2] << endl;
				cout << "Face2Area[i_zone][faceID] : "; cout << Face2Area[i_zone][faceID] << endl;
				*/
				double normale_x = Face2Norm[i_zone][faceID][0];
				double normale_y = Face2Norm[i_zone][faceID][1];
				double normale_z = Face2Norm[i_zone][faceID][2];
				int elem1 = solve.face2elem[i_zone][faceID][0];
				int elem2 = solve.face2elem[i_zone][faceID][1];
				//cout << "elem1 : "; cout << elem1; cout << "  ;  elem2 : "; cout << elem2 << endl;
				if(i_elem==elem1)
				{
					if(i_elem>elem2)
					{
						//cout << "elem2 ---> elem1 &&& i_elem==elem1" << endl;
						normale_x = -normale_x;
						normale_y = -normale_y;
						normale_z = -normale_z;
					}
				}
				else if(i_elem==elem2)
				{
					if(i_elem>elem1)
					{
						normale_x = -normale_x;
						normale_y = -normale_y;
						normale_z = -normale_z;
					}
				}
	      TroisVolume += Face2Area[i_zone][faceID]*(normale_x*r_mid_x + normale_y*r_mid_y + normale_z*r_mid_z);
	    }
	    Elem2Vol[i_zone][i_elem] = TroisVolume/3.0;
		}
	}
}

void Metric_c::SumNorm(Solver_c& solve, Reader_c& read, int choix)
{
	double permitedError = 1.0e-12;
	string state = "NA";
	choix = 1;
	if(choix == 1)
	{
		double*** Elem2sumNorm	= new double** [solve.nzone];
	  for(int i_zone=0; i_zone<solve.nzone; i_zone++)
		{
			Elem2sumNorm[i_zone]	= new double* [solve.zone2nelem[i_zone]];
			for(int i_elem=0; i_elem<solve.zone2nelem[i_zone]; i_elem++)
			{
				Elem2sumNorm[i_zone][i_elem] = new double[read.ndime];
				Elem2sumNorm[i_zone][i_elem][0] = 0.0;
				Elem2sumNorm[i_zone][i_elem][1] = 0.0;
				Elem2sumNorm[i_zone][i_elem][2] = 0.0;
				int ElemIdentifier = solve.elem2vtk[i_zone][i_elem];
		    int Nbr_of_face = vtklookup[0][ElemIdentifier][0];
		    for(int i_face=0; i_face<Nbr_of_face; i_face++)
		    {
					int faceID = solve.elem2face[i_zone][i_elem][i_face];
					double normale_x = Face2Area[i_zone][faceID]*Face2Norm[i_zone][faceID][0];
					double normale_y = Face2Area[i_zone][faceID]*Face2Norm[i_zone][faceID][1];
					double normale_z = Face2Area[i_zone][faceID]*Face2Norm[i_zone][faceID][2];
					int elem1 = solve.face2elem[i_zone][faceID][0];
					int elem2 = solve.face2elem[i_zone][faceID][1];
					if(i_elem==elem1)
					{
						if(i_elem>elem2)
						{
							normale_x = -normale_x;
							normale_y = -normale_y;
							normale_z = -normale_z;
						}
					}
					else if(i_elem==elem2)
					{
						if(i_elem>elem1)
						{
							normale_x = -normale_x;
							normale_y = -normale_y;
							normale_z = -normale_z;
						}
					}
					/*
					cout << "faceID "; cout << faceID << endl;
					cout << "normale_x "; cout << normale_x << endl;
					cout << "normale_y "; cout << normale_y << endl;
					cout << "normale_z "; cout << normale_z << endl;
					*/
					Elem2sumNorm[i_zone][i_elem][0] += normale_x;
					Elem2sumNorm[i_zone][i_elem][1] += normale_y;
					Elem2sumNorm[i_zone][i_elem][2] += normale_z;
				}
				if(Elem2sumNorm[i_zone][i_elem][0] > permitedError)
				{
					state = "bad";
				}
				else if(Elem2sumNorm[i_zone][i_elem][1] > permitedError)
				{
					state = "bad";
				}
				else if(Elem2sumNorm[i_zone][i_elem][2] > permitedError)
				{
					state = "bad";
				}
				/*
				cout << "------------------------"<<endl;
				cout << "Elem # : "; cout << i_elem << endl;
				cout << "Elem2sumNorm[i_zone][i_elem][0] : "; cout <<  Elem2sumNorm[i_zone][i_elem][0] << endl;
				cout << "Elem2sumNorm[i_zone][i_elem][1] : "; cout <<  Elem2sumNorm[i_zone][i_elem][1] << endl;
				cout << "Elem2sumNorm[i_zone][i_elem][2] : "; cout <<  Elem2sumNorm[i_zone][i_elem][2] << endl;
				*/
				/*
				cout << "------------------" <<endl;
				cout << "elem # "; cout << i_elem << endl;
				cout << "Nbr_of_face "; cout << Nbr_of_face << endl;
				cout << Elem2sumNorm[i_zone][i_elem][0] << endl;
				cout << Elem2sumNorm[i_zone][i_elem][1] << endl;
				cout << Elem2sumNorm[i_zone][i_elem][2] << endl;
				*/
			}
		}
		if(state == "NA")
		{
			state = "Sum=0 for every elements";
		}
	}
	else
	{
		double SumNormAllZones[3];
	  double SumNormLocalZone[3];
	  SumNormAllZones[0] = 0.0;
	  SumNormAllZones[1] = 0.0;
	  SumNormAllZones[2] = 0.0;
	  SumNormLocalZone[0] = 0.0;
	  SumNormLocalZone[1] = 0.0;
	  SumNormLocalZone[2] = 0.0;
	  for(int i_zone=0; i_zone<solve.nzone; i_zone++)
		{
	    for(int i_face=0; i_face<solve.zone2nface[i_zone]; i_face++)
	    {
	      SumNormLocalZone[0] += Face2Norm[i_zone][i_face][0];
	      SumNormLocalZone[1] += Face2Norm[i_zone][i_face][1];
	      SumNormLocalZone[2] += Face2Norm[i_zone][i_face][2];
	    }
	    SumNormAllZones[0] += SumNormLocalZone[0];
	    SumNormAllZones[1] += SumNormLocalZone[1];
	    SumNormAllZones[2] += SumNormLocalZone[2];
			/*
			cout << "SumNormAllZones[0] : "; cout <<  SumNormAllZones[0] << endl;
			cout << "SumNormAllZones[1] : "; cout <<  SumNormAllZones[1] << endl;
			cout << "SumNormAllZones[2] : "; cout <<  SumNormAllZones[2] << endl;
			*/
	  }
	}
	cout << "Verificaion de la somme des normales pour chaques elements : ";
	cout << state << endl;
}

void Metric_c::Elem2DeltaS(Solver_c& solve, Reader_c& read)
{
	double Gamma = 1.4;
	Elem2DeltaS_xyz	= new double** [solve.nzone];
	for(int i_zone=0; i_zone<solve.nzone; i_zone++)
	{
		Elem2DeltaS_xyz[i_zone]	= new double* [solve.zone2nelem[i_zone]];
		for(int i_elem=0; i_elem<solve.zone2nelem[i_zone]; i_elem++)
		{
			double somme_sX = 0.0;
		  double somme_sY = 0.0;
		  double somme_sZ = 0.0;
			Elem2DeltaS_xyz[i_zone][i_elem] = new double[read.ndime];
			int ElemIdentifier = solve.elem2vtk[i_zone][i_elem];
			int Nbr_of_face = vtklookup[0][ElemIdentifier][0];
			for (int i_face=0; i_face<Nbr_of_face; i_face++)
	    {
	        int faceID = solve.elem2face[i_zone][i_elem][i_face];
	        somme_sX += abs(Face2Norm[i_zone][faceID][0]*Face2Area[i_zone][faceID]);
	        somme_sY += abs(Face2Norm[i_zone][faceID][1]*Face2Area[i_zone][faceID]);
	        somme_sZ += abs(Face2Norm[i_zone][faceID][2]*Face2Area[i_zone][faceID]);
	    }
			Elem2DeltaS_xyz[i_zone][i_elem][0] = 0.5*somme_sX;
	    Elem2DeltaS_xyz[i_zone][i_elem][1] = 0.5*somme_sY;
	    Elem2DeltaS_xyz[i_zone][i_elem][2] = 0.5*somme_sZ;
		}
	}
}

void Metric_c::Face2Center(Solver_c& solve, Reader_c& read)
{
	Elem2Center	= new double** [solve.nzone];
	Face2ElemCenter	= new double** [solve.nzone];
	for(int i_zone=0; i_zone<solve.nzone; i_zone++)
	{
		Elem2Center[i_zone]	= new double* [solve.zone2nelem[i_zone]];
		Face2ElemCenter[i_zone]	= new double* [solve.zone2nface[i_zone]];
		for(int face_i=0; face_i<solve.zone2nface[i_zone]; face_i++)
		{
			Face2ElemCenter[i_zone][face_i] = new double[2];
		}
		for(int i_elem=0; i_elem<solve.zone2nelem[i_zone]; i_elem++)
		{
			Elem2Center[i_zone][i_elem] = new double[read.ndime];
			//cout << "Elemt ID : "; cout << i_elem << endl;
			int ElemIdentifier = solve.elem2vtk[i_zone][i_elem];
			int Nbr_of_face = vtklookup[0][ElemIdentifier][0];
			int Nbr_of_node = vtklookup[0][ElemIdentifier][1];
			//cout << "Nbr_of_face : "; cout << Nbr_of_face << endl;
			//cout << "Nbr_of_node : "; cout << Nbr_of_node << endl;
			double c_x = 0.0;
	    double c_y = 0.0;
	    double c_z = 0.0;
	    for(int Node_i=0; Node_i<Nbr_of_node; Node_i++)
	    {
	      int Node_ID = solve.elem2node[i_zone][i_elem][Node_i];
	      c_x += solve.zone2coord[i_zone][Node_ID][0];
	      c_y += solve.zone2coord[i_zone][Node_ID][1];
	      c_z += solve.zone2coord[i_zone][Node_ID][2];
	    }
	    Elem2Center[i_zone][i_elem][0] = c_x/Nbr_of_node;
	    Elem2Center[i_zone][i_elem][1] = c_y/Nbr_of_node;
	    Elem2Center[i_zone][i_elem][2] = c_z/Nbr_of_node;
			//cout << "Center_x : "; cout << Elem2Center[i_zone][i_elem][0] << endl;
			//cout << "Center_y : "; cout << Elem2Center[i_zone][i_elem][1] << endl;
			//cout << "Center_z : "; cout << Elem2Center[i_zone][i_elem][2] << endl;
	    for(int i_face=0; i_face<Nbr_of_face;i_face++)
	    {
	      int faceID = solve.elem2face[i_zone][i_elem][i_face];
				//cout << "faceID : "; cout << faceID << endl;
				int nnode = solve.face2Nbr_of_node[i_zone][faceID];
				//cout << "nnode : "; cout << nnode << endl;
	      double r_mid_x = 0.0;
	      double r_mid_y = 0.0;
	      double r_mid_z = 0.0;
	      for(int Node_i=0; Node_i<nnode; Node_i++)
	      {
	        r_mid_x += solve.zone2coord[i_zone][solve.face2node[i_zone][faceID][Node_i]][0];
	        r_mid_y += solve.zone2coord[i_zone][solve.face2node[i_zone][faceID][Node_i]][1];
	        r_mid_z += solve.zone2coord[i_zone][solve.face2node[i_zone][faceID][Node_i]][2];
	      }
	      r_mid_x = r_mid_x/nnode;
	      r_mid_y = r_mid_y/nnode;
	      r_mid_z = r_mid_z/nnode;
				/*
				cout << "Pour face ID : "; cout << faceID << endl;
				cout << "Center_x : "; cout << r_mid_x << endl;
				cout << "Center_y : "; cout << r_mid_y << endl;
				cout << "Center_z : "; cout << r_mid_z << endl;
				*/
				int elem1 = solve.face2elem[i_zone][faceID][0];
				int elem2 = solve.face2elem[i_zone][faceID][1];
				/*
				cout << "i_elem : "; cout << i_elem <<endl;
				cout << "faceID : "; cout << faceID <<endl;
				*/
				double delta_x = r_mid_x-Elem2Center[i_zone][i_elem][0];
				double delta_y = r_mid_y-Elem2Center[i_zone][i_elem][1];
				double delta_z = r_mid_z-Elem2Center[i_zone][i_elem][2];
	      if(i_elem == elem1)
	      {
					Face2ElemCenter[i_zone][faceID][0] = sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z);
	      }
	      else if(i_elem == elem2)
	      {
					Face2ElemCenter[i_zone][faceID][1] = sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z);
	      }
				else
				{
					cout << "ERROR" << endl;
				}
	    }
		}
  }
}
