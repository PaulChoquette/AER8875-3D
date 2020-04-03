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
void Metric_c::ComputeMetric( Reader_c& read)
{
	cout << "Metric Calculation \tSTARTING...";
	//cout << "-------------------------------------- COMPUTE Normales & Surfaces --------------------------------------"<<endl;
  	Norm_Area(read);
	//cout << "-------------------------------------- COMPUTE Volumes --------------------------------------"<<endl;
	Volume(read);
	//cout << "-------------------------------------- COMPUTE elem2deltaSxyz --------------------------------------"<<endl;
	Elem2DeltaS(read);
	//cout << "-------------------------------------- COMPUTE Face2Center --------------------------------------"<<endl;
	Face2Center(read);
  /*
  Metric_c.Face2Center(Nelem_real, Nface);
	*/
	cout << "...............DONE" << endl;
}
////////////////////////////////////////////////////////////////////////////////
void Metric_c::Norm_Area(Reader_c& read)
{
	face2norm = new double*[nface];
	face2area = new double[nface];
	for(int iface=0; iface<nface; iface++)
	{
		face2norm[iface] = new double[read.ndime];
		int nnode = face2nnofa[iface];

		if(nnode == 3)
		{
			int pt1 = face2node[iface][0];
			int pt2 = face2node[iface][1];
			int pt3 = face2node[iface][2];

			double delta_xy_a = (read.coord[pt1][0] - read.coord[pt2][0])*(read.coord[pt1][1] + read.coord[pt2][1]);
			double delta_xy_b = (read.coord[pt2][0] - read.coord[pt3][0])*(read.coord[pt2][1] + read.coord[pt3][1]);
			double delta_xy_c = (read.coord[pt3][0] - read.coord[pt1][0])*(read.coord[pt3][1] + read.coord[pt1][1]);
			double delta_yz_a = (read.coord[pt1][1] - read.coord[pt2][1])*(read.coord[pt1][2] + read.coord[pt2][2]);
			double delta_yz_b = (read.coord[pt2][1] - read.coord[pt3][1])*(read.coord[pt2][2] + read.coord[pt3][2]);
			double delta_yz_c = (read.coord[pt3][1] - read.coord[pt1][1])*(read.coord[pt3][2] + read.coord[pt1][2]);
			double delta_zx_a = (read.coord[pt1][2] - read.coord[pt2][2])*(read.coord[pt1][0] + read.coord[pt2][0]);
			double delta_zx_b = (read.coord[pt2][2] - read.coord[pt3][2])*(read.coord[pt2][0] + read.coord[pt3][0]);
			double delta_zx_c = (read.coord[pt3][2] - read.coord[pt1][2])*(read.coord[pt3][0] + read.coord[pt1][0]);
			double S_x = 0.5*(delta_yz_a + delta_yz_b + delta_yz_c);
			double S_y = 0.5*(delta_zx_a + delta_zx_b + delta_zx_c);
			double S_z = 0.5*(delta_xy_a + delta_xy_b + delta_xy_c);

			face2area[iface] = sqrt(S_x*S_x + S_y*S_y + S_z*S_z);
			face2norm[iface][0] = S_x/face2area[iface];
			face2norm[iface][1] = S_y/face2area[iface];
			face2norm[iface][2] = S_z/face2area[iface];
		}
		else if(nnode == 4)
		{
			int pt1 = face2node[iface][0];
			int pt2 = face2node[iface][1];
			int pt3 = face2node[iface][2];
			int pt4 = face2node[iface][3];
			double delta_x_a = (read.coord[pt4][0] - read.coord[pt2][0]);
			double delta_x_b = (read.coord[pt3][0] - read.coord[pt1][0]);
			double delta_y_a = (read.coord[pt4][1] - read.coord[pt2][1]);
			double delta_y_b = (read.coord[pt3][1] - read.coord[pt1][1]);
			double delta_z_a = (read.coord[pt4][2] - read.coord[pt2][2]);
			double delta_z_b = (read.coord[pt3][2] - read.coord[pt1][2]);
			/* On met des (-) ci dessous mais on sait pas trop pourquoi.... ca regle le probleme */
			double S_x = -0.5*(delta_y_a*delta_z_b - delta_z_a*delta_y_b);
			double S_y = -0.5*(delta_z_a*delta_x_b - delta_x_a*delta_z_b);
			double S_z = -0.5*(delta_x_a*delta_y_b - delta_y_a*delta_x_b);
			face2area[iface] = sqrt(S_x*S_x + S_y*S_y + S_z*S_z);
			face2norm[iface][0] = S_x/face2area[iface];
			face2norm[iface][1] = S_y/face2area[iface];
			face2norm[iface][2] = S_z/face2area[iface];
			/*
			if(iface==22)
			{
				cout << "FACE ID" << endl;
				cout << "face2norm[iface][0] : "; cout << face2norm[iface][0] << endl;
				cout << "face2norm[iface][1] : "; cout << face2norm[iface][1] << endl;
				cout << "face2norm[iface][2] : "; cout << face2norm[iface][2] << endl;
			}
			*/
		}
	}
}
void Metric_c::Volume(Reader_c& read)
{


		elem2vol	= new double [nface];
		elem2vol    = new double [nelem];
		for(int ielem=0; ielem<nelem; ielem++)
		{

			int vtk = read.elem2vtk[ielem];
			//cout << "----------------- vtk : "; cout << vtk << endl;
	    	int Nbr_of_face = vtklookup[ndime-2][vtk][0];
			//cout << "Nbr_of_face : "; cout << Nbr_of_face << endl;
	    double TroisVolume = 0.0;
	    for(int iface=0; iface<Nbr_of_face; iface++)
	    {
			int faceID = elem2face[ielem][iface];
			//cout << "Face ID : "; cout << faceID <<endl;
			int nnode = face2nnofa[faceID];
			//cout << "nnode : "; cout << nnode << endl;
	      	double r_mid_x = 0.0;
	      	double r_mid_y = 0.0;
	      	double r_mid_z = 0.0;
	      	for(int inode=0; inode<nnode; inode++)
	      	{
				int nodeID = face2node[faceID][inode];
				//cout << "nodeID : "; cout << nodeID << endl;
				r_mid_x += read.coord[nodeID][0];
				r_mid_y += read.coord[nodeID][1];
				r_mid_z += read.coord[nodeID][2];
	        }
			r_mid_x = r_mid_x/nnode;
			r_mid_y = r_mid_y/nnode;
			r_mid_z = r_mid_z/nnode;
			/*
			cout << "r_mid_x : " ;cout << r_mid_x << endl;
			cout << "r_mid_y: " ;cout << r_mid_y << endl;
			cout << "r_mid_z : " ;cout << r_mid_z << endl;
			cout << "face2norm : "; cout << face2norm[faceID][0];cout << " ; "; cout << face2norm[faceID][1];cout << " ; "; cout << face2norm[faceID][2] << endl;
			cout << "face2area[faceID] : "; cout << face2area[faceID] << endl;
			*/
			double normale_x = face2norm[faceID][0];
			double normale_y = face2norm[faceID][1];
			double normale_z = face2norm[faceID][2];
			int elem1 = face2elem[faceID][0];
			int elem2 = face2elem[faceID][1];
			//cout << "elem1 : "; cout << elem1; cout << "  ;  elem2 : "; cout << elem2 << endl;
			if(ielem==elem1)
			{
				if(ielem>elem2)
				{
					//cout << "elem2 ---> elem1 &&& ielem==elem1" << endl;
					normale_x = -normale_x;
					normale_y = -normale_y;
					normale_z = -normale_z;
				}
			}
			else if(ielem==elem2)
			{
				if(ielem>elem1)
				{
					normale_x = -normale_x;
					normale_y = -normale_y;
					normale_z = -normale_z;
				}
			}
	      	TroisVolume += face2area[faceID]*(normale_x*r_mid_x + normale_y*r_mid_y + normale_z*r_mid_z);
	    }
	   	elem2vol [ielem] = TroisVolume/3.0;
	}
}
void Metric_c::SumNorm(Reader_c& read, int choix)
{
	double permitedError = 1.0e-12;
	string state = "NA";
	choix = 1;
	if(choix == 1)
	{
		double** Elem2sumNorm;
		Elem2sumNorm	= new double* [nelem];
		for(int ielem=0; ielem<nelem; ielem++)
		{
			Elem2sumNorm[ielem] = new double[read.ndime];
			Elem2sumNorm[ielem][0] = 0.0;
			Elem2sumNorm[ielem][1] = 0.0;
			Elem2sumNorm[ielem][2] = 0.0;
			int vtk = read.elem2vtk[ielem];
			int Nbr_of_face = vtklookup[ndime-2][vtk][0];
			for(int iface=0; iface<Nbr_of_face; iface++)
			{
				int faceID = elem2face[ielem][iface];
				double normale_x = face2area[faceID]*face2norm[faceID][0];
				double normale_y = face2area[faceID]*face2norm[faceID][1];
				double normale_z = face2area[faceID]*face2norm[faceID][2];
				int elem1 = face2elem[faceID][0];
				int elem2 = face2elem[faceID][1];
				if(ielem==elem1)
				{
					if(ielem>elem2)
					{
						normale_x = -normale_x;
						normale_y = -normale_y;
						normale_z = -normale_z;
					}
				}
				else if(ielem==elem2)
				{
					if(ielem>elem1)
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
				Elem2sumNorm[ielem][0] += normale_x;
				Elem2sumNorm[ielem][1] += normale_y;
				Elem2sumNorm[ielem][2] += normale_z;
			}
			if(Elem2sumNorm[ielem][0] > permitedError)
			{
				state = "bad";
			}
			else if(Elem2sumNorm[ielem][1] > permitedError)
			{
				state = "bad";
			}
			else if(Elem2sumNorm[ielem][2] > permitedError)
			{
				state = "bad";
			}
			/*
			cout << "------------------------"<<endl;
			cout << "Elem # : "; cout << ielem << endl;
			cout << "Elem2sumNorm[ielem][0] : "; cout <<  Elem2sumNorm[ielem][0] << endl;
			cout << "Elem2sumNorm[ielem][1] : "; cout <<  Elem2sumNorm[ielem][1] << endl;
			cout << "Elem2sumNorm[ielem][2] : "; cout <<  Elem2sumNorm[ielem][2] << endl;
			*/
			/*
			cout << "------------------" <<endl;
			cout << "elem # "; cout << ielem << endl;
			cout << "Nbr_of_face "; cout << Nbr_of_face << endl;
			cout << Elem2sumNorm[ielem][0] << endl;
			cout << Elem2sumNorm[ielem][1] << endl;
			cout << Elem2sumNorm[ielem][2] << endl;
			*/
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

	for(int iface=0; iface<nface; iface++)
	{
		SumNormLocalZone[0] += face2norm[iface][0];
		SumNormLocalZone[1] += face2norm[iface][1];
		SumNormLocalZone[2] += face2norm[iface][2];
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
	cout << "Verificaion de la somme des normales pour chaques elements : ";
	cout << state << endl;
}
void Metric_c::Elem2DeltaS(Reader_c& read)
{
	double Gamma = 1.4;
	elem2deltaSxyz	= new double* [nelem];
	for(int ielem=0; ielem<nelem; ielem++)
	{
		double somme_sX = 0.0;
		double somme_sY = 0.0;
		double somme_sZ = 0.0;
		elem2deltaSxyz[ielem] = new double[read.ndime];
		int vtk = read.elem2vtk[ielem];
		int Nbr_of_face = vtklookup[ndime-2][vtk][0];
		for (int iface=0; iface<Nbr_of_face; iface++)
	{
		int faceID = elem2face[ielem][iface];
		somme_sX += fabs(face2norm[faceID][0]*face2area[faceID]);
		somme_sY += fabs(face2norm[faceID][1]*face2area[faceID]);
		somme_sZ += fabs(face2norm[faceID][2]*face2area[faceID]);
	}
		elem2deltaSxyz[ielem][0] = 0.5*somme_sX;
	elem2deltaSxyz[ielem][1] = 0.5*somme_sY;
	elem2deltaSxyz[ielem][2] = 0.5*somme_sZ;
	}

}
void Metric_c::Face2Center(Reader_c& read)
{
	elem2center	= new double* [nelem];
	face2elemCenter	= new double**[nface];
	face2midpoint = new double* [nface];

	for(int face_i=0; face_i<nface; face_i++)
	{
		face2midpoint[face_i] = new double[ndime];
		face2elemCenter[face_i] = new double*[2];
		for (int ix=0;ix<2;++ix) {
			face2elemCenter[face_i][ix] = new double[ndime];
		}
	}
	for(int ielem=0; ielem<nelem; ielem++)
	{
		elem2center[ielem] = new double[read.ndime];
		//cout << "Elemt ID : "; cout << ielem << endl;
		int vtk = read.elem2vtk[ielem];
		int Nbr_of_face = vtklookup[ndime-2][vtk][0];
		int Nbr_of_node = vtklookup[ndime-2][vtk][1];
		//cout << "Nbr_of_face : "; cout << Nbr_of_face << endl;
		//cout << "Nbr_of_node : "; cout << Nbr_of_node << endl;
		double c_x = 0.0;
		double c_y = 0.0;
		double c_z = 0.0;
		for(int inode=0; inode<Nbr_of_node; inode++)
		{
			int inodeD = read.elem2node[ielem][inode];
			c_x += read.coord[inodeD][0];
			c_y += read.coord[inodeD][1];
			c_z += read.coord[inodeD][2];
		}
		elem2center[ielem][0] = c_x/Nbr_of_node;
		elem2center[ielem][1] = c_y/Nbr_of_node;
		elem2center[ielem][2] = c_z/Nbr_of_node;
		//cout << "Center_x : "; cout << elem2center[ielem][0] << endl;
		//cout << "Center_y : "; cout << elem2center[ielem][1] << endl;
		//cout << "Center_z : "; cout << elem2center[ielem][2] << endl;
		for(int iface=0; iface<Nbr_of_face;iface++)
		{
			int faceID = elem2face[ielem][iface];
				//cout << "faceID : "; cout << faceID << endl;
				int nnode = face2nnofa[faceID];
				//cout << "nnode : "; cout << nnode << endl;
			double r_mid_x = 0.0;
			double r_mid_y = 0.0;
			double r_mid_z = 0.0;
			for(int inode=0; inode<nnode; inode++)
			{
				r_mid_x += read.coord[face2node[faceID][inode]][0];
				r_mid_y += read.coord[face2node[faceID][inode]][1];
				r_mid_z += read.coord[face2node[faceID][inode]][2];
			}
			r_mid_x = r_mid_x/nnode;
			r_mid_y = r_mid_y/nnode;
			r_mid_z = r_mid_z/nnode;

			face2midpoint[faceID][0] = r_mid_x;
			face2midpoint[faceID][1] = r_mid_y;
			face2midpoint[faceID][2] = r_mid_z;
			/*
			cout << "Pour face ID : "; cout << faceID << endl;
			cout << "Center_x : "; cout << r_mid_x << endl;
			cout << "Center_y : "; cout << r_mid_y << endl;
			cout << "Center_z : "; cout << r_mid_z << endl;
			*/
			int elem1 = face2elem[faceID][0];
			int elem2 = face2elem[faceID][1];
			/*
			cout << "ielem : "; cout << ielem <<endl;
			cout << "faceID : "; cout << faceID <<endl;
			*/
			double delta_x = r_mid_x-elem2center[ielem][0];
			double delta_y = r_mid_y-elem2center[ielem][1];
			double delta_z = r_mid_z-elem2center[ielem][2];
			if(ielem == elem1)
			{
				face2elemCenter[faceID][0][0] = delta_x;
				face2elemCenter[faceID][0][1] = delta_y;
				face2elemCenter[faceID][0][2] = delta_z;
			}
			else if(ielem == elem2)
			{
				face2elemCenter[faceID][1][0] = delta_x;
				face2elemCenter[faceID][1][1] = delta_y;
				face2elemCenter[faceID][1][2] = delta_z;
			}
			else
			{
				cout << "ERROR" << endl;
			}
		}
	}
}
