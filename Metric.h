#ifndef METRIC_H // include guard
#define METRIC_H

#include <iostream>
#include "Connect.h"
class Solver_c;
class Metric_c : public Connect_c
{
public:

	// ================================================= INITIALIZATION ====================================================
	// (Mettre les variables a initialiser ici)

	// :
	double SumNormLocalZone;
	double SumNormAllZones;
	// :
	double** Face2Area;
	double** Elem2Vol;
	// :
	double*** Face2Norm;
	double*** Elem2Center;
	double*** Face2ElemCenter;
	double*** Elem2DeltaS_xyz;
	// =========================================== FUNCTION MEMBERS ============================================
	void Compute(Solver_c& solve, Reader_c& read);
	void Norm_Area(Solver_c& solve, Reader_c& read);
	void SumNorm(Solver_c& solve, Reader_c& read, int choix);
	void Volume(Solver_c& solve, Reader_c& read);
	void Elem2DeltaS(Solver_c& solve, Reader_c& read);
	void Face2Center(Solver_c& solve, Reader_c& read);
};

#endif
