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
	double* face2area;
	double* elem2vol;
	// :
	double** face2norm;
	double** elem2center;
	double** face2elemCenter;
	double** elem2deltaSxyz;
	// =========================================== FUNCTION MEMBERS ============================================
	void ComputeMetric(Reader_c& read);
	void Norm_Area(Reader_c& read);
	void SumNorm(Reader_c& read, int choix);
	void Volume(Reader_c& read);
	void Elem2DeltaS(Reader_c& read);
	void Face2Center(Reader_c& read);
};

#endif
