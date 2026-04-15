/*-------------------------------------------------------------------------
*
* File name:      rad_filament.h
*
* Project:        RADIA
*
* Description:    Magnetic field source: filament conductor
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#ifndef __RAD_FILAMENT_H
#define __RAD_FILAMENT_H

#include "rad_geometry_3d.h"
#include "rad_transform_def.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTFlmLinCur : public radTg3d {
	radTrans NativeRotation;
	double Length;

public:
	double I;
	TVector3d StartPoint, EndPoint;

	radTFlmLinCur(const TVector3d&, const TVector3d&, double);

	radTFlmLinCur() {}

	int Type_g3d() { return 4;}

	void SetNativeRotation(const TVector3d&, double);

	void B_comp(radTField*);
	void B_intComp(radTField*);

	// SimpleEnergyComp REMOVED (Phase C, 2026-04-16, energy-based API gone)

	// Dump / DumpBin REMOVED (Phase B2b/B2c, 2026-04-15)

	// SubdivideItself* / CutItself / FindLowestAndUppestVertices REMOVED (Phase C, 2026-04-16)

	int DuplicateItself(radThg& hg, radTApplication*, char)
	{ 
		return FinishDuplication(new radTFlmLinCur(*this), hg);
	}

	int NumberOfDegOfFreedom() { return 0;}
	int SizeOfThis() { return sizeof(radTFlmLinCur);}

	int ScaleCurrent(double scaleCoef) //virtual in g3d
	{
		I *= scaleCoef; 
		return 1;
	}
};

//-------------------------------------------------------------------------

#endif
