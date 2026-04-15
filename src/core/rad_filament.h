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

	void SimpleEnergyComp(radTField* FieldPtr)
	{
		// Energy: E = -integral(I dl . A) = -I * (EndPoint - StartPoint) . A
		// Radia now uses SI units (meters) internally, matching ELF.
		// I in Amperes, dl in meters, A in T*m -> Energy in Joules
		const double ConstForJ = -1.0;
		radTFieldKey LocFieldKey; LocFieldKey.A_ = 1;
		radTField LocField(LocFieldKey, FieldPtr->CompCriterium);
		LocField.P = 0.5*(StartPoint + EndPoint);
		((radTg3d*)(FieldPtr->HandleEnergyForceTorqueCompData.rep->hSource.rep))->B_genComp(&LocField);
		FieldPtr->Energy += (ConstForJ*I)*((EndPoint - StartPoint)*LocField.A);
	}

	// Dump / DumpBin REMOVED (Phase B2b/B2c, 2026-04-15)

	int SubdivideItself(double*, radThg&, radTApplication*, radTSubdivOptions*);
	int SubdivideItselfByOneSetOfParPlanes(TVector3d&, TVector3d*, int, radThg&, radTApplication*, radTSubdivOptions* pSubdivOptions, radTvhg*) 
	{
		radTSend Send;
		if(!pSubdivOptions->SubdivideCoils) return 1;
		else { Send.ErrorMessage("Radia::Error109"); return 0;}
	}
	int SubdivideItselfByParPlanes(double*, int, radThg&, radTApplication*, radTSubdivOptions* pSubdivOptions) 
	{
		radTSend Send;
		if(!pSubdivOptions->SubdivideCoils) return 1;
		else { Send.ErrorMessage("Radia::Error109"); return 0;}
	}
	int CutItself(TVector3d* InCuttingPlane, radThg& In_hg, radTPair_int_hg&, radTPair_int_hg&, radTApplication* radPtr, radTSubdivOptions* pSubdivOptions)
	{
		radTSend Send;
		if(!pSubdivOptions->SubdivideCoils) return 1;
		else { Send.ErrorMessage("Radia::Error109"); return 0;}
	}
	int SubdivideItselfByEllipticCylinder(double*, radTCylindricSubdivSpec*, radThg&, radTApplication*, radTSubdivOptions* pSubdivOptions)
	{ 
		radTSend Send;
		if(!pSubdivOptions->SubdivideCoils) return 1;
		else { Send.ErrorMessage("Radia::Error109"); return 0;}
	}

	int FindLowestAndUppestVertices(TVector3d&, radTSubdivOptions* pSubdivOptions, TVector3d&, TVector3d&, radTrans&, char&, char& Ignore)
	{
		Ignore = 1;
		if(!pSubdivOptions->SubdivideCoils) return 1;
		else 
		{ 
			radTSend Send;
			Send.ErrorMessage("Radia::Error109"); return 0;
		}
	}

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
