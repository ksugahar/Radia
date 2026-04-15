/*-------------------------------------------------------------------------
*
* File name:      rad_arc_current.h
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 rectangular cross-section arc with azimuthal current
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#ifndef __RAD_ARC_CURRENT_H
#define __RAD_ARC_CURRENT_H

#include "rad_geometry_3d.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTArcCur : public radTg3d {
public:
	//TVector3d CentrPoint; //moved to radTg3d
	TVector3d CircleCentrPoint;
	double R_min, R_max;
	double Phi_min, Phi_max;
	double Height;
	double J_azim;
	int NumberOfSectors;

	short BasedOnPrecLevel;

	short InternalFacesAfterCut;
	char J_IsNotZero;

	radTArcCur(const TVector3d& InCPoiVect, const double* InRadiiPtr, 
			   const double* InAnglesPtr, double InHeight, 
			   double InJ_azim, int InNumberOfSectors =0, short InBasedOnPrecLevel =0)
	{
		CircleCentrPoint = InCPoiVect;
		R_min = *InRadiiPtr; R_max = *(InRadiiPtr+1);
		Phi_min = *InAnglesPtr; Phi_max = *(InAnglesPtr+1);
		Height = InHeight;
		J_azim = InJ_azim;
		NumberOfSectors = InNumberOfSectors;
		BasedOnPrecLevel = InBasedOnPrecLevel;

		ComputeCentrPoint();

		if(J_azim != 0.) J_IsNotZero = 1;
		InternalFacesAfterCut = 0;
		ConsiderOnlyWithTrans = 0;
	}
	radTArcCur()
	{ 
		J_IsNotZero = 0; InternalFacesAfterCut = 0; ConsiderOnlyWithTrans = 0;
	}

	int Type_g3d() { return 3;}

	void ComputeCentrPoint()
	{
		double Phic = 0.5*(Phi_min + Phi_max);
		double rc = (2./3.)*(R_max*R_max + R_max*R_min + R_min*R_min)/(R_max + R_min);
		CentrPoint = TVector3d(rc*cos(Phic) + CircleCentrPoint.x, rc*sin(Phic) + CircleCentrPoint.y, CircleCentrPoint.z);
	}

	void B_comp(radTField* FieldPtr)
	{
		if(FieldPtr->FieldKey.J_) J_comp(FieldPtr);
		if(!(FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_ || FieldPtr->FieldKey.A_ || FieldPtr->FieldKey.Phi_ || FieldPtr->FieldKey.PreRelax_)) return;

		// Use analytical elliptic integral method
		B_compElliptic(FieldPtr);
	}
	void B_compElliptic(radTField*);  // Analytical method using elliptic integrals

	void J_comp(radTField* FieldPtr)
	{
		const double twoPi = 2.*3.141592653589793238;
		const int maxNumTwoPi = 10;

		TVector3d vR = FieldPtr->P - CircleCentrPoint;
		double r = vR.Abs();
		if((r < R_min) || (r > R_max)) return; //point is outside
		if(::fabs(vR.z) > 0.5*Height) return; //point is outside

		double phi = Argument(vR.x, vR.y); // -Pi < phi <= Pi
		if(phi < Phi_min)
		{
			for(int i=0; i<maxNumTwoPi; i++)
			{
				phi += twoPi; if(phi >= Phi_min) break;
			}
		}
		else if(phi > Phi_max)
		{
			for(int i=0; i<maxNumTwoPi; i++)
			{
				phi -= twoPi; if(phi <= Phi_max) break;
			}
		}
		if((phi < Phi_min) || (phi > Phi_max)) return; //point is outside

		TVector3d locJ(-J_azim*sin(phi), J_azim*cos(phi), 0);
		FieldPtr->J += locJ;
	}

	// Note: B_intComp is now handled through numerical integration for line integrals.
	// The analytical elliptic integral method (B_compElliptic) is used for point field computation.

	double Volume() { return 0.5*(Phi_max - Phi_min)*(R_max*R_max - R_min*R_min)*Height;}
	void VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique);

	// SimpleEnergyComp REMOVED (Phase C, 2026-04-16, energy-based API gone)

	void Push_backCenterPointAndField(radTFieldKey*, radTVectPairOfVect3d*, radTrans* pBaseTrans, radTg3d*, radTApplication*);

	// Dump / DumpPureObjInfo / DumpBin REMOVED (Phase B2b/B2c, 2026-04-15)


	// SubdivideItself* / CutItself / FindLowestAndUppestVertices REMOVED (Phase C, 2026-04-16)
	int DuplicateItself(radThg& hg, radTApplication*, char)
	{
		return FinishDuplication(new radTArcCur(*this), hg);
	}
	int SizeOfThis() { return sizeof(radTArcCur);}

	int ItemIsNotFullyInternalAfterCut()
	{
		return (InternalFacesAfterCut == 63)? 0 : 1;
	}
	void MapFaceAsInternalAfterCut(short FaceNo)
	{
		short ExtraFaceCode;
		switch(FaceNo)
		{
			case 1:
				ExtraFaceCode = 1;
				break;
			case 2:
				ExtraFaceCode = 2;
				break;
			case 3:
				ExtraFaceCode = 4;
				break;
			case 4:
				ExtraFaceCode = 8;
				break;
			case 5:
				ExtraFaceCode = 16;
				break;
			case 6:
				ExtraFaceCode = 32;
				break;
		}
		InternalFacesAfterCut |= ExtraFaceCode;
	}
	void ListFacesInternalAfterCut(short* FacesState)
	{
		short BufNum = InternalFacesAfterCut;
		for(int k=0; k<6; k++) { *(FacesState++) = BufNum & 1; BufNum >>= 1;}
	}
	void SetFacesInternalAfterCut(short* FacesState)
	{
		for(int k=0; k<6; k++) if(*(FacesState++)) MapFaceAsInternalAfterCut(k+1);
	}

	int ScaleCurrent(double scaleCoef) //virtual in g3d
	{
		J_azim *= scaleCoef; 
		return 1;
	}
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTBackgroundFieldSource : public radTg3d {
public:
	TVector3d BackgrB;  // Background magnetic flux density B in Tesla

	radTBackgroundFieldSource(const TVector3d& InBackgrB)
	{
		BackgrB = InBackgrB;  // Input in Tesla
	}

	radTBackgroundFieldSource()
	{
		BackgrB.Zero();
	}

	void B_comp(radTField* FieldPtr)
	{
		// BackgrB is magnetic flux density B in Tesla
		if(FieldPtr->FieldKey.B_) FieldPtr->B += BackgrB;
		// Convert B to H using H = B/μ₀ where μ₀ = 4π×10⁻⁷ H/m = 1.25663706212e-6 T/(A/m)
		if(FieldPtr->FieldKey.H_) FieldPtr->H += BackgrB * 795774.715459;  // 1/μ₀ = 795774.715459 (A/m)/T
		if(FieldPtr->FieldKey.A_) 
		{
			TVector3d BufA(FieldPtr->P.z*BackgrB.y, FieldPtr->P.x*BackgrB.z, FieldPtr->P.y*BackgrB.x);
			FieldPtr->A += BufA;
		}
	}
	void B_intComp(radTField* FieldPtr)
	{
		if(FieldPtr->FieldKey.FinInt_)
		{
			TVector3d D = FieldPtr->NextP - FieldPtr->P;
			double pathLength = sqrt(D.x*D.x + D.y*D.y + D.z*D.z);
			// BackgrB is in Tesla, convert to appropriate units for integrals
			TVector3d BufIb = pathLength * BackgrB;  // Ib: Integral of B
			if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufIb;
			// Ih: Integral of H = Integral of (B/μ₀)
			TVector3d BufIh = pathLength * BackgrB * 795774.715459;  // 1/μ₀
			if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufIh;
		}
		// Infinite integral is set to zero (though, formally, its contribution is infinite)
	}

	// Dump / DumpBin REMOVED (Phase B2b/B2c, 2026-04-15)

	int DuplicateItself(radThg& hg, radTApplication*, char)
	{
		return FinishDuplication(new radTBackgroundFieldSource(*this), hg);
	}
	int SizeOfThis() { return sizeof(radTBackgroundFieldSource);}
};

//-------------------------------------------------------------------------

#endif
