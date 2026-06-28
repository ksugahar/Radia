/*-------------------------------------------------------------------------
*
* File name:      rad_rectangular_block.h
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 rectangular parallelepiped with constant magnetization 
*                 or currect density
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#ifndef __RAD_RECTANGULAR_BLOCK_H
#define __RAD_RECTANGULAR_BLOCK_H

#include "rad_geometry_3d.h"

#include "rad_polyhedron.h"
#include "rad_convergence.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

extern radTConvergRepair& radCR;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct radTParallelepSurfIntData {
	TVector3d PointOnSurface;
	int SurfBoundInd;
// Surface bound indicator: 1 - lower, 2 - upper, 3 - left, 4 - right, 5 - back, 6 - front, 0 - No
	int IntegrandLen;
	void (radTg3d::*IntegrandFunPtr)(radTField*);
	radTField Field;

	double* InnerAbsPrecAndLimitsArray;
	short* InnerElemCompNotFinished;
	TVector3d** InnerIntegVal;
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTRecCur : public radTg3dRelax {
	radTParallelepSurfIntData* SurfIntDataPtr;
public:
	TVector3d Dimensions;
	TVector3d J;
	short J_IsNotZero;
	short InternalFacesAfterCut;

	radTRecCur(const TVector3d& InCPoiVect, const TVector3d& InDimsVect, 
			   const TVector3d& InMagnVect, const TVector3d& InJ_vect, const radThg& InMaterHandle, short InJ_IsNotZero =0) 
			   : radTg3dRelax(InCPoiVect, InMagnVect, InMaterHandle)
	{
		Dimensions=InDimsVect; J=InJ_vect;
		if(InMaterHandle.rep != 0) J_IsNotZero = 0;
		InternalFacesAfterCut = 0;

		J_IsNotZero = InJ_IsNotZero;
	}
	radTRecCur(const TVector3d& InCPoiVect, const TVector3d& InDimsVect, 
			   const TVector3d& InMagnVect, 
			   const TVector3d& InJ_vect, short InJ_IsNotZero)
			   : radTg3dRelax(InMagnVect)
	{
		CentrPoint=InCPoiVect; Dimensions=InDimsVect; //Magn=InMagnVect; 
		J=InJ_vect; J_IsNotZero = InJ_IsNotZero;
		InternalFacesAfterCut = 0;
	}
	// radTRecCur(CAuxBinStrVect&, ...) REMOVED (Phase B2c, 2026-04-15)
	radTRecCur() : radTg3dRelax()
	{ 
		InternalFacesAfterCut = 0;
	}

	int Type_g3dRelax() { return 1;}
	virtual int Type_RecCur() { return 0;}

	void B_comp(radTField*);
	void B_compMultipole(radTField*, double*);
	void B_intComp(radTField*);
	void B_intUtilSpecCaseZeroVxVy(const TVector3d&, const TVector3d&, short, TMatrix3d&, TVector3d&);

	void IntOverShape(radTField* FieldPtr) 
	{
		if(FieldPtr->ShapeIntDataPtr->IntOverSurf_) IntOverSurf(FieldPtr);
		else if(FieldPtr->ShapeIntDataPtr->IntOverVol_) IntOverVol(FieldPtr);
	}
	void IntOverSurf(radTField*);
	void FunForOuterIntAtSurfInt(double, TVector3d*);
	inline void FunForInnerIntAtSurfInt(double, TVector3d*);
	void IntOverVol(radTField*) {}

	double Volume() { return Dimensions.x*Dimensions.y*Dimensions.z;}
	void VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique);

	// SimpleEnergyComp REMOVED (Phase C, 2026-04-16, energy-based API gone)
	void UniformlyDistrPoints(double* q, TVector3d& P)
	{// This is not used
		P.x = CentrPoint.x + Dimensions.x*((*(q++))-0.5); 
		P.y = CentrPoint.y + Dimensions.y*((*(q++))-0.5); 
		P.z = CentrPoint.z + Dimensions.z*((*q)-0.5);
	}

	void Push_backCenterPointAndField(radTFieldKey*, radTVectPairOfVect3d*, radTrans*, radTg3d*, radTApplication*);
	
	// Dump / DumpPureObjInfo / DumpBin / DumpBin_RecMag / DumpBinParse_RecMag REMOVED (Phase B2b/B2c, 2026-04-15)


	int DuplicateItself(radThg& hg, radTApplication*, char) 
	{
		return FinishDuplication(new radTRecCur(*this), hg);
	}

	// SubdivideItself* / CutItself / FindLowestAndUppestVertices REMOVED (Phase C, 2026-04-16)

	int SetMaterial(radThg& InMatHandle, radTApplication* ApPtr) 
	{ 
		if(!J_IsNotZero) return radTg3dRelax::SetMaterial(InMatHandle, ApPtr);
		else return 1;
	}

	int ScaleCurrent(double scaleCoef) //virtual in g3d
	{//note: if(scaleCoef == 0) this still doesn't change J_IsNotZero
		if(J_IsNotZero) 
		{
			J *= scaleCoef; return 1;
		}
		else return 0;
	}

	int NumberOfDegOfFreedom() { return (MaterHandle.rep == 0)? 0 : 3;}
	int SizeOfThis() { return sizeof(radTRecCur);}

	// Current block is never converted to a polyhedron (the magnet path builds the MMMM polyhedron
	// directly in SetRecMag). Return 1 ("handled", like radTPolyhedron) so radTGroup::ConvertToPolyhedron
	// does not treat a current block as a conversion failure (base radTg3dRelax returns 0).
	int ConvertToPolyhedron(radThg&, radTApplication*, char) { return 1;}
	// CheckVertexPtsPositionsWithRespectToPlane REMOVED (Phase C, 2026-04-16)

	void DefineRelAndAbsTol(double* RelAbsTol)
	{
		double RelZeroToler = 1.E-09;
		RelZeroToler = 500.*((RelZeroToler>radCR.RelRand)? RelZeroToler : radCR.RelRand);

		TVector3d VectToCenter = 0.5*Dimensions;
		RelAbsTol[1] = RelZeroToler*NormAbs(VectToCenter);
		RelAbsTol[0] = RelZeroToler;
	}

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
	void MapFaceAsExternal(short FaceNo)
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
		InternalFacesAfterCut &= (!ExtraFaceCode);
	}
	void ListFacesInternalAfterCut(short* FacesState)
	{
		short BufNum = InternalFacesAfterCut;
		for(int k=0; k<6; k++) { *(FacesState++) = BufNum & 1; BufNum >>= 1;}
	}
	void SetFacesInternalAfterCut(short* FacesState)
	{
		for(int k=0; k<6; k++) 
		{
			if(*(FacesState++)) MapFaceAsInternalAfterCut(k+1);
		}
	}
};

//-------------------------------------------------------------------------

inline void radTRecCur::FunForInnerIntAtSurfInt(double Arg, TVector3d* VectArray)
{
	if(SurfIntDataPtr->SurfBoundInd==1 || SurfIntDataPtr->SurfBoundInd==2 || 
	   SurfIntDataPtr->SurfBoundInd==3 || SurfIntDataPtr->SurfBoundInd==4)
		SurfIntDataPtr->PointOnSurface.x = Arg;
	else if(SurfIntDataPtr->SurfBoundInd==5 || SurfIntDataPtr->SurfBoundInd==6)
		SurfIntDataPtr->PointOnSurface.y = Arg;

	SurfIntDataPtr->Field.P = SurfIntDataPtr->PointOnSurface;
	(((radTg3d*)this)->*(SurfIntDataPtr->IntegrandFunPtr))(&(SurfIntDataPtr->Field));

	for(int i=0; i<SurfIntDataPtr->IntegrandLen; i++) 
		VectArray[i] = (SurfIntDataPtr->Field.ShapeIntDataPtr->VectArray)[i];
}

//-------------------------------------------------------------------------

#endif