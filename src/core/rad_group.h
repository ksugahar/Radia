/*-------------------------------------------------------------------------
*
* File name:      radgroup.h
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 group (/container) of magnetic field sources
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#ifndef __RADGROUP_H
#define __RADGROUP_H

#include "rad_geometry_3d.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

using radTVectOfInt = vector<int>;

//-------------------------------------------------------------------------

class radTGroup : public radTg3d {
public:
	radTmhg GroupMapOfHandlers;

	radTGroup() {}
	~radTGroup() {}

	int Type_g3d() override { return 2;}
	virtual int Type_Group() { return 0;}

	inline void AddElement(int, const radThg&);

	inline void B_comp(radTField*) override;  // Modified by P. Elleaume 8 Nov 96
	// B_genComp is NOT overridden: groups inherit radTg3d::B_genComp, which applies
	// g3dListOfTransform by transforming the observation point / field (NestedFor_B).
	// Field evaluation must never mutate shared element state -- batch evaluation
	// (rad.Fld / RadiaField CF) calls B_genComp concurrently from ParallelFor threads.
	void B_intComp(radTField* FieldPtr) { B_comp(FieldPtr);} // This is not an Error!!!

	// Dump / DumpPureObjInfo / DumpBin REMOVED (Phase B2b/B2c, 2026-04-15)


	int DuplicateItself(radThg& hg, radTApplication* radPtr, char PutNewStuffIntoGenCont) override
	{
		return DuplicateGroupStuff(new radTGroup(*this), hg, radPtr, PutNewStuffIntoGenCont);
	}
	int DuplicateGroupStuff(radTGroup*, radThg&, radTApplication*, char);
	virtual int DuplicateWithoutDuplicatingGroupStuff(radThg& hgGroup)
	{
		radTSend Send;
		radTGroup* pNewGroup = new radTGroup(*this);
		if(pNewGroup==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hgLoc(pNewGroup);
		hgGroup = hgLoc;
		return 1;
	}

	int CreateFromSym(radThg&, radTApplication*, char);

	// All Subdivide* / CutItself / FindLowestAndUppestVertices / SetUpCuttingPlanes / FlattenNestedStructure*
	// CollectNonGroupElements / JustTraverse / SimpleEnergyComp / ActualEnergyForceTorqueCompWithAdd
	// MarkFurtherSubdNeed* / SetupFurtherSubdInd / SetMessageChar / NextStepEnergyForceTorqueComp
	// ProceedNextStepEnergyForceTorqueComp / EstimateCenterPointOverRelaxables / FindEdgePointsOver*
	// SubdivideByEllipses / SubdivideByPlanesPerpToCylAx REMOVED (Phase C, 2026-04-16)

	int ConvertToPolyhedron(radThg&, radTApplication*, char);

	int SetMaterial(radThg&, radTApplication*) override;
	void SetM(TVector3d& M) override; //virtual
	inline int ScaleCurrent(double) override; //virtual in radTg3d

	inline int ItemIsNotFullyInternalAfterCut();

	inline int NumberOfDegOfFreedom() override;
	inline int SizeOfThis() override;

	void Push_backCenterPointAndField(radTFieldKey*, radTVectPairOfVect3d*, radTrans*, radTg3d*, radTApplication*);

	inline double Volume();
	inline void VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique);

	radTGroup* CreateGroupIncludingAllMembersExceptIt(const radTmhg::const_iterator&);
};

//-------------------------------------------------------------------------

inline void radTGroup::AddElement(int ElemKey, const radThg& hg)
{
	GroupMapOfHandlers[ElemKey] = hg;
	static_cast<radTg3d*>(hg.rep)->IsGroupMember = 1;
}

//-------------------------------------------------------------------------

inline void radTGroup::B_comp(radTField* FieldPtr)
{
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
		static_cast<radTg3d*>(iter->second.rep)->B_genComp(FieldPtr);
}

//-------------------------------------------------------------------------

inline int radTGroup::NumberOfDegOfFreedom()
{
	int DegFrCount = 0;
	for(const auto& pair : GroupMapOfHandlers)
		DegFrCount += static_cast<radTg3d*>(pair.second.rep)->NumberOfDegOfFreedom();
	return DegFrCount;
}

//-------------------------------------------------------------------------

inline int radTGroup::SizeOfThis()
{
	int GenSize = sizeof(*this);
	for(const auto& pair : GroupMapOfHandlers)
		GenSize += pair.second.rep->SizeOfThis();
	return GenSize;
}

//-------------------------------------------------------------------------

inline int radTGroup::ItemIsNotFullyInternalAfterCut()
{
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
		if(static_cast<radTg3d*>(iter->second.rep)->ItemIsNotFullyInternalAfterCut()) return 1;
	return 0;
}

//-------------------------------------------------------------------------

inline double radTGroup::Volume()
{
	double SumVol = 0.;
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
		//SumVol += static_cast<radTg3d*>(iter->second.rep)->Volume();
		SumVol += static_cast<radTg3d*>(iter->second.rep)->VolumeWithSym();

	return SumVol;
}

//-------------------------------------------------------------------------

inline void radTGroup::VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique)
{
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
		static_cast<radTg3d*>(iter->second.rep)->VerticesInLocFrame(OutVect, EnsureUnique);
}

//-------------------------------------------------------------------------

// radTGroup::SimpleEnergyComp REMOVED (Phase C, 2026-04-16, energy-based API gone)
// MarkFurtherSubdNeed / MarkFurtherSubdNeed1D / SetupFurtherSubdInd / SetMessageChar REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

inline int radTGroup::ScaleCurrent(double scaleCoef)
{
	int scalingWasApplied = 0;
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		if(static_cast<radTg3d*>(iter->second.rep)->ScaleCurrent(scaleCoef)) scalingWasApplied = 1;
	}
	return scalingWasApplied;
}

//-------------------------------------------------------------------------

#endif
