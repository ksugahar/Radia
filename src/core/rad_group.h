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
	void B_genComp(radTField* FieldPtr) override;  // Override for transformation propagation
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

	int SubdivideItself(double*, radThg&, radTApplication*, radTSubdivOptions*) override;
	int SubdivideItselfAsWholeInLabFrame(double*, radThg&, radTApplication*, radTSubdivOptions*);
	int SetUpCuttingPlanes(TVector3d&, double*, radTSubdivOptions*, TVector3d*);
	int FindLowestAndUppestVertices(TVector3d&, radTSubdivOptions*, TVector3d&, TVector3d&, radTrans&, char&, char&);

	int SubdivideItselfByEllipticCylinder(double*, radTCylindricSubdivSpec*, radThg&, radTApplication*, radTSubdivOptions*) override;
	int SubdivideItselfByEllipticCylinderAsWholeInLabFrame(double*, radTCylindricSubdivSpec*, radThg&, radTApplication*, radTSubdivOptions*);
	int FindEdgePointsOverPhiAndAxForCylSubd(radTCylindricSubdivSpec*, TVector3d*, double*);
	int ConvertToPolyhedron(radThg&, radTApplication*, char);
	int EstimateCenterPointOverRelaxables(TVector3d&);
	int SubdivideItselfOverAzimuth(double*, double*, radTCylindricSubdivSpec*, radThg&, radTApplication*, radTSubdivOptions*);
	int FindEdgePointsOverEllipseSet0(double*, radTCylindricSubdivSpec*, radThg, TVector3d*, double*, radTSubdivOptions*);
	int FindEdgePointsOverEllipseSet(double*, radTCylindricSubdivSpec*, radThg, TVector3d*, double*, radTSubdivOptions*);
	int SubdivideByEllipses(double*, double*, radTCylindricSubdivSpec*, radThg&, radTApplication*, radTSubdivOptions*);
	int SubdivideByPlanesPerpToCylAx(double*, double*, TVector3d*, radTCylindricSubdivSpec*, radThg&, radTApplication*, radTSubdivOptions*);
	void FlattenNestedStructureIfMessageCharIsSet(radTApplication*, char, radTSubdivOptions*);
	void JustTraverse();

	int CutItself(TVector3d*, radThg&, radTPair_int_hg&, radTPair_int_hg&, radTApplication*, radTSubdivOptions*);

	int SubdivideItselfByParPlanes(double*, int, radThg&, radTApplication*, radTSubdivOptions*) override;
	int SubdivideItselfByParPlanesAsWholeInLabFrame(double*, int, radThg&, radTApplication*, radTSubdivOptions*);

	int SubdivideItselfByOneSetOfParPlanes(TVector3d&, TVector3d*, int, radThg&, radTApplication*, radTSubdivOptions*, radTvhg*) override;

	int SetMaterial(radThg&, radTApplication*) override;
	void SetM(TVector3d& M) override; //virtual
	inline int ScaleCurrent(double) override; //virtual in radTg3d

	void FlattenNestedStructure(radTApplication*, char);
	void CollectNonGroupElements(radTmhg*, int&, radTApplication*, char);

	inline int ItemIsNotFullyInternalAfterCut();

	inline int NumberOfDegOfFreedom() override;
	inline int SizeOfThis() override;

	void Push_backCenterPointAndField(radTFieldKey*, radTVectPairOfVect3d*, radTrans*, radTg3d*, radTApplication*);

	inline double Volume();
	inline void VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique);

	inline void SimpleEnergyComp(radTField*);
	void ActualEnergyForceTorqueCompWithAdd(radTField*);
	inline void MarkFurtherSubdNeed(char, char, char);
	inline void MarkFurtherSubdNeed1D(char, char);
	int NextStepEnergyForceTorqueComp(double*, radThg&, radTField*, char&);
	int ProceedNextStepEnergyForceTorqueComp(double*, radThg&, radTField*, radTField*, char&, char);
	inline void SetupFurtherSubdInd(char);

	inline void SetMessageChar(char);

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

inline void radTGroup::SimpleEnergyComp(radTField* FieldPtr)
{
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
	{
		auto* g3dPtr = static_cast<radTg3d*>(iter->second.rep);
		if(g3dPtr->g3dListOfTransform.empty()) g3dPtr->SimpleEnergyComp(FieldPtr);
		else
		{
			g3dPtr->NestedFor_Energy(FieldPtr, g3dPtr->g3dListOfTransform.begin());
		}
	}
}

//-------------------------------------------------------------------------

inline void radTGroup::MarkFurtherSubdNeed(char SubdNeedX, char SubdNeedY, char SubdNeedZ)
{
	radTg3d::MarkFurtherSubdNeed(SubdNeedX, SubdNeedY, SubdNeedZ);

	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
		static_cast<radTg3d*>(iter->second.rep)->MarkFurtherSubdNeed(SubdNeedX, SubdNeedY, SubdNeedZ);
}

//-------------------------------------------------------------------------

inline void radTGroup::MarkFurtherSubdNeed1D(char SubdNeed, char XorYorZ)
{
	radTg3d::MarkFurtherSubdNeed1D(SubdNeed, XorYorZ);

	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
		static_cast<radTg3d*>(iter->second.rep)->MarkFurtherSubdNeed1D(SubdNeed, XorYorZ);
}

//-------------------------------------------------------------------------

inline void radTGroup::SetupFurtherSubdInd(char InSubdInd)
{
	radTg3d::SetupFurtherSubdInd(InSubdInd);

	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
		static_cast<radTg3d*>(iter->second.rep)->SetupFurtherSubdInd(InSubdInd);
}

//-------------------------------------------------------------------------

inline void radTGroup::SetMessageChar(char InMessageChar)
{
	MessageChar = InMessageChar;
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
		static_cast<radTg3d*>(iter->second.rep)->SetMessageChar(InMessageChar);
}

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
