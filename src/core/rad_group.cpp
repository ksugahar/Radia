/*-------------------------------------------------------------------------
*
* File name:      radgroup.cpp
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 group (/container) of magnetic field sources
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#include "rad_application.h"
#include "rad_group.h"
#include "rad_geometry_3d_aux.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

// radTGroup::Dump / DumpPureObjInfo / DumpBin REMOVED (Phase B2b/B2c, 2026-04-15)

//-------------------------------------------------------------------------


//-------------------------------------------------------------------------

int radTGroup::DuplicateGroupStuff(radTGroup* NewGroupPtr, radThg& hg, radTApplication* radPtr, char PutNewStuffIntoGenCont)
{
	radTSend Send;
	if(NewGroupPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
	NewGroupPtr->IsGroupMember = 0;
	NewGroupPtr->GroupMapOfHandlers.erase(NewGroupPtr->GroupMapOfHandlers.begin(), NewGroupPtr->GroupMapOfHandlers.end());
	radThg hgLoc(NewGroupPtr);

	int NewElemKey, NewStuffCounter = 0;
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin();
		iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg hgSubLoc;
		if(!((radTg3d*)(((*iter).second).rep))->DuplicateItself(hgSubLoc, radPtr, PutNewStuffIntoGenCont)) return 0;

		if(PutNewStuffIntoGenCont)
		{
			NewElemKey = radPtr->AddElementToContainer(hgSubLoc);
			NewGroupPtr->AddElement(NewElemKey, hgSubLoc);
		}
		else NewGroupPtr->AddElement(++NewStuffCounter, hgSubLoc);
	}
	hg = hgLoc; return 1;
	// This really creates new copies of all the GroupMembers
}

//-------------------------------------------------------------------------

int radTGroup::SetMaterial(radThg& InMatHandle, radTApplication* ApPtr)
{
	char PutNewStuffIntoGenCont = 1; // For Material: Maybe not necessary?
	radTMaterial* pMat = static_cast<radTMaterial*>(InMatHandle.rep);
	char EasyAxisDefinedInMat = pMat->EasyAxisDefined;

	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg hgMat = InMatHandle;
		if(!EasyAxisDefinedInMat) 
		{
			if(!pMat->DuplicateItself(hgMat, ApPtr, PutNewStuffIntoGenCont)) return 0;
			if(PutNewStuffIntoGenCont) ApPtr->AddElementToContainer(hgMat); // Maybe not necessary
		}
		if(!((radTg3d*)(((*iter).second).rep))->SetMaterial(hgMat, ApPtr)) return 0;
	}
	return 1;
}

//-------------------------------------------------------------------------

void radTGroup::SetM(TVector3d& M)
{
	for(radTmhg::const_iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		((radTg3d*)(((*iter).second).rep))->SetM(M);
	}
}

//-------------------------------------------------------------------------
// radTGroup::B_genComp override REMOVED (2026-07-10).  It propagated the
// group's transforms by push_front/restore on every CHILD's
// g3dListOfTransform, mutating shared elements during field evaluation --
// concurrent B_genComp calls from the batch ParallelFor (rad.Fld batch,
// RadiaField CF assembly) then corrupted the heap (0xC0000374/0xC0000005).
// Groups now inherit radTg3d::B_genComp (NestedFor_B), which applies the
// same transforms via observation-point/field transformation without
// touching any shared state.
//-------------------------------------------------------------------------

// radTGroup::SubdivideItself REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::SubdivideItselfAsWholeInLabFrame REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::SetUpCuttingPlanes REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::FindLowestAndUppestVertices REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::SubdivideItselfByOneSetOfParPlanes REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::CutItself REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::SubdivideItselfByParPlanes REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::SubdivideItselfByParPlanesAsWholeInLabFrame REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::FlattenNestedStructure REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::CollectNonGroupElements REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::SubdivideItselfByEllipticCylinder REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::SubdivideItselfByEllipticCylinderAsWholeInLabFrame REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::FindEdgePointsOverPhiAndAxForCylSubd REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

int radTGroup::ConvertToPolyhedron(radThg& In_hg, radTApplication* radPtr, char ReplaceOldStuff)
{
	//radTCast Cast;
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg& NewHandle = (*iter).second;
		radTGroup* pGroup = radTCast::GroupCast(static_cast<radTg3d*>(NewHandle.rep));
		if(pGroup != 0)
		{
			radThg OldHandle = NewHandle;
			if(!pGroup->ConvertToPolyhedron(NewHandle, radPtr, ReplaceOldStuff)) return 0;
		}
		else
		{
			radTg3dRelax* g3dRelaxPtr = radTCast::g3dRelaxCast(static_cast<radTg3d*>(NewHandle.rep));
			if(g3dRelaxPtr != 0)
			{
				radThg OldHandle = NewHandle;
				if(!g3dRelaxPtr->ConvertToPolyhedron(NewHandle, radPtr, ReplaceOldStuff)) return 0;

				if(ReplaceOldStuff)
				{
					radPtr->ReplaceInGlobalMap(OldHandle, NewHandle);
					if((static_cast<radTg3d*>(OldHandle.rep))->IsGroupMember) radPtr->ReplaceInAllGroups(OldHandle, NewHandle);
				}
			}
		}
	}
	return 1;
}

//-------------------------------------------------------------------------

// radTGroup::EstimateCenterPointOverRelaxables REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::SubdivideItselfOverAzimuth REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::FindEdgePointsOverEllipseSet0 REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::FindEdgePointsOverEllipseSet REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::SubdivideByEllipses REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::SubdivideByPlanesPerpToCylAx REMOVED (Phase C, 2026-04-16)
		
//-------------------------------------------------------------------------

// radTGroup::FlattenNestedStructureIfMessageCharIsSet REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::JustTraverse REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

int radTGroup::CreateFromSym(radThg& In_hg, radTApplication* radPtr, char PutNewStuffIntoGenCont)
{
	radThg hgDpl;
	if(!DuplicateWithoutDuplicatingGroupStuff(hgDpl)) return 0;
	radTGroup* pGroupDpl = (radTGroup*)(static_cast<radTg3d*>(hgDpl.rep));

	pGroupDpl->GroupMapOfHandlers.erase(pGroupDpl->GroupMapOfHandlers.begin(), pGroupDpl->GroupMapOfHandlers.end());

	int NewStuffCounter = 0;
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg hgNew = (*iter).second;
		radThg hgOld = hgNew;
		(static_cast<radTg3d*>(hgOld.rep))->CreateFromSym(hgNew, radPtr, PutNewStuffIntoGenCont);

		if(PutNewStuffIntoGenCont)
		{
			int NewElemKey = radPtr->AddElementToContainer(hgNew);
			pGroupDpl->AddElement(NewElemKey, hgNew);
		}
		else pGroupDpl->AddElement(++NewStuffCounter, hgNew);
	}

	int MaxMult = 0;
	for(radTlphg::const_iterator iterTr = g3dListOfTransform.begin(); iterTr != g3dListOfTransform.end(); ++iterTr)
		if((*iterTr).m > MaxMult) MaxMult = (*iterTr).m;

	In_hg = hgDpl;
	if(MaxMult > 1) 
	{
		radTSend Send;
		radTGroup* pGroup = new radTGroup();
		if(pGroup==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hgGroup(pGroup);

		radTrans BaseTrans;
		BaseTrans.SetupIdent();
		if(!pGroupDpl->NestedFor_CreateFromSym(pGroup, radPtr, PutNewStuffIntoGenCont, &BaseTrans, pGroupDpl->g3dListOfTransform.begin())) return 0;
		if(!pGroup->GroupMapOfHandlers.empty()) 
		{
			pGroup->IsGroupMember = pGroupDpl->IsGroupMember;
			pGroup->ConsiderOnlyWithTrans = pGroupDpl->ConsiderOnlyWithTrans;
			// pGroup->MessageChar = pGroupDpl->MessageChar; REMOVED (Phase C, 2026-04-16)

			In_hg = hgGroup;
		}
	}
	return 1;
}

//-------------------------------------------------------------------------

void radTGroup::Push_backCenterPointAndField(radTFieldKey* pFieldKey, radTVectPairOfVect3d* pVectPairOfVect3d, radTrans* pBaseTrans, radTg3d* g3dSrcPtr, radTApplication* pAppl)
{// Attention: this assumes no more than one transformation with mult. no more than 1 !!!
	radTrans* pTrans = (g3dListOfTransform.empty())? 0 : (radTrans*)((*(g3dListOfTransform.begin())).Handler_g.rep);
	radTrans TotTrans;
	if(pTrans != 0)
	{
		if(pBaseTrans != 0) 
		{
			TrProduct(pBaseTrans, pTrans, TotTrans);
			pTrans = &TotTrans;
		}
	}
	else
	{
		if(pBaseTrans != 0) pTrans = pBaseTrans;
	}
	
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		radThg hg = (*iter).second;
		radTg3d* g3dPtr = static_cast<radTg3d*>(hg.rep);
		//g3dPtr->Push_backCenterPointAndField(pFieldKey, pVectPairOfVect3d, pTrans, g3dSrcPtr);

		radThg hgDplWithoutSym; //OC061007_BNL
		char PutNewStuffIntoGenCont = 0;
		if(!g3dPtr->CreateFromSym(hgDplWithoutSym, pAppl, PutNewStuffIntoGenCont)) return;

		radTg3d* g3dDplWithoutSymPtr = static_cast<radTg3d*>(hgDplWithoutSym.rep);

		radTvhg vhFlatTransforms; //OC061007_BNL
		g3dDplWithoutSymPtr->FlattenSpaceTransforms(vhFlatTransforms);
		if(vhFlatTransforms.size() > 0)
		{
			g3dDplWithoutSymPtr->EraseAllTransformations();
			g3dDplWithoutSymPtr->AddTransform(1, vhFlatTransforms[0]);
		}

		g3dDplWithoutSymPtr->Push_backCenterPointAndField(pFieldKey, pVectPairOfVect3d, pTrans, g3dSrcPtr, pAppl);
	}
}

//-------------------------------------------------------------------------

// radTGroup::NextStepEnergyForceTorqueComp REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::ProceedNextStepEnergyForceTorqueComp REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTGroup::ActualEnergyForceTorqueCompWithAdd REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

radTGroup* radTGroup::CreateGroupIncludingAllMembersExceptIt(const radTmhg::const_iterator& it)
{
	radTGroup* OutGroupPtr = new radTGroup();

	int LocInd = -1;
	for(radTmhg::iterator iter = GroupMapOfHandlers.begin(); iter != GroupMapOfHandlers.end(); ++iter)
	{
		LocInd++;
		if(iter == it) continue;

		radThg cur_hg = (*iter).second;
		OutGroupPtr->AddElement(LocInd, cur_hg); 
	}
	return OutGroupPtr;
}

//-------------------------------------------------------------------------
