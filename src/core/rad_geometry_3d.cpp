/*-------------------------------------------------------------------------
*
* File name:      radg3d.cpp
*
* Project:        RADIA
*
* Description:    Base class for 3D objects - magnetic field sources
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#include "rad_application.h"
#include "rad_geometry_3d.h"
#include "rad_geometry_3d_aux.h"
//#include "rad_transform_def.h"
#include "gmtrans.h"
#include "rad_type_cast.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

extern radTApplication rad;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTg3d::B_genComp(radTField* FieldPtr)
{
	radTFieldKey& FieldKey = FieldPtr->FieldKey;

	if(g3dListOfTransform.empty())
	{
		if(FieldKey.Ib_ || FieldKey.Ih_) B_intComp(FieldPtr);
		else if(FieldKey.Force_) IntOverShape(FieldPtr);
		else B_comp(FieldPtr);
	}
	else
	{
		if(FieldKey.Force_) NestedFor_IntOverShape(FieldPtr, g3dListOfTransform.begin());
		else NestedFor_B(FieldPtr, g3dListOfTransform.begin());
	}
}

//-------------------------------------------------------------------------

void radTg3d::NestedFor_B(radTField* FieldPtr, const radTlphg::iterator& Iter)
{
	radTrans* TransPtr = (radTrans*)(((*Iter).Handler_g).rep);
	radTlphg::iterator LocalNextIter = Iter;
	LocalNextIter++;

	TVector3d ZeroVect(0.,0.,0.);

	short FldIntNeeded = FieldPtr->FieldKey.Ib_ || FieldPtr->FieldKey.Ih_; // Plus this string

	if((*Iter).m == 1)
	{
		radTField BufField = TransPtr->TrField_inv(*FieldPtr);
		BufField.P = TransPtr->TrPoint_inv(FieldPtr->P);
		if(FldIntNeeded) BufField.NextP = TransPtr->TrPoint_inv(FieldPtr->NextP); // Plus this string

		B_comp_Or_NestedFor(&BufField, LocalNextIter);

		BufField.P = FieldPtr->P;
		if(FldIntNeeded) BufField.NextP = FieldPtr->NextP; // Plus this string

		*FieldPtr = TransPtr->TrField(BufField);
	}
	else
	{
		radTField BufField(FieldPtr->FieldKey, FieldPtr->CompCriterium, FieldPtr->P, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
		if(FldIntNeeded) // Plus this
		{
			BufField.Ib = ZeroVect; BufField.Ih = ZeroVect; BufField.NextP = FieldPtr->NextP;
		}

		B_comp_Or_NestedFor(&BufField, LocalNextIter);
		radTField BufField1 = BufField;

		BufField = radTField(FieldPtr->FieldKey, FieldPtr->CompCriterium, TransPtr->TrPoint_inv(BufField.P), ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
		if(FldIntNeeded) // Plus this
		{
			BufField.Ib = ZeroVect; BufField.Ih = ZeroVect; 
			BufField.NextP = TransPtr->TrPoint_inv(BufField1.NextP);
		}

		int Mult = (*Iter).m;
		for(int km = 1; km < Mult-1; km++)
		{
			B_comp_Or_NestedFor(&BufField, LocalNextIter);
			BufField = TransPtr->TrField_inv(BufField);
			BufField.P = TransPtr->TrPoint_inv(BufField.P);
			if(FldIntNeeded) BufField.NextP = TransPtr->TrPoint_inv(BufField.NextP); // Plus this string
		}

		B_comp_Or_NestedFor(&BufField, LocalNextIter);
		radTField BufField2 = BufField;

		for(int km1 = 1; km1 < Mult; km1++)	BufField2 = TransPtr->TrField(BufField2);

		BufField1 += BufField2;
		*FieldPtr += BufField1;
	}
}

//-------------------------------------------------------------------------

// radTg3d::NestedFor_Energy REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

void radTg3d::NestedFor_IntOverShape(radTField* FieldPtr, const radTlphg::iterator& Iter)
{
	radTrans* TransPtr = (radTrans*)(((*Iter).Handler_g).rep);
	radTlphg::iterator LocalNextIter = Iter;
	LocalNextIter++;

	radTg3d* SourcePtr = (radTg3d*)((FieldPtr->ShapeIntDataPtr->HandleOfSource).rep);
	radTrans* InvTransPtr = new radTrans(*TransPtr); // This is to let handler work correctly
	InvTransPtr->Invert();
	radThg hInvTrans((radTg*)InvTransPtr);

	int LocLenVal = FieldPtr->ShapeIntDataPtr->IntegrandLength;
	TVector3d* LocVectArray = FieldPtr->ShapeIntDataPtr->VectArray;
	char* LocVectTypeArray = FieldPtr->ShapeIntDataPtr->VectTypeArray;
	TVector3d ZeroVect(0.,0.,0.);

	if((*Iter).m == 1)
	{
		for(int i=0; i<LocLenVal; i++)
		{
			if(LocVectTypeArray[i]=='r') LocVectArray[i] = TransPtr->TrVectField_inv(LocVectArray[i]);
			if(LocVectTypeArray[i]=='a') LocVectArray[i] = TransPtr->TrVectPoten_inv(LocVectArray[i]);
		}
		SourcePtr->AddTransform(1, hInvTrans);

		IntOverShape_Or_NestedFor(FieldPtr, LocalNextIter);

		for(int ii=0; ii<LocLenVal; ii++)
		{
			if(LocVectTypeArray[ii]=='r') LocVectArray[ii] = TransPtr->TrVectField(LocVectArray[ii]);
			if(LocVectTypeArray[ii]=='a') LocVectArray[ii] = TransPtr->TrVectPoten(LocVectArray[ii]);
		}
		SourcePtr->EraseOuterTransform();
	}
	else
	{
		radTStructForShapeInt LocStructForShapeInt1 = *(FieldPtr->ShapeIntDataPtr);
		radTStructForShapeInt LocStructForShapeInt = *(FieldPtr->ShapeIntDataPtr);
		std::vector<TVector3d> vLocVectArrayMult1(LocLenVal);
		std::vector<TVector3d> vLocVectArrayMult(LocLenVal);
		TVector3d* LocVectArrayMult1 = vLocVectArrayMult1.data();
		TVector3d* LocVectArrayMult = vLocVectArrayMult.data();
		LocStructForShapeInt1.VectArray = LocVectArrayMult1;
		LocStructForShapeInt.VectArray = LocVectArrayMult;
		for(int i=0; i<LocLenVal; i++)
		{
			LocVectArrayMult1[i] = LocVectArrayMult[i] = ZeroVect;
		}
		radTField BufField(FieldPtr->FieldKey, FieldPtr->CompCriterium, &LocStructForShapeInt1);

		IntOverShape_Or_NestedFor(&BufField, LocalNextIter);

		BufField = radTField(FieldPtr->FieldKey, FieldPtr->CompCriterium, &LocStructForShapeInt);
		SourcePtr->AddTransform(1, hInvTrans);

		int Mult = (*Iter).m;
		for(int km = 1; km < Mult-1; km++)
		{
			IntOverShape_Or_NestedFor(&BufField, LocalNextIter);

			for(int ii=0; ii<LocLenVal; ii++)
			{
				if(LocVectTypeArray[ii]=='r') LocVectArrayMult[ii] = TransPtr->TrVectField_inv(LocVectArrayMult[ii]);
				if(LocVectTypeArray[ii]=='a') LocVectArrayMult[ii] = TransPtr->TrVectPoten_inv(LocVectArrayMult[ii]);
			}
			SourcePtr->AddTransform(1, hInvTrans);
		}

		IntOverShape_Or_NestedFor(&BufField, LocalNextIter);
		for(int km1 = 1; km1 < Mult; km1++)
		{
			for(int iii=0; iii<LocLenVal; iii++)
			{
				if(LocVectTypeArray[iii]=='r') LocVectArrayMult[iii] = TransPtr->TrVectField(LocVectArrayMult[iii]);
				if(LocVectTypeArray[iii]=='a') LocVectArrayMult[iii] = TransPtr->TrVectPoten(LocVectArrayMult[iii]);
			}
			SourcePtr->EraseOuterTransform();
		}
		for(int iiii=0; iiii<LocLenVal; iiii++) LocVectArray[iiii] += LocVectArrayMult1[iiii] + LocVectArrayMult[iiii];
		// RAII: automatic cleanup
	}
}

//-------------------------------------------------------------------------

void radTg3d::B_intCompFinNum(radTField* FieldPtr)
{// This uses Newton method (n=3)
	const double IntegWeight[] = {3./8., 9./8., 9./8., 3./4.};

	TVector3d VectV = FieldPtr->NextP - FieldPtr->P;
	double Fact = sqrt(VectV.x*VectV.x + VectV.y*VectV.y + VectV.z*VectV.z);

	radTFieldKey LocFieldKey;
	short LocIb_ = FieldPtr->FieldKey.Ib_;
	short LocIh_ = FieldPtr->FieldKey.Ih_;
	LocFieldKey.B_ = LocIb_;
	LocFieldKey.H_ = LocIh_;

	radTCompCriterium LocCompCriterium;
	LocCompCriterium = FieldPtr->CompCriterium;

	TVector3d ZeroVect(0.,0.,0.), S_forB, S_forH, GenS_forB(0.,0.,0.), GenS_forH(0.,0.,0.),
			  IntForB(1.E+23, 1.E+23, 1.E+23), IntForH(1.E+23, 1.E+23, 1.E+23),
			  PrIntForB, PrIntForH;

	radTField LocField(LocFieldKey, LocCompCriterium, FieldPtr->P, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);

	double t_min = 0.;
	double Step_t, t;
	short IndForWeight, IndForPass;
	short NotFirstPass = 0;

	int AmOfPoi = 4;
	int AmOfPoi_mi_1;
	double PrecParamB_int, PrecParamH_int, PrecParamInt;
	PrecParamB_int = PrecParamH_int = 0.;
	PrecParamInt = 1.E+23; 

	while(PrecParamInt > FieldPtr->CompCriterium.AbsPrecB_int)
	{
		AmOfPoi_mi_1 = AmOfPoi - 1;
		Step_t = 1./AmOfPoi_mi_1;
		t = t_min;

		PrIntForB = IntForB; PrIntForH = IntForH;
		S_forB = ZeroVect; S_forH = ZeroVect;

		IndForWeight = IndForPass = 0;

		for(int i=0; i<AmOfPoi; i++)
		{
			if(IndForPass==2) IndForPass = 0;
			if(IndForWeight==4) IndForWeight = 1;
			if(NotFirstPass && (IndForPass==0)) goto BottomOfThisLoop;
			if(i==AmOfPoi_mi_1) IndForWeight = 0;

			LocField.P = FieldPtr->P + (t * VectV);
			if(LocIb_) LocField.B = ZeroVect;
			if(LocIh_) LocField.H = ZeroVect;

			B_comp(&LocField);
						
			if(LocIb_) S_forB += IntegWeight[IndForWeight] * LocField.B;
			if(LocIh_) S_forH += IntegWeight[IndForWeight] * LocField.H;

BottomOfThisLoop:
			IndForPass++; IndForWeight++;
			t += Step_t;
		}

		if(LocIb_)
		{
			GenS_forB += S_forB; 
			IntForB = Step_t * GenS_forB;
			PrecParamB_int = Fact * Max( Max( Abs(IntForB.x-PrIntForB.x), Abs(IntForB.y-PrIntForB.y)), Abs(IntForB.z-PrIntForB.z));
		}
		if(LocIh_)
		{
			GenS_forH += S_forH; 
			IntForH = Step_t * GenS_forH;
			PrecParamH_int = Fact * Max( Max( Abs(IntForH.x-PrIntForH.x), Abs(IntForH.y-PrIntForH.y)), Abs(IntForH.z-PrIntForH.z));
		}

		PrecParamInt = Max(PrecParamB_int, PrecParamH_int);
		AmOfPoi = AmOfPoi_mi_1 * 2 + 1;
		NotFirstPass = 1;
	}

	if(LocIb_) FieldPtr->Ib += Fact * IntForB;
	if(LocIh_) FieldPtr->Ih += Fact * IntForH;
}

//-------------------------------------------------------------------------

void radTg3d::NormStressTensor(radTField* FieldPtr)
{
	TVector3d ZeroVect(0.,0.,0.);
	FieldPtr->FieldKey.Force_= 0; 
	short PrevB_= FieldPtr->FieldKey.B_; FieldPtr->FieldKey.B_= 1;
	FieldPtr->B = ZeroVect;

	(static_cast<radTg3d*>(FieldPtr->ShapeIntDataPtr->HandleOfSource.rep))->B_genComp(FieldPtr);

	// Maxwell stress tensor: T = (1/mu_0) * (B B - 0.5 |B|^2 I)
	// Radia now uses SI units (meters) internally, matching ELF.
	// ConForStrTensInSI = 1/mu_0 = 1/(4*pi*1e-7) [m/H]
	const double ConForStrTensInSI = 1.0/(4*3.14159265358979*1.E-07);
	TVector3d LocB = FieldPtr->B;

	//Out normal projection of the Maxwell Stress Tensor
	*(FieldPtr->ShapeIntDataPtr->VectArray) =
		ConForStrTensInSI*((LocB*FieldPtr->ShapeIntDataPtr->Normal)*LocB
		-(0.5*(LocB*LocB))*FieldPtr->ShapeIntDataPtr->Normal);

	FieldPtr->FieldKey.Force_= 1; FieldPtr->FieldKey.B_= PrevB_;
}

//-------------------------------------------------------------------------

// radTg3d::DumpTransApplied REMOVED (Phase B2b, 2026-04-15)

//-------------------------------------------------------------------------

// radTg3d::ActualEnergyForceTorqueComp REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTg3d::EnergyForceTorqueComp REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTg3d::CheckIfMoreEnrFrcTrqCompNeededAndUpdate REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTg3d::ActualEnergyForceTorqueCompWithAdd REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTg3d::ProceedNextStepEnergyForceTorqueComp REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTg3d::NextStepEnergyForceTorqueComp REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTg3d::EnergyForceTorqueCompAutoDestSubd REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTg3d::CheckAxesExchangeForSubdInLabFrame REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTg3d::TransferSubdivisionStructToLocalFrame REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTg3d::FindEllipticCoordOfPoint REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTg3d::EstimateLengthAlongEllipse REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

int radTg3d::CreateFromSym(radThg& In_hg, radTApplication* radPtr, char PutNewStuffIntoGenCont)
{
	int MaxMult = 0;
	for(radTlphg::const_iterator iter = g3dListOfTransform.begin(); iter != g3dListOfTransform.end(); ++iter)
		if((*iter).m > MaxMult) MaxMult = (*iter).m;

	if(MaxMult <= 1) 
	//if((MaxMult <= 1) && (g3dListOfTransform.size() <= 1)) //OC061007_BNL
	{
		return DuplicateItself(In_hg, radPtr, PutNewStuffIntoGenCont);
	}
	else
	{
		radTSend Send;
		radTGroup* pGroup = new radTGroup();
		if(pGroup==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hgGroup(pGroup);

		radTrans BaseTrans;
		BaseTrans.SetupIdent();
		if(!NestedFor_CreateFromSym(pGroup, radPtr, PutNewStuffIntoGenCont, &BaseTrans, g3dListOfTransform.begin())) return 0;

		if(!pGroup->GroupMapOfHandlers.empty()) 
		{
			pGroup->IsGroupMember = IsGroupMember;
			pGroup->ConsiderOnlyWithTrans = ConsiderOnlyWithTrans;
			// pGroup->MessageChar = MessageChar; REMOVED (Phase C, 2026-04-16)
			In_hg = hgGroup;
		}
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTg3d::NestedFor_CreateFromSym(radTGroup* pGroup, radTApplication* radPtr, char PutNewStuffIntoGenCont, radTrans* BaseTransPtr, const radTlphg::iterator& Iter)
{
	radTrans* TransPtr = (radTrans*)(((*Iter).Handler_g).rep);
	radTlphg::iterator LocalNextIter = Iter;
	LocalNextIter++;

	radTrans LocTotTrans = *BaseTransPtr;

	if((*Iter).m == 1)
	{
		LocTotTrans = Product(LocTotTrans, *TransPtr);
		if(!CreateAndAddToGroupOrNestedFor(pGroup, radPtr, PutNewStuffIntoGenCont, &LocTotTrans, LocalNextIter)) return 0;
	}
	else
	{
		if(!CreateAndAddToGroupOrNestedFor(pGroup, radPtr, PutNewStuffIntoGenCont, &LocTotTrans, LocalNextIter)) return 0;
		int Mult = (*Iter).m;
		for(int km = 1; km < Mult; km++)
		{
			LocTotTrans = Product(LocTotTrans, *TransPtr);
			if(!CreateAndAddToGroupOrNestedFor(pGroup, radPtr, PutNewStuffIntoGenCont, &LocTotTrans, LocalNextIter)) return 0;
		}
	}
	return 1;
}

//-------------------------------------------------------------------------

int radTg3d::CreateAndAddToGroupOrNestedFor(radTGroup* pGroup, radTApplication* radPtr, char PutNewStuffIntoGenCont, radTrans* BaseTransPtr, const radTlphg::iterator& Iter)
{
	if(Iter == g3dListOfTransform.end())
	{
		radThg hgNew;
		if(!DuplicateItself(hgNew, radPtr, PutNewStuffIntoGenCont)) return 0;

		radTg3d* g3dDplPtr = static_cast<radTg3d*>(hgNew.rep);
		g3dDplPtr->EraseAllTransformations();

		double RelTol = 1.E-12;		
		if(!BaseTransPtr->IsIdent(RelTol))
		{
			radTSend Send;
			radTrans* pNewTrans = new radTrans(*BaseTransPtr);
			if(pNewTrans == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
			radThg hgTrans(pNewTrans);
			int LocMult = 1;
			g3dDplPtr->AddTransform(LocMult, hgTrans);
		}

		if(PutNewStuffIntoGenCont)
		{
			pGroup->AddElement(radPtr->AddElementToContainer(hgNew), hgNew);
		}
		else
		{
			int LocKey = (int)(pGroup->GroupMapOfHandlers.size()) + 1;
			pGroup->AddElement(LocKey, hgNew);
		}
		return 1;
	}
	else return NestedFor_CreateFromSym(pGroup, radPtr, PutNewStuffIntoGenCont, BaseTransPtr, Iter);
}

//-------------------------------------------------------------------------

double radTg3d::VolumeWithSym()
{
	int TotMult = TotalMultiplicity();
	return TotMult*Volume();
}

//-------------------------------------------------------------------------

int radTg3d::TotalMultiplicity()
{
	int TotMult = 1;
	if(g3dListOfTransform.empty()) return TotMult;

	for(radTlphg::reverse_iterator iter = g3dListOfTransform.rbegin(); iter != g3dListOfTransform.rend(); ++iter)
	{
		radTrans* pTrans = (radTrans*)(((*iter).Handler_g).rep);
		TotMult *= (*iter).m;
	}

	return TotMult;
}

//-------------------------------------------------------------------------

void radTg3d::FlattenSpaceTransforms(radTvhg& FlatTransforms)
{
	if(g3dListOfTransform.empty()) return;

	//radThg ihg(new radIdentTrans());
	radTrans *pOrigIdentTrf = new radTrans(); //OC061007
	pOrigIdentTrf->SetupIdent();
	radThg ihg(pOrigIdentTrf);
	FlatTransforms.push_back(ihg);

	for(radTlphg::reverse_iterator iter = g3dListOfTransform.rbegin(); iter != g3dListOfTransform.rend(); ++iter)
	{
		radTrans* pTrans = (radTrans*)(((*iter).Handler_g).rep);
		int mult = (*iter).m;

		int CurFlatTranSize = (int)(FlatTransforms.size());
		
		if(mult == 1)
		{
			for(int k=0; k<CurFlatTranSize; k++)
			{
				radTrans* pCurFlatTrans = (radTrans*)(FlatTransforms[k].rep);
				*pCurFlatTrans = Product(*pTrans, *pCurFlatTrans); //multiply flat trans from left
			}
			continue;
		}
		
		radTvhg AuxDuplVect;
		for(int k=0; k<CurFlatTranSize; k++)
		{
			radTrans* pCurFlatTrans = (radTrans*)(FlatTransforms[k].rep);
			radThg LocHg(new radTrans(*pCurFlatTrans));
			AuxDuplVect.push_back(LocHg);
		}

		for(int j=1; j<mult; j++)
		{
			for(int k=0; k<CurFlatTranSize; k++)
			{
				radTrans* pCurFlatTrans = (radTrans*)(AuxDuplVect[k].rep);
				*pCurFlatTrans = Product(*pTrans, *pCurFlatTrans); //multiply from left

				radThg LocHg(new radTrans(*pCurFlatTrans));
				FlatTransforms.push_back(LocHg);
			}
		}

		AuxDuplVect.erase(AuxDuplVect.begin(), AuxDuplVect.end());
	}
}

//-------------------------------------------------------------------------
// radTg3d::DumpBinParse_g3d / radTg3dRelax::DumpBinParse_g3dRelax / DumpMaterApplied REMOVED (Phase B2b/B2c, 2026-04-15)
//-------------------------------------------------------------------------

void radTg3dRelax::Push_backCenterPointAndField(radTFieldKey* pFieldKey, radTVectPairOfVect3d* pVectPairOfVect3d, radTrans* pBaseTrans, radTg3d* g3dSrcPtr, radTApplication* pAppl)
{// Attention: this assumes no more than one transformation with mult. no more than 1 !!!
	TVector3d CP = CentrPoint;
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

	if(pTrans != 0) CP = pTrans->TrPoint(CP);
	radTPairOfVect3d Pair(CP);

	if(pFieldKey->J_) return;
	else if(pFieldKey->M_) Pair.V2 = (pTrans == 0)? Magn : pTrans->TrVectField(Magn);
	else
	{
		radTCompCriterium CompCriterium;
		TVector3d ZeroVect(0.,0.,0.);
		radTField Field(*pFieldKey, CompCriterium, CP, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
		g3dSrcPtr->B_genComp(&Field);
		Pair.V2 = (pFieldKey->B_)? Field.B : ((pFieldKey->H_)? Field.H : ((pFieldKey->A_)? Field.A : ZeroVect));
	}
	pVectPairOfVect3d->push_back(Pair);
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTStructForShapeInt::radTStructForShapeInt(const radTStructForShapeInt& InStr) 
{
	HandleOfSource = InStr.HandleOfSource;
	HandleOfShape = InStr.HandleOfShape;
	IntegrandLength = InStr.IntegrandLength;
	Normal = InStr.Normal;
	VectArray = InStr.VectArray;
	VectTypeArray = InStr.VectTypeArray;
	IntegrandFunPtr = InStr.IntegrandFunPtr;
	IntOverLine_ = InStr.IntOverLine_;
	IntOverSurf_ = InStr.IntOverSurf_;
	IntOverVol_ = InStr.IntOverVol_;
	AbsPrecArray = InStr.AbsPrecArray;
}

//-------------------------------------------------------------------------
