/*-------------------------------------------------------------------------
*
* File name:      radapl3.cpp
*
* Project:        RADIA
*
* Description:    Wrapping RADIA application function calls
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#include "rad_application.h"
#include "rad_geometry_3d_aux.h"
#include "rad_operation_names.h"
#include "rad_interaction.h"
#include "rad_field_unified.h"

#include <math.h>
#include <string.h>
#include <vector>

#include "rad_parallel.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTApplication::ComputeField(int ElemKey, char* FieldChar, double* StObsPoi, long lenStObsPoi,
								   double* FiObsPoi, long lenFiObsPoi, int Np, char* ShowArgFlag, double StrtArg)
{
	radTField* FieldArray = nullptr;
	double* ArgArray = nullptr;
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return;}
		radTFieldKey FieldKey;
		if(!ValidateFieldChar(FieldChar, &FieldKey)) return;
		TVector3d StObsPoiVect, FiObsPoiVect;
		if(!ValidateVector3d(StObsPoi, lenStObsPoi, &StObsPoiVect)) return;
		if(!ValidateVector3d(FiObsPoi, lenFiObsPoi, &FiObsPoiVect)) return;

		// Setup IMA context if IMA was used in the last solve AND we're computing field for that model
		// CRITICAL: Only set IMA context if ElemKey matches the cached model (m_cached_obj_key)
		bool imaWasSet = false;
		if(m_cached_interact_key > 0 && m_cached_obj_key == ElemKey)
		{
			radTInteraction* pIntrc = GetInteractionByKey(m_cached_interact_key);
			if(pIntrc && pIntrc->IsIMAEnabled())
			{
				RadIMAFieldContext::Set(
					pIntrc->GetIMASymmetry(),
					pIntrc->GetIMASignX(),
					pIntrc->GetIMASignY(),
					pIntrc->GetIMASignZ()
				);
				imaWasSet = true;
			}
		}

		short ArgumentNeeded = 0;
		if(!strcmp(ShowArgFlag, "arg")) ArgumentNeeded = 1;
		else if(strcmp(ShowArgFlag, "noarg")) { Send.ErrorMessage("Radia::Error034"); return;}

		if(Np==1 && (FiObsPoiVect.x < 1.E+22) && (FiObsPoiVect.y < 1.E+22) && (FiObsPoiVect.z < 1.E+22)) Np = 101; // New Default

		TVector3d ZeroVect(0.,0.,0.);
		TVector3d ObsPoiVect = StObsPoiVect;

		radTField Field(FieldKey, CompCriterium, ObsPoiVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
		g3dPtr->B_genComp(&Field);

		std::vector<radTField> vFieldArray;
		std::vector<double> vArgArray;
		if(Np>1)
		{
			vFieldArray.resize(Np);
			FieldArray = vFieldArray.data();
			if(ArgumentNeeded)
			{
				vArgArray.resize(Np);
				ArgArray = vArgArray.data();
			}

			FieldArray[0] = Field;
			TVector3d TranslVect = (1./double(Np-1))*(FiObsPoiVect-StObsPoiVect);
			double StepArg;
			if(ArgumentNeeded)
			{
				ArgArray[0] = StrtArg;
				StepArg	= sqrt(TranslVect.x*TranslVect.x + TranslVect.y*TranslVect.y + TranslVect.z*TranslVect.z);
			}

		ngcore::ParallelFor(ngcore::IntRange(1, Np), [&](size_t i) {
			TVector3d LocalObsPoiVect = StObsPoiVect + double(i) * TranslVect;
				FieldArray[i] = radTField(FieldKey, CompCriterium, LocalObsPoiVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
				g3dPtr->B_genComp(&(FieldArray[i]));
				if(ArgumentNeeded) ArgArray[i] = StrtArg + double(i) * StepArg;
			});
		}
		else FieldArray = &Field;

		if(SendingIsRequired) Send.OutFieldCompRes(FieldChar, FieldArray, ArgArray, Np);

		// Clear IMA context after computation
		if(imaWasSet) RadIMAFieldContext::Clear();
		// RAII: vFieldArray and vArgArray cleaned up automatically
	}
	catch(...)
	{
		// Clear IMA context on error
		RadIMAFieldContext::Clear();
		// RAII: vFieldArray and vArgArray cleaned up automatically
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeField(int ElemKey, char* FieldChar, radTVectorOfVector3d& VectorOfVector3d, radTVectInputCell& VectInputCell)
{
	radTField* FieldArray = 0;
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return;}
		radTFieldKey FieldKey;
		if(!ValidateFieldChar(FieldChar, &FieldKey)) return;

		// Setup IMA context if IMA was used in the last solve AND we're computing field for that model
		// CRITICAL: Only set IMA context if ElemKey matches the cached model (m_cached_obj_key)
		bool imaWasSet = false;
		if(m_cached_interact_key > 0 && m_cached_obj_key == ElemKey)
		{
			radTInteraction* pIntrc = GetInteractionByKey(m_cached_interact_key);
			if(pIntrc && pIntrc->IsIMAEnabled())
			{
				RadIMAFieldContext::Set(
					pIntrc->GetIMASymmetry(),
					pIntrc->GetIMASignX(),
					pIntrc->GetIMASignY(),
					pIntrc->GetIMASignZ()
				);
				imaWasSet = true;
			}
		}

		long Np = (long)(VectorOfVector3d.size());
		std::vector<radTField> vFieldArray(Np);
		FieldArray = vFieldArray.data();
		radTField* tField = FieldArray;

		TVector3d ZeroVect(0.,0.,0.);
		ngcore::ParallelFor(ngcore::IntRange(Np), [&](size_t i) {
			FieldArray[i] = radTField(FieldKey, CompCriterium, VectorOfVector3d[i], ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
			g3dPtr->B_genComp(&(FieldArray[i]));
		});

		if(SendingIsRequired) OutFieldCompRes(FieldChar, FieldArray, Np, VectInputCell);

		// Clear IMA context after computation
		if(imaWasSet) RadIMAFieldContext::Clear();
		// RAII: vFieldArray cleaned up automatically
	}
	catch(...)
	{
		// Clear IMA context on error
		RadIMAFieldContext::Clear();
		// RAII: vFieldArray cleaned up automatically
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeField(int ElemKey, char* FieldChar, double** Points, long Np)
{
	radTField *FieldArray = 0;
	double *arFldVals = 0, *arFldValsRecv = 0; //OC02012020
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return;}
		radTFieldKey FieldKey;
		if(!ValidateFieldChar(FieldChar, &FieldKey)) return;

		// Setup IMA context if IMA was used in the last solve AND we're computing field for that model
		// CRITICAL: Only set IMA context if ElemKey matches the cached model (m_cached_obj_key)
		bool imaWasSet = false;
		if(m_cached_interact_key > 0 && m_cached_obj_key == ElemKey)
		{
			radTInteraction* pIntrc = GetInteractionByKey(m_cached_interact_key);
			if(pIntrc && pIntrc->IsIMAEnabled())
			{
				RadIMAFieldContext::Set(
					pIntrc->GetIMASymmetry(),
					pIntrc->GetIMASignX(),
					pIntrc->GetIMASignY(),
					pIntrc->GetIMASignZ()
				);
				imaWasSet = true;
			}
		}

		std::vector<radTField> vFieldArray;
		if(m_nProcMPI < 2) //OC01012020
		{
			vFieldArray.resize(Np);
			FieldArray = vFieldArray.data();
			radTField* tField = FieldArray;

			TVector3d ZeroVect(0.,0.,0.), v;
			double **tPoints = Points;
			ngcore::ParallelFor(ngcore::IntRange(Np), [&](size_t i) {
				double *t = Points[i];
				TVector3d v; v.x = *(t++); v.y = *(t++); v.z = *t;

				FieldArray[i] = radTField(FieldKey, CompCriterium, v, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
				g3dPtr->B_genComp(&(FieldArray[i]));
			});
			if(SendingIsRequired) OutFieldCompRes(FieldChar, FieldArray, Np);
			// RAII: vFieldArray cleaned up automatically
		}

		// Clear IMA context after computation
		if(imaWasSet) RadIMAFieldContext::Clear();
	}
	catch(...)
	{
		// Clear IMA context on error
		RadIMAFieldContext::Clear();
		// RAII: vFieldArray cleaned up automatically
		if(arFldVals != 0) delete[] arFldVals; //OC02012020
		if(arFldValsRecv != 0) delete[] arFldValsRecv; //OC02012020
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeFieldInt(int ElemKey, char* FinOrInfChar, char* FieldIntChar, double* StPoi, long lenStPoi, double* FiPoi, long lenFiPoi)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return;}
		TVector3d StPoiVect;
		if(!ValidateVector3d(StPoi, lenStPoi, &StPoiVect)) return;
		TVector3d FiPoiVect;
		if(!ValidateVector3d(FiPoi, lenFiPoi, &FiPoiVect)) return;

		radTFieldKey FieldKey;
		if(!ValidateFieldIntChar(FieldIntChar, FinOrInfChar, &FieldKey)) return;

		TVector3d ZeroVect(0.,0.,0.);
		radTField Field(FieldKey, CompCriterium, StPoiVect, FiPoiVect, ZeroVect, ZeroVect);

		g3dPtr->B_genComp(&Field);

		if(SendingIsRequired) Send.OutFieldIntCompRes(FieldIntChar, &Field);
	}
	catch(...) 
	{ 
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeFieldForce(int SourceElemKey, int ShapeElemKey)
{
	try
	{
		radThg hSource;
		if(!ValidateElemKey(SourceElemKey, hSource)) return;
		radTg3d* SourcePtr = Cast.g3dCast(hSource.rep);
		if(SourcePtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		radThg hShape;
		if(!ValidateElemKey(ShapeElemKey, hShape)) return;
		radTg3d* ShapePtr = Cast.g3dCast(hShape.rep);
		if(ShapePtr==0) { Send.ErrorMessage("Radia::Error003"); return;}
		radTRectangle* RectanglePtr = Cast.RectangleCast(ShapePtr); 
		if(RectanglePtr==0)
		{
			radTg3dRelax* g3dRelaxPtr = Cast.g3dRelaxCast(ShapePtr);
			if(g3dRelaxPtr!=0)
			{
				radTRecCur* RecMagPtr = Cast.RecCurCast(g3dRelaxPtr);
				if(RecMagPtr==0) { Send.ErrorMessage("Radia::Error036"); return;}
			}
			else { Send.ErrorMessage("Radia::Error036"); return;}
			// Modify this later (incl. Error message), as integration methods for other primitives are ready
			// What about Group Shape?
		}

		radTFieldKey FieldKey;
		FieldKey.Force_= 1;

		radTStructForShapeInt ShapeIntData;
		ShapeIntData.HandleOfSource = hSource;
		ShapeIntData.HandleOfShape = hShape;
		ShapeIntData.IntegrandLength = 1; // Number of elements in TVector3d* to be integrated over a Shape
		ShapeIntData.IntegrandFunPtr = &radTg3d::NormStressTensor;
		ShapeIntData.IntOverLine_= ShapeIntData.IntOverVol_= 0;
		ShapeIntData.IntOverSurf_= 1;
		ShapeIntData.AbsPrecArray = &(CompCriterium.AbsPrecForce);
		TVector3d LocForce(0.,0.,0.);
		ShapeIntData.VectArray = &LocForce;
		char LocForceType = 'r'; // 'a' - axial, 'r' - regular
		ShapeIntData.VectTypeArray = &LocForceType;

		TVector3d ZeroVect(0.,0.,0.);
		radTField Field(FieldKey, CompCriterium, &ShapeIntData /*, ZeroVect */);

		ShapePtr->B_genComp(&Field);

		if(SendingIsRequired) Send.Vector3d(&LocForce);
	}
	catch(...) 
	{ 
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

// ComputeFieldEnergy / ComputeFieldForceThroughEnergy / ComputeFieldTorqueThroughEnergy REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------
/**
void radTApplication::OutFieldForceOrTorqueThroughEnergyCompRes(char* ForceComponID, TVector3d& Vect, char ID)
{// This is only for Force and Torque!
	char* BufChar = ForceComponID;
	char* EqEmptyStr = (ID=='f')? "FxFyFz" : "TxTyTz";

	char SmallID = ID;
	char CapitalID = (SmallID=='f')? 'F' : 'T';

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			char* BufChar_pl_1 = BufChar+1;
			if((((*BufChar==CapitalID) || (*BufChar==SmallID)) && 
			   (*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z')) ||
			   (*BufChar == 'X') || (*BufChar == 'x') ||
			   (*BufChar == 'Y') || (*BufChar == 'y') ||
			   (*BufChar == 'Z') || (*BufChar == 'z')) ItemCount++;
			BufChar++;
		}
		BufChar = ForceComponID;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 3;
	}
	if(ItemCount > 1) Send.InitOutList(ItemCount);

	while (*BufChar != '\0') 
	{
		if((*(BufChar)==CapitalID) || (*(BufChar)==SmallID))
		{
			char* BufChar_pl_1 = BufChar+1;
			if((*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z')) Send.Vector3d(&Vect);
		}
		else if((*(BufChar)=='X') || (*(BufChar)=='x')) Send.Double(Vect.x);
		else if((*(BufChar)=='Y') || (*(BufChar)=='y')) Send.Double(Vect.y);
		else if((*(BufChar)=='Z') || (*(BufChar)=='z')) Send.Double(Vect.z);
		BufChar++;
	}
}
**/
// NOTE: ComputeParticleTrajectory, ComputeFocusingPotential, ComputeSecondOrderKickPer,
void radTApplication::ComputeShimSignature(int ElemKey, char* FldID, double* Disp, double* StPoi, double* FiPoi, int Np, double* Vi)
{
	radTField **arr_pField = 0, *arr_resField = 0;
	int TwoNp = Np << 1;
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return;
		radTg3d* SourcePtr = Cast.g3dCast(hg.rep);
		if(SourcePtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		radTFieldKey FieldKey;
		if(FldID == 0) { Send.ErrorMessage("Radia::Error096"); return;}
		if((*FldID == 'i') || (*FldID == 'I'))
		{
			//if(!ValidateFieldIntChar(FldID, "inf", &FieldKey)) 
			if(!ValidateFieldIntChar(FldID, (char*)"inf", &FieldKey)) //OC01052013 
			{
				Send.ErrorMessage("Radia::Error096"); return;
			}
		}
		else if(!ValidateFieldChar(FldID, &FieldKey, false))
		{
			Send.ErrorMessage("Radia::Error096"); return;
		}

		TVector3d vDisp;
		if(!ValidateVector3d(Disp, 3, &vDisp)) return;
		TVector3d vStPoi;
		if(!ValidateVector3d(StPoi, 3, &vStPoi)) return;
		TVector3d vFiPoi;
		if(!ValidateVector3d(FiPoi, 3, &vFiPoi)) return;
		TVector3d vVi;
		if(!ValidateVector3d(Vi, 3, &vVi)) return;

		if(Np <= 0) { Send.ErrorMessage("Radia::Error079"); return;}

		TVector3d ZeroVect(0.,0.,0.);
		std::vector<radTField*> vArr_pField(TwoNp, nullptr);
		arr_pField = vArr_pField.data();
		std::vector<radTField> vArr_resField(Np);
		arr_resField = vArr_resField.data();

		TVector3d vTransl = ((Np > 1)? (1./double(Np - 1)) : 1)*(vFiPoi - vStPoi);

		short prevSendingIsRequired = SendingIsRequired;
		SendingIsRequired = false;

		int indTranslat = SetTranslation(Disp, 3);
		int indSrcDisp = DuplicateElement_g3d(ElemKey, 0, 0, 0);
		indSrcDisp = ApplySymmetry(indSrcDisp, indTranslat, 1);
		radThg hgSrcDisp;
		if(!ValidateElemKey(indSrcDisp, hgSrcDisp)) return;
		radTg3d* SourceDispPtr = Cast.g3dCast(hgSrcDisp.rep);
		if(SourceDispPtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		radTg3d *arrSourceDispPtr[] = {SourcePtr, SourceDispPtr};

		radTField **tField = arr_pField;
		for(int k=0; k<2; k++)
		{
			radTg3d *curSourceDispPtr = arrSourceDispPtr[k];
			TVector3d vP1 = vStPoi;
			for(int i=0; i<Np; i++)
			{
				if(FieldKey.Ib_	|| FieldKey.Ih_)
				{
					TVector3d vP2 = vP1 + vVi;
					*tField = new radTField(FieldKey, CompCriterium, vP1, vP2, ZeroVect, ZeroVect);
				}
				else
				{
					*tField = new radTField(FieldKey, CompCriterium, vP1, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
				}

				curSourceDispPtr->B_genComp(*tField);

				if(k == 1) arr_resField[i] = **tField - *(arr_pField[i]);

				tField++;
				vP1 += vTransl;
			}
			//vStPoi += vDisp;
			curSourceDispPtr++;
		}

		DeleteElement(indSrcDisp);
		DeleteElement(indTranslat);

		SendingIsRequired = prevSendingIsRequired;

		if(SendingIsRequired)
		{
			if(FieldKey.Ib_	|| FieldKey.Ih_)
			{
				Send.OutFieldIntCompRes(FldID, arr_resField, nullptr, Np);
				//Send.OutFieldIntCompRes(FldID, &Field);
			}
			else
			{
				Send.OutFieldCompRes(FldID, arr_resField, nullptr, Np);
			}
		}

		if(arr_pField != 0) 
		{
			for(int i=0; i<TwoNp; i++) if(arr_pField[i] != 0) delete arr_pField[i];
			// RAII: vArr_pField cleaned up automatically
		}
		// RAII: vArr_resField cleaned up automatically
	}
	catch(...) 
	{ 
		if(arr_pField != 0) 
		{
			for(int i=0; i<TwoNp; i++) if(arr_pField[i] != 0) delete arr_pField[i];
			// RAII: vArr_pField cleaned up automatically
		}
		// RAII: vArr_resField cleaned up automatically
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

// FldUnits is deprecated. Radia always uses meters.
// These functions are kept as no-ops for backward compatibility.

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTApplication::ReplaceInAllGroups(radThg& OldHandle, radThg& NewHandle)
{
	radTg3d* g3dOldPtr = static_cast<radTg3d*>(OldHandle.rep);

	try
	{
		for(radTmhg::iterator GenIter = GlobalMapOfHandlers.begin(); GenIter != GlobalMapOfHandlers.end(); ++GenIter)
		{
			radTg* gP = ((*GenIter).second).rep;
			radTg3d* g3dP = Cast.g3dCast(gP); 
			if(g3dP != 0) 
			{
				radTGroup* GroupP = Cast.GroupCast(g3dP); 
				if(GroupP != 0) 
					for(radTmhg::iterator GroupIter = GroupP->GroupMapOfHandlers.begin();
						GroupIter != GroupP->GroupMapOfHandlers.end(); ++GroupIter)
						if(((*GroupIter).second).rep == g3dOldPtr) (*GroupIter).second = NewHandle;
			}
		}
	}
	catch(...) 
	{ 
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ReplaceInGlobalMap(radThg& OldHandle, radThg& NewHandle)
{
	radTg3d* g3dOldPtr = static_cast<radTg3d*>(OldHandle.rep);

	try
	{
		for(radTmhg::iterator GenIter = GlobalMapOfHandlers.begin(); GenIter != GlobalMapOfHandlers.end(); ++GenIter)
		{
			if(((*GenIter).second).rep == g3dOldPtr) 
			{
				(*GenIter).second = NewHandle;
			}
		}
	}
	catch(...) 
	{ 
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

// SubdivideElement_g3d / SubdivideElement_g3dByParPlanes / CutElement_g3d REMOVED (Phase C, 2026-04-16)
// FieldCompMethForSubdividedRecMag / SetLocMgnInSbdRecMag REMOVED (Phase C, 2026-04-16, radTSubdividedRecMag gone)

//-------------------------------------------------------------------------

int radTApplication::DuplicateElement_g3d(int ElemKey, const char** OptionNames, const char** OptionValues, int OptionCount)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); 
		radTrans* TransPtr = 0; //OC270307: added duplication of transformations
		if(g3dPtr==0) 
		{
			TransPtr = Cast.TransCast(hg.rep); 
			if(TransPtr == 0) 
			{
				Send.ErrorMessage("Radia::Error003"); return 0;
			}
		}

		if(TransPtr != 0) //OC270307: added duplication of transformations
		{
			radTrans *pNewTrans = new radTrans(*TransPtr);
			radThg hg(pNewTrans);
			int NewElemKey = AddElementToContainer(hg);
			if(SendingIsRequired) Send.Int(NewElemKey);
			return NewElemKey;
		}

		radTOptionNames OptNam;
		const char** BufNameString = OptionNames;
		const char** BufValString = OptionValues;

		char ReleaseSym = 0;
		for(int i=0; i<OptionCount; i++)
		{
			if(!strcmp(*BufNameString, OptNam.FreeSym))
			{
				if(!strcmp(*BufValString, (OptNam.FreeSymValues)[0])) ReleaseSym = 0;
				else if(!strcmp(*BufValString, (OptNam.FreeSymValues)[1])) ReleaseSym = 1;
				else if(!strcmp(*BufValString, (OptNam.FreeSymValues)[2])) ReleaseSym = 0;
				else if(!strcmp(*BufValString, (OptNam.FreeSymValues)[3])) ReleaseSym = 1;
				else { Send.ErrorMessage("Radia::Error062"); return 0;}
			}
			else { Send.ErrorMessage("Radia::Error062"); return 0;}
			BufNameString++; BufValString++;
		}

		char PutNewStuffIntoGenCont = 1;
		if(ReleaseSym)
		{
			if(!g3dPtr->CreateFromSym(hg, this, PutNewStuffIntoGenCont)) return 0;
		}
		else
		{
			if(!g3dPtr->DuplicateItself(hg, this, PutNewStuffIntoGenCont)) return 0;
		}

		int NewElemKey = AddElementToContainer(hg);

		if(SendingIsRequired) Send.Int(NewElemKey);
		return NewElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::CreateFromObj_g3dWithSym(int ElemKey)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); 
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		char PutNewStuffIntoGenCont = 1;
		if(!g3dPtr->CreateFromSym(hg, this, PutNewStuffIntoGenCont)) return 0;

		if(hg.rep != g3dPtr)
		{
			int NewElemKey = AddElementToContainer(hg);

			if(SendingIsRequired) Send.Int(NewElemKey);
			return NewElemKey;
		}
		else
		{
			if(SendingIsRequired) Send.Int(ElemKey);
			return ElemKey;
		}
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::DeleteAllElements(int DeletionMethNo)
{
	try
	{
		if(DeletionMethNo == 1) GlobalMapOfHandlers.erase(GlobalMapOfHandlers.begin(), GlobalMapOfHandlers.end());
		else if(DeletionMethNo == 2)
		{
			radTmhg::iterator IterStartDel = GlobalMapOfHandlers.end(); --IterStartDel;
			for(radTmhg::iterator iter = IterStartDel; iter != GlobalMapOfHandlers.begin(); --iter)
				GlobalMapOfHandlers.erase(iter);
			GlobalMapOfHandlers.erase(GlobalMapOfHandlers.begin());
		}
		GlobalUniqueMapKey = 1;

		// Reset solve cache - interaction object keys are now invalid
		m_cached_interact_key = 0;
		m_cached_obj_key = 0;

		// Clear element caches (stale after element deletion)
		RadFieldUnified::ClearAllCaches();

		if(SendingIsRequired) Send.Int(0);
		return 1;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::ComputeGeometricalVolume(int ElemKey)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); 
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		if(SendingIsRequired) Send.Double(g3dPtr->VolumeWithSym());
		return 1;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ReturnInput(double Input, int NumTimes)
{
	try
	{
		// RAII: Use std::vector for automatic cleanup (exception-safe)
		std::vector<double> vOutArray(NumTimes);
		double* OutArray = vOutArray.data();

		for(int i=0; i<NumTimes; ++i)
		{
			OutArray[i] = Input;
		}
		Send.DoubleList(OutArray, NumTimes);
		// RAII: automatic cleanup
	}
	catch(...)
	{
		// RAII: automatic cleanup even on exception
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetMemAllocMethForIntrctMatr(char* TotOrParts)
{
	short InMemAllocForIntrctMatrTotAtOnce = 0;
	if((!strcmp(TotOrParts, "tot")) || (!strcmp(TotOrParts, "Tot")) || (!strcmp(TotOrParts, "TOT"))) InMemAllocForIntrctMatrTotAtOnce = 1;
	else if((!strcmp(TotOrParts, "parts")) || (!strcmp(TotOrParts, "Parts")) || (!strcmp(TotOrParts, "PARTS"))) InMemAllocForIntrctMatrTotAtOnce = 0;
	else { Send.ErrorMessage("Radia::Error046"); return 0;}

	MemAllocForIntrctMatrTotAtOnce = InMemAllocForIntrctMatrTotAtOnce;

	if(SendingIsRequired) Send.Int(1);
	return 1;
}

//-------------------------------------------------------------------------
