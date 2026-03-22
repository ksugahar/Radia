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
				radTRecMag* RecMagPtr = Cast.RecMagCast(g3dRelaxPtr);
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

void radTApplication::ComputeFieldEnergy(int DestElemKey, int SourceElemKey, int* SubdivArray, long lenSubdivArray)
{
	radTStructForEnergyForceTorqueComp* StructForEnergyForceTorqueCompPtr = nullptr;
	try
	{
		radThg hDest;
		if(!ValidateElemKey(DestElemKey, hDest)) return;
		radTg3d* DestPtr = Cast.g3dCast(hDest.rep);
		if(DestPtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		radThg hSource;
		if(!ValidateElemKey(SourceElemKey, hSource)) return;
		radTg3d* SourcePtr = Cast.g3dCast(hSource.rep);
		if(SourcePtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		radTFieldKey FieldKey;
		FieldKey.Energy_ = 1;

		if((lenSubdivArray != 3) || (SubdivArray[0] < 0) || (SubdivArray[1] < 0) || (SubdivArray[2] < 0))
		{
			Send.ErrorMessage("Radia::Error021"); return;
		}
		//double ActualSubdivisionArray[] = {SubdivArray[0], 1., SubdivArray[1], 1., SubdivArray[2], 1.};
		double ActualSubdivisionArray[] = {(double)SubdivArray[0], 1., (double)SubdivArray[1], 1., (double)SubdivArray[2], 1.}; //OC101015

		StructForEnergyForceTorqueCompPtr = new radTStructForEnergyForceTorqueComp();
		StructForEnergyForceTorqueCompPtr->hSource = hSource;
		StructForEnergyForceTorqueCompPtr->hDest = hDest;
		StructForEnergyForceTorqueCompPtr->radPtr = this;
		StructForEnergyForceTorqueCompPtr->DestSubdivArray = ActualSubdivisionArray;
		StructForEnergyForceTorqueCompPtr->AutoDestSubdivision = CheckForAutoDestSubdivision(ActualSubdivisionArray);
		radTHandleStructForEnergyForceTorqueComp HandleStructForEnergyForceTorqueComp(StructForEnergyForceTorqueCompPtr);
		StructForEnergyForceTorqueCompPtr = nullptr;  // Ownership transferred to handle

		radTField Field(FieldKey, CompCriterium, HandleStructForEnergyForceTorqueComp);
		DestPtr->EnergyForceTorqueComp(&Field);

		if(Field.HandleEnergyForceTorqueCompData.rep->SomethingIsWrong) return;
		if(SendingIsRequired) Send.Double(Field.Energy);
	}
	catch(...)
	{
		if(StructForEnergyForceTorqueCompPtr) delete StructForEnergyForceTorqueCompPtr;  // Clean up if exception before handle ownership transfer
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeFieldForceThroughEnergy(int DestElemKey, int SourceElemKey, char* ForceComponID, int* SubdivArray, long lenSubdivArray)
{
	radTStructForEnergyForceTorqueComp* StructForEnergyForceTorqueCompPtr = nullptr;
	try
	{
		radThg hDest;
		if(!ValidateElemKey(DestElemKey, hDest)) return;
		radTg3d* DestPtr = Cast.g3dCast(hDest.rep);
		if(DestPtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		radThg hSource;
		if(!ValidateElemKey(SourceElemKey, hSource)) return;
		radTg3d* SourcePtr = Cast.g3dCast(hSource.rep);
		if(SourcePtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		if(!ValidateForceChar(ForceComponID)) return;

		radTFieldKey FieldKey;
		FieldKey.ForceEnr_ = 1;

		if((lenSubdivArray != 3) || (SubdivArray[0] < 0) || (SubdivArray[1] < 0) || (SubdivArray[2] < 0))
		{
			Send.ErrorMessage("Radia::Error021"); return;
		}
		//double ActualSubdivisionArray[] = {SubdivArray[0], 1., SubdivArray[1], 1., SubdivArray[2], 1.};
		double ActualSubdivisionArray[] = {(double)SubdivArray[0], 1., (double)SubdivArray[1], 1., (double)SubdivArray[2], 1.}; //OC101015

		StructForEnergyForceTorqueCompPtr = new radTStructForEnergyForceTorqueComp();
		StructForEnergyForceTorqueCompPtr->hSource = hSource;
		StructForEnergyForceTorqueCompPtr->hDest = hDest;
		StructForEnergyForceTorqueCompPtr->radPtr = this;
		StructForEnergyForceTorqueCompPtr->DestSubdivArray = ActualSubdivisionArray;
		StructForEnergyForceTorqueCompPtr->AutoDestSubdivision = CheckForAutoDestSubdivision(ActualSubdivisionArray);
		radTHandleStructForEnergyForceTorqueComp HandleStructForEnergyForceTorqueComp(StructForEnergyForceTorqueCompPtr);
		StructForEnergyForceTorqueCompPtr = nullptr;  // Ownership transferred to handle

		radTField Field(FieldKey, CompCriterium, HandleStructForEnergyForceTorqueComp);
		DestPtr->EnergyForceTorqueComp(&Field);

		if(Field.HandleEnergyForceTorqueCompData.rep->SomethingIsWrong) return;
		if(SendingIsRequired) Send.OutFieldForceOrTorqueThroughEnergyCompRes(ForceComponID, Field.Force, 'f');
	}
	catch(...)
	{
		if(StructForEnergyForceTorqueCompPtr) delete StructForEnergyForceTorqueCompPtr;  // Clean up if exception before handle ownership transfer
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeFieldTorqueThroughEnergy(int DestElemKey, int SourceElemKey, char* TorqueComponID, int* SubdivArray, long lenSubdivArray, double* TorqueCenPo, long lenTorqueCenPo)
{
	radTStructForEnergyForceTorqueComp* StructForEnergyForceTorqueCompPtr = nullptr;
	try
	{
		radThg hDest;
		if(!ValidateElemKey(DestElemKey, hDest)) return;
		radTg3d* DestPtr = Cast.g3dCast(hDest.rep);
		if(DestPtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		radThg hSource;
		if(!ValidateElemKey(SourceElemKey, hSource)) return;
		radTg3d* SourcePtr = Cast.g3dCast(hSource.rep);
		if(SourcePtr==0) { Send.ErrorMessage("Radia::Error003"); return;}

		if(!ValidateTorqueChar(TorqueComponID)) return;

		radTFieldKey FieldKey;
		FieldKey.Torque_ = 1;

		if((lenSubdivArray != 3) || (SubdivArray[0] < 0) || (SubdivArray[1] < 0) || (SubdivArray[2] < 0))
		{
			Send.ErrorMessage("Radia::Error021"); return;
		}
		//double ActualSubdivisionArray[] = {SubdivArray[0], 1., SubdivArray[1], 1., SubdivArray[2], 1.};
		double ActualSubdivisionArray[] = {(double)SubdivArray[0], 1., (double)SubdivArray[1], 1., (double)SubdivArray[2], 1.}; //OC101015

		StructForEnergyForceTorqueCompPtr = new radTStructForEnergyForceTorqueComp();
		StructForEnergyForceTorqueCompPtr->hSource = hSource;
		StructForEnergyForceTorqueCompPtr->hDest = hDest;
		StructForEnergyForceTorqueCompPtr->radPtr = this;
		StructForEnergyForceTorqueCompPtr->DestSubdivArray = ActualSubdivisionArray;
		StructForEnergyForceTorqueCompPtr->AutoDestSubdivision = CheckForAutoDestSubdivision(ActualSubdivisionArray);
		radTHandleStructForEnergyForceTorqueComp HandleStructForEnergyForceTorqueComp(StructForEnergyForceTorqueCompPtr);
		StructForEnergyForceTorqueCompPtr = nullptr;  // Ownership transferred to handle

		radTField Field(FieldKey, CompCriterium, HandleStructForEnergyForceTorqueComp);
		if(!ValidateVector3d(TorqueCenPo, lenTorqueCenPo, &(Field.P))) return;

		DestPtr->EnergyForceTorqueComp(&Field);

		if(Field.HandleEnergyForceTorqueCompData.rep->SomethingIsWrong) return;
		if(SendingIsRequired) Send.OutFieldForceOrTorqueThroughEnergyCompRes(TorqueComponID, Field.Torque, 't');
	}
	catch(...)
	{
		if(StructForEnergyForceTorqueCompPtr) delete StructForEnergyForceTorqueCompPtr;  // Clean up if exception before handle ownership transfer
		Initialize(); return;
	}
}

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
// ComputeSecondOrderKick all REMOVED (2026-03-22). Use CERN Xsuite/Xtrack.
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

int radTApplication::SubdivideElement_g3d(int ElemKey, double* SubdivArray, long lenSubdivArray, char TypeExtraSpec, double* ExtraSpec, long lenExtraSpec, const char** OptionNames, const char** OptionValues, int OptionCount)
{
	radTCylindricSubdivSpec* pCylindricSubdivSpec = 0;

	try
	{
		radTmhg::iterator iter = GlobalMapOfHandlers.find(ElemKey);
		if(iter == GlobalMapOfHandlers.end()) { Send.ErrorMessage("Radia::Error002"); return 0;}
		radTg3d* g3dPtr = Cast.g3dCast((*iter).second.rep); 
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		if(lenSubdivArray != 6) { Send.ErrorMessage("Radia::Error052"); return 0;}
		if((SubdivArray[0] <= 0) || (SubdivArray[1] <= 0) || (SubdivArray[2] <= 0)
			|| (SubdivArray[3] <= 0) || (SubdivArray[4] <= 0) || (SubdivArray[5] <= 0))
		{
			Send.ErrorMessage("Radia::Error053"); return 0;
		}

		char SubdivisionFrame = 0; // 0- Local; 1- Laboratory, each Group Member separately; 2- Laboratory, all Group as whole;
		char SubdivisionParamCode = 0; // 0- kx,ky,kz are subdiv. numbers; 1- kx,ky,kz are average sizes of pieces;
		char SubdivideCoils = 0; // 0- do not subdivide coils; 1- subdivide coils;

		radTOptionNames OptNam;
		const char** BufNameString = OptionNames;
		const char** BufValString = OptionValues;

		for(int i=0; i<OptionCount; i++)
		{
			if(!strcmp(*BufNameString, OptNam.Frame))
			{
				if(!strcmp(*BufValString, (OptNam.FrameValues)[0])) SubdivisionFrame = 0;
				else if(!strcmp(*BufValString, (OptNam.FrameValues)[1])) SubdivisionFrame = 1;
				else if(!strcmp(*BufValString, (OptNam.FrameValues)[2])) SubdivisionFrame = 2;
				else { Send.ErrorMessage("Radia::Error062"); return 0;}
			}
			else if(!strcmp(*BufNameString, OptNam.SubdParamCode))
			{
				if(!strcmp(*BufValString, (OptNam.SubdParamCodeValues)[0])) SubdivisionParamCode = 0;
				else if(!strcmp(*BufValString, (OptNam.SubdParamCodeValues)[1])) SubdivisionParamCode = 1;
				else { Send.ErrorMessage("Radia::Error062"); return 0;}
			}
			else if(!strcmp(*BufNameString, OptNam.SubdCoils))
			{
				if((!strcmp(*BufValString, (OptNam.SubdCoilsValues)[0])) || (!strcmp(*BufValString, (OptNam.SubdCoilsValues)[2]))) SubdivideCoils = 0;
				else if((!strcmp(*BufValString, (OptNam.SubdCoilsValues)[1])) || (!strcmp(*BufValString, (OptNam.SubdCoilsValues)[3]))) SubdivideCoils = 1;
				else { Send.ErrorMessage("Radia::Error062"); return 0;}
			}
			else { Send.ErrorMessage("Radia::Error062"); return 0;}
			BufNameString++; BufValString++;
		}

		if(SubdivisionParamCode == 0)
		{// Important: setting q=1 if k=1
			const double SubdZeroTol = 1.E-12;
			for(int kk=0; kk<3; kk++)
			{
				int TwoKk = kk*2;
				if(SubdivArray[TwoKk] < 1.) SubdivArray[TwoKk] = 1.;
				if(fabs(SubdivArray[TwoKk]-1.) < SubdZeroTol) SubdivArray[TwoKk+1] = 1.;
			}
		}

		double NewSubdivArray[15];
		for(int ii=0; ii<lenSubdivArray; ii++) NewSubdivArray[ii] = SubdivArray[ii];

		radTSubdivOptions SubdivOptions;
		SubdivOptions.SubdivisionFrame = SubdivisionFrame;
		SubdivOptions.SubdivisionParamCode = SubdivisionParamCode;
		SubdivOptions.SubdivideCoils = SubdivideCoils;
		SubdivOptions.PutNewStuffIntoGenCont = 1;
		SubdivOptions.ReplaceOldStuff = 1;

		char CylindricSubdivision = (TypeExtraSpec == 1)? 1 : 0;
		if(CylindricSubdivision)
		{
			if((ExtraSpec[3]==0.) && (ExtraSpec[4]==0.) && (ExtraSpec[5]==0.)) { Send.ErrorMessage("Radia::Error066"); return 0;}
			if((lenExtraSpec > 6) && (ExtraSpec[9]<=0.)) { Send.ErrorMessage("Radia::Error067"); return 0;}

			pCylindricSubdivSpec = new radTCylindricSubdivSpec(ExtraSpec, lenExtraSpec);
			if(pCylindricSubdivSpec == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

			SubdivOptions.MethForRadialSegmAtEllCylSubd = 0; // Make accessible from outside if necessary
		}

		char SubdivisionByParPlanes = (TypeExtraSpec == 2)? 1 : 0;
		int AmOfDir;
		if(SubdivisionByParPlanes)
		{
			AmOfDir = int((lenExtraSpec + 1.E-06)/3.);
			if(AmOfDir < 1) { Send.ErrorMessage("Radia::Error069"); return 0;}

			for(int i=0; i<AmOfDir; i++)
			{
				int im5 = i*5;
				int im2 = i*2;
				int im3 = i*3;
				NewSubdivArray[im5] = ExtraSpec[im3];
				NewSubdivArray[im5+1] = ExtraSpec[im3+1];
				NewSubdivArray[im5+2] = ExtraSpec[im3+2];

				NewSubdivArray[im5+3] = SubdivArray[im2];
				NewSubdivArray[im5+4] = SubdivArray[im2+1];

				if((NewSubdivArray[im5]==0.) && (NewSubdivArray[im5+1]==0.) && (NewSubdivArray[im5+2]==0.)) { Send.ErrorMessage("Radia::Error069"); return 0;}
				if((NewSubdivArray[im5+3]<=0.) || (NewSubdivArray[im5+4]<=0.)) { Send.ErrorMessage("Radia::Error053"); return 0;}
			}
			//double NewSubdivArray[] = {n1x,n1y,n1z,k1x,q1x, n2x,n2y,n2z,k2x,q2x, n3x,n3y,n3z,k3x,q3x};
		}

		radThg& hgNew = (*iter).second;
		radThg hgOld = hgNew;

		if(CylindricSubdivision)
		{
			if(!g3dPtr->SubdivideItselfByEllipticCylinder(NewSubdivArray, pCylindricSubdivSpec, hgNew, this, &SubdivOptions)) return 0;
		}
		else if(SubdivisionByParPlanes)
		{
			if(!g3dPtr->SubdivideItselfByParPlanes(NewSubdivArray, AmOfDir, hgNew, this, &SubdivOptions)) return 0;
		}
		else if(!g3dPtr->SubdivideItself(NewSubdivArray, hgNew, this, &SubdivOptions)) return 0;

		if(g3dPtr->IsGroupMember) ReplaceInAllGroups(hgOld, hgNew);
		if(SendingIsRequired) Send.Int(ElemKey);

		if(pCylindricSubdivSpec != 0) { delete pCylindricSubdivSpec; pCylindricSubdivSpec = 0;}
		return ElemKey;
	}
	catch(...)
	{
		if(pCylindricSubdivSpec != 0) delete pCylindricSubdivSpec;
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SubdivideElement_g3dByParPlanes(int ElemKey, double* SubdivArray, int AmOfSubdivDirections, const char* LocOrLabFrame)
{
	try
	{
		radTmhg::iterator iter = GlobalMapOfHandlers.find(ElemKey);
		if(iter == GlobalMapOfHandlers.end()) { Send.ErrorMessage("Radia::Error002"); return 0;}
		radTg3d* g3dPtr = Cast.g3dCast((*iter).second.rep); 
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		char SubdivisionFrame = 0; // 0- Local; 1- Laboratory;
		char SubdivisionParamCode = 0; // 0- kx,ky,kz are subdiv. numbers; 1- kx,ky,kz are average sizes of pieces;
		char SubdivideCoils = 0; // 0- do not subdivide coils; 1- subdivide coils;

		if((!strcmp(LocOrLabFrame, "loc")) || (!strcmp(LocOrLabFrame, "Loc")) || (!strcmp(LocOrLabFrame, "LOC"))) SubdivisionFrame = 0;
		else if((!strcmp(LocOrLabFrame, "lab")) || (!strcmp(LocOrLabFrame, "Lab")) || (!strcmp(LocOrLabFrame, "LAB"))) SubdivisionFrame = 1;
		else { Send.ErrorMessage("Radia::Error054"); return 0;}

		if(SubdivisionParamCode == 0)
		{// Important: setting q=1 if k=1
			const double SubdZeroTol = 1.E-12;
			for(int kk=0; kk<AmOfSubdivDirections; kk++)
			{
				int BufInd = 5*kk+4;
				if(fabs(SubdivArray[BufInd] - 1.) < SubdZeroTol) SubdivArray[BufInd] = 1.;
			}
		}

		radTSubdivOptions SubdivOptions;
		SubdivOptions.SubdivisionFrame = SubdivisionFrame;
		SubdivOptions.SubdivisionParamCode = SubdivisionParamCode;
		SubdivOptions.SubdivideCoils = SubdivideCoils;
		SubdivOptions.PutNewStuffIntoGenCont = 1;

		radThg& hgNew = (*iter).second;
		radThg hgOld = hgNew;

		int SubdOK = g3dPtr->SubdivideItselfByParPlanes(SubdivArray, AmOfSubdivDirections, hgNew, this, &SubdivOptions);
		if(!SubdOK) return 0;
		if(g3dPtr->IsGroupMember) ReplaceInAllGroups(hgOld, hgNew);

		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::CutElement_g3d(int ElemKey, double* PointOnPlane, long lenPointOnPlane, double* PlaneNormal, long lenPlaneNormal, const char** OptionNames, const char** OptionValues, int OptionCount)
{
	try
	{
		radTmhg::iterator iter = GlobalMapOfHandlers.find(ElemKey);
		if(iter == GlobalMapOfHandlers.end()) { Send.ErrorMessage("Radia::Error002"); return 0;}
		radTg3d* g3dPtr = Cast.g3dCast((*iter).second.rep); 
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		TVector3d PointOnPlaneVect, PlaneNormalVect;
		if(!ValidateVector3d(PointOnPlane, lenPointOnPlane, &PointOnPlaneVect)) return 0;
		if(!ValidateVector3d(PlaneNormal, lenPlaneNormal, &PlaneNormalVect)) return 0;

		TVector3d PlaneSpecification[] = {PointOnPlaneVect, PlaneNormalVect};

		radTOptionNames OptNam;
		const char* OptNamesToFind[] = {OptNam.Frame, OptNam.SubdCoils};
		char OptValsFoundParsed[] = {0, 0};
		char &SubdivisionFrame = OptValsFoundParsed[0]; // 0- Local; 1- LabTot; 2- Lab;
		char &CutCoils = OptValsFoundParsed[1]; // 0- No; 1- Yes;

		if(!OptNam.findParseOptionValues(OptionNames, OptionValues, OptionCount, OptNamesToFind, 2, OptValsFoundParsed, 0, 0))
		{
			Send.ErrorMessage("Radia::Error062"); return 0;
		}

		//char SubdivisionFrame = 0; // 0- Local; 1- Laboratory;
		//const char** BufNameString = OptionNames;
		//const char** BufValString = OptionValues;
		//for(int i=0; i<OptionCount; i++)
		//{
		//	if(!strcmp(*BufNameString, OptNam.Frame))
		//	{
		//		if(!strcmp(*BufValString, (OptNam.FrameValues)[0])) SubdivisionFrame = 0;
		//		else if(!strcmp(*BufValString, (OptNam.FrameValues)[1])) SubdivisionFrame = 1;
		//		else if(!strcmp(*BufValString, (OptNam.FrameValues)[2])) SubdivisionFrame = 2;
		//		else { Send.ErrorMessage("Radia::Error062"); return 0;}
		//	}
		//	else { Send.ErrorMessage("Radia::Error062"); return 0;}
		//	BufNameString++; BufValString++;
		//}

		radTSubdivOptions SubdivOptions;
		SubdivOptions.SubdivisionFrame = SubdivisionFrame;
		SubdivOptions.SubdivisionParamCode = 0; // Is not used by Cut
		SubdivOptions.SubdivideCoils = CutCoils; //0;
		SubdivOptions.PutNewStuffIntoGenCont = 1;
		SubdivOptions.ReplaceOldStuff = 0;
		SubdivOptions.SeparatePiecesAtCutting = 1;
		SubdivOptions.MapInternalFacesAfterCut = 0;

		radThg hgOld = (*iter).second;
		radThg hgNewTot = hgOld;
		radTPair_int_hg NewLowerPair_int_hg, NewUpperPair_int_hg;

		if(!g3dPtr->CutItself(PlaneSpecification, hgNewTot, NewLowerPair_int_hg, NewUpperPair_int_hg, this, &SubdivOptions)) return 0;

		int LowerElemKey = 0, UpperElemKey = 0;
		if(NewLowerPair_int_hg.Handler_g.rep != 0)
		{
			LowerElemKey = NewLowerPair_int_hg.m;
			if(LowerElemKey == 0) LowerElemKey = AddElementToContainer(NewLowerPair_int_hg.Handler_g);
		}
		if(NewUpperPair_int_hg.Handler_g.rep != 0) 
		{
			UpperElemKey = NewUpperPair_int_hg.m;
			if(UpperElemKey == 0) UpperElemKey = AddElementToContainer(NewUpperPair_int_hg.Handler_g);
		}

		if(SendingIsRequired) 
		{
			//char LowerElemKeyIsNotZero = (LowerElemKey != 0);
			//char UpperElemKeyIsNotZero = (UpperElemKey != 0);
			//int NumberOfElem = (LowerElemKeyIsNotZero && UpperElemKeyIsNotZero)? 2 : 1;
			int NumberOfElem = 0; //OC290908

			int OutArray[2];
			int *pOutArray = OutArray;
			if(LowerElemKey != 0) { *(pOutArray++) = LowerElemKey; NumberOfElem++;}
			if(UpperElemKey != 0) { *pOutArray = UpperElemKey; NumberOfElem++;}

			if(NumberOfElem <= 0) { *pOutArray = ElemKey; NumberOfElem++;} //OC290908

			Send.IntList(OutArray, NumberOfElem);
		}
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::FieldCompMethForSubdividedRecMag(int ElemKey, int InFldCmpMeth, int SubLevel)
{
	radThg hg;
	if(!ValidateElemKey(ElemKey, hg)) return 0;
	radTg3d* g3dP = Cast.g3dCast(hg.rep); if(g3dP==0) { Send.ErrorMessage("Radia::Error003"); return 0;}
	radTGroup* GroupP = Cast.GroupCast(g3dP); if(GroupP==0) { Send.ErrorMessage("Radia::Error004"); return 0;}
	radTSubdividedRecMag* SubdividedRecMagP = Cast.SubdividedRecMagCast(GroupP); 
	if(SubdividedRecMagP==0) { Send.ErrorMessage("Radia::Error037"); return 0;}

	int Probe = SubdividedRecMagP->SetupFldCmpData((short)InFldCmpMeth, SubLevel);
	if(Probe==0) { Send.ErrorMessage("Radia::Error040"); return 0;}
	else if(Probe==-38) { Send.ErrorMessage("Radia::Error038"); return 0;}

	if(SendingIsRequired) Send.Int(ElemKey);
	return ElemKey;
}

//-------------------------------------------------------------------------

int radTApplication::SetLocMgnInSbdRecMag(int ElemKey, TVector3d* ArrayOfVectIndx, TVector3d* ArrayOfMagn, int Len)
{
	radThg hg;
	if(!ValidateElemKey(ElemKey, hg)) return 0;
	radTg3d* g3dP = Cast.g3dCast(hg.rep); if(g3dP==0) { Send.ErrorMessage("Radia::Error003"); return 0;}
	radTGroup* GroupP = Cast.GroupCast(g3dP); if(GroupP==0) { Send.ErrorMessage("Radia::Error004"); return 0;}
	radTSubdividedRecMag* SubdividedRecMagP = Cast.SubdividedRecMagCast(GroupP); 
	if(SubdividedRecMagP==0) { Send.ErrorMessage("Radia::Error037"); return 0;}
	
	for(int i=0; i<Len; i++)
	{
		int ix = int(ArrayOfVectIndx[i].x) - 1;
		if((ix >= int(SubdividedRecMagP->kx)) || (ix < 0)) { Send.ErrorMessage("Radia::Error039"); return 0;}
		int iy = int(ArrayOfVectIndx[i].y) - 1;
		if((iy >= int(SubdividedRecMagP->ky)) || (iy < 0)) { Send.ErrorMessage("Radia::Error039"); return 0;}
		int iz = int(ArrayOfVectIndx[i].z) - 1;
		if((iz >= int(SubdividedRecMagP->kz)) || (iz < 0)) { Send.ErrorMessage("Radia::Error039"); return 0;}

		int SubElemNo = (ix*int(SubdividedRecMagP->ky) + iy)*int(SubdividedRecMagP->kz) + iz;

		SubdividedRecMagP->SetLocalMagn(SubElemNo, ArrayOfMagn[i]);
	}

	if(SendingIsRequired) Send.Int(ElemKey);
	return ElemKey;
}

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
		CopyDrawAttr(ElemKey, NewElemKey);

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
			CopyDrawAttr(ElemKey, NewElemKey);

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

		MapOfDrawAttr.erase(MapOfDrawAttr.begin(), MapOfDrawAttr.end());

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

int radTApplication::ComputeGeometricalLimits(int ElemKey)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); 
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		double LimArr[6];
		g3dPtr->Limits(0, LimArr);

		if(SendingIsRequired) Send.DoubleList(LimArr, 6);
		return 1;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------
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
