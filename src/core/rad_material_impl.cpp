/*-------------------------------------------------------------------------
*
* File name:      radapl2.cpp
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
#include "rad_interaction.h"
#include "rad_material_def.h"
#include "rad_relaxation_methods.h"
#include "rad_particle_trajectory.h"
#include "rad_operation_names.h"
#include "rad_material_aux.h"
#include "rad_point_classify.h"
#include "gmvbstr.h"

#include <math.h>
#include <string.h>
#include <cstring>
#include <cstdio>
#include <chrono>
#include <cmath>

#include "rad_parallel.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

int radTApplication::SetLinearMaterial(double* KsiArray, long lenKsiArray, double* RemMagnArray, long lenRemMagnArray)
{
	radTMaterial* MaterPtr = nullptr;
	try
	{
		if(lenKsiArray != 2)
		{
			Send.ErrorMessage("Radia::Error022"); return 0;
		}
		TVector3d RemMagnVect;
		char EasyAxisDefined;
		if(lenRemMagnArray == 3)
		{
			EasyAxisDefined = 1;

			if(!ValidateVector3d(RemMagnArray, lenRemMagnArray, &RemMagnVect)) return 0;
			if((RemMagnVect.x==0) && (RemMagnVect.y==0) && (RemMagnVect.z==0) && (KsiArray[0]!=KsiArray[1]))
			{ Send.ErrorMessage("Radia::Error023"); return 0;}
		}
		else if(lenRemMagnArray == 1)
		{
			EasyAxisDefined = 0;
			RemMagnVect.x = *RemMagnArray;
		}

		MaterPtr = new radTLinearAnisotropMaterial(KsiArray, RemMagnVect, EasyAxisDefined);
		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		radThg hg(MaterPtr);
		MaterPtr = nullptr;  // Ownership transferred to radThg
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		if(MaterPtr) delete MaterPtr;  // Clean up if exception before radThg ownership transfer
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetLinearIsotropicMaterial(double Ksi)
{
	radTLinearIsotropMaterial* MaterPtr = nullptr;
	try
	{
		MaterPtr = new radTLinearIsotropMaterial(Ksi);
		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		radThg hg(MaterPtr);
		MaterPtr = nullptr;  // Ownership transferred to radThg
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		if(MaterPtr) delete MaterPtr;  // Clean up if exception before radThg ownership transfer
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetLinearAnisotropicMaterial(double* KsiArray, long lenKsiArray, double* EasyAxisArray, long lenEasyAxisArray)
{
	radTLinearAnisotropMaterial* MaterPtr = nullptr;
	try
	{
		if(lenKsiArray != 2)
		{
			Send.ErrorMessage("Radia::Error022"); return 0;
		}
		if(lenEasyAxisArray != 3)
		{
			Send.ErrorMessage("Radia::Error102"); return 0;  // Easy axis must have 3 components
		}

		TVector3d EasyAxisVect;
		if(!ValidateVector3d(EasyAxisArray, lenEasyAxisArray, &EasyAxisVect)) return 0;

		// Normalize easy axis
		double AbsEasyAxis = sqrt(EasyAxisVect.x*EasyAxisVect.x + EasyAxisVect.y*EasyAxisVect.y + EasyAxisVect.z*EasyAxisVect.z);
		if(AbsEasyAxis < 1.E-10)
		{
			Send.ErrorMessage("Radia::Error103"); return 0;  // Easy axis cannot be zero vector
		}
		TVector3d RemMagnVect = (1.0/AbsEasyAxis) * EasyAxisVect;  // Use normalized easy axis as RemMagn

		MaterPtr = new radTLinearAnisotropMaterial(KsiArray, RemMagnVect, 1);  // EasyAxisDefined = 1
		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		radThg hg(MaterPtr);
		MaterPtr = nullptr;  // Ownership transferred to radThg
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		if(MaterPtr) delete MaterPtr;  // Clean up if exception before radThg ownership transfer
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetPermanentMagnet(double Br, double Hc, double* MagAxisArray, long lenMagAxisArray)
{
	radTPermanentMagnet* MaterPtr = nullptr;
	try
	{
		if(lenMagAxisArray != 3)
		{
			Send.ErrorMessage("Radia::Error104"); return 0;  // Magnetization axis must have 3 components
		}

		TVector3d MagAxisVect;
		if(!ValidateVector3d(MagAxisArray, lenMagAxisArray, &MagAxisVect)) return 0;

		double AbsMagAxis = sqrt(MagAxisVect.x*MagAxisVect.x + MagAxisVect.y*MagAxisVect.y + MagAxisVect.z*MagAxisVect.z);
		if(AbsMagAxis < 1.E-10)
		{
			Send.ErrorMessage("Radia::Error105"); return 0;  // Magnetization axis cannot be zero vector
		}

		if(Br <= 0 || Hc <= 0)
		{
			Send.ErrorMessage("Radia::Error106"); return 0;  // Br and Hc must be positive
		}

		MaterPtr = new radTPermanentMagnet(Br, Hc, MagAxisVect);
		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		radThg hg(MaterPtr);
		MaterPtr = nullptr;  // Ownership transferred to radThg
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		if(MaterPtr) delete MaterPtr;  // Clean up if exception before radThg ownership transfer
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetMagFixed(double* MagnArray, long lenMagnArray)
{
	radTMagFixed* MaterPtr = nullptr;
	try
	{
		if(lenMagnArray != 3)
		{
			Send.ErrorMessage("Radia::Error104"); return 0;  // Magnetization must have 3 components
		}

		TVector3d MagnVect;
		if(!ValidateVector3d(MagnArray, lenMagnArray, &MagnVect)) return 0;

		MaterPtr = new radTMagFixed(MagnVect);
		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		radThg hg(MaterPtr);
		MaterPtr = nullptr;  // Ownership transferred to radThg
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		if(MaterPtr) delete MaterPtr;
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetMagLinear(double Br, double Hc, double* MagAxisArray, long lenMagAxisArray)
{
	radTMagLinear* MaterPtr = nullptr;
	try
	{
		if(lenMagAxisArray != 3)
		{
			Send.ErrorMessage("Radia::Error104"); return 0;  // Magnetization axis must have 3 components
		}

		TVector3d MagAxisVect;
		if(!ValidateVector3d(MagAxisArray, lenMagAxisArray, &MagAxisVect)) return 0;

		double AbsMagAxis = sqrt(MagAxisVect.x*MagAxisVect.x + MagAxisVect.y*MagAxisVect.y + MagAxisVect.z*MagAxisVect.z);
		if(AbsMagAxis < 1.E-10)
		{
			Send.ErrorMessage("Radia::Error105"); return 0;  // Magnetization axis cannot be zero vector
		}

		if(Br <= 0 || Hc <= 0)
		{
			Send.ErrorMessage("Radia::Error106"); return 0;  // Br and Hc must be positive
		}

		MaterPtr = new radTMagLinear(Br, Hc, MagAxisVect);
		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		radThg hg(MaterPtr);
		MaterPtr = nullptr;  // Ownership transferred to radThg
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		if(MaterPtr) delete MaterPtr;
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetMagCurve(double* pCurveData, int numPoints, double* MagAxisArray, long lenMagAxisArray)
{
	radTMagCurve* MaterPtr = nullptr;
	try
	{
		if(lenMagAxisArray != 3)
		{
			Send.ErrorMessage("Radia::Error104"); return 0;
		}

		TVector3d MagAxisVect;
		if(!ValidateVector3d(MagAxisArray, lenMagAxisArray, &MagAxisVect)) return 0;

		double AbsMagAxis = sqrt(MagAxisVect.x*MagAxisVect.x + MagAxisVect.y*MagAxisVect.y + MagAxisVect.z*MagAxisVect.z);
		if(AbsMagAxis < 1.E-10)
		{
			Send.ErrorMessage("Radia::Error105"); return 0;
		}

		if(numPoints < 2)
		{
			Send.ErrorMessage("Radia::Error024"); return 0;  // Need at least 2 points
		}

		// Convert flat array to TVector2d array
		std::vector<TVector2d> vCurve(numPoints);
		for(int i = 0; i < numPoints; i++)
		{
			vCurve[i].x = pCurveData[2*i];     // H value
			vCurve[i].y = pCurveData[2*i + 1]; // B value
		}

		MaterPtr = new radTMagCurve(vCurve.data(), numPoints, MagAxisVect);
		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		radThg hg(MaterPtr);
		MaterPtr = nullptr;
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		if(MaterPtr) delete MaterPtr;
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetNonlinearIsotropMaterial(double* Ms, long lenMs, double* ks, long len_ks)
{
	radTNonlinearIsotropMaterial* MaterPtr = nullptr;
	try
	{
		if((lenMs != len_ks) || (lenMs > 3)) { Send.ErrorMessage("Radia::Error024"); return 0;}

		MaterPtr = new radTNonlinearIsotropMaterial(Ms, ks, (int)len_ks);
		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hg(MaterPtr);
		MaterPtr = nullptr;  // Ownership transferred to radThg
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		if(MaterPtr) delete MaterPtr;  // Clean up if exception before radThg ownership transfer
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetNonlinearIsotropMaterial(TVector2d* ArrayHB, int LenArrayHB)
{
	radTNonlinearIsotropMaterial* MaterPtr = nullptr;
	try
	{
		if(!ValidateIsotropMaterDescrByPoints(ArrayHB, LenArrayHB)) return 0;

		MaterPtr = new radTNonlinearIsotropMaterial(ArrayHB, LenArrayHB);
		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hg(MaterPtr);
		MaterPtr = nullptr;  // Ownership transferred to radThg
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		if(MaterPtr) delete MaterPtr;  // Clean up if exception before radThg ownership transfer
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::ValidateIsotropMaterDescrByPoints(TVector2d* ArrayHB, int LenArrayHB)
{
	try
	{
		if((ArrayHB->x < 0.) || (ArrayHB->y < 0.)) { Send.ErrorMessage("Radia::Error071"); return 0;}

		TVector2d* tArrayHB = ArrayHB;
		double Hprev = -1, Bprev = -1;
		for(int i=0; i<LenArrayHB; i++)
		{
			if(tArrayHB->x < Hprev) { Send.ErrorMessage("Radia::Error071"); return 0;}
			if(tArrayHB->y < 0.95*Bprev) { Send.WarningMessage("Radia::Warning014");}

			Hprev = tArrayHB->x; Bprev = tArrayHB->y;
			tArrayHB++;
		}
		return 1;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetNonlinearLaminatedMaterial(TVector2d* ArrayOfPoints2d, int lenArrayOfPoints2d, double PackFactor, double* dN)
{
	if(lenArrayOfPoints2d <= 0) { Send.ErrorMessage("Radia::Error000"); return 0;}
	if((PackFactor <= 0) || (PackFactor > 1)) { Send.ErrorMessage("Radia::Error074"); return 0;}

	radTMaterial *MaterPtr = nullptr;

	try
	{
		if(lenArrayOfPoints2d <= 3)
		{
			double Ms[] = {0,0,0};
			double Ks[] = {0,0,0};
			int lenMs = lenArrayOfPoints2d;
			for(int i=0; i<lenMs; i++)
			{
				Ks[i] = ArrayOfPoints2d[i].x;
				Ms[i] = ArrayOfPoints2d[i].y;
			}

			if((PackFactor <= 0.) || (PackFactor >= 1.)) MaterPtr = new radTNonlinearIsotropMaterial(Ms, Ks, lenMs);
			else MaterPtr = new radTNonlinearLaminatedMaterial(Ms, Ks, lenMs, PackFactor, dN);
		}
		else
		{
			if(!ValidateIsotropMaterDescrByPoints(ArrayOfPoints2d, lenArrayOfPoints2d)) return 0;

			if((PackFactor <= 0.) || (PackFactor >= 1.)) MaterPtr = new radTNonlinearIsotropMaterial(ArrayOfPoints2d, lenArrayOfPoints2d);
			else MaterPtr = new radTNonlinearLaminatedMaterial(ArrayOfPoints2d, lenArrayOfPoints2d, PackFactor, dN);
		}

		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error000"); return 0;}
		radThg hg(MaterPtr);
		MaterPtr = nullptr;  // Ownership transferred to radThg
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		if(MaterPtr) delete MaterPtr;  // Clean up if exception before radThg ownership transfer
		Send.ErrorMessage("Radia::Error075");
		return 0;
	}
}

//-------------------------------------------------------------------------

//int radTApplication::SetNonlinearAnisotropMaterial(double** Ksi, double** Ms, double* Hc, char* DependenceIsNonlinear)
int radTApplication::SetNonlinearAnisotropMaterial(double** Ksi, double** Ms, double* Hc, int lenHc, char* DependenceIsNonlinear)
{
	radTMaterial* MaterPtr = nullptr;
	try
	{
		char MaterialIsIsotropic = 1;
		if(DependenceIsNonlinear[0] || DependenceIsNonlinear[1])
		{
			if(!(DependenceIsNonlinear[0] && DependenceIsNonlinear[1])) MaterialIsIsotropic = 0;
			else if((lenHc == 2) && (!(Hc[0]==0. && Hc[1]==0.))) MaterialIsIsotropic = 0;
			else if((lenHc == 4) && ((Hc[0]!=0.) || (Hc[1]!=0.) || (Hc[2]!=0.) || (Hc[3]!=0.))) MaterialIsIsotropic = 0;
			else
			{
				double *KsiPar = Ksi[0], *KsiPer = Ksi[1], *MsPar = Ms[0], *MsPer = Ms[1];
				for(int i=0; i<3; i++)
				{
					if(!((*(KsiPar++)==*(KsiPer++)) && (*(MsPar++)==*(MsPer++)))) { MaterialIsIsotropic = 0; break;}
				}
				if((*KsiPar != 0.) || (*KsiPer != 0.)) MaterialIsIsotropic = 0;
			}
			if(MaterialIsIsotropic) MaterPtr = new radTNonlinearIsotropMaterial(Ms[0], Ksi[0], 3);
			else 
			{
				double Hci[4];
				if(lenHc <= 2) 
				{
					for(int i=0; i<lenHc; i++) Hci[i] = Hc[0];
				}
				else
				{
					for(int i=0; i<lenHc; i++) Hci[i] = Hc[i];
				}
				//MaterPtr = new radTNonlinearAnisotropMaterial(Ksi, Ms, Hc, DependenceIsNonlinear);
				MaterPtr = new radTNonlinearAnisotropMaterial(Ksi, Ms, Hci, DependenceIsNonlinear);
			}
		}
		else { Send.ErrorMessage("Radia::Error060"); return 0;}

		if(DependenceIsNonlinear[1] && (Hc[1] != 0.) && (lenHc == 2)) { Send.ErrorMessage("Radia::Error061"); return 0;}

		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hg(MaterPtr);
		MaterPtr = nullptr;  // Ownership transferred to radThg
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		if(MaterPtr) delete MaterPtr;  // Clean up if exception before radThg ownership transfer
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetNonlinearAnisotropMaterial0(double* pDataPar, int lenDataPar, double* pDataPer, int lenDataPer)
{
	radTMaterial* MaterPtr = nullptr;
	try
	{
		bool MaterialIsIsotropic = true;
		char DependenceIsNonlinear[] = {(lenDataPar > 1), (lenDataPer > 1)};

		double KsiPar[4], KsiPer[4], MsPar[3], MsPer[3], Hci[4];
		double *tKsiPar = KsiPar, *tKsiPer = KsiPer, *tMsPar = MsPar, *tMsPer = MsPer, *tHci = Hci;
		for(int j=0; j<3; j++)
		{
			*(tKsiPar++) = 0;
			*(tKsiPer++) = 0;
			*(tMsPar++) = 0;
			*(tMsPer++) = 0;
			*(tHci++) = 0;
		}
		*tKsiPar = 0; *tKsiPer = 0; *tHci = 0;

		double *Ksi[] = {KsiPar, KsiPer}, *Ms[] = {MsPar, MsPer};

		if(lenDataPar == 11) //{ksi1,ms1,hc1,ksi2,ms2,hc2,ksi3,ms3,hc3,ksi0,hc0}
		{
			double *tKsi = Ksi[0], *tMs = Ms[0], *tHci = Hci;
			double *tDataPar = pDataPar;
			for(int i=0; i<3; i++) 
			{
				*(tKsi++) = *(tDataPar++);
				*(tMs++) = *(tDataPar++);
				*(tHci++) = *(tDataPar++);
			}
			*tKsi = *(tDataPar++);
			*tHci = *tDataPar;
		}
		else if(lenDataPar == 8) //{ksi1,ms1,ksi2,ms2,ksi3,ms3,ksi0,hc}
		{
			double *tKsi = Ksi[0], *tMs = Ms[0];
			double *tDataPar = pDataPar;
			for(int i=0; i<3; i++) 
			{
				*(tKsi++) = *(tDataPar++);
				*(tMs++) = *(tDataPar++);
			}
			*tKsi = *(tDataPar++);
			Hci[0] = Hci[1] = Hci[2] = Hci[3] = *tDataPar;
		}
		else if(lenDataPar == 1) //{ksi0}
		{
			double *tKsi = Ksi[0]; //, *tMs = Ms[0];
			double *tDataPar = pDataPar;
			//for(int i=0; i<3; i++) 
			//{
			//	*(tKsi++) = 0;
			//	*(tMs++) = 0;
			//}
			*tKsi = *tDataPar;
			Hci[0] = Hci[1] = Hci[2] = Hci[3] = 0;
		}

		if(lenDataPer == 7) //{ksi1,ms1,ksi2,ms2,ksi3,ms3,ksi0}
		{
			double *tKsi = Ksi[1], *tMs = Ms[1];
			double *tDataPer = pDataPer;
			for(int i=0; i<3; i++) 
			{
				*(tKsi++) = *(tDataPer++);
				*(tMs++) = *(tDataPer++);
			}
			*tKsi = *tDataPer;
		}
		else if(lenDataPer == 1) //{ksi0}
		{
			double *tKsi = Ksi[1]; //, *tMs = Ms[1];
			double *tDataPer = pDataPer;
			//for(int i=0; i<3; i++) 
			//{
			//	*(tKsi++) = 0;
			//	*(tMs++) = 0;
			//}
			*tKsi = *tDataPer;
		}

		if(DependenceIsNonlinear[0] || DependenceIsNonlinear[1])
		{
			if(!(DependenceIsNonlinear[0] && DependenceIsNonlinear[1])) MaterialIsIsotropic = false;
			else if(!(Hci[0]==0. && Hci[1]==0. && Hci[2]==0. && Hci[3]==0.)) MaterialIsIsotropic = false;
			else
			{
				double *KsiPar = Ksi[0], *KsiPer = Ksi[1], *MsPar = Ms[0], *MsPer = Ms[1];
				for(int i=0; i<3; i++)
				{
					if(!((*(KsiPar++)==*(KsiPer++)) && (*(MsPar++)==*(MsPer++)))) { MaterialIsIsotropic = false; break;}
				}
				if((*KsiPar != 0.) || (*KsiPer != 0.)) MaterialIsIsotropic = false;
			}
			if(MaterialIsIsotropic) MaterPtr = new radTNonlinearIsotropMaterial(Ms[0], Ksi[0], 3);
			else 
			{
				MaterPtr = new radTNonlinearAnisotropMaterial(Ksi, Ms, Hci, DependenceIsNonlinear);
			}
		}
		else { Send.ErrorMessage("Radia::Error060"); return 0;}

		//if(DependenceIsNonlinear[1] && (Hc[1] != 0.)) { Send.ErrorMessage("Radia::Error061"); return 0;}

		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hg(MaterPtr);
		MaterPtr = nullptr;  // Ownership transferred to radThg
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		if(MaterPtr) delete MaterPtr;  // Clean up if exception before radThg ownership transfer
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::ApplyMaterial(int g3dElemKey, int MaterElemKey)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(g3dElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); 
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		radTg3dRelax* g3dRelaxPtr = Cast.g3dRelaxCast(g3dPtr); 
		if(g3dRelaxPtr==0) 
		{
			radTGroup* GroupPtr = Cast.GroupCast(g3dPtr);
			if(GroupPtr==0) { Send.ErrorMessage("Radia::Error015"); return 0;}
		}

		if(!ValidateElemKey(MaterElemKey, hg)) return 0;
		radTMaterial* MaterPtr = Cast.MaterCast(hg.rep);
		if(MaterPtr==0) { Send.ErrorMessage("Radia::Error016"); return 0;}

		char PutNewStuffIntoGenCont = 1;
		if(!MaterPtr->EasyAxisDefined) 
		{
			if(!MaterPtr->DuplicateItself(hg, this, PutNewStuffIntoGenCont)) return 0;
			AddElementToContainer(hg); // Maybe not necessary
		}

		if(!g3dPtr->SetMaterial(hg, this)) return 0;

		if(SendingIsRequired) Send.Int(g3dElemKey);
		return g3dElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------
/**
void radTApplication::OutMagnetizCompRes(char* MagnChar, TVector3d& M_vect)
{
	char* BufChar = MagnChar;
	char* EqEmptyStr = "MxMyMz";

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			char* BufChar_pl_1 = BufChar+1;
			if((((*BufChar == 'M') || (*BufChar == 'm')) && 
			   (*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z')) ||
			   (*BufChar == 'X') || (*BufChar == 'x') ||
			   (*BufChar == 'Y') || (*BufChar == 'y') ||
			   (*BufChar == 'Z') || (*BufChar == 'z')) ItemCount++;
			BufChar++;
		}
		BufChar = MagnChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 3;
	}

	if(ItemCount > 1) Send.InitOutList(ItemCount);

	while (*BufChar != '\0') 
	{
		if((*(BufChar)=='M') || (*(BufChar)=='m'))
		{
			char* BufChar_pl_1 = BufChar+1;
			if((*(BufChar_pl_1)!='x') && (*(BufChar_pl_1)!='X') &&
			   (*(BufChar_pl_1)!='y') && (*(BufChar_pl_1)!='Y') &&
			   (*(BufChar_pl_1)!='z') && (*(BufChar_pl_1)!='Z')) Send.Vector3d(&M_vect);
		}
		else if((*(BufChar)=='X') || (*(BufChar)=='x')) Send.Double(M_vect.x);
		else if((*(BufChar)=='Y') || (*(BufChar)=='y')) Send.Double(M_vect.y);
		else if((*(BufChar)=='Z') || (*(BufChar)=='z')) Send.Double(M_vect.z);

		BufChar++;
	}
}
**/
//-------------------------------------------------------------------------

void radTApplication::ComputeMvsH(int g3dRelaxOrMaterElemKey, char* MagnChar, double* H, long lenH)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(g3dRelaxOrMaterElemKey, hg)) return;

		radTMaterial* MaterPtr = nullptr;

		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr!=nullptr)
		{
			radTg3dRelax* g3dRelaxPtr = Cast.g3dRelaxCast(g3dPtr);
			if(g3dRelaxPtr!=nullptr)
			{
				MaterPtr = static_cast<radTMaterial*>(g3dRelaxPtr->MaterHandle.rep);
				if(MaterPtr==nullptr) { Send.ErrorMessage("Radia::Error027"); return;}
			}
			else
			{
				radTGroup* GroupPtr = Cast.GroupCast(g3dPtr);
				if(GroupPtr!=nullptr)
				{
					radTg3dRelax* g3dSubdRelaxPtr = nullptr;

					radTSubdividedRecMag* SubdividedRecMagPtr = Cast.SubdividedRecMagCast(GroupPtr);
					if(SubdividedRecMagPtr!=nullptr) g3dSubdRelaxPtr = static_cast<radTg3dRelax*>(SubdividedRecMagPtr);
					else
					{
						radTSubdividedExtrPolygon* SubdividedExtrPolygonPtr = Cast.SubdExtrPolygonCastFromGroup(GroupPtr);
						if(SubdividedExtrPolygonPtr!=nullptr) g3dSubdRelaxPtr = static_cast<radTg3dRelax*>(SubdividedExtrPolygonPtr);
						else
						{
							radTSubdividedPolyhedron* SubdividedPolyhedronPtr = Cast.SubdPolyhedronCastFromGroup(GroupPtr);
							if(SubdividedPolyhedronPtr!=nullptr) g3dSubdRelaxPtr = static_cast<radTg3dRelax*>(SubdividedPolyhedronPtr);
						}
					}
					if(g3dSubdRelaxPtr!=nullptr)
					{
						MaterPtr = static_cast<radTMaterial*>(g3dSubdRelaxPtr->MaterHandle.rep);
						if(MaterPtr==nullptr) { Send.ErrorMessage("Radia::Error027"); return;}
					}
				}
			}
		}
		if(MaterPtr==nullptr)
		{
			MaterPtr = Cast.MaterCast(hg.rep);
			if(MaterPtr==nullptr) { Send.ErrorMessage("Radia::Error025"); return;}
		}

		if(!ValidateMagnChar(MagnChar)) return;
		TVector3d H_vect;
		if(!ValidateVector3d(H, lenH, &H_vect)) return;

		TVector3d M_vect = MaterPtr->M(H_vect);
		if(SendingIsRequired) Send.OutMagnetizCompRes(MagnChar, M_vect);
	}
	catch(...)
	{
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------

int radTApplication::MatHysSaveState(int MaterElemKey, double* pState, int* pLen)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(MaterElemKey, hg)) return -1;

		radTMaterial* MaterPtr = Cast.MaterCast(hg.rep);
		if(MaterPtr==nullptr) { Send.ErrorMessage("Radia::Error025"); return -1; }

		radTHysteresisMaterial* HysPtr = dynamic_cast<radTHysteresisMaterial*>(MaterPtr);
		if(HysPtr==nullptr) { Send.ErrorMessage("Radia::Error025"); return -1; }

		*pLen = HysPtr->GetStateSize();
		if(pState != nullptr) HysPtr->SaveStateToArray(pState);
		return 0;
	}
	catch(...) { Initialize(); return -1; }
}

int radTApplication::MatHysRestoreState(int MaterElemKey, const double* pState, int Len)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(MaterElemKey, hg)) return -1;

		radTMaterial* MaterPtr = Cast.MaterCast(hg.rep);
		if(MaterPtr==nullptr) { Send.ErrorMessage("Radia::Error025"); return -1; }

		radTHysteresisMaterial* HysPtr = dynamic_cast<radTHysteresisMaterial*>(MaterPtr);
		if(HysPtr==nullptr) { Send.ErrorMessage("Radia::Error025"); return -1; }

		if(Len != HysPtr->GetStateSize()) { Send.ErrorMessage("Radia::Error000"); return -1; }
		HysPtr->RestoreStateFromArray(pState);
		return 0;
	}
	catch(...) { Initialize(); return -1; }
}

int radTApplication::MatHysCommitState(int MaterElemKey)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(MaterElemKey, hg)) return -1;

		radTMaterial* MaterPtr = Cast.MaterCast(hg.rep);
		if(MaterPtr==nullptr) { Send.ErrorMessage("Radia::Error025"); return -1; }

		radTHysteresisMaterial* HysPtr = dynamic_cast<radTHysteresisMaterial*>(MaterPtr);
		if(HysPtr==nullptr) { Send.ErrorMessage("Radia::Error025"); return -1; }

		HysPtr->CommitState();
		return 0;
	}
	catch(...) { Initialize(); return -1; }
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//void radTApplication::DumpElem(int ElemKey)
void radTApplication::DumpElem(int* arKeys, int nElem, const char* strFormat, bool arKeysAllocInMathLink)
{
	try
	{
		if((strFormat == 0) || (strcmp(strFormat, "asc") == 0) || (strcmp(strFormat, "ascii") == 0) || (strcmp(strFormat, "ASC") == 0) || (strcmp(strFormat, "ASCII") == 0))
		{
			ostringstream OutDumpStream;
			int nElem_mi_1 = nElem - 1;
			for(int i=0; i<nElem; i++)
			{
				int elemKey = arKeys[i];
				radThg hg;
				if(!ValidateElemKey(elemKey, hg)) return;
				
				OutDumpStream << "Index " << elemKey << ": ";
				(hg.rep)->Dump(OutDumpStream);

				if(i < nElem_mi_1) OutDumpStream << endl;
			}
			
			OutDumpStream << ends;
			Send.String(OutDumpStream.str().c_str());
		}
		else if((strcmp(strFormat, "bin") == 0) || (strcmp(strFormat, "binary") == 0) || (strcmp(strFormat, "BIN") == 0) || (strcmp(strFormat, "BINARY") == 0))
		{
			//CAuxBinStr oStr;
			CAuxBinStrVect oStr;

			oStr << (char)arKeysAllocInMathLink; //indicates whether at parsing a list of elements should be expected

			oStr << nElem; //number of "directly listed" elements, to start with
			for(int j=0; j<nElem; j++) //first writing keys of directly-dumped elements 
			{
				oStr << arKeys[j];
			}
			oStr << nElem; //number of elements (it will be corrected later with actual number of elements)

			//radTmhg locElMap;
			//int elemCount=0;
			vector<int> vElemKeysOut;
			for(int i=0; i<nElem; i++)
			{
				int elemKey = arKeys[i];

				int indExist = CAuxParse::FindElemInd(elemKey, vElemKeysOut);
				if(indExist < 0)
				{//to avoid duplication of objects in output byte string
					radThg hg;
					if(!ValidateElemKey(elemKey, hg)) return;
					//(hg.rep)->DumpBin(oStr, locElMap, hg);
					(hg.rep)->DumpBin(oStr, vElemKeysOut, GlobalMapOfHandlers, GlobalUniqueMapKey, elemKey);
				}
			}

			//int elemCount = (int)locElMap.size();
			//oStr.setFromPos(0, elemCount);
			int elemCount = (int)vElemKeysOut.size();
			//oStr.setFromPos((long)((nElem + 1)*(sizeof(int))), elemCount);
			oStr.setFromPos((long)((nElem + 1)*(sizeof(int)) + 1), elemCount); //OC060713

			//Saving "Drawing Attributes" of objects
			long drwAttrOfst = oStr.getCurOfst();
			int nDrwAttrFound = 0;
			oStr << nDrwAttrFound;
			for(int j=0; j<elemCount; j++)
			{
				int elemKey = vElemKeysOut[j];
				radTMapOfDrawAttr::const_iterator itDrw = MapOfDrawAttr.find(elemKey);
				if(itDrw != MapOfDrawAttr.end())
				{
					const radTDrawAttr &drwAttr = itDrw->second;
					oStr << elemKey;
					//Members of radTDrawAttr
					//double Red, Green, Blue; 
					oStr << drwAttr.RGB_col.Red << drwAttr.RGB_col.Green << drwAttr.RGB_col.Blue; 
					//double LineThickness;
					oStr << drwAttr.LineThickness;
					nDrwAttrFound++;
				}
			}
			if(nDrwAttrFound > 0) oStr.setFromPos(drwAttrOfst, nDrwAttrFound);

			Send.ByteString(reinterpret_cast<const unsigned char*>(oStr.data()), (long)oStr.size());
		}
		else 
		{
			Send.ErrorMessage("Radia::Error000");
		}
	}
	catch(...)
	{
		Initialize();
	}
}

//-------------------------------------------------------------------------

int radTApplication::DumpElemParse(const unsigned char *bstr, int bstrLen)
{
	if((bstr == 0) || (bstrLen <= 0)) { Send.ErrorMessage("Radia::Error000"); return 0;}
	int *arDirElemOldKeys = 0;
	try
	{
		CAuxBinStrVect inStr(bstr, bstrLen);

		char listIsExpectedInOutput = 0;
		inStr >> listIsExpectedInOutput; //indicates whether at parsing a list of elements should be expected

		int nElemDir = 0;
		inStr >> nElemDir;

		std::vector<int> vArDirElemOldKeys(nElemDir);
		arDirElemOldKeys = vArDirElemOldKeys.data();
		int *t_arDirElemOldKeys = arDirElemOldKeys;
		for(int j=0; j<nElemDir; j++)
		{
			inStr >> (*(t_arDirElemOldKeys++));
		}

		int nElemTot = 0;
		inStr >> nElemTot;

		int oldKey;
		char cType1, cType2, cType3, cType4, cType5;
		map<int, int> mKeysOldNew;
		vector<int> vElemKeys;
		
		radTg3d g3d;
		radTArcCur arcCur;
		radTFlmLinCur flmCur;
		radTBackgroundFieldSource bkgFldSrc;
		radTGroup grp;
		radTg3dRelax g3dRelax;
		radTRecMag recMag;
		radTSubdividedRecMag sbdRecMag;
		radTExtrPolygon extPgn;
		radTSubdividedExtrPolygon sbdExtPgn;
		radTPolyhedron polyhdr;
		radTSubdividedPolyhedron sbdPolyhdr;
		radTrans tr;
		radTMaterial mat;
		radTLinearAnisotropMaterial matLinAniso;
		radTLinearIsotropMaterial matLinIso;
		radTNonlinearIsotropMaterial matNonLinIso;
		radTNonlinearAnisotropMaterial matNonLinAniso;
		radTNonlinearLaminatedMaterial matNonLinLam;
		radTInteraction intrc;

		for(int i=0; i<nElemTot; i++)
		{
			inStr >> oldKey;

			inStr >> cType1;
			inStr >> cType2;
			inStr >> cType3;
			inStr >> cType4;
			inStr >> cType5;

			radThg hg;

			if(cType1 == g3d.Type_g())
			{
				if(cType2 == g3dRelax.Type_g3d())
				{
					if(cType3 == recMag.Type_g3dRelax())
					{//Instantiate RecMag
						hg = radThg(new radTRecMag(inStr, mKeysOldNew, GlobalMapOfHandlers));
					}
					else if(cType3 == extPgn.Type_g3dRelax())
					{//Instantiate ExtrPolygon
						hg = radThg(new radTExtrPolygon(inStr, mKeysOldNew, GlobalMapOfHandlers));
					}
					else if(cType3 == polyhdr.Type_g3dRelax())
					{//Instantiate Polyhedron
						hg = radThg(new radTPolyhedron(inStr, mKeysOldNew, GlobalMapOfHandlers));
					}
				}
				else if(cType2 == grp.Type_g3d())
				{//Instantiate Group
					if(cType3 == grp.Type_Group())
					{//Instantiate Group
						hg = radThg(new radTGroup(inStr, mKeysOldNew, GlobalMapOfHandlers));
					}
					else if(cType3 == sbdRecMag.Type_Group())
					{
						hg = radThg(static_cast<radTGroup*>(new radTSubdividedRecMag(inStr, mKeysOldNew, GlobalMapOfHandlers)));
					}
					else if(cType3 == sbdExtPgn.Type_Group())
					{//Instantiate Subdivided ExtrPolygon
						hg = radThg(static_cast<radTGroup*>(new radTSubdividedExtrPolygon(inStr, mKeysOldNew, GlobalMapOfHandlers)));
					}
					else if(cType3 == sbdPolyhdr.Type_Group())
					{//Instantiate Subdivided Polyhedron
						hg = radThg(static_cast<radTGroup*>(new radTSubdividedPolyhedron(inStr, mKeysOldNew, GlobalMapOfHandlers)));
					}
				}
				else if(cType2 == arcCur.Type_g3d())
				{//Instantiate ArcCur
					hg = radThg(new radTArcCur(inStr, mKeysOldNew, GlobalMapOfHandlers));
				}
				else if(cType2 == flmCur.Type_g3d())
				{//Instantiate FlmLinCur
					hg = radThg(new radTFlmLinCur(inStr, mKeysOldNew, GlobalMapOfHandlers));
				}
				else if(cType2 == bkgFldSrc.Type_g3d())
				{//Instantiate BackgroundFieldSource
					hg = radThg(new radTBackgroundFieldSource(inStr, mKeysOldNew, GlobalMapOfHandlers));
				}
			}
			else if(cType1 == tr.Type_g())
			{//Instantiate Transformation
				hg = radThg(new radTrans(inStr));
			}
			else if(cType1 == mat.Type_g())
			{
				if(cType2 == matLinIso.Type_Material())
				{//Instantiate Linear Isotropic Material
					hg = radThg(new radTLinearIsotropMaterial(inStr));
				}
				else if(cType2 == matLinAniso.Type_Material())
				{//Instantiate Linear Anisotropic Material
					hg = radThg(new radTLinearAnisotropMaterial(inStr));
				}
				else if(cType2 == matNonLinIso.Type_Material())
				{//Instantiate Non-Linear Isotropic Material
					hg = radThg(new radTNonlinearIsotropMaterial(inStr));
				}
				else if(cType2 == matNonLinAniso.Type_Material())
				{//Instantiate Non-Linear Anisotropic Material
					if(cType3 == matNonLinAniso.Type_NonlinearAnisotropMaterial())
					{
						hg = radThg(new radTNonlinearAnisotropMaterial(inStr));
					}
					else if(cType3 == matNonLinLam.Type_NonlinearAnisotropMaterial())
					{
						hg = radThg(new radTNonlinearLaminatedMaterial(inStr));
					}
				}
			}
			else if(cType1 == intrc.Type_g())
			{//Instantiate Interaction Matrix
				hg = radThg(new radTInteraction(inStr, mKeysOldNew, GlobalMapOfHandlers));
			}

			int elemKey = AddElementToContainer(hg);
			vElemKeys.push_back(elemKey);
			mKeysOldNew[oldKey] = elemKey;
		}

		//Drawing Attributes
		int nDrwAttrFound = 0;
		inStr >> nDrwAttrFound;
		//double red, green, blue, lineThick;
		for(int id=0; id<nDrwAttrFound; id++)
		{
			int oldElemKey=0;
			inStr >> oldElemKey;

			radTDrawAttr DrawAttr;
			inStr >> DrawAttr.RGB_col.Red;
			inStr >> DrawAttr.RGB_col.Green;
			inStr >> DrawAttr.RGB_col.Blue;
			inStr >> DrawAttr.LineThickness;

			int newElemKey=0;
			map<int, int>::const_iterator itOldNewKey = mKeysOldNew.find(oldElemKey);
			if(itOldNewKey != mKeysOldNew.end())
			{
				MapOfDrawAttr[itOldNewKey->second] = DrawAttr;
			}
		}

		int res = 0;
		//int trueNumElems = (int)(vElemKeys.size());
		//if(trueNumElems == 1)
		//if(nElemDir == 1)
		if((nElemDir == 1) && (!listIsExpectedInOutput))
		{
			//int elemKey = vElemKeys[0];
			map<int, int>::const_iterator itKeyOldNew = mKeysOldNew.find(*arDirElemOldKeys);
			int elemKey = 0;
			if(itKeyOldNew != mKeysOldNew.end())
			{
				elemKey = itKeyOldNew->second;
			}
			if(SendingIsRequired) Send.Int(elemKey);
			//return elemKey;
			res = elemKey;
		}
		//else if(trueNumElems > 1)
		//else if(nElemDir > 1)
		else if(listIsExpectedInOutput)
		{
			if(SendingIsRequired)
			{
				for(int j=0; j<nElemDir; j++)
				{
					int oldElemKey = arDirElemOldKeys[j];
					map<int, int>::const_iterator itKeyOldNew = mKeysOldNew.find(oldElemKey);
					int elemKey = 0;
					if(itKeyOldNew != mKeysOldNew.end())
					{
						elemKey = itKeyOldNew->second;
					}
					arDirElemOldKeys[j] = elemKey;
				}
				Send.IntList(arDirElemOldKeys, nElemDir);
			}
			//return 1;
			res = 1;
		}

		// RAII: vArDirElemOldKeys cleaned up automatically
		return res;
	}
	catch(...)
	{
		Send.ErrorMessage("Radia::Error202");
		Initialize();
		// RAII: vArDirElemOldKeys cleaned up automatically
		return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::RetrieveElemKey(const radTg* IngPtr)
{
	try
	{
		int ElemKey = 0;
		for(radTmhg::iterator GenIter = GlobalMapOfHandlers.begin();
			GenIter != GlobalMapOfHandlers.end(); ++GenIter)
			if(((*GenIter).second).rep == IngPtr) { ElemKey = (*GenIter).first; break;}

			return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

bool radTApplication::UnsafeGetElemByKey(int ElemKey, radThg& outHandle)
{
	radTmhg::iterator iter = GlobalMapOfHandlers.find(ElemKey);
	if (iter == GlobalMapOfHandlers.end()) {
		return false;
	}
	outHandle = iter->second;
	return true;
}

//-------------------------------------------------------------------------

void radTApplication::GenDump()
{
	try
	{
//#ifdef __GNUC__
//		ostrstream OutDumpStream;
//#else
		ostringstream OutDumpStream; // Porting
//#endif

		OutDumpStream << "rad: Currently in memory:\n";
		int AmOfElem = (int)(GlobalMapOfHandlers.size());
		if(AmOfElem > 0)
		{
			for(radTmhg::const_iterator iter = GlobalMapOfHandlers.begin();
				iter != GlobalMapOfHandlers.end(); ++iter)
			{
				OutDumpStream << "  Elem. No.:" << (*iter).first << endl;
				(((*iter).second).rep)->Dump(OutDumpStream, 1);
			}
			OutDumpStream << ends;

//#ifdef __GNUC__
//			Send.String(OutDumpStream.str());
//#else
			Send.String(OutDumpStream.str().c_str()); // Porting
//#endif
		}
		else Send.ErrorMessage("Radia::Error100");
	}
	catch(...)
	{
		Initialize(); return;
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

// Graphics3D functions removed: ApplyDrawAttrToElem_g3d, RemoveDrawAttrFromElem_g3d,
// GraphicsForElem_g3d, GraphicsForElem_g3d_VTK, GraphicsForAll_g3d

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

int radTApplication::PreRelax(int ElemKey, int SrcElemKey, char skipDenseMatrix)
{
	radThg hg;
	if(!ValidateElemKey(ElemKey, hg)) return 0;
	radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
	if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

	radThg hgMoreExtSrc;
	if(SrcElemKey!=0)
	{
		if(!ValidateElemKey(SrcElemKey, hgMoreExtSrc)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hgMoreExtSrc.rep);
		if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}
	}

	radTInteraction* InteractionPtr = nullptr;
	try
	{
		char AllocateExtraArray = 1; //OC300504
		char KeepTransData = 1; //OC240408 to enable update after scaling of currents

		InteractionPtr = new radTInteraction(hg, hgMoreExtSrc, CompCriterium, MemAllocForIntrctMatrTotAtOnce, AllocateExtraArray, KeepTransData, m_rankMPI, m_nProcMPI, skipDenseMatrix); //OC08012020 + skipDenseMatrix
		//radTInteraction* InteractionPtr = new radTInteraction(hg, hgMoreExtSrc, CompCriterium, MemAllocForIntrctMatrTotAtOnce, AllocateExtraArray, KeepTransData);

		if(InteractionPtr->SomethingIsWrong)
		{
			delete InteractionPtr;
			InteractionPtr = nullptr;
			return 0;
		} // The message has already been sent
		else if(!(InteractionPtr->NotEmpty())) { delete InteractionPtr; InteractionPtr = nullptr; Send.ErrorMessage("Radia::Error102"); return 0;}
		else
		{
			radThg InteractHandle(InteractionPtr);
			InteractionPtr = nullptr;  // Ownership transferred to radThg
			int InteractElemKey = AddElementToContainer(InteractHandle);
			if(SendingIsRequired) Send.Int(InteractElemKey);
			return InteractElemKey;
		}
	}
	catch (...)
	{
		if(InteractionPtr) delete InteractionPtr;
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetRelaxSubInterval(int InteractElemKey, int StartNo, int FinNo, int RelaxTogether)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(InteractElemKey, hg)) return 0;
		radTInteraction* InteractPtr = Cast.InteractCast(hg.rep);
		if(InteractPtr==0) { Send.ErrorMessage("Radia::Error017"); return 0;}

		TRelaxSubIntervalID SubIntervalID = (RelaxTogether != 0) ?
			TRelaxSubIntervalID::RelaxTogether : TRelaxSubIntervalID::RelaxApart;

		InteractPtr->AddRelaxSubInterval(StartNo, FinNo, SubIntervalID);

		return 1;
	}
	catch (...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ShowInteractMatrix(int InteractElemKey)
{
	radThg hg;
	if(!ValidateElemKey(InteractElemKey, hg)) return;
	radTInteraction* InteractPtr = Cast.InteractCast(hg.rep);
	if(InteractPtr==0) { Send.ErrorMessage("Radia::Error017"); return;}

	InteractPtr->ShowInteractMatrix();
}

//-------------------------------------------------------------------------

int radTApplication::GetInteractMatrix(int InteractElemKey, double* pMatrix, int* pDOF)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(InteractElemKey, hg)) return 0;
		radTInteraction* InteractPtr = Cast.InteractCast(hg.rep);
		if(InteractPtr==0) { Send.ErrorMessage("Radia::Error017"); return 0;}

		// Get total DOF and matrix data
		int totalDOF = InteractPtr->GetTotalDOF();
		*pDOF = totalDOF;

		if(pMatrix != nullptr && totalDOF > 0)
		{
			const double* matrixData = InteractPtr->GetFlatInteractMatrix();
			if(matrixData != nullptr)
			{
				// Copy matrix data (row-major format: A[target][source])
				long matrixSize = (long)totalDOF * (long)totalDOF;
				std::memcpy(pMatrix, matrixData, matrixSize * sizeof(double));
			}
			else
			{
				// Matrix not built - return zeros
				long matrixSize = (long)totalDOF * (long)totalDOF;
				std::memset(pMatrix, 0, matrixSize * sizeof(double));
			}
		}

		return 1;
	}
	catch (...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

void radTApplication::ShowInteractVector(int InteractElemKey, char* FieldVectID)
{
	radThg hg;
	if(!ValidateElemKey(InteractElemKey, hg)) return;
	radTInteraction* InteractPtr = Cast.InteractCast(hg.rep); 
	if(InteractPtr==0) { Send.ErrorMessage("Radia::Error017"); return;}

	if(!strcmp(FieldVectID, "ext")) InteractPtr->ShowInteractVector('E');
	else if(!strcmp(FieldVectID, "tot")) InteractPtr->ShowInteractVector('T');
	else if(!strcmp(FieldVectID, "mag")) InteractPtr->ShowInteractVector('M');
	else { Send.ErrorMessage("Radia::Error020"); return;}
}

//-------------------------------------------------------------------------

int radTApplication::MakeManualRelax(int InteractElemKey, int MethNo, int IterNumber, double RelaxParam)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(InteractElemKey, hg)) return 0;
		radTInteraction* InteractPtr = Cast.InteractCast(hg.rep); 
		if(InteractPtr==0) { Send.ErrorMessage("Radia::Error017"); return 0;}

		// Valid methods: LU (0), BICGSTAB (1), HACAPK (2)
		{
			bool validMethod = (MethNo == RadSolverMethod::LU || MethNo == RadSolverMethod::BICGSTAB);
#ifdef RADIA_USE_HACAPK
			validMethod = validMethod || (MethNo == RadSolverMethod::BICGSTAB_HMATRIX);
#endif
			if(!validMethod) { Send.ErrorMessage("Radia::Error028"); return 0;}
		}
		if(IterNumber<0) { Send.ErrorMessage("Radia::Error019"); return 0;}
		if((RelaxParam<0.) || (RelaxParam>1.)) { Send.ErrorMessage("Radia::Error018"); return 0;}

		switch(MethNo)
		{
				case RadSolverMethod::LU:
			{
				// LU direct solver - solve linear system directly
				radTRelaxationMethNo_0 RelaxMethNo_0(InteractPtr);
				// For direct solver, precision is not used; just solve once
				RelaxMethNo_0.AutoRelax(1.0e-6, 1);  // Dummy precision, single "iteration"
			}
			break;
		case RadSolverMethod::BICGSTAB:
			{
				// BiCGSTAB iterative solver (default)
				radTRelaxationMethNo_1 RelaxMethNo_1(InteractPtr);
				RelaxMethNo_1.AutoRelax(RelaxParam, IterNumber);
			}
			break;
#ifdef RADIA_USE_HACAPK
		case RadSolverMethod::BICGSTAB_HMATRIX:
			{
				// BiCGSTAB with H-matrix (HACApK ACA+)
				radTRelaxationMethNo_2 RelaxMethNo_2(InteractPtr);

				// Set HACApK parameters from application settings
				RadHACApKParams hacapk_params;
				hacapk_params.aca_eps = m_hacapk_eps;
				hacapk_params.leaf_size = m_hacapk_leaf_size;
				hacapk_params.eta = m_hacapk_eta;
				hacapk_params.print_level = 0;  // Silent (set to 1 for debug output)
				RelaxMethNo_2.SetHACApKParams(hacapk_params);

				RelaxMethNo_2.AutoRelax(RelaxParam, IterNumber);

				// Store statistics for later retrieval
				const RadHACApKStats& stats = RelaxMethNo_2.GetHMatrixStats();
				m_hacapk_n_lowrank = stats.n_lowrank;
				m_hacapk_n_dense = stats.n_dense;
				m_hacapk_max_rank = stats.max_rank;
				m_hacapk_n_leaves = stats.n_leaves;
				m_hacapk_n_dof = stats.n_dof;
				m_hacapk_compression = stats.compression;
				m_hacapk_build_time = stats.build_time;
				m_hacapk_memory_mb = stats.memory_mb;
				m_hacapk_dense_memory_mb = stats.dense_memory_mb;
				m_hacapk_stats_valid = true;
			}
			break;
#endif
		}

		int lenRelaxStatusParamArray = 3;
		double RelaxStatusParamArray[3];
		InteractPtr->OutRelaxStatusParam(RelaxStatusParamArray);

		if(SendingIsRequired) Send.DoubleList(RelaxStatusParamArray, lenRelaxStatusParamArray);
		return InteractElemKey;
	}
	catch (...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::MakeAutoRelax(int InteractElemKey, double PrecOnMagnetiz, int MaxIterNumber, int MethNo, const char** arOptionNames, const char** arOptionValues, int numOptions)
{
	// Initialize solve statistics
	m_solve_stats_valid = false;
	m_solve_t_matrix_build = 0.0;
	m_solve_t_linear_solve = 0.0;
	m_solve_t_lu_decomp = 0.0;  // Reset LU decomposition time
	m_solve_linear_iterations = 0;
	m_solve_nonl_iterations = 0;
	m_solve_num_threads = radia::GetNumThreads();  // Record while TaskManager is active

	try
	{
		radThg hg;
		if(!ValidateElemKey(InteractElemKey, hg)) return 0;
		radTInteraction* InteractPtr = Cast.InteractCast(hg.rep);
		if(InteractPtr==0) { Send.ErrorMessage("Radia::Error017"); return 0;}

		int ActualIterNum = 0; //OC02012020
		int lenRelaxStatusParamArray = 3;
		double RelaxStatusParamArray[3];

		if(m_rankMPI <= 0) //OC02012020
		{
			if(PrecOnMagnetiz <= 0.) { Send.ErrorMessage("Radia::Error030"); return 0; }
			if(MaxIterNumber <= 0) { Send.ErrorMessage("Radia::Error031"); return 0; }

			// Valid methods: LU (0), BICGSTAB (1), HACAPK (2)
			{
				bool validMethod = (MethNo == RadSolverMethod::LU || MethNo == RadSolverMethod::BICGSTAB);
#ifdef RADIA_USE_HACAPK
				validMethod = validMethod || (MethNo == RadSolverMethod::BICGSTAB_HMATRIX);
#endif
				if(!validMethod) { Send.ErrorMessage("Radia::Error041"); return 0; }
			}

			radTOptionNames OptNam;
			const char** BufNameString = arOptionNames;
			const char** BufValString = arOptionValues;
			char MagnResetIsNotNeeded = 0;
			for(int i=0; i<numOptions; i++)
			{
				if(!strcmp(*BufNameString, OptNam.ZeroM))
				{
					if(!strcmp(*BufValString, (OptNam.ZeroM_Values)[0])) MagnResetIsNotNeeded = 1; //no
					else if(!strcmp(*BufValString, (OptNam.ZeroM_Values)[1])) MagnResetIsNotNeeded = 0; //yes
					else if(!strcmp(*BufValString, (OptNam.ZeroM_Values)[2])) MagnResetIsNotNeeded = 1; //false
					else if(!strcmp(*BufValString, (OptNam.ZeroM_Values)[3])) MagnResetIsNotNeeded = 0; //true
					else { Send.ErrorMessage("Radia::Error062"); return 0; }
				}
				else { Send.ErrorMessage("Radia::Error062"); return 0; }
				BufNameString++; BufValString++;
			}

			//int ActualIterNum = 0;

			// B-input solvers for hysteresis (energy or play)
			// Auto-activates when b_input_newton=True or b_input_hantila=True
			// and all materials are radTHysteresisMaterial
			if(m_b_input_newton || m_b_input_hantila)
			{
				// Check if all elements have hysteresis materials
				bool all_hysteresis = true;
				int n_elem = InteractPtr->GetAmOfMainElem();
				for(int i = 0; i < n_elem && all_hysteresis; i++)
				{
					radTg3dRelax* g3d = InteractPtr->GetElement(i);
					radTMaterial* mat = (radTMaterial*)(g3d->MaterHandle.rep);
					if(dynamic_cast<radTHysteresisMaterial*>(mat) == nullptr)
						all_hysteresis = false;
				}

				if(all_hysteresis && n_elem > 0)
				{
					radTRelaxationMethNo_0 RelaxMethNo_0(InteractPtr);

					if(m_b_input_hantila)
					{
						// Hybrid: Newton warmup (3 iter) -> Hantila refinement
						// Newton gets close quickly (quadratic convergence),
						// then Hantila refines with constant LHS (O(N^2) per iter)
						int newton_warmup = 3;
						int newton_iters = RelaxMethNo_0.AutoRelax_BInput_Newton(
							PrecOnMagnetiz, newton_warmup, MagnResetIsNotNeeded);

						// If Newton already converged, skip Hantila
						double rsp[3]; InteractPtr->OutRelaxStatusParam(rsp);
						double newton_residual = rsp[0];  // MisfitM
						if(newton_iters < newton_warmup || newton_residual >= PrecOnMagnetiz)
						{
							// Newton didn't fully converge -> refine with Hantila
							ActualIterNum = newton_iters + RelaxMethNo_0.AutoRelax_BInput_Hantila(
								PrecOnMagnetiz, MaxIterNumber,
								m_hantila_alpha, m_hantila_relax,
								/*MagnResetIsNotNeeded=*/1);
						}
						else
						{
							// Newton converged within warmup -> use Newton result
							ActualIterNum = newton_iters;
						}
					}
					else
					{
						// Use B-input Newton (full Jacobian each iteration)
						ActualIterNum = RelaxMethNo_0.AutoRelax_BInput_Newton(
							PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
					}

					// Skip standard solver
					goto binput_done;
				}
			}

			switch(MethNo)
			{
						case RadSolverMethod::LU:
			{
				radTRelaxationMethNo_0 RelaxMethNo_0(InteractPtr);
				ActualIterNum = RelaxMethNo_0.AutoRelax(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
			}
			break;
			case RadSolverMethod::BICGSTAB:
			{
				// BiCGSTAB iterative solver
				radTRelaxationMethNo_1 RelaxMethNo_1(InteractPtr);
				ActualIterNum = RelaxMethNo_1.AutoRelax(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
			}
			break;
#ifdef RADIA_USE_HACAPK
			case RadSolverMethod::BICGSTAB_HMATRIX:
			{
				// BiCGSTAB with H-matrix (HACApK ACA+)
				radTRelaxationMethNo_2 RelaxMethNo_2(InteractPtr);

				// Set HACApK parameters from application settings
				RadHACApKParams hacapk_params;
				hacapk_params.aca_eps = m_hacapk_eps;
				hacapk_params.leaf_size = m_hacapk_leaf_size;
				hacapk_params.eta = m_hacapk_eta;
				hacapk_params.print_level = 0;  // Silent (set to 1 for debug output)
				RelaxMethNo_2.SetHACApKParams(hacapk_params);

				ActualIterNum = RelaxMethNo_2.AutoRelax(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);

				// Store statistics for later retrieval
				const RadHACApKStats& stats = RelaxMethNo_2.GetHMatrixStats();
				m_hacapk_n_lowrank = stats.n_lowrank;
				m_hacapk_n_dense = stats.n_dense;
				m_hacapk_max_rank = stats.max_rank;
				m_hacapk_n_leaves = stats.n_leaves;
				m_hacapk_n_dof = stats.n_dof;
				m_hacapk_compression = stats.compression;
				m_hacapk_build_time = stats.build_time;
				m_hacapk_memory_mb = stats.memory_mb;
				m_hacapk_dense_memory_mb = stats.dense_memory_mb;
				m_hacapk_stats_valid = true;
			}
			break;
#endif
			}

			binput_done:
			InteractPtr->OutRelaxStatusParam(RelaxStatusParamArray);
		}

		// Set solve statistics
		m_solve_nonl_iterations = ActualIterNum;
		m_solve_stats_valid = true;

		if(ActualIterNum >= MaxIterNumber) { Send.WarningMessage("Radia::Warning015");}
		if(SendingIsRequired) Send.OutRelaxResultsInfo(RelaxStatusParamArray, lenRelaxStatusParamArray, ActualIterNum);

		return ActualIterNum;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::UpdateSourcesForRelax(int InteractElemKey)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(InteractElemKey, hg)) return 0;
		radTInteraction* InteractPtr = Cast.InteractCast(hg.rep); 
		if(InteractPtr==0) { Send.ErrorMessage("Radia::Error017"); return 0;}

		InteractPtr->UpdateExternalField();

		if(SendingIsRequired) Send.Int(InteractElemKey);
		return InteractElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::SolveGen(int ObjKey, double PrecOnMagnetiz, int MaxIterNumber, int MethNo, const char* image)
{
	// Methods 6-7 have been removed (deprecated)
	// All methods now go through the standard PreRelax + MakeAutoRelax path
	long ActualIterNum = 0;
	try
	{
		short PrevSendingIsRequired = SendingIsRequired;
		SendingIsRequired = 0;

		// For HACApK (method 2), skip dense matrix construction
		// HACApK builds its own data structures
		char skipDenseMatrix = 0;
#ifdef RADIA_USE_HACAPK
		if(MethNo == RadSolverMethod::BICGSTAB_HMATRIX)
		{
			skipDenseMatrix = 1;
		}
#endif

		// Check if we can reuse cached interaction matrix
		// The interaction matrix N only depends on geometry, not on chi (material)
		// This avoids the expensive O(N^2) matrix construction on repeated Solve() calls
		int InteractElemKey = 0;
		double t_matrix_build = 0.0;

		// Normalize image string for comparison
		std::string imageSpec = (image != nullptr) ? image : "";

		bool cacheValid = (m_cached_interact_key > 0 && m_cached_obj_key == ObjKey);

		// Invalidate cache if image spec changed
		if(cacheValid && imageSpec != m_cached_image_spec)
		{
			cacheValid = false;
			m_cached_interact_key = 0;
			m_cached_obj_key = 0;
			m_cached_image_spec.clear();
		}

		// For HACApK, cache invalidation is handled differently (H-matrix is rebuilt internally)
		// Also invalidate cache if solver method changed (different matrix type might be needed)
#ifdef RADIA_USE_HACAPK
		if(cacheValid && MethNo == RadSolverMethod::BICGSTAB_HMATRIX)
		{
			// HACApK builds its own H-matrix, may need to rebuild
			// For now, keep the cache - HACApK handles its own caching
		}
#endif

		if(cacheValid)
		{
			// Reuse cached interaction object
			InteractElemKey = m_cached_interact_key;
			t_matrix_build = 0.0;  // No rebuild needed

			// Validate the cached key is still valid in GlobalMapOfHandlers
			radThg hg;
			if(!ValidateElemKey(InteractElemKey, hg))
			{
				// Cache is stale, need to rebuild
				cacheValid = false;
				m_cached_interact_key = 0;
				m_cached_obj_key = 0;
				m_cached_image_spec.clear();
			}
			else
			{
				// Re-evaluate external field (ObjBckg callback may have changed)
				// The interaction matrix N is geometry-only and can be cached,
				// but external fields must be refreshed for hysteresis stepping
				radTInteraction* pIntrc = dynamic_cast<radTInteraction*>(hg.rep);
				if(pIntrc != nullptr)
					pIntrc->SetupExternFieldArray();
			}
		}

		if(!cacheValid)
		{
			// Time matrix construction (PreRelax includes SetupInteractMatrix)
			auto t_prerelax_start = std::chrono::high_resolution_clock::now();
			InteractElemKey = PreRelax(ObjKey, 0, skipDenseMatrix);
			auto t_prerelax_end = std::chrono::high_resolution_clock::now();
			t_matrix_build = std::chrono::duration<double>(t_prerelax_end - t_prerelax_start).count();
			if(InteractElemKey <= 0) return 0;

			// Apply IMA symmetry if specified
			if(!imageSpec.empty())
			{
				// Get interaction object
				radThg hg;
				if(ValidateElemKey(InteractElemKey, hg))
				{
					radTInteraction* pIntrc = dynamic_cast<radTInteraction*>(hg.rep);
					if(pIntrc != nullptr)
					{
						// Parse and apply IMA symmetry
						// For HACApK, skip dense IMA matrix (kernel handles IMA)
						ApplyIMASymmetryToInteraction(pIntrc, imageSpec.c_str(), skipDenseMatrix != 0);
					}
				}
			}

			// Cache the interaction key for future reuse
			m_cached_interact_key = InteractElemKey;
			m_cached_obj_key = ObjKey;
			m_cached_image_spec = imageSpec;
		}

		SendingIsRequired = PrevSendingIsRequired;

		try
		{
			ActualIterNum = MakeAutoRelax(InteractElemKey, PrecOnMagnetiz, MaxIterNumber, MethNo);
			// Store matrix build time (MakeAutoRelax resets m_solve_t_matrix_build to 0)
			m_solve_t_matrix_build = t_matrix_build;
		}
		catch(...)
		{
			SendingIsRequired = 0;
			// Don't delete cached interaction on error - keep for potential retry
			throw 0;
		}

		// DON'T delete the interaction element - keep it cached for next Solve() call
		// The cache will be invalidated when geometry changes or UtiDelAll() is called
		// PrevSendingIsRequired = SendingIsRequired; SendingIsRequired = 0;
		// DeleteElement(InteractElemKey);
		// SendingIsRequired = PrevSendingIsRequired;
	}
	catch(...) { Initialize(); return 0;}
	return ActualIterNum;
}

//-------------------------------------------------------------------------

int radTApplication::SolveGenNonl(int ObjKey, double PrecOnMagnetiz, int MaxIterNumber, int MethNo, int NonlMethod, const char* image)
{
	// Set nonlinear method before calling SolveGen
	// 0 = mucal1 (chi-change), 1 = mucal2 (B-change/Newton)
	int OldNonlMethod = this->NonlinearMethod;
	this->NonlinearMethod = NonlMethod;

	int result = SolveGen(ObjKey, PrecOnMagnetiz, MaxIterNumber, MethNo, image);

	// Restore previous setting
	this->NonlinearMethod = OldNonlMethod;

	return result;
}

//-------------------------------------------------------------------------

int radTApplication::BuildMatrix(int ObjKey, const char* image)
{
	// Build interaction matrix without solving
	// This allows users to inspect the matrix before solving
	int InteractElemKey = 0;
	try
	{
		short PrevSendingIsRequired = SendingIsRequired;
		SendingIsRequired = 0;

		// Normalize image string for comparison
		std::string imageSpec = (image != nullptr) ? image : "";

		// Invalidate cache if object or image spec changed
		bool cacheValid = (m_cached_interact_key > 0 && m_cached_obj_key == ObjKey);
		if(cacheValid && imageSpec != m_cached_image_spec)
		{
			cacheValid = false;
			m_cached_interact_key = 0;
			m_cached_obj_key = 0;
			m_cached_image_spec.clear();
		}

		if(cacheValid)
		{
			// Reuse cached interaction object
			InteractElemKey = m_cached_interact_key;

			// Validate the cached key is still valid in GlobalMapOfHandlers
			radThg hg;
			if(!ValidateElemKey(InteractElemKey, hg))
			{
				// Cache is stale, need to rebuild
				cacheValid = false;
				m_cached_interact_key = 0;
				m_cached_obj_key = 0;
				m_cached_image_spec.clear();
			}
		}

		if(!cacheValid)
		{
			// Build interaction matrix (PreRelax includes SetupInteractMatrix)
			InteractElemKey = PreRelax(ObjKey, 0, 0);
			if(InteractElemKey <= 0)
			{
				SendingIsRequired = PrevSendingIsRequired;
				return 0;
			}

			// Apply IMA symmetry if specified
			if(!imageSpec.empty())
			{
				// Get interaction object
				radThg hg;
				if(ValidateElemKey(InteractElemKey, hg))
				{
					radTInteraction* pIntrc = dynamic_cast<radTInteraction*>(hg.rep);
					if(pIntrc != nullptr)
					{
						// Parse and apply IMA symmetry
						ApplyIMASymmetryToInteraction(pIntrc, imageSpec.c_str());
					}
				}
			}

			// Cache the interaction key for future reuse
			m_cached_interact_key = InteractElemKey;
			m_cached_obj_key = ObjKey;
			m_cached_image_spec = imageSpec;
		}

		SendingIsRequired = PrevSendingIsRequired;
	}
	catch(...)
	{
		Initialize();
		return 0;
	}
	return InteractElemKey;
}

//-------------------------------------------------------------------------
// ApplyIMASymmetryToInteraction: Parse image string and apply IMA symmetry
// Format: "+x", "-z", "+x-z", "+x+y-z", etc.
// Each axis is prefixed by + (symmetric) or - (antisymmetric)
// No prefix defaults to + (symmetric)
// Examples:
//   "+x" -> X mirror with symmetric BC
//   "-z" -> Z mirror with antisymmetric BC
//   "+x-z" -> X and Z mirrors, X symmetric, Z antisymmetric (ELF quarter model)
//-------------------------------------------------------------------------
bool radTApplication::ApplyIMASymmetryToInteraction(radTInteraction* pIntrc, const char* imageSpec, bool skipDenseIMA)
{
	if(pIntrc == nullptr || imageSpec == nullptr || imageSpec[0] == '\0')
	{
		return false;
	}

	std::string spec(imageSpec);

	// Parse the image specification string
	int symmetryFlags = radTInteraction::IMA_NONE;
	int signX = 1, signY = 1, signZ = 1;  // Default to symmetric (+)

	size_t pos = 0;
	while(pos < spec.length())
	{
		// Skip whitespace
		while(pos < spec.length() && (spec[pos] == ' ' || spec[pos] == '\t'))
			pos++;

		if(pos >= spec.length()) break;

		// Check for sign prefix
		int currentSign = 1;  // Default to symmetric
		if(spec[pos] == '+')
		{
			currentSign = 1;
			pos++;
		}
		else if(spec[pos] == '-')
		{
			currentSign = -1;
			pos++;
		}

		if(pos >= spec.length()) break;

		// Read axis
		char axis = tolower(spec[pos]);
		pos++;

		switch(axis)
		{
		case 'x':
			symmetryFlags |= radTInteraction::IMA_X;
			signX = currentSign;
			break;
		case 'y':
			symmetryFlags |= radTInteraction::IMA_Y;
			signY = currentSign;
			break;
		case 'z':
			symmetryFlags |= radTInteraction::IMA_Z;
			signZ = currentSign;
			break;
		default:
			// Invalid axis - skip
			break;
		}
	}

	if(symmetryFlags == radTInteraction::IMA_NONE)
	{
		return false;  // No valid symmetry specified
	}

	// Apply IMA symmetry with per-axis signs
	int numElements = pIntrc->SetIMASymmetry(symmetryFlags, signX, signY, signZ);

	// Setup IMA subsystem (element reduction + optional dense matrix)
	// For HACApK: skipDenseIMA=true skips dense matrix, kernel handles IMA
	if(pIntrc->IsIMAEnabled())
	{
		pIntrc->SetupInteractMatrix_IMA(skipDenseIMA);
	}

	return (numElements > 0);
}

//-------------------------------------------------------------------------
/**
void radTApplication::OutFieldCompRes(char* FieldChar, radTField* FieldArray, double* ArgArray, int Np)
{
	char* BufChar = FieldChar;
	char* EqEmptyStr = "BHAM";

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			if((*BufChar == 'B') || (*BufChar == 'b') || 
			   (*BufChar == 'H') || (*BufChar == 'h') ||
			   (*BufChar == 'A') || (*BufChar == 'a') ||
			   (*BufChar == 'M') || (*BufChar == 'm') ||
			   (*BufChar == 'P') || (*BufChar == 'p')) ItemCount++;
			BufChar++;
		}
		BufChar = FieldChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 4;
	}
	char* ActualInitCharPtr = BufChar;

	if(Np > 1) Send.InitOutList(Np);

	radTField* FieldPtr = FieldArray;
	for(int i=0; i<Np; i++)
	{
		if(ArgArray != nullptr) // Argument Needed
		{
			Send.InitOutList(2);
			Send.Double(ArgArray[i]);
		}

		if(ItemCount > 1) Send.InitOutList(ItemCount);
		while (*BufChar != '\0') 
		{
			char* BufChar_p_1 = BufChar+1;
			if(*(BufChar)=='B' || *(BufChar)=='b')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->B.x);
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->B.y);
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->B.z);
				else Send.Vector3d(&(FieldPtr->B));
			}
			else if(*(BufChar)=='H' || *(BufChar)=='h')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->H.x);
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->H.y);
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->H.z);
				else Send.Vector3d(&(FieldPtr->H));
			}
			else if(*(BufChar)=='A' || *(BufChar)=='a')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->A.x);
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->A.y);
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->A.z);
				else Send.Vector3d(&(FieldPtr->A));
			}
			else if(*(BufChar)=='M' || *(BufChar)=='m')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->M.x);
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->M.y);
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->M.z);
				else Send.Vector3d(&(FieldPtr->M));
			}
			else if(*(BufChar)=='P' || *(BufChar)=='p')	Send.Double(FieldPtr->Phi);
			BufChar++;
		}
		FieldPtr++;
		BufChar = ActualInitCharPtr;
	}
}
**/
//-------------------------------------------------------------------------

void radTApplication::OutFieldCompRes(char* FieldChar, radTField* FieldArray, long Np)
{
	char* BufChar = FieldChar;
	//char* EqEmptyStr = "BHAMJ";
	char EqEmptyStr[] = "BHAMJ"; //OC01052013

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			if((*BufChar == 'B') || (*BufChar == 'b') || 
			   (*BufChar == 'H') || (*BufChar == 'h') ||
			   (*BufChar == 'A') || (*BufChar == 'a') ||
			   (*BufChar == 'M') || (*BufChar == 'm') ||
			   (*BufChar == 'J') || (*BufChar == 'j') ||
			   (*BufChar == 'P') || (*BufChar == 'p') ||
			   (*BufChar == 'Q') || (*BufChar == 'q')) ItemCount++;
			BufChar++;
		}
		BufChar = FieldChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 4;
	}
	char* ActualInitCharPtr = BufChar;

	// RAII: Use std::vector for automatic cleanup
	std::vector<double> vTotOutArray(Np*12);
	double *TotOutArray = vTotOutArray.data();

	int InnerCount=0;
	// Magnetic constant (permeability of free space)
	const double Mu0 = 4. * 3.1415926535897932 * 1.e-7; // T*m/A
	
	radTField* FieldPtr = FieldArray;
	double *t = TotOutArray;
	for(int i=0; i<Np; i++)
	{
		InnerCount = 0;
		while (*BufChar != '\0') 
		{
			char* BufChar_p_1 = BufChar+1;
			if(*(BufChar)=='B' || *(BufChar)=='b')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->B.x * Mu0; InnerCount++; } // Fix: Convert A/m to Tesla
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->B.y * Mu0; InnerCount++; } // Fix: Convert A/m to Tesla
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->B.z * Mu0; InnerCount++; } // Fix: Convert A/m to Tesla
				else { *(t++) = FieldPtr->B.x * Mu0; *(t++) = FieldPtr->B.y * Mu0; *(t++) = FieldPtr->B.z * Mu0; InnerCount += 3;} // Fix: Convert A/m to Tesla
			}
			else if(*(BufChar)=='H' || *(BufChar)=='h')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->H.x; InnerCount++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->H.y; InnerCount++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->H.z; InnerCount++;}
				else { *(t++) = FieldPtr->H.x; *(t++) = FieldPtr->H.y; *(t++) = FieldPtr->H.z; InnerCount += 3;}
			}
			else if(*(BufChar)=='A' || *(BufChar)=='a')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->A.x; InnerCount++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->A.y; InnerCount++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->A.z; InnerCount++;}
				else { *(t++) = FieldPtr->A.x; *(t++) = FieldPtr->A.y; *(t++) = FieldPtr->A.z; InnerCount += 3;}
			}
			else if(*(BufChar)=='M' || *(BufChar)=='m')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->M.x; InnerCount++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->M.y; InnerCount++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->M.z; InnerCount++;}
				else { *(t++) = FieldPtr->M.x; *(t++) = FieldPtr->M.y; *(t++) = FieldPtr->M.z; InnerCount += 3;}
			}
			else if(*(BufChar)=='J' || *(BufChar)=='j')
			{
				if(*BufChar_p_1=='x' || *BufChar_p_1=='X') { *(t++) = FieldPtr->J.x; InnerCount++;}
				else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') { *(t++) = FieldPtr->J.y; InnerCount++;}
				else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') { *(t++) = FieldPtr->J.z; InnerCount++;}
				else { *(t++) = FieldPtr->J.x; *(t++) = FieldPtr->J.y; *(t++) = FieldPtr->J.z; InnerCount += 3;}
			}
			else if(*(BufChar)=='P' || *(BufChar)=='p')
			{
				// Scalar potential phi_m (units: Ampere)
				*(t++) = FieldPtr->Phi; InnerCount++;
			}
			else if(*(BufChar)=='Q' || *(BufChar)=='q') //OC161005
			{
				*(t++) = FieldPtr->B.x * Mu0; *(t++) = FieldPtr->B.y * Mu0; *(t++) = FieldPtr->B.z * Mu0; InnerCount += 3; // Fix: Convert A/m to Tesla
				*(t++) = FieldPtr->H.x; *(t++) = FieldPtr->H.y; *(t++) = FieldPtr->H.z; InnerCount += 3;
				*(t++) = FieldPtr->A.x; *(t++) = FieldPtr->A.y; *(t++) = FieldPtr->A.z; InnerCount += 3;
			}
			
			BufChar++;
		}
		FieldPtr++;
		BufChar = ActualInitCharPtr;
	}
	int Dims[] = { InnerCount, Np};
	Send.MultiDimArrayOfDouble(TotOutArray, Dims, 2);
	// RAII: automatic cleanup
}

//-------------------------------------------------------------------------

void radTApplication::OutFieldCompRes(char* FieldChar, radTField* FieldArray, long Np, radTVectInputCell& VectInputCell)
{
	char* BufChar = FieldChar;
	//char* EqEmptyStr = "BHAMJ";
	char EqEmptyStr[] = "BHAMJ"; //OC01052013

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			if((*BufChar == 'B') || (*BufChar == 'b') || 
			   (*BufChar == 'H') || (*BufChar == 'h') ||
			   (*BufChar == 'A') || (*BufChar == 'a') ||
			   (*BufChar == 'M') || (*BufChar == 'm') ||
			   (*BufChar == 'J') || (*BufChar == 'j') ||
			   (*BufChar == 'P') || (*BufChar == 'p')) ItemCount++;
			BufChar++;
		}
		BufChar = FieldChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 4;
	}
	char* ActualInitCharPtr = BufChar;

	radTField* tField = FieldArray;
	for(radTVectInputCell::iterator iterCell = VectInputCell.begin(); iterCell != VectInputCell.end(); ++iterCell)
	{
		radTInputCell& Cell = *iterCell;
		if(Cell.Type == 'L')
		{
			Send.InitOutList(Cell.AuxNum);
		}
		else if(Cell.Type == 'P')
		{
			ParseAndSendOneFieldValue(tField++, BufChar, ItemCount);
		}
	}
}

//-------------------------------------------------------------------------

void radTApplication::ParseAndSendOneFieldValue(radTField* FieldPtr, char* BufChar, int AmOfItem)
{
	if(AmOfItem > 1) Send.InitOutList(AmOfItem);
	while(*BufChar != '\0') 
	{
		char* BufChar_p_1 = BufChar+1;
		if(*(BufChar)=='B' || *(BufChar)=='b')
		{
			if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->B.x);
			else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->B.y);
			else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->B.z);
			else Send.Vector3d(&(FieldPtr->B));
		}
		else if(*(BufChar)=='H' || *(BufChar)=='h')
		{
			if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->H.x);
			else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->H.y);
			else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->H.z);
			else Send.Vector3d(&(FieldPtr->H));
		}
		else if(*(BufChar)=='A' || *(BufChar)=='a')
		{
			if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->A.x);
			else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->A.y);
			else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->A.z);
			else Send.Vector3d(&(FieldPtr->A));
		}
		else if(*(BufChar)=='M' || *(BufChar)=='m')
		{
			if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->M.x);
			else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->M.y);
			else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->M.z);
			else Send.Vector3d(&(FieldPtr->M));
		}
		else if(*(BufChar)=='J' || *(BufChar)=='j')
		{
			if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->J.x);
			else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->J.y);
			else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->J.z);
			else Send.Vector3d(&(FieldPtr->J));
		}
		else if(*(BufChar)=='P' || *(BufChar)=='p')	Send.Double(FieldPtr->Phi);
		BufChar++;
	}
}

//-------------------------------------------------------------------------
/**
void radTApplication::OutFieldIntCompRes(char* FieldIntChar, radTField* FieldPtr)
{
	char* BufChar = FieldIntChar;
	char* BufCharPrev = nullptr;
	char* EqEmptyStr = "Ib";

	short I_used = 0;
	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			if(((*BufChar == 'B') || (*BufChar == 'b') || 
			    (*BufChar == 'H') || (*BufChar == 'h')) ||
			   (((*BufChar == 'X') || (*BufChar == 'x') ||
			     (*BufChar == 'Y') || (*BufChar == 'y') ||
				 (*BufChar == 'Z') || (*BufChar == 'z')) &&
				(*BufCharPrev != 'B') && (*BufCharPrev != 'b') &&
				(*BufCharPrev != 'H') && (*BufCharPrev != 'h'))) ItemCount++;

			if((*BufChar == 'I') || (*BufChar == 'i')) I_used = 1;
			BufCharPrev = BufChar;
			BufChar++;
		}
		BufChar = FieldIntChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 1;
	}
	if(I_used && (ItemCount == 0))
	{
		BufChar = EqEmptyStr;
		ItemCount = 1;
	}

	if(ItemCount > 1) Send.InitOutList(ItemCount);

	while (*BufChar != '\0') 
	{
		char* BufChar_pl_1 = BufChar+1;
		char* BufChar_mi_1 = BufChar-1;

		if((*BufChar =='I') || (*BufChar == 'i'))
		{
			if((*BufChar_pl_1 == 'X') || (*BufChar_pl_1 == 'x')) Send.Double(FieldPtr->Ib.x);
			else if((*BufChar_pl_1 == 'Y') || (*BufChar_pl_1 == 'y')) Send.Double(FieldPtr->Ib.y);
			else if((*BufChar_pl_1 == 'Z') || (*BufChar_pl_1 == 'z')) Send.Double(FieldPtr->Ib.z);
			else if((*BufChar_pl_1 != 'B') && (*BufChar_pl_1 != 'b') &&
					(*BufChar_pl_1 != 'H') && (*BufChar_pl_1 != 'h') &&
					(*BufChar_pl_1 != 'X') && (*BufChar_pl_1 != 'x') &&
					(*BufChar_pl_1 != 'Y') && (*BufChar_pl_1 != 'y') &&
					(*BufChar_pl_1 != 'Z') && (*BufChar_pl_1 != 'z')) { Send.Vector3d(&(FieldPtr->Ib));	break;}
		}
		else if((*BufChar == 'B') || (*BufChar == 'b'))
		{
			if((*BufChar_pl_1 == 'X') || (*BufChar_pl_1 == 'x')) Send.Double(FieldPtr->Ib.x);
			else if((*BufChar_pl_1 == 'Y') || (*BufChar_pl_1 == 'y')) Send.Double(FieldPtr->Ib.y);
			else if((*BufChar_pl_1 == 'Z') || (*BufChar_pl_1 == 'z')) Send.Double(FieldPtr->Ib.z);
			else Send.Vector3d(&(FieldPtr->Ib));
		}
		else if((*BufChar == 'H') || (*BufChar == 'h'))
		{
			if((*BufChar_pl_1 == 'X') || (*BufChar_pl_1 == 'x')) Send.Double(FieldPtr->Ih.x);
			else if((*BufChar_pl_1 == 'Y') || (*BufChar_pl_1 == 'y')) Send.Double(FieldPtr->Ih.y);
			else if((*BufChar_pl_1 == 'Z') || (*BufChar_pl_1 == 'z')) Send.Double(FieldPtr->Ih.z);
			else Send.Vector3d(&(FieldPtr->Ih));
		}
		else if(((*BufChar == 'X') || (*BufChar == 'x')) &&
				(*BufChar_mi_1 != 'I') && (*BufChar_mi_1 != 'i') &&
				(*BufChar_mi_1 != 'B') && (*BufChar_mi_1 != 'b') &&
 				(*BufChar_mi_1 != 'H') && (*BufChar_mi_1 != 'h')) Send.Double(FieldPtr->Ib.x);
		else if(((*BufChar == 'Y') || (*BufChar == 'y')) &&
				(*BufChar_mi_1 != 'I') && (*BufChar_mi_1 != 'i') &&
				(*BufChar_mi_1 != 'B') && (*BufChar_mi_1 != 'b') &&
 				(*BufChar_mi_1 != 'H') && (*BufChar_mi_1 != 'h')) Send.Double(FieldPtr->Ib.y);
		else if(((*BufChar == 'Z') || (*BufChar == 'z')) &&
				(*BufChar_mi_1 != 'I') && (*BufChar_mi_1 != 'i') &&
				(*BufChar_mi_1 != 'B') && (*BufChar_mi_1 != 'b') &&
 				(*BufChar_mi_1 != 'H') && (*BufChar_mi_1 != 'h')) Send.Double(FieldPtr->Ib.z);
		BufChar++;
	}
}
**/
//-------------------------------------------------------------------------

void radTApplication::OutFieldEnergyForceCompRes(char* FieldChar, radTField* FieldPtr)
{
	char* BufChar = FieldChar;
	//char* EqEmptyStr = "EFT";
	char EqEmptyStr[] = "EFT"; //OC01052013

	int ItemCount = 0;
	if(*BufChar != '\0')
	{
		while (*BufChar != '\0') 
		{
			if((*BufChar == 'E') || (*BufChar == 'e') || 
			   (*BufChar == 'F') || (*BufChar == 'f') ||
			   (*BufChar == 'T') || (*BufChar == 't')) ItemCount++;
			BufChar++;
		}
		BufChar = FieldChar;
	}
	else
	{
		BufChar = EqEmptyStr;
		ItemCount = 3;
	}
	char* ActualInitCharPtr = BufChar;

	if(ItemCount > 1) Send.InitOutList(ItemCount);
	while (*BufChar != '\0') 
	{
		char* BufChar_p_1 = BufChar+1;

		if(*(BufChar)=='E' || *(BufChar)=='e') Send.Double(FieldPtr->Energy);
		else if(*(BufChar)=='F' || *(BufChar)=='f')
		{
			if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->Force.x);
			else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->Force.y);
			else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->Force.z);
			else Send.Vector3d(&(FieldPtr->Force));
		}
		else if(*(BufChar)=='T' || *(BufChar)=='t')
		{
			if(*BufChar_p_1=='x' || *BufChar_p_1=='X') Send.Double(FieldPtr->Torque.x);
			else if(*BufChar_p_1=='y' || *BufChar_p_1=='Y') Send.Double(FieldPtr->Torque.y);
			else if(*BufChar_p_1=='z' || *BufChar_p_1=='Z') Send.Double(FieldPtr->Torque.z);
			else Send.Vector3d(&(FieldPtr->Torque));
		}
		BufChar++;
	}
}

//-------------------------------------------------------------------------

void radTApplication::OutCenFieldCompRes(radTVectPairOfVect3d* pVectPairOfVect3d)
{
	int AmOfPoints = (int)(pVectPairOfVect3d->size());
	radTSend Send;
	Send.InitOutList(AmOfPoints);

	for(int i=0; i<AmOfPoints; i++)
	{
		Send.InitOutList(2);
		radTPairOfVect3d& Pair = (*pVectPairOfVect3d)[i];
		Send.Vector3d(&(Pair.V1));
		Send.Vector3d(&(Pair.V2));
	}
}

//-------------------------------------------------------------------------

int radTApplication::SetCompPrecisions(const char** ValNames, double* Values, int ValCount)
{
	try
	{
		radTOptionNames OptionNames;
		const char** BufString = ValNames;
		double* Ptr = Values;
		for(int i=0; i<ValCount; i++)
		{
			if(!strcmp(*BufString, OptionNames.B)) CompCriterium.AbsPrecB = *Ptr;
			else if(!strcmp(*BufString, OptionNames.A)) CompCriterium.AbsPrecA = *Ptr;
			else if(!strcmp(*BufString, OptionNames.BInt)) CompCriterium.AbsPrecB_int = *Ptr;
			else if(!strcmp(*BufString, OptionNames.Force)) CompCriterium.AbsPrecForce = *Ptr;
			else if(!strcmp(*BufString, OptionNames.Torque)) CompCriterium.AbsPrecTorque = *Ptr;
			else if(!strcmp(*BufString, OptionNames.Energy)) CompCriterium.AbsPrecEnergy = *Ptr;
			else if(!strcmp(*BufString, OptionNames.Coord)) CompCriterium.AbsPrecTrjCoord = *Ptr;
			else if(!strcmp(*BufString, OptionNames.Angle)) CompCriterium.AbsPrecTrjAngle = *Ptr;
			else { Send.ErrorMessage("Radia::Error057"); return 0;}
			BufString++; Ptr++;
		}
		if(SendingIsRequired) Send.Int(1);
		return 1;
	}
	catch(...) { Initialize(); return 0;}
}

//-------------------------------------------------------------------------

int radTApplication::SetCompCriterium(double InAbsPrecB, double InAbsPrecA, double InAbsPrecB_int, double InAbsPrecFrc, double InAbsPrecTrjCoord, double InAbsPrecTrjAngle)
{
	short InBasedOnPrecFlag = 0;

	try
	{
		CompCriterium.BasedOnPrecLevel = InBasedOnPrecFlag; 
		CompCriterium.AbsPrecB = InAbsPrecB;
		CompCriterium.AbsPrecA = InAbsPrecA;
		CompCriterium.AbsPrecB_int = InAbsPrecB_int;
		CompCriterium.AbsPrecForce = InAbsPrecFrc;
		CompCriterium.AbsPrecTrjCoord = InAbsPrecTrjCoord;
		CompCriterium.AbsPrecTrjAngle = InAbsPrecTrjAngle;

		if(SendingIsRequired) Send.Int(1);
		return 1;
	}
	catch(...) { Initialize(); return 0;}
}

//-------------------------------------------------------------------------

int radTApplication::SetMltplThresh(double* InMltplThresh) // Maybe to be removed later
{
	for(int i=0; i<4; i++) CompCriterium.MltplThresh[i] = InMltplThresh[i]*InMltplThresh[i];
	if(SendingIsRequired) Send.Int(1);
	return 1;
}

//-------------------------------------------------------------------------

int radTApplication::SetTolForConvergence(double AbsRandMagnitude, double RelRandMagnitude, double ZeroRandMagnitude)
{
	CnRep.SwitchActOnDoubles(1, fabs(AbsRandMagnitude), fabs(RelRandMagnitude), fabs(ZeroRandMagnitude));
	if(SendingIsRequired) Send.Int(1);
	return 1;
}

//-------------------------------------------------------------------------

int radTApplication::RandomizationOnOrOff(char* OnOrOff)
{
	char SwitchOn;
	if((!strcmp(OnOrOff, "on")) || (!strcmp(OnOrOff, "On")) || (!strcmp(OnOrOff, "ON"))) SwitchOn = 1;
	else if((!strcmp(OnOrOff, "off")) || (!strcmp(OnOrOff, "Off")) || (!strcmp(OnOrOff, "OFF"))) SwitchOn = 0;
	else { Send.ErrorMessage("Radia::Error043"); return 0;}

	if(SwitchOn) CnRep = CnRepAux;
	else
	{
		CnRepAux = CnRep;
		CnRep.AbsRand = CnRep.RelRand = CnRep.ZeroRand = 0.;
	}
	if(SendingIsRequired) Send.Int(int(SwitchOn));
	return 1;
}

//-------------------------------------------------------------------------

#ifdef RADIA_USE_HACAPK
void radTApplication::SetHACApKParams(double eps, int leaf_size, double eta)
{
	// Update only if positive value provided (-1 means keep current)
	// This allows SetHMatrixEpsilon to change only eps
	if(eps > 0) m_hacapk_eps = eps;
	if(leaf_size > 0) m_hacapk_leaf_size = leaf_size;
	if(eta > 0) m_hacapk_eta = eta;
	// Invalidate previous statistics when parameters change
	m_hacapk_stats_valid = false;

	if(SendingIsRequired) Send.Int(1);
}

void radTApplication::GetHACApKStats(double* dOut, int* nOut)
{
	if(!m_hacapk_stats_valid)
	{
		*nOut = 0;
		return;
	}

	// H-matrix structure statistics
	dOut[0] = (double)m_hacapk_n_lowrank;
	dOut[1] = (double)m_hacapk_n_dense;
	dOut[2] = (double)m_hacapk_max_rank;
	dOut[3] = (double)m_hacapk_n_leaves;
	dOut[4] = (double)m_hacapk_n_dof;
	dOut[5] = m_hacapk_compression;
	dOut[6] = m_hacapk_build_time;

	// Timing statistics (ELF-compatible)
	dOut[7] = m_timing_hmatrix_build;   // H-matrix construction time [s]
	dOut[8] = m_timing_linear_solve;    // Total BiCGSTAB solve time [s]
	dOut[9] = (double)m_linear_iterations;  // Total BiCGSTAB iterations

	// Memory statistics (ELF-compatible)
	dOut[10] = m_hacapk_memory_mb;       // H-matrix memory [MB]
	dOut[11] = m_hacapk_dense_memory_mb; // Dense matrix memory [MB]

	*nOut = 12;
}
#endif

//-------------------------------------------------------------------------

//-------------------------------------------------------------------------

void radTApplication::GetSolveStats(double* dOut, int* nOut)
{
	if(!m_solve_stats_valid)
	{
		*nOut = 0;
		return;
	}

	// Solve statistics
	dOut[0] = m_solve_t_matrix_build;   // Interaction matrix build time [s]
	dOut[1] = m_solve_t_linear_solve;   // Total linear solver time [s]
	dOut[2] = (double)m_solve_linear_iterations;  // Total linear iterations
	dOut[3] = (double)m_solve_nonl_iterations;    // Total nonlinear iterations

	// Parallelism status
	dOut[4] = 1.0;  // TaskManager enabled
	dOut[5] = (double)m_solve_num_threads;

	// LU decomposition time (Method 0 only)
	dOut[6] = m_solve_t_lu_decomp;

#ifdef RADIA_USE_HACAPK
	// H-matrix timing (Method 2 only) - ELF-compatible
	dOut[7] = m_timing_hmatrix_build;   // H-matrix construction time [s]
	*nOut = 8;
#else
	*nOut = 7;
#endif
}

//-------------------------------------------------------------------------

void radTApplication::ClassifyPoints(int* classification, int* nearest_elem, int n_points,
                                     double* points, int container_handle, double near_threshold)
{
	try
	{
		// Validate handle
		radThg hg;
		if(!ValidateElemKey(container_handle, hg))
		{
			Send.ErrorMessage("Radia::Error003");
			return;
		}

		// Use the solid angle implementation from rad_point_classify
		RadPointClassify::ClassifyPointsFromHandle(
			n_points,
			points,
			container_handle,
			near_threshold,
			classification,
			nearest_elem
		);
	}
	catch(...)
	{
		Send.ErrorMessage("Radia::Error000");
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeFieldBatch(double* B_out, double* H_out, int n_points,
                                        double* points, int container_handle)
{
	try
	{
		// Validate handle and get 3D object pointer
		radThg hg;
		if(!ValidateElemKey(container_handle, hg))
		{
			Send.ErrorMessage("Radia::Error003");
			return;
		}
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr == 0)
		{
			Send.ErrorMessage("Radia::Error003");
			return;
		}

		// Setup IMA context if IMA was used in the last solve AND we're computing field for that model
		// CRITICAL: Only set IMA context if container_handle matches the cached model (m_cached_obj_key)
		// Otherwise, we'd incorrectly add mirror contributions for models that weren't solved with IMA
		bool imaWasSet = false;

		if(m_cached_interact_key > 0 && m_cached_obj_key == container_handle)
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

		// Initialize output arrays to zero (if provided)
		if(B_out) std::memset(B_out, 0, n_points * 3 * sizeof(double));
		if(H_out) std::memset(H_out, 0, n_points * 3 * sizeof(double));

		// Magnetic constant (permeability of free space)
		// B field is stored internally in A/m (same as H), needs conversion to Tesla
		const double Mu0 = 4. * 3.1415926535897932 * 1.e-7; // T*m/A

		// Create field key for B and H computation
		radTFieldKey FieldKey;
		FieldKey.B_ = true;
		FieldKey.H_ = true;

		TVector3d ZeroVect(0., 0., 0.);

		// TaskManager parallelization is safe here because:
		// 1. Each iteration creates its own thread-local radTField object
		// 2. Each iteration writes to different output array indices
		// 3. B_genComp() only reads from the g3dPtr object (no writes)
		// 4. The FieldKey and ZeroVect are copied by value into each thread's radTField
		ngcore::ParallelFor(ngcore::IntRange(n_points), [&](size_t i) {
			TVector3d pt;
			pt.x = points[i * 3 + 0];
			pt.y = points[i * 3 + 1];
			pt.z = points[i * 3 + 2];

			// Create thread-local field object for each point
			radTFieldKey localFieldKey;
			localFieldKey.B_ = true;
			localFieldKey.H_ = true;
			TVector3d localZeroVect(0., 0., 0.);

			radTField Field(localFieldKey, CompCriterium, pt, localZeroVect, localZeroVect, localZeroVect, localZeroVect, 0.);
			g3dPtr->B_genComp(&Field);

			// Convert B from A/m to Tesla (multiply by Mu0)
			if(B_out) {
				B_out[i * 3 + 0] = Field.B.x * Mu0;
				B_out[i * 3 + 1] = Field.B.y * Mu0;
				B_out[i * 3 + 2] = Field.B.z * Mu0;
			}

			// H is already in A/m
			if(H_out) {
				H_out[i * 3 + 0] = Field.H.x;
				H_out[i * 3 + 1] = Field.H.y;
				H_out[i * 3 + 2] = Field.H.z;
			}
		});

		// Clear IMA context after computation
		if(imaWasSet) RadIMAFieldContext::Clear();
	}
	catch(...)
	{
		// Clear IMA context on error
		RadIMAFieldContext::Clear();
		Send.ErrorMessage("Radia::Error000");
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeScalarPotentialBatch(double* phi_out, int n_points,
                                                  double* points, int container_handle)
{
	try
	{
		// Validate handle and get 3D object pointer
		radThg hg;
		if(!ValidateElemKey(container_handle, hg))
		{
			Send.ErrorMessage("Radia::Error003");
			return;
		}
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr == 0)
		{
			Send.ErrorMessage("Radia::Error003");
			return;
		}

		// Setup IMA context if IMA was used in the last solve AND we're computing field for that model
		// CRITICAL: Only set IMA context if container_handle matches the cached model (m_cached_obj_key)
		// Otherwise, we'd incorrectly add mirror contributions for models that weren't solved with IMA
		bool imaWasSet = false;
		if(m_cached_interact_key > 0 && m_cached_obj_key == container_handle)
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

		// Initialize output array to zero
		std::memset(phi_out, 0, n_points * sizeof(double));

		// Create field key for Phi computation
		radTFieldKey FieldKey;
		FieldKey.Phi_ = true;

		TVector3d ZeroVect(0., 0., 0.);

		// Compute scalar potential at each point (TaskManager parallelized)
		ngcore::ParallelFor(ngcore::IntRange(n_points), [&](size_t i) {
			TVector3d pt;
			pt.x = points[i * 3 + 0];
			pt.y = points[i * 3 + 1];
			pt.z = points[i * 3 + 2];

			radTFieldKey localFieldKey;
			localFieldKey.Phi_ = true;
			TVector3d localZeroVect(0., 0., 0.);

			radTField Field(localFieldKey, CompCriterium, pt, localZeroVect, localZeroVect, localZeroVect, localZeroVect, 0.);
			g3dPtr->B_genComp(&Field);

			phi_out[i] = Field.Phi;
		});

		// Clear IMA context after computation
		if(imaWasSet) RadIMAFieldContext::Clear();
	}
	catch(...)
	{
		// Clear IMA context on error
		RadIMAFieldContext::Clear();
		Send.ErrorMessage("Radia::Error000");
	}
}

//-------------------------------------------------------------------------

void radTApplication::ComputeVectorPotentialBatch(double* A_out, int n_points,
                                                  double* points, int container_handle)
{
	try
	{
		// Validate handle and get 3D object pointer
		radThg hg;
		if(!ValidateElemKey(container_handle, hg))
		{
			Send.ErrorMessage("Radia::Error003");
			return;
		}
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep);
		if(g3dPtr == 0)
		{
			Send.ErrorMessage("Radia::Error003");
			return;
		}

		// Initialize output array to zero
		std::memset(A_out, 0, n_points * 3 * sizeof(double));

		// Create field key for A computation
		radTFieldKey FieldKey;
		FieldKey.A_ = true;

		TVector3d ZeroVect(0., 0., 0.);

		// Compute vector potential at each point (TaskManager parallelized)
		ngcore::ParallelFor(ngcore::IntRange(n_points), [&](size_t i) {
			TVector3d pt;
			pt.x = points[i * 3 + 0];
			pt.y = points[i * 3 + 1];
			pt.z = points[i * 3 + 2];

			radTFieldKey localFieldKey;
			localFieldKey.A_ = true;
			TVector3d localZeroVect(0., 0., 0.);

			radTField Field(localFieldKey, CompCriterium, pt, localZeroVect, localZeroVect, localZeroVect, localZeroVect, 0.);
			g3dPtr->B_genComp(&Field);

			A_out[i * 3 + 0] = Field.A.x;
			A_out[i * 3 + 1] = Field.A.y;
			A_out[i * 3 + 2] = Field.A.z;
		});
	}
	catch(...)
	{
		Send.ErrorMessage("Radia::Error000");
	}
}

//-------------------------------------------------------------------------

radTInteraction* radTApplication::GetInteractionByKey(int interactKey)
{
	radThg hg;
	if(!ValidateElemKey(interactKey, hg))
	{
		return nullptr;
	}
	return Cast.InteractCast(hg.rep);
}

//=========================================================================
// radTEnergyHysteresisMaterial implementation
// Energy-based vector hysteresis (Francois-Lavet / Egger formulation)
//=========================================================================

static const double EHYST_MU_0 = 4.0 * 3.14159265358979323846 * 1.0e-7;
static const double EHYST_NU_0 = 1.0 / EHYST_MU_0;

//-------------------------------------------------------------------------
// Table interpolation helper: binary search + linear interpolation
//-------------------------------------------------------------------------

static double InterpolateTable(const std::vector<double>& x,
                                const std::vector<double>& y,
                                int n, double xi)
{
	if(n <= 0) return 0.0;
	if(xi <= x[0]) return y[0];
	if(xi >= x[n-1]) return y[n-1];  // clamp (flat extrapolation)
	// Binary search
	int lo = 0, hi = n - 1;
	while(hi - lo > 1)
	{
		int mid = (lo + hi) / 2;
		if(x[mid] <= xi) lo = mid; else hi = mid;
	}
	double t = (xi - x[lo]) / (x[hi] - x[lo]);
	return y[lo] + t * (y[hi] - y[lo]);
}

//-------------------------------------------------------------------------
// Precompute U (integral) and df (derivative) tables from r, f data
//-------------------------------------------------------------------------

void radTEnergyHysteresisMaterial::PrecomputeOperatorTable(OperatorTable& tab)
{
	int n = tab.n;
	tab.U.resize(n);
	tab.df.resize(n);

	// U[0] = 0, trapezoidal integration: U[i] = U[i-1] + 0.5*(f[i-1]+f[i])*(r[i]-r[i-1])
	tab.U[0] = 0.0;
	for(int i = 1; i < n; i++)
		tab.U[i] = tab.U[i-1] + 0.5 * (tab.f[i-1] + tab.f[i]) * (tab.r[i] - tab.r[i-1]);

	// df: finite differences
	if(n == 1)
	{
		tab.df[0] = 0.0;
	}
	else if(n == 2)
	{
		double d = (tab.f[1] - tab.f[0]) / (tab.r[1] - tab.r[0]);
		tab.df[0] = d;
		tab.df[1] = d;
	}
	else
	{
		// Forward difference for first point
		tab.df[0] = (tab.f[1] - tab.f[0]) / (tab.r[1] - tab.r[0]);
		// Central differences for interior
		for(int i = 1; i < n - 1; i++)
			tab.df[i] = (tab.f[i+1] - tab.f[i-1]) / (tab.r[i+1] - tab.r[i-1]);
		// Backward difference for last point
		tab.df[n-1] = (tab.f[n-1] - tab.f[n-2]) / (tab.r[n-1] - tab.r[n-2]);
	}
}

//-------------------------------------------------------------------------
// Constructor from table-based shape functions
//-------------------------------------------------------------------------

radTEnergyHysteresisMaterial::radTEnergyHysteresisMaterial(
	int K, const double* chi,
	const std::vector<std::vector<double>>& r_tables,
	const std::vector<std::vector<double>>& f_tables,
	double eps)
	: m_K(K), m_eps(eps), m_last_chi(0), m_last_chi_d(0), m_has_result(false)
{
	m_chi.resize(K);
	m_tables.resize(K);
	for(int k = 0; k < K; k++)
	{
		m_chi[k] = chi[k];
		auto& tab = m_tables[k];
		int n = (int)r_tables[k].size();
		tab.n = n;
		tab.r = r_tables[k];
		tab.f = f_tables[k];
		tab.r_max = (n > 0) ? tab.r[n-1] : 0.0;
		PrecomputeOperatorTable(tab);
	}

	TVector3d zero(0, 0, 0);
	m_Jk_prev.resize(K, zero);
	m_Jk_pinning.resize(K, zero);
	m_Jk_current.resize(K, zero);
	m_last_H = zero;
	m_last_B = zero;
}

//-------------------------------------------------------------------------
// Serialization constructor
//-------------------------------------------------------------------------

radTEnergyHysteresisMaterial::radTEnergyHysteresisMaterial(CAuxBinStrVect& inStr)
	: m_last_chi(0), m_last_chi_d(0), m_has_result(false)
{
	DumpBinParse_Material(inStr);
	inStr >> m_K;
	inStr >> m_eps;

	m_chi.resize(m_K);
	m_tables.resize(m_K);
	for(int k = 0; k < m_K; k++)
	{
		inStr >> m_chi[k];
		auto& tab = m_tables[k];
		inStr >> tab.n;
		tab.r.resize(tab.n);
		tab.f.resize(tab.n);
		for(int i = 0; i < tab.n; i++) inStr >> tab.r[i];
		for(int i = 0; i < tab.n; i++) inStr >> tab.f[i];
		tab.r_max = (tab.n > 0) ? tab.r[tab.n-1] : 0.0;
		PrecomputeOperatorTable(tab);  // Recompute U and df
	}

	TVector3d zero(0, 0, 0);
	m_Jk_prev.resize(m_K, zero);
	m_Jk_pinning.resize(m_K, zero);
	m_Jk_current.resize(m_K, zero);

	// Deserialize state
	for(int k = 0; k < m_K; k++) inStr >> m_Jk_prev[k];
	m_Jk_pinning = m_Jk_prev;
	m_last_H = zero;
	m_last_B = zero;
}

//-------------------------------------------------------------------------
// Serialization
//-------------------------------------------------------------------------

void radTEnergyHysteresisMaterial::DumpBin(CAuxBinStrVect& oStr,
	std::vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers,
	int& gUniqueMapKey, int elemKey)
{
	vElemKeysOut.push_back(elemKey);
	int MatType = Type_Material();
	oStr << MatType;
	DumpBin_Material(oStr);

	oStr << m_K;
	oStr << m_eps;
	for(int k = 0; k < m_K; k++)
	{
		oStr << m_chi[k];
		oStr << m_tables[k].n;
		for(int i = 0; i < m_tables[k].n; i++) oStr << m_tables[k].r[i];
		for(int i = 0; i < m_tables[k].n; i++) oStr << m_tables[k].f[i];
	}
	for(int k = 0; k < m_K; k++) oStr << m_Jk_prev[k];
}

//-------------------------------------------------------------------------
// Internal energy U_k(|J|) and derivatives (table interpolation)
//-------------------------------------------------------------------------

double radTEnergyHysteresisMaterial::Uk(int k, double J_mag) const
{
	if(J_mag < 1e-30) return 0.0;
	return InterpolateTable(m_tables[k].r, m_tables[k].U, m_tables[k].n, J_mag);
}

double radTEnergyHysteresisMaterial::dUk(int k, double J_mag) const
{
	if(J_mag < 1e-30) return 0.0;
	return InterpolateTable(m_tables[k].r, m_tables[k].f, m_tables[k].n, J_mag);
}

double radTEnergyHysteresisMaterial::d2Uk(int k, double J_mag) const
{
	return InterpolateTable(m_tables[k].r, m_tables[k].df, m_tables[k].n, J_mag);
}

//-------------------------------------------------------------------------
// Vector gradient and Hessian of U_k
//-------------------------------------------------------------------------

TVector3d radTEnergyHysteresisMaterial::GradUk(int k, const TVector3d& J) const
{
	double J_mag = NormAbs(J);
	if(J_mag < 1e-30) return TVector3d(0, 0, 0);
	return (dUk(k, J_mag) / J_mag) * J;
}

TMatrix3d radTEnergyHysteresisMaterial::HessUk(int k, const TVector3d& J) const
{
	double J_mag = NormAbs(J);
	if(J_mag < 1e-30) return d2Uk(k, 0.0) * Eye();

	TVector3d e = (1.0 / J_mag) * J;
	double du = dUk(k, J_mag);
	double d2u = d2Uk(k, J_mag);
	double du_over_r = du / J_mag;

	// H_U = (d2u - du/|J|) * e*e^T + (du/|J|) * I
	return (d2u - du_over_r) * OuterProduct(e, e) + du_over_r * Eye();
}

//-------------------------------------------------------------------------
// Regularized norm |x|_eps and derivatives
//-------------------------------------------------------------------------

double radTEnergyHysteresisMaterial::NormEps(const TVector3d& x) const
{
	return sqrt(x * x + m_eps);  // x*x is dot product via operator*
}

TVector3d radTEnergyHysteresisMaterial::GradNormEps(const TVector3d& x) const
{
	return (1.0 / NormEps(x)) * x;
}

TMatrix3d radTEnergyHysteresisMaterial::HessNormEps(const TVector3d& x) const
{
	double n = NormEps(x);
	double inv_n = 1.0 / n;
	// (I - x*x^T / n^2) / n
	return inv_n * (Eye() - (1.0 / (n * n)) * OuterProduct(x, x));
}

//-------------------------------------------------------------------------
// Objective functions
//-------------------------------------------------------------------------

double radTEnergyHysteresisMaterial::ObjectiveForwardK(
	int k, const TVector3d& J, const TVector3d& H) const
{
	TVector3d diff = J - m_Jk_prev[k];
	return Uk(k, NormAbs(J)) - (H * J) + m_chi[k] * NormEps(diff);
}

double radTEnergyHysteresisMaterial::ObjectiveForward(
	const TVector3d& B, const std::vector<TVector3d>& Jk_list) const
{
	TVector3d J_sum(0, 0, 0);
	for(int k = 0; k < m_K; k++) J_sum += Jk_list[k];

	TVector3d diff_B = B - J_sum;
	double G = 0.5 * EHYST_NU_0 * (diff_B * diff_B);

	for(int k = 0; k < m_K; k++)
	{
		TVector3d diff = Jk_list[k] - m_Jk_prev[k];
		G += Uk(k, NormAbs(Jk_list[k]));
		G += m_chi[k] * NormEps(diff);
	}
	return G;
}

//-------------------------------------------------------------------------
// Newton solver for single J_k (used by Inverse: H -> B)
//-------------------------------------------------------------------------

TVector3d radTEnergyHysteresisMaterial::SolveInverseK(int k, const TVector3d& H) const
{
	TVector3d Jk = m_Jk_prev[k];
	const int max_iter = 30;
	const double tol = 1e-12;

	for(int it = 0; it < max_iter; it++)
	{
		TVector3d diff = Jk - m_Jk_prev[k];
		TVector3d grad = GradUk(k, Jk) - H + m_chi[k] * GradNormEps(diff);
		double grad_norm = NormAbs(grad);
		if(grad_norm < tol) break;

		TMatrix3d hess = HessUk(k, Jk) + m_chi[k] * HessNormEps(diff);
		TVector3d dJ = Solve3x3(hess, (-1.0) * grad);

		// Armijo backtracking
		double tau = 1.0;
		double F_curr = ObjectiveForwardK(k, Jk, H);
		double dir_deriv = grad * dJ;  // directional derivative
		for(int ls = 0; ls < 20; ls++)
		{
			double F_new = ObjectiveForwardK(k, Jk + tau * dJ, H);
			if(F_new <= F_curr + 0.1 * tau * dir_deriv) break;
			tau *= 0.5;
		}
		Jk += tau * dJ;
		// Clamp to saturation bound
		double Jmax = 0.9999 * m_tables[k].r_max;
		if(Jk.x > Jmax) Jk.x = Jmax;
		else if(Jk.x < -Jmax) Jk.x = -Jmax;
		if(Jk.y > Jmax) Jk.y = Jmax;
		else if(Jk.y < -Jmax) Jk.y = -Jmax;
		if(Jk.z > Jmax) Jk.z = Jmax;
		else if(Jk.z < -Jmax) Jk.z = -Jmax;
	}
	return Jk;
}

//-------------------------------------------------------------------------
// Inverse operator: H -> B (energy-based approximation of B-input Play)
// Each J_k solved independently -- no Schur needed
//-------------------------------------------------------------------------

TVector3d radTEnergyHysteresisMaterial::Inverse(const TVector3d& H)
{
	// Save pinning reference before update (for Jacobian)
	for(int k = 0; k < m_K; k++)
		m_Jk_pinning[k] = m_Jk_prev[k];

	TVector3d J_total(0, 0, 0);
	for(int k = 0; k < m_K; k++)
	{
		TVector3d Jk = SolveInverseK(k, H);
		m_Jk_current[k] = Jk;
		J_total += Jk;
	}

	TVector3d B = EHYST_MU_0 * H + J_total;

	// Update state
	for(int k = 0; k < m_K; k++)
		m_Jk_prev[k] = m_Jk_current[k];

	// Cache result
	m_last_H = H;
	m_last_B = B;
	m_has_result = true;

	return B;
}

//-------------------------------------------------------------------------
// Forward operator: B -> H (natural for B-input Play, Schur complement)
//-------------------------------------------------------------------------

TVector3d radTEnergyHysteresisMaterial::Forward(const TVector3d& B)
{
	// Save pinning reference before update
	for(int k = 0; k < m_K; k++)
		m_Jk_pinning[k] = m_Jk_prev[k];

	std::vector<TVector3d> Jk_list(m_K);
	for(int k = 0; k < m_K; k++)
		Jk_list[k] = m_Jk_prev[k];

	const int max_iter_inv = 100;
	TMatrix3d I3 = Eye();

	for(int n_iter = 0; n_iter < max_iter_inv; n_iter++)
	{
		// Compute common residual
		TVector3d J_sum(0, 0, 0);
		for(int k = 0; k < m_K; k++) J_sum += Jk_list[k];
		TVector3d residual_common = (-EHYST_NU_0) * (B - J_sum);

		// Compute gradients and Hessians
		std::vector<TVector3d> grad_list(m_K);
		std::vector<TMatrix3d> hess_list(m_K);
		double max_grad = 0.0;

		for(int k = 0; k < m_K; k++)
		{
			TVector3d diff = Jk_list[k] - m_Jk_prev[k];
			grad_list[k] = residual_common + GradUk(k, Jk_list[k])
			               + m_chi[k] * GradNormEps(diff);
			hess_list[k] = EHYST_NU_0 * I3 + HessUk(k, Jk_list[k])
			               + m_chi[k] * HessNormEps(diff);

			double gn = NormAbs(grad_list[k]);
			if(gn > max_grad) max_grad = gn;
		}

		if(max_grad < 1e-12) break;

		// Schur complement
		std::vector<TMatrix3d> hk_priv_inv(m_K);
		for(int k = 0; k < m_K; k++)
		{
			TMatrix3d hk_priv = hess_list[k] - EHYST_NU_0 * I3;
			hk_priv_inv[k] = Matrix3d_inv(hk_priv);
		}

		TMatrix3d S = I3;
		TVector3d rhs(0, 0, 0);
		for(int k = 0; k < m_K; k++)
		{
			S += EHYST_NU_0 * (hk_priv_inv[k]);  // Use operator* for matrix*matrix
			rhs -= hk_priv_inv[k] * grad_list[k];
		}

		TVector3d delta = Solve3x3(S, rhs);

		// Compute step for each J_k
		std::vector<TVector3d> dJk_list(m_K);
		double dir_deriv = 0.0;
		for(int k = 0; k < m_K; k++)
		{
			dJk_list[k] = hk_priv_inv[k] * ((-1.0) * grad_list[k] - EHYST_NU_0 * delta);
			dir_deriv += grad_list[k] * dJk_list[k];
		}

		// Armijo backtracking on total objective
		double G_curr = ObjectiveForward(B, Jk_list);
		double tau = 1.0;
		for(int ls = 0; ls < 30; ls++)
		{
			std::vector<TVector3d> Jk_trial(m_K);
			for(int k = 0; k < m_K; k++)
				Jk_trial[k] = Jk_list[k] + tau * dJk_list[k];

			double G_trial = ObjectiveForward(B, Jk_trial);
			if(G_trial <= G_curr + 0.1 * tau * dir_deriv) break;
			tau *= 0.5;
		}

		for(int k = 0; k < m_K; k++)
		{
			Jk_list[k] += tau * dJk_list[k];
			// Clamp each component to saturation bound (matches Python L-BFGS-B bounds)
			double Jmax = 0.9999 * m_tables[k].r_max;
			if(Jk_list[k].x > Jmax) Jk_list[k].x = Jmax;
			else if(Jk_list[k].x < -Jmax) Jk_list[k].x = -Jmax;
			if(Jk_list[k].y > Jmax) Jk_list[k].y = Jmax;
			else if(Jk_list[k].y < -Jmax) Jk_list[k].y = -Jmax;
			if(Jk_list[k].z > Jmax) Jk_list[k].z = Jmax;
			else if(Jk_list[k].z < -Jmax) Jk_list[k].z = -Jmax;
		}
	}

	// Compute H = nu_0 * (B - sum J_k)
	TVector3d J_total(0, 0, 0);
	for(int k = 0; k < m_K; k++) J_total += Jk_list[k];
	TVector3d H = EHYST_NU_0 * (B - J_total);

	// Update state
	for(int k = 0; k < m_K; k++)
	{
		m_Jk_current[k] = Jk_list[k];
		m_Jk_prev[k] = Jk_list[k];
	}

	m_last_H = H;
	m_last_B = B;
	m_has_result = true;

	return H;
}

//-------------------------------------------------------------------------
// Jacobian dB/dH
//-------------------------------------------------------------------------

void radTEnergyHysteresisMaterial::ComputeJacobian(TMatrix3d& dBdH, double& chi_d) const
{
	dBdH = EHYST_MU_0 * Eye();
	for(int k = 0; k < m_K; k++)
	{
		TVector3d Jk = m_Jk_current[k];
		TVector3d diff = Jk - m_Jk_pinning[k];
		TMatrix3d hess = HessUk(k, Jk) + m_chi[k] * HessNormEps(diff);
		dBdH += Matrix3d_inv(hess);
	}

	// Scalar differential chi = trace(dBdH) / (3 * mu_0) - 1
	double trace = dBdH.Str0.x + dBdH.Str1.y + dBdH.Str2.z;
	chi_d = trace / (3.0 * EHYST_MU_0) - 1.0;
}

//-------------------------------------------------------------------------
// M(H) virtual override
//-------------------------------------------------------------------------

TVector3d radTEnergyHysteresisMaterial::M(const TVector3d& H)
{
	TVector3d B = Inverse(H);  // H -> B
	// M = B/mu_0 - H
	return (1.0 / EHYST_MU_0) * B - H;
}

//-------------------------------------------------------------------------
// DefineInstantKsiTensor virtual override
//-------------------------------------------------------------------------

void radTEnergyHysteresisMaterial::DefineInstantKsiTensor(
	const TVector3d& H, TMatrix3d& KsiTensor, TVector3d& Mr)
{
	TVector3d B = Inverse(H);  // H -> B
	double H_mag = NormAbs(H);
	double B_mag = NormAbs(B);

	double chi_scalar = 0.0;
	if(H_mag > 1e-30)
		chi_scalar = B_mag / (EHYST_MU_0 * H_mag) - 1.0;

	KsiTensor = chi_scalar * Eye();
	Mr = TVector3d(0, 0, 0);  // No remanence in this formulation
}

//-------------------------------------------------------------------------
// Solver integration methods
//-------------------------------------------------------------------------

double radTEnergyHysteresisMaterial::ComputeChiFromH(const TVector3d& H)
{
	TVector3d B = Inverse(H);  // H -> B
	double H_mag = NormAbs(H);
	double B_mag = NormAbs(B);

	if(H_mag < 1e-30) return GetInitialChi_ELF_Style();
	m_last_chi = B_mag / (EHYST_MU_0 * H_mag) - 1.0;
	return m_last_chi;
}

double radTEnergyHysteresisMaterial::ComputeChiDualMethod(
	double H_mag, double mu_old, double relax)
{
	if(H_mag < 1e-30) return GetInitialChi_ELF_Style();

	// Reconstruct H vector using last known direction
	TVector3d H;
	double last_H_mag = NormAbs(m_last_H);
	if(last_H_mag > 1e-30)
		H = (H_mag / last_H_mag) * m_last_H;
	else
		H = TVector3d(H_mag, 0, 0);  // Default direction

	TVector3d B = Inverse(H);  // H -> B
	double B_mag = NormAbs(B);

	m_last_chi = B_mag / (EHYST_MU_0 * H_mag) - 1.0;
	if(m_last_chi < 0) m_last_chi = 0;
	return m_last_chi;
}

double radTEnergyHysteresisMaterial::ComputeDifferentialChi(double H_mag)
{
	TMatrix3d dBdH;
	double chi_d;
	ComputeJacobian(dBdH, chi_d);
	m_last_chi_d = chi_d;
	return chi_d;
}

//-------------------------------------------------------------------------
// Application method
//-------------------------------------------------------------------------

int radTApplication::SetEnergyHysteresisMaterial(
	int K, const double* chi,
	const std::vector<std::vector<double>>& r_tables,
	const std::vector<std::vector<double>>& f_tables,
	double eps)
{
	radTEnergyHysteresisMaterial* MaterPtr = nullptr;
	try
	{
		MaterPtr = new radTEnergyHysteresisMaterial(K, chi, r_tables, f_tables, eps);
		if(MaterPtr == nullptr) { Send.ErrorMessage("Radia::Error900"); return 0; }

		radThg hg(MaterPtr);
		MaterPtr = nullptr;
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		if(MaterPtr) delete MaterPtr;
		Initialize();
		return 0;
	}
}

//=========================================================================
// radTPlayHysteresisMaterial implementation
//=========================================================================

static const double PLAY_MU_0 = 4.0 * 3.14159265358979323846 * 1.0e-7;

void radTPlayHysteresisMaterial::PrecomputePlayTable(PlayTable& tab)
{
	int n = tab.n;
	tab.df.resize(n);
	if(n == 1) {
		tab.df[0] = 0.0;
	} else if(n == 2) {
		double d = (tab.f[1] - tab.f[0]) / (tab.r[1] - tab.r[0]);
		tab.df[0] = d;
		tab.df[1] = d;
	} else {
		tab.df[0] = (tab.f[1] - tab.f[0]) / (tab.r[1] - tab.r[0]);
		for(int i = 1; i < n - 1; i++)
			tab.df[i] = (tab.f[i+1] - tab.f[i-1]) / (tab.r[i+1] - tab.r[i-1]);
		tab.df[n-1] = (tab.f[n-1] - tab.f[n-2]) / (tab.r[n-1] - tab.r[n-2]);
	}
}

void radTPlayHysteresisMaterial::ComputeMonotoneLimits()
{
	// Scan virgin initial curve (1D, all p_k start at 0) to find where dH/dB <= 0.
	// On virgin ascending: p_k = max(0, B - eta_k), dp_k/dB = 1 if B > eta_k.
	// H(B) = sum_{k: B > eta_k} f_k(B - eta_k)
	// dH/dB = sum_{k: B > eta_k} f'_k(B - eta_k)
	double B_sat = GetBsaturation();
	if(B_sat < 1e-10) { m_B_mono_max = 1.0; m_H_mono_max = 0; return; }

	int n_scan = 500;
	double dB_step = B_sat / n_scan;
	m_B_mono_max = B_sat;
	m_H_mono_max = 0;

	double H_prev = 0;
	for(int i = 1; i <= n_scan; i++)
	{
		double B_val = i * dB_step;
		double H_val = 0;
		double dHdB_val = 0;

		for(int k = 0; k < m_K; k++)
		{
			double pk_val = B_val - m_eta[k];
			if(pk_val <= 1e-30) continue;  // stuck operator
			H_val += fk(k, pk_val);
			dHdB_val += dfk(k, pk_val);
		}

		if(dHdB_val <= 0 && i > 1)
		{
			// Non-monotone detected: previous point is the safe limit
			m_B_mono_max = (i - 1) * dB_step;
			m_H_mono_max = H_prev;
			return;
		}
		H_prev = H_val;
	}
	// No non-monotonicity found: full range is monotone
	m_H_mono_max = H_prev;
}

radTPlayHysteresisMaterial::radTPlayHysteresisMaterial(
	int K, const double* eta,
	const std::vector<std::vector<double>>& r_tables,
	const std::vector<std::vector<double>>& f_tables)
	: m_K(K), m_B_mono_max(1.0), m_H_mono_max(0),
	  m_last_chi(0), m_last_chi_d(0), m_has_result(false)
{
	m_eta.resize(K);
	m_tables.resize(K);
	for(int k = 0; k < K; k++)
	{
		m_eta[k] = eta[k];
		auto& tab = m_tables[k];
		int n = (int)r_tables[k].size();
		tab.n = n;
		tab.r = r_tables[k];
		tab.f = f_tables[k];
		tab.r_max = (n > 0) ? tab.r[n-1] : 0.0;
		PrecomputePlayTable(tab);
	}

	TVector3d zero(0, 0, 0);
	m_pk_prev.resize(K, zero);
	m_pk_pinning.resize(K, zero);
	m_pk_current.resize(K, zero);
	m_last_H = zero;
	m_last_B = zero;
	m_last_dHdB = Eye();

	ComputeMonotoneLimits();
}

double radTPlayHysteresisMaterial::fk(int k, double r_mag) const
{
	return InterpolateTable(m_tables[k].r, m_tables[k].f, m_tables[k].n, r_mag);
}

double radTPlayHysteresisMaterial::dfk(int k, double r_mag) const
{
	return InterpolateTable(m_tables[k].r, m_tables[k].df, m_tables[k].n, r_mag);
}

//-------------------------------------------------------------------------
// Forward: B -> H, O(K) direct evaluation
// Also computes and caches analytical Jacobian dH/dB
//-------------------------------------------------------------------------
TVector3d radTPlayHysteresisMaterial::Forward(const TVector3d& B)
{
	TVector3d H(0, 0, 0);
	TMatrix3d dHdB(TVector3d(0,0,0), TVector3d(0,0,0), TVector3d(0,0,0));

	for(int k = 0; k < m_K; k++)
	{
		double eta_k = m_eta[k];
		TVector3d q = B - m_pk_pinning[k];
		double q_mag = sqrt(q.x*q.x + q.y*q.y + q.z*q.z);

		TVector3d pk;
		TMatrix3d dpk_dB(TVector3d(0,0,0), TVector3d(0,0,0), TVector3d(0,0,0));

		if(eta_k < 1e-30)
		{
			// eta=0: p_k always tracks B exactly
			pk = B;
			dpk_dB = Eye();
		}
		else if(q_mag <= eta_k)
		{
			// Stuck: p_k doesn't move
			pk = m_pk_pinning[k];
			// dpk_dB = 0 (already zeros)
		}
		else
		{
			// Following: p_k = B - eta_k * q / |q|
			double ratio = eta_k / q_mag;
			pk.x = B.x - ratio * q.x;
			pk.y = B.y - ratio * q.y;
			pk.z = B.z - ratio * q.z;

			// dp/dB = (1 - eta/|q|)*I + (eta/|q|) * q_hat * q_hat^T
			TVector3d q_hat(q.x/q_mag, q.y/q_mag, q.z/q_mag);
			double c1 = 1.0 - ratio;
			dpk_dB = c1 * Eye() + ratio * OuterProduct(q_hat, q_hat);
		}
		m_pk_current[k] = pk;

		double pk_mag = sqrt(pk.x*pk.x + pk.y*pk.y + pk.z*pk.z);
		if(pk_mag > 1e-30)
		{
			double f_val = fk(k, pk_mag);
			double df_val = dfk(k, pk_mag);
			TVector3d pk_hat(pk.x/pk_mag, pk.y/pk_mag, pk.z/pk_mag);

			// H contribution: f_k(|p_k|) * p_k / |p_k|
			H.x += f_val * pk_hat.x;
			H.y += f_val * pk_hat.y;
			H.z += f_val * pk_hat.z;

			// dh_k/dp_k = df_k * p_hat*p_hat^T + (f_k/|p_k|)*(I - p_hat*p_hat^T)
			double f_over_r = f_val / pk_mag;
			TMatrix3d pp = OuterProduct(pk_hat, pk_hat);
			TMatrix3d I_pp = Eye() - pp;
			TMatrix3d dhk_dpk = df_val * pp + f_over_r * I_pp;

			// dH/dB += dh_k/dp_k * dp_k/dB (chain rule)
			dHdB = dHdB + dhk_dpk * dpk_dB;
		}
	}

	m_last_dHdB = dHdB;
	m_last_H = H;
	m_last_B = B;
	m_has_result = true;
	return H;
}

//-------------------------------------------------------------------------
// Inverse: H -> B, Newton iteration with analytical Jacobian
//-------------------------------------------------------------------------
TVector3d radTPlayHysteresisMaterial::Inverse(const TVector3d& H_target)
{
	// Ensure pinning is set from committed state
	for(int k = 0; k < m_K; k++)
		m_pk_pinning[k] = m_pk_prev[k];

	double H_target_mag = sqrt(H_target.x*H_target.x + H_target.y*H_target.y + H_target.z*H_target.z);

	// Monotone constraint: B clamped to [0, m_B_mono_max] where dH/dB > 0.
	// For H > m_H_mono_max (beyond model range), return B_mono_max directly.
	double B_max = m_B_mono_max;
	double max_dB_step = 0.05 * B_max;  // max |dB| per Newton iteration

	// Early return: H_target beyond monotone range -> saturated at B_mono_max
	if(H_target_mag > m_H_mono_max && m_H_mono_max > 1e-10)
	{
		TVector3d B_sat_dir;
		if(H_target_mag > 1e-20)
		{
			double sc = m_B_mono_max / H_target_mag;
			B_sat_dir.x = sc * H_target.x;
			B_sat_dir.y = sc * H_target.y;
			B_sat_dir.z = sc * H_target.z;
		}
		else
		{
			B_sat_dir.x = B_sat_dir.y = B_sat_dir.z = 0;
		}
		Forward(B_sat_dir);
		// Auto-commit
		for(int k = 0; k < m_K; k++)
			m_pk_prev[k] = m_pk_current[k];
		return m_last_B;
	}

	// Initial guess
	TVector3d B;
	if(m_has_result && H_target_mag > 1e-20)
	{
		B = m_last_B;
	}
	else
	{
		// B = mu_0 * (1 + chi_init) * H
		double chi_init = GetInitialChi_ELF_Style();
		if(chi_init < 1.0) chi_init = 1.0;
		double mu_r = 1.0 + chi_init;
		B.x = PLAY_MU_0 * mu_r * H_target.x;
		B.y = PLAY_MU_0 * mu_r * H_target.y;
		B.z = PLAY_MU_0 * mu_r * H_target.z;
	}

	// Clamp initial guess
	double B_init_mag = sqrt(B.x*B.x + B.y*B.y + B.z*B.z);
	if(B_init_mag > B_max)
	{
		double scale = B_max / B_init_mag;
		B.x *= scale; B.y *= scale; B.z *= scale;
	}

	if(H_target_mag < 1e-20)
	{
		// H=0: B=0 for virgin state, otherwise use Forward at B=0
		if(!m_has_result)
		{
			B.x = B.y = B.z = 0;
			Forward(B);
			// Auto-commit
			for(int k = 0; k < m_K; k++)
				m_pk_prev[k] = m_pk_current[k];
			return m_last_B;
		}
	}

	const int max_iter = 100;
	const double tol = 1e-10;

	// Track best solution (lowest residual) for fallback
	TVector3d B_best = B;
	double best_res = 1e30;

	for(int it = 0; it < max_iter; it++)
	{
		TVector3d H_comp = Forward(B);
		TVector3d res(H_comp.x - H_target.x, H_comp.y - H_target.y, H_comp.z - H_target.z);
		double res_norm = sqrt(res.x*res.x + res.y*res.y + res.z*res.z);

		if(res_norm < best_res)
		{
			best_res = res_norm;
			B_best = B;
		}
		if(res_norm < tol) break;

		double trace_dHdB = m_last_dHdB.Str0.x + m_last_dHdB.Str1.y + m_last_dHdB.Str2.z;

		// Non-monotone detection: if dH/dB becomes negative despite B_mono_max clamp
		// (can happen on non-virgin branches), retreat toward lower |B|
		if(trace_dHdB < -1e-6)
		{
			double B_mag = sqrt(B.x*B.x + B.y*B.y + B.z*B.z);
			double retreat = 0.05 * B_max;
			if(B_mag > retreat)
			{
				double sc = (B_mag - retreat) / B_mag;
				B.x *= sc; B.y *= sc; B.z *= sc;
			}
			continue;
		}

		double dHdB_scale = fabs(trace_dHdB) / 3.0;

		TVector3d dB;
		if(dHdB_scale > 1e-6)
		{
			// Newton step: dB = -inv(dH/dB) * residual
			TMatrix3d dHdB_reg = m_last_dHdB;
			double reg = 1e-8 * dHdB_scale;
			dHdB_reg.Str0.x += reg;
			dHdB_reg.Str1.y += reg;
			dHdB_reg.Str2.z += reg;

			TMatrix3d dBdH = Matrix3d_inv(dHdB_reg);
			dB = (-1.0) * (dBdH * res);
		}
		else
		{
			// Jacobian too small: steepest descent
			double step = (H_target_mag > 1.0) ? PLAY_MU_0 * H_target_mag * 0.05 : 0.005;
			if(H_target_mag > 1e-20)
			{
				TVector3d H_hat(H_target.x/H_target_mag, H_target.y/H_target_mag, H_target.z/H_target_mag);
				double res_proj = res.x*H_hat.x + res.y*H_hat.y + res.z*H_hat.z;
				double sign = (res_proj > 0) ? -1.0 : 1.0;
				dB.x = sign * step * H_hat.x;
				dB.y = sign * step * H_hat.y;
				dB.z = sign * step * H_hat.z;
			}
			else
			{
				dB.x = (res.x > 0) ? -step : step;
				dB.y = (res.y > 0) ? -step : step;
				dB.z = (res.z > 0) ? -step : step;
			}
		}

		// NaN check
		if(dB.x != dB.x || dB.y != dB.y || dB.z != dB.z)
		{
			double chi = GetInitialChi_ELF_Style();
			if(chi < 1.0) chi = 1.0;
			B.x = PLAY_MU_0 * (1.0 + chi) * H_target.x;
			B.y = PLAY_MU_0 * (1.0 + chi) * H_target.y;
			B.z = PLAY_MU_0 * (1.0 + chi) * H_target.z;
			continue;
		}

		// Trust region: limit step magnitude
		double dB_mag = sqrt(dB.x*dB.x + dB.y*dB.y + dB.z*dB.z);
		if(dB_mag > max_dB_step)
		{
			double sc = max_dB_step / dB_mag;
			dB.x *= sc; dB.y *= sc; dB.z *= sc;
		}

		// Armijo backtracking with B-magnitude clamping
		double tau = 1.0;
		bool accepted = false;
		for(int ls = 0; ls < 20; ls++)
		{
			TVector3d B_trial(B.x + tau*dB.x, B.y + tau*dB.y, B.z + tau*dB.z);

			// Clamp trial B magnitude to B_max
			double B_trial_mag = sqrt(B_trial.x*B_trial.x + B_trial.y*B_trial.y + B_trial.z*B_trial.z);
			if(B_trial_mag > B_max)
			{
				double sc = B_max / B_trial_mag;
				B_trial.x *= sc; B_trial.y *= sc; B_trial.z *= sc;
			}

			TVector3d H_trial = Forward(B_trial);
			TVector3d res_trial(H_trial.x - H_target.x, H_trial.y - H_target.y, H_trial.z - H_target.z);
			double res_trial_norm = sqrt(res_trial.x*res_trial.x + res_trial.y*res_trial.y + res_trial.z*res_trial.z);
			if(res_trial_norm < res_norm)
			{
				B = B_trial;
				accepted = true;
				break;
			}
			tau *= 0.5;
		}
		if(!accepted)
		{
			// All backtracking failed: take tiny step
			B.x += 0.001 * dB.x;
			B.y += 0.001 * dB.y;
			B.z += 0.001 * dB.z;
		}
	}

	// If Newton didn't converge well, use best B found
	{
		TVector3d H_final = Forward(B);
		TVector3d res_final(H_final.x - H_target.x, H_final.y - H_target.y, H_final.z - H_target.z);
		double res_final_norm = sqrt(res_final.x*res_final.x + res_final.y*res_final.y + res_final.z*res_final.z);
		if(res_final_norm > best_res * 1.01)
		{
			Forward(B_best);  // restore best state
		}
	}

	// Auto-commit state: advance play operators for next step
	for(int k = 0; k < m_K; k++)
		m_pk_prev[k] = m_pk_current[k];

	return m_last_B;
}

//-------------------------------------------------------------------------
// ComputeJacobian: returns dB/dH = inv(dH/dB), compatible with solver
//-------------------------------------------------------------------------
void radTPlayHysteresisMaterial::ComputeJacobian(TMatrix3d& dBdH, double& chi_d) const
{
	if(!m_has_result)
	{
		dBdH = PLAY_MU_0 * Eye();
		chi_d = 0.0;
		return;
	}
	dBdH = Matrix3d_inv(m_last_dHdB);
	double trace = dBdH.Str0.x + dBdH.Str1.y + dBdH.Str2.z;
	chi_d = trace / (3.0 * PLAY_MU_0) - 1.0;
}

//-------------------------------------------------------------------------
// Material interface methods (same pattern as energy model)
//-------------------------------------------------------------------------
TVector3d radTPlayHysteresisMaterial::M(const TVector3d& H)
{
	TVector3d B = Inverse(H);
	double inv_mu0 = 1.0 / PLAY_MU_0;
	return TVector3d(B.x * inv_mu0 - H.x, B.y * inv_mu0 - H.y, B.z * inv_mu0 - H.z);
}

void radTPlayHysteresisMaterial::DefineInstantKsiTensor(
	const TVector3d& H, TMatrix3d& KsiTensor, TVector3d& Mr)
{
	TVector3d Mv = M(H);
	Mr = Mv;

	TMatrix3d dBdH;
	double chi_d;
	ComputeJacobian(dBdH, chi_d);

	// KsiTensor = dM/dH = dB/dH / mu_0 - I
	double inv_mu0 = 1.0 / PLAY_MU_0;
	KsiTensor = inv_mu0 * dBdH - Eye();
}

double radTPlayHysteresisMaterial::ComputeChiFromH(const TVector3d& H)
{
	TVector3d B = Inverse(H);
	double H_mag = sqrt(H.x*H.x + H.y*H.y + H.z*H.z);
	if(H_mag < 1e-30) return GetInitialChi_ELF_Style();
	double B_mag = sqrt(B.x*B.x + B.y*B.y + B.z*B.z);
	return B_mag / (PLAY_MU_0 * H_mag) - 1.0;
}

double radTPlayHysteresisMaterial::ComputeChiDualMethod(double H_mag, double mu_old, double relax)
{
	TVector3d Hv(H_mag, 0, 0);
	return ComputeChiFromH(Hv);
}

double radTPlayHysteresisMaterial::ComputeDifferentialChi(double H_mag)
{
	// Save state so the synthetic 1D Inverse doesn't clobber real 3D state
	std::vector<TVector3d> saved_prev = m_pk_prev;
	std::vector<TVector3d> saved_pinning = m_pk_pinning;
	std::vector<TVector3d> saved_current = m_pk_current;
	TVector3d saved_H = m_last_H, saved_B = m_last_B;
	TMatrix3d saved_dHdB = m_last_dHdB;
	double saved_chi = m_last_chi, saved_chi_d = m_last_chi_d;
	bool saved_has = m_has_result;

	TVector3d Hv(H_mag, 0, 0);
	Inverse(Hv);
	TMatrix3d dBdH;
	double chi_d;
	ComputeJacobian(dBdH, chi_d);

	// Restore state
	m_pk_prev = saved_prev;
	m_pk_pinning = saved_pinning;
	m_pk_current = saved_current;
	m_last_H = saved_H; m_last_B = saved_B;
	m_last_dHdB = saved_dHdB;
	m_last_chi = saved_chi; m_last_chi_d = saved_chi_d;
	m_has_result = saved_has;
	return chi_d;
}

//-------------------------------------------------------------------------
// Serialization
//-------------------------------------------------------------------------
void radTPlayHysteresisMaterial::DumpBin(CAuxBinStrVect& oStr,
	std::vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers,
	int& gUniqueMapKey, int elemKey)
{
	vElemKeysOut.push_back(elemKey);
	int matType = Type_Material();
	oStr << matType;
	oStr << m_K;
	for(int k = 0; k < m_K; k++) oStr << m_eta[k];
	for(int k = 0; k < m_K; k++)
	{
		oStr << m_tables[k].n;
		for(int i = 0; i < m_tables[k].n; i++) oStr << m_tables[k].r[i];
		for(int i = 0; i < m_tables[k].n; i++) oStr << m_tables[k].f[i];
	}
	// Save state
	for(int k = 0; k < m_K; k++)
	{
		oStr << m_pk_prev[k].x << m_pk_prev[k].y << m_pk_prev[k].z;
	}
}

radTPlayHysteresisMaterial::radTPlayHysteresisMaterial(CAuxBinStrVect& inStr)
	: m_B_mono_max(1.0), m_H_mono_max(0),
	  m_last_chi(0), m_last_chi_d(0), m_has_result(false)
{
	inStr >> m_K;
	m_eta.resize(m_K);
	for(int k = 0; k < m_K; k++) inStr >> m_eta[k];
	m_tables.resize(m_K);
	for(int k = 0; k < m_K; k++)
	{
		int n;
		inStr >> n;
		m_tables[k].n = n;
		m_tables[k].r.resize(n);
		m_tables[k].f.resize(n);
		for(int i = 0; i < n; i++) inStr >> m_tables[k].r[i];
		for(int i = 0; i < n; i++) inStr >> m_tables[k].f[i];
		m_tables[k].r_max = (n > 0) ? m_tables[k].r[n-1] : 0.0;
		PrecomputePlayTable(m_tables[k]);
	}
	TVector3d zero(0, 0, 0);
	m_pk_prev.resize(m_K, zero);
	m_pk_pinning.resize(m_K, zero);
	m_pk_current.resize(m_K, zero);
	for(int k = 0; k < m_K; k++)
	{
		inStr >> m_pk_prev[k].x >> m_pk_prev[k].y >> m_pk_prev[k].z;
		m_pk_pinning[k] = m_pk_prev[k];
	}
	m_last_H = zero;
	m_last_B = zero;
	m_last_dHdB = Eye();

	ComputeMonotoneLimits();
}

//-------------------------------------------------------------------------
// Application method: create and register play hysteresis material
//-------------------------------------------------------------------------
int radTApplication::SetPlayHysteresisMaterial(
	int K, const double* eta,
	const std::vector<std::vector<double>>& r_tables,
	const std::vector<std::vector<double>>& f_tables)
{
	radTPlayHysteresisMaterial* MaterPtr = nullptr;
	try
	{
		MaterPtr = new radTPlayHysteresisMaterial(K, eta, r_tables, f_tables);
		if(MaterPtr == nullptr) { Send.ErrorMessage("Radia::Error900"); return 0; }

		radThg hg(MaterPtr);
		MaterPtr = nullptr;
		int ElemKey = AddElementToContainer(hg);
		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		if(MaterPtr) delete MaterPtr;
		Initialize();
		return 0;
	}
}

//-------------------------------------------------------------------------
