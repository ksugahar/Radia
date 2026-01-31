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

#ifdef _OPENMP
#include <omp.h>
#endif


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
		double Hprev = -1, Mprev = -1;
		for(int i=0; i<LenArrayHB; i++)
		{
			if(tArrayHB->x < Hprev) { Send.ErrorMessage("Radia::Error071"); return 0;}
			if(tArrayHB->y < 0.95*Mprev) { Send.WarningMessage("Radia::Warning014");}

			Hprev = tArrayHB->x; Mprev = tArrayHB->y;
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

int radTApplication::ApplyDrawAttrToElem_g3d(int ElemKey, double* RGB_col, long lenRGB_col, double InLineThickness)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}
		TVector3d RGB_colVect;
		if(!ValidateVector3d(RGB_col, lenRGB_col, &RGB_colVect)) return 0;
		// May be not necessary?
		if((RGB_col[0]<0.) || (RGB_col[1]<0.) || (RGB_col[2]<0.) || 
		   (RGB_col[0]>1.) || (RGB_col[1]>1.) || (RGB_col[2]>1.))
		{
			Send.ErrorMessage("Radia::Error008"); return 0;
		}

		radTMapOfDrawAttr::iterator iter = MapOfDrawAttr.find(ElemKey);
		if(!(iter == MapOfDrawAttr.end())) MapOfDrawAttr.erase(iter);

		radRGB ColRGB(RGB_colVect.x, RGB_colVect.y, RGB_colVect.z);
		radTDrawAttr DrawAttr;
		DrawAttr.RGB_col = ColRGB;
		DrawAttr.LineThickness = (InLineThickness<0)? 0.001 : InLineThickness;

		MapOfDrawAttr[ElemKey] = DrawAttr;

		if(SendingIsRequired) Send.Int(ElemKey);
		return ElemKey;
	}
	catch(...)
	{
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::RemoveDrawAttrFromElem_g3d(int ElemKey)
{
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		radTMapOfDrawAttr::iterator iter = MapOfDrawAttr.find(ElemKey);
		if(iter == MapOfDrawAttr.end()) 
		{
			Send.ErrorMessage("Radia::Error013");
			return 0;
		}
		else
		{
			MapOfDrawAttr.erase(iter);
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

int radTApplication::GraphicsForElem_g3d(int ElemKey, int InShowSymmetryChilds, const char** arOptionNames, const char** arOptionValues, int numOptions)
{
	radTg3dGraphPresent* g3dGraphPresentPtr = nullptr;
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		radTOptionNames OptNam;
		const char* OptNamesToFind[] = {OptNam.Debug};
		char OptValsFoundParsed[] = {0};
		char &doDebug = OptValsFoundParsed[0]; // 0- No; 1- Yes;
		if(!OptNam.findParseOptionValues(arOptionNames, arOptionValues, numOptions, OptNamesToFind, 1, OptValsFoundParsed, 0, 0))
		{
			Send.ErrorMessage("Radia::Error062"); return 0;
		}

		radGraphPresOptions InGraphPresOptions((char)InShowSymmetryChilds, doDebug);

		Send.GenInitDraw();

		g3dGraphPresentPtr = g3dPtr->CreateGraphPresent();

		g3dGraphPresentPtr->SetGraphPresOptions(InGraphPresOptions);
		g3dGraphPresentPtr->MapOfDrawAttrPtr = &MapOfDrawAttr;
		g3dGraphPresentPtr->RetrieveDrawAttr(ElemKey);

		g3dGraphPresentPtr->GenDraw();

		delete g3dGraphPresentPtr;
		g3dGraphPresentPtr = nullptr;
		return ElemKey;
	}
	catch(...)
	{
		if(g3dGraphPresentPtr) delete g3dGraphPresentPtr;
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

int radTApplication::GraphicsForElem_g3d_VTK(int ElemKey, const char** OptionNames, const char** OptionValues, int OptionCount) //OC04112019 (from R. Nagler's radTApplication::GoObjGeometry)
{
	radTg3dGraphPresent* g3dGraphPresentPtr = nullptr;
	try
	{
		radThg hg;
		if(!ValidateElemKey(ElemKey, hg)) return 0;
		radTg3d* g3dPtr = Cast.g3dCast(hg.rep); if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); return 0;}

		//bool SendingWasAlreadyDone = false;

		char OptBits[4];
		char& DoShowLines = OptBits[0];
		char& DoShowFaces = OptBits[1];
		char& DoShowFrameAxes = OptBits[2];
		char& DoShowSymChilds = OptBits[3];
		if(!DecodeViewingOptions(OptionNames, OptionValues, OptionCount, OptBits)) return 0;

		radGraphPresOptions InGraphPresOptions(DoShowSymChilds);
		g3dGraphPresentPtr = g3dPtr->CreateGraphPresent();

		char DrawFacilityInd = 2; // VTK export facility index
		g3dGraphPresentPtr->DrawFacilityInd = DrawFacilityInd;

		radTg3dGraphPresent::Send = Send;

		g3dGraphPresentPtr->SetGraphPresOptionsExt(InGraphPresOptions, DoShowLines, DoShowFaces);
		g3dGraphPresentPtr->MapOfDrawAttrPtr = &MapOfDrawAttr;
		g3dGraphPresentPtr->RetrieveDrawAttr(ElemKey);

		g3dGraphPresentPtr->GenDraw();
		if(DoShowFrameAxes) g3dGraphPresentPtr->DrawFrameLines();

		int keyGeomData = (radTg3dGraphPresent::Send).GeomDataToBuffer();

		delete g3dGraphPresentPtr;
		g3dGraphPresentPtr = nullptr;
		return keyGeomData;
	}
	catch(...)
	{
		if(g3dGraphPresentPtr) delete g3dGraphPresentPtr;
		Initialize(); return 0;
	}
}

//-------------------------------------------------------------------------

void radTApplication::GraphicsForAll_g3d(int InShowSymmetryChilds)
{
	radTg3dGraphPresent* g3dGraphPresentPtr = nullptr;
	try
	{
		int TotalElem = (int)(GlobalMapOfHandlers.size());
		// RAII: Use std::vector for automatic cleanup
		std::vector<radTg3d*> vG3dPtrPtr(TotalElem);
		std::vector<int> vKeyPtr(TotalElem);
		radTg3d** g3dPtrPtr = vG3dPtrPtr.data();
		int* KeyPtr = vKeyPtr.data();

		radGraphPresOptions InGraphPresOptions((char)InShowSymmetryChilds);

		int g3dPresElemCount = 0;
		for(radTmhg::const_iterator iter = GlobalMapOfHandlers.begin();
			iter != GlobalMapOfHandlers.end(); ++iter)
		{
			radTg* gPtr = ((*iter).second).rep;
			radTg3d g3d;
			if(gPtr->Type_g() == g3d.Type_g())
				if(!static_cast<radTg3d*>(gPtr)->IsGroupMember)
				{
					g3dPtrPtr[g3dPresElemCount] = static_cast<radTg3d*>(gPtr);
					KeyPtr[g3dPresElemCount++] = (*iter).first;
				}
		}
		if(g3dPresElemCount != 0)
		{
			Send.GenInitDraw();

			Send.InitOutList(g3dPresElemCount);
			for(int i = 0; i < g3dPresElemCount; i++)
			{
				g3dGraphPresentPtr = g3dPtrPtr[i]->CreateGraphPresent();

				g3dGraphPresentPtr->SetGraphPresOptions(InGraphPresOptions);
				g3dGraphPresentPtr->MapOfDrawAttrPtr = &MapOfDrawAttr;
				g3dGraphPresentPtr->RetrieveDrawAttr(KeyPtr[i]);
				g3dGraphPresentPtr->GenDraw();
				delete g3dGraphPresentPtr;
				g3dGraphPresentPtr = nullptr;
			}
		}
		else Send.ErrorMessage("Radia::Error101");
		// RAII: automatic cleanup (also fixes missing delete[] KeyPtr!)
	}
	catch(...)
	{
		if(g3dGraphPresentPtr) delete g3dGraphPresentPtr;
		Initialize(); return;
	}
}

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
				// Copy matrix data (column-major format)
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

		// Valid methods: LU (0), BICGSTAB (1)
		// Note: BICGSTAB_HMATRIX (2) requires RADIA_USE_HACAPK to be defined
#ifdef RADIA_USE_HACAPK
		if(MethNo != RadSolverMethod::LU && MethNo != RadSolverMethod::BICGSTAB && MethNo != RadSolverMethod::BICGSTAB_HMATRIX) { Send.ErrorMessage("Radia::Error028"); return 0;}
#else
		if(MethNo != RadSolverMethod::LU && MethNo != RadSolverMethod::BICGSTAB) { Send.ErrorMessage("Radia::Error028"); return 0;}
#endif
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
				hacapk_params.print_level = 1;  // Standard output
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

			// Valid methods for AutoRelax: LU (0), BICGSTAB (1)
			// Note: BICGSTAB_HMATRIX (2) requires RADIA_USE_HACAPK to be defined
#ifdef RADIA_USE_HACAPK
			if(MethNo != RadSolverMethod::LU && MethNo != RadSolverMethod::BICGSTAB && MethNo != RadSolverMethod::BICGSTAB_HMATRIX) { Send.ErrorMessage("Radia::Error041"); return 0; }
#else
			if(MethNo != RadSolverMethod::LU && MethNo != RadSolverMethod::BICGSTAB) { Send.ErrorMessage("Radia::Error041"); return 0; }
#endif

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
				hacapk_params.print_level = 1;  // Standard output
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

		// For HACApK solver (method 2), skip dense matrix construction
		// HACApK builds its own H-matrix, dense matrix is unnecessary overhead
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
bool radTApplication::ApplyIMASymmetryToInteraction(radTInteraction* pIntrc, const char* imageSpec)
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

	// Build IMA matrix
	if(pIntrc->IsIMAEnabled())
	{
		pIntrc->SetupInteractMatrix_IMA();
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
				// Reference: ELF_MAGIC implementation (src/dll/m_fmm3d.f90)
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
	Send.InitOutList(AmOfPoints, 0);

	for(int i=0; i<AmOfPoints; i++)
	{
		Send.InitOutList(2, 0);
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

	// OpenMP status
#ifdef _OPENMP
	dOut[4] = 1.0;  // OpenMP enabled
	dOut[5] = (double)omp_get_max_threads();
#else
	dOut[4] = 0.0;  // OpenMP disabled
	dOut[5] = 1.0;
#endif

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

		// Convert input points from user units to internal SI units (meters)
		// Radia v1.4.3+: Internal unit system is SI (meters), matching ELF
		// m_lengthUnitScale is 1.0 for m (default), 0.001 for mm
		double scale = m_lengthUnitScale;

		// Create a copy of points in meters (internal SI units)
		std::vector<double> points_internal(n_points * 3);
		for(int i = 0; i < n_points * 3; ++i)
		{
			points_internal[i] = points[i] * scale;
		}

		// Use the solid angle implementation from rad_point_classify
		RadPointClassify::ClassifyPointsFromHandle(
			n_points,
			points_internal.data(),
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
                                        double* points, int container_handle, int method)
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

		// Initialize output arrays to zero (if provided)
		if(B_out) std::memset(B_out, 0, n_points * 3 * sizeof(double));
		if(H_out) std::memset(H_out, 0, n_points * 3 * sizeof(double));

		// Magnetic constant (permeability of free space)
		// B field is stored internally in A/m (same as H), needs conversion to Tesla
		const double Mu0 = 4. * 3.1415926535897932 * 1.e-7; // T*m/A

		// Apply unit scaling: convert user units to internal mm
		// Same as rad_c_interface.cpp Field() function
		double scale = GetLengthUnitScale();

		// Create field key for B and H computation
		radTFieldKey FieldKey;
		FieldKey.B_ = true;
		FieldKey.H_ = true;

		TVector3d ZeroVect(0., 0., 0.);

		// Compute field at each point using OpenMP parallelization
		// method: 0 = direct, 1 = FMM (not yet implemented)
		//
		// OpenMP parallelization is safe here because:
		// 1. Each iteration creates its own thread-local radTField object
		// 2. Each iteration writes to different output array indices
		// 3. B_genComp() only reads from the g3dPtr object (no writes)
		// 4. The FieldKey and ZeroVect are copied by value into each thread's radTField
		//
		// Note: Exception handling inside parallel region is tricky - we avoid throwing
		// exceptions inside B_genComp by ensuring g3dPtr is valid before the loop.
		#ifdef _OPENMP
		#pragma omp parallel for schedule(static) if(n_points > 100)
		#endif
		for(int i = 0; i < n_points; ++i)
		{
			TVector3d pt;
			pt.x = points[i * 3 + 0] * scale;
			pt.y = points[i * 3 + 1] * scale;
			pt.z = points[i * 3 + 2] * scale;

			// Create thread-local field object for each point
			// All members are either values (not pointers/references) or initialized here
			radTFieldKey localFieldKey;
			localFieldKey.B_ = true;
			localFieldKey.H_ = true;
			TVector3d localZeroVect(0., 0., 0.);

			// Constructor: (FieldKey, CompCriterium, P, B, H, A, M, J=0. -> TVector3d(0,0,0))
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
		}
	}
	catch(...)
	{
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

		// Initialize output array to zero
		std::memset(phi_out, 0, n_points * sizeof(double));

		// Apply unit scaling: convert user units to internal mm
		double scale = GetLengthUnitScale();

		// Create field key for Phi computation
		radTFieldKey FieldKey;
		FieldKey.Phi_ = true;

		TVector3d ZeroVect(0., 0., 0.);

		// Compute scalar potential at each point
		// OpenMP parallelization disabled due to Intel OpenMP deadlock issues
		for(int i = 0; i < n_points; ++i)
		{
			TVector3d pt;
			pt.x = points[i * 3 + 0] * scale;
			pt.y = points[i * 3 + 1] * scale;
			pt.z = points[i * 3 + 2] * scale;

			// Create thread-local field key and vectors
			radTFieldKey localFieldKey;
			localFieldKey.Phi_ = true;
			TVector3d localZeroVect(0., 0., 0.);

			radTField Field(localFieldKey, CompCriterium, pt, localZeroVect, localZeroVect, localZeroVect, localZeroVect, 0.);
			g3dPtr->B_genComp(&Field);

			// Phi is dimensionless (in internal units), convert based on length unit
			// phi_m = (1/4pi) * integral(M . r / r^3) dV
			// Units: [phi_m] = A (magnetic scalar potential)
			phi_out[i] = Field.Phi;
		}
	}
	catch(...)
	{
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

		// Apply unit scaling: convert user units to internal mm
		double scale = GetLengthUnitScale();

		// Create field key for A computation
		radTFieldKey FieldKey;
		FieldKey.A_ = true;

		TVector3d ZeroVect(0., 0., 0.);

		// Compute vector potential at each point
		// OpenMP parallelization disabled due to Intel OpenMP deadlock issues
		for(int i = 0; i < n_points; ++i)
		{
			TVector3d pt;
			pt.x = points[i * 3 + 0] * scale;
			pt.y = points[i * 3 + 1] * scale;
			pt.z = points[i * 3 + 2] * scale;

			// Create thread-local field key and vectors
			radTFieldKey localFieldKey;
			localFieldKey.A_ = true;
			TVector3d localZeroVect(0., 0., 0.);

			radTField Field(localFieldKey, CompCriterium, pt, localZeroVect, localZeroVect, localZeroVect, localZeroVect, 0.);
			g3dPtr->B_genComp(&Field);

			// A has units of T*m (Tesla-meter)
			// Radia internal A is stored in consistent units
			A_out[i * 3 + 0] = Field.A.x;
			A_out[i * 3 + 1] = Field.A.y;
			A_out[i * 3 + 2] = Field.A.z;
		}
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

//-------------------------------------------------------------------------
