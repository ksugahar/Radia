/*-------------------------------------------------------------------------
*
* File name:      radmater.h
*
* Project:        RADIA
*
* Description:    Material relations and auxiliary functions for relaxation
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#ifndef __RADMATER_H
#define __RADMATER_H

#include "rad_serialization.h"
#include "gmvect.h"
#include "rad_geometry_base.h"

#include <math.h>
#include <sstream>

//-------------------------------------------------------------------------

class radTMaterial;

//-------------------------------------------------------------------------

class radTg3dRelax;

//-------------------------------------------------------------------------

class radTMaterial : public radTg {
public:
	TVector3d RemMagn; // Don't make it private nor protected
	char EasyAxisDefined;

	radTMaterial(const TVector3d& InRemMagn, char InEasyAxisDefined)
	{
		RemMagn = InRemMagn; EasyAxisDefined = InEasyAxisDefined;
	}
	radTMaterial() { EasyAxisDefined = 0;}

	int Type_g() { return 3;}
	virtual int Type_Material() { return 0;}

	virtual TVector3d M(const TVector3d& H) { return 0.*H;}
	virtual void DefineInstantKsiTensor(const TVector3d&, TMatrix3d&, TVector3d&) {}
	virtual void MultMatrByInstKsiAndMr(const TVector3d&, const TMatrix3d&, TMatrix3d&, TVector3d&) {}

	virtual void FindNewH(TVector3d&, const TMatrix3d&, const TVector3d&, double) {}
	//virtual void FindNewH(TVector3d&, const TMatrix3d&, const TVector3d&, double, radTg3dRelax*) {}
	//virtual void FindNewH(TVector3d&, const TMatrix3d&, const TVector3d&, double, radTg3dRelax*, void* p=0) {} //OC140103

	virtual int FinishSetup(TVector3d&) { return 1;}

	int FinishDuplication(radTMaterial* MatPtr, radThg& hg)
	{
		radTSend Send;
		if(MatPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		radThg hgLoc(MatPtr); hg = hgLoc; return 1;
	}

	// Dump / DumpBin_Material / DumpBinParse_Material REMOVED (Phase B2b/B2c, 2026-04-15)

	//static void SteerNewH(TVector3d& PrevH, TVector3d& InstantH, void* pvAuxRelax); //OC140103
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTLinearAnisotropMaterial : public radTMaterial {
	double KsiPar, KsiPerp;
	TMatrix3d KsiTensor;

public:
	radTLinearAnisotropMaterial(const double* InKsiArray, const TVector3d& InRemMagn, char InEasyAxisDefined)
		: radTMaterial(InRemMagn, InEasyAxisDefined)
	{
		KsiPar = InKsiArray[0]; KsiPerp = InKsiArray[1];
		if(InEasyAxisDefined) SetupKsiTensor();
	}

	radTLinearAnisotropMaterial(double InKsiPar, double InKsiPerp, double InMr)
	{
		KsiPar = InKsiPar; KsiPerp = InKsiPerp;
		RemMagn.x = InMr; RemMagn.y = RemMagn.z = 0.;
	}

/*
	radTLinearAnisotropMaterial(int IndDB, const TVector3d& InRemMagn, char InEasyAxisDefined)
		: radTMaterial(InRemMagn, InEasyAxisDefined)
	{
		radTLinearAnisotropMaterialDB *pMat = (radTLinearAnisotropMaterialDB*)(MaterDB[IndDB].rep);

		KsiPar = pMat->KsiPar; KsiPerp = pMat->KsiPerp;
		if(InEasyAxisDefined) SetupKsiTensor();
		else
		{
			RemMagn.x = pMat->Mr;
			RemMagn.y = RemMagn.z = 0.;
		}
	}
*/

	// radTLinearAnisotropMaterial(CAuxBinStrVect&) REMOVED (Phase B2c, 2026-04-15)

	radTLinearAnisotropMaterial() {}

	int Type_Material() { return 1;}
	
	inline void SetupKsiTensor();
	TVector3d M(const TVector3d& H) { return KsiTensor*H;}  // Pure linear response (no remanence)
	void DefineInstantKsiTensor(const TVector3d& InstantH, TMatrix3d& InstantKsiTensor, TVector3d& InstantMr)
	{
		InstantKsiTensor = KsiTensor; InstantMr = TVector3d(0,0,0);  // No remanence for linear material
	}
	void MultMatrByInstKsiAndMr(const TVector3d& InstantH, const TMatrix3d& Matr, TMatrix3d& MultByKsi, TVector3d& MultByMr)
	{
		MultByKsi = Matr * KsiTensor; MultByMr = TVector3d(0,0,0);  // No remanence for linear material
	}
	//void FindNewH(TVector3d& H, const TMatrix3d& Matr, const TVector3d& H_Ext, double DesiredPrecOnMagnetizE2, radTg3dRelax* pMag, void* p=0) //OC140103
	void FindNewH(TVector3d& H, const TMatrix3d& Matr, const TVector3d& H_Ext, double DesiredPrecOnMagnetizE2)
	{
		TVector3d ESt1(1.,0.,0.), ESt2(0.,1.,0.), ESt3(0.,0.,1.);
		TMatrix3d E(ESt1, ESt2, ESt3);
		TMatrix3d BufMatr = E - Matr*KsiTensor;
		TMatrix3d InvBufMatr;
		Matrix3d_inv(BufMatr, InvBufMatr);
		H = InvBufMatr*H_Ext;  // No remanence term for linear material
	}

	int FinishSetup(TVector3d& Magn)
	{
		if(!EasyAxisDefined)
		{
			double AbsLocMagn = sqrt(Magn.x*Magn.x + Magn.y*Magn.y + Magn.z*Magn.z);
		
			radTSend Send;
			const double AbsTol = 1.E-10;
			if(AbsLocMagn < AbsTol) { Send.ErrorMessage("Radia::Error107"); return 0;}

			// For pure linear material (no remanence): RemMagn should stay zero
			// Only rescale if RemMagn was explicitly set to non-zero value
			double AbsRemMagn = sqrt(RemMagn.x*RemMagn.x + RemMagn.y*RemMagn.y + RemMagn.z*RemMagn.z);
			if(AbsRemMagn > AbsTol) {
				// RemMagn is non-zero - rescale to align with Magn direction
				RemMagn = (RemMagn.x/AbsLocMagn)*Magn;
			}
			// If RemMagn is zero, leave it zero (pure linear material)
			SetupKsiTensor();
			EasyAxisDefined = 1;
		}
		return 1;
	}

	int DuplicateItself(radThg& hg, radTApplication*, char) 
	{
		return FinishDuplication(new radTLinearAnisotropMaterial(*this), hg);
	}

	// Dump / DumpBin REMOVED (Phase B2b/B2c, 2026-04-15)

	int SizeOfThis() { return sizeof(radTLinearAnisotropMaterial);}
};

//-------------------------------------------------------------------------

inline void radTLinearAnisotropMaterial::SetupKsiTensor()
{
	double AbsRemMagn = sqrt(RemMagn.x*RemMagn.x + RemMagn.y*RemMagn.y + RemMagn.z*RemMagn.z);
	TVector3d L = (1./AbsRemMagn)*RemMagn;
	double DeltaKsi = KsiPar-KsiPerp;
	double LxLx, LyLy, LzLz;
	LxLx=L.x*L.x; LyLy=L.y*L.y; LzLz=L.z*L.z;
	TVector3d Str0(KsiPar*LxLx+KsiPerp*(LyLy+LzLz), DeltaKsi*L.x*L.y, DeltaKsi*L.x*L.z);
	TVector3d Str1(Str0.y, KsiPar*LyLy+KsiPerp*(LxLx+LzLz), DeltaKsi*L.y*L.z);
	TVector3d Str2(Str0.z, Str1.z, KsiPar*LzLz+KsiPerp*(LxLx+LyLy));
	KsiTensor.Str0 = Str0; KsiTensor.Str1 = Str1; KsiTensor.Str2 = Str2;
}

//-------------------------------------------------------------------------

// radTLinearAnisotropMaterial::Dump REMOVED (Phase B2b, 2026-04-15)

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTLinearIsotropMaterial : public radTMaterial {
	double Ksi;

public:
	radTLinearIsotropMaterial(double InKsi) { Ksi = InKsi;}
	radTLinearIsotropMaterial(const double* InKsiArray, const TVector3d& InRemMagn, char InEasyAxisDefined) 
		: radTMaterial(InRemMagn, InEasyAxisDefined) { Ksi = InKsiArray[0];}
	
	// radTLinearIsotropMaterial(CAuxBinStrVect&) REMOVED (Phase B2c, 2026-04-15)

	radTLinearIsotropMaterial() {}

	int Type_Material() { return 2;}

	TVector3d M(const TVector3d& H) { return Ksi * H;} // Pure linear response (no remanence)
	void DefineInstantKsiTensor(const TVector3d&, TMatrix3d&, TVector3d&);
	void MultMatrByInstKsiAndMr(const TVector3d&, const TMatrix3d& Matr, TMatrix3d& MultByKsi, TVector3d& MultByMr)
	{
		MultByKsi = Ksi * Matr; MultByMr = TVector3d(0,0,0);  // No remanence for linear material
	}
	//void FindNewH(TVector3d& H, const TMatrix3d& Matr, const TVector3d& H_Ext, double DesiredPrecOnMagnetizE2, radTg3dRelax* pMag, void* p=0) //OC140103
	void FindNewH(TVector3d& H, const TMatrix3d& Matr, const TVector3d& H_Ext, double DesiredPrecOnMagnetizE2) //OC140103
	{
		TVector3d ESt1(1.,0.,0.), ESt2(0.,1.,0.), ESt3(0.,0.,1.);
		TMatrix3d E(ESt1, ESt2, ESt3);
		TMatrix3d BufMatr = E - Ksi*Matr;
		TMatrix3d InvBufMatr;
		Matrix3d_inv(BufMatr, InvBufMatr);
		H = InvBufMatr*H_Ext;  // No remanence term for linear material
	}

	int DuplicateItself(radThg& hg, radTApplication*, char) 
	{
		return FinishDuplication(new radTLinearIsotropMaterial(*this), hg);
	}

	// Dump / DumpBin REMOVED (Phase B2b/B2c, 2026-04-15)

	int SizeOfThis() { return sizeof(radTLinearIsotropMaterial);}
};

//-------------------------------------------------------------------------

inline void radTLinearIsotropMaterial::DefineInstantKsiTensor(const TVector3d& InstantH, TMatrix3d& InstantKsiTensor, TVector3d& InstantMr)
{
	TVector3d E_Str0(1.,0.,0.), E_Str1(0.,1.,0.), E_Str2(0.,0.,1.);
	TMatrix3d E(E_Str0, E_Str1, E_Str2);
	InstantKsiTensor = Ksi*E; InstantMr = TVector3d(0,0,0);  // No remanence for linear material
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTNonlinearIsotropMaterial : public radTMaterial {
	double Ms[3], ks[3];
	int lenMs_ks;

	// B-H curve storage (ELF-compatible approach)
	// gArrayHB stores (H, B) pairs directly, NOT (H, M)
	std::vector<TVector2d> vgArrayHB;
	TVector2d* gArrayHB;          // gArrayHB[i].x = H, gArrayHB[i].y = B (not M!)
	std::vector<double> vgdBdH;   // derivatives dB/dH
	double* gdBdH;
	int gLenArrayHB;

	double gMaxKsi;

public:
	radTNonlinearIsotropMaterial(const double* InMsArray, const double* In_ksArray, int In_lenMs_ks)	
	{ 
		gMaxKsi = 0;
		lenMs_ks = In_lenMs_ks;
		for(int i=0; i<lenMs_ks; i++) { Ms[i]=InMsArray[i]; ks[i]=In_ksArray[i]; gMaxKsi+=ks[i];}

		gArrayHB = 0; gLenArrayHB = 0; gdBdH = 0;
	}
	radTNonlinearIsotropMaterial(TVector2d* InArrayHB, int InLenArrayHB)
	{
		gArrayHB = 0; gdBdH = 0; gLenArrayHB = 0;
		double ZeroTol = 1e-10;
		char PrependZero = 0;
		if((InArrayHB->x > ZeroTol) && (InArrayHB->y > ZeroTol))
		{
			InLenArrayHB++; PrependZero = 1;
		}

		gLenArrayHB = InLenArrayHB;
		AllocateArrays(InLenArrayHB);
		CopyArrayHB(gArrayHB, InArrayHB, InLenArrayHB, PrependZero);
		// Compute dB/dH derivatives for B-H curve
		Compute_dBdH(gArrayHB, gdBdH, gLenArrayHB, gMaxKsi);
	}
	// radTNonlinearIsotropMaterial(CAuxBinStrVect&) REMOVED (Phase B2c, 2026-04-15)

	radTNonlinearIsotropMaterial() { gArrayHB = 0; gLenArrayHB = 0; gdBdH = 0; gMaxKsi = 0;}
	~radTNonlinearIsotropMaterial() { DeallocateArrays();}

	int Type_Material() { return 3;}

	TVector3d M(const TVector3d& H);
	void DefineInstantKsiTensor(const TVector3d&, TMatrix3d&, TVector3d&);
	void MultMatrByInstKsiAndMr(const TVector3d&, const TMatrix3d&, TMatrix3d&, TVector3d&);
	//void FindNewH(TVector3d&, const TMatrix3d&, const TVector3d&, double, radTg3dRelax*, void*); //OC140103
	void FindNewH(TVector3d&, const TMatrix3d&, const TVector3d&, double); //OC140103

	// Dump / DumpBin REMOVED (Phase B2b/B2c, 2026-04-15)

	// Compute dB/dH derivatives for B-H curve (renamed from Compute_dMdH)
	static void Compute_dBdH(TVector2d* ArrayHB, double* dBdH, int LenArrayHB, double& MaxKsi);

	// Check and correct dB/dH derivatives (renamed from CheckAndCorrect_dMdH)
	static void CheckAndCorrect_dBdH(TVector2d* ArrayHB, double* dBdH, int LenArrayHB);

	static double Derivative5(TVector2d* f, int PoIndx);
	static double Derivative3(TVector2d* f, int PoIndx);

	// ELF-compatible B(H) interpolation: returns M = B/mu_0 - H
	static double AbsMvsAbsH_Interpol(double AbsH, TVector2d* ArrayHB, double* dBdH, int LenArrayHB);

	// ELF-compatible inverse: find H given M
	static double AbsHvsAbsM_Interpol(double AbsM, TVector2d* ArrayHB, double* dBdH, int LenArrayHB);

	// ELF Method 2: H+B sum interpolation (CGS normalized)
	// Given hb_sum = H/H_scale + B/B_scale, find H and B on BH curve
	static void InterpolateBH_HBSum(double hb_sum_target, TVector2d* ArrayHB, int LenArrayHB,
	                                double H_scale, double B_scale, double& H_out, double& B_out);

	// ELF dual-method chi update: returns optimal chi using both methods
	// Method 1: Standard mu = B(H)/(mu_0*H)
	// Method 2: H+B sum interpolation
	// Selects method with smaller |mu_new - mu_old|
	// relax: under-relaxation parameter (0.0 = full step, >0 = under-relaxation)
	double ComputeChiDualMethod(double H_mag, double mu_old, double relax = 0.0) const;

	// Compute differential susceptibility chi_d = (dB/dH)/mu_0 - 1 at given H
	// Used for Newton-Raphson nonlinear iteration (tangent linearization)
	double ComputeDifferentialChi(double H_mag) const;

	// Public wrapper for inverse B-H lookup (used by relaxation methods)
	double GetHfromM(double AbsM) const
	{
		if(gLenArrayHB > 0 && gArrayHB != nullptr && gdBdH != nullptr)
		{
			return AbsHvsAbsM_Interpol(AbsM, gArrayHB, gdBdH, gLenArrayHB);
		}
		// Fallback for analytical formula: use Newton iteration
		return GetHfromM_Analytical(AbsM);
	}

	double GetHfromM_Analytical(double AbsM) const;

	// ELF mucal0 style: get initial chi from 2nd point of B-H curve
	// Returns chi = B2/(mu0*H2) - 1 where (H2, B2) is the 2nd data point
	// Returns -1 if B-H table is not available or has fewer than 2 points
	double GetInitialChi_ELF_Style() const
	{
		if(gLenArrayHB >= 2 && gArrayHB != nullptr)
		{
			double H2 = gArrayHB[1].x;  // H value at 2nd point (index 1)
			double B2 = gArrayHB[1].y;  // B value at 2nd point
			const double MU_0_local = 4.0 * 3.14159265358979323846 * 1.0e-7;
			if(H2 > 1.0e-10)
			{
				// chi = B/(mu0*H) - 1 = mu_r - 1
				return B2 / (MU_0_local * H2) - 1.0;
			}
		}
		return -1.0;  // Indicate failure
	}

	// Get B saturation value from BH curve (ELF-compatible)
	// ELF uses: B_sat = B(last) - H(last)
	// Note: In Radia's BH curve, gArrayHB[i].x = H, gArrayHB[i].y = B
	// ELF subtracts H from B for normalization (though dimensionally unusual)
	double GetBsaturation() const
	{
		if(gLenArrayHB >= 1 && gArrayHB != nullptr)
		{
			// ELF-compatible: B_sat = B(last) - H(last)
			double B_last = gArrayHB[gLenArrayHB - 1].y;
			double H_last = gArrayHB[gLenArrayHB - 1].x;
			double B_sat = B_last - H_last;
			if(B_sat < 1.0e-10) return 1.0;  // Fallback
			return B_sat;
		}
		return 1.0;  // Fallback: 1 Tesla
	}

	// ELF-compatible: returns M and dM/dH from B-H curve
	static void AbsMvsAbsH_FuncAndDer_Interpol(double AbsH, TVector2d* ArrayHB, double* dBdH, int LenArrayHB, double& f, double& fDer);
	void DefineScalarM(double AbsInstantH, double& f, double& InstKsi);
	void DefineScalarM_dMdH(double AbsInstantH, double& f, double& dfdH);

	double FuncNewAbsH(double AbsH, const TMatrix3d& Matr, const TVector3d& H_Ext);
	double FuncToZero(double AbsH, const TMatrix3d& Matr, const TVector3d& H_Ext) 
	{
		return FuncNewAbsH(AbsH, Matr, H_Ext) - AbsH;
	}

	int DuplicateItself(radThg& hg, radTApplication*, char)
	{// Add more if new members!
		radTSend Send;
		radTNonlinearIsotropMaterial* pNewMater = 0;
		if((gArrayHB != 0) && (gdBdH != 0) && (gLenArrayHB != 0))
		{
			pNewMater = new radTNonlinearIsotropMaterial();
			if(pNewMater == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

			if(!pNewMater->AllocateArrays(gLenArrayHB)) { Send.ErrorMessage("Radia::Error900"); return 0;}
			CopyArrayHB(pNewMater->gArrayHB, gArrayHB, gLenArrayHB, 0);
			CopyArray_dBdH(pNewMater->gdBdH, gdBdH, gLenArrayHB);
		}
		else pNewMater = new radTNonlinearIsotropMaterial(*this);
		if(pNewMater == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}

		return FinishDuplication(pNewMater, hg);
	}

	int SizeOfThis() { return sizeof(radTNonlinearIsotropMaterial);}

	void FindNewH_FromKsi(TVector3d& InstantH, const TMatrix3d& Matr, const TVector3d& H_Ext, double gInstKsi) 
	{
		TMatrix3d E; E.Str0.x = 1.; E.Str1.y = 1.; E.Str2.z = 1.;
		TMatrix3d BufMatr = gInstKsi*Matr;
		BufMatr = E - BufMatr;
		TMatrix3d InvBufMatr;
		Matrix3d_inv(BufMatr, InvBufMatr);
		InstantH = InvBufMatr*H_Ext;
	}

	int AllocateArrays(int InLenArrayHB)
	{
		DeallocateArrays();
		gLenArrayHB = InLenArrayHB;

		vgArrayHB.resize(gLenArrayHB);
		gArrayHB = vgArrayHB.data();

		vgdBdH.resize(gLenArrayHB);
		gdBdH = vgdBdH.data();
		return 1;
	}
	void DeallocateArrays()
	{
		// RAII: vgArrayHB and vgdBdH cleaned up automatically
		vgArrayHB.clear();
		gArrayHB = 0;
		vgdBdH.clear();
		gdBdH = 0;
		gLenArrayHB = 0;
	}

	static void CopyArrayHB(TVector2d* Dst, TVector2d* Src, int InLenArrayHB, char PrependZero)
	{
		if((Dst == 0) || (Src == 0) || (InLenArrayHB <= 0)) return;

		TVector2d *tArrayHB = Dst, *tInArrayHB = Src;
		if(PrependZero)
		{
			tArrayHB->x = 0.; (tArrayHB++)->y = 0.;
			InLenArrayHB--;
		}
		for(int i=0; i<InLenArrayHB; i++) *(tArrayHB++) = *(tInArrayHB++);
	}
	void CopyArray_dBdH(double* Dst, double* Src, int InLenArrayHB)
	{
		if((Dst == 0) || (Src == 0) || (InLenArrayHB <= 0)) return;

		double *tDst = Dst, *tSrc = Src;
		for(int i=0; i<InLenArrayHB; i++) *(tDst++) = *(tSrc++);
	}
	
	static void CubPln(double Step, double f1, double f2, double fpr1, double fpr2, double* aa)
	{
		double InvStep = 1./Step;
		double f1mf2_d_s1ms2 = (f2 - f1)*InvStep;
		*(aa++) = f1;
		*(aa++) = fpr1;
		*(aa++) = (3.*f1mf2_d_s1ms2 - 2.*fpr1 - fpr2)*InvStep;
		*aa = (-2.*f1mf2_d_s1ms2 + fpr1 + fpr2)*InvStep*InvStep;
	}
};

//-------------------------------------------------------------------------

inline TVector3d radTNonlinearIsotropMaterial::M(const TVector3d& H)
{
	double AbsH = sqrt(H.x*H.x + H.y*H.y + H.z*H.z);
	double AbsM = 0.;
	if(gLenArrayHB == 0)
	{
		// Analytical formula (tanh model)
		for(int i=0; i<lenMs_ks; i++)
			if(Ms[i]!=0.) AbsM += Ms[i]*tanh(ks[i]*AbsH/Ms[i]);
	}
	else
	{
		// ELF-compatible: B-H curve interpolation, returns M = B/mu_0 - H
		AbsM = AbsMvsAbsH_Interpol(AbsH, gArrayHB, gdBdH, gLenArrayHB);
	}

	if(AbsH!=0) return (AbsM/AbsH)*H + RemMagn;
	else return RemMagn;
}

//-------------------------------------------------------------------------

inline void radTNonlinearIsotropMaterial::MultMatrByInstKsiAndMr(const TVector3d& InstantH, const TMatrix3d& Matr, TMatrix3d& MultByKsi, TVector3d& MultByMr)
{
/**
	double AbsInstantH = sqrt(InstantH.x*InstantH.x + InstantH.y*InstantH.y + InstantH.z*InstantH.z);
	double Der, f, InstKsi;
	Der = f = InstKsi = 0.;

	if(gLenArrayHB == 0)
	{
		if(AbsInstantH==0.)
		{
			for(int j=0; j<lenMs_ks; j++) Der += ks[j];
			InstKsi = Der;
		}
		else
		{
			for(int i=0; i<lenMs_ks; i++) if(Ms[i]!=0.) f += Ms[i]*tanh(ks[i]*AbsInstantH/Ms[i]);
			InstKsi = f/AbsInstantH;
		}
	}
	else
	{
		if(AbsInstantH == 0.) InstKsi = *gdMdH;
		//else InstKsi = AbsMvsAbsH_Interpol(AbsInstantH)/AbsInstantH;
		else InstKsi = AbsMvsAbsH_Interpol(AbsInstantH, gArrayHB, gdMdH, gLenArrayHB)/AbsInstantH;
	}
	MultByKsi = InstKsi*Matr; MultByMr = Matr*RemMagn;
**/

	TMatrix3d InstantKsiTensor;
	TVector3d InstantMr;
	DefineInstantKsiTensor(InstantH, InstantKsiTensor, InstantMr);

	MultByKsi = Matr*InstantKsiTensor; 
	MultByMr = Matr*InstantMr;
}

//-------------------------------------------------------------------------

// radTNonlinearIsotropMaterial::Dump / DumpBin REMOVED (Phase B2b/B2c, 2026-04-15)

//-------------------------------------------------------------------------

/**
 * Permanent magnet material with demagnetization curve (Br/Hc model)
 *
 * Magnetic behavior:
 *   B = Br + μ₀·μ_rec·H_parallel
 *   M = B/μ₀ - H = (Br/μ₀) + (μ_rec - 1)·H_parallel
 *
 * where:
 *   Br = residual flux density [T]
 *   Hc = coercivity [A/m]
 *   μ_rec = Br/(μ₀·Hc) = recoil permeability
 *   H_parallel = component of H along easy axis
 */
class radTPermanentMagnet : public radTMaterial {
	double Br;           // Residual flux density [T]
	double Hc;           // Coercivity [A/m]
	double mu_rec;       // Recoil permeability (calculated from Br/Hc)
	TVector3d MagAxis;   // Easy magnetization axis (normalized)

	static constexpr double mu_0 = 1.25663706212e-6;  // Permeability of free space [T/(A/m)]

public:
	radTPermanentMagnet(double InBr, double InHc, const TVector3d& InMagAxis)
		: radTMaterial(TVector3d(0,0,0), 1)  // EasyAxisDefined = 1
	{
		Br = InBr;
		Hc = InHc;
		MagAxis = InMagAxis;

		// Normalize easy axis
		double AbsMagAxis = sqrt(MagAxis.x*MagAxis.x + MagAxis.y*MagAxis.y + MagAxis.z*MagAxis.z);
		if(AbsMagAxis > 0) MagAxis = (1.0/AbsMagAxis) * MagAxis;

		// Calculate recoil permeability: μ_rec = Br / (μ₀·Hc)
		if(Hc > 0) mu_rec = Br / (mu_0 * Hc);
		else mu_rec = 1.0;  // Default to vacuum permeability if Hc = 0

		// Set remanent magnetization Mr = Br/μ₀
		RemMagn = (Br / mu_0) * MagAxis;
	}

	// radTPermanentMagnet(CAuxBinStrVect&) REMOVED (Phase B2c, 2026-04-15)

	radTPermanentMagnet() {}

	int Type_Material() { return 101;}  // New type ID for permanent magnet

	TVector3d M(const TVector3d& H)
	{
		// M = Mr + (μ_rec - 1)·H_parallel·MagAxis
		double H_parallel = H.x*MagAxis.x + H.y*MagAxis.y + H.z*MagAxis.z;  // Dot product
		TVector3d M_induced = (mu_rec - 1.0) * H_parallel * MagAxis;
		return RemMagn + M_induced;
	}

	void DefineInstantKsiTensor(const TVector3d& InstantH, TMatrix3d& InstantKsiTensor, TVector3d& InstantMr)
	{
		// Susceptibility tensor for permanent magnet with demagnetization
		// χ_parallel = μ_rec - 1 along easy axis
		double ksi_par = mu_rec - 1.0;
		double ksi_perp = 0.0;  // No response perpendicular to easy axis

		// Construct anisotropic susceptibility tensor
		TVector3d L = MagAxis;
		double LxLx = L.x*L.x, LyLy = L.y*L.y, LzLz = L.z*L.z;
		double DeltaKsi = ksi_par - ksi_perp;

		TVector3d Str0(ksi_par*LxLx + ksi_perp*(LyLy+LzLz), DeltaKsi*L.x*L.y, DeltaKsi*L.x*L.z);
		TVector3d Str1(Str0.y, ksi_par*LyLy + ksi_perp*(LxLx+LzLz), DeltaKsi*L.y*L.z);
		TVector3d Str2(Str0.z, Str1.z, ksi_par*LzLz + ksi_perp*(LxLx+LyLy));

		InstantKsiTensor.Str0 = Str0;
		InstantKsiTensor.Str1 = Str1;
		InstantKsiTensor.Str2 = Str2;
		InstantMr = RemMagn;
	}

	void MultMatrByInstKsiAndMr(const TVector3d& InstantH, const TMatrix3d& Matr, TMatrix3d& MultByKsi, TVector3d& MultByMr)
	{
		TMatrix3d KsiTensor;
		TVector3d InstantMr;
		DefineInstantKsiTensor(InstantH, KsiTensor, InstantMr);

		MultByKsi = Matr * KsiTensor;
		MultByMr = Matr * RemMagn;
	}

	void FindNewH(TVector3d& H, const TMatrix3d& Matr, const TVector3d& H_Ext, double DesiredPrecOnMagnetizE2)
	{
		// Similar to anisotropic material, but using permanent magnet susceptibility
		TMatrix3d KsiTensor;
		TVector3d InstantMr;
		DefineInstantKsiTensor(H, KsiTensor, InstantMr);

		TVector3d ESt1(1.,0.,0.), ESt2(0.,1.,0.), ESt3(0.,0.,1.);
		TMatrix3d E(ESt1, ESt2, ESt3);
		TMatrix3d BufMatr = E - Matr*KsiTensor;
		TMatrix3d InvBufMatr;
		Matrix3d_inv(BufMatr, InvBufMatr);
		H = InvBufMatr*(H_Ext + Matr*RemMagn);
	}

	int DuplicateItself(radThg& hg, radTApplication*, char)
	{
		return FinishDuplication(new radTPermanentMagnet(*this), hg);
	}

	// Dump / DumpBin REMOVED (Phase B2b/B2c, 2026-04-15)

	int SizeOfThis() { return sizeof(radTPermanentMagnet);}
};

//-------------------------------------------------------------------------

//=========================================================================
// PERMANENT MAGNET MATERIALS
//=========================================================================
//
// Three types of permanent magnet materials are provided:
//
// 1. radTMagFixed (Type 100) - Fixed magnetization, no demagnetization
//    API: rad.MatMagFixed([Mx, My, Mz])
//    Use case: Ideal permanent magnets, fast simulation
//
// 2. radTMagLinear (Type 102) - Linear demagnetization curve (Br/Hc model)
//    API: rad.MatMagLinear(Br, Hc, [ex, ey, ez])  [NOT YET IMPLEMENTED]
//    Use case: NdFeB, SmCo with linear recoil behavior
//
// 3. radTMagCurve (Type 103) - User-defined demagnetization curve
//    API: rad.MatMagCurve([[H1,B1], [H2,B2], ...], [ex, ey, ez])  [NOT YET IMPLEMENTED]
//    Use case: Ferrite, AlNiCo with nonlinear demagnetization
//
// Note: radTPermanentMagnet (Type 101) is the legacy Br/Hc implementation
//       and is retained for backward compatibility.
//=========================================================================

/**
 * Fixed magnetization material (no demagnetization, no relaxation)
 *
 * This material represents an ideal permanent magnet with fixed magnetization
 * that does not change during Solve. The element acts as a pure field source.
 *
 * Magnetic behavior:
 *   M = Mr (constant, set at creation time)
 *   chi = 0 (no susceptibility - magnetization is fixed)
 *
 * Use case:
 *   - Permanent magnets without demagnetization effects
 *   - Elements that should not participate in self-consistent iteration
 *   - Fast simulation when demagnetization is negligible
 *
 * API: rad.MatMagFixed([Mx, My, Mz])
 *   Mx, My, Mz: Magnetization components [A/m]
 *
 * Note: Elements with this material contribute to the external field
 * of other elements but their magnetization is not updated by Solve.
 */
class radTMagFixed : public radTMaterial {
	TVector3d FixedMagn;  // Fixed magnetization [A/m]

public:
	radTMagFixed(const TVector3d& InMagn)
		: radTMaterial(InMagn, 1)  // EasyAxisDefined = 1 (direction defined by Magn)
	{
		FixedMagn = InMagn;
		RemMagn = InMagn;  // RemMagn is the remanent magnetization
	}

	// radTMagFixed(CAuxBinStrVect&) REMOVED (Phase B2c, 2026-04-15)

	radTMagFixed() {}

	int Type_Material() { return 100;}  // Type ID for fixed magnet (non-relaxable)

	// M(H) = Mr (constant, independent of H)
	TVector3d M(const TVector3d& H)
	{
		return FixedMagn;  // Fixed magnetization, does not depend on H
	}

	// Zero susceptibility - magnetization is fixed
	void DefineInstantKsiTensor(const TVector3d& InstantH, TMatrix3d& InstantKsiTensor, TVector3d& InstantMr)
	{
		// chi = 0 (no induced magnetization from H)
		InstantKsiTensor.Str0 = TVector3d(0, 0, 0);
		InstantKsiTensor.Str1 = TVector3d(0, 0, 0);
		InstantKsiTensor.Str2 = TVector3d(0, 0, 0);
		InstantMr = FixedMagn;
	}

	void MultMatrByInstKsiAndMr(const TVector3d& InstantH, const TMatrix3d& Matr, TMatrix3d& MultByKsi, TVector3d& MultByMr)
	{
		// chi = 0, so Matr * chi = 0
		MultByKsi.Str0 = TVector3d(0, 0, 0);
		MultByKsi.Str1 = TVector3d(0, 0, 0);
		MultByKsi.Str2 = TVector3d(0, 0, 0);
		MultByMr = Matr * FixedMagn;
	}

	void FindNewH(TVector3d& H, const TMatrix3d& Matr, const TVector3d& H_Ext, double DesiredPrecOnMagnetizE2)
	{
		// For fixed magnet, H is determined by external field and own remanence
		// With chi = 0, the equation simplifies to: H = H_Ext + Matr * Mr
		H = H_Ext + Matr * FixedMagn;
	}

	int DuplicateItself(radThg& hg, radTApplication*, char)
	{
		return FinishDuplication(new radTMagFixed(*this), hg);
	}

	// Dump / DumpBin REMOVED (Phase B2b/B2c, 2026-04-15)

	int SizeOfThis() { return sizeof(radTMagFixed);}

	// Accessor for the fixed magnetization
	const TVector3d& GetFixedMagn() const { return FixedMagn; }
};

// Backward compatibility alias
typedef radTMagFixed radTFixedMagnet;

//-------------------------------------------------------------------------

/**
 * Linear demagnetization permanent magnet material (Br/Hc model)
 *
 * Magnetic behavior along easy axis:
 *   B = Br + mu_0 * mu_rec * H_parallel
 *   M = B/mu_0 - H = (Br/mu_0) + (mu_rec - 1) * H_parallel
 *
 * where:
 *   Br = residual flux density [T]
 *   Hc = coercivity [A/m]
 *   mu_rec = Br/(mu_0*Hc) = recoil permeability
 *   H_parallel = component of H along easy axis
 *
 * API: rad.MatMagLinear(Br, Hc, [ex, ey, ez])  [NOT YET IMPLEMENTED - uses MatMagFixed]
 *   Br: Residual flux density [T]
 *   Hc: Coercivity [A/m]
 *   [ex, ey, ez]: Easy axis direction (will be normalized)
 *
 * Use case: NdFeB, SmCo permanent magnets with linear recoil behavior
 *
 * Note: Current implementation falls back to radTMagFixed (no demagnetization).
 *       Full implementation will be added in future version.
 */
class radTMagLinear : public radTMaterial {
	double Br;           // Residual flux density [T]
	double Hc;           // Coercivity [A/m]
	double mu_rec;       // Recoil permeability (calculated from Br/Hc)
	TVector3d MagAxis;   // Easy magnetization axis (normalized)

	static constexpr double mu_0 = 1.25663706212e-6;  // Permeability of free space [T/(A/m)]

public:
	radTMagLinear(double InBr, double InHc, const TVector3d& InMagAxis)
		: radTMaterial(TVector3d(0,0,0), 1)  // EasyAxisDefined = 1
	{
		Br = InBr;
		Hc = InHc;
		MagAxis = InMagAxis;

		// Normalize easy axis
		double AbsMagAxis = sqrt(MagAxis.x*MagAxis.x + MagAxis.y*MagAxis.y + MagAxis.z*MagAxis.z);
		if(AbsMagAxis > 0) MagAxis = (1.0/AbsMagAxis) * MagAxis;

		// Calculate recoil permeability: mu_rec = Br / (mu_0*Hc)
		if(Hc > 0) mu_rec = Br / (mu_0 * Hc);
		else mu_rec = 1.0;  // Default to vacuum permeability if Hc = 0

		// Set remanent magnetization Mr = Br/mu_0
		RemMagn = (Br / mu_0) * MagAxis;
	}

	// radTMagLinear(CAuxBinStrVect&) REMOVED (Phase B2c, 2026-04-15)

	radTMagLinear() {}

	int Type_Material() { return 102;}  // Type ID for linear demagnetization PM

	// Current implementation: Fixed magnetization (no demagnetization)
	// TODO: Implement full linear demagnetization model
	TVector3d M(const TVector3d& H)
	{
		// For now, return fixed magnetization (same as radTMagFixed)
		return RemMagn;

		// Full implementation would be:
		// double H_parallel = H.x*MagAxis.x + H.y*MagAxis.y + H.z*MagAxis.z;
		// TVector3d M_induced = (mu_rec - 1.0) * H_parallel * MagAxis;
		// return RemMagn + M_induced;
	}

	void DefineInstantKsiTensor(const TVector3d& InstantH, TMatrix3d& InstantKsiTensor, TVector3d& InstantMr)
	{
		// Current implementation: chi = 0 (fixed magnetization)
		InstantKsiTensor.Str0 = TVector3d(0, 0, 0);
		InstantKsiTensor.Str1 = TVector3d(0, 0, 0);
		InstantKsiTensor.Str2 = TVector3d(0, 0, 0);
		InstantMr = RemMagn;
	}

	void MultMatrByInstKsiAndMr(const TVector3d& InstantH, const TMatrix3d& Matr, TMatrix3d& MultByKsi, TVector3d& MultByMr)
	{
		MultByKsi.Str0 = TVector3d(0, 0, 0);
		MultByKsi.Str1 = TVector3d(0, 0, 0);
		MultByKsi.Str2 = TVector3d(0, 0, 0);
		MultByMr = Matr * RemMagn;
	}

	void FindNewH(TVector3d& H, const TMatrix3d& Matr, const TVector3d& H_Ext, double DesiredPrecOnMagnetizE2)
	{
		H = H_Ext + Matr * RemMagn;
	}

	int DuplicateItself(radThg& hg, radTApplication*, char)
	{
		return FinishDuplication(new radTMagLinear(*this), hg);
	}

	// Dump / DumpBin REMOVED (Phase B2b/B2c, 2026-04-15)

	int SizeOfThis() { return sizeof(radTMagLinear);}

	// Accessors
	double GetBr() const { return Br; }
	double GetHc() const { return Hc; }
	double GetMuRec() const { return mu_rec; }
	const TVector3d& GetMagAxis() const { return MagAxis; }
};

//-------------------------------------------------------------------------

/**
 * User-defined demagnetization curve permanent magnet material
 *
 * Magnetic behavior:
 *   B = f(H) where f is defined by user-provided B-H curve
 *   M = B/mu_0 - H
 *
 * API: rad.MatMagCurve([[H1,B1], [H2,B2], ...], [ex, ey, ez])  [NOT YET IMPLEMENTED - uses MatMagFixed]
 *   [[H1,B1], [H2,B2], ...]: B-H curve data points (H in A/m, B in Tesla)
 *   [ex, ey, ez]: Easy axis direction (will be normalized)
 *
 * Use case: Ferrite, AlNiCo permanent magnets with nonlinear demagnetization curves
 *
 * Note: Current implementation falls back to radTMagFixed (no demagnetization).
 *       Full implementation will be added in future version.
 */
class radTMagCurve : public radTMaterial {
	std::vector<TVector2d> vDemagCurve;  // Demagnetization curve [[H1,B1], [H2,B2], ...]
	TVector2d* DemagCurve;               // Pointer to curve data
	int LenDemagCurve;                   // Number of points in curve
	TVector3d MagAxis;                   // Easy magnetization axis (normalized)

	static constexpr double mu_0 = 1.25663706212e-6;  // Permeability of free space [T/(A/m)]

public:
	radTMagCurve(TVector2d* InCurve, int InLenCurve, const TVector3d& InMagAxis)
		: radTMaterial(TVector3d(0,0,0), 1)
	{
		MagAxis = InMagAxis;

		// Normalize easy axis
		double AbsMagAxis = sqrt(MagAxis.x*MagAxis.x + MagAxis.y*MagAxis.y + MagAxis.z*MagAxis.z);
		if(AbsMagAxis > 0) MagAxis = (1.0/AbsMagAxis) * MagAxis;

		// Copy curve data
		LenDemagCurve = InLenCurve;
		vDemagCurve.resize(LenDemagCurve);
		DemagCurve = vDemagCurve.data();
		for(int i = 0; i < LenDemagCurve; i++)
		{
			DemagCurve[i] = InCurve[i];
		}

		// Set remanent magnetization from Br (B at H=0)
		// Find Br by interpolating to H=0
		double Br = 0.0;
		if(LenDemagCurve >= 2)
		{
			// Assume first point is at H=0 or interpolate
			if(DemagCurve[0].x <= 0.0)
			{
				Br = DemagCurve[0].y;
			}
			else
			{
				// Linear interpolation to H=0
				double H1 = DemagCurve[0].x, B1 = DemagCurve[0].y;
				double H2 = DemagCurve[1].x, B2 = DemagCurve[1].y;
				if(H2 != H1) Br = B1 + (B2 - B1) * (0.0 - H1) / (H2 - H1);
				else Br = B1;
			}
		}
		RemMagn = (Br / mu_0) * MagAxis;
	}

	// radTMagCurve(CAuxBinStrVect&) REMOVED (Phase B2c, 2026-04-15)

	radTMagCurve() { DemagCurve = nullptr; LenDemagCurve = 0; }

	int Type_Material() { return 103;}  // Type ID for user-defined demagnetization curve PM

	// Current implementation: Fixed magnetization (no demagnetization)
	// TODO: Implement full demagnetization curve model
	TVector3d M(const TVector3d& H)
	{
		return RemMagn;
	}

	void DefineInstantKsiTensor(const TVector3d& InstantH, TMatrix3d& InstantKsiTensor, TVector3d& InstantMr)
	{
		InstantKsiTensor.Str0 = TVector3d(0, 0, 0);
		InstantKsiTensor.Str1 = TVector3d(0, 0, 0);
		InstantKsiTensor.Str2 = TVector3d(0, 0, 0);
		InstantMr = RemMagn;
	}

	void MultMatrByInstKsiAndMr(const TVector3d& InstantH, const TMatrix3d& Matr, TMatrix3d& MultByKsi, TVector3d& MultByMr)
	{
		MultByKsi.Str0 = TVector3d(0, 0, 0);
		MultByKsi.Str1 = TVector3d(0, 0, 0);
		MultByKsi.Str2 = TVector3d(0, 0, 0);
		MultByMr = Matr * RemMagn;
	}

	void FindNewH(TVector3d& H, const TMatrix3d& Matr, const TVector3d& H_Ext, double DesiredPrecOnMagnetizE2)
	{
		H = H_Ext + Matr * RemMagn;
	}

	radTMagCurve(const radTMagCurve& other)
		: radTMaterial(other), vDemagCurve(other.vDemagCurve),
		  LenDemagCurve(other.LenDemagCurve), MagAxis(other.MagAxis)
	{
		DemagCurve = vDemagCurve.data();
	}

	int DuplicateItself(radThg& hg, radTApplication*, char)
	{
		return FinishDuplication(new radTMagCurve(*this), hg);
	}

	// Dump / DumpBin REMOVED (Phase B2b/B2c, 2026-04-15)

	int SizeOfThis() { return sizeof(radTMagCurve);}

	// Accessors
	int GetCurveLength() const { return LenDemagCurve; }
	const TVector3d& GetMagAxis() const { return MagAxis; }
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
// Energy-based vector hysteresis material (Francois-Lavet / Egger formulation)
// Type_Material() = 5
//
// Forward operator: B -> H, O(K) via Schur complement (natural for B-input Play)
// Inverse operator: H -> B, O(K) per evaluation (each J_k independent)
//
// Reference:
//   Egger, Engertsberger, Schafelner: "Efficient evaluation of forward and
//   inverse energy-based magnetic hysteresis operators", MAGCON-25-07-0171
//-------------------------------------------------------------------------

// Per-operator shape function table (replaces analytical A_s*tan formula)
//-------------------------------------------------------------------------
// Abstract base class for hysteresis materials
// Provides common interface for Energy-based and Play-based models.
//-------------------------------------------------------------------------
class radTHysteresisMaterial : public radTMaterial {
public:
	virtual ~radTHysteresisMaterial() {}

	// Core operators (B-input naming)
	virtual TVector3d Forward(const TVector3d& B) = 0;    // B -> H
	virtual TVector3d Inverse(const TVector3d& H) = 0;    // H -> B
	virtual void ComputeJacobian(TMatrix3d& dBdH, double& chi_d) const = 0;

	// State management
	virtual void SaveState(std::vector<TVector3d>& saved) const = 0;
	virtual void RestoreState(const std::vector<TVector3d>& saved) = 0;
	virtual void ResetState() = 0;
	virtual void CommitState() = 0;

	// Array-based state (for Python API)
	virtual int GetStateSize() const = 0;
	virtual void SaveStateToArray(double* pState) const = 0;
	virtual void RestoreStateFromArray(const double* pState) = 0;

	// Chi computation (for relaxation solver)
	virtual double ComputeChiFromH(const TVector3d& H) = 0;
	virtual double ComputeDifferentialChi(double H_mag) = 0;

	// B-input chi: chi = |M|/|H| where H = Forward(B), M = B/mu0 - H.
	// Implemented once on the base via the virtual Forward(B) so it works for
	// BOTH the Play model (Type 6, O(K) direct Forward) and the Energy model
	// (Type 5, which is a Forward-equivalent subclass).  NOTE: Forward(B)
	// advances the material's m_pk_current/m_last_* state -- the moment B-input
	// Picard saves/restores per-element state around each call so each element
	// keeps its own play trajectory (see UpdateChiAndCheckConvergence).
	//
	// Mirrors the clamping/return conventions of ComputeChiFromH (which the
	// Play/Energy override as B = Inverse(H); chi = |B|/(mu0*|H|) - 1).  Here
	// we keep chi >= 0 and clamp to a large finite ceiling so the moment system
	// matrix diagonal (1/chi) never blows up.
	virtual double ComputeChiFromB(const TVector3d& B)
	{
		const double MU_0_local = 4.0 * 3.14159265358979323846 * 1.0e-7;
		const double INV_MU_0_local = 1.0 / MU_0_local;
		TVector3d H = Forward(B);                 // B -> H (advances play state)
		double Hn = sqrt(H.x*H.x + H.y*H.y + H.z*H.z);
		// M = B/mu0 - H
		TVector3d M(B.x*INV_MU_0_local - H.x,
		            B.y*INV_MU_0_local - H.y,
		            B.z*INV_MU_0_local - H.z);
		double Mn = sqrt(M.x*M.x + M.y*M.y + M.z*M.z);
		if(Hn < 1.0e-30) return GetInitialChi_ELF_Style();
		double chi = Mn / Hn;
		if(chi < 0.0) chi = 0.0;
		if(chi > 1.0e10) chi = 1.0e10;
		return chi;
	}

	// Initialization helpers
	virtual double GetInitialChi_ELF_Style() const = 0;
	virtual double GetBsaturation() const = 0;
};

// radTEnergyHysteresisMaterial (Type 5) is defined as a thin subclass of
// radTPlayHysteresisMaterial below, AFTER Type 6 is fully declared.

//-------------------------------------------------------------------------
// B-input Play hysteresis material
//   Direct play model: H = sum f_k(|p_k|) * p_k/|p_k|
//   where p_k = play(B, eta_k, p_k_prev)
//   Forward(B->H) is O(K), no Newton iteration needed.
//   Inverse(H->B) uses Newton with analytical Jacobian.
//-------------------------------------------------------------------------

// Per-operator shape function table for play model
struct PlayTable {
	std::vector<double> r;   // |p| grid points [0, r_max], monotonically increasing
	std::vector<double> f;   // f_k(r) shape function H contribution (can be negative)
	std::vector<double> df;  // f_k'(r) derivative (for Jacobian)
	int n;                   // Number of grid points
	double r_max;            // Maximum |p| = r[n-1]
};

class radTPlayHysteresisMaterial : public radTHysteresisMaterial {

	// Model parameters (immutable after construction)
	int m_K;                              // Number of play operators
	std::vector<PlayTable> m_tables;      // Per-operator shape function tables
	std::vector<double> m_eta;            // Play thresholds [Tesla]

	// Internal state: play operator outputs
	std::vector<TVector3d> m_pk_prev;     // Committed play outputs [K]
	std::vector<TVector3d> m_pk_pinning;  // Pinning reference for current step [K]
	std::vector<TVector3d> m_pk_current;  // Last Forward result [K]

	// Monotone limit: largest B where dH/dB > 0 on virgin curve
	double m_B_mono_max;                  // B at monotone boundary [Tesla]
	double m_H_mono_max;                  // H(B_mono_max) [A/m]

	// Energy-based decomposition: H = nu_rev*B + H_irr(B)
	double m_nu_rev;                      // Reversible reluctivity [A/m/T]

	// Cached results from last Forward() call
	TMatrix3d m_last_dHdB;                // Analytical Jacobian dH/dB
	TVector3d m_last_H;
	TVector3d m_last_B;
	double m_last_chi;
	double m_last_chi_d;
	bool m_has_result;

public:
	// Constructor from table-based shape functions
	radTPlayHysteresisMaterial(int K, const double* eta,
	                            const std::vector<std::vector<double>>& r_tables,
	                            const std::vector<std::vector<double>>& f_tables);

	// Serialization constructor radTPlayHysteresisMaterial(CAuxBinStrVect&) REMOVED (Phase B2c, 2026-04-15)

	// Default constructor
	radTPlayHysteresisMaterial()
		: m_K(0), m_B_mono_max(1.0), m_H_mono_max(0), m_nu_rev(0),
		  m_last_chi(0), m_last_chi_d(0), m_has_result(false) {}

	// Precompute monotone limit on virgin curve
	void ComputeMonotoneLimits();

	int Type_Material() { return 6; }

	// radTMaterial virtual overrides
	TVector3d M(const TVector3d& H);
	void DefineInstantKsiTensor(const TVector3d& H, TMatrix3d& KsiTensor, TVector3d& Mr);

	// Solver integration
	double ComputeChiFromH(const TVector3d& H) override;
	double ComputeChiDualMethod(double H_mag, double mu_old, double relax = 0.0);
	double ComputeDifferentialChi(double H_mag) override;

	double GetInitialChi_ELF_Style() const override
	{
		// Play model: H = sum_k f_k(|p_k|) * p_hat
		// For small B (all operators following): dH/dB = sum_k f'_k(0)
		// So chi = mu_r - 1 = 1/(mu_0 * dH/dB) - 1
		if(m_K <= 0 || m_tables.empty()) return 1.0;
		const double MU_0_local = 4.0 * 3.14159265358979323846 * 1.0e-7;
		double dHdB_sum = 0.0;
		for(int k = 0; k < m_K; k++)
		{
			const auto& tab = m_tables[k];
			if(tab.n >= 2 && tab.r[1] > 1e-30)
				dHdB_sum += (tab.f[1] - tab.f[0]) / tab.r[1];
		}
		if(fabs(dHdB_sum) < 1e-10) return 1.0;
		double chi_init = 1.0 / (MU_0_local * fabs(dHdB_sum)) - 1.0;
		if(chi_init < 1.0) chi_init = 1.0;
		return chi_init;
	}

	double GetBsaturation() const override
	{
		double rmax = 0.0;
		for(int k = 0; k < m_K; k++)
			if(m_tables[k].r_max > rmax) rmax = m_tables[k].r_max;
		if(rmax < 1e-10) return 1.0;
		return rmax;
	}

	// Core operators
	TVector3d Forward(const TVector3d& B) override;    // B -> H, O(K) direct
	TVector3d Inverse(const TVector3d& H) override;    // H -> B, Newton with analytical Jacobian

	void ComputeJacobian(TMatrix3d& dBdH, double& chi_d) const override;

	// Energy-based decomposition: H = nu_rev*B + H_irr(B)
	void ComputeNuRev();                               // Called after construction
	double GetNuRev() const { return m_nu_rev; }       // For Hantila solver
	TVector3d Irreversible(const TVector3d& B);        // H_irr = Forward(B) - nu_rev*B

	// State management
	void SaveState(std::vector<TVector3d>& saved) const override
	{
		saved.resize(m_K);
		for(int k = 0; k < m_K; k++) saved[k] = m_pk_prev[k];
	}
	void RestoreState(const std::vector<TVector3d>& saved) override
	{
		for(int k = 0; k < m_K; k++)
		{
			m_pk_prev[k] = saved[k];
			m_pk_pinning[k] = saved[k];
		}
		m_has_result = false;
	}
	void ResetState() override
	{
		TVector3d zero(0, 0, 0);
		for(int k = 0; k < m_K; k++)
		{
			m_pk_prev[k] = zero;
			m_pk_pinning[k] = zero;
			m_pk_current[k] = zero;
		}
		m_has_result = false;
	}
	void CommitState() override
	{
		for(int k = 0; k < m_K; k++)
		{
			m_pk_prev[k] = m_pk_current[k];
			m_pk_pinning[k] = m_pk_current[k];
		}
	}

	// State size = K*9 (Play states) + 7 (warm-start cache: m_last_B[3] +
	// m_last_H[3] + m_has_result flag).
	//
	// Bug fix 2026-05-02: prior format omitted the warm-start cache so
	// RestoreStateFromArray was NOT state-equivalent to "having driven the
	// material to that point". On hard branches (post-saturation
	// crossings), the next Inverse cold-started instead of warmstarting
	// from m_last_B → Newton landed in a different basin. Cross-validation
	// observed a 2.83e6 A/m sign flip in M (Type 6) at H = (2000, 500,
	// -300) post-restore vs the equivalent in-line drive sequence.
	int GetStateSize() const override { return m_K * 9 + 7; }
	void SaveStateToArray(double* pState) const override
	{
		int idx = 0;
		for(int k = 0; k < m_K; k++)
		{
			pState[idx++] = m_pk_prev[k].x; pState[idx++] = m_pk_prev[k].y; pState[idx++] = m_pk_prev[k].z;
			pState[idx++] = m_pk_pinning[k].x; pState[idx++] = m_pk_pinning[k].y; pState[idx++] = m_pk_pinning[k].z;
			pState[idx++] = m_pk_current[k].x; pState[idx++] = m_pk_current[k].y; pState[idx++] = m_pk_current[k].z;
		}
		pState[idx++] = m_last_B.x; pState[idx++] = m_last_B.y; pState[idx++] = m_last_B.z;
		pState[idx++] = m_last_H.x; pState[idx++] = m_last_H.y; pState[idx++] = m_last_H.z;
		pState[idx++] = m_has_result ? 1.0 : 0.0;
	}
	void RestoreStateFromArray(const double* pState) override
	{
		int idx = 0;
		for(int k = 0; k < m_K; k++)
		{
			m_pk_prev[k].x = pState[idx++]; m_pk_prev[k].y = pState[idx++]; m_pk_prev[k].z = pState[idx++];
			m_pk_pinning[k].x = pState[idx++]; m_pk_pinning[k].y = pState[idx++]; m_pk_pinning[k].z = pState[idx++];
			m_pk_current[k].x = pState[idx++]; m_pk_current[k].y = pState[idx++]; m_pk_current[k].z = pState[idx++];
		}
		m_last_B.x = pState[idx++]; m_last_B.y = pState[idx++]; m_last_B.z = pState[idx++];
		m_last_H.x = pState[idx++]; m_last_H.y = pState[idx++]; m_last_H.z = pState[idx++];
		m_has_result = (pState[idx++] >= 0.5);
	}

	// DumpBin REMOVED (Phase B2c, 2026-04-15)

	int DuplicateItself(radThg& hg, radTApplication*, char)
	{
		radTSend Send;
		radTPlayHysteresisMaterial* pNew = new radTPlayHysteresisMaterial();
		if(pNew == 0) { Send.ErrorMessage("Radia::Error900"); return 0; }
		*pNew = *this;
		return FinishDuplication(pNew, hg);
	}

	int SizeOfThis() { return sizeof(radTPlayHysteresisMaterial); }

private:
	// Table interpolation (same binary search as energy model)
	double fk(int k, double r_mag) const;   // f_k(r) shape function
	double dfk(int k, double r_mag) const;  // f_k'(r) derivative

	// Precompute df table from f data
	static void PrecomputePlayTable(PlayTable& tab);

	// Helpers (reuse from energy model)
	static TMatrix3d OuterProduct(const TVector3d& a, const TVector3d& b)
	{
		return TMatrix3d(
			TVector3d(a.x*b.x, a.x*b.y, a.x*b.z),
			TVector3d(a.y*b.x, a.y*b.y, a.y*b.z),
			TVector3d(a.z*b.x, a.z*b.y, a.z*b.z));
	}
	static TMatrix3d Eye()
	{
		return TMatrix3d(
			TVector3d(1, 0, 0),
			TVector3d(0, 1, 0),
			TVector3d(0, 0, 1));
	}
};

//-------------------------------------------------------------------------
// Energy-based vector hysteresis material (Type 5)
//
// Refactor 2026-05-XX: Type 5 is now a thin subclass of Type 6
// (radTPlayHysteresisMaterial). The original Schur-complement Newton
// implementation (Henrotte–Egger formulation with per-hysteron J_k as
// independent variables) is structurally incompatible with the
// Potter B-input shape functions identified by Hane–Sugahara: f_k < 0
// for k >= 1 makes the per-hysteron Hessian h^priv_k = f'_k
// sign-indefinite, the Sherman–Morrison Schur scalar can vanish or
// change sign, and the Newton descent direction reverses.
//
// The current implementation uses Type 6's 3D Newton with analytical
// Jacobian on F(B) = nu_rev*B + sum_k g_k(|p_k|) p_k/|p_k| - H = 0,
// which preserves the energy formulation interpretation
// W(B) = (1/2) nu_rev |B|^2 + sum_k G_k(|p_k|)
// where g_0 = f_0 - nu_rev*r and g_k = f_k for k >= 1. The forward
// map is algebraically identical to direct Play (H-axis congruency
// preserved exactly), and nu_rev is exposed for Hantila one-time-LU
// coupling in MMM/MSC nonlinear iteration.
//
// API: MatEnergyHysteresis(K, chi, tables, eps). The chi parameter is
// numerically equal to the Play threshold eta_k by Hane–Sugahara
// identification convention; it is forwarded to the Type 6 base
// constructor as eta. The eps parameter is retained for ABI
// compatibility but is no longer used (the Bergqvist regularization
// |x|_eps had no role in the rev/irrev separated formulation).
//-------------------------------------------------------------------------

class radTEnergyHysteresisMaterial : public radTPlayHysteresisMaterial {
public:
	// Constructor: chi -> eta (numerically identical), eps unused.
	radTEnergyHysteresisMaterial(int K, const double* chi,
	                             const std::vector<std::vector<double>>& r_tables,
	                             const std::vector<std::vector<double>>& f_tables,
	                             double /*eps*/)
		: radTPlayHysteresisMaterial(K, chi, r_tables, f_tables) {}

	// Default constructor
	radTEnergyHysteresisMaterial() : radTPlayHysteresisMaterial() {}

	// Distinct material type ID (5) — distinguishes from Type 6 (Play)
	int Type_Material() { return 5; }

	// DuplicateItself must produce another Type 5 instance, not Type 6
	int DuplicateItself(radThg& hg, radTApplication*, char)
	{
		radTSend Send;
		radTEnergyHysteresisMaterial* pNew = new radTEnergyHysteresisMaterial();
		if(pNew == 0) { Send.ErrorMessage("Radia::Error900"); return 0; }
		*pNew = *this;
		return FinishDuplication(pNew, hg);
	}

	int SizeOfThis() { return sizeof(radTEnergyHysteresisMaterial); }

	// All other methods (M, DefineInstantKsiTensor, ComputeChiFromH,
	// ComputeDifferentialChi, GetInitialChi_ELF_Style, GetBsaturation,
	// Forward, Inverse, ComputeJacobian, SaveState, RestoreState,
	// ResetState, CommitState, GetStateSize, SaveStateToArray,
	// RestoreStateFromArray) are inherited from radTPlayHysteresisMaterial.
};

//-------------------------------------------------------------------------

#endif
