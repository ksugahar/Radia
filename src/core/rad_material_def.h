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

	void Dump(std::ostream& o, int ShortSign =0) // Porting
	{
		radTg::Dump(o);
		o << "Magnetic material: ";
	}

	void DumpBin_Material(CAuxBinStrVect& oStr)
	{
		//static radTHMatDBVect MaterDB; //no need to dump static members
		//TVector3d RemMagn; // Don't make it private nor protected
		oStr << RemMagn;
		
		//char EasyAxisDefined;
		oStr << EasyAxisDefined;
	}

	void DumpBinParse_Material(CAuxBinStrVect& inStr)
	{
		//TVector3d RemMagn; // Don't make it private nor protected
		inStr >> RemMagn;

		//char EasyAxisDefined;
		inStr >> EasyAxisDefined;
	}

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

	radTLinearAnisotropMaterial(CAuxBinStrVect& inStr) //, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
	{//Instantiates from string according to DumpBin
		DumpBinParse_Material(inStr);

		//double KsiPar, KsiPerp;
		inStr >> KsiPar;
		inStr >> KsiPerp;

		//TMatrix3d KsiTensor;
		inStr >> KsiTensor;
	}

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

	inline void Dump(std::ostream& o, int ShortSign =0); // Porting

	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
	//void DumpBin(CAuxBinStrVect& oStr, radTmhg& mEl, radThg& hg)
	{
		//int newKey = (int)mEl.size() + 1;
		//mEl[newKey] = hg;
		//Start dumping this object
		//oStr << newKey;

		vElemKeysOut.push_back(elemKey);
		oStr << elemKey;

		//Next 5 bytes define/encode element type:
		oStr << (char)Type_g();
		oStr << (char)Type_Material();
		oStr << (char)0;
		oStr << (char)0;
		oStr << (char)0;

		//Members of radTMaterial
		DumpBin_Material(oStr);

		//double KsiPar, KsiPerp;
		oStr << KsiPar << KsiPerp;

		//TMatrix3d KsiTensor;
		oStr << KsiTensor;
	}

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

inline void radTLinearAnisotropMaterial::Dump(std::ostream& o, int ShortSign) //Porting
{
	radTMaterial::Dump(o);
	o << "Linear anisotropic";

	if(ShortSign==1) return;
	o << endl;
	o << "   {ksipar,ksiper}= {" << KsiPar << ',' << KsiPerp << "}" << endl;

	if(EasyAxisDefined)
		o << "   {mrx,mry,mrz}= {" << RemMagn.x << ',' << RemMagn.y << ',' << RemMagn.z << "}";
	else
		o << "   mr= " << RemMagn.x;

	o << endl;
	o << "   Memory occupied: " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTLinearIsotropMaterial : public radTMaterial {
	double Ksi;

public:
	radTLinearIsotropMaterial(double InKsi) { Ksi = InKsi;}
	radTLinearIsotropMaterial(const double* InKsiArray, const TVector3d& InRemMagn, char InEasyAxisDefined) 
		: radTMaterial(InRemMagn, InEasyAxisDefined) { Ksi = InKsiArray[0];}
	
	radTLinearIsotropMaterial(CAuxBinStrVect& inStr) //, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
	{//Instantiates from string according to DumpBin
		DumpBinParse_Material(inStr);
		//double Ksi;
		inStr << Ksi;
	}

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

	void Dump(std::ostream& o, int ShortSign =0) // Porting
	//inline void radTLinearIsotropMaterial::Dump(std::ostream& o, int ShortSign) // Porting
	{
		radTMaterial::Dump(o);
		o << "Linear isotropic";

		if(ShortSign==1) return;
		o << endl;
		o << "   ksi= " << Ksi;

		o << endl;
		o << "   Memory occupied: " << SizeOfThis() << " bytes";
	}

	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
	{
		vElemKeysOut.push_back(elemKey);
		oStr << elemKey;

		//Next 5 bytes define/encode element type:
		oStr << (char)Type_g();
		oStr << (char)Type_Material();
		oStr << (char)0;
		oStr << (char)0;
		oStr << (char)0;

		//Members of radTMaterial
		DumpBin_Material(oStr);

		//double Ksi;
		oStr << Ksi;
	}

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
	radTNonlinearIsotropMaterial(CAuxBinStrVect& inStr) //, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
	{//Instantiates from string according to DumpBin
		DumpBinParse_Material(inStr);

		//double Ms[3];
		inStr >> Ms[0]; inStr >> Ms[1]; inStr >> Ms[2];

		//double ks[3];
		inStr >> ks[0]; inStr >> ks[1]; inStr >> ks[2];

		//int lenMs_ks;
		inStr >> lenMs_ks;

		//int gLenArrayHB;
		inStr >> gLenArrayHB;

		//TVector2d* gArrayHB;
		gArrayHB = 0;
		char cTest=0;
		inStr >> cTest;
		if(cTest > 0)
		{
			vgArrayHB.resize(gLenArrayHB);
			gArrayHB = vgArrayHB.data();
			TVector2d *t_gArrayHB = gArrayHB;
			for(int i=0; i<gLenArrayHB; i++) inStr >> (*(t_gArrayHB++));
		}
		//double* gdBdH;
		gdBdH = 0;
		inStr >> cTest;
		if(cTest > 0)
		{
			vgdBdH.resize(gLenArrayHB);
			gdBdH = vgdBdH.data();
			double *t_gdBdH = gdBdH;
			for(int i=0; i<gLenArrayHB; i++) inStr >> (*(t_gdBdH++));
		}

		//double gMaxKsi;
		inStr >> gMaxKsi;
	}

	radTNonlinearIsotropMaterial() { gArrayHB = 0; gLenArrayHB = 0; gdBdH = 0; gMaxKsi = 0;}
	~radTNonlinearIsotropMaterial() { DeallocateArrays();}

	int Type_Material() { return 3;}

	TVector3d M(const TVector3d& H);
	void DefineInstantKsiTensor(const TVector3d&, TMatrix3d&, TVector3d&);
	void MultMatrByInstKsiAndMr(const TVector3d&, const TMatrix3d&, TMatrix3d&, TVector3d&);
	//void FindNewH(TVector3d&, const TMatrix3d&, const TVector3d&, double, radTg3dRelax*, void*); //OC140103
	void FindNewH(TVector3d&, const TMatrix3d&, const TVector3d&, double); //OC140103

	inline void Dump(std::ostream& o, int ShortSign =0);
	inline void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey);

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

inline void radTNonlinearIsotropMaterial::Dump(std::ostream& o, int ShortSign) // Porting
{
	radTMaterial::Dump(o);
	o << "Nonlinear isotropic";

	if(ShortSign==1) return;

	o << endl;
	if((gArrayHB == 0) || (gLenArrayHB == 0))
	{
		o << "   {ms1,ms2,ms3}= {" << Ms[0] << ',' << Ms[1] << ',' << Ms[2] << "}" << endl;
		o << "   {ks1,ks2,ks3}= {" << ks[0] << ',' << ks[1] << ',' << ks[2] << "}";
	}
	else
	{
		o << "   M(H) defined by table of values";
	}

	o << endl;
	o << "   Memory occupied: " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------

inline void radTNonlinearIsotropMaterial::DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
{
	vElemKeysOut.push_back(elemKey);
	oStr << elemKey;

	//Next 5 bytes define/encode element type:
	oStr << (char)Type_g();
	oStr << (char)Type_Material();
	oStr << (char)0;
	oStr << (char)0;
	oStr << (char)0;

	//Members of radTMaterial
	DumpBin_Material(oStr);

	//double Ms[3];
	oStr << Ms[0] << Ms[1] << Ms[2];

	//double ks[3];
	oStr << ks[0] << ks[1] << ks[2];

	//int lenMs_ks;
	oStr << lenMs_ks;

	//int gLenArrayHB;
	oStr << gLenArrayHB;

	//TVector2d* gArrayHB;
	if((gLenArrayHB > 0) && (gArrayHB != 0))
	{
		oStr << (char)1;
		TVector2d *t_gArrayHB = gArrayHB;
		for(int i=0; i<gLenArrayHB; i++) oStr << (*(t_gArrayHB++));
	}
	else oStr << (char)0;

	//double* gdBdH;
	if((gLenArrayHB > 0) && (gdBdH != 0))
	{
		oStr << (char)1;
		double *t_gdBdH = gdBdH;
		for(int i=0; i<gLenArrayHB; i++) oStr << (*(t_gdBdH++));
	}
	else oStr << (char)0;

	//double gMaxKsi;
	oStr << gMaxKsi;
}

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

	radTPermanentMagnet(CAuxBinStrVect& inStr)
	{
		// Instantiates from string according to DumpBin
		DumpBinParse_Material(inStr);

		inStr >> Br;
		inStr >> Hc;
		inStr >> mu_rec;
		inStr >> MagAxis;
	}

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

	void Dump(std::ostream& o, int ShortSign =0)
	{
		radTMaterial::Dump(o);
		o << "Permanent magnet (Br/Hc)";

		if(ShortSign==1) return;
		o << endl;
		o << "   Br= " << Br << " T, Hc= " << Hc << " A/m, mu_rec= " << mu_rec << endl;
		o << "   Easy axis: {" << MagAxis.x << ',' << MagAxis.y << ',' << MagAxis.z << "}" << endl;
		o << "   Remanent magnetization: {" << RemMagn.x << ',' << RemMagn.y << ',' << RemMagn.z << "} A/m" << endl;
		o << "   Memory occupied: " << SizeOfThis() << " bytes";
	}

	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
	{
		vElemKeysOut.push_back(elemKey);
		oStr << elemKey;

		// Next 5 bytes define/encode element type
		oStr << (char)Type_g();
		oStr << (char)Type_Material();
		oStr << (char)0;
		oStr << (char)0;
		oStr << (char)0;

		// Members of radTMaterial
		DumpBin_Material(oStr);

		// Members of radTPermanentMagnet
		oStr << Br;
		oStr << Hc;
		oStr << mu_rec;
		oStr << MagAxis;
	}

	int SizeOfThis() { return sizeof(radTPermanentMagnet);}
};

//-------------------------------------------------------------------------

#endif
