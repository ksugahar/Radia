/*-------------------------------------------------------------------------
*
* File name:      radmtra1.h
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

#ifndef __RADMTRA1_H
#define __RADMTRA1_H

#include "rad_material_def.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTNonlinearAnisotropMaterial : public radTMaterial {

protected:
	double Ksi[2][4], Ms[2][3], Hci[4];
	char DependenceIsNonlinear[2];

	std::vector<TVector2d> vgArrayHB_Par, vgArrayHB_Perp;
	TVector2d *gArrayHB_Par, *gArrayHB_Perp;
	int gLenArrayHB_Par, gLenArrayHB_Perp;
	std::vector<double> vgdBdH_Par, vgdBdH_Perp;
	double *gdBdH_Par, *gdBdH_Perp;
	double gMaxKsi_Par, gMaxKsi_Perp;

	TVector3d UnitEasyAxisVect;

public:

	radTNonlinearAnisotropMaterial(double** InKsi, double** InMs, double* InHci, char* InDependenceIsNonlinear)
	{
		gLenArrayHB_Par = gLenArrayHB_Perp = 0;
		gArrayHB_Par = gArrayHB_Perp = 0;
		gdBdH_Par = gdBdH_Perp = 0;

		EasyAxisDefined = 0;
		for(int k=0; k<2; k++)
		{
			DependenceIsNonlinear[k] = InDependenceIsNonlinear[k];
			if(DependenceIsNonlinear[k])
			{
				double *CurInKsi = InKsi[k], *CurInMs = InMs[k], *CurInHci = InHci;
				double *CurKsi = Ksi[k], *CurMs = Ms[k], *CurHci = Hci;
				for(int i=0; i<3; i++)
				{
					*(CurKsi++) = *(CurInKsi++); *(CurMs++) = *(CurInMs++); 
					if(k == 0) *(CurHci++) = *(CurInHci++);
				}
				*CurKsi = *CurInKsi;
				if(k == 0) *CurHci = *CurInHci;
				//Hc[k] = InHc[k];
			}
			else *(Ksi[k]) = *(InKsi[k]);
		}
	}

	radTNonlinearAnisotropMaterial(CAuxBinStrVect& inStr)
	{//Instantiates from string according to DumpBin
		DumpBinParse_Material(inStr);
		//Members of radTNonlinearAnisotropMaterial
		DumpBinParse_NonlinearAnisotropMaterial(inStr);
	}

	radTNonlinearAnisotropMaterial() 
	{
		gLenArrayHB_Par = gLenArrayHB_Perp = 0;
		gArrayHB_Par = gArrayHB_Perp = 0;
		gdBdH_Par = gdBdH_Perp = 0;
	}

	int Type_Material() { return 4;}
	virtual int Type_NonlinearAnisotropMaterial() { return 0;}

	inline TVector3d M(const TVector3d& H);
	inline void DefineInstantKsiTensor(const TVector3d&, TMatrix3d&, TVector3d&);
	inline void MultMatrByInstKsiAndMr(const TVector3d&, const TMatrix3d&, TMatrix3d&, TVector3d&);
	//inline void FindNewH(TVector3d&, const TMatrix3d&, const TVector3d&, double, radTg3dRelax*, void*); //OC140103
	inline void FindNewH(TVector3d&, const TMatrix3d&, const TVector3d&, double);

	inline void SetKsiTensor(double, double, TMatrix3d&);
	inline double ScalarInstantKsi(double, char, TVector3d&);
	inline double ScalarInstantKsiAsForIsotropic(double, char);

	int FinishSetup(TVector3d& Magn)
	{
		if(!EasyAxisDefined)
		{
			double AbsLocMagn = sqrt(Magn.x*Magn.x + Magn.y*Magn.y + Magn.z*Magn.z);

			radTSend Send;
			const double AbsTol = 1.E-10;
			if(AbsLocMagn < AbsTol) { Send.ErrorMessage("Radia::Error107"); return 0;}

			UnitEasyAxisVect = (1./AbsLocMagn)*Magn;
			EasyAxisDefined = 1;
		}
		return 1;
	}

	int DuplicateItself(radThg& hg, radTApplication*, char) 
	{
		return FinishDuplication(new radTNonlinearAnisotropMaterial(*this), hg);
	}

	inline void Dump(std::ostream& o, int ShortSign =0);
	inline void DumpBin_NonlinearAnisotropMaterial(CAuxBinStrVect& oStr);
	inline void DumpBinParse_NonlinearAnisotropMaterial(CAuxBinStrVect& inStr);
	inline void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey);

	int SizeOfThis() { return sizeof(radTNonlinearAnisotropMaterial);}

	int AllocateArrays(int InLenArrayHB_Par, int InLenArrayHB_Perp)
	{
		DeallocateArrays();
		if(InLenArrayHB_Par > 0)
		{
			vgArrayHB_Par.resize(InLenArrayHB_Par);
			gArrayHB_Par = vgArrayHB_Par.data();
			vgdBdH_Par.resize(InLenArrayHB_Par);
			gdBdH_Par = vgdBdH_Par.data();
			gLenArrayHB_Par = InLenArrayHB_Par;
		}
		if(InLenArrayHB_Perp > 0)
		{
			vgArrayHB_Perp.resize(InLenArrayHB_Perp);
			gArrayHB_Perp = vgArrayHB_Perp.data();
			vgdBdH_Perp.resize(InLenArrayHB_Perp);
			gdBdH_Perp = vgdBdH_Perp.data();
			gLenArrayHB_Perp = InLenArrayHB_Perp;
		}
		return 1;
	}
	void DeallocateArrays()
	{
		// RAII: vgArrayHB_Par and vgdBdH_Par cleaned up automatically
		vgArrayHB_Par.clear();
		gArrayHB_Par = 0;
		vgdBdH_Par.clear();
		gdBdH_Par = 0;
		gLenArrayHB_Par = 0;

		// RAII: vgArrayHB_Perp and vgdBdH_Perp cleaned up automatically
		vgArrayHB_Perp.clear();
		gArrayHB_Perp = 0;
		vgdBdH_Perp.clear(); 
		gdBdH_Perp = 0;
		gLenArrayHB_Perp = 0;
	}

	void SetupEasyAxisIfPossible(double *pN)
	{
		EasyAxisDefined = 0;
		if(pN != 0)
		{
			double Ne2 = pN[0]*pN[0] + pN[1]*pN[1] + pN[2]*pN[2];
			if(Ne2 > 0)
			{
				double InvAbsN = 1./sqrt(Ne2);
				UnitEasyAxisVect.x = InvAbsN*pN[0];
				UnitEasyAxisVect.y = InvAbsN*pN[1];
				UnitEasyAxisVect.z = InvAbsN*pN[2];
				EasyAxisDefined = 1;
			}
		}
	}
};

//-------------------------------------------------------------------------

inline void radTNonlinearAnisotropMaterial::SetKsiTensor(double KsiPar, double KsiPerp, TMatrix3d& KsiTensor)
{
	TVector3d& L = UnitEasyAxisVect;
	double DeltaKsi = KsiPar-KsiPerp;
	double LxLx, LyLy, LzLz;
	LxLx=L.x*L.x; LyLy=L.y*L.y; LzLz=L.z*L.z;
	TVector3d Str0(KsiPar*LxLx+KsiPerp*(LyLy+LzLz), DeltaKsi*L.x*L.y, DeltaKsi*L.x*L.z);
	TVector3d Str1(Str0.y, KsiPar*LyLy+KsiPerp*(LxLx+LzLz), DeltaKsi*L.y*L.z);
	TVector3d Str2(Str0.z, Str1.z, KsiPar*LzLz+KsiPerp*(LxLx+LyLy));
	KsiTensor.Str0 = Str0; KsiTensor.Str1 = Str1; KsiTensor.Str2 = Str2;
}

//-------------------------------------------------------------------------

inline TVector3d radTNonlinearAnisotropMaterial::M(const TVector3d& InstantH)
{
	TVector3d M_Par, M_Per;
	if(DependenceIsNonlinear[0])
	{
		double Hpar = UnitEasyAxisVect*InstantH;
		double *ks = Ksi[0], *ms = Ms[0], *hc = Hci;  //hc = *Hc;
		//double ScalarH_mi_hc = Hpar - hc, f = 0.;
		double f = 0.;

		for(int j=0; j<3; j++) 
		{
			double ScalarH_mi_hci = Hpar - (*hc);

			if((*ms)!=0.) f += (*ms)*tanh((*ks)*ScalarH_mi_hci/(*ms));
			ks++; ms++; hc++;
		}
		//f += (Ksi[0][3])*ScalarH_mi_hc;
		f += (Ksi[0][3])*(Hpar - (*hc));
		M_Par = f*UnitEasyAxisVect;
	}
	else
	{
		M_Par = ((*(Ksi[0]))*(UnitEasyAxisVect*InstantH))*UnitEasyAxisVect;
	}

	TVector3d H_Per = InstantH - (UnitEasyAxisVect*InstantH)*UnitEasyAxisVect;
	if(DependenceIsNonlinear[1])
	{
		double AbsH_Per = sqrt(H_Per.x*H_Per.x + H_Per.y*H_Per.y + H_Per.z*H_Per.z);
		const double AbsTol = 1.E-13;
		if(AbsH_Per < AbsTol) { M_Per.x = M_Per.y = M_Per.z = 0.;}
		else
		{
			double *ks = Ksi[1], *ms = Ms[1], f = 0.;
			for(int j=0; j<3; j++) 
			{
				if((*ms)!=0.) f += (*ms)*tanh((*ks)*AbsH_Per/(*ms));
				ks++; ms++;
			}
			f += (Ksi[1][3])*AbsH_Per;
			M_Per = (f/AbsH_Per)*H_Per;
		}
	}
	else
	{
		M_Per = (*(Ksi[1]))*H_Per;
	}

	return M_Par + M_Per;
}

//-------------------------------------------------------------------------

inline double radTNonlinearAnisotropMaterial::ScalarInstantKsi(double ScalarH, char ParOrPerp, TVector3d& InstMr) // 0 - Par, 1 - Perp
{
	double f=0., InstKsi=0., ScalInstMr=0.;
	if(ParOrPerp == 0) //Par
	{
		//double hc = Hc[ParOrPerp];
		//if(hc != 0.)
		//{
		double *ks = Ksi[0], *ms = Ms[0], *hci = Hci;
		//double ScalarH_mi_hc = ScalarH - hc;
		for(int j=0; j<3; j++) 
		{
			if((*ms)!=0.) 
			{
				double ScalarH_mi_hc = ScalarH - *(hci);
				//double Sech = 1./cosh(ScalarH_mi_hc*(*ks)); //OC040213 (commented-out)
				//InstKsi += (*ks)*(*ms)*Sech*Sech;
				//f += (*ms)*tanh((*ks)*ScalarH_mi_hc/(*ms));

				double arg = (*ks)*ScalarH_mi_hc/(*ms); //OC040213
				double Sech = 1./cosh(arg);
				InstKsi += (*ks)*Sech*Sech;
				f += (*ms)*tanh(arg);
			}
			ks++; ms++; hci++;
		}
		InstKsi += *ks;
		//f += (*ks)*ScalarH_mi_hc;
		f += (*ks)*(ScalarH - *(hci));

		ScalInstMr = f - InstKsi*ScalarH;
		InstMr = ScalInstMr*UnitEasyAxisVect;
		return InstKsi;
		//}
		//else return ScalarInstantKsiAsForIsotropic(ScalarH, ParOrPerp);
	}
	else return ScalarInstantKsiAsForIsotropic(ScalarH, ParOrPerp);
}

//-------------------------------------------------------------------------

inline double radTNonlinearAnisotropMaterial::ScalarInstantKsiAsForIsotropic(double ScalarH, char ParOrPerp) // 0 - Par, 1 - Perp
{
	double *ks = Ksi[ParOrPerp], *ms = Ms[ParOrPerp];
	double f=0., InstKsi=0.;
	const double AbsTol = 1.E-13;

	double AbsInstantH = ::fabs(ScalarH);
	if(gLenArrayHB_Par <= 0)
	{
		//if(AbsInstantH < AbsTol) //OC
		if(AbsInstantH == 0.) //OC
		{
			for(int j=0; j<4; j++) InstKsi += *(ks++);
		}
		else
		{
			for(int i=0; i<3; i++)
			{
				if((*ms)!=0.) f += (*ms)*tanh((*ks)*AbsInstantH/(*ms));
				ks++; ms++;
			}
			f += (*ks)*AbsInstantH;
			InstKsi = f/AbsInstantH;
		}
	}
	else
	{
		if(ParOrPerp == 0) //Par
		{
			if(AbsInstantH == 0.) InstKsi = *gdBdH_Par;
			else InstKsi = radTNonlinearIsotropMaterial::AbsMvsAbsH_Interpol(AbsInstantH, gArrayHB_Par, gdBdH_Par, gLenArrayHB_Par)/AbsInstantH;
		}
		else //Perp
		{
			if(AbsInstantH == 0.) InstKsi = *gdBdH_Perp;
			else InstKsi = radTNonlinearIsotropMaterial::AbsMvsAbsH_Interpol(AbsInstantH, gArrayHB_Perp, gdBdH_Perp, gLenArrayHB_Perp)/AbsInstantH;
		}
	}

	return InstKsi;
}

//-------------------------------------------------------------------------

inline void radTNonlinearAnisotropMaterial::DefineInstantKsiTensor(
	const TVector3d& InstantH, TMatrix3d& InstantKsiTensor, TVector3d& InstantMr)
{
	TVector3d H_Per, DummyVect(0.,0.,0.);
	double Hpar, AbsHper, InstantKsiPar, InstantKsiPerp;

	if(DependenceIsNonlinear[0])
	{
		Hpar = UnitEasyAxisVect*InstantH;
		InstantKsiPar = ScalarInstantKsi(Hpar, 0, InstantMr);
	}
	else InstantKsiPar = *(Ksi[0]);
	if(DependenceIsNonlinear[1])
	{
		H_Per = InstantH - Hpar*UnitEasyAxisVect;
		AbsHper = sqrt(H_Per.x*H_Per.x + H_Per.y*H_Per.y + H_Per.z*H_Per.z);
		InstantKsiPerp = ScalarInstantKsi(AbsHper, 1, DummyVect);
	}
	else InstantKsiPerp = *(Ksi[1]);

	SetKsiTensor(InstantKsiPar, InstantKsiPerp, InstantKsiTensor);
}

//-------------------------------------------------------------------------

inline void radTNonlinearAnisotropMaterial::MultMatrByInstKsiAndMr(const TVector3d& InstantH, const TMatrix3d& Matr, TMatrix3d& MultByKsi, TVector3d& MultByMr)
{
	TVector3d H = InstantH;
	if(DependenceIsNonlinear[0])
	{
		if(InstantH.x==0. && InstantH.y==0. && InstantH.z==0.)
		{
			//H = Hc[0]*UnitEasyAxisVect;
			H = Hci[0]*UnitEasyAxisVect; //????
		}
	}

	TMatrix3d InstantKsiTensor;
	DefineInstantKsiTensor(H, InstantKsiTensor, RemMagn);
	MultByKsi = Matr * InstantKsiTensor; MultByMr = Matr * RemMagn;
}

//-------------------------------------------------------------------------

//inline void radTNonlinearAnisotropMaterial::FindNewH(TVector3d& H, const TMatrix3d& Matr, const TVector3d& H_Ext, double, radTg3dRelax* pMag, void* pvAuxRelax) //OC140103
inline void radTNonlinearAnisotropMaterial::FindNewH(TVector3d& H, const TMatrix3d& Matr, const TVector3d& H_Ext, double)
{
	TVector3d ESt1(1.,0.,0.), ESt2(0.,1.,0.), ESt3(0.,0.,1.);
	TMatrix3d E(ESt1, ESt2, ESt3), InstantKsiTensor, BufMatr, InvBufMatr;

	TVector3d PrevH = H; //OC140103

	DefineInstantKsiTensor(H, InstantKsiTensor, RemMagn);

	BufMatr = E - Matr*InstantKsiTensor;
	Matrix3d_inv(BufMatr, InvBufMatr);
	H = InvBufMatr*(H_Ext + Matr*RemMagn);

	//radTMaterial::SteerNewH(PrevH, H, pvAuxRelax); //OC140103
	//insert the above into most unstable materials
}

//-------------------------------------------------------------------------

inline void radTNonlinearAnisotropMaterial::Dump(std::ostream& o, int ShortSign)
{
	radTMaterial::Dump(o);
	o << "Nonlinear Anisotropic";
	if(ShortSign==1) return;

	o << endl;
	o << "   Parallel to the easy magnetization axis:" << endl;
	if(DependenceIsNonlinear[0])
	{
		//o << "      {ksi1,ms1}= {" << Ksi[0][0] << ',' << Ms[0][0] << "}" << endl;
		//o << "      {ksi2,ms2}= {" << Ksi[0][1] << ',' << Ms[0][1] << "}" << endl;
		//o << "      {ksi3,ms3}= {" << Ksi[0][2] << ',' << Ms[0][2] << "}" << endl;
		//o << "      ksi0= " << Ksi[0][3] << endl;
		//o << "      hc= " << Hc[0] << endl;

		o << "      {ksi1,ms1,hc1}= {" << Ksi[0][0] << ',' << Ms[0][0] << ',' << Hci[0] << "}" << endl;
		o << "      {ksi2,ms2,hc2}= {" << Ksi[0][1] << ',' << Ms[0][1] << ',' << Hci[1] << "}" << endl;
		o << "      {ksi3,ms3,hc3}= {" << Ksi[0][2] << ',' << Ms[0][2] << ',' << Hci[2] << "}" << endl;
		o << "      {ksi0,hc0}= {" << Ksi[0][3] << ',' << Hci[3] << "}" << endl;
	}
	else
	{
		o << "      ksi= " << Ksi[0][0] << endl;
	}

	o << "   Perpendicular to the easy magnetization axis:" << endl;
	if(DependenceIsNonlinear[1])
	{
		o << "      {ksi1,ms1}= {" << Ksi[1][0] << ',' << Ms[1][0] << "}" << endl;
		o << "      {ksi2,ms2}= {" << Ksi[1][1] << ',' << Ms[1][1] << "}" << endl;
		o << "      {ksi3,ms3}= {" << Ksi[1][2] << ',' << Ms[1][2] << "}" << endl;
		o << "      ksi0= " << Ksi[1][3] << endl;
		//o << "      hc= " << Hc[1] << endl;
	}
	else
	{
		o << "      ksi= " << Ksi[1][0] << endl;
	}

	o << "   Easy magnetization axis: ";
	if(EasyAxisDefined) o << "{" << UnitEasyAxisVect.x << "," << UnitEasyAxisVect.y << "," << UnitEasyAxisVect.z << "}" << endl;
	else o << "not defined" << endl;

	o << "   Memory occupied: " << SizeOfThis() << " bytes";
}

//-------------------------------------------------------------------------

inline void radTNonlinearAnisotropMaterial::DumpBin_NonlinearAnisotropMaterial(CAuxBinStrVect& oStr)
{
	//double Ksi[2][4];
	int i, j;
	for(i=0; i<2; i++)
		for(j=0; j<3; j++)
			oStr << Ksi[i][j];

	//double Ms[2][3];
	for(i=0; i<2; i++)
		for(j=0; j<3; j++)
			oStr << Ms[i][j];

	//double Hci[4];
	for(i=0; i<4; i++) oStr << Hci[i];

	//char DependenceIsNonlinear[2];
	oStr << DependenceIsNonlinear[0] << DependenceIsNonlinear[1];

	//int gLenArrayHB_Par
	oStr << gLenArrayHB_Par;
	//TVector2d *gArrayHB_Par
	if((gLenArrayHB_Par) && (gArrayHB_Par != 0))
	{
		oStr << (char)1;
		TVector2d *t_gArrayHB_Par = gArrayHB_Par;
		for(int i=0; i<gLenArrayHB_Par; i++) oStr << (*(t_gArrayHB_Par++));
	}
	else oStr << (char)0;
	//double *gdBdH_Par
	if((gLenArrayHB_Par) && (gdBdH_Par != 0))
	{
		oStr << (char)1;
		double *t_gdBdH_Par = gdBdH_Par;
		for(int i=0; i<gLenArrayHB_Par; i++) oStr << (*(t_gdBdH_Par++));
	}
	else oStr << (char)0;

	//int gLenArrayHB_Perp;
	oStr << gLenArrayHB_Par;
	//TVector2d *gArrayHB_Perp;
	if((gLenArrayHB_Perp) && (gArrayHB_Perp != 0))
	{
		oStr << (char)1;
		TVector2d *t_gArrayHB_Perp = gArrayHB_Perp;
		for(int i=0; i<gLenArrayHB_Perp; i++) oStr << (*(t_gArrayHB_Perp++));
	}
	else oStr << (char)0;
	//double *gdBdH_Perp;
	if((gLenArrayHB_Perp) && (gdBdH_Perp != 0))
	{
		oStr << (char)1;
		double *t_gdBdH_Perp = gdBdH_Perp;
		for(int i=0; i<gLenArrayHB_Perp; i++) oStr << (*(t_gdBdH_Perp++));
	}
	else oStr << (char)0;

	//double gMaxKsi_Par, gMaxKsi_Perp;
	oStr << gMaxKsi_Par << gMaxKsi_Perp;

	//TVector3d UnitEasyAxisVect;
	oStr << UnitEasyAxisVect;
}

//-------------------------------------------------------------------------

inline void radTNonlinearAnisotropMaterial::DumpBinParse_NonlinearAnisotropMaterial(CAuxBinStrVect& inStr)
{
	//double Ksi[2][4];
	int i, j;
	for(i=0; i<2; i++)
		for(j=0; j<3; j++)
			inStr >> Ksi[i][j];

	//double Ms[2][3];
	for(i=0; i<2; i++)
		for(j=0; j<3; j++)
			inStr >> Ms[i][j];

	//double Hci[4];
	for(i=0; i<4; i++) inStr >> Hci[i];

	//char DependenceIsNonlinear[2];
	inStr >> DependenceIsNonlinear[0];
	inStr >> DependenceIsNonlinear[1];

	//int gLenArrayHB_Par
	inStr >> gLenArrayHB_Par;
	char cTest = 0;
	inStr >> cTest;
	//TVector2d *gArrayHB_Par
	if(cTest > 0)
	{
		vgArrayHB_Par.resize(gLenArrayHB_Par);
		gArrayHB_Par = vgArrayHB_Par.data();
		TVector2d *t_gArrayHB_Par = gArrayHB_Par;
		for(int i=0; i<gLenArrayHB_Par; i++) inStr >> (*(t_gArrayHB_Par++));
	}
	//double *gdBdH_Par
	inStr >> cTest;
	if(cTest > 0)
	{
		vgdBdH_Par.resize(gLenArrayHB_Par);
		gdBdH_Par = vgdBdH_Par.data();
		double *t_gdBdH_Par = gdBdH_Par;
		for(int i=0; i<gLenArrayHB_Par; i++) inStr >> (*(t_gdBdH_Par++));
	}

	//int gLenArrayHB_Perp;
	inStr >> gLenArrayHB_Perp;
	//TVector2d *gArrayHB_Perp;
	inStr >> cTest;
	if(cTest > 0)
	{
		vgArrayHB_Perp.resize(gLenArrayHB_Perp);
		gArrayHB_Perp = vgArrayHB_Perp.data();
		TVector2d *t_gArrayHB_Perp = gArrayHB_Perp;
		for(int i=0; i<gLenArrayHB_Perp; i++) inStr >> (*(t_gArrayHB_Perp++));
	}
	//double *gdBdH_Perp;
	inStr >> cTest;
	if(cTest > 0)
	{
		vgdBdH_Perp.resize(gLenArrayHB_Perp);
		gdBdH_Perp = vgdBdH_Perp.data();
		double *t_gdBdH_Perp = gdBdH_Perp;
		for(int i=0; i<gLenArrayHB_Perp; i++) inStr >> (*(t_gdBdH_Perp++));
	}

	//double gMaxKsi_Par, gMaxKsi_Perp;
	inStr >> gMaxKsi_Par;
	inStr >> gMaxKsi_Perp;

	//TVector3d UnitEasyAxisVect;
	inStr >> UnitEasyAxisVect;
}

//-------------------------------------------------------------------------

inline void radTNonlinearAnisotropMaterial::DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
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

	//Members of radTNonlinearAnisotropMaterial
	DumpBin_NonlinearAnisotropMaterial(oStr);
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTNonlinearLaminatedMaterial : public radTNonlinearAnisotropMaterial {
	
	double gPackFactor;

public:

	radTNonlinearLaminatedMaterial(double* InMs, double* InKsi, int lenMs, double InPackFactor, double* pN) : radTNonlinearAnisotropMaterial()
	{// This sets up isotropic data only. The anisotropic data is derived from isotropic ones at run time.
	 // The isotropic data are stored in Ksi[0], Ms[0] arrays of radTNonlinearAnisotropMaterial

		if((InKsi == 0) || (InMs == 0) || (lenMs <= 0) || (lenMs > 3) || (InPackFactor <= 0) || (InPackFactor > 1)) throw 0;

		if(lenMs > 3) lenMs = 3;
		double *CurInKsi = InKsi, *CurInMs = InMs;
		double *CurKsi = Ksi[0], *CurMs = Ms[0];

		CurKsi[0] = 0; CurKsi[1] = 0; CurKsi[2] = 0; CurKsi[3] = 0;
		for(int i=0; i<3; i++)
		{
			if(i < lenMs) { *(CurKsi++) = *(CurInKsi++); *(CurMs++) = *(CurInMs++);}
			else { *(CurKsi++) = 0.; *(CurMs++) = 0.;}
		}
		DependenceIsNonlinear[0] = 1;

		gPackFactor = InPackFactor;
		SetupEasyAxisIfPossible(pN);
	}

	radTNonlinearLaminatedMaterial(TVector2d* InArrayHB, int InLenArrayHB, double InPackFactor, double* pN)
	{// This sets up isotropic data only. The anisotropic data is derived from isotropic ones at run time.
	 // The isotropic data are stored in .._Par arrays of radTNonlinearAnisotropMaterial

		if((InArrayHB == 0) || (InLenArrayHB <= 0) || (InPackFactor <= 0) || (InPackFactor > 1)) throw 0;

		gLenArrayHB_Par = gLenArrayHB_Perp = 0;
		gArrayHB_Par = gArrayHB_Perp = 0;
		gdBdH_Par = gdBdH_Perp = 0;
		gMaxKsi_Par = gMaxKsi_Perp = 0;

		double ZeroTol = 1e-10;

		char PrependZero = 0;
		if((InArrayHB->x > ZeroTol) && (InArrayHB->y > ZeroTol))
		{
			InLenArrayHB++; PrependZero = 1;
		}

		gLenArrayHB_Par = InLenArrayHB;
		AllocateArrays(gLenArrayHB_Par, 0);
		radTNonlinearIsotropMaterial::CopyArrayHB(gArrayHB_Par, InArrayHB, InLenArrayHB, PrependZero);
		radTNonlinearIsotropMaterial::Compute_dBdH(gArrayHB_Par, gdBdH_Par, gLenArrayHB_Par, gMaxKsi_Par);

		gPackFactor = InPackFactor;
		SetupEasyAxisIfPossible(pN);
	}

	radTNonlinearLaminatedMaterial(CAuxBinStrVect& inStr)
	{//Instantiates from string according to DumpBin
		DumpBinParse_Material(inStr);
		//Members of radTNonlinearAnisotropMaterial
		DumpBinParse_NonlinearAnisotropMaterial(inStr);

		//double gPackFactor;
		inStr >> gPackFactor;
	}

	radTNonlinearLaminatedMaterial() : radTNonlinearAnisotropMaterial() {}

	//int Type_Material() { return 4;} //same as radTNonlinearAnisotropMaterial
	//double ScalarInstantKsiAsForIsotropic(double, char); //to correct in radTNonlinearAnisotropMaterial
	//double ScalarInstantKsi(double, char, TVector3d&); //don't need
	int Type_NonlinearAnisotropMaterial() { return 1;}

	//void FindNewH(TVector3d& H, const TMatrix3d& Matr, const TVector3d& H_Ext, double, radTg3dRelax* pMag, void* pvAuxRelax) //OC140103
	void FindNewH(TVector3d& H, const TMatrix3d& Matr, const TVector3d& H_Ext, double)
	{
		TVector3d PrevH = H; //OC140103

		TVector3d ESt1(1.,0.,0.), ESt2(0.,1.,0.), ESt3(0.,0.,1.);
		TMatrix3d E(ESt1, ESt2, ESt3), InstantKsiTensor, BufMatr, InvBufMatr;

		const int AmOfExtraCycles = 8; //to tune //OC140103

		for(int i=0; i<AmOfExtraCycles; i++) //OC140103
		{
			DefineInstantKsiTensor(H, InstantKsiTensor, RemMagn);
			BufMatr = E - Matr*InstantKsiTensor;
			Matrix3d_inv(BufMatr, InvBufMatr);
			H = InvBufMatr*(H_Ext + Matr*RemMagn);
		}

		//radTMaterial::SteerNewH(PrevH, H, pvAuxRelax); //OC140103
		//insert the above into most unstable materials
	}

	TVector3d M(const TVector3d& InstantH)
	{
		TVector3d M_Norm(0,0,0), M_Tang(0,0,0);
		const double AbsTol = 1.E-13;

		double AbsH_Norm = UnitEasyAxisVect*InstantH;
		TVector3d H_Norm = AbsH_Norm*UnitEasyAxisVect;
		TVector3d H_Tang = InstantH - H_Norm;

		//double AbsH_Tang = sqrt(H_Tang.x*H_Tang.x + H_Tang.y*H_Tang.y + H_Tang.z*H_Tang.z);

		//char AbsH_NormIsZero = (AbsH_Norm < AbsTol);
		//char AbsH_TangIsZero = (AbsH_Tang < AbsTol);
		//if(AbsH_NormIsZero && AbsH_TangIsZero) { return M_Norm;} //zero

		//double Ksi_Tang = 0, Ksi_Norm = 0;
		//if(!AbsH_TangIsZero)
		//{
		//  Ksi_Tang = gPackFactor*ScalarInstantKsiAsForIsotropic(AbsH_Tang, 0); // 0 - Par, 1 - Perp
		//	M_Tang = Ksi_Tang*H_Tang;
		//}
		//if(!AbsH_NormIsZero)
		//{
		//  double KsiIso_Norm = ScalarInstantKsiAsForIsotropic(AbsH_Norm, 0); // 0 - Par, 1 - Perp
		//  Ksi_Norm = gPackFactor*KsiIso_Norm/(1 + (1 - gPackFactor)*KsiIso_Norm);
		//	M_Norm = Ksi_Norm*H_Norm;
		//}

		double AbsH = sqrt(InstantH.x*InstantH.x + InstantH.y*InstantH.y + InstantH.z*InstantH.z);
		if(AbsH < AbsTol) { return M_Norm;} //zero
		double KsiIso = ScalarInstantKsiAsForIsotropic(AbsH, 0); // 0 - Par, 1 - Perp

		double Ksi_Tang = gPackFactor*KsiIso;
		M_Tang = Ksi_Tang*H_Tang;
		double Ksi_Norm = gPackFactor*KsiIso/(1. + (1. - gPackFactor)*KsiIso);
		M_Norm = Ksi_Norm*H_Norm;

		return (M_Norm + M_Tang);
	}

	void DefineInstantKsiTensor(const TVector3d& InstantH, TMatrix3d& InstantKsiTensor, TVector3d& InstantMr)
	{
		InstantMr.x = 0; InstantMr.y = 0; InstantMr.z = 0;

		//TVector3d M_Norm(0,0,0), M_Tang(0,0,0);
		//const double AbsTol = 1.E-13;

		//double AbsH_Norm = UnitEasyAxisVect*InstantH;
		//TVector3d H_Norm = AbsH_Norm*UnitEasyAxisVect;
		//TVector3d H_Tang = InstantH - H_Norm;
		//double AbsH_Tang = sqrt(H_Tang.x*H_Tang.x + H_Tang.y*H_Tang.y + H_Tang.z*H_Tang.z);

		//double Ksi_Tang = gPackFactor*ScalarInstantKsiAsForIsotropic(AbsH_Tang, 0); // 0 - Par, 1 - Perp
		//double KsiIso_Norm = ScalarInstantKsiAsForIsotropic(AbsH_Norm, 0); // 0 - Par, 1 - Perp
		//double Ksi_Norm = gPackFactor*KsiIso_Norm/(1 + (1 - gPackFactor)*KsiIso_Norm);

		double AbsH = sqrt(InstantH.x*InstantH.x + InstantH.y*InstantH.y + InstantH.z*InstantH.z);
		//if(AbsH < AbsTol) { return M_Norm;} //zero
		double KsiIso = ScalarInstantKsiAsForIsotropic(AbsH, 0); // 0 - Par, 1 - Perp

		double Ksi_Tang = gPackFactor*KsiIso;
		//M_Tang = Ksi_Tang*H_Tang;
		double Ksi_Norm = gPackFactor*KsiIso/(1. + (1. - gPackFactor)*KsiIso);
		//M_Norm = Ksi_Norm*H_Norm;

		SetKsiTensor(Ksi_Norm, Ksi_Tang, InstantKsiTensor);
	}

	void MultMatrByInstKsiAndMr(const TVector3d& InstantH, const TMatrix3d& Matr, TMatrix3d& MultByKsi, TVector3d& MultByMr)
	{
		TVector3d H = InstantH;

		TMatrix3d InstantKsiTensor;
		DefineInstantKsiTensor(H, InstantKsiTensor, RemMagn);
		MultByKsi = Matr * InstantKsiTensor; 
		MultByMr.x = 0; MultByMr.y = 0; MultByMr.z = 0; 
	}

	int DuplicateItself(radThg& hg, radTApplication*, char) 
	{
		return FinishDuplication(new radTNonlinearLaminatedMaterial(*this), hg);
	}
	int SizeOfThis() { return sizeof(radTNonlinearLaminatedMaterial);}

	void Dump(std::ostream& o, int ShortSign)
	{
		radTMaterial::Dump(o);
		o << "Nonlinear Laminated";
		if(ShortSign==1) return;

		o << endl;
		o << "   Packing factor: " << gPackFactor << endl;

		o << "   Isotropic parameters:" << endl;

		if((gArrayHB_Par == 0) || (gLenArrayHB_Par == 0))
		{
			o << "      {ksi1,ms1}= {" << Ksi[0][0] << ',' << Ms[0][0] << "}" << endl;
			o << "      {ksi2,ms2}= {" << Ksi[0][1] << ',' << Ms[0][1] << "}" << endl;
			o << "      {ksi3,ms3}= {" << Ksi[0][2] << ',' << Ms[0][2] << "}" << endl;
		}
		else o << "      M(H) defined by table of values" << endl;

		o << "   Vector normal to the lamination planes: ";
		if(EasyAxisDefined) o << "{" << UnitEasyAxisVect.x << "," << UnitEasyAxisVect.y << "," << UnitEasyAxisVect.z << "}" << endl;
		else o << "not defined" << endl;
		
		o << "   Memory occupied: " << SizeOfThis() << " bytes";
	}

	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, radTmhg& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
	{
		vElemKeysOut.push_back(elemKey);
		oStr << elemKey;

		//Next 5 bytes define/encode element type:
		oStr << (char)Type_g();
		oStr << (char)Type_Material();
		oStr << (char)Type_NonlinearAnisotropMaterial();
		oStr << (char)0;
		oStr << (char)0;

		//Members of radTMaterial
		DumpBin_Material(oStr);
	
		//Members of radTNonlinearAnisotropMaterial
		DumpBin_NonlinearAnisotropMaterial(oStr);

		//double gPackFactor;
		oStr << gPackFactor;
	}
};

//-------------------------------------------------------------------------

#endif
