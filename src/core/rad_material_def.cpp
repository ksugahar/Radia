/*-------------------------------------------------------------------------
*
* File name:      radmater.cpp
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

#include "rad_material_def.h"
#include "rad_material_aux.h"
#include "rad_geometry_3d.h"
#include "auxparse.h"
#include "rad_relaxation_methods.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
/**
void radTMaterial::SteerNewH(TVector3d& PrevH, TVector3d& InstantH, void* pvAuxRelax) //OC140103
{// may modify InstantH
	const float MinRelaxPar = (float)0.1; //(float)0.0001; //(float)0.24; //to tune
	const float RelaxParReduceCoef = (float)0.9; //to tune
	const double RelDifE2 = 0.99; //1E-05; //to tune
	const int AmOfBadPassesToAllow = 1; //1; //to tune

	if(pvAuxRelax == 0) return;

	//double AbsInstantH_E2 = InstantH.x*InstantH.x + InstantH.y*InstantH.y + InstantH.z*InstantH.z;
	double dx = InstantH.x - PrevH.x, dy = InstantH.y - PrevH.y, dz = InstantH.z - PrevH.z;
	double AbsDeltH_E2 = dx*dx + dy*dy + dz*dz;

	radTRelaxAuxData *pCurRelaxAuxData = (radTRelaxAuxData*)pvAuxRelax;

	float &PrevAbsDeltH_E2 = pCurRelaxAuxData->AbsDeltH;
	float &RelaxPar = pCurRelaxAuxData->RelaxPar;
	int &BadPassCount = pCurRelaxAuxData->BadPassCounts;

	//if(AbsDeltH_E2 >= PrevAbsDeltH_E2 - AbsInstantH_E2*RelMinDifE2)
	if(AbsDeltH_E2 >= RelDifE2*PrevAbsDeltH_E2)
	{
		if(BadPassCount >= AmOfBadPassesToAllow)
		{
			RelaxPar *= RelaxParReduceCoef;
			if(RelaxPar < MinRelaxPar) RelaxPar = MinRelaxPar;
			BadPassCount = 0;
		}
		else 
		{
			BadPassCount++;
		}
	}
	else
	{
		BadPassCount = 0;
	}
	PrevAbsDeltH_E2 = (float)AbsDeltH_E2;
	InstantH = RelaxPar*InstantH + (1. - RelaxPar)*PrevH;
}
**/
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//void radTNonlinearIsotropMaterial::FindNewH(TVector3d& InstantH, const TMatrix3d& Matr, const TVector3d& H_Ext, double DesiredPrecOnMagnE2, radTg3dRelax* pMag, void* p) //OC140103
void radTNonlinearIsotropMaterial::FindNewH(TVector3d& InstantH, const TMatrix3d& Matr, const TVector3d& H_Ext, double DesiredPrecOnMagnE2) //OC140103
{
	static const double mu_0 = 4.0e-7 * 3.14159265358979323846;  // 4*pi*1e-7
	const double AbsHZeroTol = 1.E-10;
	const double AbsMZeroTol = 1.E-10;
	const int MaxIterToFindH = 15; //50; //15;

	TMatrix3d E; E.Str0.x = 1.; E.Str1.y = 1.; E.Str2.z = 1.;
	TMatrix3d BufMatr, InvBufMatr;

	double MisfitM = 3;
	double AbsInstantH = sqrt(InstantH.x*InstantH.x + InstantH.y*InstantH.y + InstantH.z*InstantH.z);
	double f, InstKsi=0.;

	for(int i=0; i<MaxIterToFindH; i++)
	{
		f=0.;
		if(gLenArrayHB == 0)
		{
			// Analytical formula (tanh model)
			if(AbsInstantH <= AbsHZeroTol)
			{
				for(int j=0; j<lenMs_ks; j++) InstKsi += ks[j];
				f = 0;
			}
			else
			{
				for(int i=0; i<lenMs_ks; i++) if(Ms[i]!=0.) f += Ms[i]*tanh(ks[i]*AbsInstantH/Ms[i]);
				InstKsi = f/AbsInstantH;
			}
		}
		else
		{
			// ELF-compatible: B-H curve interpolation
			if(AbsInstantH <= AbsHZeroTol)
			{
				// ELF method: At H=0, use 2nd B-H point for initial chi
				// chi = B[1] / (mu_0 * H[1]) - 1 (mucal0 style)
				if(gLenArrayHB >= 2 && gArrayHB[1].x > 1.0e-15)
				{
					InstKsi = gArrayHB[1].y / (mu_0 * gArrayHB[1].x) - 1.0;
				}
				else
				{
					InstKsi = gdBdH[0] / mu_0 - 1.0;  // fallback
				}
				f = 0;
			}
			else
			{
				// AbsMvsAbsH_Interpol now returns M = B/mu_0 - H from B-H curve
				f = AbsMvsAbsH_Interpol(AbsInstantH, gArrayHB, gdBdH, gLenArrayHB);
				InstKsi = f/AbsInstantH;
			}
		}

		BufMatr = E - InstKsi*Matr;
		Matrix3d_inv(BufMatr, InvBufMatr);
		InstantH = InvBufMatr*H_Ext;

		double NewAbsInstantH = sqrt(InstantH.x*InstantH.x + InstantH.y*InstantH.y + InstantH.z*InstantH.z);
		double NewMisfitM = f - NewAbsInstantH*InstKsi;

		double AbsMisfitM = ::fabs(MisfitM), AbsNewMisfitM = ::fabs(NewMisfitM);

		double ProbNew = AbsMisfitM;
		if(NewMisfitM*MisfitM > 0) ProbNew += 0.5*AbsNewMisfitM;
		double ProbOld = AbsNewMisfitM;
		double Alpha = ProbNew/(ProbNew + ProbOld);
		AbsInstantH = Alpha*NewAbsInstantH + (1 - Alpha)*AbsInstantH;

		MisfitM = f - AbsInstantH*InstKsi;

/**
		TVector3d CurM = M(InstantH);
		f = CurM.Abs();

		TMatrix3d InstantKsiTensor;
		TVector3d InstMr, PrevInstantH = InstantH;
		DefineInstantKsiTensor(InstantH, InstantKsiTensor, InstMr);
		BufMatr = E - Matr*InstantKsiTensor;
		Matrix3d_inv(BufMatr, InvBufMatr);
		InstantH = InvBufMatr*(H_Ext + Matr*InstMr);
		TVector3d NewM = InstantKsiTensor*InstantH + InstMr;
		double NewMisfitM = f - NewM.Abs();

		double NewAbsInstantH = sqrt(InstantH.x*InstantH.x + InstantH.y*InstantH.y + InstantH.z*InstantH.z);
		double AbsMisfitM = ::fabs(MisfitM), AbsNewMisfitM = ::fabs(NewMisfitM);

		if(AbsNewMisfitM < AbsMZeroTol) break;

		double ProbNew = AbsMisfitM;
		if(NewMisfitM*MisfitM > 0) ProbNew += 0.5*AbsNewMisfitM;
		double ProbOld = AbsNewMisfitM;
		double Alpha = ProbNew/(ProbNew + ProbOld);
		InstantH = Alpha*InstantH + (1. - Alpha)*PrevInstantH;
		AbsInstantH = InstantH.Abs();

		NewM = InstantKsiTensor*InstantH + InstMr;
		MisfitM = f - NewM.Abs();
**/

		if(MisfitM*MisfitM <= DesiredPrecOnMagnE2) break;
	}
}

//-------------------------------------------------------------------------

double radTNonlinearIsotropMaterial::FuncNewAbsH(double AbsH, const TMatrix3d& Matr, const TVector3d& H_Ext)
{
	static const double mu_0 = 4.0e-7 * 3.14159265358979323846;  // 4*pi*1e-7
	const double AbsHZeroTol = 1.E-10;
	double f=0., InstKsi=0.;

	if(gLenArrayHB == 0)
	{
		// Analytical formula (tanh model)
		if(AbsH <= AbsHZeroTol)
		{
			for(int j=0; j<lenMs_ks; j++) InstKsi += ks[j];
			f = 0;
		}
		else
		{
			for(int i=0; i<lenMs_ks; i++) if(Ms[i]!=0.) f += Ms[i]*tanh(ks[i]*AbsH/Ms[i]);
			InstKsi = f/AbsH;
		}
	}
	else
	{
		// ELF-compatible: B-H curve interpolation
		if(AbsH <= AbsHZeroTol)
		{
			// ELF method: At H=0, use 2nd B-H point for initial chi
			// chi = B[1] / (mu_0 * H[1]) - 1 (mucal0 style)
			if(gLenArrayHB >= 2 && gArrayHB[1].x > 1.0e-15)
			{
				InstKsi = gArrayHB[1].y / (mu_0 * gArrayHB[1].x) - 1.0;
			}
			else
			{
				InstKsi = gdBdH[0] / mu_0 - 1.0;  // fallback
			}
			f = 0;
		}
		else
		{
			// AbsMvsAbsH_Interpol now returns M = B/mu_0 - H from B-H curve
			f = AbsMvsAbsH_Interpol(AbsH, gArrayHB, gdBdH, gLenArrayHB);
			InstKsi = f/AbsH;
		}
	}
	TMatrix3d E; E.Str0.x = 1.; E.Str1.y = 1.; E.Str2.z = 1.;
	TMatrix3d BufMatr, InvBufMatr;
	BufMatr = E - InstKsi*Matr;
	Matrix3d_inv(BufMatr, InvBufMatr);
	TVector3d AuxVectInstH = InvBufMatr*H_Ext;
	return AuxVectInstH.Abs();
}

//-------------------------------------------------------------------------

void radTNonlinearIsotropMaterial::DefineInstantKsiTensor(const TVector3d& InstantH, TMatrix3d& InstKsi, TVector3d& InstMr)
{
	double H0 = sqrt(InstantH.x*InstantH.x + InstantH.y*InstantH.y + InstantH.z*InstantH.z);

	double f=0, dfdH=0;
	DefineScalarM_dMdH(H0, f, dfdH);

	double AbsZeroTolH = 1.e-10;

	if(::fabs(H0) < AbsZeroTolH)
	{
		InstKsi.Str0.x = dfdH; InstKsi.Str0.y = 0; InstKsi.Str0.z = 0; 
		InstKsi.Str1.x = 0; InstKsi.Str1.y = dfdH; InstKsi.Str1.z = 0; 
		InstKsi.Str2.x = 0; InstKsi.Str2.y = 0; InstKsi.Str2.z = dfdH;
		InstMr = RemMagn;
		return;
	}

	double H0e3 = H0*H0*H0;
	double InvH0e3 = 1./H0e3;
	double Hx0 = InstantH.x, Hy0 = InstantH.y, Hz0 = InstantH.z;
	double Hx0e2 = Hx0*Hx0, Hy0e2 = Hy0*Hy0, Hz0e2 = Hz0*Hz0;

	double H0_dfdH_mi_f = H0*dfdH - f;
	double BufOffDiag = InvH0e3*H0_dfdH_mi_f;
	double H0_dfdH = H0*dfdH;
	double f_mi_H0_dfdH_d_H0 = -H0_dfdH_mi_f/H0;

	InstKsi.Str0.x = InvH0e3*((Hy0e2 + Hz0e2)*f + Hx0e2*H0_dfdH);
	InstKsi.Str0.y = BufOffDiag*Hx0*Hy0; //InvH0e3*Hx0*Hy0*(H0*dfdH - f);
	InstKsi.Str0.z = BufOffDiag*Hx0*Hz0; //InvH0e3*Hx0*Hz0*(H0*dfdH - f);
	InstKsi.Str1.x = InstKsi.Str0.y;
	InstKsi.Str1.y = InvH0e3*((Hx0e2 + Hz0e2)*f + Hy0e2*H0_dfdH);
	InstKsi.Str1.z = BufOffDiag*Hy0*Hz0; //InvH0e3*Hy0*Hz0*(H0*dfdH - f);
	InstKsi.Str2.x = InstKsi.Str0.z;
	InstKsi.Str2.y = InstKsi.Str1.z;
	InstKsi.Str2.z = InvH0e3*((Hx0e2 + Hy0e2)*f + Hz0e2*H0_dfdH);

	InstMr.x = Hx0*f_mi_H0_dfdH_d_H0 + RemMagn.x;
	InstMr.y = Hy0*f_mi_H0_dfdH_d_H0 + RemMagn.y;
	InstMr.z = Hz0*f_mi_H0_dfdH_d_H0 + RemMagn.z;

/**
	double Der, f, dInstKsi;
	Der = f = dInstKsi = 0.;

	if(gLenArrayHB == 0)
	{
		if(H0 == 0.)
		{
			for(int j=0; j<lenMs_ks; j++) Der += ks[j];
			dInstKsi = Der;
		}
		else
		{
			for(int i=0; i<lenMs_ks; i++) if(Ms[i]!=0.) f += Ms[i]*tanh(ks[i]*H0/Ms[i]);
			dInstKsi = f/H0;
		}
	}
	else
	{
		if(H0 == 0.) dInstKsi = *gdMdH;
		else dInstKsi = AbsMvsAbsH_Interpol(H0, gArrayHB, gdMdH, gLenArrayHB)/H0;
	}

	TVector3d Ksi_Str0(dInstKsi,0.,0.), Ksi_Str1(0.,dInstKsi,0.), Ksi_Str2(0.,0.,dInstKsi);
	InstKsi.Str0 = Ksi_Str0; InstKsi.Str1 = Ksi_Str1; InstKsi.Str2 = Ksi_Str2;
	InstMr = RemMagn;
**/
}

//-------------------------------------------------------------------------

void radTNonlinearIsotropMaterial::DefineScalarM_dMdH(double AbsH, double& f, double& dfdH)
{
	static const double mu_0 = 4.0e-7 * 3.14159265358979323846;  // 4*pi*1e-7
	double AbsZeroTolH = 1.e-10;
	if(AbsH < 0.) AbsH = 0.;

	f = dfdH = 0.;
	if(gLenArrayHB == 0)
	{
		// Analytical formula (tanh model)
		for(int i=0; i<lenMs_ks; i++)
		{
			double ms_i = Ms[i], ks_i = ks[i];
			if(ms_i == 0.) continue;

			double arg = ks_i*AbsH/ms_i;
			f += ms_i*tanh(arg);

			double cosh_arg = cosh(arg);
			dfdH += ks_i/(cosh_arg*cosh_arg);
		}
	}
	else
	{
		// ELF-compatible: B-H curve interpolation
		if(AbsH < AbsZeroTolH)
		{
			f = 0;
			// ELF method: use 2nd B-H point for initial chi
			if(gLenArrayHB >= 2 && gArrayHB[1].x > 1.0e-15)
			{
				dfdH = gArrayHB[1].y / (mu_0 * gArrayHB[1].x) - 1.0;
			}
			else
			{
				dfdH = gdBdH[0] / mu_0 - 1.0;
			}
			return;
		}
		AbsMvsAbsH_FuncAndDer_Interpol(AbsH, gArrayHB, gdBdH, gLenArrayHB, f, dfdH);
	}
}

//-------------------------------------------------------------------------

void radTNonlinearIsotropMaterial::DefineScalarM(double AbsInstantH, double& f, double& InstKsi)
{
	static const double mu_0 = 4.0e-7 * 3.14159265358979323846;  // 4*pi*1e-7
	const double AbsHZeroTol = 1.E-09;
	if(gLenArrayHB == 0)
	{
		// Analytical formula (tanh model)
		if(AbsInstantH <= AbsHZeroTol)
		{
			for(int j=0; j<lenMs_ks; j++) InstKsi += ks[j];
			f = 0;
		}
		else
		{
			for(int i=0; i<lenMs_ks; i++) if(Ms[i]!=0.) f += Ms[i]*tanh(ks[i]*AbsInstantH/Ms[i]);
			InstKsi = f/AbsInstantH;
		}
	}
	else
	{
		// ELF-compatible: B-H curve interpolation
		if(AbsInstantH <= AbsHZeroTol)
		{
			// ELF method: At H=0, use 2nd B-H point for initial chi
			// chi = B[1] / (mu_0 * H[1]) - 1 (mucal0 style)
			if(gLenArrayHB >= 2 && gArrayHB[1].x > 1.0e-15)
			{
				InstKsi = gArrayHB[1].y / (mu_0 * gArrayHB[1].x) - 1.0;
			}
			else
			{
				InstKsi = gdBdH[0] / mu_0 - 1.0;  // fallback
			}
			f = 0;
		}
		else
		{
			// AbsMvsAbsH_Interpol now returns M = B/mu_0 - H from B-H curve
			f = AbsMvsAbsH_Interpol(AbsInstantH, gArrayHB, gdBdH, gLenArrayHB);
			InstKsi = f/AbsInstantH;
		}
	}
}

//-------------------------------------------------------------------------

double radTNonlinearIsotropMaterial::Derivative5(TVector2d* f, int PoIndx)
{
	double x0 = f[0].x, x1 = f[1].x, x2 = f[2].x, x3 = f[3].x, x4 = f[4].x;
	x1 -= x0; x2 -= x0; x3 -= x0; x4 -= x0; 
	double f0 = f[0].y, f1 = f[1].y, f2 = f[2].y, f3 = f[3].y, f4 = f[4].y;

	char IsIncreasing = ((x1 > 0.) && (x2 > x1) && (x3 > x2) && (x4 > x3)) && 
						((f0 <= f1) && (f1 <= f2) && (f2 <= f3) && (f3 <= f4));

	double x1mix2 = (x1-x2), x1mix3 = (x1-x3), x1mix4 = (x1-x4), x2mix3 = (x2-x3), x2mix4 = (x2-x4), x3mix4 = (x3-x4);
	double x1e2 = x1*x1, x2e2 = x2*x2, x3e2 = x3*x3, x4e2 = x4*x4;
	if(PoIndx==0)
	{
		double x1mix3e2 = x1mix3*x1mix3;
		double Der = (f4*x1e2*x1mix2*x2e2*x1mix3*x2mix3*x3e2 + 
					(-f3*x1e2*x1mix2*x2e2*x1mix4*x2mix4 + 
					x3e2*(f2*x1e2*x1mix3*x1mix4 - f1*x2e2*x2mix3*x2mix4)*
					x3mix4)*x4e2 - f0*x1mix2*x1mix3*x2mix3*x1mix4*
					x2mix4*x3mix4*(x1*x2*x3 + x2*x3*x4 + x1*(x2 + x3)*x4))/
					(x1*x1mix2*x2*x1mix3*x2mix3*x3*x1mix4*x2mix4*x3mix4*x4);
		if(IsIncreasing && (Der < 0.)) Der = 0.;
		return Der;
	}
	else if(PoIndx==1)
	{
		double x1e3 = x1e2*x1;
		double x1mix2e2 = x1mix2*x1mix2, x1mix3e2 = x1mix3*x1mix3, x1mix4e2 = x1mix4*x1mix4;
		double x1mix3e4 = x1mix3e2*x1mix3e2;
		double Der = (-f4*x1e2*x1mix2e2*x2*x1mix3e2*x2mix3*x3 + 
					f0*x1mix2e2*x1mix3e2*x2mix3*x1mix4e2*x2mix4*x3mix4 + 
					x4*(f3*x1e2*x1mix2e2*x2*x1mix4e2*x2mix4 - 
					x3*x3mix4*(f2*x1e2*x1mix3e2*x1mix4e2 + f1*x2*x2mix3*x2mix4*
					(-4*x1e3 + x2*x3*x4 + 3*x1e2*(x2 + x3 + x4) - 
					2*x1*(x3*x4 + x2*(x3 + x4))))))/
					(x1*x1mix2*x2*x1mix3*x2mix3*x3*x1mix4*x2mix4*x3mix4*x4);
		if(IsIncreasing && (Der < 0.)) Der = 0.;
		return Der;
	}
	else if(PoIndx==2)
	{
		double x1mix2e2 = x1mix2*x1mix2, x2mix3e2 = x2mix3*x2mix3, x1mix3e2 = x1mix3*x1mix3, x2mix4e2 = x2mix4*x2mix4;
		double Der = (f4*x1*x1mix2e2*x2e2*x1mix3*x2mix3e2*x3 - 
					f0*x1mix2e2*x1mix3*x2mix3e2*x1mix4*x2mix4e2*x3mix4 + 
					x4*(-f3*x1*x1mix2e2*x2e2*x1mix4*x2mix4e2 + 
					x3*x3mix4*(f1*x2e2*x2mix3e2*x2mix4e2 + f2*x1*x1mix3*x1mix4*
					(x1*(3*x2e2 + x3*x4 - 2*x2*(x3 + x4)) + 
					x2*(-4*x2e2 - 2*x3*x4 + 3*x2*(x3 + x4))))))/
					(x1*x1mix2*x2*x1mix3*x2mix3*x3*x1mix4*x2mix4*x3mix4*x4);
		if(IsIncreasing && (Der < 0.)) Der = 0.;
		return Der;
	}
	else if(PoIndx==3)
	{
		double x1mix3e2 = x1mix3*x1mix3, x2mix3e2 = x2mix3*x2mix3, x3mix4e2 = x3mix4*x3mix4;
		double x1mix3e3 = x1mix3e2*x1mix3;
		double Der = (-f4*x1*x1mix2*x2*x1mix3e2*x2mix3e2*x3e2 + 
					f0*x1mix2*x1mix3e2*x2mix3e2*x1mix4*x2mix4*x3mix4e2 + 
					x4*(f3*x1*x1mix2*x2*x1mix4*x2mix4*
					(x3*(2*x1*x2 - 3*(x1 + x2)*x3 + 4*x3e2) - 
					(x1*x2 - 2*(x1 + x2)*x3 + 3*x3e2)*x4) + 
					x3e2*x3mix4e2*(f2*x1*x1mix3e2*x1mix4 + 
					f1*x2*x2mix3e2*(-x2 + x4))))/
					(x1*x1mix2*x2*x1mix3*x2mix3*x3*x1mix4*x2mix4*x3mix4*x4);
		if(IsIncreasing && (Der < 0.)) Der = 0.;
		return Der;
	}
	else if(PoIndx==4)
	{
		double x4e3 = x4e2*x4;
		double x2mix4e2 = x2mix4*x2mix4, x3mix4e2 = x3mix4*x3mix4, x1mix3e2 = x1mix3*x1mix3, x1mix4e2 = x1mix4*x1mix4;
		double x1mix3e3 = x1mix3e2*x1mix3;
		double Der = (-f0*x1mix2*x1mix3*x2mix3*x1mix4e2*x2mix4e2*
					x3mix4e2 + (f3*x1*x1mix2*x2*x1mix4e2*x2mix4e2 + 
					x3*(-f2*x1*x1mix3*x1mix4e2 + f1*x2*x2mix3*x2mix4e2)*
					x3mix4e2)*x4e2 + f4*x1*x1mix2*x2*x1mix3*x2mix3*x3*
					(x1*x2*x3 - 2*(x2*x3 + x1*(x2 + x3))*x4 + 3*(x1 + x2 + x3)*x4e2 - 4*x4e3))/
					(x1*x1mix2*x2*x1mix3*x2mix3*x3*x1mix4*x2mix4*x3mix4*x4);
		if(IsIncreasing && (Der < 0.)) Der = 0.;
		return Der;
	}
	else return 0.;
}

//-------------------------------------------------------------------------

double radTNonlinearIsotropMaterial::Derivative3(TVector2d* f, int PoIndx)
{
	double x0 = f[0].x, x1 = f[1].x, x2 = f[2].x;
	x1 -= x0; x2 -= x0;
	double f0 = f[0].y, f1 = f[1].y, f2 = f[2].y;

	char IsIncreasing = ((x1 > 0.) && (x2 > x1)) && ((f0 <= f1) && (f1 <= f2));

	double x1e2 = x1*x1, x2e2 = x2*x2;
	double x1mix2 = (x1-x2);
	if(PoIndx == 0)
	{
		double Der = (f0*x1e2 - f2*x1e2 - f0*x2e2 + f1*x2e2)/(-x1*x2*x1mix2);
		if(IsIncreasing && (Der < 0.)) Der = 0.;
		return Der;
	}
	else if(PoIndx == 1)
	{
		double Der = (-f2*x1e2 + f0*x1mix2*x1mix2 + f1*(2*x1 - x2)*x2)/(x1*x1mix2*x2);
		if(IsIncreasing && (Der < 0.)) Der = 0.;
		return Der;
	}
	else if(PoIndx == 2)
	{
		double Der = (f2*x1*(x1 - 2*x2) - f0*x1mix2*x1mix2 + f1*x2e2)/(x1*x1mix2*x1mix2*x2);
		if(IsIncreasing && (Der < 0.)) Der = 0.;
		return Der;
	}
	else return 0.;
}

//-------------------------------------------------------------------------

// Compute dB/dH derivatives for B-H curve (ELF-compatible)
// Now processes B-H curve directly (ArrayHB stores H,B pairs)
void radTNonlinearIsotropMaterial::Compute_dBdH(TVector2d* ArrayHB, double* dBdH, int LenArrayHB, double& MaxKsi)
{
	static const double mu_0 = 4.0e-7 * 3.14159265358979323846;  // 4*pi*1e-7

	if((ArrayHB == 0) || (dBdH == 0) || (LenArrayHB <= 0)) return;
	double *tdBdH = dBdH;
	TVector2d *tHB = ArrayHB;

	MaxKsi = 0;

	int LenArrayHB_mi_4 = LenArrayHB - 4;
	int LenArrayHB_mi_1 = LenArrayHB - 1;

	*tdBdH = Derivative3(tHB, 0);
	double ksi = *tdBdH / mu_0 - 1.0;  // chi = dB/dH / mu_0 - 1
	if(MaxKsi < ksi) MaxKsi = ksi;
	tdBdH++;

	*tdBdH = Derivative3(tHB, 1);
	ksi = *tdBdH / mu_0 - 1.0;
	if(MaxKsi < ksi) MaxKsi = ksi;
	tdBdH++; tHB++;

	*tdBdH = Derivative3(tHB, 1);
	ksi = *tdBdH / mu_0 - 1.0;
	if(MaxKsi < ksi) MaxKsi = ksi;
	tdBdH++;

	for(int i=3; i<LenArrayHB_mi_4; i++)
	{
		*tdBdH = Derivative5(tHB, 2);
		ksi = *tdBdH / mu_0 - 1.0;
		if(MaxKsi < ksi) MaxKsi = ksi;
		tdBdH++; tHB++;
	}

	if(LenArrayHB >= 7)
	{
		*tdBdH = Derivative5(tHB, 2);
		ksi = *tdBdH / mu_0 - 1.0;
		if(MaxKsi < ksi) MaxKsi = ksi;
		tdBdH++; tHB += 2;

		*tdBdH = Derivative3(tHB, 1);
		ksi = *tdBdH / mu_0 - 1.0;
		if(MaxKsi < ksi) MaxKsi = ksi;
		tdBdH++; tHB++;

		*tdBdH = Derivative3(tHB, 1);
		ksi = *tdBdH / mu_0 - 1.0;
		if(MaxKsi < ksi) MaxKsi = ksi;
		tdBdH++; tHB++;

		*tdBdH = ((tHB+1)->y - tHB->y)/((tHB+1)->x - tHB->x);
		ksi = *tdBdH / mu_0 - 1.0;
		if(MaxKsi < ksi) MaxKsi = ksi;
	}
	else if(LenArrayHB >= 6)
	{
		*tdBdH = Derivative5(tHB, 2);
		ksi = *tdBdH / mu_0 - 1.0;
		if(MaxKsi < ksi) MaxKsi = ksi;
		tdBdH++; tHB += 2;

		*tdBdH = Derivative3(tHB, 1);
		ksi = *tdBdH / mu_0 - 1.0;
		if(MaxKsi < ksi) MaxKsi = ksi;
		tdBdH++; tHB++;

		*tdBdH = ((tHB+1)->y - tHB->y)/((tHB+1)->x - tHB->x);
		ksi = *tdBdH / mu_0 - 1.0;
		if(MaxKsi < ksi) MaxKsi = ksi;
	}
	else if(LenArrayHB >= 5)
	{
		*tdBdH = Derivative3(tHB, 1);
		ksi = *tdBdH / mu_0 - 1.0;
		if(MaxKsi < ksi) MaxKsi = ksi;
		tdBdH++; tHB++;

		*tdBdH = ((tHB+1)->y - tHB->y)/((tHB+1)->x - tHB->x);
		ksi = *tdBdH / mu_0 - 1.0;
		if(MaxKsi < ksi) MaxKsi = ksi;
	}
	else if(LenArrayHB >= 4)
	{
		*tdBdH = ((tHB+1)->y - tHB->y)/((tHB+1)->x - tHB->x);
		ksi = *tdBdH / mu_0 - 1.0;
		if(MaxKsi < ksi) MaxKsi = ksi;
	}

	CheckAndCorrect_dBdH(ArrayHB, dBdH, LenArrayHB);
}

//-------------------------------------------------------------------------

// Check and correct dB/dH derivatives for B-H curve (ELF-compatible)
void radTNonlinearIsotropMaterial::CheckAndCorrect_dBdH(TVector2d* ArrayHB, double* dBdH, int LenArrayHB)
{
	if((ArrayHB == 0) || (dBdH == 0) || (LenArrayHB <= 0)) return;

	int LenArrayHB_mi_1 = LenArrayHB - 1;
	double *tdBdH = dBdH + LenArrayHB_mi_1;
	TVector2d *tHB = ArrayHB + LenArrayHB_mi_1;

	for(int i=0; i<LenArrayHB_mi_1; i++)
	{
		tdBdH--; tHB--;
		double x2 = (tHB+1)->x - tHB->x;
		double f1 = tHB->y, f2 = (tHB+1)->y;  // B values
		double fd1 = *tdBdH, fd2 = *(tdBdH+1);

		double x2e2 = x2*x2, f1mif2 = f1 - f2;
		double x2e3 = x2e2*x2;

		double a0 = f1;
		double a1 = fd1;
		double a2 = -((3*f1mif2 + (2*fd1 + fd2)*x2)/(x2e2));
		double a3 = -((-2*f1mif2 - (fd1 + fd2)*x2)/x2e3);

		double D = a2*a2 - 3.*a1*a3;
		double xc1 = 1.E+23, xc2 = 1.E+23;
		if(a3 != 0.)
		{
			if(D >= 0.)
			{
				double R = sqrt(D), Buf = 1./(3.*a3);
				xc1 = (-a2 - R)*Buf;
				xc2 = (-a2 + R)*Buf;
			}
		}

		char IsIncreasing = ((f1 < f2) && (fd1 > 0.) && (fd2 > 0.));
		char IsDecreasing = ((f1 > f2) && (fd1 < 0.) && (fd2 < 0.));
		char IsMonotone = (IsIncreasing || IsDecreasing);

		char xc1IsInside = 0, xc2IsInside = 0;
		if(a3 != 0.)
		{
			xc1IsInside = ((0. < xc1) && (xc1 < x2));
			xc2IsInside = ((0. < xc2) && (xc2 < x2));
		}

		char CorrectionNeeded = (IsMonotone && (D > 0.) && (xc1IsInside || xc2IsInside));
		if(CorrectionNeeded)
		{
			double LocD = -3*fd2*x2e3*(4*f1mif2 + fd2*x2);
			if(LocD >= 0.)
			{
				double LocR = sqrt(LocD), LocBuf1 = 1./(2*x2e2), LocBuf2 = -x2*(6*f1mif2 + fd2*x2);
				double fd1a = (LocBuf2 - LocR)*LocBuf1, fd1b = (LocBuf2 + LocR)*LocBuf1;
				double ra = fabs(fd1 - fd1a), rb = fabs(fd1 - fd1b);
				*tdBdH = (ra < rb)? fd1a : fd1b;
			}
		}
	}
}

//-------------------------------------------------------------------------

// ELF-compatible B(H) interpolation approach
// ArrayHB now stores (H, B) pairs directly, NOT (H, M)
// Returns M = B/mu_0 - H computed from interpolated B value
double radTNonlinearIsotropMaterial::AbsMvsAbsH_Interpol(double AbsH, TVector2d* ArrayHB, double* dBdH, int LenArrayHB)
{
	static const double mu_0 = 4.0e-7 * 3.14159265358979323846;  // 4*pi*1e-7

	// Use LINEAR interpolation for B(H) - same as ELF
	TVector2d *tHB = ArrayHB;  // Renamed to clarify: stores (H, B) not (H, M)
	int Indx = 0;
	for(int i=0; i<LenArrayHB; i++)
	{
		if(tHB->x > AbsH) break;
		tHB++; Indx++;
	}
	tHB--; Indx--;

	double B;  // Interpolated B value

	if(Indx < 0)
	{
		// Below first point: linear extrapolation from origin
		// For H=0, B=0, so B/H = initial slope
		if(AbsH <= 1e-10) return 0.0;  // Avoid division by zero
		B = ArrayHB->y * AbsH / ArrayHB->x;  // Linear from origin
	}
	else if(Indx >= (LenArrayHB - 1))
	{
		// Above last point: use last point's value (saturation)
		TVector2d *pHB = ArrayHB + (LenArrayHB - 1);
		B = pHB->y;  // Saturated B value
	}
	else
	{
		// Linear interpolation: B = B1 + (B2 - B1) * (H - H1) / (H2 - H1)
		double H1 = tHB->x, H2 = (tHB + 1)->x;
		double B1 = tHB->y, B2 = (tHB + 1)->y;
		double t = (AbsH - H1) / (H2 - H1);
		B = B1 + t * (B2 - B1);
	}

	// Convert B to M: M = B/mu_0 - H
	double M = B / mu_0 - AbsH;
	return M;
}

//-------------------------------------------------------------------------

// ELF-compatible inverse interpolation: find H given M
// ArrayHB stores (H, B) pairs, need to find H such that B/mu_0 - H = M
double radTNonlinearIsotropMaterial::AbsHvsAbsM_Interpol(double AbsM, TVector2d* ArrayHB, double* dBdH, int LenArrayHB)
{
	static const double mu_0 = 4.0e-7 * 3.14159265358979323846;  // 4*pi*1e-7

	// Convert B-H array to M values for search
	// M = B/mu_0 - H for each point
	TVector2d *tHB = ArrayHB;
	int Indx = 0;

	// Find interval where M falls
	for(int i=0; i<LenArrayHB; i++)
	{
		double M_i = tHB->y / mu_0 - tHB->x;  // M = B/mu_0 - H
		if(M_i > AbsM) break;
		tHB++; Indx++;
	}
	tHB--; Indx--;

	if(Indx < 0) return ArrayHB->x;
	if(Indx >= (LenArrayHB - 1)) return ArrayHB[LenArrayHB - 1].x;

	// Linear interpolation in H
	double H1 = tHB->x, H2 = (tHB + 1)->x;
	double B1 = tHB->y, B2 = (tHB + 1)->y;
	double M1 = B1 / mu_0 - H1;
	double M2 = B2 / mu_0 - H2;
	double t = (AbsM - M1) / (M2 - M1);
	return H1 + t * (H2 - H1);
}

//-------------------------------------------------------------------------

// ELF Method 2: H+B sum interpolation (CGS normalized)
// Given hb_sum = H/H_scale + B/B_scale, find H and B on BH curve
// Note: ArrayHB stores (H, M) pairs - we compute B = mu_0 * (H + M)
void radTNonlinearIsotropMaterial::InterpolateBH_HBSum(double hb_sum_target, TVector2d* ArrayHB, int LenArrayHB,
                                                        double H_scale, double B_scale, double& H_out, double& B_out)
{
	static const double mu_0 = 4.0e-7 * 3.14159265358979323846;

	// Handle edge cases
	if(LenArrayHB < 2 || ArrayHB == nullptr)
	{
		// Linear material fallback
		H_out = hb_sum_target * H_scale / 2.0;
		B_out = mu_0 * H_out;
		return;
	}

	// Compute B from H, M: B = mu_0 * (H + M)
	// Note: gArrayHB stores (H, B) pairs, not (H, M)!
	// So just use B directly without conversion
	auto getB = [](double H, double B) { return B; };  // B is stored directly

	// Compute normalized H+B for first point
	double B0 = getB(ArrayHB[0].x, ArrayHB[0].y);
	double sum1 = ArrayHB[0].x / H_scale + B0 / B_scale;

	// Below first point: extrapolate from origin
	if(hb_sum_target <= sum1)
	{
		if(sum1 > 1.0e-10)
		{
			double t = hb_sum_target / sum1;
			H_out = t * ArrayHB[0].x;
			B_out = t * B0;
		}
		else
		{
			H_out = 0.0;
			B_out = 0.0;
		}
		return;
	}

	// Compute normalized H+B for last point
	int n = LenArrayHB;
	double Bn = getB(ArrayHB[n-1].x, ArrayHB[n-1].y);
	double sum2 = ArrayHB[n-1].x / H_scale + Bn / B_scale;

	// Above last point: extrapolate with last slope
	if(hb_sum_target >= sum2)
	{
		if(n >= 2)
		{
			double h1 = ArrayHB[n-2].x, h2 = ArrayHB[n-1].x;
			double b1 = getB(h1, ArrayHB[n-2].y);
			double b2 = getB(h2, ArrayHB[n-1].y);
			double sum1_last = h1 / H_scale + b1 / B_scale;
			double hs = h2 - h1;
			double bs = b2 - b1;
			double hbs = hs / H_scale + bs / B_scale;
			if(std::fabs(hbs) > 1.0e-10)
			{
				double hd = hb_sum_target - sum1_last;
				double hi = hd / hbs;
				H_out = hs * hi + h1;
				B_out = bs * hi + b1;
			}
			else
			{
				H_out = h2;
				B_out = b2;
			}
		}
		else
		{
			H_out = ArrayHB[n-1].x;
			B_out = Bn;
		}
		return;
	}

	// Find interval for interpolation
	for(int i = 0; i < n-1; i++)
	{
		double Bi = getB(ArrayHB[i].x, ArrayHB[i].y);
		double Bip1 = getB(ArrayHB[i+1].x, ArrayHB[i+1].y);
		double hb_sum_i = ArrayHB[i].x / H_scale + Bi / B_scale;
		double hb_sum_ip1 = ArrayHB[i+1].x / H_scale + Bip1 / B_scale;

		if(hb_sum_target >= hb_sum_i && hb_sum_target < hb_sum_ip1)
		{
			double h1 = ArrayHB[i].x, h2 = ArrayHB[i+1].x;
			double b1 = Bi, b2 = Bip1;
			double hs = h2 - h1;
			double bs = b2 - b1;
			double hbs = hs / H_scale + bs / B_scale;
			if(std::fabs(hbs) > 1.0e-10)
			{
				double hd = hb_sum_target - hb_sum_i;
				double hi = hd / hbs;
				H_out = hs * hi + h1;
				B_out = bs * hi + b1;
			}
			else
			{
				H_out = h1;
				B_out = b1;
			}
			return;
		}
	}

	// Fallback
	H_out = hb_sum_target * H_scale / 2.0;
	B_out = mu_0 * H_out;
}

//-------------------------------------------------------------------------

// ELF dual-method chi update
// Method 1: Standard mu = B(H)/(mu_0*H)
// Method 2: H+B sum interpolation
// Selects method with smaller |mu_new - mu_old|
double radTNonlinearIsotropMaterial::ComputeChiDualMethod(double H_mag, double mu_old) const
{
	static const double mu_0 = 4.0e-7 * 3.14159265358979323846;

	// CGS conversion factors (ELF uses CGS internally)
	const double H_scale = 79.577;   // A/m per Oe
	const double B_scale = 1.0e-4;   // T per Gauss

	if(gLenArrayHB < 2 || gArrayHB == nullptr)
	{
		// Fallback for analytical or undefined material
		return mu_old - 1.0;
	}

	double B_interp, dBdH_unused;

	// Method 1: Standard H-based interpolation: mu = B(H) / (mu_0 * H)
	double un1;
	if(H_mag > 1.0e-10)
	{
		// Get B from BH curve
		double M = AbsMvsAbsH_Interpol(H_mag, gArrayHB, gdBdH, gLenArrayHB);
		B_interp = mu_0 * (M + H_mag);  // B = mu_0 * (M + H)
		un1 = B_interp / (mu_0 * H_mag);
	}
	else
	{
		// Use 2nd BH point for initial chi (ELF mucal0 style)
		// gArrayHB stores (H, B) pairs directly
		double H1 = gArrayHB[1].x;
		double B1 = gArrayHB[1].y;  // This is B, not M!
		un1 = B1 / (mu_0 * H1);  // mu_r = B / (mu_0 * H)
	}

	// ELF Method 1 correction: hn = H + (bo - bn) / uo; unn = bn / hn; un = MAX(un, unn)
	if(mu_old > 1.0 && H_mag > 1.0e-10)
	{
		double B_old = mu_0 * mu_old * H_mag;
		double hn = H_mag + (B_old - B_interp) / (mu_0 * mu_old);
		if(hn > 1.0e-10)
		{
			double unn = B_interp / (mu_0 * hn);
			if(unn > un1) un1 = unn;
		}
	}

	// Method 2: H+B sum interpolation
	double B_old = mu_0 * mu_old * H_mag;
	double hb_sum = H_mag / H_scale + B_old / B_scale;
	double H_new, B_new;
	InterpolateBH_HBSum(hb_sum, gArrayHB, gLenArrayHB, H_scale, B_scale, H_new, B_new);

	double un2;
	if(H_new > 1.0e-10)
	{
		un2 = B_new / (mu_0 * H_new);
	}
	else
	{
		un2 = un1;  // Fallback to Method 1
	}

	// Select method with smaller change from mu_old
	double us1 = std::fabs(un1 - mu_old);
	double us2 = std::fabs(un2 - mu_old);
	double un = (us1 <= us2) ? un1 : un2;

	// Ensure mu_r >= 1 (chi >= 0)
	if(un < 1.0) un = 1.0;

	double chi_new = un - 1.0;
	if(chi_new < 1.0e-6) chi_new = 1.0e-6;

	return chi_new;
}

//-------------------------------------------------------------------------

// Get H from M using Newton iteration for analytical formula (tanh model)
double radTNonlinearIsotropMaterial::GetHfromM_Analytical(double AbsM) const
{
	// For tanh model: M = sum(ms_i * tanh(ks_i * H / ms_i))
	// Use Newton iteration to find H given M
	const int MaxIter = 50;
	const double Tol = 1.0e-8;

	double H = AbsM / gMaxKsi;  // Initial guess using max susceptibility
	if(H < 1.0e-6) H = 1.0e-6;

	for(int iter = 0; iter < MaxIter; iter++)
	{
		double M_calc = 0.0;
		double dMdH = 0.0;
		for(int i = 0; i < lenMs_ks; i++)
		{
			if(Ms[i] == 0.) continue;
			double arg = ks[i] * H / Ms[i];
			double th = tanh(arg);
			M_calc += Ms[i] * th;
			double ch = cosh(arg);
			dMdH += ks[i] / (ch * ch);
		}

		double dM = M_calc - AbsM;
		if(fabs(dM) < Tol) break;
		if(fabs(dMdH) < 1.0e-15) break;

		H = H - dM / dMdH;
		if(H < 0) H = 1.0e-6;
	}

	return H;
}

//-------------------------------------------------------------------------

// ELF-compatible: returns M and dM/dH from B-H curve
// ArrayHB stores (H, B) pairs, computes M = B/mu_0 - H and dM/dH = dB/dH/mu_0 - 1
void radTNonlinearIsotropMaterial::AbsMvsAbsH_FuncAndDer_Interpol(double AbsH, TVector2d* ArrayHB, double* dBdH, int LenArrayHB, double& f, double& fDer)
{
	static const double mu_0 = 4.0e-7 * 3.14159265358979323846;  // 4*pi*1e-7

	TVector2d *tHB = ArrayHB;
	int Indx = 0;
	for(int i=0; i<LenArrayHB; i++)
	{
		if(tHB->x > AbsH) break;
		tHB++; Indx++;
	}
	tHB--; Indx--;

	double B, dB_dH;

	if(Indx < 0)
	{
		// Below first point: linear extrapolation from origin
		if(AbsH <= 1e-10) {
			f = 0.0;
			fDer = dBdH[0] / mu_0 - 1.0;
			return;
		}
		B = ArrayHB->y * AbsH / ArrayHB->x;
		dB_dH = ArrayHB->y / ArrayHB->x;  // Initial slope
	}
	else if(Indx >= (LenArrayHB - 1))
	{
		// Above last point: saturation
		B = ArrayHB[LenArrayHB - 1].y;
		dB_dH = 0.0;  // Saturated
	}
	else
	{
		// Linear interpolation: B = B1 + (B2 - B1) * (H - H1) / (H2 - H1)
		double H1 = tHB->x, H2 = (tHB + 1)->x;
		double B1 = tHB->y, B2 = (tHB + 1)->y;
		double t = (AbsH - H1) / (H2 - H1);
		B = B1 + t * (B2 - B1);
		dB_dH = (B2 - B1) / (H2 - H1);
	}

	// Convert to M and dM/dH
	// M = B/mu_0 - H
	// dM/dH = dB/dH / mu_0 - 1
	f = B / mu_0 - AbsH;
	fDer = dB_dH / mu_0 - 1.0;
}

//-------------------------------------------------------------------------
