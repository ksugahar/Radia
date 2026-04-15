/*-------------------------------------------------------------------------
*
* File name:      radflm.cpp
*
* Project:        RADIA
*
* Description:    Magnetic field source: filament conductor
*
* Author(s):      Oleg Chubar, Pascal Elleaume
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#include "rad_application.h"
#include "rad_filament.h"

#include <math.h>
#include <sstream>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTFlmLinCur::radTFlmLinCur(const TVector3d& InStartPoint, const TVector3d& InEndPoint, double InI)
{
	I = InI; StartPoint = InStartPoint; EndPoint = InEndPoint;

	TVector3d LinVect = EndPoint - StartPoint;
	double SqLength = LinVect.x*LinVect.x + LinVect.y*LinVect.y + LinVect.z*LinVect.z;
	Length = sqrt(SqLength);

	if(LinVect.y==0. && LinVect.z==0.)
	{
		TVector3d St0(1.,0.,0.);
		TVector3d St1(0.,1.,0.);
		TVector3d St2(0.,0.,1.);
		TMatrix3d M(St0, St1, St2);
		TVector3d ZeroVect(0.,0.,0.);
		NativeRotation = radTrans(M, M, ZeroVect, 1., 1.); // Identity
		return;
	}

	TVector3d LinVectProto(Length, 0., 0.);
	TVector3d RotAx(LinVectProto.y*LinVect.z-LinVectProto.z*LinVect.y,
					LinVectProto.z*LinVect.x-LinVectProto.x*LinVect.z,
					LinVectProto.x*LinVect.y-LinVectProto.y*LinVect.x);
	double cosAngle = (LinVectProto * LinVect)/SqLength;
	double Angle = acos(cosAngle);
	SetNativeRotation(RotAx, Angle);

	TVector3d TestVect = NativeRotation.TrBiPoint(LinVectProto);
	const double SmallPositive = 1.E-10;
	const double TwoPi = 2.*3.141592653589793238;
	if((Abs(TestVect.x-LinVect.x)/Length > SmallPositive) ||
	   (Abs(TestVect.y-LinVect.y)/Length > SmallPositive) ||
	   (Abs(TestVect.z-LinVect.z)/Length > SmallPositive)) SetNativeRotation(RotAx, TwoPi-Angle);
}

//-------------------------------------------------------------------------

void radTFlmLinCur::SetNativeRotation(const TVector3d& InAxVect, double Angle)
{
	double NormFact = 1./sqrt(InAxVect.x*InAxVect.x+InAxVect.y*InAxVect.y+InAxVect.z*InAxVect.z);
	TVector3d AxVect = NormFact*InAxVect;
	double VxVx, VyVy, VzVz;
	VxVx=AxVect.x*AxVect.x; VyVy=AxVect.y*AxVect.y; VzVz=AxVect.z*AxVect.z;

	double cosAng, sinAng, One_m_cos;
	cosAng = cos(Angle); sinAng = sin(Angle); One_m_cos = 1. - cosAng;
	double One_m_cosVxVy, One_m_cosVxVz, One_m_cosVyVz, sinVx, sinVy, sinVz;
	One_m_cosVxVy = One_m_cos*AxVect.x*AxVect.y;
	One_m_cosVxVz = One_m_cos*AxVect.x*AxVect.z;
	One_m_cosVyVz = One_m_cos*AxVect.y*AxVect.z;
	sinVx = sinAng*AxVect.x; sinVy = sinAng*AxVect.y; sinVz = sinAng*AxVect.z;

	TVector3d St0(VxVx+cosAng*(VyVy+VzVz), One_m_cosVxVy-sinVz, One_m_cosVxVz+sinVy);
	TVector3d St1(One_m_cosVxVy+sinVz, VyVy+cosAng*(VxVx+VzVz), One_m_cosVyVz-sinVx);
	TVector3d St2(One_m_cosVxVz-sinVy, One_m_cosVyVz+sinVx, VzVz+cosAng*(VxVx+VyVy));
	TMatrix3d M(St0, St1, St2);
	TVector3d St00(1.-St0.x, -St0.y, -St0.z);
	TVector3d St01(-St1.x, 1.-St1.y, -St1.z);
	TVector3d St02(-St2.x, -St2.y, 1.-St2.z);
	TMatrix3d M0(St00, St01, St02);

	NativeRotation = radTrans(M, M0*StartPoint, 1., 1.);
}

//-------------------------------------------------------------------------

void radTFlmLinCur::B_comp(radTField* FieldPtr)
{
	TVector3d BufP = NativeRotation.TrPoint_inv(FieldPtr->P);
	TVector3d V0 = StartPoint - BufP;
	double x1 = V0.x + Length;

	const double SmallPositive = 1.E-23;
	double y0y0_p_z0z0 = V0.y*V0.y + V0.z*V0.z;
	if(y0y0_p_z0z0 == 0.) y0y0_p_z0z0 = SmallPositive;

	double SqRt0 = sqrt(V0.x*V0.x + y0y0_p_z0z0);
	double SqRt1 = sqrt(x1*x1 + y0y0_p_z0z0);

	double ConI = 1.E-04 * I;
	double ComMult;

	if(FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_)
	{
		ComMult = ConI*(V0.x/SqRt0 - x1/SqRt1)/y0y0_p_z0z0;
		TVector3d BufB(0., -ComMult*V0.z, ComMult*V0.y);
		BufB = NativeRotation.TrVectField(BufB);
		if(FieldPtr->FieldKey.B_) FieldPtr->B += BufB;
		if(FieldPtr->FieldKey.H_) FieldPtr->H += BufB;
	}
	if(FieldPtr->FieldKey.A_)
	{
		double V0x_p_SqRt0 = V0.x+SqRt0;
		if(V0x_p_SqRt0==0.) V0x_p_SqRt0 = SmallPositive;
		ComMult = ConI*log((x1+SqRt1)/V0x_p_SqRt0);
		TVector3d BufA(ComMult, 0., 0.);
		FieldPtr->A += NativeRotation.TrVectPoten(BufA);
	}
	if(FieldPtr->FieldKey.Phi_)
	{
		// Magnetic scalar potential from a line current segment
		// For a line current along x from x=0 to x=L (in local coordinates),
		// the solid angle subtended at point (x0, y0, z0) is computed as:
		//
		// The scalar potential of a finite line segment is:
		//   Phi = (I / 4*pi) * Omega
		// where Omega is the "solid angle" of the ribbon connecting the
		// line segment to the observation point.
		//
		// For numerical stability, we use:
		//   Phi = (I / 4*pi) * [atan2(y*z, r0*|x0|) - atan2(y*z, r1*|x1|)]
		//       when the point is not on the current axis
		//
		// Reference: J.D. Jackson, "Classical Electrodynamics", 3rd ed.
		//
		// Note: In the local frame, the wire is along x-axis from x=0 to x=Length.
		// V0 = StartPoint - P gives vector from P to StartPoint (x=0).
		// V0.x is negative of x-coordinate of P relative to start.
		// So x0 = -V0.x (distance from P to start along x)
		// x1 = -(V0.x + Length) = -x1_local

		const double Pi = 3.141592653589793238;
		double rho2 = y0y0_p_z0z0;  // y^2 + z^2

		if(rho2 > SmallPositive) {
			// General case: point is not on the wire axis
			// Use the solid angle formula for a line segment
			//
			// The solid angle subtended by a line segment from point A to B
			// at observation point P is related to the angles at P.
			//
			// For a wire along x-axis, the scalar potential is:
			//   Phi = (I / 4*pi) * [atan(x0*y / (rho*r0)) - atan(x1*y / (rho*r1))]
			// where rho = sqrt(y^2 + z^2), r0 = |P - Start|, r1 = |P - End|
			//
			// Simplified form using signed angles:
			double x0 = -V0.x;  // x-coordinate relative to start point
			double x1_loc = x0 - Length;  // x-coordinate relative to end point
			double rho = sqrt(rho2);

			// Angles: atan2(y, x) gives angle from x-axis
			// For solid angle, we need the angle in the perpendicular plane
			double angle0 = atan2(V0.y * x0, rho * SqRt0);
			double angle1 = atan2(V0.y * x1_loc, rho * SqRt1);

			double Omega = angle0 - angle1;
			double Phi = I * Omega / (4.0 * Pi);
			FieldPtr->Phi += Phi;
		}
		// If on the wire axis (rho = 0), the potential is undefined/multivalued
		// We leave Phi unchanged (effectively 0 contribution)
	}
}

//-------------------------------------------------------------------------

void radTFlmLinCur::B_intComp(radTField* FieldPtr)
{
	if(FieldPtr->FieldKey.FinInt_) { B_intCompFinNum(FieldPtr); return;}

// An analytical algorithm for infinite Field Integrals:
	TVector3d BufP = NativeRotation.TrPoint_inv(FieldPtr->P);
	TVector3d BufNextP = NativeRotation.TrPoint_inv(FieldPtr->NextP);

	TVector3d VV = StartPoint - BufP;
	double x2 = VV.x + Length;

	TVector3d v = BufNextP - BufP;
	double Mod_v = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
	v.x /= Mod_v; v.y /= Mod_v; v.z /= Mod_v;

	double ConI = 2.E-04 * I;
	TVector3d BufIntB(0.,0.,0.);
	double ComMult;

	const double SpecCaseZeroToler = 1.E-12;
	double vyvy_p_vzvz = v.y*v.y + v.z*v.z;
	if(vyvy_p_vzvz < SpecCaseZeroToler)
	{
		ComMult = ConI*(x2-VV.x)/(VV.z*VV.z + VV.y*VV.y);
		BufIntB.y = ComMult*VV.z;
		BufIntB.z = -ComMult*VV.y;
		goto FinalDefinitionOfFieldIntegrals;
	}
	{
		const double Pi = 3.141592653589793238;
		const double SmallestZeroToler = 1.E-14;
		double vzY_mi_vyZ = v.z*VV.y - v.y*VV.z; if(vzY_mi_vyZ==0.) vzY_mi_vyZ = SmallestZeroToler;
		double vxZ = v.x*VV.z;
		double vzX1_mi_vxZ = v.z*VV.x - vxZ;
		double vzX2_mi_vxZ = v.z*x2 - vxZ;
		double vxY = v.x*VV.y;
		double vxY_mi_vyX1 = vxY - v.y*VV.x;
		double vxY_mi_vyX2 = vxY - v.y*x2;
		double vxvyY_p_vxvzZ = v.x*(v.y*VV.y+v.z*VV.z);
		double vzYmivyZ_mu_vzYmivyZ = vzY_mi_vyZ*vzY_mi_vyZ;

		double PiMult = 0.;
		double F = (atan(TransAtans((vyvy_p_vzvz*VV.x-vxvyY_p_vxvzZ)/vzY_mi_vyZ, -(vyvy_p_vzvz*x2-vxvyY_p_vxvzZ)/vzY_mi_vyZ, PiMult)) + Pi*PiMult)/vyvy_p_vzvz;
		double G = 0.5*v.x*log((vzYmivyZ_mu_vzYmivyZ + vzX1_mi_vxZ*vzX1_mi_vxZ + vxY_mi_vyX1*vxY_mi_vyX1)/(vzYmivyZ_mu_vzYmivyZ + vzX2_mi_vxZ*vzX2_mi_vxZ + vxY_mi_vyX2*vxY_mi_vyX2))/vyvy_p_vzvz;

		BufIntB.y = ConI*(v.y*F+v.z*G);
		BufIntB.z = ConI*(v.z*F-v.y*G);
	}
FinalDefinitionOfFieldIntegrals:
	BufIntB = NativeRotation.TrVectField(BufIntB);
	if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufIntB;
	if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufIntB;
}

//-------------------------------------------------------------------------


//-------------------------------------------------------------------------

// radTFlmLinCur::Dump / DumpBin REMOVED (Phase B2b/B2c, 2026-04-15)

//-------------------------------------------------------------------------

// radTFlmLinCur::SubdivideItself REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------
