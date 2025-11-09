/*-------------------------------------------------------------------------
*
* File name:      radplnr2.cpp
*
* Project:        RADIA
*
* Description:    Auxiliary 2D objects
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#include "rad_planar_2d.h"
#include "radappl.h"
#include "radpoly_analytical.h"

//-------------------------------------------------------------------------

extern radTConvergRepair& radCR;
extern radTYield radYield;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTPolygon::B_comp(radTField* FieldPtr)
{
// Orientation: The Polygon normal parallel to vertical ort !!!
// NOTE: This function uses the analytical formula from radpoly_analytical.cpp
//       for improved performance and simplicity.

	const double PI = 3.14159265358979;
	const double ConstForH = 1./4./PI;

	// Check for multitasking
	if(radYield.Check()==0) return;

	// Get observation point
	TVector3d& ObsPo = FieldPtr->P;

	// Check which field components are needed
	short A_CompNeeded = FieldPtr->FieldKey.A_;
	short B_orH_CompNeeded = FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_ || FieldPtr->FieldKey.PreRelax_;

	if(!A_CompNeeded && !B_orH_CompNeeded) return;

	// Handle singularity: if observation point is on polygon plane
	double z = CoordZ - ObsPo.z;
	TVector3d ObsPoOriginal = ObsPo;  // Save original
	if(z == 0.)
	{
		// Shift observation point slightly
		double AbsRandZ = radCR.AbsRandMagnitude(CoordZ);
		if(AbsRandZ == 0.) AbsRandZ = 1.e-15;
		z = AbsRandZ;
		ObsPo.z -= AbsRandZ;
	}

	// ========================================================================
	// Use analytical formula for field calculation
	// ========================================================================

	// Setup local coordinate system (polygon is in XY plane at z=CoordZ)
	TVector3d AA(1, 0, 0);  // Local X-axis
	TVector3d BB(0, 1, 0);  // Local Y-axis
	TVector3d CC(0, 0, 1);  // Normal (Z-axis)
	TVector3d YY(0, 0, CoordZ);  // Reference point on polygon plane

	// Prepare observation point in single-element vector
	std::vector<TVector3d> obs_points(1, ObsPo);
	std::vector<TVector3d> field_result(1, TVector3d(0, 0, 0));

	// Magnetic charge density (from Magn.z)
	double W = ConstForH * Magn.z;

	// Call analytical formula
	RadAnalyticalFieldFromPolygonCharge(
		AA, BB, CC, YY,
		EdgePointsVector,  // Polygon vertices in 2D
		obs_points,
		field_result,
		W,
		1,  // Element index (for error reporting)
		AmOfEdgePoints
	);

	// Extract result
	TVector3d& H_field = field_result[0];

	// ========================================================================
	// Apply results to FieldPtr
	// ========================================================================

	if(B_orH_CompNeeded)
	{
		if(FieldPtr->FieldKey.PreRelax_)
		{
			// Special handling for relaxation
			TVector3d St0(0., 0., -H_field.x);
			TVector3d St1(0., 0., -H_field.y);
			TVector3d St2(0., 0., -H_field.z);
			FieldPtr->B += St0;
			FieldPtr->H += St1;
			FieldPtr->A += St2;
		}
		else
		{
			// Normal case: add to H field
			FieldPtr->H += H_field;
		}
	}

	if(A_CompNeeded)
	{
		// Vector potential calculation
		// A = -z * H_z in the direction perpendicular to Magn
		double AS = -z * H_field.z;
		TVector3d BufA(-Magn.y, Magn.x, 0.);
		FieldPtr->A += AS * BufA;
	}

	// Restore original observation point
	ObsPo = ObsPoOriginal;
}

//-------------------------------------------------------------------------

void radTPolygon::B_intComp(radTField* FieldPtr)
{
// Orientation: The Polygon normal parallel to vertical ort !!!
// Field Integrals are correct only if line does not cross the prism body.
// This was noticed to have the "divide by zero" problem if (Vx*Vx+Vy*Vy==0)||(Vx*Vx+Vz*Vz==0)||(Vy*Vy+Vz*Vz==0)

	if(FieldPtr->FieldKey.FinInt_) { B_intCompFinNum(FieldPtr); return;}

	const double PI = 3.14159265358979;
	const double ConstForH = 1./4./PI;
	const double Max_k = 1.E+08; // To skip segments in general field int. computation loop.
	const double ZeroToler = 1.E-06; // This is to switch to Special Cases
	const double SmallestRelTolerV = 1.E-12; // Relative tolerance to repair trapping V.i to zero in general case

	TVector3d V = FieldPtr->NextP - FieldPtr->P;
	double InvAbsV = 1./sqrt(V.x*V.x + V.y*V.y + V.z*V.z);
	V = InvAbsV*V;

	short InitVxIsZero = (Abs(V.x) <= ZeroToler);
	short InitVyIsZero = (Abs(V.y) <= ZeroToler);
	short InitVzIsZero = (Abs(V.z) <= ZeroToler);

// This is the attempt to avoid "divide by zero" problem
	TSpecCaseID SpecCaseID = TSpecCaseID::NoSpecCase;
	if(InitVxIsZero && InitVyIsZero) SpecCaseID = TSpecCaseID::ZeroVxVy;
	else if(InitVxIsZero && InitVzIsZero) SpecCaseID = TSpecCaseID::ZeroVxVz;
	else if(InitVyIsZero && InitVzIsZero) SpecCaseID = TSpecCaseID::ZeroVyVz;
	if(SpecCaseID==TSpecCaseID::ZeroVxVy || SpecCaseID==TSpecCaseID::ZeroVxVz || SpecCaseID==TSpecCaseID::ZeroVyVz) { B_intCompSpecCases(FieldPtr, SpecCaseID); return;}

	double AbsRandX = radCR.AbsRandMagnitude(CentrPoint.x);
	double AbsRandY = radCR.AbsRandMagnitude(CentrPoint.y);
	double AbsRandZ = radCR.AbsRandMagnitude(CoordZ);

	double Vx = V.x; if(Vx==0.) Vx = SmallestRelTolerV;
	double Vy = V.y; if(Vy==0.) Vy = SmallestRelTolerV;
	double Vz = V.z; if(Vz==0.) Vz = SmallestRelTolerV;

	double Vxe2 = Vx*Vx, Vye2 = Vy*Vy, Vze2 = Vz*Vz;
	double Vxe2pVye2 = Vxe2+Vye2, Vxe2pVze2 = Vxe2+Vze2, Vye2pVze2 = Vye2+Vze2;
	double Vye2pVze2Vy = Vye2pVze2*Vy, Vx1pVye2 = Vx*(1.+Vye2);
	double VxVy = Vx*Vy, VyVz = Vy*Vz;

	double p2dVxe2pVze2 = 2./Vxe2pVze2;
	TVector3d& StPo = FieldPtr->P;

	double z = CoordZ - StPo.z;
// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
	if(z==0.) z = AbsRandZ;
	double Vxz = Vx*z, Vyz = Vy*z, VyVzz = VyVz*z, Vzz = Vz*z;

	radTVect2dVect::iterator BaseIter = EdgePointsVector.begin();
	int AmOfEdgePoints_mi_1 = AmOfEdgePoints - 1;

	TVector2d& FirstPo2d = *BaseIter;
	TVector2d First2d(FirstPo2d.x - StPo.x, FirstPo2d.y - StPo.y);
	TVector2d Vect2dToAdd(-StPo.x, -StPo.y);

	double x1 = First2d.x;
	double y1 = First2d.y;
// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
	if(x1==0.) x1 = AbsRandY;
	if(y1==0.) y1 = AbsRandZ;
	double x2, y2;

	double Sx=0., Sy=0., Sz=0.;

	for(int i=0; i<AmOfEdgePoints; i++)
	{
		++BaseIter;
		if(i!=AmOfEdgePoints_mi_1)
		{
			TVector2d& BufPo2d = *BaseIter;
			x2 = BufPo2d.x + Vect2dToAdd.x;
			y2 = BufPo2d.y + Vect2dToAdd.y;
		}
		else
		{
			x2 = First2d.x;
			y2 = First2d.y;
		}
		// Artificial shift of an observation point a bit right of the block's border
		// if the point is exactly on the boarder (to avoid "divide by zero" error):
		if(x2==0.) x2 = AbsRandX;
		if(y2==0.) y2 = AbsRandY;

		double x2mx1 = x2-x1;
		double y2my1 = y2-y1;
		double abs_x2mx1 = Abs(x2mx1), abs_y2my1 = Abs(y2my1);

		if(abs_x2mx1*Max_k > abs_y2my1)
		{
			double k = y2my1/x2mx1, b = y1 - k*x1;
			double ke2 = k*k;
			double ke2p1 = ke2 + 1.;
			double kVx = k*Vx, kVy = k*Vy, kVz = k*Vz;
			double kVxmVy = kVx - Vy, kVypVx = kVy + Vx;
			double kVypVxVz = kVypVx*Vz;
			double kVxmVye2p1pke2Vze2 = kVxmVy*kVxmVy + ke2p1*Vze2;
			double bVx = b*Vx, bVz = b*Vz;
			double kVxe2pVze2mVxVy = k*Vxe2pVze2 - VxVy;
			double bkVxe2pVze2mVxVy = b*kVxe2pVze2mVxVy;
			double kVxmVyx1 = kVxmVy*x1, kVxmVyx2 = kVxmVy*x2;
			double kVxmVye2p1pke2Vze2x1 = kVxmVye2p1pke2Vze2*x1, kVxmVye2p1pke2Vze2x2 = kVxmVye2p1pke2Vze2*x2;
			double kVypVxVzz = kVypVxVz*z, kVxmVyz = kVxmVy*z;
			double bVxpkVxmVyx1 = bVx + kVxmVyx1, bVxpkVxmVyx2 = bVx + kVxmVyx2;
			double Vzx1 = Vz*x1, Vzx2 = Vz*x2;
			double Vzx1mVxz = Vzx1-Vxz, Vzx2mVxz = Vzx2-Vxz;
			double kVzx1 = kVz*x1, kVzx2 = kVz*x2;
			double bVzpkVzx1mVyz = bVz + kVzx1 - Vyz, bVzpkVzx2mVyz = bVz + kVzx2 - Vyz;
			double bVzpkVxmVyz = bVz + kVxmVyz; 
			double bVxpkVxmVyx1e2 = bVxpkVxmVyx1*bVxpkVxmVyx1, bVxpkVxmVyx2e2 = bVxpkVxmVyx2*bVxpkVxmVyx2;
			double bVxe2pVze2 = b*Vxe2pVze2;
			double kVxe2pVze2mVxVyx1 = kVxe2pVze2mVxVy*x1, kVxe2pVze2mVxVyx2 = kVxe2pVze2mVxVy*x2;

			double PiMult1=0.;
			double SumAtans1 = atan(TransAtans(-(bkVxe2pVze2mVxVy + kVxmVye2p1pke2Vze2x1 - kVypVxVzz)/bVzpkVxmVyz, (bkVxe2pVze2mVxVy + kVxmVye2p1pke2Vze2x2 - kVypVxVzz)/bVzpkVxmVyz, PiMult1));
			SumAtans1 += PiMult1*PI;

			double AtanX1 = atan((bVxe2pVze2 + kVxe2pVze2mVxVyx1 - VyVzz)/(Vxz - Vzx1));
			double AtanX2 = atan((bVxe2pVze2 + kVxe2pVze2mVxVyx2 - VyVzz)/(Vxz - Vzx2));

			double LogX1 = log(bVxpkVxmVyx1e2 + bVzpkVzx1mVyz*bVzpkVzx1mVyz + Vzx1mVxz*Vzx1mVxz);
			double LogX2 = log(bVxpkVxmVyx2e2 + bVzpkVzx2mVyz*bVzpkVzx2mVyz + Vzx2mVxz*Vzx2mVxz);

			double kVypVxVzVz = kVypVxVz*Vz;
			double kVxVy = k*VxVy;
			double BufLogMult1 = (kVxVy-Vye2pVze2)*bVxe2pVze2 + Vzz*(Vye2pVze2Vy-k*Vx1pVye2);

			double BufLogMult2 = (kVypVxVzz - b*kVxe2pVze2mVxVy)/kVxmVye2p1pke2Vze2;
			double BufLogMult3 = kVypVxVz*bVxe2pVze2 + (Vx*kVxmVy - Vy*kVypVxVzVz)*z;

			double p2dVxe2pVze2dkVxmVye2p1pke2Vze2 = p2dVxe2pVze2/kVxmVye2p1pke2Vze2;

			Sx += p2dVxe2pVze2dkVxmVye2p1pke2Vze2*((kVz*bVxe2pVze2 + (VxVy*kVxmVy - kVypVxVzVz)*z)*SumAtans1 + Vz*(kVxmVye2p1pke2Vze2x2*AtanX2 - kVxmVye2p1pke2Vze2x1*AtanX1) 
				+ 0.5*((VxVy*kVxmVye2p1pke2Vze2x2 + BufLogMult1)*LogX2 - (VxVy*kVxmVye2p1pke2Vze2x1 + BufLogMult1)*LogX1));

			Sy += -2.*(bVzpkVxmVyz/kVxmVye2p1pke2Vze2)*SumAtans1 + (BufLogMult2-x2)*LogX2 - (BufLogMult2-x1)*LogX1;

			Sz += p2dVxe2pVze2dkVxmVye2p1pke2Vze2*(-(b*kVxmVy*Vxe2pVze2 + (Vye2-Vxe2-2.*kVxVy)*Vzz)*SumAtans1 - Vx*(kVxmVye2p1pke2Vze2x2*AtanX2 - kVxmVye2p1pke2Vze2x1*AtanX1)
				+ 0.5*((VyVz*kVxmVye2p1pke2Vze2x2 + BufLogMult3)*LogX2 - (VyVz*kVxmVye2p1pke2Vze2x1 + BufLogMult3)*LogX1));
		}
		x1 = x2; y1 = y2;
	}
	TVector3d BufIh(Sx, Sy, Sz);
	double MultIh = -ConstForH*Magn.z; // Check if "-" is necessary
	BufIh = MultIh*BufIh;

	if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufIh;
	if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufIh;
}

//-------------------------------------------------------------------------

void radTPolygon::B_intCompSpecCases(radTField* FieldPtr, const TSpecCaseID& SpecCaseID)
{
	const double PI = 3.14159265358979;
	const double ConstForH = 1./4./PI;

	const double Max_k = 1.E+08;

	TVector3d& StPo = FieldPtr->P;

	double AbsRandX = radCR.AbsRandMagnitude(CentrPoint.x);
	double AbsRandY = radCR.AbsRandMagnitude(CentrPoint.y);
	double AbsRandZ = radCR.AbsRandMagnitude(CoordZ);

	double z = CoordZ - StPo.z;
// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
	if(z==0.) z = AbsRandZ;
	double ze2 = z*z;

	radTVect2dVect::iterator BaseIter = EdgePointsVector.begin();
	int AmOfEdgePoints_mi_1 = AmOfEdgePoints - 1;

	TVector2d& FirstPo2d = *BaseIter;
	TVector2d First2d(FirstPo2d.x - StPo.x, FirstPo2d.y - StPo.y);
	TVector2d Vect2dToAdd(-StPo.x, -StPo.y);

	double x1 = First2d.x;
	double y1 = First2d.y;
// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
	if(x1==0.) x1 = AbsRandX;
	if(y1==0.) y1 = AbsRandY;
	double x2, y2;
	double x1e2 = x1*x1, x2e2, y1e2 = y1*y1, y2e2;

	double Sx=0., Sy=0., Sz=0.;
	double PiMult1=0.;

	for(int i=0; i<AmOfEdgePoints; i++)
	{
		++BaseIter;
		if(i!=AmOfEdgePoints_mi_1)
		{
			TVector2d& BufPo2d = *BaseIter;
			x2 = BufPo2d.x + Vect2dToAdd.x;
			y2 = BufPo2d.y + Vect2dToAdd.y;
		}
		else
		{
			x2 = First2d.x;
			y2 = First2d.y;
		}
		// Artificial shift of an observation point a bit right of the block's border
		// if the point is exactly on the boarder (to avoid "divide by zero" error):
		if(x2==0.) x2 = AbsRandX;
		if(y2==0.) y2 = AbsRandY;
		x2e2 = x2*x2; y2e2 = y2*y2;

		double x2mx1 = x2-x1;
		double y2my1 = y2-y1;
		double abs_x2mx1 = Abs(x2mx1), abs_y2my1 = Abs(y2my1);

		if(abs_x2mx1*Max_k > abs_y2my1)
		{
			if(SpecCaseID==TSpecCaseID::ZeroVxVy)
			{
				double k = y2my1/x2mx1, b = y1 - k*x1;
				double ke2 = k*k, bk = b*k;
				double ke2p1 = ke2 + 1.;

				double AtanX1 = atan(k+b/x1);
				double AtanX2 = atan(k+b/x2);
				double SumAtans1 = atan(TransAtans((bk+ke2p1*x2)/b, -(bk+ke2p1*x1)/b, PiMult1));

				SumAtans1 += PiMult1*PI;
				double Log1 = log(x1e2 + y1e2);
				double Log2 = log(x2e2 + y2e2);
				double SumLogs1 = Log2 - Log1;
				double bdke2p1 = b/ke2p1;
				double bkdke2p1 = bdke2p1*k;

				Sx += -2.*((x2*AtanX2 - x1*AtanX1) - bkdke2p1*SumAtans1) - bdke2p1*SumLogs1;
				Sy += -2.*bdke2p1*SumAtans1 - ((bkdke2p1+x2)*Log2 - (bkdke2p1+x1)*Log1);
			}
			else if(SpecCaseID==TSpecCaseID::ZeroVxVz)
			{
				double k = y2my1/x2mx1, b = y1 - k*x1;
				double kz = k*z;
				double SumAtans1 = atan(TransAtans(x2/z, -x1/z, PiMult1));
				SumAtans1 += PiMult1*PI;
				double SumLogs1 = log((x2e2 + ze2)/(x1e2 + ze2));

				Sx += 2.*(-y2my1 + kz*SumAtans1) - b*SumLogs1;
				Sz += -2.*b*SumAtans1 - kz*SumLogs1;
			}
		}
		if(abs_y2my1*Max_k > abs_x2mx1)
		{
			if(SpecCaseID==TSpecCaseID::ZeroVyVz)
			{
				double k1 = x2mx1/y2my1, b1 = x1 - k1*y1;
				double k1z = k1*z;
				double SumAtans1 = atan(TransAtans(y2/z, -y1/z, PiMult1));
				SumAtans1 += PiMult1*PI;
				double SumLogs1 = log((y2e2 + ze2)/(y1e2 + ze2));

				Sy -= 2.*(-x2mx1 + k1z*SumAtans1) - b1*SumLogs1;
				Sz -= -2.*b1*SumAtans1 - k1z*SumLogs1;
			}
		}
		x1 = x2; y1 = y2;
		x1e2 = x2e2; y1e2 = y2e2;
	}
	TVector3d BufIh(Sx, Sy, Sz);
	double MultIh = -ConstForH*Magn.z; // Check if "-" is necessary
	BufIh = MultIh*BufIh;

	if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufIh;
	if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufIh;
}

//-------------------------------------------------------------------------

//void radTPolygon::B_comp_frJ(radTField* pField)
//{
//	radTVect2dVect::iterator baseIter = EdgePointsVector.begin();
//
//	for(int i=0; i<AmOfEdgePoints; i++)
//	{
//		++baseIter;
//	}
//}

//-------------------------------------------------------------------------
