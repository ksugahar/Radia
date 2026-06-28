/*-------------------------------------------------------------------------
*
* File name:      radrec.cpp
*
* Project:        RADIA
*
* Description:    Magnetic field source:
*                 rectangular parallelepiped with constant magnetization 
*                 or currect density
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

//-------------------------------------------------------------------------
// Implementation of class radTRecCur - a class of objects of rectangular
// parallelipipedic shape capable to generate magnetic field.
// RecMag is derived from radTg3d.
//-------------------------------------------------------------------------

#include "rad_rectangular_block.h"
#include "rad_group.h"
#include "rad_application.h"
#include "rad_geometry_3d_aux.h"

#ifndef _INC_MATH
#include <math.h>
#endif

//#ifdef __GNUC__
//#include <strstream.h>
//#else
#include <sstream>
//#endif

//-------------------------------------------------------------------------

extern radTYield radYield;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTRecCur::B_comp(radTField* FieldPtr)
{
	const double ZeroToler = 1.E-20;

	short J_is_Zero = 1;
	if(J.x!=0. || J.y!=0. || J.z!=0.) J_is_Zero = 0;
	TVector3d P_min_CenPo = FieldPtr->P - CentrPoint;

// [far-field magnetization multipole shortcut + its precompute REMOVED 2026-06-28: radTRecCur is current-only]

	if(radYield.Check()==0) return; // To allow multitasking on Mac: consider better places for this

	TVector3d T0(0,0,0), T1(0,0,0), T(0,0,0), S(1.,1.,1.);
	const double dConst1 = 0.5;
	TVector3d HalfDim = dConst1*Dimensions;

	struct { double x[2],y[2],z[2];} BfSt;
	for(int ii=0; ii<=1; ii++)
	{
		int Eps=ii*2-1;
		BfSt.x[ii] = -P_min_CenPo.x+Eps*HalfDim.x;
		BfSt.y[ii] = -P_min_CenPo.y+Eps*HalfDim.y;
		BfSt.z[ii] = -P_min_CenPo.z+Eps*HalfDim.z;

// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
		if(BfSt.x[ii]==0.) BfSt.x[ii] = radCR.AbsRandMagnitude(HalfDim.x);
		if(BfSt.y[ii]==0.) BfSt.y[ii] = radCR.AbsRandMagnitude(HalfDim.y);
		if(BfSt.z[ii]==0.) BfSt.z[ii] = radCR.AbsRandMagnitude(HalfDim.z);
	}

	double x0 = BfSt.x[0], x1 = BfSt.x[1];
	double y0 = BfSt.y[0], y1 = BfSt.y[1];
	double z0 = BfSt.z[0], z1 = BfSt.z[1];

// [FieldKey.M_ query REMOVED 2026-06-28: a current block has no magnetization]
	if(FieldPtr->FieldKey.J_) //OC061008
	{
		if(!J_is_Zero) //to drop this test?
		{
			if((x0*x1<0) && (y0*y1<0) && (z0*z1<0)) FieldPtr->J += J;
		}
	}
	if((!FieldPtr->FieldKey.B_) && (!FieldPtr->FieldKey.H_) && (!FieldPtr->FieldKey.A_) && (!FieldPtr->FieldKey.Phi_) && (!FieldPtr->FieldKey.PreRelax_)) return;

	double x0e2 = x0*x0, x1e2 = x1*x1;
	double y0e2 = y0*y0, y1e2 = y1*y1;
	double z0e2 = z0*z0, z1e2 = z1*z1;

	double D000 = sqrt(x0e2+y0e2+z0e2);	
	double D100 = sqrt(x1e2+y0e2+z0e2);	
	double D010 = sqrt(x0e2+y1e2+z0e2);	
	double D110 = sqrt(x1e2+y1e2+z0e2);	
	double D001 = sqrt(x0e2+y0e2+z1e2);	
	double D101 = sqrt(x1e2+y0e2+z1e2);	
	double D011 = sqrt(x0e2+y1e2+z1e2);	
	double D111 = sqrt(x1e2+y1e2+z1e2);	

	const double Pi = 3.141592653589793238;
	double PiMult1, PiMult2, PiMult3;
	PiMult1 = PiMult2 = PiMult3 = 0.;

	T0.x = atan(TransAtans(TransAtans(y0*z0/(x0*D000), -y0*z1/(x0*D001), PiMult1), 
						   TransAtans(-y1*z0/(x0*D010), y1*z1/(x0*D011), PiMult2), PiMult3))+Pi*(PiMult1+PiMult2+PiMult3);
	T1.x = atan(TransAtans(TransAtans(-y0*z0/(x1*D100), y0*z1/(x1*D101), PiMult1), 
						   TransAtans(y1*z0/(x1*D110), -y1*z1/(x1*D111), PiMult2), PiMult3))+Pi*(PiMult1+PiMult2+PiMult3);
	T0.y = atan(TransAtans(TransAtans(x0*z0/(y0*D000), -x0*z1/(y0*D001), PiMult1), 
						   TransAtans(-x1*z0/(y0*D100), x1*z1/(y0*D101), PiMult2), PiMult3))+Pi*(PiMult1+PiMult2+PiMult3);
	T1.y = atan(TransAtans(TransAtans(-x0*z0/(y1*D010), x0*z1/(y1*D011), PiMult1), 
						   TransAtans(x1*z0/(y1*D110), -x1*z1/(y1*D111), PiMult2), PiMult3))+Pi*(PiMult1+PiMult2+PiMult3);
	T0.z = atan(TransAtans(TransAtans(x0*y0/(z0*D000), -x1*y0/(z0*D100), PiMult1), 
						   TransAtans(-x0*y1/(z0*D010), x1*y1/(z0*D110), PiMult2), PiMult3))+Pi*(PiMult1+PiMult2+PiMult3);
	T1.z = atan(TransAtans(TransAtans(-x0*y0/(z1*D001), x1*y0/(z1*D101), PiMult1), 
						   TransAtans(x0*y1/(z1*D011), -x1*y1/(z1*D111), PiMult2), PiMult3))+Pi*(PiMult1+PiMult2+PiMult3);

	double AbsRandD000 = 10.*radCR.AbsRandMagnitude(D000);
	double AbsRandD010 = 10.*radCR.AbsRandMagnitude(D010);
	double AbsRandD001 = 10.*radCR.AbsRandMagnitude(D001);
	double AbsRandD011 = 10.*radCR.AbsRandMagnitude(D011);
	double AbsRandD100 = 10.*radCR.AbsRandMagnitude(D100);
	double AbsRandD110 = 10.*radCR.AbsRandMagnitude(D110);
	double AbsRandD101 = 10.*radCR.AbsRandMagnitude(D101);
	double AbsRandD111 = 10.*radCR.AbsRandMagnitude(D111);

	double z0plD100 = z0+D100; if(z0plD100 < AbsRandD100) z0plD100 = 0.5*(x1e2 + y0e2)/Abs(z0);
	double z1plD101 = z1+D101; if(z1plD101 < AbsRandD101) z1plD101 = 0.5*(x1e2 + y0e2)/Abs(z1);
	double z1plD001 = z1+D001; if(z1plD001 < AbsRandD001) z1plD001 = 0.5*(x0e2 + y0e2)/Abs(z1);
	double z0plD000 = z0+D000; if(z0plD000 < AbsRandD000) z0plD000 = 0.5*(x0e2 + y0e2)/Abs(z0);
	double z0plD010 = z0+D010; if(z0plD010 < AbsRandD010) z0plD010 = 0.5*(x0e2 + y1e2)/Abs(z0);
	double z1plD011 = z1+D011; if(z1plD011 < AbsRandD011) z1plD011 = 0.5*(x0e2 + y1e2)/Abs(z1);
	double z1plD111 = z1+D111; if(z1plD111 < AbsRandD111) z1plD111 = 0.5*(x1e2 + y1e2)/Abs(z1);
	double z0plD110 = z0+D110; if(z0plD110 < AbsRandD110) z0plD110 = 0.5*(x1e2 + y1e2)/Abs(z0);

	double y0plD100 = y0+D100; if(y0plD100 < AbsRandD100) y0plD100 = 0.5*(x1e2 + z0e2)/Abs(y0);
	double y1plD110 = y1+D110; if(y1plD110 < AbsRandD110) y1plD110 = 0.5*(x1e2 + z0e2)/Abs(y1);
	double y1plD010 = y1+D010; if(y1plD010 < AbsRandD010) y1plD010 = 0.5*(x0e2 + z0e2)/Abs(y1);
	double y0plD000 = y0+D000; if(y0plD000 < AbsRandD000) y0plD000 = 0.5*(x0e2 + z0e2)/Abs(y0);
	double y0plD001 = y0+D001; if(y0plD001 < AbsRandD001) y0plD001 = 0.5*(x0e2 + z1e2)/Abs(y0);
	double y1plD011 = y1+D011; if(y1plD011 < AbsRandD011) y1plD011 = 0.5*(x0e2 + z1e2)/Abs(y1);
	double y1plD111 = y1+D111; if(y1plD111 < AbsRandD111) y1plD111 = 0.5*(x1e2 + z1e2)/Abs(y1);
	double y0plD101 = y0+D101; if(y0plD101 < AbsRandD101) y0plD101 = 0.5*(x1e2 + z1e2)/Abs(y0);

	double x0plD010 = x0+D010; if(x0plD010 < AbsRandD010) x0plD010 = 0.5*(y1e2 + z0e2)/Abs(x0);
	double x1plD110 = x1+D110; if(x1plD110 < AbsRandD110) x1plD110 = 0.5*(y1e2 + z0e2)/Abs(x1);
	double x1plD100 = x1+D100; if(x1plD100 < AbsRandD100) x1plD100 = 0.5*(y0e2 + z0e2)/Abs(x1);
	double x0plD000 = x0+D000; if(x0plD000 < AbsRandD000) x0plD000 = 0.5*(y0e2 + z0e2)/Abs(x0);
	double x0plD001 = x0+D001; if(x0plD001 < AbsRandD001) x0plD001 = 0.5*(y0e2 + z1e2)/Abs(x0);
	double x1plD101 = x1+D101; if(x1plD101 < AbsRandD101) x1plD101 = 0.5*(y0e2 + z1e2)/Abs(x1);
	double x1plD111 = x1+D111; if(x1plD111 < AbsRandD111) x1plD111 = 0.5*(y1e2 + z1e2)/Abs(x1);
	double x0plD011 = x0+D011; if(x0plD011 < AbsRandD011) x0plD011 = 0.5*(y1e2 + z1e2)/Abs(x0);

	double z0plD100_di_z1plD101 = z0plD100/z1plD101;
	double z1plD001_di_z0plD000 = z1plD001/z0plD000;
	double z0plD010_di_z1plD011 = z0plD010/z1plD011;
	double z1plD111_di_z0plD110 = z1plD111/z0plD110;
	double y0plD100_di_y1plD110 = y0plD100/y1plD110;
	double y1plD010_di_y0plD000 = y1plD010/y0plD000;
	double y0plD001_di_y1plD011 = y0plD001/y1plD011;
	double y1plD111_di_y0plD101 = y1plD111/y0plD101;
	double x0plD010_di_x1plD110 = x0plD010/x1plD110;
	double x1plD100_di_x0plD000 = x1plD100/x0plD000;
	double x0plD001_di_x1plD101 = x0plD001/x1plD101;
	double x1plD111_di_x0plD011 = x1plD111/x0plD011;

	const double dConst2 = 1./4./Pi;  // For Phi: H = -grad(Phi), so Phi uses 1/(4*pi)
	// For A: B = curl(A), so A uses mu_0/(4*pi) = 1e-7 H/m
	// Radia uses SI units (B in Tesla, H in A/m), so A should be in T*m
	const double dConstA = 1.0e-7;  // mu_0/(4*pi)
	// ConstForJ: Biot-Savart constant for current density
	// B = (mu_0/4*pi) * integral(J x r / r^3) dV
	// For SI units (J in A/m^2, B in T): ConstForJ = mu_0/(4*pi) = 1e-7
	const double ConstForJ = 1.0e-7;

	double ln_z0plD100_di_z1plD101, ln_z1plD001_di_z0plD000, ln_z0plD010_di_z1plD011, ln_z1plD111_di_z0plD110,
		   ln_y0plD100_di_y1plD110, ln_y1plD010_di_y0plD000, ln_y0plD001_di_y1plD011, ln_y1plD111_di_y0plD101,
		   ln_x0plD010_di_x1plD110, ln_x1plD100_di_x0plD000, ln_x0plD001_di_x1plD101, ln_x1plD111_di_x0plD011;

	if(FieldPtr->FieldKey.A_ || FieldPtr->FieldKey.Phi_ || !J_is_Zero)
	{
		ln_z0plD100_di_z1plD101 = log(z0plD100_di_z1plD101);
		ln_z1plD001_di_z0plD000 = log(z1plD001_di_z0plD000);
		ln_z0plD010_di_z1plD011 = log(z0plD010_di_z1plD011);
		ln_z1plD111_di_z0plD110 = log(z1plD111_di_z0plD110);
		ln_y0plD100_di_y1plD110 = log(y0plD100_di_y1plD110);
		ln_y1plD010_di_y0plD000 = log(y1plD010_di_y0plD000);
		ln_y0plD001_di_y1plD011 = log(y0plD001_di_y1plD011);
		ln_y1plD111_di_y0plD101 = log(y1plD111_di_y0plD101);
		ln_x0plD010_di_x1plD110 = log(x0plD010_di_x1plD110);
		ln_x1plD100_di_x0plD000 = log(x1plD100_di_x0plD000);
		ln_x0plD001_di_x1plD101 = log(x0plD001_di_x1plD101);
		ln_x1plD111_di_x0plD011 = log(x1plD111_di_x0plD011);

		TVector3d BufVect(x0*T0.x + x1*T1.x
						 +y0*(ln_z0plD100_di_z1plD101+ln_z1plD001_di_z0plD000)
						 +y1*(ln_z0plD010_di_z1plD011+ln_z1plD111_di_z0plD110)
						 +z0*(ln_y0plD100_di_y1plD110+ln_y1plD010_di_y0plD000)
						 +z1*(ln_y0plD001_di_y1plD011+ln_y1plD111_di_y0plD101),
						  y0*T0.y + y1*T1.y
						 +x0*(ln_z0plD010_di_z1plD011+ln_z1plD001_di_z0plD000)
						 +x1*(ln_z0plD100_di_z1plD101+ln_z1plD111_di_z0plD110)
						 +z0*(ln_x0plD010_di_x1plD110+ln_x1plD100_di_x0plD000)
						 +z1*(ln_x0plD001_di_x1plD101+ln_x1plD111_di_x0plD011),
						  z0*T0.z + z1*T1.z
						 +y0*(ln_x0plD001_di_x1plD101+ln_x1plD100_di_x0plD000)
						 +y1*(ln_x0plD010_di_x1plD110+ln_x1plD111_di_x0plD011)
						 +x0*(ln_y0plD001_di_y1plD011+ln_y1plD010_di_y0plD000)
						 +x1*(ln_y0plD100_di_y1plD110+ln_y1plD111_di_y0plD101));
		{
			if(FieldPtr->FieldKey.A_)
			{
				// Vector potential A from uniform current density J:
				// A = (mu_0 / 4*pi) * J * integral(1/|r-r'|) dV'
				//
				// For a rectangular parallelepiped, the integral has the form (Durand):
				// I = -sum over 8 corners of (+/-) * f(X,Y,Z)
				// where f(X,Y,Z) = X*Y*asinh(Z/sqrt(X^2+Y^2)) + Y*Z*asinh(X/sqrt(Y^2+Z^2)) + Z*X*asinh(Y/sqrt(Z^2+X^2))
				//                  - X^2/2*atan(YZ/XR) - Y^2/2*atan(ZX/YR) - Z^2/2*atan(XY/ZR)
				// with R = sqrt(X^2 + Y^2 + Z^2)
				//
				// Note: asinh(z/sqrt(x^2+y^2)) = ln((z+R)/sqrt(x^2+y^2)), which differs from ln(|z|+R)
				//
				// Corners are at (x0,y0,z0), (x1,y0,z0), (x0,y1,z0), (x1,y1,z0),
				//                (x0,y0,z1), (x1,y0,z1), (x0,y1,z1), (x1,y1,z1)
				// with signs: +, -, -, +, -, +, +, -

				// Small epsilon to avoid division by zero
				const double eps = 1e-30;

				// Helper function macro for asinh terms
				// asinh(z/sqrt(x^2+y^2)) where we need to handle the sqrt carefully
				#define ASINH_TERM(x,y,z) ((x)*(y)*asinh((z)/(sqrt((x)*(x)+(y)*(y)+eps))))
				#define ATAN_TERM(x,y,z,R) ((x)*(x)/2.*atan((y)*(z)/((x)*(R)+eps)))

				// Compute f(X,Y,Z) for each corner using asinh form
				// f = X*Y*asinh(Z/sqrt(X^2+Y^2)) + Y*Z*asinh(X/sqrt(Y^2+Z^2)) + Z*X*asinh(Y/sqrt(Z^2+X^2))
				//     - X^2/2*atan(YZ/XR) - Y^2/2*atan(ZX/YR) - Z^2/2*atan(XY/ZR)

				// Corner 000: (x0, y0, z0), sign = +1
				double f000 = ASINH_TERM(x0,y0,z0) + ASINH_TERM(y0,z0,x0) + ASINH_TERM(z0,x0,y0)
				            - ATAN_TERM(x0,y0,z0,D000) - ATAN_TERM(y0,z0,x0,D000) - ATAN_TERM(z0,x0,y0,D000);

				// Corner 100: (x1, y0, z0), sign = -1
				double f100 = ASINH_TERM(x1,y0,z0) + ASINH_TERM(y0,z0,x1) + ASINH_TERM(z0,x1,y0)
				            - ATAN_TERM(x1,y0,z0,D100) - ATAN_TERM(y0,z0,x1,D100) - ATAN_TERM(z0,x1,y0,D100);

				// Corner 010: (x0, y1, z0), sign = -1
				double f010 = ASINH_TERM(x0,y1,z0) + ASINH_TERM(y1,z0,x0) + ASINH_TERM(z0,x0,y1)
				            - ATAN_TERM(x0,y1,z0,D010) - ATAN_TERM(y1,z0,x0,D010) - ATAN_TERM(z0,x0,y1,D010);

				// Corner 110: (x1, y1, z0), sign = +1
				double f110 = ASINH_TERM(x1,y1,z0) + ASINH_TERM(y1,z0,x1) + ASINH_TERM(z0,x1,y1)
				            - ATAN_TERM(x1,y1,z0,D110) - ATAN_TERM(y1,z0,x1,D110) - ATAN_TERM(z0,x1,y1,D110);

				// Corner 001: (x0, y0, z1), sign = -1
				double f001 = ASINH_TERM(x0,y0,z1) + ASINH_TERM(y0,z1,x0) + ASINH_TERM(z1,x0,y0)
				            - ATAN_TERM(x0,y0,z1,D001) - ATAN_TERM(y0,z1,x0,D001) - ATAN_TERM(z1,x0,y0,D001);

				// Corner 101: (x1, y0, z1), sign = +1
				double f101 = ASINH_TERM(x1,y0,z1) + ASINH_TERM(y0,z1,x1) + ASINH_TERM(z1,x1,y0)
				            - ATAN_TERM(x1,y0,z1,D101) - ATAN_TERM(y0,z1,x1,D101) - ATAN_TERM(z1,x1,y0,D101);

				// Corner 011: (x0, y1, z1), sign = +1
				double f011 = ASINH_TERM(x0,y1,z1) + ASINH_TERM(y1,z1,x0) + ASINH_TERM(z1,x0,y1)
				            - ATAN_TERM(x0,y1,z1,D011) - ATAN_TERM(y1,z1,x0,D011) - ATAN_TERM(z1,x0,y1,D011);

				// Corner 111: (x1, y1, z1), sign = -1
				double f111 = ASINH_TERM(x1,y1,z1) + ASINH_TERM(y1,z1,x1) + ASINH_TERM(z1,x1,y1)
				            - ATAN_TERM(x1,y1,z1,D111) - ATAN_TERM(y1,z1,x1,D111) - ATAN_TERM(z1,x1,y1,D111);

				#undef ASINH_TERM
				#undef ATAN_TERM

				// Sum with alternating signs, then negate (from derivation matching scipy)
				// I = -(f000 - f100 - f010 + f110 - f001 + f101 + f011 - f111)
				double scalarIntegral = -(f000 - f100 - f010 + f110 - f001 + f101 + f011 - f111);

				// A = (mu_0 / 4*pi) * J * scalarIntegral
				FieldPtr->A += ConstForJ * scalarIntegral * J;
			}
			if(FieldPtr->FieldKey.B_ || FieldPtr->FieldKey.H_)
			{
				// B = (mu_0/4*pi) * (J x BufVect) in Tesla
				// But Radia stores B internally as B/mu_0 (in A/m), and OutFieldCompRes
				// multiplies by mu_0 on output. So we need to divide by mu_0 here.
				// ConstForJ = mu_0/(4*pi) = 1e-7, and 1/mu_0 = 1/(4*pi*1e-7) = 1/(ConstForJ * 4)
				// So: B_internal = (1/mu_0) * ConstForJ * (J x BufVect) = (J x BufVect) / (4*pi)
				const double dConst2 = 1./4./Pi;  // 1/(4*pi)
				TVector3d BufForB(J.y*BufVect.z-J.z*BufVect.y,
								  J.z*BufVect.x-J.x*BufVect.z,
								  J.x*BufVect.y-J.y*BufVect.x);
				FieldPtr->B += dConst2*BufForB;
				FieldPtr->H += dConst2*BufForB;
			}
		}
	}
	
// [magnetization demag-tensor B/H + PreRelax block REMOVED 2026-06-28: current-only]
}

//-------------------------------------------------------------------------

// radTRecCur::B_compMultipole (magnetization multipole expansion) REMOVED 2026-06-28 (current-only)

//-------------------------------------------------------------------------

void radTRecCur::B_intComp(radTField* FieldPtr)
{
	if(FieldPtr->FieldKey.FinInt_) { B_intCompFinNum(FieldPtr); return;}

// An analytical algorithm for infinite Field Integrals:
	TVector3d CenPo_mi_StPo = CentrPoint - FieldPtr->P;
	TVector3d HalfDim = 0.5 * Dimensions;
	TVector3d P1 = CenPo_mi_StPo - HalfDim;
	TVector3d P2 = CenPo_mi_StPo + HalfDim;
	TVector3d V = FieldPtr->NextP - FieldPtr->P;

	double ModV = sqrt(V.x*V.x + V.y*V.y + V.z*V.z);
	V.x /= ModV;  V.y /= ModV;  V.z /= ModV;

	const double Pi = 3.141592653589793238;
	const double ZeroToler = 1.E-06; // Relative tolerance to switch to special cases
	const double SmallestRelTolerV = 1.E-12; // Relative tolerance to repair trapping V.i to zero at general case

	double AbsRandX = radCR.AbsRandMagnitude(CentrPoint.x);
	double AbsRandY = radCR.AbsRandMagnitude(CentrPoint.y);
	double AbsRandZ = radCR.AbsRandMagnitude(CentrPoint.z);

	TMatrix3d F;
	TVector3d G;

	short J_is_Zero = 1;
	if(J.x!=0. || J.y!=0. || J.z!=0.) J_is_Zero = 0;

// Tests for special cases:
	double AbsVx = Abs(V.x);
	double AbsVy = Abs(V.y);
	double AbsVz = Abs(V.z);
	if(AbsVx<ZeroToler && AbsVy<ZeroToler)
	{
// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
		if(P1.x==0.) P1.x = AbsRandX;
		if(P1.y==0.) P1.y = AbsRandY;
		if(P2.x==0.) P2.x = AbsRandX;
		if(P2.y==0.) P2.y = AbsRandY;

		B_intUtilSpecCaseZeroVxVy(P1, P2, J_is_Zero, F, G);
		goto FinalDefinitionOfFieldIntegrals;
	}
	if(AbsVx<ZeroToler && AbsVz<ZeroToler) 
	{
// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
		if(P1.x==0.) P1.x = AbsRandX;
		if(P1.z==0.) P1.z = AbsRandZ;
		if(P2.x==0.) P2.x = AbsRandX;
		if(P2.z==0.) P2.z = AbsRandZ;

		TVector3d LocP1(P1.x, P1.z, P1.y), LocP2(P2.x, P2.z, P2.y), LocG;
		TMatrix3d LocF;
		B_intUtilSpecCaseZeroVxVy(LocP1, LocP2, J_is_Zero, LocF, LocG);
		TVector3d& F_str0 = F.Str0;
		TVector3d& F_str1 = F.Str1;
		TVector3d& F_str2 = F.Str2;
		TVector3d& LocF_str0 = LocF.Str0;
		TVector3d& LocF_str1 = LocF.Str1;
		TVector3d& LocF_str2 = LocF.Str2;
		F_str0.x = LocF_str0.y; F_str0.y = LocF_str0.x; F_str0.z = LocF_str0.z;
		F_str1.x = LocF_str2.y; F_str1.y = LocF_str2.x; F_str1.z = LocF_str2.z;
		F_str2.x = LocF_str1.y; F_str2.y = LocF_str1.x; F_str2.z = LocF_str1.z;
		G.x = LocG.x; G.y = LocG.z; G.z = LocG.y;
		goto FinalDefinitionOfFieldIntegrals;
	}
	if(AbsVy<ZeroToler && AbsVz<ZeroToler) 
	{
// Artificial shift of an observation point a bit right of the block's border
// if the point is exactly on the boarder (to avoid "divide by zero" error):
		if(P1.y==0.) P1.y = AbsRandY;
		if(P1.z==0.) P1.z = AbsRandZ;
		if(P2.y==0.) P2.y = AbsRandY;
		if(P2.z==0.) P2.z = AbsRandZ;

		TVector3d LocP1(P1.z, P1.y, P1.x), LocP2(P2.z, P2.y, P2.x), LocG;
		TMatrix3d LocF;
		B_intUtilSpecCaseZeroVxVy(LocP1, LocP2, J_is_Zero, LocF, LocG);
		TVector3d& F_str0 = F.Str0;
		TVector3d& F_str1 = F.Str1;
		TVector3d& F_str2 = F.Str2;
		TVector3d& LocF_str0 = LocF.Str0;
		TVector3d& LocF_str1 = LocF.Str1;
		TVector3d& LocF_str2 = LocF.Str2;
		F_str0.x = LocF_str2.z; F_str0.y = LocF_str2.y; F_str0.z = LocF_str2.x;
		F_str1.x = LocF_str1.z; F_str1.y = LocF_str1.y; F_str1.z = LocF_str1.x;
		F_str2.x = LocF_str0.z; F_str2.y = LocF_str0.y; F_str2.z = LocF_str0.x;
		G.x = LocG.z; G.y = LocG.y; G.z = LocG.x;
		goto FinalDefinitionOfFieldIntegrals;
	}

	{
// Prevent trapping each of V.x,V.y,V.z to Zero separately???
		if(AbsVx<SmallestRelTolerV) V.x = SmallestRelTolerV;
		if(AbsVy<SmallestRelTolerV) V.y = SmallestRelTolerV;
		if(AbsVz<SmallestRelTolerV) V.z = SmallestRelTolerV;
//GeneralCaseStart:
		double vxvx = V.x*V.x;
		double vyvy = V.y*V.y;
		double vzvz = V.z*V.z;
		double vxvx_p_vyvy = vxvx + vyvy;
		double vxvx_p_vzvz = vxvx + vzvz;
		double vyvy_p_vzvz = vyvy + vzvz;
		double vzvz_m_vyvy = vzvz - vyvy;
		double vzvz_m_vxvx = vzvz - vxvx;
		double vxvx_m_vyvy = vxvx - vyvy;
		double vyvy_p_vzvz_p_2vxvx = vyvy_p_vzvz + 2.*vxvx;
		double vxvx_p_vzvz_p_2vyvy = vxvx_p_vzvz + 2.*vyvy;
		double vxvx_p_vyvy_p_2vzvz = vxvx_p_vyvy + 2.*vzvz;
		double One_d_vxvxpvzvz = 1./vxvx_p_vzvz;
		double One_d_vxvxpvyvy = 1./vxvx_p_vyvy;
		double One_d_vyvypvzvz = 1./vyvy_p_vzvz;
		double One_d_vxvxpvzvz_d_vyvypvzvz = One_d_vxvxpvzvz*One_d_vyvypvzvz;
		double One_d_vxvxpvyvy_d_vxvxpvzvz = One_d_vxvxpvzvz*One_d_vxvxpvyvy;
		double One_d_vxvxpvyvy_d_vyvypvzvz = One_d_vxvxpvyvy*One_d_vyvypvzvz;
		double vxx1 = V.x*P1.x;
		double vxy1 = V.x*P1.y;
		double vxz1 = V.x*P1.z;
		double vyx1 = V.y*P1.x;
		double vyy1 = V.y*P1.y;
		double vyz1 = V.y*P1.z;
		double vzx1 = V.z*P1.x;
		double vzy1 = V.z*P1.y;
		double vzz1 = V.z*P1.z;
		double vxx2 = V.x*P2.x;
		double vxy2 = V.x*P2.y;
		double vxz2 = V.x*P2.z;
		double vyx2 = V.y*P2.x;
		double vyy2 = V.y*P2.y;
		double vyz2 = V.y*P2.z;
		double vzx2 = V.z*P2.x;
		double vzy2 = V.z*P2.y;
		double vzz2 = V.z*P2.z;   		  // Remove this ???
		double vzy1_m_vyz1 = vzy1 - vyz1; if(vzy1_m_vyz1==0.) vzy1_m_vyz1 = SmallestRelTolerV*AbsRandY;
		double vzx1_m_vxz1 = vzx1 - vxz1; if(vzx1_m_vxz1==0.) vzx1_m_vxz1 = SmallestRelTolerV*AbsRandX;
		double vxy1_m_vyx1 = vxy1 - vyx1; if(vxy1_m_vyx1==0.) vxy1_m_vyx1 = SmallestRelTolerV*AbsRandY;
		double vzy2_m_vyz1 = vzy2 - vyz1; if(vzy2_m_vyz1==0.) vzy2_m_vyz1 = SmallestRelTolerV*AbsRandY;
		double vzx1_m_vxz2 = vzx1 - vxz2; if(vzx1_m_vxz2==0.) vzx1_m_vxz2 = SmallestRelTolerV*AbsRandX;
		double vxy1_m_vyx2 = vxy1 - vyx2; if(vxy1_m_vyx2==0.) vxy1_m_vyx2 = SmallestRelTolerV*AbsRandY;
		double vzy1_m_vyz2 = vzy1 - vyz2; if(vzy1_m_vyz2==0.) vzy1_m_vyz2 = SmallestRelTolerV*AbsRandY;
		double vzx2_m_vxz1 = vzx2 - vxz1; if(vzx2_m_vxz1==0.) vzx2_m_vxz1 = SmallestRelTolerV*AbsRandX;
		double vxy2_m_vyx1 = vxy2 - vyx1; if(vxy2_m_vyx1==0.) vxy2_m_vyx1 = SmallestRelTolerV*AbsRandY;
		double vzy2_m_vyz2 = vzy2 - vyz2; if(vzy2_m_vyz2==0.) vzy2_m_vyz2 = SmallestRelTolerV*AbsRandY;
		double vzx2_m_vxz2 = vzx2 - vxz2; if(vzx2_m_vxz2==0.) vzx2_m_vxz2 = SmallestRelTolerV*AbsRandX;
		double vxy2_m_vyx2 = vxy2 - vyx2; if(vxy2_m_vyx2==0.) vxy2_m_vyx2 = SmallestRelTolerV*AbsRandY;
		double ArgAtanX111 = (V.z*vzx1_m_vxz1 - V.y*vxy1_m_vyx1)/vzy1_m_vyz1;
		double ArgAtanX211 = (V.z*vzx2_m_vxz1 - V.y*vxy1_m_vyx2)/vzy1_m_vyz1;
		double ArgAtanX121 = (V.z*vzx1_m_vxz1 - V.y*vxy2_m_vyx1)/vzy2_m_vyz1;
		double ArgAtanX221 = (V.z*vzx2_m_vxz1 - V.y*vxy2_m_vyx2)/vzy2_m_vyz1;
		double ArgAtanX112 = (V.z*vzx1_m_vxz2 - V.y*vxy1_m_vyx1)/vzy1_m_vyz2;
		double ArgAtanX212 = (V.z*vzx2_m_vxz2 - V.y*vxy1_m_vyx2)/vzy1_m_vyz2;
		double ArgAtanX122 = (V.z*vzx1_m_vxz2 - V.y*vxy2_m_vyx1)/vzy2_m_vyz2;
		double ArgAtanX222 = (V.z*vzx2_m_vxz2 - V.y*vxy2_m_vyx2)/vzy2_m_vyz2;
		double ArgAtanY111 = (V.z*vzy1_m_vyz1 + V.x*vxy1_m_vyx1)/vzx1_m_vxz1;
		double ArgAtanY211 = (V.z*vzy1_m_vyz1 + V.x*vxy1_m_vyx2)/vzx2_m_vxz1;
		double ArgAtanY121 = (V.z*vzy2_m_vyz1 + V.x*vxy2_m_vyx1)/vzx1_m_vxz1;
		double ArgAtanY221 = (V.z*vzy2_m_vyz1 + V.x*vxy2_m_vyx2)/vzx2_m_vxz1;
		double ArgAtanY112 = (V.z*vzy1_m_vyz2 + V.x*vxy1_m_vyx1)/vzx1_m_vxz2;
		double ArgAtanY212 = (V.z*vzy1_m_vyz2 + V.x*vxy1_m_vyx2)/vzx2_m_vxz2;
		double ArgAtanY122 = (V.z*vzy2_m_vyz2 + V.x*vxy2_m_vyx1)/vzx1_m_vxz2;
		double ArgAtanY222 = (V.z*vzy2_m_vyz2 + V.x*vxy2_m_vyx2)/vzx2_m_vxz2;
		double ArgAtanZ111 = (V.x*vzx1_m_vxz1 + V.y*vzy1_m_vyz1)/vxy1_m_vyx1;
		double ArgAtanZ211 = (V.x*vzx2_m_vxz1 + V.y*vzy1_m_vyz1)/vxy1_m_vyx2;
		double ArgAtanZ121 = (V.x*vzx1_m_vxz1 + V.y*vzy2_m_vyz1)/vxy2_m_vyx1;
		double ArgAtanZ221 = (V.x*vzx2_m_vxz1 + V.y*vzy2_m_vyz1)/vxy2_m_vyx2;
		double ArgAtanZ112 = (V.x*vzx1_m_vxz2 + V.y*vzy1_m_vyz2)/vxy1_m_vyx1;
		double ArgAtanZ212 = (V.x*vzx2_m_vxz2 + V.y*vzy1_m_vyz2)/vxy1_m_vyx2;
		double ArgAtanZ122 = (V.x*vzx1_m_vxz2 + V.y*vzy2_m_vyz2)/vxy2_m_vyx1;
		double ArgAtanZ222 = (V.x*vzx2_m_vxz2 + V.y*vzy2_m_vyz2)/vxy2_m_vyx2;
		double PiMult = 0.;
		double SumAtanXy1z1 = atan(TransAtans(ArgAtanX111, -ArgAtanX211, PiMult)) + Pi*PiMult;
		double SumAtanXy2z1 = atan(TransAtans(ArgAtanX121, -ArgAtanX221, PiMult)) + Pi*PiMult;
		double SumAtanXy1z2 = atan(TransAtans(ArgAtanX112, -ArgAtanX212, PiMult)) + Pi*PiMult;
		double SumAtanXy2z2 = atan(TransAtans(ArgAtanX122, -ArgAtanX222, PiMult)) + Pi*PiMult;
		double SumAtanYx1z1 = atan(TransAtans(ArgAtanY111, -ArgAtanY121, PiMult)) + Pi*PiMult;
		double SumAtanYx2z1 = atan(TransAtans(ArgAtanY211, -ArgAtanY221, PiMult)) + Pi*PiMult;
		double SumAtanYx1z2 = atan(TransAtans(ArgAtanY112, -ArgAtanY122, PiMult)) + Pi*PiMult;
		double SumAtanYx2z2 = atan(TransAtans(ArgAtanY212, -ArgAtanY222, PiMult)) + Pi*PiMult;
		double SumAtanZx1y1 = atan(TransAtans(ArgAtanZ111, -ArgAtanZ112, PiMult)) + Pi*PiMult;
		double SumAtanZx2y1 = atan(TransAtans(ArgAtanZ211, -ArgAtanZ212, PiMult)) + Pi*PiMult;
		double SumAtanZx1y2 = atan(TransAtans(ArgAtanZ121, -ArgAtanZ122, PiMult)) + Pi*PiMult;
		double SumAtanZx2y2 = atan(TransAtans(ArgAtanZ221, -ArgAtanZ222, PiMult)) + Pi*PiMult;
		double vxy1mvyx1_mu_vxy1mvyx1 = vxy1_m_vyx1*vxy1_m_vyx1;
		double vxy1mvyx2_mu_vxy1mvyx2 = vxy1_m_vyx2*vxy1_m_vyx2;
		double vxy2mvyx1_mu_vxy2mvyx1 = vxy2_m_vyx1*vxy2_m_vyx1;
		double vxy2mvyx2_mu_vxy2mvyx2 = vxy2_m_vyx2*vxy2_m_vyx2;
		double vzy1mvyz1_mu_vzy1mvyz1 = vzy1_m_vyz1*vzy1_m_vyz1;
		double vzy2mvyz1_mu_vzy2mvyz1 = vzy2_m_vyz1*vzy2_m_vyz1;
		double vzy1mvyz2_mu_vzy1mvyz2 = vzy1_m_vyz2*vzy1_m_vyz2;
		double vzy2mvyz2_mu_vzy2mvyz2 = vzy2_m_vyz2*vzy2_m_vyz2;
		double vzx1mvxz1_mu_vzx1mvxz1 = vzx1_m_vxz1*vzx1_m_vxz1;
		double vzx2mvxz1_mu_vzx2mvxz1 = vzx2_m_vxz1*vzx2_m_vxz1;
		double vzx1mvxz2_mu_vzx1mvxz2 = vzx1_m_vxz2*vzx1_m_vxz2;
		double vzx2mvxz2_mu_vzx2mvxz2 = vzx2_m_vxz2*vzx2_m_vxz2;
		double Log111 = log(vxy1mvyx1_mu_vxy1mvyx1 + vzy1mvyz1_mu_vzy1mvyz1 + vzx1mvxz1_mu_vzx1mvxz1);
		double Log211 = log(vxy1mvyx2_mu_vxy1mvyx2 + vzy1mvyz1_mu_vzy1mvyz1 + vzx2mvxz1_mu_vzx2mvxz1);
		double Log121 = log(vxy2mvyx1_mu_vxy2mvyx1 + vzy2mvyz1_mu_vzy2mvyz1 + vzx1mvxz1_mu_vzx1mvxz1);
		double Log221 = log(vxy2mvyx2_mu_vxy2mvyx2 + vzy2mvyz1_mu_vzy2mvyz1 + vzx2mvxz1_mu_vzx2mvxz1);
		double Log112 = log(vxy1mvyx1_mu_vxy1mvyx1 + vzy1mvyz2_mu_vzy1mvyz2 + vzx1mvxz2_mu_vzx1mvxz2);
		double Log212 = log(vxy1mvyx2_mu_vxy1mvyx2 + vzy1mvyz2_mu_vzy1mvyz2 + vzx2mvxz2_mu_vzx2mvxz2);
		double Log122 = log(vxy2mvyx1_mu_vxy2mvyx1 + vzy2mvyz2_mu_vzy2mvyz2 + vzx1mvxz2_mu_vzx1mvxz2);
		double Log222 = log(vxy2mvyx2_mu_vxy2mvyx2 + vzy2mvyz2_mu_vzy2mvyz2 + vzx2mvxz2_mu_vzx2mvxz2);
	
		if(J_is_Zero)
		{
			double CommonLogTermZ = 0.5*(P1.z*(Log111-Log211+Log221-Log121) + P2.z*(Log212-Log112+Log122-Log222));
			double CommonLogTermX = 0.5*(P1.x*(Log111-Log121+Log122-Log112) + P2.x*(Log221-Log211+Log212-Log222));
			double CommonLogTermY = 0.5*(P1.y*(Log111-Log211+Log212-Log112) + P2.y*(Log221-Log121+Log122-Log222));
/* Fxyx*/	F.Str0.x = One_d_vxvxpvzvz*(vzx1_m_vxz1*SumAtanYx1z1 - vzx2_m_vxz1*SumAtanYx2z1 - vzx1_m_vxz2*SumAtanYx1z2 + vzx2_m_vxz2*SumAtanYx2z2)
					 + 0.5*V.y*One_d_vxvxpvzvz*((vzz1+vxx1)*(Log121-Log111) + (vzz1+vxx2)*(Log211-Log221) + (vzz2+vxx1)*(Log112-Log122) + (vzz2+vxx2)*(Log222-Log212))
					 + CommonLogTermY;
/* Fxzx*/	F.Str0.y = One_d_vxvxpvyvy*(-vxy1_m_vyx1*SumAtanZx1y1 + vxy1_m_vyx2*SumAtanZx2y1 + vxy2_m_vyx1*SumAtanZx1y2 - vxy2_m_vyx2*SumAtanZx2y2)
					 + 0.5*V.z*One_d_vxvxpvyvy*((vyy1+vxx1)*(Log112-Log111) + (vyy1+vxx2)*(Log211-Log212) + (vyy2+vxx1)*(Log121-Log122) + (vyy2+vxx2)*(Log222-Log221))
					 + CommonLogTermZ;
/* -Fyz0*/	F.Str0.z = One_d_vxvxpvyvy_d_vxvxpvzvz*((vzvz_m_vyvy*vxx1+vxvx_p_vzvz*vyy1)*SumAtanZx1y1 - (vzvz_m_vyvy*vxx1+vxvx_p_vzvz*vyy2)*SumAtanZx1y2 - (vzvz_m_vyvy*vxx2+vxvx_p_vzvz*vyy1)*SumAtanZx2y1 + (vzvz_m_vyvy*vxx2+vxvx_p_vzvz*vyy2)*SumAtanZx2y2)
					 + V.z*One_d_vxvxpvzvz*(P1.z*(SumAtanYx1z1-SumAtanYx2z1) + P2.z*(SumAtanYx2z2-SumAtanYx1z2))
					 - 0.5*V.z*One_d_vxvxpvyvy_d_vxvxpvzvz*((vyvy_p_vzvz_p_2vxvx*vyx1-vxvx_p_vzvz*vxy1)*(Log112-Log111) + (vyvy_p_vzvz_p_2vxvx*vyx1-vxvx_p_vzvz*vxy2)*(Log121-Log122) + (vyvy_p_vzvz_p_2vxvx*vyx2-vxvx_p_vzvz*vxy1)*(Log211-Log212) + (vyvy_p_vzvz_p_2vxvx*vyx2-vxvx_p_vzvz*vxy2)*(Log222-Log221))
					 - V.x*V.y*One_d_vxvxpvzvz*CommonLogTermZ
					 - (Pi*vzvz*One_d_vxvxpvzvz/V.x)*(P1.x*Step(vzx1/V.x-P1.z)*Step(P2.z-vzx1/V.x)*(Sign(vxy1_m_vyx1)-Sign(vxy2_m_vyx1)) + P2.x*Step(vzx2/V.x-P1.z)*Step(P2.z-vzx2/V.x)*(Sign(vxy2_m_vyx2)-Sign(vxy1_m_vyx2)));
/* Fxyy*/	F.Str1.x = One_d_vyvypvzvz*(vzy1_m_vyz1*SumAtanXy1z1 - vzy2_m_vyz1*SumAtanXy2z1 - vzy1_m_vyz2*SumAtanXy1z2 + vzy2_m_vyz2*SumAtanXy2z2)
					 + 0.5*V.x*One_d_vyvypvzvz*((vzz1+vyy1)*(Log211-Log111) + (vzz1+vyy2)*(Log121-Log221) + (vzz2+vyy1)*(Log112-Log212) + (vzz2+vyy2)*(Log222-Log122))
					 + CommonLogTermX;
/* -Fxz0*/	F.Str1.y =-One_d_vxvxpvyvy_d_vyvypvzvz*((vzvz_m_vxvx*vyy1+vyvy_p_vzvz*vxx1)*SumAtanZx1y1 - (vzvz_m_vxvx*vyy1+vyvy_p_vzvz*vxx2)*SumAtanZx2y1 - (vzvz_m_vxvx*vyy2+vyvy_p_vzvz*vxx1)*SumAtanZx1y2 + (vzvz_m_vxvx*vyy2+vyvy_p_vzvz*vxx2)*SumAtanZx2y2)
					 + V.z*One_d_vyvypvzvz*(P1.z*(SumAtanXy1z1-SumAtanXy2z1) + P2.z*(SumAtanXy2z2-SumAtanXy1z2))
					 - 0.5*V.z*One_d_vxvxpvyvy_d_vyvypvzvz*((vxvx_p_vzvz_p_2vyvy*vxy1-vyvy_p_vzvz*vyx1)*(Log112-Log111) + (vxvx_p_vzvz_p_2vyvy*vxy1-vyvy_p_vzvz*vyx2)*(Log211-Log212) + (vxvx_p_vzvz_p_2vyvy*vxy2-vyvy_p_vzvz*vyx1)*(Log121-Log122) + (vxvx_p_vzvz_p_2vyvy*vxy2-vyvy_p_vzvz*vyx2)*(Log222-Log221))
					 - V.x*V.y*One_d_vyvypvzvz*CommonLogTermZ
					 - (Pi*vzvz*One_d_vyvypvzvz/V.y)*(P1.y*Step(vzy1/V.y-P1.z)*Step(P2.z-vzy1/V.y)*(Sign(vxy1_m_vyx2)-Sign(vxy1_m_vyx1)) + P2.y*Step(vzy2/V.y-P1.z)*Step(P2.z-vzy2/V.y)*(Sign(vxy2_m_vyx1)-Sign(vxy2_m_vyx2)));
/* Fyzy*/	F.Str1.z = F.Str0.y;
/* -Fxy0*/	F.Str2.x =-One_d_vxvxpvzvz_d_vyvypvzvz*((vxvx_m_vyvy*vzz1+vxvx_p_vzvz*vyy1)*SumAtanXy1z1 - (vxvx_m_vyvy*vzz1+vxvx_p_vzvz*vyy2)*SumAtanXy2z1 - (vxvx_m_vyvy*vzz2+vxvx_p_vzvz*vyy1)*SumAtanXy1z2 + (vxvx_m_vyvy*vzz2+vxvx_p_vzvz*vyy2)*SumAtanXy2z2)
					 - V.x*One_d_vxvxpvzvz*(P1.x*(SumAtanYx1z1-SumAtanYx1z2) + P2.x*(SumAtanYx2z2-SumAtanYx2z1))
					 - 0.5*V.x*One_d_vxvxpvzvz_d_vyvypvzvz*((vxvx_p_vyvy_p_2vzvz*vyz1-vxvx_p_vzvz*vzy1)*(Log211-Log111) + (vxvx_p_vyvy_p_2vzvz*vyz1-vxvx_p_vzvz*vzy2)*(Log121-Log221) + (vxvx_p_vyvy_p_2vzvz*vyz2-vxvx_p_vzvz*vzy1)*(Log112-Log212) + (vxvx_p_vyvy_p_2vzvz*vyz2-vxvx_p_vzvz*vzy2)*(Log222-Log122))
					 - V.y*V.z*One_d_vxvxpvzvz*CommonLogTermX
					 - (Pi*vxvx*One_d_vxvxpvzvz/V.z)*(P1.z*Step(vxz1/V.z-P1.x)*Step(P2.x-vxz1/V.z)*(Sign(vzy1_m_vyz1)-Sign(vzy2_m_vyz1)) + P2.z*Step(vxz2/V.z-P1.x)*Step(P2.x-vxz2/V.z)*(Sign(vzy2_m_vyz2)-Sign(vzy1_m_vyz2)));
/* Fxzz*/	F.Str2.y = F.Str1.x;
/* Fyzz*/	F.Str2.z = F.Str0.x;
		}
		else
		{
			double vy_d_vxvxpvzvz_mu_vyvypvzvzp2vxvx = V.y*One_d_vxvxpvzvz*vyvy_p_vzvz_p_2vxvx;
			double vx_d_vyvypvzvz_mu_vxvxpvzvzp2vyvy = V.x*One_d_vyvypvzvz*vxvx_p_vzvz_p_2vyvy;
			double vy_d_vxvxpvzvz_mu_vxvxpvyvyp2vzvz = V.y*One_d_vxvxpvzvz*vxvx_p_vyvy_p_2vzvz;
			double x1x1 = P1.x*P1.x;
			double x2x2 = P2.x*P2.x;
			double y1y1 = P1.y*P1.y;
			double y2y2 = P2.y*P2.y;
			double z1z1 = P1.z*P1.z;
			double z2z2 = P2.z*P2.z;
			double TwoGx, TwoGy, TwoGz;

			TwoGx =-One_d_vxvxpvyvy_d_vxvxpvzvz*((vxvx_p_vzvz*(vxy1_m_vyx1-vyx1)*P1.y-vzvz_m_vyvy*vxx1*P1.x)*SumAtanZx1y1 - (vxvx_p_vzvz*(vxy2_m_vyx1-vyx1)*P2.y-vzvz_m_vyvy*vxx1*P1.x)*SumAtanZx1y2 - (vxvx_p_vzvz*(vxy1_m_vyx2-vyx2)*P1.y-vzvz_m_vyvy*vxx2*P2.x)*SumAtanZx2y1 + (vxvx_p_vzvz*(vxy2_m_vyx2-vyx2)*P2.y-vzvz_m_vyvy*vxx2*P2.x)*SumAtanZx2y2)
				  + One_d_vxvxpvzvz*(P1.z*(vzx1_m_vxz1+vzx1)*SumAtanYx1z1 - P2.z*(vzx1_m_vxz2+vzx1)*SumAtanYx1z2 - P1.z*(vzx2_m_vxz1+vzx2)*SumAtanYx2z1 + P2.z*(vzx2_m_vxz2+vzx2)*SumAtanYx2z2)
				  + 0.5*V.y*One_d_vxvxpvzvz*((2.*vxx1+vzz1)*P1.z*(Log121-Log111) + (2.*vxx1+vzz2)*P2.z*(Log112-Log122) + (2.*vxx2+vzz1)*P1.z*(Log211-Log221) + (2.*vxx2+vzz2)*P2.z*(Log222-Log212))
				  + 0.5*V.z*One_d_vxvxpvyvy*((vy_d_vxvxpvzvz_mu_vyvypvzvzp2vxvx*x1x1-(2.*vxx1+vyy1)*P1.y)*(Log111-Log112) + (-vy_d_vxvxpvzvz_mu_vyvypvzvzp2vxvx*x1x1+(2.*vxx1+vyy2)*P2.y)*(Log121-Log122) + (vy_d_vxvxpvzvz_mu_vyvypvzvzp2vxvx*x2x2-(2.*vxx2+vyy1)*P1.y)*(Log212-Log211) + (-vy_d_vxvxpvzvz_mu_vyvypvzvzp2vxvx*x2x2+(2.*vxx2+vyy2)*P2.y)*(Log222-Log221))
				  + P1.y*P1.z*(Log111-Log211) + P1.y*P2.z*(Log212-Log112) + P2.y*P1.z*(Log221-Log121) + P2.y*P2.z*(Log122-Log222)
				  + (Pi*V.z*V.z*One_d_vxvxpvzvz/V.x)*(x1x1*Step(vzx1/V.x-P1.z)*Step(P2.z-vzx1/V.x)*(Sign(vxy2_m_vyx1)-Sign(vxy1_m_vyx1)) + x2x2*Step(vzx2/V.x-P1.z)*Step(P2.z-vzx2/V.x)*(Sign(vxy1_m_vyx2)-Sign(vxy2_m_vyx2)));
			TwoGy = One_d_vxvxpvyvy_d_vyvypvzvz*(-(vyvy_p_vzvz*(vxy1_m_vyx1+vxy1)*P1.x+vzvz_m_vxvx*vyy1*P1.y)*SumAtanZx1y1 + (vyvy_p_vzvz*(vxy1_m_vyx2+vxy1)*P2.x+vzvz_m_vxvx*vyy1*P1.y)*SumAtanZx2y1 + (vyvy_p_vzvz*(vxy2_m_vyx1+vxy2)*P1.x+vzvz_m_vxvx*vyy2*P2.y)*SumAtanZx1y2 - (vyvy_p_vzvz*(vxy2_m_vyx2+vxy2)*P2.x+vzvz_m_vxvx*vyy2*P2.y)*SumAtanZx2y2)
				  + One_d_vyvypvzvz*(P1.z*(vzy1_m_vyz1+vzy1)*SumAtanXy1z1 - P2.z*(vzy1_m_vyz2+vzy1)*SumAtanXy1z2 - P1.z*(vzy2_m_vyz1+vzy2)*SumAtanXy2z1 + P2.z*(vzy2_m_vyz2+vzy2)*SumAtanXy2z2)
				  + 0.5*V.x*One_d_vyvypvzvz*((2.*vyy1+vzz1)*P1.z*(Log211-Log111) + (2.*vyy1+vzz2)*P2.z*(Log112-Log212) + (2.*vyy2+vzz1)*P1.z*(Log121-Log221) + (2.*vyy2+vzz2)*P2.z*(Log222-Log122))
				  + 0.5*V.z*One_d_vxvxpvyvy*((vx_d_vyvypvzvz_mu_vxvxpvzvzp2vyvy*y1y1-(2.*vyy1+vxx1)*P1.x)*(Log111-Log112) + (-vx_d_vyvypvzvz_mu_vxvxpvzvzp2vyvy*y1y1+(2.*vyy1+vxx2)*P2.x)*(Log211-Log212) + (vx_d_vyvypvzvz_mu_vxvxpvzvzp2vyvy*y2y2-(2.*vyy2+vxx1)*P1.x)*(Log122-Log121) + (-vx_d_vyvypvzvz_mu_vxvxpvzvzp2vyvy*y2y2+(2.*vyy2+vxx2)*P2.x)*(Log222-Log221))
				  + P1.x*P1.z*(Log111-Log121) + P1.x*P2.z*(Log122-Log112) + P2.x*P1.z*(Log221-Log211) + P2.x*P2.z*(Log212-Log222)
				  + (Pi*V.z*V.z*One_d_vyvypvzvz/V.y)*(y1y1*Step(vzy1/V.y-P1.z)*Step(P2.z-vzy1/V.y)*(Sign(vxy1_m_vyx1)-Sign(vxy1_m_vyx2)) + y2y2*Step(vzy2/V.y-P1.z)*Step(P2.z-vzy2/V.y)*(Sign(vxy2_m_vyx2)-Sign(vxy2_m_vyx1)));
			TwoGz = One_d_vxvxpvzvz_d_vyvypvzvz*((vxvx_p_vzvz*(vzy1_m_vyz1-vyz1)*P1.y-vxvx_m_vyvy*vzz1*P1.z)*SumAtanXy1z1 - (vxvx_p_vzvz*(vzy2_m_vyz1-vyz1)*P2.y-vxvx_m_vyvy*vzz1*P1.z)*SumAtanXy2z1 - (vxvx_p_vzvz*(vzy1_m_vyz2-vyz2)*P1.y-vxvx_m_vyvy*vzz2*P2.z)*SumAtanXy1z2 + (vxvx_p_vzvz*(vzy2_m_vyz2-vyz2)*P2.y-vxvx_m_vyvy*vzz2*P2.z)*SumAtanXy2z2)
				  + One_d_vxvxpvzvz*(P1.x*(vzx1_m_vxz1-vxz1)*SumAtanYx1z1 - P2.x*(vzx2_m_vxz1-vxz1)*SumAtanYx2z1 - P1.x*(vzx1_m_vxz2-vxz2)*SumAtanYx1z2 + P2.x*(vzx2_m_vxz2-vxz2)*SumAtanYx2z2)
				  + 0.5*V.y*One_d_vxvxpvzvz*((2.*vzz1+vxx1)*P1.x*(Log121-Log111) + (2.*vzz1+vxx2)*P2.x*(Log211-Log221) + (2.*vzz2+vxx1)*P1.x*(Log112-Log122) + (2.*vzz2+vxx2)*P2.x*(Log222-Log212))
				  + 0.5*V.x*One_d_vyvypvzvz*((vy_d_vxvxpvzvz_mu_vxvxpvyvyp2vzvz*z1z1-(2.*vzz1+vyy1)*P1.y)*(Log111-Log211) + (-vy_d_vxvxpvzvz_mu_vxvxpvyvyp2vzvz*z1z1+(2.*vzz1+vyy2)*P2.y)*(Log121-Log221) + (vy_d_vxvxpvzvz_mu_vxvxpvyvyp2vzvz*z2z2-(2.*vzz2+vyy1)*P1.y)*(Log212-Log112) + (-vy_d_vxvxpvzvz_mu_vxvxpvyvyp2vzvz*z2z2+(2.*vzz2+vyy2)*P2.y)*(Log222-Log122))
				  + P1.x*P1.y*(Log111-Log112) + P2.x*P1.y*(Log212-Log211) + P1.x*P2.y*(Log122-Log121) + P2.x*P2.y*(Log221-Log222)
				  + (Pi*V.x*V.x*One_d_vxvxpvzvz/V.z)*(z1z1*Step(vxz1/V.z-P1.x)*Step(P2.x-vxz1/V.z)*(Sign(vzy2_m_vyz1)-Sign(vzy1_m_vyz1)) + z2z2*Step(vxz2/V.z-P1.x)*Step(P2.x-vxz2/V.z)*(-Sign(vzy2_m_vyz2)+Sign(vzy1_m_vyz2)));
			G.x = 0.5*TwoGx; G.y = 0.5*TwoGy; G.z = 0.5*TwoGz;
		}
	}

FinalDefinitionOfFieldIntegrals:
	if(J_is_Zero)
	{
		const double ConstForM = 1./2./Pi;
		TVector3d ConByM = ConstForM * Magn;
		if(FieldPtr->FieldKey.Ib_)
		{
			TVector3d BufIb(-(F.Str2.x+F.Str1.y)*ConByM.x + F.Str1.z*ConByM.y + F.Str2.z*ConByM.z,
							 F.Str0.y*ConByM.x - (F.Str2.x+F.Str0.z)*ConByM.y + F.Str2.y*ConByM.z,
							 F.Str0.x*ConByM.x + F.Str1.x*ConByM.y - (F.Str1.y+F.Str0.z)*ConByM.z);
			FieldPtr->Ib += BufIb;
		}
		if(FieldPtr->FieldKey.Ih_)
		{
			TVector3d BufIh(F.Str0.z*ConByM.x + F.Str0.y*ConByM.y + F.Str0.x*ConByM.z, 
							F.Str1.z*ConByM.x + F.Str1.y*ConByM.y + F.Str1.x*ConByM.z,
							F.Str2.z*ConByM.x + F.Str2.y*ConByM.y + F.Str2.x*ConByM.z);
			FieldPtr->Ih += BufIh;
		}
	}
	else
	{
		const double ConstForJ = 0.0002;
		TVector3d BufI((J.y*G.z-J.z*G.y)*ConstForJ, (J.z*G.x-J.x*G.z)*ConstForJ, (J.x*G.y-J.y*G.x)*ConstForJ);
		if(FieldPtr->FieldKey.Ib_) FieldPtr->Ib += BufI;
		if(FieldPtr->FieldKey.Ih_) FieldPtr->Ih += BufI;
	}
}

//-------------------------------------------------------------------------

void radTRecCur::B_intUtilSpecCaseZeroVxVy(const TVector3d& P1, const TVector3d& P2, short J_is_Zero, TMatrix3d& F, TVector3d& G)
{
	double z2_m_z1 = P2.z - P1.z;

	const double Pi = 3.141592653589793238;
	double PiMult1, PiMult2, PiMult3;
	PiMult1 = PiMult2 = PiMult3 = 0.;

	double x1x1 = P1.x*P1.x;
	double y1y1 = P1.y*P1.y;
	double x2x2 = P2.x*P2.x;
	double y2y2 = P2.y*P2.y;

	if(J_is_Zero)
	{
		F.Str0.x = F.Str1.x = F.Str2.x = F.Str2.y = F.Str2.z = 0.;
		double SumAtan1 = atan(TransAtans(TransAtans(-P2.x/P1.y, P1.x/P1.y, PiMult1), TransAtans(P2.x/P2.y, -P1.x/P2.y, PiMult2), PiMult3)) 
			            + Pi*(PiMult1+PiMult2+PiMult3);
		F.Str1.y = -z2_m_z1*SumAtan1;
		double SumAtan2 = atan(TransAtans(TransAtans(-P2.y/P1.x, P1.y/P1.x, PiMult1), TransAtans(P2.y/P2.x, -P1.y/P2.x, PiMult2), PiMult3)) 
			            + Pi*(PiMult1+PiMult2+PiMult3);
		F.Str0.z = -z2_m_z1*SumAtan2;
		F.Str0.y = F.Str1.z = 0.5*z2_m_z1*log(((x1x1+y2y2)*(x2x2+y1y1))/((x1x1+y1y1)*(x2x2+y2y2)));
	}
	else
	{
		double SumAtan_x1 = atan(TransAtans(-P1.y/P1.x, P2.y/P1.x, PiMult1)) + Pi*PiMult1;
		double SumAtan_x2 = atan(TransAtans(-P2.y/P2.x, P1.y/P2.x, PiMult1)) + Pi*PiMult1;
		G.x = z2_m_z1*(P1.x*SumAtan_x1 + P2.x*SumAtan_x2 + 0.5*(P1.y*log((x2x2+y1y1)/(x1x1+y1y1)) + P2.y*log((x1x1+y2y2)/(x2x2+y2y2))));
		double SumAtan_y1 = atan(TransAtans(-P1.x/P1.y, P2.x/P1.y, PiMult1)) + Pi*PiMult1;
		double SumAtan_y2 = atan(TransAtans(-P2.x/P2.y, P1.x/P2.y, PiMult1)) + Pi*PiMult1;
		G.y = z2_m_z1*(P1.y*SumAtan_y1 + P2.y*SumAtan_y2 + 0.5*(P1.x*log((x1x1+y2y2)/(x1x1+y1y1)) + P2.x*log((x2x2+y1y1)/(x2x2+y2y2))));
		G.z = 0.;
	}
}

//-------------------------------------------------------------------------

void radTRecCur::FunForOuterIntAtSurfInt(double Arg, TVector3d* VectArray)
{
	const double PrecEnhFact = 1.; // Don't make it >1 : it's dangerous for convergence of outer itegral !
	double* OuterIntPrecArray = SurfIntDataPtr->Field.ShapeIntDataPtr->AbsPrecArray;

	double SmallPositive = 1.E-10;

	if(SurfIntDataPtr->SurfBoundInd==1 || SurfIntDataPtr->SurfBoundInd==2)
	{
		for(int i=0; i<SurfIntDataPtr->IntegrandLen; i++)
			(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[i] = PrecEnhFact*OuterIntPrecArray[i]/Dimensions.y;
		(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[SurfIntDataPtr->IntegrandLen] = CentrPoint.x - 0.5*Dimensions.x + SmallPositive;
		(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[SurfIntDataPtr->IntegrandLen + 1] = CentrPoint.x + 0.5*Dimensions.x;
		SurfIntDataPtr->PointOnSurface.y = Arg;
	}
	else if(SurfIntDataPtr->SurfBoundInd==3 || SurfIntDataPtr->SurfBoundInd==4)
	{
		for(int i=0; i<SurfIntDataPtr->IntegrandLen; i++)
			(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[i] = PrecEnhFact*OuterIntPrecArray[i]/Dimensions.z;
		(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[SurfIntDataPtr->IntegrandLen] = CentrPoint.x - 0.5*Dimensions.x + SmallPositive;
		(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[SurfIntDataPtr->IntegrandLen + 1] = CentrPoint.x + 0.5*Dimensions.x;
		SurfIntDataPtr->PointOnSurface.z = Arg;
	}
	else if(SurfIntDataPtr->SurfBoundInd==5 || SurfIntDataPtr->SurfBoundInd==6)
	{
		for(int i=0; i<SurfIntDataPtr->IntegrandLen; i++)
			(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[i] = PrecEnhFact*OuterIntPrecArray[i]/Dimensions.z;
		(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[SurfIntDataPtr->IntegrandLen] = CentrPoint.y - 0.5*Dimensions.y + SmallPositive;
		(SurfIntDataPtr->InnerAbsPrecAndLimitsArray)[SurfIntDataPtr->IntegrandLen + 1] = CentrPoint.y + 0.5*Dimensions.y;
		SurfIntDataPtr->PointOnSurface.z = Arg;
	}
	FormalOneFoldInteg(this, &radTRecCur::FunForInnerIntAtSurfInt, SurfIntDataPtr->IntegrandLen, 
					   SurfIntDataPtr->InnerAbsPrecAndLimitsArray, 
					   SurfIntDataPtr->InnerElemCompNotFinished, SurfIntDataPtr->InnerIntegVal);

	for(int i=0; i<SurfIntDataPtr->IntegrandLen; i++) VectArray[i] = ((SurfIntDataPtr->InnerIntegVal)[0])[i];
}

//-------------------------------------------------------------------------

void radTRecCur::IntOverSurf(radTField* FieldPtr)
{
	int LenVal = FieldPtr->ShapeIntDataPtr->IntegrandLength;
	int LenValp2 = LenVal+2;
	TVector3d ZeroVect(0.,0.,0.);

	TVector3d* InnerIntegVal[6];
	TVector3d* OuterIntegVal[6];

	std::vector<std::vector<TVector3d>> vInnerIntegValStorage(6);
	std::vector<std::vector<TVector3d>> vOuterIntegValStorage(6);
	std::vector<short> vInnerElemCompNotFinished(LenVal);
	std::vector<short> vOuterElemCompNotFinished(LenVal);
	std::vector<double> vInnerAbsPrecAndLimitsArray(LenValp2);
	std::vector<double> vOuterAbsPrecAndLimitsArray(LenValp2);
	std::vector<TVector3d> vLocalVectArray(LenVal);

	int j;
	for(j=0; j<6; j++)
	{
		vInnerIntegValStorage[j].resize(LenVal);
		vOuterIntegValStorage[j].resize(LenVal);
		InnerIntegVal[j] = vInnerIntegValStorage[j].data();
		OuterIntegVal[j] = vOuterIntegValStorage[j].data();
	}
	short* InnerElemCompNotFinished = vInnerElemCompNotFinished.data();
	short* OuterElemCompNotFinished = vOuterElemCompNotFinished.data();
	double* InnerAbsPrecAndLimitsArray = vInnerAbsPrecAndLimitsArray.data();
	double* OuterAbsPrecAndLimitsArray = vOuterAbsPrecAndLimitsArray.data();
	TVector3d* LocalVectArray = vLocalVectArray.data();
	SurfIntDataPtr = new radTParallelepSurfIntData();

	SurfIntDataPtr->IntegrandLen = LenVal;
	SurfIntDataPtr->IntegrandFunPtr = FieldPtr->ShapeIntDataPtr->IntegrandFunPtr;
	SurfIntDataPtr->InnerAbsPrecAndLimitsArray = InnerAbsPrecAndLimitsArray;
	SurfIntDataPtr->InnerElemCompNotFinished = InnerElemCompNotFinished;
	SurfIntDataPtr->InnerIntegVal = InnerIntegVal;
	
	SurfIntDataPtr->Field = *FieldPtr;
	radTStructForShapeInt LocShapeIntData = *(FieldPtr->ShapeIntDataPtr);

	TVector3d* InputFieldPtrVectArrayPtr = FieldPtr->ShapeIntDataPtr->VectArray;

	LocShapeIntData.VectArray = LocalVectArray;
	SurfIntDataPtr->Field.ShapeIntDataPtr = &LocShapeIntData;

	int i;
	for(i=0; i<LenVal; i++)
	{
		OuterAbsPrecAndLimitsArray[i] = (FieldPtr->ShapeIntDataPtr->AbsPrecArray)[i];
	}
	TVector3d HalfDim = 0.5*Dimensions;

	TVector3d* OutVectArray = InputFieldPtrVectArrayPtr;

	double SmallPositive = 1.E-10;

//For lower and upper bounds
	OuterAbsPrecAndLimitsArray[LenVal] = CentrPoint.y - HalfDim.y + SmallPositive;
	OuterAbsPrecAndLimitsArray[LenVal+1] = CentrPoint.y + HalfDim.y;
//Integration over lower bound
	SurfIntDataPtr->SurfBoundInd = 1;
	SurfIntDataPtr->PointOnSurface.z = CentrPoint.z - HalfDim.z + SmallPositive;
	SurfIntDataPtr->Field.ShapeIntDataPtr->Normal = TVector3d(0.,0.,-1.);
	FormalOneFoldInteg(this, &radTRecCur::FunForOuterIntAtSurfInt, LenVal, OuterAbsPrecAndLimitsArray, OuterElemCompNotFinished, OuterIntegVal);
	for(i=0; i<LenVal; i++) OutVectArray[i] += (OuterIntegVal[0])[i];
//Integration over upper bound
	SurfIntDataPtr->SurfBoundInd = 2;
	SurfIntDataPtr->PointOnSurface.z = CentrPoint.z + HalfDim.z + SmallPositive;
	SurfIntDataPtr->Field.ShapeIntDataPtr->Normal = TVector3d(0.,0.,1.);
	FormalOneFoldInteg(this, &radTRecCur::FunForOuterIntAtSurfInt, LenVal, OuterAbsPrecAndLimitsArray, OuterElemCompNotFinished, OuterIntegVal);
	for(i=0; i<LenVal; i++) OutVectArray[i] += (OuterIntegVal[0])[i];

//For left, right, back and front bounds
	OuterAbsPrecAndLimitsArray[LenVal] = CentrPoint.z - HalfDim.z + SmallPositive;
	OuterAbsPrecAndLimitsArray[LenVal+1] = CentrPoint.z + HalfDim.z;
//Integration over left bound
	SurfIntDataPtr->SurfBoundInd = 3;
	SurfIntDataPtr->PointOnSurface.y = CentrPoint.y - HalfDim.y + SmallPositive;
	SurfIntDataPtr->Field.ShapeIntDataPtr->Normal = TVector3d(0.,-1.,0.);
	FormalOneFoldInteg(this, &radTRecCur::FunForOuterIntAtSurfInt, LenVal, OuterAbsPrecAndLimitsArray, OuterElemCompNotFinished, OuterIntegVal);
	for(i=0; i<LenVal; i++) OutVectArray[i] += (OuterIntegVal[0])[i];
//Integration over right bound
	SurfIntDataPtr->SurfBoundInd = 4;
	SurfIntDataPtr->PointOnSurface.y = CentrPoint.y + HalfDim.y + SmallPositive;
	SurfIntDataPtr->Field.ShapeIntDataPtr->Normal = TVector3d(0.,1.,0.);
	FormalOneFoldInteg(this, &radTRecCur::FunForOuterIntAtSurfInt, LenVal, OuterAbsPrecAndLimitsArray, OuterElemCompNotFinished, OuterIntegVal);
	for(i=0; i<LenVal; i++) OutVectArray[i] += (OuterIntegVal[0])[i];
//Integration over back bound
	SurfIntDataPtr->SurfBoundInd = 5;
	SurfIntDataPtr->PointOnSurface.x = CentrPoint.x - HalfDim.x + SmallPositive;
	SurfIntDataPtr->Field.ShapeIntDataPtr->Normal = TVector3d(-1.,0.,0.);
	FormalOneFoldInteg(this, &radTRecCur::FunForOuterIntAtSurfInt, LenVal, OuterAbsPrecAndLimitsArray, OuterElemCompNotFinished, OuterIntegVal);
	for(i=0; i<LenVal; i++) OutVectArray[i] += (OuterIntegVal[0])[i];
//Integration over right bound
	SurfIntDataPtr->SurfBoundInd = 6;
	SurfIntDataPtr->PointOnSurface.x = CentrPoint.x + HalfDim.x + SmallPositive;
	SurfIntDataPtr->Field.ShapeIntDataPtr->Normal = TVector3d(1.,0.,0.);
	FormalOneFoldInteg(this, &radTRecCur::FunForOuterIntAtSurfInt, LenVal, OuterAbsPrecAndLimitsArray, OuterElemCompNotFinished, OuterIntegVal);
	for(i=0; i<LenVal; i++) OutVectArray[i] += (OuterIntegVal[0])[i];

// Automatic cleanup via RAII
	delete SurfIntDataPtr;
}

//-------------------------------------------------------------------------


//-------------------------------------------------------------------------

// radTRecCur::SubdivideItself REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTRecCur::Dump / DumpPureObjInfo REMOVED (Phase B2b, 2026-04-15)

//-------------------------------------------------------------------------

// DumpBin / DumpBin_RecMag / DumpBinParse_RecMag REMOVED (Phase B2c, 2026-04-15)

//-------------------------------------------------------------------------

// radTRecCur::ConvertToPolyhedron body REMOVED 2026-06-28: a current block is never converted
// (the magnet path now builds the MMMM polyhedron directly in SetRecMag). The override is now an
// inline { return 1; } in rad_rectangular_block.h (1 = "handled", so radTGroup conversion succeeds).

//-------------------------------------------------------------------------

// radTRecCur::SubdivideItselfByOneSetOfParPlanes REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTRecCur::SubdivideItselfByPlanesParToFace REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTRecCur::CutItself REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTRecCur::FindLowestAndUppestVertices REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTRecCur::CheckVertexPtsPositionsWithRespectToPlane REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

void radTRecCur::Push_backCenterPointAndField(radTFieldKey* pFieldKey, radTVectPairOfVect3d* pVectPairOfVect3d, radTrans* pBaseTrans, radTg3d* g3dSrcPtr, radTApplication* pAppl)
{// Attention: this assumes no more than one transformation with mult. no more than 1 !!!
	TVector3d CP = CentrPoint;
	radTrans* pTrans = (g3dListOfTransform.empty())? 0 : (radTrans*)((*(g3dListOfTransform.begin())).Handler_g.rep);

	radTrans TotTrans;
	if(pTrans != 0)
	{
		if(pBaseTrans != 0) 
		{
			TrProduct(pBaseTrans, pTrans, TotTrans);
			pTrans = &TotTrans;
		}
	}
	else
	{
		if(pBaseTrans != 0) pTrans = pBaseTrans;
	}
	
	if(pTrans != 0) CP = pTrans->TrPoint(CP);
	radTPairOfVect3d Pair(CP);

	if(pFieldKey->M_)
	{
		if(J_IsNotZero) return;
		Pair.V2 = (pTrans == 0)? Magn : pTrans->TrVectField(Magn);
	}
	else if(pFieldKey->J_)
	{
		if(!J_IsNotZero) return;
		Pair.V2 = (pTrans == 0)? J : pTrans->TrVectField(J);
	}
	else
	{
		radTCompCriterium CompCriterium;
		TVector3d ZeroVect(0.,0.,0.);
		radTField Field(*pFieldKey, CompCriterium, CP, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
		g3dSrcPtr->B_genComp(&Field);
		Pair.V2 = (pFieldKey->B_)? Field.B : ((pFieldKey->H_)? Field.H : ((pFieldKey->A_)? Field.A : ZeroVect));
	}
	pVectPairOfVect3d->push_back(Pair);
}

//-------------------------------------------------------------------------

void radTRecCur::VerticesInLocFrame(radTVectorOfVector3d& OutVect, bool EnsureUnique)
{
	TVector3d HalfDim = 0.5*Dimensions;
	double xMin = CentrPoint.x - HalfDim.x, xMax = CentrPoint.x + HalfDim.x;
	double yMin = CentrPoint.y - HalfDim.y, yMax = CentrPoint.y + HalfDim.y;
	double zMin = CentrPoint.z - HalfDim.z, zMax = CentrPoint.z + HalfDim.z;

	TVector3d P1(xMin, yMin, zMin); OutVect.push_back(P1);
	TVector3d P2(xMax, yMin, zMin); OutVect.push_back(P2);
	TVector3d P3(xMin, yMax, zMin); OutVect.push_back(P3);
	TVector3d P4(xMax, yMax, zMin); OutVect.push_back(P4);

	TVector3d P5(xMin, yMin, zMax); OutVect.push_back(P5);
	TVector3d P6(xMax, yMin, zMax); OutVect.push_back(P6);
	TVector3d P7(xMin, yMax, zMax); OutVect.push_back(P7);
	TVector3d P8(xMax, yMax, zMax); OutVect.push_back(P8);
}

//-------------------------------------------------------------------------
