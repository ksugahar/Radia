/*-------------------------------------------------------------------------
*
* File name:      radpoly_analytical.cpp
*
* Project:        RADIA
*
* Description:    Analytical field formula from polygon magnetic charges
*
* Author(s):      Radia Development Team
*
* First release:  2025-11-08
*
*-------------------------------------------------------------------------*/

#include "radpoly_analytical.h"
#include <cmath>
#include <algorithm>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

//-------------------------------------------------------------------------
// Helper functions
//-------------------------------------------------------------------------

inline double SQR3(double x, double y, double z)
{
	return std::sqrt(x*x + y*y + z*z);
}

inline double SQ2(double x, double y)
{
	return std::sqrt(x*x + y*y);
}

//-------------------------------------------------------------------------
/**
 * Compute field from polygon magnetic charge using analytical formula
 *
 * @param AA Local coordinate system X-axis (unit vector)
 * @param BB Local coordinate system Y-axis (unit vector)
 * @param CC Local coordinate system Z-axis (normal vector)
 * @param YY Reference point on polygon (3D)
 * @param XY Polygon vertices in local 2D coordinates (KAdo points)
 * @param XX Observation points in global 3D coordinates (MXX points)
 * @param FGH Output: magnetic field at each observation point (3 x MXX)
 * @param W Magnetic charge density (weight)
 * @param MXX Number of observation points
 * @param NII Element index (for error reporting)
 * @param KAdo Number of polygon vertices (3 or 4)
 *
 * @note This function uses the analytical formula:
 *       H = (1/4π) * ∮ σ * dΩ
 *       where dΩ is the solid angle subtended by each edge
 */
void RadAnalyticalFieldFromPolygonCharge(
	const TVector3d& AA,
	const TVector3d& BB,
	const TVector3d& CC,
	const TVector3d& YY,
	const std::vector<TVector2d>& XY,  // Polygon vertices in 2D local coords
	const std::vector<TVector3d>& XX,  // Observation points in 3D
	std::vector<TVector3d>& FGH,       // Output: field at each point
	double W,                           // Magnetic charge density
	int NII,                            // Element index
	int KAdo)                           // Number of vertices (3 or 4)
{
	const double ONE = 1.0;
	const double ZER = 0.0;
	const double EPS = 1.0e-20;
	const double BIG = 1.0e20;

	int MXX = static_cast<int>(XX.size());
	if(MXX == 0) return;

	// Ensure output array is sized correctly
	if(FGH.size() != static_cast<size_t>(MXX)) {
		FGH.resize(MXX, TVector3d(0, 0, 0));
	}

	// Compute edge properties
	std::vector<double> DS(4);   // Edge lengths
	std::vector<double> AM(4);   // Slopes (dy/dx)
	std::vector<double> SM(4);   // sqrt(1 + slope^2)
	std::vector<double> XD(4);   // Edge direction X
	std::vector<double> YD(4);   // Edge direction Y

	double EPSG = 0.0;
	double ZONE = (KAdo == 3) ? ZER : ONE;

	// Compute edge parameters
	for(int J = 0; J < KAdo; J++) {
		int L = (J + 1) % KAdo;  // Next vertex

		double XS1 = XY[L].x - XY[J].x;
		double XS2 = XY[L].y - XY[J].y;

		if(std::abs(XS1) < EPS) {
			// Vertical edge - handle specially
			std::cerr << "[Warning] RadAnalyticalFieldFromPolygonCharge: vertical edge in element "
			          << NII << std::endl;
			XS1 = EPS;
		}

		DS[J] = SQ2(XS2, XS1);
		AM[J] = XS2 / XS1;
		SM[J] = std::sqrt(AM[J]*AM[J] + ONE);
		XD[J] = -XS1 / DS[J];
		YD[J] =  XS2 / DS[J];

		EPSG = std::max(EPSG, DS[J]);
	}

	EPSG = EPSG * 1.0e-12;  // Tolerance for z=0 check

	// For triangle, set 4th edge to dummy values
	if(KAdo == 3) {
		DS[3] = ONE;
		AM[3] = ZER;
		SM[3] = BIG;
		XD[3] = ZER;
		YD[3] = ZER;
	}

	// Main loop over observation points (can be parallelized)
	#pragma omp parallel for if(MXX > 100)
	for(int I = 0; I < MXX; I++) {
		// Transform observation point to local coordinates
		TVector3d DD = XX[I] - YY;

		double EE1 = DD.x*AA.x + DD.y*AA.y + DD.z*AA.z;  // X in local frame
		double EE2 = DD.x*BB.x + DD.y*BB.y + DD.z*BB.z;  // Y in local frame
		double EE3 = DD.x*CC.x + DD.y*CC.y + DD.z*CC.z;  // Z in local frame (height)

		// Compute distances from observation point to vertices
		std::vector<double> X(4), Y(4), H(4), E(4), R(4);

		for(int J = 0; J < 4; J++) {
			X[J] = EE1 - XY[J].x;
			Y[J] = EE2 - XY[J].y;
			H[J] = Y[J] * X[J];
		}

		double Z = EE3;
		double Z2 = Z * Z;

		for(int J = 0; J < 4; J++) {
			E[J] = Z2 + X[J]*X[J];
			R[J] = SQR3(X[J], Y[J], Z);
		}

		// Compute edge contributions
		std::vector<double> RM(4), RP(4), RR(4), AL(4);

		for(int J = 0; J < 4; J++) {
			int JP1 = (J + 1) % 4;

			RM[J] = R[J] + R[JP1] - DS[J];
			RP[J] = R[J] + R[JP1] + DS[J];
			RR[J] = std::max(RM[J] / RP[J], EPS);
			AL[J] = std::log(RR[J]);
		}

		// Compute field components in local frame
		double HH1 = W * (-YD[0]*AL[0] - YD[1]*AL[1] - YD[2]*AL[2] - YD[3]*AL[3]);
		double HH2 = W * (-XD[0]*AL[0] - XD[1]*AL[1] - XD[2]*AL[2] - XD[3]*AL[3]);
		double HH3 = ZER;

		// Z-component (solid angle contribution)
		if(std::abs(Z) > EPSG) {
			std::vector<double> ZR(4), AT(4), BT(4);

			for(int J = 0; J < 4; J++) {
				ZR[J] = Z * R[J];
			}

			// Compute arctan terms
			for(int J = 0; J < 4; J++) {
				int JP1 = (J + 1) % 4;

				AT[J] = (AM[J]*E[J] - H[J]) / ZR[J];
				BT[J] = (AM[J]*E[JP1] - H[JP1]) / ZR[JP1];
			}

			// Special handling for triangle (4th term)
			AT[3] = AT[3] * ZONE;
			BT[3] = BT[3] * ZONE;

			HH3 = W * (-std::atan(AT[0]) - std::atan(AT[1]) - std::atan(AT[2]) - std::atan(AT[3])
			          + std::atan(BT[0]) + std::atan(BT[1]) + std::atan(BT[2]) + std::atan(BT[3]));
		}

		// Transform field back to global coordinates
		FGH[I].x += HH1*AA.x + HH2*BB.x + HH3*CC.x;
		FGH[I].y += HH1*AA.y + HH2*BB.y + HH3*CC.y;
		FGH[I].z += HH1*AA.z + HH2*BB.z + HH3*CC.z;
	}
}

//-------------------------------------------------------------------------
/**
 * Compute field from triangular magnetic charge
 *
 * Convenience wrapper for 3-vertex polygons
 */
void RadAnalyticalFieldFromTriangleCharge(
	const TVector3d& AA,
	const TVector3d& BB,
	const TVector3d& CC,
	const TVector3d& YY,
	const TVector2d& V1,
	const TVector2d& V2,
	const TVector2d& V3,
	const std::vector<TVector3d>& XX,
	std::vector<TVector3d>& FGH,
	double W,
	int NII)
{
	std::vector<TVector2d> XY = {V1, V2, V3};
	RadAnalyticalFieldFromPolygonCharge(AA, BB, CC, YY, XY, XX, FGH, W, NII, 3);
}

//-------------------------------------------------------------------------
/**
 * Compute field from quadrilateral magnetic charge
 *
 * Convenience wrapper for 4-vertex polygons
 */
void RadAnalyticalFieldFromQuadCharge(
	const TVector3d& AA,
	const TVector3d& BB,
	const TVector3d& CC,
	const TVector3d& YY,
	const TVector2d& V1,
	const TVector2d& V2,
	const TVector2d& V3,
	const TVector2d& V4,
	const std::vector<TVector3d>& XX,
	std::vector<TVector3d>& FGH,
	double W,
	int NII)
{
	std::vector<TVector2d> XY = {V1, V2, V3, V4};
	RadAnalyticalFieldFromPolygonCharge(AA, BB, CC, YY, XY, XX, FGH, W, NII, 4);
}
