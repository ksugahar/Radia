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

#include "rad_poly_analytical.h"
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
 *       H = (1/4pi) * integral sigma * dOmega
 *       where dOmega is the solid angle subtended by each edge
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
			// Vertical edge - handle specially by replacing with small non-zero value
			// This is a numerical stability measure for near-vertical edges in tetrahedral meshes
			// Warning suppressed to avoid console spam (vertical edges are common in Netgen meshes)
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
		// Initialize to zero; only use KAdo vertices from XY (triangles have 3, quads have 4)
		std::vector<double> X(4, 0.0), Y(4, 0.0), H(4, 0.0), E(4, 0.0), R(4, 0.0);

		for(int J = 0; J < KAdo; J++) {
			X[J] = EE1 - XY[J].x;
			Y[J] = EE2 - XY[J].y;
			H[J] = Y[J] * X[J];
		}
		// For triangles (KAdo=3), set element [3] = element [0] for edge connectivity
		// This is needed because the edge loop (J=0,1,2,3) uses JP1=(J+1)%4
		// Edge 3 connects vertex[3] to vertex[0], but for triangles there is no vertex[3]
		if(KAdo == 3) {
			X[3] = X[0];
			Y[3] = Y[0];
			H[3] = H[0];
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

//-------------------------------------------------------------------------
/**
 * Compute H field from triangular face using GLOBAL 3D coordinates directly
 *
 * This function computes the magnetic field from a triangular face with
 * surface charge density sigma = M dot n, using global coordinates only.
 * No local coordinate transformation is needed.
 *
 * @param V0, V1, V2  Triangle vertices in GLOBAL 3D coordinates
 * @param M           Magnetization vector (in global coordinates)
 * @param obsPoint    Observation point in GLOBAL 3D coordinates
 * @return            H field at observation point (in global coordinates)
 *
 * Formula: H = (sigma/4pi) * Omega  where Omega is the solid angle
 *
 * Reference: Lindholm 1984, "Three-Dimensional Magnetostatic Fields from
 *            Point-Matched Integral Equations with Linearly Varying
 *            Scalar Sources on a Triangular Grid"
 */
TVector3d RadFieldFromTriangleFaceGlobal(
	const TVector3d& V0,
	const TVector3d& V1,
	const TVector3d& V2,
	const TVector3d& M,
	const TVector3d& obsPoint)
{
	const double PI4 = 4.0 * 3.14159265358979323846;
	const double EPS = 1.0e-15;

	// Compute face normal from cross product (V1-V0) x (V2-V0)
	TVector3d e1 = V1 - V0;
	TVector3d e2 = V2 - V0;
	TVector3d normal;
	normal.x = e1.y * e2.z - e1.z * e2.y;
	normal.y = e1.z * e2.x - e1.x * e2.z;
	normal.z = e1.x * e2.y - e1.y * e2.x;

	double normalLen = std::sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
	if(normalLen < EPS) {
		return TVector3d(0., 0., 0.);  // Degenerate triangle
	}

	// Normalize
	normal.x /= normalLen;
	normal.y /= normalLen;
	normal.z /= normalLen;

	// Surface charge density: sigma = M dot n
	double sigma = M.x * normal.x + M.y * normal.y + M.z * normal.z;
	if(std::abs(sigma) < EPS) {
		return TVector3d(0., 0., 0.);  // No surface charge
	}

	// Vectors from each vertex to observation point (IMPORTANT: direction matters for sign!)
	TVector3d r0 = obsPoint - V0;
	TVector3d r1 = obsPoint - V1;
	TVector3d r2 = obsPoint - V2;

	double R0 = std::sqrt(r0.x*r0.x + r0.y*r0.y + r0.z*r0.z);
	double R1 = std::sqrt(r1.x*r1.x + r1.y*r1.y + r1.z*r1.z);
	double R2 = std::sqrt(r2.x*r2.x + r2.y*r2.y + r2.z*r2.z);

	// Check for observation point at vertex
	if(R0 < EPS || R1 < EPS || R2 < EPS) {
		// Point is at or very near a vertex - use regularized value
		return TVector3d(0., 0., 0.);
	}

	// ====================================================================
	// Compute solid angle using Van Oosterom & Strackee formula
	// Omega = 2 * atan2(numerator, denominator)
	// numerator = r0 . (r1 x r2)
	// denominator = R0*R1*R2 + (r0.r1)*R2 + (r0.r2)*R1 + (r1.r2)*R0
	// ====================================================================

	// Cross product r1 x r2
	TVector3d r1xr2;
	r1xr2.x = r1.y * r2.z - r1.z * r2.y;
	r1xr2.y = r1.z * r2.x - r1.x * r2.z;
	r1xr2.z = r1.x * r2.y - r1.y * r2.x;

	// Triple product: r0 . (r1 x r2)
	double tripleProduct = r0.x * r1xr2.x + r0.y * r1xr2.y + r0.z * r1xr2.z;

	// Dot products
	double r0_r1 = r0.x * r1.x + r0.y * r1.y + r0.z * r1.z;
	double r0_r2 = r0.x * r2.x + r0.y * r2.y + r0.z * r2.z;
	double r1_r2 = r1.x * r2.x + r1.y * r2.y + r1.z * r2.z;

	double denom = R0*R1*R2 + r0_r1*R2 + r0_r2*R1 + r1_r2*R0;

	double Omega = 2.0 * std::atan2(tripleProduct, denom);

	// ====================================================================
	// H field contribution from solid angle
	// For a uniformly charged triangle, H = -(sigma/4pi) * grad(Omega)
	// For observation point NOT on the triangle plane:
	//   H = (sigma/4pi) * Omega * n  (along normal direction)
	// This is the contribution to H from surface charge.
	// ====================================================================

	// However, the correct formula for H from surface charge is more complex.
	// Let's use the edge contribution formula (Lindholm approach):
	//
	// For each edge of the triangle, compute:
	//   contribution = (sigma/4pi) * ln((R_start + R_end + L)/(R_start + R_end - L)) * (edge x rho)
	// where L is edge length and rho is perpendicular from obs point to edge

	// Simpler approach: use the gradient of scalar potential
	// phi = -(sigma/4pi) * Omega
	// H = -grad(phi) = (sigma/4pi) * grad(Omega)
	//
	// For a planar triangle, grad(Omega) has components both in-plane and out-of-plane.
	// The closed-form expression is:
	//   grad(Omega) = sum over edges of: (1/|r_cross|^2) * (r_cross) * angle_contribution
	//
	// A simpler result from the paper:
	//   H_n = (sigma/4pi) * Omega  (normal component)
	// But we need the full 3D vector.

	// Use the edge sum formula from ELF:
	TVector3d H(0., 0., 0.);

	// Edge 0: V0 -> V1
	{
		TVector3d edge = V1 - V0;
		double L = std::sqrt(edge.x*edge.x + edge.y*edge.y + edge.z*edge.z);
		if(L > EPS) {
			double Rsum = R0 + R1;
			if(std::abs(Rsum - L) > EPS && std::abs(Rsum + L) > EPS) {
				double logArg = (Rsum + L) / (Rsum - L);
				if(logArg > EPS) {
					double logVal = std::log(logArg);
					// Cross product: edge x r0 gives perpendicular direction
					TVector3d edgeCrossR0;
					edgeCrossR0.x = edge.y * r0.z - edge.z * r0.y;
					edgeCrossR0.y = edge.z * r0.x - edge.x * r0.z;
					edgeCrossR0.z = edge.x * r0.y - edge.y * r0.x;
					double factor = logVal / (L * L);
					H.x += factor * edgeCrossR0.x;
					H.y += factor * edgeCrossR0.y;
					H.z += factor * edgeCrossR0.z;
				}
			}
		}
	}

	// Edge 1: V1 -> V2
	{
		TVector3d edge = V2 - V1;
		double L = std::sqrt(edge.x*edge.x + edge.y*edge.y + edge.z*edge.z);
		if(L > EPS) {
			double Rsum = R1 + R2;
			if(std::abs(Rsum - L) > EPS && std::abs(Rsum + L) > EPS) {
				double logArg = (Rsum + L) / (Rsum - L);
				if(logArg > EPS) {
					double logVal = std::log(logArg);
					TVector3d edgeCrossR1;
					edgeCrossR1.x = edge.y * r1.z - edge.z * r1.y;
					edgeCrossR1.y = edge.z * r1.x - edge.x * r1.z;
					edgeCrossR1.z = edge.x * r1.y - edge.y * r1.x;
					double factor = logVal / (L * L);
					H.x += factor * edgeCrossR1.x;
					H.y += factor * edgeCrossR1.y;
					H.z += factor * edgeCrossR1.z;
				}
			}
		}
	}

	// Edge 2: V2 -> V0
	{
		TVector3d edge = V0 - V2;
		double L = std::sqrt(edge.x*edge.x + edge.y*edge.y + edge.z*edge.z);
		if(L > EPS) {
			double Rsum = R2 + R0;
			if(std::abs(Rsum - L) > EPS && std::abs(Rsum + L) > EPS) {
				double logArg = (Rsum + L) / (Rsum - L);
				if(logArg > EPS) {
					double logVal = std::log(logArg);
					TVector3d edgeCrossR2;
					edgeCrossR2.x = edge.y * r2.z - edge.z * r2.y;
					edgeCrossR2.y = edge.z * r2.x - edge.x * r2.z;
					edgeCrossR2.z = edge.x * r2.y - edge.y * r2.x;
					double factor = logVal / (L * L);
					H.x += factor * edgeCrossR2.x;
					H.y += factor * edgeCrossR2.y;
					H.z += factor * edgeCrossR2.z;
				}
			}
		}
	}

	// Add normal component from solid angle
	double HnormalMag = Omega;
	H.x += HnormalMag * normal.x;
	H.y += HnormalMag * normal.y;
	H.z += HnormalMag * normal.z;

	// Scale by sigma/(4*pi)
	double scale = sigma / PI4;
	H.x *= scale;
	H.y *= scale;
	H.z *= scale;

	return H;
}
