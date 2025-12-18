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
 *
 * FIXED 2025-12-06: Rewritten to use RadAnalyticalFieldFromPolygonCharge
 * which is the proven analytical formula used by the standard polygon method.
 *
 * FIXED 2025-12-06: Added outward normal check. If the computed normal points
 * inward (toward element centroid), the surface charge sign is negated to
 * ensure correct field direction.
 *
 * FIXED 2025-12-07: Basis vector construction uses triple cross product method
 * to ensure consistent orthonormal basis across all triangles.
 *
 * @param V0, V1, V2    Triangle vertices in GLOBAL 3D coordinates
 * @param M             Magnetization vector (in global coordinates)
 * @param obsPoint      Observation point in GLOBAL 3D coordinates
 * @param elemCentroid  Element (tetrahedron) centroid for outward normal check
 * @return              H field at observation point (in global coordinates)
 */
TVector3d RadFieldFromTriangleFaceGlobal(
	const TVector3d& V0,
	const TVector3d& V1,
	const TVector3d& V2,
	const TVector3d& M,
	const TVector3d& obsPoint,
	const TVector3d& elemCentroid)
{
	// =========================================================================
	// Rewritten to use RadAnalyticalFieldFromPolygonCharge for correctness
	// This ensures consistency with the working standard polygon method.
	// =========================================================================
	const double PI = 3.14159265358979323846;
	const double ConstForH = 1.0 / (4.0 * PI);
	const double EPS = 1.0e-15;

	// Compute face edges
	TVector3d e1 = V1 - V0;
	TVector3d e2 = V2 - V0;

	// Compute face normal from cross product (V1-V0) x (V2-V0)
	TVector3d normal;
	normal.x = e1.y * e2.z - e1.z * e2.y;
	normal.y = e1.z * e2.x - e1.x * e2.z;
	normal.z = e1.x * e2.y - e1.y * e2.x;

	double normalLen = std::sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
	if(normalLen < EPS) {
		return TVector3d(0., 0., 0.);  // Degenerate triangle
	}

	// Normalize to get CC (local Z-axis = face normal)
	TVector3d CC;
	CC.x = normal.x / normalLen;
	CC.y = normal.y / normalLen;
	CC.z = normal.z / normalLen;

	// Surface charge density: sigma = M dot n
	// Note: The sign of sigma will be corrected below based on normal direction
	double sigma = M.x * CC.x + M.y * CC.y + M.z * CC.z;

	// =========================================================================
	// Outward normal check:
	// Compute face center and check if normal points outward from element centroid
	// If normal . (faceCenter - elemCentroid) < 0, the normal points inward
	// In that case, we negate the surface charge (NOT the normal itself).
	// The coordinate system stays as-is, but the charge sign is corrected.
	// =========================================================================
	TVector3d faceCenter;
	faceCenter.x = (V0.x + V1.x + V2.x) / 3.0;
	faceCenter.y = (V0.y + V1.y + V2.y) / 3.0;
	faceCenter.z = (V0.z + V1.z + V2.z) / 3.0;

	TVector3d outwardVec;
	outwardVec.x = faceCenter.x - elemCentroid.x;
	outwardVec.y = faceCenter.y - elemCentroid.y;
	outwardVec.z = faceCenter.z - elemCentroid.z;

	double dotProduct = CC.x * outwardVec.x + CC.y * outwardVec.y + CC.z * outwardVec.z;

	// If normal points inward, negate the surface charge
	// (do not flip the coordinate system)
	if(dotProduct < 0.0) {
		sigma = -sigma;
	}
	if(std::abs(sigma) < EPS) {
		return TVector3d(0., 0., 0.);  // No surface charge
	}

	// =========================================================================
	// Build local coordinate system using triple cross product method
	//
	// Basis construction:
	//   vert_rel(2) = v2 - v1 (edge from v1 to v2)
	//   vert_rel(3) = v3 - v1 (edge from v1 to v3)
	//   basis_a = vert_rel(3)
	//   basis_b = vert_rel(3) - vert_rel(2) * 0.5
	//   Then triple cross product to orthonormalize:
	//     basis_c = basis_a x basis_b (normalized)
	//     basis_a = basis_b x basis_c (normalized)
	//     basis_b = basis_c x basis_a (normalized)
	// =========================================================================

	// Formulation: basis_a = e2, basis_b = e2 - e1*0.5
	TVector3d basis_a = e2;
	TVector3d basis_b;
	basis_b.x = e2.x - e1.x * 0.5;
	basis_b.y = e2.y - e1.y * 0.5;
	basis_b.z = e2.z - e1.z * 0.5;

	// Triple cross product #1: basis_c = basis_a x basis_b (normalized)
	TVector3d basis_c;
	basis_c.x = basis_a.y * basis_b.z - basis_a.z * basis_b.y;
	basis_c.y = basis_a.z * basis_b.x - basis_a.x * basis_b.z;
	basis_c.z = basis_a.x * basis_b.y - basis_a.y * basis_b.x;
	double cLen = std::sqrt(basis_c.x*basis_c.x + basis_c.y*basis_c.y + basis_c.z*basis_c.z);
	if(cLen < EPS) {
		return TVector3d(0., 0., 0.);  // Degenerate triangle
	}
	basis_c.x /= cLen; basis_c.y /= cLen; basis_c.z /= cLen;

	// Triple cross product #2: basis_a = basis_b x basis_c (normalized)
	TVector3d new_basis_a;
	new_basis_a.x = basis_b.y * basis_c.z - basis_b.z * basis_c.y;
	new_basis_a.y = basis_b.z * basis_c.x - basis_b.x * basis_c.z;
	new_basis_a.z = basis_b.x * basis_c.y - basis_b.y * basis_c.x;
	double aLen = std::sqrt(new_basis_a.x*new_basis_a.x + new_basis_a.y*new_basis_a.y + new_basis_a.z*new_basis_a.z);
	if(aLen < EPS) {
		return TVector3d(0., 0., 0.);  // Degenerate triangle
	}
	new_basis_a.x /= aLen; new_basis_a.y /= aLen; new_basis_a.z /= aLen;
	basis_a = new_basis_a;

	// Triple cross product #3: basis_b = basis_c x basis_a (normalized)
	TVector3d new_basis_b;
	new_basis_b.x = basis_c.y * basis_a.z - basis_c.z * basis_a.y;
	new_basis_b.y = basis_c.z * basis_a.x - basis_c.x * basis_a.z;
	new_basis_b.z = basis_c.x * basis_a.y - basis_c.y * basis_a.x;
	double bLen = std::sqrt(new_basis_b.x*new_basis_b.x + new_basis_b.y*new_basis_b.y + new_basis_b.z*new_basis_b.z);
	if(bLen < EPS) {
		return TVector3d(0., 0., 0.);  // Degenerate triangle
	}
	new_basis_b.x /= bLen; new_basis_b.y /= bLen; new_basis_b.z /= bLen;
	basis_b = new_basis_b;

	// AA = local X-axis, BB = local Y-axis
	TVector3d AA = basis_a;
	TVector3d BB = basis_b;

	// Reference point YY = V0 (first vertex)
	TVector3d YY = V0;

	// Convert triangle vertices to 2D local coordinates
	TVector2d local_V0(0.0, 0.0);

	TVector3d d1 = V1 - V0;
	TVector2d local_V1(
		d1.x*AA.x + d1.y*AA.y + d1.z*AA.z,
		d1.x*BB.x + d1.y*BB.y + d1.z*BB.z
	);

	TVector3d d2 = V2 - V0;
	TVector2d local_V2(
		d2.x*AA.x + d2.y*AA.y + d2.z*AA.z,
		d2.x*BB.x + d2.y*BB.y + d2.z*BB.z
	);

	// Prepare observation point and output vectors
	std::vector<TVector3d> obs_points(1, obsPoint);
	std::vector<TVector3d> field_result(1, TVector3d(0., 0., 0.));

	// Weight: sigma / (4*pi)
	double W = ConstForH * sigma;

	// Create 2D vertex array for the triangle
	std::vector<TVector2d> XY = {local_V0, local_V1, local_V2};

	// Call the proven analytical formula
	RadAnalyticalFieldFromPolygonCharge(
		AA, BB, basis_c, YY,
		XY,
		obs_points,
		field_result,
		W,
		1,   // Element index (for error reporting)
		3    // Triangle has 3 vertices
	);

	return field_result[0];
}
