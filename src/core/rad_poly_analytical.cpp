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
#include "rad_constants.h"  // Unified mathematical/physical constants

#include "rad_parallel.h"

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

	// Main loop over observation points (parallelized via TaskManager)
	ngcore::ParallelFor(ngcore::IntRange(MXX), [&](size_t I) {
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
	});
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
	// Uses RadConst::INV_FOUR_PI for 1/(4*pi) factor.
	// =========================================================================
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
	double W = RadConst::INV_FOUR_PI * sigma;

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

//-------------------------------------------------------------------------
// ELF-style optimized functions for HACApK 3DOF tetrahedra
//-------------------------------------------------------------------------

void RadComputeTriangleFaceBasis(
	const TVector3d& V0,
	const TVector3d& V1,
	const TVector3d& V2,
	const TVector3d& elemCentroid,
	RadTriangleFaceBasis& basis)
{
	const double EPS = 1.0e-15;

	basis.valid = false;
	basis.charge_sign = 1.0;
	basis.face_origin = V0;

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
		return;  // Degenerate triangle
	}

	// Normalize to get basis_c (local Z-axis = face normal)
	basis.basis_c.x = normal.x / normalLen;
	basis.basis_c.y = normal.y / normalLen;
	basis.basis_c.z = normal.z / normalLen;

	// Check if normal points outward from element centroid
	TVector3d faceCenter;
	faceCenter.x = (V0.x + V1.x + V2.x) / 3.0;
	faceCenter.y = (V0.y + V1.y + V2.y) / 3.0;
	faceCenter.z = (V0.z + V1.z + V2.z) / 3.0;

	TVector3d outwardVec;
	outwardVec.x = faceCenter.x - elemCentroid.x;
	outwardVec.y = faceCenter.y - elemCentroid.y;
	outwardVec.z = faceCenter.z - elemCentroid.z;

	double dotProduct = basis.basis_c.x * outwardVec.x +
	                    basis.basis_c.y * outwardVec.y +
	                    basis.basis_c.z * outwardVec.z;

	if(dotProduct < 0.0) {
		basis.charge_sign = -1.0;  // Normal points inward, negate charge
	}

	// Build local coordinate system using triple cross product method
	// Same as RadFieldFromTriangleFaceGlobal for consistency
	TVector3d tmp_basis_a = e2;
	TVector3d tmp_basis_b;
	tmp_basis_b.x = e2.x - e1.x * 0.5;
	tmp_basis_b.y = e2.y - e1.y * 0.5;
	tmp_basis_b.z = e2.z - e1.z * 0.5;

	// Triple cross product #1: basis_c = basis_a x basis_b (normalized)
	TVector3d tmp_c;
	tmp_c.x = tmp_basis_a.y * tmp_basis_b.z - tmp_basis_a.z * tmp_basis_b.y;
	tmp_c.y = tmp_basis_a.z * tmp_basis_b.x - tmp_basis_a.x * tmp_basis_b.z;
	tmp_c.z = tmp_basis_a.x * tmp_basis_b.y - tmp_basis_a.y * tmp_basis_b.x;
	double cLen = std::sqrt(tmp_c.x*tmp_c.x + tmp_c.y*tmp_c.y + tmp_c.z*tmp_c.z);
	if(cLen < EPS) return;
	tmp_c.x /= cLen; tmp_c.y /= cLen; tmp_c.z /= cLen;

	// Triple cross product #2: basis_a = basis_b x basis_c (normalized)
	TVector3d new_a;
	new_a.x = tmp_basis_b.y * tmp_c.z - tmp_basis_b.z * tmp_c.y;
	new_a.y = tmp_basis_b.z * tmp_c.x - tmp_basis_b.x * tmp_c.z;
	new_a.z = tmp_basis_b.x * tmp_c.y - tmp_basis_b.y * tmp_c.x;
	double aLen = std::sqrt(new_a.x*new_a.x + new_a.y*new_a.y + new_a.z*new_a.z);
	if(aLen < EPS) return;
	basis.basis_a.x = new_a.x / aLen;
	basis.basis_a.y = new_a.y / aLen;
	basis.basis_a.z = new_a.z / aLen;

	// Triple cross product #3: basis_b = basis_c x basis_a (normalized)
	TVector3d new_b;
	new_b.x = tmp_c.y * basis.basis_a.z - tmp_c.z * basis.basis_a.y;
	new_b.y = tmp_c.z * basis.basis_a.x - tmp_c.x * basis.basis_a.z;
	new_b.z = tmp_c.x * basis.basis_a.y - tmp_c.y * basis.basis_a.x;
	double bLen = std::sqrt(new_b.x*new_b.x + new_b.y*new_b.y + new_b.z*new_b.z);
	if(bLen < EPS) return;
	basis.basis_b.x = new_b.x / bLen;
	basis.basis_b.y = new_b.y / bLen;
	basis.basis_b.z = new_b.z / bLen;

	// Use tmp_c as the local Z-axis (may differ from basis.basis_c in sign)
	basis.basis_c = tmp_c;

	// Convert vertices to local 2D coordinates
	basis.local_v[0] = TVector2d(0.0, 0.0);

	TVector3d d1 = V1 - V0;
	basis.local_v[1].x = d1.x*basis.basis_a.x + d1.y*basis.basis_a.y + d1.z*basis.basis_a.z;
	basis.local_v[1].y = d1.x*basis.basis_b.x + d1.y*basis.basis_b.y + d1.z*basis.basis_b.z;

	TVector3d d2 = V2 - V0;
	basis.local_v[2].x = d2.x*basis.basis_a.x + d2.y*basis.basis_a.y + d2.z*basis.basis_a.z;
	basis.local_v[2].y = d2.x*basis.basis_b.x + d2.y*basis.basis_b.y + d2.z*basis.basis_b.z;

	// ELF optimization: pre-compute edge parameters (saves recomputation for each M direction)
	const double ONE = 1.0;
	const double ZER = 0.0;
	const double EDGE_EPS = 1.0e-20;

	basis.EPSG = 0.0;

	for(int J = 0; J < 3; J++) {
		int L = (J + 1) % 3;
		double XS1 = basis.local_v[L].x - basis.local_v[J].x;
		double XS2 = basis.local_v[L].y - basis.local_v[J].y;

		if(std::abs(XS1) < EDGE_EPS) XS1 = EDGE_EPS;

		basis.DS[J] = std::sqrt(XS2*XS2 + XS1*XS1);
		basis.AM[J] = XS2 / XS1;
		basis.XD[J] = -XS1 / basis.DS[J];
		basis.YD[J] = XS2 / basis.DS[J];

		basis.EPSG = std::max(basis.EPSG, basis.DS[J]);
	}

	basis.EPSG = basis.EPSG * 1.0e-12;

	// Triangle: set 4th edge to dummy values
	basis.DS[3] = ONE;
	basis.AM[3] = ZER;
	basis.XD[3] = ZER;
	basis.YD[3] = ZER;

	basis.valid = true;
}

TVector3d RadFieldFromTriangleFaceWithBasis(
	const RadTriangleFaceBasis& basis,
	double sigma,
	const TVector3d& obsPoint)
{
	if(!basis.valid || std::abs(sigma) < 1.0e-15) {
		return TVector3d(0., 0., 0.);
	}

	// ELF-style optimized field computation using pre-computed edge parameters
	// No edge parameter recomputation - uses basis.DS, basis.AM, basis.XD, basis.YD
	const double EPS = 1.0e-20;

	const TVector3d& AA = basis.basis_a;
	const TVector3d& BB = basis.basis_b;
	const TVector3d& CC = basis.basis_c;
	const TVector3d& YY = basis.face_origin;

	// Weight: sigma / (4*pi)
	double W = RadConst::INV_FOUR_PI * sigma;

	// Transform observation point to local coordinates
	TVector3d DD = obsPoint - YY;
	double EE1 = DD.x*AA.x + DD.y*AA.y + DD.z*AA.z;
	double EE2 = DD.x*BB.x + DD.y*BB.y + DD.z*BB.z;
	double EE3 = DD.x*CC.x + DD.y*CC.y + DD.z*CC.z;

	// Distances from observation point to vertices
	double X[4], Y[4], H[4], E[4], R[4];

	for(int J = 0; J < 3; J++) {
		X[J] = EE1 - basis.local_v[J].x;
		Y[J] = EE2 - basis.local_v[J].y;
		H[J] = Y[J] * X[J];
	}
	X[3] = X[0]; Y[3] = Y[0]; H[3] = H[0];

	double Z = EE3;
	double Z2 = Z * Z;

	for(int J = 0; J < 4; J++) {
		E[J] = Z2 + X[J]*X[J];
		R[J] = std::sqrt(X[J]*X[J] + Y[J]*Y[J] + Z*Z);
	}

	// Edge contributions (use pre-computed DS from basis)
	double RM[4], RP[4], RR[4], AL[4];

	for(int J = 0; J < 4; J++) {
		int JP1 = (J + 1) % 4;
		RM[J] = R[J] + R[JP1] - basis.DS[J];
		RP[J] = R[J] + R[JP1] + basis.DS[J];
		RR[J] = std::max(RM[J] / RP[J], EPS);
		AL[J] = std::log(RR[J]);
	}

	// Field components in local frame (use pre-computed XD, YD from basis)
	double HH1 = W * (-basis.YD[0]*AL[0] - basis.YD[1]*AL[1] - basis.YD[2]*AL[2] - basis.YD[3]*AL[3]);
	double HH2 = W * (-basis.XD[0]*AL[0] - basis.XD[1]*AL[1] - basis.XD[2]*AL[2] - basis.XD[3]*AL[3]);
	double HH3 = 0.0;

	// Z-component (solid angle) - use pre-computed EPSG and AM from basis
	if(std::abs(Z) > basis.EPSG) {
		double ZR[4], AT[4], BT[4];

		for(int J = 0; J < 4; J++) {
			ZR[J] = Z * R[J];
		}

		for(int J = 0; J < 4; J++) {
			int JP1 = (J + 1) % 4;
			AT[J] = (basis.AM[J]*E[J] - H[J]) / ZR[J];
			BT[J] = (basis.AM[J]*E[JP1] - H[JP1]) / ZR[JP1];
		}

		// Triangle: 4th term = 0
		AT[3] = 0.0;
		BT[3] = 0.0;

		HH3 = W * (-std::atan(AT[0]) - std::atan(AT[1]) - std::atan(AT[2]) - std::atan(AT[3])
		          + std::atan(BT[0]) + std::atan(BT[1]) + std::atan(BT[2]) + std::atan(BT[3]));
	}

	// Transform field back to global coordinates
	TVector3d result;
	result.x = HH1*AA.x + HH2*BB.x + HH3*CC.x;
	result.y = HH1*AA.y + HH2*BB.y + HH3*CC.y;
	result.z = HH1*AA.z + HH2*BB.z + HH3*CC.z;

	return result;
}

//-------------------------------------------------------------------------
/**
 * Compute vector potential A from a triangular face using face integration
 *
 * This function computes the vector potential from a uniformly magnetized
 * volume element using face integration (NOT dipole approximation).
 *
 * Formula (matching radTRecCur::B_comp):
 *   A(r) = (1/4pi) * M x BufVect
 *
 * where BufVect = integral_S n/|r-r'| dS is the surface integral vector.
 *
 * Note: This follows Radia's internal convention which does NOT include
 * the mu_0 factor. The relation B = curl(A) uses Radia's internal units.
 *
 * @param V0, V1, V2     Triangle vertices in GLOBAL 3D coordinates
 * @param M              Magnetization vector (in global coordinates)
 * @param obsPoint       Observation point in GLOBAL 3D coordinates
 * @param elemCentroid   Element centroid for outward normal check
 * @return               A field at observation point (in global coordinates)
 */
TVector3d RadVectorPotentialFromTriangleFaceGlobal(
	const TVector3d& V0,
	const TVector3d& V1,
	const TVector3d& V2,
	const TVector3d& M,
	const TVector3d& obsPoint,
	const TVector3d& elemCentroid)
{
	const double EPS = 1.0e-15;
	// Vector potential: A = (mu_0/4pi) * (M x BufVect)
	// Radia uses SI units (B in Tesla, H in A/m), so A should be in T*m
	const double dConst2 = 1.0e-7;  // mu_0/(4*pi)

	// Compute face edges
	TVector3d e1 = V1 - V0;
	TVector3d e2 = V2 - V0;

	// Compute face normal from cross product
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

	// Check if normal points outward
	TVector3d faceCenter;
	faceCenter.x = (V0.x + V1.x + V2.x) / 3.0;
	faceCenter.y = (V0.y + V1.y + V2.y) / 3.0;
	faceCenter.z = (V0.z + V1.z + V2.z) / 3.0;

	TVector3d outwardVec = faceCenter - elemCentroid;
	double dotProduct = normal.x * outwardVec.x + normal.y * outwardVec.y + normal.z * outwardVec.z;
	double signFactor = (dotProduct < 0.0) ? -1.0 : 1.0;

	// Build local coordinate system using face edge (same as FieldFromChargedTriangle)
	// basis_a = e1 normalized (local X)
	// basis_c = face normal (local Z)
	// basis_b = basis_c x basis_a (local Y)

	TVector3d basis_c = normal;
	if(signFactor < 0.0) {
		// Flip normal to point outward
		basis_c.x = -basis_c.x;
		basis_c.y = -basis_c.y;
		basis_c.z = -basis_c.z;
	}

	TVector3d basis_a = e1;
	double aLen = std::sqrt(basis_a.x*basis_a.x + basis_a.y*basis_a.y + basis_a.z*basis_a.z);
	if(aLen < EPS) return TVector3d(0., 0., 0.);
	basis_a.x /= aLen; basis_a.y /= aLen; basis_a.z /= aLen;

	TVector3d basis_b;
	basis_b.x = basis_c.y * basis_a.z - basis_c.z * basis_a.y;
	basis_b.y = basis_c.z * basis_a.x - basis_c.x * basis_a.z;
	basis_b.z = basis_c.x * basis_a.y - basis_c.y * basis_a.x;
	double bLen = std::sqrt(basis_b.x*basis_b.x + basis_b.y*basis_b.y + basis_b.z*basis_b.z);
	if(bLen < EPS) return TVector3d(0., 0., 0.);
	basis_b.x /= bLen; basis_b.y /= bLen; basis_b.z /= bLen;

	TVector3d AA = basis_a;  // Local X
	TVector3d BB = basis_b;  // Local Y
	TVector3d CC = basis_c;  // Local Z (outward normal)

	// Transform vertices to local 2D coordinates (V0 = origin)
	double xy0_x = 0.0, xy0_y = 0.0;
	double xy1_x = e1.x*AA.x + e1.y*AA.y + e1.z*AA.z;
	double xy1_y = e1.x*BB.x + e1.y*BB.y + e1.z*BB.z;
	double xy2_x = e2.x*AA.x + e2.y*AA.y + e2.z*AA.z;
	double xy2_y = e2.x*BB.x + e2.y*BB.y + e2.z*BB.z;

	// Transform observation point to local coordinates
	TVector3d DD = obsPoint - V0;
	double X = DD.x*AA.x + DD.y*AA.y + DD.z*AA.z;
	double Y = DD.x*BB.x + DD.y*BB.y + DD.z*BB.z;
	double Z = DD.x*CC.x + DD.y*CC.y + DD.z*CC.z;

	// Distance from observation point to each vertex
	double x0 = X - xy0_x, y0 = Y - xy0_y;
	double x1 = X - xy1_x, y1 = Y - xy1_y;
	double x2 = X - xy2_x, y2 = Y - xy2_y;

	double R0 = std::sqrt(x0*x0 + y0*y0 + Z*Z);
	double R1 = std::sqrt(x1*x1 + y1*y1 + Z*Z);
	double R2 = std::sqrt(x2*x2 + y2*y2 + Z*Z);

	// Compute I_scalar = integral_face 1/|r-r'| dS using Wilton formula
	// Reference: Wilton, Rao, Glisson, IEEE Trans. AP, Vol. 32, No. 3, March 1984
	//
	// For each edge from vertex A to vertex B:
	//   I_edge = P0 * log((R_plus + s_plus)/(R_minus + s_minus))
	//          - |z| * [atan2(P0*s_plus, R0^2 + |z|*R_plus)
	//                 - atan2(P0*s_minus, R0^2 + |z|*R_minus)]
	// where:
	//   P0 = perpendicular distance from obs projection to edge (positive INWARD)
	//   R0 = sqrt(P0^2 + z^2)
	//   s_plus/s_minus = projections along edge direction
	//   R_plus/R_minus = 3D distances to edge vertices

	double I_scalar = 0.0;
	double absZ = std::abs(Z);
	double z2 = Z * Z;

	// Define 2D vertices array for edge loop
	double verts_x[3] = {xy0_x, xy1_x, xy2_x};
	double verts_y[3] = {xy0_y, xy1_y, xy2_y};
	double R_verts[3] = {R0, R1, R2};

	for(int i = 0; i < 3; i++) {
		int j = (i + 1) % 3;

		// Edge from vertex i to vertex j
		double edge_x = verts_x[j] - verts_x[i];
		double edge_y = verts_y[j] - verts_y[i];
		double L = std::sqrt(edge_x*edge_x + edge_y*edge_y);
		if(L < EPS) continue;

		// Unit tangent along edge
		double tx = edge_x / L;
		double ty = edge_y / L;

		// Inward normal (LEFT perpendicular for CCW vertices)
		double mx = -ty;
		double my = tx;

		// P0 = perpendicular distance from obs to edge line (positive = inward side)
		double P0 = (X - verts_x[i]) * mx + (Y - verts_y[i]) * my;

		// R0_edge = perpendicular distance from obs point to edge LINE (in 3D)
		double R0_edge = std::sqrt(P0*P0 + z2);

		// s_plus = projection of (obs - vertex_i) onto edge direction
		double s_plus = (X - verts_x[i]) * tx + (Y - verts_y[i]) * ty;

		// s_minus = projection of (obs - vertex_j) onto edge direction
		double s_minus = (X - verts_x[j]) * tx + (Y - verts_y[j]) * ty;

		// R_plus = |obs - vertex_i| in 3D, R_minus = |obs - vertex_j| in 3D
		double R_plus = R_verts[i];
		double R_minus = R_verts[j];

		// Log term: P0 * log((R_plus + s_plus) / (R_minus + s_minus))
		double log_term = 0.0;
		double num = R_plus + s_plus;
		double den = R_minus + s_minus;
		if(std::abs(num) > EPS && std::abs(den) > EPS) {
			log_term = P0 * std::log(std::abs(num / den));
		}

		// Arctan term: |z| * [atan2(...) - atan2(...)]
		double atan_term = 0.0;
		if(absZ > EPS && std::abs(R0_edge) > EPS) {
			double R0_sq = R0_edge * R0_edge;
			double atan_plus = std::atan2(P0 * s_plus, R0_sq + absZ * R_plus);
			double atan_minus = std::atan2(P0 * s_minus, R0_sq + absZ * R_minus);
			atan_term = absZ * (atan_plus - atan_minus);
		}

		// Edge contribution
		I_scalar += log_term - atan_term;
	}

	// BufVect = n * I_scalar (contribution from this face)
	// where n is the outward normal (basis_c after sign correction)
	TVector3d BufVect;
	BufVect.x = CC.x * I_scalar;
	BufVect.y = CC.y * I_scalar;
	BufVect.z = CC.z * I_scalar;

	// A = (1/4pi) * (M x BufVect) - matching radTRecCur::B_comp formula
	TVector3d MxBuf;
	MxBuf.x = M.y * BufVect.z - M.z * BufVect.y;
	MxBuf.y = M.z * BufVect.x - M.x * BufVect.z;
	MxBuf.z = M.x * BufVect.y - M.y * BufVect.x;

	TVector3d A;
	A.x = dConst2 * MxBuf.x;
	A.y = dConst2 * MxBuf.y;
	A.z = dConst2 * MxBuf.z;

	return A;
}

//-------------------------------------------------------------------------
/**
 * Compute scalar potential Phi from a triangular face using face integration
 *
 * Same as RadVectorPotentialFromTriangleFaceGlobal but computes:
 *   Phi = (1/4pi) * (M dot BufVect)
 * instead of:
 *   A = (1/4pi) * (M x BufVect)
 */
double RadScalarPotentialFromTriangleFaceGlobal(
	const TVector3d& V0,
	const TVector3d& V1,
	const TVector3d& V2,
	const TVector3d& M,
	const TVector3d& obsPoint,
	const TVector3d& elemCentroid)
{
	const double EPS = 1.0e-15;
	const double dConst2 = RadConst::INV_FOUR_PI;  // 1/(4*pi) - matching radTRecCur

	// Compute face edges
	TVector3d e1 = V1 - V0;
	TVector3d e2 = V2 - V0;

	// Compute face normal from cross product
	TVector3d normal;
	normal.x = e1.y * e2.z - e1.z * e2.y;
	normal.y = e1.z * e2.x - e1.x * e2.z;
	normal.z = e1.x * e2.y - e1.y * e2.x;

	double normalLen = std::sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
	if(normalLen < EPS) {
		return 0.0;  // Degenerate triangle
	}

	// Normalize
	normal.x /= normalLen;
	normal.y /= normalLen;
	normal.z /= normalLen;

	// Check if normal points outward
	TVector3d faceCenter;
	faceCenter.x = (V0.x + V1.x + V2.x) / 3.0;
	faceCenter.y = (V0.y + V1.y + V2.y) / 3.0;
	faceCenter.z = (V0.z + V1.z + V2.z) / 3.0;

	TVector3d outwardVec = faceCenter - elemCentroid;
	double dotProduct = normal.x * outwardVec.x + normal.y * outwardVec.y + normal.z * outwardVec.z;
	double signFactor = (dotProduct < 0.0) ? -1.0 : 1.0;

	// Build local coordinate system using face edge (same as FieldFromChargedTriangle)
	TVector3d basis_c = normal;
	if(signFactor < 0.0) {
		// Flip normal to point outward
		basis_c.x = -basis_c.x;
		basis_c.y = -basis_c.y;
		basis_c.z = -basis_c.z;
	}

	TVector3d basis_a = e1;
	double aLen = std::sqrt(basis_a.x*basis_a.x + basis_a.y*basis_a.y + basis_a.z*basis_a.z);
	if(aLen < EPS) return 0.0;
	basis_a.x /= aLen; basis_a.y /= aLen; basis_a.z /= aLen;

	TVector3d basis_b;
	basis_b.x = basis_c.y * basis_a.z - basis_c.z * basis_a.y;
	basis_b.y = basis_c.z * basis_a.x - basis_c.x * basis_a.z;
	basis_b.z = basis_c.x * basis_a.y - basis_c.y * basis_a.x;
	double bLen = std::sqrt(basis_b.x*basis_b.x + basis_b.y*basis_b.y + basis_b.z*basis_b.z);
	if(bLen < EPS) return 0.0;
	basis_b.x /= bLen; basis_b.y /= bLen; basis_b.z /= bLen;

	TVector3d AA = basis_a;  // Local X
	TVector3d BB = basis_b;  // Local Y
	TVector3d CC = basis_c;  // Local Z (outward normal)

	// Transform vertices to local 2D coordinates (V0 = origin)
	double xy0_x = 0.0, xy0_y = 0.0;
	double xy1_x = e1.x*AA.x + e1.y*AA.y + e1.z*AA.z;
	double xy1_y = e1.x*BB.x + e1.y*BB.y + e1.z*BB.z;
	double xy2_x = e2.x*AA.x + e2.y*AA.y + e2.z*AA.z;
	double xy2_y = e2.x*BB.x + e2.y*BB.y + e2.z*BB.z;

	// Transform observation point to local coordinates
	TVector3d DD = obsPoint - V0;
	double X = DD.x*AA.x + DD.y*AA.y + DD.z*AA.z;
	double Y = DD.x*BB.x + DD.y*BB.y + DD.z*BB.z;
	double Z = DD.x*CC.x + DD.y*CC.y + DD.z*CC.z;

	// Distance from observation point to each vertex
	double x0 = X - xy0_x, y0 = Y - xy0_y;
	double x1 = X - xy1_x, y1 = Y - xy1_y;
	double x2 = X - xy2_x, y2 = Y - xy2_y;

	double R0 = std::sqrt(x0*x0 + y0*y0 + Z*Z);
	double R1 = std::sqrt(x1*x1 + y1*y1 + Z*Z);
	double R2 = std::sqrt(x2*x2 + y2*y2 + Z*Z);

	// Compute I_scalar = integral_face 1/|r-r'| dS using Wilton formula
	// Reference: Wilton, Rao, Glisson, IEEE Trans. AP, Vol. 32, No. 3, March 1984
	//
	// For each edge from vertex A to vertex B:
	//   I_edge = P0 * log((R_plus + s_plus)/(R_minus + s_minus))
	//          - |z| * [atan2(P0*s_plus, R0^2 + |z|*R_plus)
	//                 - atan2(P0*s_minus, R0^2 + |z|*R_minus)]
	// where:
	//   P0 = perpendicular distance from obs projection to edge (positive INWARD)
	//   R0 = sqrt(P0^2 + z^2)
	//   s_plus/s_minus = projections along edge direction
	//   R_plus/R_minus = 3D distances to edge vertices

	double I_scalar = 0.0;
	double absZ = std::abs(Z);
	double z2 = Z * Z;

	// Define 2D vertices array for edge loop
	double verts_x[3] = {xy0_x, xy1_x, xy2_x};
	double verts_y[3] = {xy0_y, xy1_y, xy2_y};
	double R_verts[3] = {R0, R1, R2};

	for(int i = 0; i < 3; i++) {
		int j = (i + 1) % 3;

		// Edge from vertex i to vertex j
		double edge_x = verts_x[j] - verts_x[i];
		double edge_y = verts_y[j] - verts_y[i];
		double L = std::sqrt(edge_x*edge_x + edge_y*edge_y);
		if(L < EPS) continue;

		// Unit tangent along edge
		double tx = edge_x / L;
		double ty = edge_y / L;

		// Inward normal (LEFT perpendicular for CCW vertices)
		double mx = -ty;
		double my = tx;

		// P0 = perpendicular distance from obs to edge line (positive = inward side)
		double P0 = (X - verts_x[i]) * mx + (Y - verts_y[i]) * my;

		// R0 = perpendicular distance from obs point to edge LINE (in 3D)
		double R0_edge = std::sqrt(P0*P0 + z2);

		// s_plus = projection of (obs - vertex_i) onto edge direction
		double s_plus = (X - verts_x[i]) * tx + (Y - verts_y[i]) * ty;

		// s_minus = projection of (obs - vertex_j) onto edge direction
		double s_minus = (X - verts_x[j]) * tx + (Y - verts_y[j]) * ty;

		// R_plus = |obs - vertex_i| in 3D, R_minus = |obs - vertex_j| in 3D
		double R_plus = R_verts[i];
		double R_minus = R_verts[j];

		// Log term: P0 * log((R_plus + s_plus) / (R_minus + s_minus))
		double log_term = 0.0;
		double num = R_plus + s_plus;
		double den = R_minus + s_minus;
		if(std::abs(num) > EPS && std::abs(den) > EPS) {
			log_term = P0 * std::log(std::abs(num / den));
		}

		// Arctan term: |z| * [atan2(...) - atan2(...)]
		double atan_term = 0.0;
		if(absZ > EPS && std::abs(R0_edge) > EPS) {
			double R0_sq = R0_edge * R0_edge;
			double atan_plus = std::atan2(P0 * s_plus, R0_sq + absZ * R_plus);
			double atan_minus = std::atan2(P0 * s_minus, R0_sq + absZ * R_minus);
			atan_term = absZ * (atan_plus - atan_minus);
		}

		// Edge contribution
		I_scalar += log_term - atan_term;
	}

	// BufVect = n * I_scalar (contribution from this face)
	TVector3d BufVect;
	BufVect.x = CC.x * I_scalar;
	BufVect.y = CC.y * I_scalar;
	BufVect.z = CC.z * I_scalar;

	// Phi = (1/4pi) * (M dot BufVect) - same formula as radTRecCur but with dot product
	double M_dot_BufVect = M.x * BufVect.x + M.y * BufVect.y + M.z * BufVect.z;

	return dConst2 * M_dot_BufVect;
}

void RadDemagTensorFromTetrahedron(
	const TVector3d* face_vertices,
	const TVector3d& elemCentroid,
	const TVector3d& obsPoint,
	double* N_mat)
{
	// Initialize output to zero
	for(int i = 0; i < 9; i++) {
		N_mat[i] = 0.0;
	}

	// Unit magnetization vectors
	TVector3d unit_M[3] = {
		TVector3d(1.0, 0.0, 0.0),
		TVector3d(1.0, 0.0, 0.0),  // Will be set properly below
		TVector3d(1.0, 0.0, 0.0)
	};
	unit_M[0] = TVector3d(1.0, 0.0, 0.0);
	unit_M[1] = TVector3d(0.0, 1.0, 0.0);
	unit_M[2] = TVector3d(0.0, 0.0, 1.0);

	// Process each of the 4 triangular faces
	for(int f = 0; f < 4; f++) {
		// Get face vertices
		const TVector3d& V0 = face_vertices[f * 3 + 0];
		const TVector3d& V1 = face_vertices[f * 3 + 1];
		const TVector3d& V2 = face_vertices[f * 3 + 2];

		// Compute face basis ONCE per face (ELF optimization)
		RadTriangleFaceBasis basis;
		RadComputeTriangleFaceBasis(V0, V1, V2, elemCentroid, basis);

		if(!basis.valid) continue;

		// Compute surface charge for each magnetization direction
		// sigma_j = M_j dot n (with sign correction)
		// The face normal is basis.basis_c, but we need the ORIGINAL normal
		// for computing M dot n. Recompute it.
		TVector3d e1 = V1 - V0;
		TVector3d e2 = V2 - V0;
		TVector3d faceNormal;
		faceNormal.x = e1.y * e2.z - e1.z * e2.y;
		faceNormal.y = e1.z * e2.x - e1.x * e2.z;
		faceNormal.z = e1.x * e2.y - e1.y * e2.x;
		double nLen = std::sqrt(faceNormal.x*faceNormal.x +
		                        faceNormal.y*faceNormal.y +
		                        faceNormal.z*faceNormal.z);
		if(nLen > 1.0e-15) {
			faceNormal.x /= nLen;
			faceNormal.y /= nLen;
			faceNormal.z /= nLen;
		}

		// Apply charge sign correction (for inward-pointing normals)
		double signFactor = basis.charge_sign;

		// Compute H field for each unit magnetization direction
		for(int j = 0; j < 3; j++) {
			// Surface charge: sigma = M dot n * sign_correction
			double sigma = (unit_M[j].x * faceNormal.x +
			               unit_M[j].y * faceNormal.y +
			               unit_M[j].z * faceNormal.z) * signFactor;

			if(std::abs(sigma) < 1.0e-15) continue;

			// Compute field using pre-computed basis
			TVector3d H = RadFieldFromTriangleFaceWithBasis(basis, sigma, obsPoint);

			// Accumulate into tensor: N_mat[i*3+j] = dH_i / dM_j
			N_mat[0*3 + j] += H.x;
			N_mat[1*3 + j] += H.y;
			N_mat[2*3 + j] += H.z;
		}
	}
}
