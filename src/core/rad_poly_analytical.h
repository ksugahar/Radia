/*-------------------------------------------------------------------------
*
* File name:      radpoly_analytical.h
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

#ifndef __RADPOLY_ANALYTICAL_H
#define __RADPOLY_ANALYTICAL_H

#include "gmvect.h"
#include <vector>

//-------------------------------------------------------------------------
/**
 * Compute field from polygon magnetic charge using analytical formula
 *
 * This function implements an analytical formula for computing the
 * magnetic field from a planar polygon with uniform magnetic charge density.
 *
 * The formula is exact and based on solid angle integration:
 *   H = (sigma/4pi) integral dOmega
 *
 * where sigma is the magnetic charge density and dOmega is the solid angle
 * subtended by each edge of the polygon.
 *
 * @param AA Local X-axis (tangent to polygon plane)
 * @param BB Local Y-axis (tangent to polygon plane)
 * @param CC Local Z-axis (normal to polygon plane)
 * @param YY Reference point on polygon (typically center or first vertex)
 * @param XY Polygon vertices in local 2D coordinates
 * @param XX Observation points in global 3D coordinates
 * @param FGH Output: accumulated magnetic field at each observation point
 * @param W Magnetic charge density (weight)
 * @param NII Element index (for error reporting)
 * @param KAdo Number of polygon vertices (3 for triangle, 4 for quad)
 */
void RadAnalyticalFieldFromPolygonCharge(
	const TVector3d& AA,
	const TVector3d& BB,
	const TVector3d& CC,
	const TVector3d& YY,
	const std::vector<TVector2d>& XY,
	const std::vector<TVector3d>& XX,
	std::vector<TVector3d>& FGH,
	double W,
	int NII,
	int KAdo);

/**
 * Compute field from triangular magnetic charge
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
	int NII);

/**
 * Compute field from quadrilateral magnetic charge
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
	int NII);

/**
 * Compute H field from triangular face using GLOBAL 3D coordinates directly
 *
 * This function computes the magnetic field from a triangular face with
 * surface charge density sigma = M dot n, using global coordinates only.
 * No local coordinate transformation is needed.
 *
 * IMPORTANT: The function ensures the normal vector points OUTWARD from the
 * element centroid. If the computed normal points inward, the surface charge
 * sign is negated.
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
	const TVector3d& elemCentroid);

/**
 * Compute vector potential A from a triangular face using face integration
 *
 * This function computes the vector potential from a uniformly magnetized
 * volume element using face integration (NOT dipole approximation).
 *
 * Formula:
 *   A(r) = (mu_0 / 4*pi) * integral_S (M x n) / |r-r'| dS
 *
 * The integral over each face uses the Wilton et al. formula with edge
 * logarithms and solid angle terms.
 *
 * @param V0, V1, V2    Triangle vertices in GLOBAL 3D coordinates
 * @param M             Magnetization vector (in global coordinates)
 * @param obsPoint      Observation point in GLOBAL 3D coordinates
 * @param elemCentroid  Element centroid for outward normal check
 * @return              A field at observation point (in global coordinates)
 */
TVector3d RadVectorPotentialFromTriangleFaceGlobal(
	const TVector3d& V0,
	const TVector3d& V1,
	const TVector3d& V2,
	const TVector3d& M,
	const TVector3d& obsPoint,
	const TVector3d& elemCentroid);

/**
 * Compute scalar potential Phi from a triangular face using face integration
 *
 * This function computes the scalar potential from a uniformly magnetized
 * volume element using face integration (NOT dipole approximation).
 *
 * Formula:
 *   Phi(r) = (1 / 4*pi) * integral_S (M dot n) / |r-r'| dS
 *          = (1 / 4*pi) * M dot BufVect
 *
 * where BufVect = n * integral_face(1/|r-r'|) dS
 *
 * The integral over each face uses the Wilton et al. formula with edge
 * logarithms and solid angle terms.
 *
 * @param V0, V1, V2    Triangle vertices in GLOBAL 3D coordinates
 * @param M             Magnetization vector (in global coordinates)
 * @param obsPoint      Observation point in GLOBAL 3D coordinates
 * @param elemCentroid  Element centroid for outward normal check
 * @return              Phi field at observation point (scalar)
 */
double RadScalarPotentialFromTriangleFaceGlobal(
	const TVector3d& V0,
	const TVector3d& V1,
	const TVector3d& V2,
	const TVector3d& M,
	const TVector3d& obsPoint,
	const TVector3d& elemCentroid);

//-------------------------------------------------------------------------
// ELF-style optimized functions for HACApK 3DOF tetrahedra
// These functions compute the full 3x3 demagnetization tensor in one pass
// by reusing coordinate transformations across all 3 magnetization directions.
//-------------------------------------------------------------------------

/**
 * Pre-computed face basis for fast 3x3 tensor computation
 * Stores all geometric data needed for field computation
 *
 * ELF-style optimization: pre-compute edge parameters (DS, AM, XD, YD)
 * to avoid recomputing them for each magnetization direction.
 */
struct RadTriangleFaceBasis {
	TVector3d basis_a;        // Local X-axis (normalized)
	TVector3d basis_b;        // Local Y-axis (normalized)
	TVector3d basis_c;        // Local Z-axis = face normal (normalized)
	TVector3d face_origin;    // Reference point (V0)
	TVector2d local_v[3];     // Vertices in local 2D coordinates
	double charge_sign;       // +1 if outward normal, -1 if inward
	bool valid;               // True if face is non-degenerate

	// Pre-computed edge parameters (ELF optimization)
	double DS[4];             // Edge lengths
	double AM[4];             // Edge slopes (dy/dx)
	double XD[4];             // Edge direction X
	double YD[4];             // Edge direction Y
	double EPSG;              // Tolerance for z=0 check
};

/**
 * Compute face basis for a triangular face (ELF-style)
 *
 * This function computes all geometric data needed for field computation
 * from a triangular face. The basis is computed ONCE per face and can be
 * reused for all 3 magnetization directions.
 *
 * @param V0, V1, V2    Triangle vertices in GLOBAL 3D coordinates
 * @param elemCentroid  Element centroid for outward normal determination
 * @param basis         Output: pre-computed face basis
 */
void RadComputeTriangleFaceBasis(
	const TVector3d& V0,
	const TVector3d& V1,
	const TVector3d& V2,
	const TVector3d& elemCentroid,
	RadTriangleFaceBasis& basis);

/**
 * Compute H field from a triangular face using pre-computed basis (ELF-style)
 *
 * This function computes the H field at an observation point from a
 * triangular face with given surface charge density, using a pre-computed
 * face basis. This is faster than RadFieldFromTriangleFaceGlobal when
 * computing fields for multiple charge densities (e.g., 3x3 tensor).
 *
 * @param basis         Pre-computed face basis from RadComputeTriangleFaceBasis
 * @param sigma         Surface charge density (M dot n, with sign correction)
 * @param obsPoint      Observation point in GLOBAL 3D coordinates
 * @return              H field at observation point (in global coordinates)
 */
TVector3d RadFieldFromTriangleFaceWithBasis(
	const RadTriangleFaceBasis& basis,
	double sigma,
	const TVector3d& obsPoint);

/**
 * Compute full 3x3 demagnetization tensor from a tetrahedron (ELF-style optimized)
 *
 * This function computes the 3x3 demagnetization tensor N such that:
 *   H = -N * M / (4*pi)
 *
 * It is optimized to compute all 9 tensor elements in one pass by:
 * 1. Computing face basis ONCE per face (4 faces total)
 * 2. Computing field for all 3 magnetization directions using the same basis
 *
 * This is 3x faster than calling RadFieldFromTriangleFaceGlobal 12 times.
 *
 * @param face_vertices Array of 4 faces, each with 3 vertices (12 TVector3d total)
 *                      face_vertices[f*3+v] = vertex v of face f
 * @param elemCentroid  Tetrahedron centroid for outward normal determination
 * @param obsPoint      Observation point (typically element center)
 * @param N_mat         Output: 3x3 demagnetization tensor in row-major order
 *                      N_mat[i*3+j] = dH_i/dM_j
 */
void RadDemagTensorFromTetrahedron(
	const TVector3d* face_vertices,  // 4 faces * 3 vertices = 12 vertices
	const TVector3d& elemCentroid,
	const TVector3d& obsPoint,
	double* N_mat);

#endif
