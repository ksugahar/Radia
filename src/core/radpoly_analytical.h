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

#endif
