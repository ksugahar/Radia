#ifndef RAD_FIELD_UNIFIED_H
#define RAD_FIELD_UNIFIED_H

/**
 * @file rad_field_unified.h
 * @brief Unified Field Computation Module
 *
 * Provides a single entry point for field computation that handles:
 * 1. Point classification (inside/outside/near magnetic elements)
 * 2. FMM acceleration (ExaFMM dipole approximation)
 * 3. Direct computation (B_genComp)
 *
 * This module is used by:
 * - rad.Fld() and rad.FldBatch()
 * - rad_particle_trajectory.cpp (beam tracking)
 * - radia_ngsolve.cpp (RadiaField CoefficientFunction)
 *
 * Design Policy:
 * - Points INSIDE magnetic elements return internal magnetization value
 * - Points OUTSIDE use direct computation or FMM
 * - Points NEAR surface are flagged for potential accuracy issues
 *
 * @version 1.6.0
 * @date 2026-01-14
 */

#include "gmvect.h"
#include "rad_point_classify.h"
#include <vector>
#include <cstdint>

// Forward declarations
class radTg3d;
struct radTField;

namespace RadFieldUnified {

/**
 * @brief Field computation method selection
 */
enum ComputeMethod {
    METHOD_AUTO = 0,       // Automatic selection based on problem size
    METHOD_DIRECT = 1,     // Direct computation (B_genComp)
    METHOD_FMM = 2,        // FMM dipole approximation
    METHOD_ADAPTIVE = 3    // Adaptive: FMM for far, direct for near
};

/**
 * @brief Field type specification
 */
enum FieldType {
    FIELD_B = 0,    // Magnetic flux density B [Tesla]
    FIELD_H = 1,    // Magnetic field H [A/m]
    FIELD_A = 2,    // Vector potential A [T*m]
    FIELD_PHI = 3,  // Scalar potential Phi [A]
    FIELD_M = 4     // Magnetization M [A/m] (inside elements only)
};

/**
 * @brief Point status after classification
 */
enum PointStatus {
    STATUS_OUTSIDE = 0,   // Point is outside all elements (normal computation)
    STATUS_INSIDE = 1,    // Point is inside an element (use internal M)
    STATUS_NEAR = 2,      // Point is near surface (accuracy warning)
    STATUS_ERROR = 3      // Computation error
};

/**
 * @brief Result structure for single point field computation
 */
struct FieldResult {
    double Bx, By, Bz;     // B field [Tesla]
    double Hx, Hy, Hz;     // H field [A/m]
    PointStatus status;     // Point classification result
    int element_id;         // ID of containing/nearest element (-1 if none)
};

/**
 * @brief Configuration for unified field computation
 */
struct ComputeConfig {
    ComputeMethod method;          // Computation method
    double fmm_eps;                // FMM tolerance (0 = disabled)
    double near_threshold;         // Distance threshold for "near" classification
    bool check_inside;             // Enable inside/outside classification
    bool return_internal_field;    // If inside, return internal M as field
    bool warn_inside;              // Issue warning when point is inside

    // Default configuration
    ComputeConfig()
        : method(METHOD_AUTO)
        , fmm_eps(0.0)
        , near_threshold(1e-6)
        , check_inside(true)
        , return_internal_field(true)
        , warn_inside(false)
    {}
};

/**
 * @brief Compute field at a single point
 *
 * Main entry point for single-point field computation.
 *
 * @param g3dPtr      Pointer to Radia 3D object (magnetic source)
 * @param point       Evaluation point [x, y, z] in current units
 * @param field_type  Type of field to compute (B, H, A, Phi)
 * @param config      Computation configuration
 * @return            FieldResult with field values and status
 */
FieldResult ComputeFieldSingle(
    radTg3d* g3dPtr,
    const TVector3d& point,
    FieldType field_type,
    const ComputeConfig& config = ComputeConfig()
);

/**
 * @brief Compute field at multiple points (batch)
 *
 * Efficient batch computation with OpenMP parallelization.
 *
 * @param g3dPtr      Pointer to Radia 3D object (magnetic source)
 * @param points      Array of evaluation points [x0,y0,z0, x1,y1,z1, ...]
 * @param n_points    Number of points
 * @param field_type  Type of field to compute (B, H, A, Phi)
 * @param config      Computation configuration
 * @param B_out       Output: B field values [Bx0,By0,Bz0, Bx1,By1,Bz1, ...] (may be nullptr)
 * @param H_out       Output: H field values [Hx0,Hy0,Hz0, Hx1,Hy1,Hz1, ...] (may be nullptr)
 * @param status_out  Output: Point status array (may be nullptr)
 */
void ComputeFieldBatch(
    radTg3d* g3dPtr,
    const double* points,
    int n_points,
    FieldType field_type,
    const ComputeConfig& config,
    double* B_out,
    double* H_out,
    PointStatus* status_out = nullptr
);

/**
 * @brief Compute field for particle trajectory integration
 *
 * Specialized function for trajectory computation.
 * Returns B field at the evaluation point, with inside/outside handling.
 *
 * @param g3dPtr      Pointer to Radia 3D object (magnetic source)
 * @param point       Evaluation point [x, y, z]
 * @param B_out       Output: B field [Bx, By, Bz]
 * @return            true if point is valid for trajectory, false if inside magnet
 */
bool ComputeFieldForTrajectory(
    radTg3d* g3dPtr,
    const TVector3d& point,
    TVector3d& B_out
);

/**
 * @brief Check if a point is inside any magnetic element
 *
 * Uses solid angle method for accurate inside/outside determination.
 *
 * @param g3dPtr      Pointer to Radia 3D object
 * @param point       Point to check
 * @param element_id  Output: ID of containing element (-1 if outside)
 * @return            true if point is inside an element
 */
bool IsPointInsideMagnet(
    radTg3d* g3dPtr,
    const TVector3d& point,
    int& element_id
);

/**
 * @brief Get magnetization at a point inside an element
 *
 * For points inside magnetic elements, returns the local magnetization.
 *
 * @param g3dPtr      Pointer to Radia 3D object
 * @param point       Point inside element
 * @param element_id  ID of containing element
 * @param M_out       Output: Magnetization [Mx, My, Mz] in A/m
 * @return            true if successful
 */
bool GetMagnetizationAtPoint(
    radTg3d* g3dPtr,
    const TVector3d& point,
    int element_id,
    TVector3d& M_out
);

/**
 * @brief Initialize FMM data structures for a Radia object
 *
 * Extracts dipole data from all elements for FMM computation.
 * Call this once before batch FMM computations.
 *
 * @param container_handle  Radia container handle
 * @return                  true if initialization successful
 */
bool InitializeFMM(int container_handle);

/**
 * @brief Release FMM data structures
 *
 * @param container_handle  Radia container handle
 */
void ReleaseFMM(int container_handle);

/**
 * @brief Get element data for inside/outside classification
 *
 * Builds element data structures from Radia container.
 *
 * @param container_handle  Radia container handle
 * @param elements          Output: Element data for classification
 * @return                  true if successful
 */
bool BuildElementData(
    int container_handle,
    std::vector<RadPointClassify::ElementData>& elements
);

} // namespace RadFieldUnified

#endif // RAD_FIELD_UNIFIED_H
