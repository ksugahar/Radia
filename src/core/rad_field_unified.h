#ifndef RAD_FIELD_UNIFIED_H
#define RAD_FIELD_UNIFIED_H

/**
 * @file rad_field_unified.h
 * @brief Unified Field Computation Module
 *
 * Provides a single entry point for field computation that handles:
 * 1. Point classification (inside/outside/near magnetic elements)
 * 2. Direct computation (B_genComp)
 *
 * This module is used by:
 * - rad.Fld() (unified single-point and batch)
 * - rad_particle_trajectory.cpp (beam tracking)
 * - radia_ngsolve.cpp (RadiaField CoefficientFunction)
 *
 * Design Policy:
 * - Points INSIDE magnetic elements return internal magnetization value
 * - Points OUTSIDE use direct computation
 * - Points NEAR surface are flagged for potential accuracy issues
 *
 * @version 1.6.0
 * @date 2026-01-14
 */

#include "gmvect.h"
#include "rad_point_classify.h"
#include <vector>
#include <cstdint>
#include <complex>

// Forward declarations
class radTg3d;
struct radTField;

namespace RadFieldUnified {

/**
 * @brief Field computation method selection
 */
enum ComputeMethod {
    METHOD_AUTO = 0,       // Automatic selection based on problem size
    METHOD_DIRECT = 1      // Direct computation (B_genComp)
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
    double near_threshold;         // Distance threshold for "near" classification
    bool check_inside;             // Enable inside/outside classification
    bool return_internal_field;    // If inside, return internal M as field
    bool warn_inside;              // Issue warning when point is inside

    // Default configuration
    ComputeConfig()
        : method(METHOD_AUTO)
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
 * @brief Clear all element caches
 *
 * Must be called when elements are deleted (e.g., UtiDelAll)
 * to prevent stale cache data.
 */
void ClearAllCaches();

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

// ============================================================================
// Complex Field Computation
// ============================================================================

/**
 * @brief Complex field result for AC/frequency domain analysis
 */
struct ComplexFieldResult {
    std::complex<double> Bx, By, Bz;  // B field [Tesla]
    std::complex<double> Hx, Hy, Hz;  // H field [A/m]
    std::complex<double> Ax, Ay, Az;  // Vector potential [T*m]
    PointStatus status;                // Point classification
    int element_id;                    // Containing element (-1 if outside)
};

/**
 * @brief Compute complex B field at a single point
 *
 * For static (DC) analysis, imaginary parts are zero.
 * For AC analysis with sinusoidal magnetization, this handles phase.
 *
 * @param g3dPtr      Pointer to Radia 3D object
 * @param point       Evaluation point [x, y, z]
 * @param M_complex   Optional: complex magnetization values for each element
 *                    Format: [Mx0_re, My0_re, Mz0_re, Mx0_im, My0_im, Mz0_im, ...]
 *                    If nullptr, use static magnetization from elements
 * @param n_elements  Number of elements (only used if M_complex is provided)
 * @param config      Computation configuration
 * @return            ComplexFieldResult with B, H field values
 */
ComplexFieldResult ComputeComplexFieldSingle(
    radTg3d* g3dPtr,
    const TVector3d& point,
    const std::complex<double>* M_complex = nullptr,
    int n_elements = 0,
    const ComputeConfig& config = ComputeConfig()
);

/**
 * @brief Batch compute complex B field at multiple points (OpenMP parallelized)
 *
 * Used by frequency-domain field workflows for efficient field evaluation.
 *
 * @param g3dPtr      Pointer to Radia 3D object
 * @param points      Array of evaluation points [x0,y0,z0, x1,y1,z1, ...]
 * @param n_points    Number of points
 * @param M_complex   Optional: complex magnetization for each element
 * @param n_elements  Number of elements (for M_complex array)
 * @param config      Computation configuration
 * @param B_out       Output: complex B field [Bx0,By0,Bz0, Bx1,By1,Bz1, ...]
 * @param status_out  Output: point status array (may be nullptr)
 */
void ComputeComplexFieldBatch(
    radTg3d* g3dPtr,
    const double* points,
    int n_points,
    const std::complex<double>* M_complex,
    int n_elements,
    const ComputeConfig& config,
    std::complex<double>* B_out,
    PointStatus* status_out = nullptr
);

/**
 * @brief Element face data for surface-charge (face-charge) integration
 *
 * Stores face geometry for accurate near-field computation.
 * For tetrahedra: 4 triangular faces
 * For hexahedra: 6 quadrilateral faces (split into 2 triangles each)
 */
struct ElementFaceData {
    int n_faces;                              // Number of faces (4 for tet, 6 for hex)
    std::vector<std::vector<TVector3d>> face_vertices;  // Vertices per face
    TVector3d centroid;                       // Element centroid
    double characteristic_size;               // Element size for distance check
};

/**
 * @brief Compute B field from compact magnetization sources.
 *
 * ADAPTIVE METHOD: Uses face integration for near-field, dipole for far-field.
 * - Near field (r < 3 * element_size): face integration (high accuracy)
 * - Far field (r >= 3 * element_size): Dipole approximation (fast, <5% error)
 *
 * @param point       Evaluation point
 * @param centers     Element centers [x0,y0,z0, x1,y1,z1, ...]
 * @param volumes     Element volumes
 * @param M_complex   Complex magnetization [Mx0,My0,Mz0, Mx1,My1,Mz1, ...]
 * @param n_elements  Number of elements
 * @param B_out       Output: complex B field [Bx, By, Bz]
 */
void ComputeBFromMagnetization(
    const TVector3d& point,
    const double* centers,
    const double* volumes,
    const std::complex<double>* M_complex,
    int n_elements,
    std::complex<double>* B_out
);

/**
 * @brief Compute B field with error-controlled adaptive method
 *
 * Switches between face integration and dipole approximation based on
 * TARGET ERROR specification. The dipole approximation error scales as:
 *   error ~ (element_size / distance)^3
 *
 * Distance threshold is computed from target error:
 *   r_threshold = element_size / target_error^(1/3)
 *
 * Examples:
 *   target_error = 0.01 (1%)  -> r > 4.6 * element_size uses dipole
 *   target_error = 0.05 (5%)  -> r > 2.7 * element_size uses dipole
 *   target_error = 0.10 (10%) -> r > 2.2 * element_size uses dipole
 *
 * Special values:
 *   target_error = 0.0  -> Always use face integration (exact, slow)
 *   target_error = 1.0  -> Always use dipole (fast, approximate)
 *
 * @param point          Evaluation point
 * @param face_data      Element face data array
 * @param M_complex      Complex magnetization [Mx0,My0,Mz0, Mx1,My1,Mz1, ...]
 * @param n_elements     Number of elements
 * @param target_error   Target relative error (0.0 to 1.0, default: 0.05 = 5%)
 * @param B_out          Output: complex B field [Bx, By, Bz]
 */
void ComputeBFromMagnetizationAdaptive(
    const TVector3d& point,
    const ElementFaceData* face_data,
    const std::complex<double>* M_complex,
    int n_elements,
    double target_error,
    std::complex<double>* B_out
);

/**
 * @brief Convert target error to distance threshold factor
 *
 * Based on dipole approximation error: error ~ (a/r)^3
 * -> r/a = 1 / error^(1/3)
 *
 * @param target_error  Target relative error (0.0 to 1.0)
 * @return              Distance threshold as multiple of element size
 */
inline double ErrorToThresholdFactor(double target_error) {
    if (target_error <= 0.0) return 1e10;  // Always face integration
    if (target_error >= 1.0) return 0.0;   // Always dipole
    return 1.0 / std::pow(target_error, 1.0/3.0);
}

/**
 * @brief Build element face data from Radia 3D object
 *
 * Extracts face geometry from all elements for near-field integration.
 *
 * @param g3dPtr     Pointer to Radia 3D object
 * @param face_data  Output: face data array (one per element)
 * @return           Number of elements processed
 */
int BuildElementFaceData(
    radTg3d* g3dPtr,
    std::vector<ElementFaceData>& face_data
);

/**
 * @brief Check if point is inside any element (with element data cache)
 *
 * Shared by static and frequency-domain field routines for inside/outside
 * classification. Uses solid angle method for accurate determination.
 *
 * @param point           Point to check
 * @param elements        Pre-built element data (from BuildElementData)
 * @param containing_elem Output: index of containing element (-1 if outside)
 * @return                true if point is inside any element
 */
bool IsPointInsideAnyElement(
    const TVector3d& point,
    const std::vector<RadPointClassify::ElementData>& elements,
    int& containing_elem
);

/**
 * @brief Get magnetization at point inside an element
 *
 * For points inside magnetic elements, returns the element's magnetization.
 * When provided, this uses the complex solution magnetization.
 *
 * @param containing_elem Index of containing element
 * @param M_complex       Complex magnetization array (nullptr for static)
 * @param n_elements      Number of elements
 * @param M_out           Output: magnetization [Mx, My, Mz] (complex)
 * @return                true if successful
 */
bool GetMagnetizationInElement(
    int containing_elem,
    const std::complex<double>* M_complex,
    int n_elements,
    std::complex<double>* M_out
);

// ============================================================================
// PEEC Conductor Field Computation (Biot-Savart)
// ============================================================================

/**
 * @brief Conductor segment data for Biot-Savart field computation
 *
 * Represents a filament segment carrying current.
 * Used for Biot-Savart field computation.
 */
struct ConductorSegment {
    TVector3d start;    // Segment start point
    TVector3d end;      // Segment end point
    TVector3d center;   // Segment center
    TVector3d tangent;  // Unit tangent vector (direction of current flow)
    double length;      // Segment length
};

/**
 * @brief Compute B field from conductor current using Biot-Savart
 *
 * Uses Laplace kernel: B = (mu0/4pi) * I * integral{ dl x r / r^3 }
 *
 *
 * @param point        Evaluation point
 * @param segments     Array of conductor segments
 * @param n_segments   Number of segments
 * @param current      Complex current [A]
 * @param B_out        Output: complex B field [Bx, By, Bz]
 */
void ComputeBFromConductor(
    const TVector3d& point,
    const ConductorSegment* segments,
    int n_segments,
    std::complex<double> current,
    std::complex<double>* B_out
);

/**
 * @brief Batch compute B field from conductor
 *
 * @param points       Evaluation points [x0,y0,z0, ...]
 * @param n_points     Number of points
 * @param segments     Conductor segments
 * @param n_segments   Number of segments
 * @param current      Complex current [A]
 * @param B_out        Output: complex B field [Bx0,By0,Bz0, ...]
 */
void ComputeBFromConductorBatch(
    const double* points,
    int n_points,
    const ConductorSegment* segments,
    int n_segments,
    std::complex<double> current,
    std::complex<double>* B_out
);

/**
 * @brief Combined conductor and magnetization field computation.
 *
 * Computes total B field from both conductor current and magnetization.
 * B_total = B_conductor + B_magnet
 *
 * Combines conductor Biot-Savart and magnet B_genComp contributions.
 *
 * @param point        Evaluation point
 * @param segments     Conductor segments (may be nullptr)
 * @param n_segments   Number of conductor segments
 * @param current      Conductor current [A]
 * @param mag_centers  Magnet element centers (may be nullptr)
 * @param mag_volumes  Magnet element volumes
 * @param M_complex    Complex magnetization
 * @param n_mag_elems  Number of magnet elements
 * @param B_out        Output: total complex B field [Bx, By, Bz]
 */
void ComputeCombinedField(
    const TVector3d& point,
    const ConductorSegment* segments,
    int n_segments,
    std::complex<double> current,
    const double* mag_centers,
    const double* mag_volumes,
    const std::complex<double>* M_complex,
    int n_mag_elems,
    std::complex<double>* B_out
);

/**
 * @brief Batch combined conductor and magnetization field computation.
 *
 * @param points       Evaluation points
 * @param n_points     Number of points
 * @param segments     Conductor segments
 * @param n_segments   Number of conductor segments
 * @param current      Conductor current [A]
 * @param mag_centers  Magnet element centers
 * @param mag_volumes  Magnet element volumes
 * @param M_complex    Complex magnetization
 * @param n_mag_elems  Number of magnet elements
 * @param B_out        Output: total complex B field
 */
void ComputeCombinedFieldBatch(
    const double* points,
    int n_points,
    const ConductorSegment* segments,
    int n_segments,
    std::complex<double> current,
    const double* mag_centers,
    const double* mag_volumes,
    const std::complex<double>* M_complex,
    int n_mag_elems,
    std::complex<double>* B_out
);

/**
 * @brief Extract conductor segments from radTConductor
 *
 * Converts conductor wire path to segment array for unified field computation.
 *
 * @param wirePath     Wire centerline path (from radTConductor)
 * @param n_path_pts   Number of path points
 * @param segments     Output: segment array (caller must allocate n_path_pts-1)
 * @return             Number of segments created
 */
int ExtractConductorSegments(
    const TVector3d* wirePath,
    int n_path_pts,
    ConductorSegment* segments
);

} // namespace RadFieldUnified

#endif // RAD_FIELD_UNIFIED_H
