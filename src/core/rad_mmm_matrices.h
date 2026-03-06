/*
 * rad_mmm_matrices.h
 *
 * MMM (Magnetic Moment Method) standalone matrix construction
 *
 * This module provides Radia-object-independent matrix construction
 * for use with NGSolve mesh integration and external solvers.
 *
 * Part of Radia project - pybind11 module: mmm_core
 */

#ifndef RAD_MMM_MATRICES_H
#define RAD_MMM_MATRICES_H

#include <vector>
#include <array>
#include <complex>
#include <cmath>

namespace radia {
namespace mmm {

//=========================================================================
// 3D Vector (standalone, no gmvectf.h dependency)
//=========================================================================

struct Vec3d {
    double x, y, z;

    Vec3d() : x(0), y(0), z(0) {}
    Vec3d(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3d operator+(const Vec3d& v) const { return Vec3d(x+v.x, y+v.y, z+v.z); }
    Vec3d operator-(const Vec3d& v) const { return Vec3d(x-v.x, y-v.y, z-v.z); }
    Vec3d operator*(double s) const { return Vec3d(x*s, y*s, z*s); }
    Vec3d operator/(double s) const { return Vec3d(x/s, y/s, z/s); }

    double dot(const Vec3d& v) const { return x*v.x + y*v.y + z*v.z; }
    Vec3d cross(const Vec3d& v) const {
        return Vec3d(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
    }
    double norm() const { return std::sqrt(x*x + y*y + z*z); }
    Vec3d normalized() const {
        double n = norm();
        return (n > 1e-30) ? (*this / n) : Vec3d(0, 0, 0);
    }
};

//=========================================================================
// Triangle Face Data
//=========================================================================

struct TriangleFace {
    Vec3d v0, v1, v2;    // Vertices
    Vec3d normal;         // Outward normal
    double area;          // Face area

    // Pre-computed edge vectors for solid angle calculation
    Vec3d edge01, edge02;

    void ComputeDerivedQuantities();
};

//=========================================================================
// Material Types
//=========================================================================

/**
 * MMMMaterialType: Magnetic material behavior type
 *
 * Hierarchy of material models:
 *   1. FixedPM: Fixed remanent magnetization, no solve needed
 *   2. DemagPM: PM with demagnetization curve (future)
 *   3. LinearIsotropic: Constant chi (default soft material)
 *   4. NonlinearIsotropic: chi(|H|) from B-H curve
 *   5. LinearAnisotropic: chi tensor (future)
 *   6. NonlinearAnisotropic: chi tensor(|H|) (future)
 *   7. Preisach: Hysteresis model (future)
 */
enum class MMMMaterialType {
    LinearIsotropic,       // Default: constant chi
    NonlinearIsotropic,    // chi(|H|) from B-H curve
    FixedPM,               // Fixed remanent magnetization
    DemagPM,               // PM with demagnetization (future)
    LinearAnisotropic,     // chi tensor (future)
    NonlinearAnisotropic,  // chi tensor(|H|) (future)
    Preisach               // Hysteresis (future)
};

//=========================================================================
// MMM Element Types
//=========================================================================

/**
 * MMMElement: Magnetic element geometry data
 *
 * Supports:
 *   - Tetrahedron (4 faces, 3 DOF: Mx, My, Mz)
 *   - Hexahedron (6 faces, 6 DOF: sigma per face - MSC method)
 *   - Permanent magnet elements (fixed Mr, 0 DOF in solve)
 */
struct MMMElement {
    Vec3d center;           // Element centroid [m]
    double volume;          // Element volume [m^3]
    int dof;                // Degrees of freedom (3 or 6, 0 for fixed PM)

    // Face data (4 for tetra, 6 for hexa)
    // Each face stores triangles (1 for tetra faces, 2 for quad hexa faces)
    std::vector<TriangleFace> triangles;

    // For hexahedra: evaluation points per face (Yano-Sugahara midpoints)
    std::vector<Vec3d> eval_points;      // 6 points for hexa
    std::vector<Vec3d> face_normals;     // 6 normals for hexa
    std::vector<double> face_areas;      // 6 areas for hexa

    // Material type and properties
    MMMMaterialType material_type;  // Material behavior type
    Vec3d Mr;                       // Remanent magnetization [A/m] (for PM elements)

    // Convenience accessor
    bool IsPermanentMagnet() const {
        return material_type == MMMMaterialType::FixedPM ||
               material_type == MMMMaterialType::DemagPM;
    }

    MMMElement() : volume(0), dof(0), material_type(MMMMaterialType::LinearIsotropic) {}
};

//=========================================================================
// MMM Matrices Result
//=========================================================================

/**
 * MMMMatrices: Interaction matrix and metadata
 */
struct MMMMatrices {
    std::vector<double> N;          // Interaction matrix [totalDOF x totalDOF] (column-major for LAPACK)
    std::vector<int> dof_offset;    // DOF offset per element [n_elem + 1]
    int total_dof;                  // Total degrees of freedom
    int n_elements;                 // Number of elements

    // Material type per element
    std::vector<MMMMaterialType> material_types;  // Material type [n_elements]
    std::vector<double> Mr;         // Remanent magnetization [n_elements * 3]
    int n_pm_elements;              // Number of PM elements (FixedPM + DemagPM)
    int n_soft_elements;            // Number of soft magnetic elements

    // Convenience check
    bool IsPM(int elem_idx) const {
        return material_types[elem_idx] == MMMMaterialType::FixedPM ||
               material_types[elem_idx] == MMMMaterialType::DemagPM;
    }

    MMMMatrices() : total_dof(0), n_elements(0), n_pm_elements(0), n_soft_elements(0) {}
};

//=========================================================================
// MMMBuilder: Interaction Matrix Builder
//=========================================================================

/**
 * MMMBuilder: Constructs interaction matrix from element geometry
 *
 * Usage:
 *   MMMBuilder builder;
 *   builder.AddTetrahedron(vertices);  // 4 vertices
 *   builder.AddHexahedron(vertices);   // 8 vertices
 *   MMMMatrices result = builder.Build();
 */
class MMMBuilder {
public:
    MMMBuilder();
    ~MMMBuilder();

    /**
     * Add tetrahedron element (4 vertices, 3 DOF)
     * @param vertices: 4 vertices in order [v0, v1, v2, v3]
     *                  Face ordering follows standard tetrahedron topology
     */
    void AddTetrahedron(const std::vector<Vec3d>& vertices);

    /**
     * Add tetrahedron element from flat array
     * @param vertices: 12 doubles [x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3]
     */
    void AddTetrahedron(const double* vertices);

    /**
     * Add hexahedron element (8 vertices, 6 DOF MSC)
     * @param vertices: 8 vertices in standard hex order
     *                  [0-3: bottom face CCW, 4-7: top face CCW]
     */
    void AddHexahedron(const std::vector<Vec3d>& vertices);

    /**
     * Add hexahedron element from flat array
     * @param vertices: 24 doubles [x0,y0,z0, ..., x7,y7,z7]
     */
    void AddHexahedron(const double* vertices);

    /**
     * Add permanent magnet tetrahedron (fixed magnetization, 0 DOF in solve)
     * @param vertices: 4 vertices
     * @param Mr: Remanent magnetization [Mx, My, Mz] in A/m
     */
    void AddPermanentMagnetTetrahedron(const std::vector<Vec3d>& vertices, const Vec3d& Mr);
    void AddPermanentMagnetTetrahedron(const double* vertices, const double* Mr);

    /**
     * Add permanent magnet hexahedron (fixed magnetization, 0 DOF in solve)
     * @param vertices: 8 vertices
     * @param Mr: Remanent magnetization [Mx, My, Mz] in A/m
     */
    void AddPermanentMagnetHexahedron(const std::vector<Vec3d>& vertices, const Vec3d& Mr);
    void AddPermanentMagnetHexahedron(const double* vertices, const double* Mr);

    /**
     * Build interaction matrix
     * @return MMMMatrices containing N matrix and metadata
     */
    MMMMatrices Build();

    // Accessors
    int NumElements() const { return static_cast<int>(elements_.size()); }
    int TotalDOF() const;
    const std::vector<MMMElement>& GetElements() const { return elements_; }

    void Clear();

private:
    std::vector<MMMElement> elements_;

    // Block computation methods
    void Compute3x3Block(int elem_i, int elem_j, double* N_block) const;
    void Compute6x6Block(int elem_i, int elem_j, double* K_block) const;
    void Compute3x6Block(int tetra_idx, int hexa_idx, double* N_block) const;
    void Compute6x3Block(int hexa_idx, int tetra_idx, double* N_block) const;

    // Field computation helpers
    Vec3d FieldFromChargedTriangle(const TriangleFace& tri, const Vec3d& obs, double sigma) const;
    double SolidAngleTriangle(const Vec3d& obs, const Vec3d& v0, const Vec3d& v1, const Vec3d& v2) const;
};

//=========================================================================
// MMMSolver: Linear System Solver
//=========================================================================

/**
 * MMMSolver: Solves (I/chi - N) * M = H_ext
 */
class MMMSolver {
public:
    MMMSolver();
    ~MMMSolver();

    /**
     * Set interaction matrix
     */
    void SetMatrix(const MMMMatrices& matrices);

    /**
     * Set interaction matrix from external numpy arrays
     * @param N: Interaction matrix (column-major) [totalDOF x totalDOF]
     * @param dof_offset: DOF offset array [n_elem + 1]
     * @param total_dof: Total DOF count
     */
    void SetMatrix(const double* N, const int* dof_offset, int total_dof, int n_elem);

    /**
     * Solve using LU decomposition (LAPACK dgesv)
     * @param inv_chi: 1/chi values [totalDOF] or [n_elem] (broadcast to DOF)
     * @param H_ext: External field [totalDOF]
     * @param M: Output magnetization [totalDOF]
     * @param chi_per_element: If true, inv_chi is [n_elem] and broadcast
     * @return 0 on success, non-zero on failure
     */
    int SolveLU(const std::vector<double>& inv_chi,
                const std::vector<double>& H_ext,
                std::vector<double>& M,
                bool chi_per_element = false);

    /**
     * Solve using BiCGSTAB iteration
     * @return Number of iterations (-1 on failure)
     */
    int SolveBiCGSTAB(const std::vector<double>& inv_chi,
                       const std::vector<double>& H_ext,
                       std::vector<double>& M,
                       double tol = 1e-8,
                       int max_iter = 1000,
                       bool chi_per_element = false);

    // Accessors
    int TotalDOF() const { return total_dof_; }
    int NumElements() const { return n_elements_; }
    bool IsReady() const { return matrix_set_; }

private:
    std::vector<double> N_;          // Interaction matrix (column-major)
    std::vector<int> dof_offset_;    // DOF offsets
    int total_dof_;
    int n_elements_;
    bool matrix_set_;

    // Build system matrix A = diag(inv_chi) - N
    void BuildSystemMatrix(std::vector<double>& A,
                           const std::vector<double>& inv_chi,
                           bool chi_per_element) const;

    // Matrix-vector product: y = A * x
    void MatVec(const std::vector<double>& A,
                const std::vector<double>& x,
                std::vector<double>& y) const;
};

//=========================================================================
// MMMFieldComputer: Field Computation from Magnetization
//=========================================================================

/**
 * MMMFieldComputer: Computes B/H fields from solved magnetization
 *
 * Stores element geometry and computes field at observation points.
 */
class MMMFieldComputer {
public:
    MMMFieldComputer();
    ~MMMFieldComputer();

    /**
     * Set element geometry from MMMBuilder
     * @param elements: Element geometry data
     */
    void SetElements(const std::vector<MMMElement>& elements);

    /**
     * Set material data for field computation
     * @param material_types: Material type per element [n_elements]
     * @param Mr: Remanent magnetization [n_elements * 3]
     */
    void SetMaterialData(const std::vector<MMMMaterialType>& material_types,
                         const std::vector<double>& Mr);

    /**
     * Compute B field at observation points
     * @param M: Solved magnetization [total_dof] (for soft magnetic elements)
     * @param obs_points: Observation points [n_points x 3]
     * @param n_points: Number of observation points
     * @return B field [n_points x 3] in Tesla
     *
     * Note: For PM elements, field is computed from fixed Mr.
     *       For soft elements, field is computed from solved M.
     */
    std::vector<double> ComputeBField(
        const std::vector<double>& M,
        const std::vector<double>& obs_points,
        int n_points) const;

    /**
     * Compute H field at observation points
     * @param M: Solved magnetization [total_dof]
     * @param obs_points: Observation points [n_points x 3]
     * @param n_points: Number of observation points
     * @return H field [n_points x 3] in A/m
     */
    std::vector<double> ComputeHField(
        const std::vector<double>& M,
        const std::vector<double>& obs_points,
        int n_points) const;

    /**
     * Get element centers for diagnostics
     * @return Element centers [n_elements x 3]
     */
    std::vector<double> GetElementCenters() const;

    /**
     * Get element volumes
     * @return Element volumes [n_elements]
     */
    std::vector<double> GetElementVolumes() const;

    int NumElements() const { return static_cast<int>(elements_.size()); }
    int TotalDOF() const;

private:
    std::vector<MMMElement> elements_;
    std::vector<int> dof_offset_;

    // Material data
    std::vector<MMMMaterialType> material_types_;
    std::vector<double> Mr_;  // [n_elements * 3]

    // Field from dipole at element center
    Vec3d DipoleField(const Vec3d& M, const Vec3d& elem_center, double volume,
                      const Vec3d& obs) const;

    // Field from surface charges (MSC)
    Vec3d SurfaceChargeField(const MMMElement& elem, const double* sigma,
                              const Vec3d& obs) const;
};

//=========================================================================
// Nonlinear Iteration Helpers
//=========================================================================

/**
 * Compute chi from B-H curve using linear interpolation
 *
 * @param M_flat: Magnetization [n_elem * 3] (3 DOF per element)
 * @param H_flat: Local H field [n_elem * 3]
 * @param bh_curve: B-H curve [[H0,B0], [H1,B1], ...] sorted by H ascending
 * @param n_elem: Number of elements
 * @return chi values [n_elem]
 */
std::vector<double> ComputeChiFromBH(
    const std::vector<double>& M_flat,
    const std::vector<double>& H_flat,
    const std::vector<std::array<double, 2>>& bh_curve,
    int n_elem);

/**
 * Check convergence based on B-field change
 *
 * @param B_new: New B field [n_elem * 3]
 * @param B_old: Old B field [n_elem * 3]
 * @param B_sat: Saturation B for relative error
 * @param n_elem: Number of elements
 * @return Maximum relative change
 */
double CheckConvergence(
    const std::vector<double>& B_new,
    const std::vector<double>& B_old,
    double B_sat,
    int n_elem);

//=========================================================================
// MMMNonlinearSolver: Nonlinear Material Solver
//=========================================================================

/**
 * Result structure for nonlinear solve
 */
struct MMMNonlinearResult {
    std::vector<double> M;       // Final magnetization [total_dof]
    std::vector<double> chi;     // Final chi per element [n_elem]
    int iterations;              // Number of iterations used
    double final_error;          // Final convergence error
    bool converged;              // True if converged
};

/**
 * MMMNonlinearSolver: Newton-Raphson iteration for nonlinear materials
 *
 * Solves the nonlinear system where chi = chi(|H|) depends on the local H field.
 * Uses fixed-point iteration with optional relaxation.
 *
 * Algorithm:
 *   1. Initialize chi with initial guess (default: linear mu_r)
 *   2. Solve linear system: (I/chi - N) * M = H_ext
 *   3. Compute local H = H_ext + N * M
 *   4. Update chi from B-H curve: chi = chi(|H|)
 *   5. Check convergence
 *   6. Repeat until converged or max iterations
 */
class MMMNonlinearSolver {
public:
    MMMNonlinearSolver();
    ~MMMNonlinearSolver();

    /**
     * Set interaction matrix
     */
    void SetMatrix(const double* N, const int* dof_offset, int total_dof, int n_elem);

    /**
     * Set B-H curve for nonlinear material
     * @param bh_curve: [[H0,B0], [H1,B1], ...] sorted by H ascending
     *                  H in A/m, B in Tesla
     */
    void SetBHCurve(const std::vector<std::array<double, 2>>& bh_curve);

    /**
     * Set solver parameters
     * @param tol: Convergence tolerance (relative B change)
     * @param max_iter: Maximum iterations
     * @param relax: Relaxation factor (0 < relax <= 1, 1=no relaxation)
     * @param linear_solver: 0=LU, 1=BiCGSTAB
     */
    void SetParameters(double tol = 1e-4, int max_iter = 100,
                       double relax = 1.0, int linear_solver = 0);

    /**
     * Set initial guess for chi (optional)
     * @param chi_init: Initial chi per element [n_elem]
     *                  If not set, uses chi from first point on B-H curve
     */
    void SetInitialChi(const std::vector<double>& chi_init);

    /**
     * Solve nonlinear system
     * @param H_ext: External field (3-component vector for uniform field,
     *               or [total_dof] for non-uniform)
     * @return MMMNonlinearResult with solution and convergence info
     */
    MMMNonlinearResult Solve(const std::vector<double>& H_ext);

    /**
     * Get local H field at element centers after solve
     * @return H_local [n_elem * 3]
     */
    std::vector<double> GetLocalH() const { return H_local_; }

    /**
     * Get B field at element centers after solve
     * @return B_local [n_elem * 3]
     */
    std::vector<double> GetLocalB() const;

    // Accessors
    int TotalDOF() const { return total_dof_; }
    int NumElements() const { return n_elements_; }

private:
    // Matrix data
    std::vector<double> N_;
    std::vector<int> dof_offset_;
    int total_dof_;
    int n_elements_;

    // B-H curve
    std::vector<std::array<double, 2>> bh_curve_;

    // Solver parameters
    double tol_;
    int max_iter_;
    double relax_;
    int linear_solver_;

    // Initial chi
    std::vector<double> chi_init_;
    bool has_chi_init_;

    // State after solve
    std::vector<double> H_local_;  // Local H field [n_elem * 3]
    std::vector<double> M_;        // Magnetization [total_dof]
    std::vector<double> chi_;      // Chi per element [n_elem]

    // Internal methods
    double InterpolateBH(double H_mag) const;
    double ComputeChiFromH(double H_mag) const;
    void ComputeLocalH(const std::vector<double>& H_ext,
                       const std::vector<double>& M);
    std::vector<double> ExpandH_ext(const std::vector<double>& H_ext) const;
};

//=========================================================================
// Physical Constants
//=========================================================================

namespace constants {
    constexpr double MU_0 = 4.0 * 3.14159265358979323846 * 1e-7;  // H/m
    constexpr double INV_FOUR_PI = 1.0 / (4.0 * 3.14159265358979323846);
    constexpr double MU_0_OVER_FOUR_PI = 1e-7;  // H/m
}

}  // namespace mmm
}  // namespace radia

#endif  // RAD_MMM_MATRICES_H
