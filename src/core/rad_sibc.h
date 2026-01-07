/*
 * rad_sibc.h
 *
 * Nonlocal Surface Impedance Boundary Condition (SIBC) for conductive magnetic materials
 *
 * Based on:
 * S. Bilicz, Z. Badics, J. Pávó, "Wide-band nonlocal impedance boundary condition
 * model for high-conductivity regions in integral equation framework",
 * ISEM 2023 (International Symposium on Electromagnetic Fields in Mechatronics,
 * Electrical and Electronic Engineering), December 2023.
 *
 * Key concept: For materials with both conductivity (σ) and permeability (μr),
 * the standard local SIBC (Et = Zs * K) is inaccurate at low frequencies.
 * The nonlocal SIBC uses 2D cross-section FEM to capture the full frequency-dependent
 * behavior from DC to resonance.
 *
 * Architecture:
 * 1. Extract 2D cross-section of conductor
 * 2. Solve 2D Helmholtz equation for E-field in cross-section
 * 3. Use boundary values to construct nonlocal impedance operator Z{.}
 * 4. Apply SIBC: Et = Z{n × K} in surface integral equation
 *
 * Part of Radia project
 */

#ifndef RAD_SIBC_H
#define RAD_SIBC_H

#include "rad_geom_types.h"
#include <vector>
#include <complex>
#include <memory>
#include <functional>

namespace radia {

// Complex type alias
using Complex = std::complex<double>;

/**
 * @brief 2D cross-section mesh for SIBC calculation
 *
 * Represents the interior cross-section of a conductor.
 * Used for solving the 2D Helmholtz equation.
 */
struct CrossSection2D {
    // Node coordinates (x, y in local 2D system)
    std::vector<double> nodeX;
    std::vector<double> nodeY;

    // Triangle elements (3 node indices per triangle)
    std::vector<std::array<int, 3>> triangles;

    // Boundary nodes (ordered around perimeter)
    std::vector<int> boundaryNodes;

    // Cross-section type
    enum class Type {
        Rectangular,
        Circular,
        Custom
    } type;

    // Dimensions (for standard shapes)
    double width;      // For rectangular
    double height;     // For rectangular
    double radius;     // For circular

    // Material properties
    double conductivity;   // σ [S/m]
    double permeability;   // μ = μ0 * μr [H/m]

    CrossSection2D()
        : type(Type::Rectangular)
        , width(0)
        , height(0)
        , radius(0)
        , conductivity(1e7)
        , permeability(4e-7 * 3.14159265358979)  // μ0
    {}
};

/**
 * @brief 2D FEM solver for Helmholtz equation in cross-section
 *
 * Solves: ∇²E - jωμσE = 0 in Ω (cross-section interior)
 * with BC: ∂E/∂n = -jωμHv on Γ (boundary)
 *
 * where Hv is the tangential H-field at the boundary.
 */
class radTCrossSection2DFEM {
public:
    radTCrossSection2DFEM();
    ~radTCrossSection2DFEM();

    /**
     * @brief Set cross-section geometry
     */
    void SetCrossSection(const CrossSection2D& cs);

    /**
     * @brief Generate mesh for rectangular cross-section
     * @param width Rectangle width
     * @param height Rectangle height
     * @param nX Elements in x direction
     * @param nY Elements in y direction
     */
    void MeshRectangle(double width, double height, int nX = 10, int nY = 10);

    /**
     * @brief Generate mesh for circular cross-section
     * @param radius Circle radius
     * @param nRadial Radial divisions
     * @param nAngular Angular divisions
     */
    void MeshCircle(double radius, int nRadial = 8, int nAngular = 24);

    /**
     * @brief Set operating frequency
     */
    void SetFrequency(double frequency);

    /**
     * @brief Set material properties
     * @param conductivity σ [S/m]
     * @param relativePermeability μr
     */
    void SetMaterial(double conductivity, double relativePermeability);

    /**
     * @brief Assemble FEM system matrix
     *
     * Builds: [K - jωμσM] where
     * K is the stiffness matrix from ∇²
     * M is the mass matrix
     */
    void AssembleSystem();

    /**
     * @brief Solve for given boundary H-field
     *
     * @param Hv_boundary Tangential H-field at boundary nodes [A/m]
     * @param E_interior Output E-field at all nodes [V/m]
     */
    void Solve(const std::vector<Complex>& Hv_boundary,
               std::vector<Complex>& E_interior);

    /**
     * @brief Get E-field at boundary (for SIBC)
     */
    void GetBoundaryE(std::vector<Complex>& E_boundary) const;

    /**
     * @brief Compute impedance matrix Z{.}
     *
     * Maps boundary H-field to boundary E-field: E_boundary = Z * H_boundary
     *
     * @param Z_matrix Output impedance matrix (n_boundary × n_boundary)
     */
    void ComputeImpedanceMatrix(std::vector<Complex>& Z_matrix);

    /**
     * @brief Get skin depth at current frequency
     */
    double GetSkinDepth() const { return skinDepth_; }

    /**
     * @brief Get number of boundary nodes
     */
    int NumBoundaryNodes() const;

    /**
     * @brief Get number of interior nodes
     */
    int NumInteriorNodes() const;

private:
    CrossSection2D cs_;
    double frequency_;
    double omega_;
    double skinDepth_;

    // FEM matrices
    std::vector<Complex> stiffnessMatrix_;  // K (sparse, stored as dense for small problems)
    std::vector<Complex> massMatrix_;        // M
    std::vector<Complex> systemMatrix_;      // K - jωμσM

    // Factorized system (for efficient multiple solves)
    std::vector<Complex> systemLU_;
    std::vector<int> pivots_;
    bool factorized_;

    // Solution storage
    std::vector<Complex> solution_;

    // Interior/boundary node mapping
    std::vector<int> interiorNodes_;
    std::vector<int> boundaryNodeMap_;  // Global to boundary index

    void AssembleElementMatrix(int elemIdx,
                                std::vector<Complex>& Ke,
                                std::vector<Complex>& Me);

    void ApplyBoundaryConditions();
    void FactorizeSystem();
};

/**
 * @brief Nonlocal SIBC operator for surface integral equation
 *
 * Implements the nonlocal surface impedance boundary condition:
 * Et = Z{n × K}
 *
 * where Z{.} is a nonlocal operator computed from 2D cross-section FEM.
 */
class radTNonlocalSIBC {
public:
    radTNonlocalSIBC();
    ~radTNonlocalSIBC();

    /**
     * @brief Set wire/conductor parameters
     * @param crossSection Cross-section type ("rectangular" or "circular")
     * @param width Width or diameter
     * @param height Height (for rectangular)
     */
    void SetConductorGeometry(const std::string& crossSection,
                               double width, double height = 0);

    /**
     * @brief Set material properties
     * @param conductivity σ [S/m]
     * @param relativePermeability μr
     */
    void SetMaterial(double conductivity, double relativePermeability);

    /**
     * @brief Set operating frequency
     */
    void SetFrequency(double frequency);

    /**
     * @brief Compute SIBC impedance matrix
     *
     * This matrix relates surface current K to tangential E-field:
     * Et = Z_sibc * K
     *
     * For nonlocal SIBC, this is NOT a diagonal matrix.
     *
     * @param n_panels Number of surface panels along wire length
     * @param Z_sibc Output impedance matrix (n_panels × n_panels)
     */
    void ComputeImpedanceMatrix(int n_panels, std::vector<Complex>& Z_sibc);

    /**
     * @brief Apply SIBC to surface current vector
     *
     * Computes Et = Z{K} using precomputed impedance matrix.
     *
     * @param K Surface current at panels [A/m]
     * @param Et Output tangential E-field [V/m]
     */
    void Apply(const std::vector<Complex>& K, std::vector<Complex>& Et);

    /**
     * @brief Get local surface impedance (for comparison)
     *
     * Zs = (1+j) / (σδ) where δ = sqrt(2/(ωμσ)) is skin depth
     */
    Complex GetLocalSurfaceImpedance() const;

    /**
     * @brief Check if local SIBC is sufficient
     *
     * Returns true if skin depth << conductor dimension,
     * meaning local SIBC is accurate and nonlocal is not needed.
     */
    bool IsLocalSIBCSufficient() const;

    /**
     * @brief Get internal resistance per unit length [Ω/m]
     */
    double GetDCResistance() const;

    /**
     * @brief Get internal inductance per unit length [H/m]
     */
    double GetInternalInductance() const;

private:
    // Cross-section parameters
    std::string crossSectionType_;
    double width_;
    double height_;

    // Material
    double conductivity_;
    double relPermeability_;

    // Frequency
    double frequency_;
    double omega_;
    double skinDepth_;

    // 2D FEM solver
    std::unique_ptr<radTCrossSection2DFEM> fem2D_;

    // Precomputed impedance matrix
    std::vector<Complex> Z_sibc_;
    int n_panels_;
    bool impedanceComputed_;

    void UpdateSkinDepth();
};

/**
 * @brief SIBC surface element for conductive magnetic material
 *
 * Represents a surface panel of a conductive magnetic object.
 * Used in the surface integral equation with SIBC.
 */
class radTSIBCElement {
public:
    radTSIBCElement();
    ~radTSIBCElement();

    /**
     * @brief Set element geometry
     * @param center Panel center
     * @param normal Outward normal
     * @param area Panel area
     * @param tangent1 First tangent direction
     * @param tangent2 Second tangent direction
     */
    void SetGeometry(const TVector3d& center,
                     const TVector3d& normal,
                     double area,
                     const TVector3d& tangent1,
                     const TVector3d& tangent2);

    /**
     * @brief Set SIBC operator
     */
    void SetSIBC(std::shared_ptr<radTNonlocalSIBC> sibc);

    /**
     * @brief Get surface current density K
     */
    const std::vector<Complex>& GetSurfaceCurrent() const { return K_; }
    std::vector<Complex>& GetSurfaceCurrent() { return K_; }

    /**
     * @brief Get tangential E-field from SIBC
     */
    void ComputeTangentialE(std::vector<Complex>& Et) const;

    // Geometry accessors
    const TVector3d& Center() const { return center_; }
    const TVector3d& Normal() const { return normal_; }
    double Area() const { return area_; }

private:
    TVector3d center_;
    TVector3d normal_;
    TVector3d tangent1_;
    TVector3d tangent2_;
    double area_;

    std::shared_ptr<radTNonlocalSIBC> sibc_;
    std::vector<Complex> K_;   // Surface current density [A/m]
};

/**
 * @brief SIBC solver for conductive magnetic materials
 *
 * Solves surface integral equation with nonlocal SIBC:
 *
 * EFIE with SIBC:
 * n × (E_inc + E_scat) = n × Z{n × K}
 *
 * where Z{.} is the nonlocal impedance operator from 2D FEM.
 */
class radTSIBCSolver {
public:
    radTSIBCSolver();
    ~radTSIBCSolver();

    /**
     * @brief Add SIBC element
     */
    void AddElement(std::shared_ptr<radTSIBCElement> element);

    /**
     * @brief Set operating frequency
     */
    void SetFrequency(double frequency);

    /**
     * @brief Set material properties
     */
    void SetMaterial(double conductivity, double relPermeability);

    /**
     * @brief Solve for surface currents
     */
    void Solve();

    /**
     * @brief Compute impedance at port
     */
    Complex ComputeImpedance();

    /**
     * @brief Frequency sweep
     */
    std::vector<Complex> ImpedanceSweep(const std::vector<double>& frequencies);

    /**
     * @brief Enable comparison with local SIBC
     *
     * When enabled, also computes solution with local SIBC for validation.
     */
    void EnableLocalComparison(bool enable) { compareWithLocal_ = enable; }

    /**
     * @brief Get relative error between local and nonlocal SIBC
     *
     * Only valid after Solve() with EnableLocalComparison(true).
     */
    double GetLocalVsNonlocalError() const { return localVsNonlocalError_; }

private:
    std::vector<std::shared_ptr<radTSIBCElement>> elements_;
    double frequency_;
    double conductivity_;
    double relPermeability_;

    // SIBC operator
    std::shared_ptr<radTNonlocalSIBC> sibc_;

    // System matrix and solution
    std::vector<Complex> systemMatrix_;
    std::vector<Complex> rhs_;
    std::vector<Complex> solution_;

    // Comparison
    bool compareWithLocal_;
    double localVsNonlocalError_;

    void BuildSystemMatrix();
    void SolveLinearSystem();
    void ExtractSolution();
};

} // namespace radia

#endif // RAD_SIBC_H
