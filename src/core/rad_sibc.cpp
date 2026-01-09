/*
 * rad_sibc.cpp
 *
 * Nonlocal Surface Impedance Boundary Condition (SIBC) implementation
 *
 * Based on:
 * S. Bilicz, Z. Badics, J. Pávó, "Wide-band nonlocal impedance boundary condition
 * model for high-conductivity regions in integral equation framework",
 * ISEM 2023.
 *
 * Part of Radia project
 */

#include "rad_sibc.h"
#include "rad_constants.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace radia {

// ============================================================================
// radTCrossSection2DFEM implementation
// ============================================================================

radTCrossSection2DFEM::radTCrossSection2DFEM()
    : frequency_(0)
    , omega_(0)
    , skinDepth_(1e30)
    , factorized_(false)
{
}

radTCrossSection2DFEM::~radTCrossSection2DFEM() {
}

void radTCrossSection2DFEM::SetCrossSection(const CrossSection2D& cs) {
    cs_ = cs;
    factorized_ = false;
}

void radTCrossSection2DFEM::MeshRectangle(double width, double height, int nX, int nY) {
    cs_.type = CrossSection2D::Type::Rectangular;
    cs_.width = width;
    cs_.height = height;

    cs_.nodeX.clear();
    cs_.nodeY.clear();
    cs_.triangles.clear();
    cs_.boundaryNodes.clear();

    // Generate nodes
    double dx = width / nX;
    double dy = height / nY;

    for (int j = 0; j <= nY; ++j) {
        for (int i = 0; i <= nX; ++i) {
            cs_.nodeX.push_back(-width/2 + i * dx);
            cs_.nodeY.push_back(-height/2 + j * dy);
        }
    }

    // Generate triangles (2 per rectangle)
    for (int j = 0; j < nY; ++j) {
        for (int i = 0; i < nX; ++i) {
            int n00 = j * (nX + 1) + i;
            int n10 = n00 + 1;
            int n01 = n00 + (nX + 1);
            int n11 = n01 + 1;

            // Lower-left triangle
            cs_.triangles.push_back(std::array<int, 3>{{n00, n10, n01}});
            // Upper-right triangle
            cs_.triangles.push_back(std::array<int, 3>{{n10, n11, n01}});
        }
    }

    // Identify boundary nodes (ordered counter-clockwise)
    // Bottom edge
    for (int i = 0; i <= nX; ++i) {
        cs_.boundaryNodes.push_back(i);
    }
    // Right edge (skip first, already added)
    for (int j = 1; j <= nY; ++j) {
        cs_.boundaryNodes.push_back(j * (nX + 1) + nX);
    }
    // Top edge (reverse, skip first)
    for (int i = nX - 1; i >= 0; --i) {
        cs_.boundaryNodes.push_back(nY * (nX + 1) + i);
    }
    // Left edge (reverse, skip first and last)
    for (int j = nY - 1; j >= 1; --j) {
        cs_.boundaryNodes.push_back(j * (nX + 1));
    }

    // Build boundary node map
    int nNodes = static_cast<int>(cs_.nodeX.size());
    boundaryNodeMap_.resize(nNodes, -1);
    for (size_t i = 0; i < cs_.boundaryNodes.size(); ++i) {
        boundaryNodeMap_[cs_.boundaryNodes[i]] = static_cast<int>(i);
    }

    // Interior nodes
    interiorNodes_.clear();
    for (int n = 0; n < nNodes; ++n) {
        if (boundaryNodeMap_[n] < 0) {
            interiorNodes_.push_back(n);
        }
    }

    factorized_ = false;
}

void radTCrossSection2DFEM::MeshCircle(double radius, int nRadial, int nAngular) {
    cs_.type = CrossSection2D::Type::Circular;
    cs_.radius = radius;

    cs_.nodeX.clear();
    cs_.nodeY.clear();
    cs_.triangles.clear();
    cs_.boundaryNodes.clear();

    // Center node
    cs_.nodeX.push_back(0);
    cs_.nodeY.push_back(0);

    // Radial rings
    for (int ir = 1; ir <= nRadial; ++ir) {
        double r = radius * ir / nRadial;
        for (int ia = 0; ia < nAngular; ++ia) {
            double theta = 2.0 * RadConst::PI * ia / nAngular;
            cs_.nodeX.push_back(r * std::cos(theta));
            cs_.nodeY.push_back(r * std::sin(theta));
        }
    }

    // Triangles: center to first ring
    for (int ia = 0; ia < nAngular; ++ia) {
        int n0 = 0;  // Center
        int n1 = 1 + ia;
        int n2 = 1 + (ia + 1) % nAngular;
        cs_.triangles.push_back(std::array<int, 3>{{n0, n1, n2}});
    }

    // Triangles: ring to ring
    for (int ir = 1; ir < nRadial; ++ir) {
        int ringStart1 = 1 + (ir - 1) * nAngular;
        int ringStart2 = 1 + ir * nAngular;

        for (int ia = 0; ia < nAngular; ++ia) {
            int n00 = ringStart1 + ia;
            int n01 = ringStart1 + (ia + 1) % nAngular;
            int n10 = ringStart2 + ia;
            int n11 = ringStart2 + (ia + 1) % nAngular;

            cs_.triangles.push_back(std::array<int, 3>{{n00, n10, n01}});
            cs_.triangles.push_back(std::array<int, 3>{{n10, n11, n01}});
        }
    }

    // Boundary nodes (outer ring)
    int outerRingStart = 1 + (nRadial - 1) * nAngular;
    for (int ia = 0; ia < nAngular; ++ia) {
        cs_.boundaryNodes.push_back(outerRingStart + ia);
    }

    // Build boundary node map
    int nNodes = static_cast<int>(cs_.nodeX.size());
    boundaryNodeMap_.resize(nNodes, -1);
    for (size_t i = 0; i < cs_.boundaryNodes.size(); ++i) {
        boundaryNodeMap_[cs_.boundaryNodes[i]] = static_cast<int>(i);
    }

    // Interior nodes
    interiorNodes_.clear();
    for (int n = 0; n < nNodes; ++n) {
        if (boundaryNodeMap_[n] < 0) {
            interiorNodes_.push_back(n);
        }
    }

    factorized_ = false;
}

void radTCrossSection2DFEM::SetFrequency(double frequency) {
    frequency_ = frequency;
    omega_ = 2.0 * RadConst::PI * frequency;

    // Update skin depth using |μ| for estimation
    if (cs_.conductivity > 0 && omega_ > 0) {
        double mu_abs = cs_.GetAbsolutePermeability();
        skinDepth_ = std::sqrt(2.0 / (omega_ * mu_abs * cs_.conductivity));
    } else {
        skinDepth_ = 1e30;
    }

    factorized_ = false;
}

void radTCrossSection2DFEM::SetMaterial(double conductivity, double relativePermeability) {
    cs_.conductivity = conductivity;
    cs_.permeability_real = MU_0 * relativePermeability;
    cs_.permeability_imag = 0.0;  // No magnetic loss
    factorized_ = false;
}

void radTCrossSection2DFEM::SetComplexMaterial(double conductivity,
                                                double mu_prime_r,
                                                double mu_double_prime_r) {
    cs_.conductivity = conductivity;
    cs_.permeability_real = MU_0 * mu_prime_r;
    cs_.permeability_imag = MU_0 * mu_double_prime_r;
    factorized_ = false;
}

void radTCrossSection2DFEM::AssembleSystem() {
    int nNodes = static_cast<int>(cs_.nodeX.size());
    int nElem = static_cast<int>(cs_.triangles.size());

    // Initialize matrices
    stiffnessMatrix_.assign(nNodes * nNodes, Complex(0, 0));
    massMatrix_.assign(nNodes * nNodes, Complex(0, 0));

    // Element assembly
    std::vector<Complex> Ke(9), Me(9);

    for (int e = 0; e < nElem; ++e) {
        AssembleElementMatrix(e, Ke, Me);

        const auto& tri = cs_.triangles[e];

        // Add to global matrices
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                int gi = tri[i];
                int gj = tri[j];
                stiffnessMatrix_[gi * nNodes + gj] += Ke[i * 3 + j];
                massMatrix_[gi * nNodes + gj] += Me[i * 3 + j];
            }
        }
    }

    // Form system matrix: K - jωμσM
    // For complex permeability μ = μ' - jμ":
    //   factor = -jω(μ' - jμ")σ = -jωμ'σ - ωμ"σ
    //          = ωμ"σ - jωμ'σ
    // This shows that:
    //   - Real part (ωμ"σ): Additional dissipation from magnetic loss
    //   - Imaginary part (-ωμ'σ): Reactive response from energy storage
    Complex mu_complex = cs_.GetComplexPermeability();  // μ' - jμ"
    Complex factor = Complex(0, -omega_) * mu_complex * cs_.conductivity;

    systemMatrix_.resize(nNodes * nNodes);

    for (int i = 0; i < nNodes * nNodes; ++i) {
        systemMatrix_[i] = stiffnessMatrix_[i] + factor * massMatrix_[i];
    }

    factorized_ = false;
}

void radTCrossSection2DFEM::AssembleElementMatrix(int elemIdx,
                                                   std::vector<Complex>& Ke,
                                                   std::vector<Complex>& Me) {
    const auto& tri = cs_.triangles[elemIdx];

    // Get node coordinates
    double x[3], y[3];
    for (int i = 0; i < 3; ++i) {
        x[i] = cs_.nodeX[tri[i]];
        y[i] = cs_.nodeY[tri[i]];
    }

    // Triangle area
    double area = 0.5 * std::abs((x[1] - x[0]) * (y[2] - y[0]) -
                                  (x[2] - x[0]) * (y[1] - y[0]));

    // Shape function gradients (constant over triangle)
    double b[3], c[3];  // grad N_i = (b_i, c_i) / (2*area)
    b[0] = y[1] - y[2];  b[1] = y[2] - y[0];  b[2] = y[0] - y[1];
    c[0] = x[2] - x[1];  c[1] = x[0] - x[2];  c[2] = x[1] - x[0];

    // Stiffness matrix: K_ij = integral{ grad(N_i) . grad(N_j) } dA
    // For linear triangles: K_ij = (b_i*b_j + c_i*c_j) / (4*area)
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Ke[i * 3 + j] = Complex((b[i] * b[j] + c[i] * c[j]) / (4.0 * area), 0);
        }
    }

    // Mass matrix: M_ij = integral{ N_i * N_j } dA
    // For linear triangles: M_ij = area * (1 + delta_ij) / 12
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double mij = (i == j) ? (area / 6.0) : (area / 12.0);
            Me[i * 3 + j] = Complex(mij, 0);
        }
    }
}

void radTCrossSection2DFEM::Solve(const std::vector<Complex>& Hv_boundary,
                                   std::vector<Complex>& E_interior) {
    if (!factorized_) {
        FactorizeSystem();
    }

    int nNodes = static_cast<int>(cs_.nodeX.size());
    int nBoundary = static_cast<int>(cs_.boundaryNodes.size());

    // Build RHS from boundary condition: ∂E/∂n = -jωμHv
    // For complex μ = μ' - jμ":
    //   bcFactor = -jω(μ' - jμ") = -jωμ' - ωμ"
    //            = -ωμ" - jωμ'
    std::vector<Complex> rhs(nNodes, Complex(0, 0));
    Complex mu_complex = cs_.GetComplexPermeability();
    Complex bcFactor = Complex(0, -omega_) * mu_complex;

    // Apply Neumann BC as surface integral (simplified for boundary nodes)
    for (int ib = 0; ib < nBoundary; ++ib) {
        int node = cs_.boundaryNodes[ib];
        rhs[node] = bcFactor * Hv_boundary[ib];
    }

    // Solve system
    E_interior = rhs;

    // Back-substitution with pre-factorized LU
    int n = nNodes;
    #ifdef USE_MKL
    // Use MKL LAPACK zgetrs
    LAPACKE_zgetrs(LAPACK_COL_MAJOR, 'N', n, 1,
                   reinterpret_cast<const lapack_complex_double*>(systemLU_.data()),
                   n, pivots_.data(),
                   reinterpret_cast<lapack_complex_double*>(E_interior.data()), n);
    #else
    // Fallback: forward/back substitution with LU
    // Forward substitution (L * y = P * b)
    for (int i = 0; i < n; ++i) {
        // Apply permutation
        if (pivots_[i] != i) {
            std::swap(E_interior[i], E_interior[pivots_[i]]);
        }
        // Forward solve
        for (int j = 0; j < i; ++j) {
            E_interior[i] -= systemLU_[j * n + i] * E_interior[j];
        }
    }
    // Back substitution (U * x = y)
    for (int i = n - 1; i >= 0; --i) {
        for (int j = i + 1; j < n; ++j) {
            E_interior[i] -= systemLU_[j * n + i] * E_interior[j];
        }
        E_interior[i] /= systemLU_[i * n + i];
    }
    #endif

    solution_ = E_interior;
}

void radTCrossSection2DFEM::GetBoundaryE(std::vector<Complex>& E_boundary) const {
    int nBoundary = static_cast<int>(cs_.boundaryNodes.size());
    E_boundary.resize(nBoundary);

    for (int ib = 0; ib < nBoundary; ++ib) {
        int node = cs_.boundaryNodes[ib];
        E_boundary[ib] = solution_[node];
    }
}

void radTCrossSection2DFEM::ComputeImpedanceMatrix(std::vector<Complex>& Z_matrix) {
    // Z maps boundary H to boundary E: E_boundary = Z * H_boundary
    // We solve for each unit H-field at boundary nodes

    int nBoundary = NumBoundaryNodes();
    Z_matrix.resize(nBoundary * nBoundary);

    std::vector<Complex> Hv_unit(nBoundary, Complex(0, 0));
    std::vector<Complex> E_sol;
    std::vector<Complex> E_boundary;

    for (int j = 0; j < nBoundary; ++j) {
        // Set unit H at boundary node j
        std::fill(Hv_unit.begin(), Hv_unit.end(), Complex(0, 0));
        Hv_unit[j] = Complex(1, 0);

        // Solve
        Solve(Hv_unit, E_sol);
        GetBoundaryE(E_boundary);

        // Column j of Z_matrix
        for (int i = 0; i < nBoundary; ++i) {
            Z_matrix[j * nBoundary + i] = E_boundary[i];
        }
    }
}

int radTCrossSection2DFEM::NumBoundaryNodes() const {
    return static_cast<int>(cs_.boundaryNodes.size());
}

int radTCrossSection2DFEM::NumInteriorNodes() const {
    return static_cast<int>(interiorNodes_.size());
}

void radTCrossSection2DFEM::FactorizeSystem() {
    int n = static_cast<int>(cs_.nodeX.size());
    systemLU_ = systemMatrix_;
    pivots_.resize(n);

    #ifdef USE_MKL
    lapack_int info = LAPACKE_zgetrf(LAPACK_COL_MAJOR, n, n,
                                      reinterpret_cast<lapack_complex_double*>(systemLU_.data()),
                                      n, pivots_.data());
    if (info != 0) {
        throw std::runtime_error("SIBC FEM: LU factorization failed");
    }
    #else
    // Fallback: simple LU with partial pivoting
    for (int k = 0; k < n; ++k) {
        // Find pivot
        int maxRow = k;
        double maxVal = std::abs(systemLU_[k * n + k]);
        for (int i = k + 1; i < n; ++i) {
            double val = std::abs(systemLU_[k * n + i]);
            if (val > maxVal) {
                maxVal = val;
                maxRow = i;
            }
        }
        pivots_[k] = maxRow;

        // Swap rows
        if (maxRow != k) {
            for (int j = 0; j < n; ++j) {
                std::swap(systemLU_[j * n + k], systemLU_[j * n + maxRow]);
            }
        }

        // Eliminate
        Complex pivot = systemLU_[k * n + k];
        if (std::abs(pivot) < 1e-15) {
            throw std::runtime_error("SIBC FEM: Singular matrix");
        }

        for (int i = k + 1; i < n; ++i) {
            Complex factor = systemLU_[k * n + i] / pivot;
            systemLU_[k * n + i] = factor;  // Store L factor
            for (int j = k + 1; j < n; ++j) {
                systemLU_[j * n + i] -= factor * systemLU_[j * n + k];
            }
        }
    }
    #endif

    factorized_ = true;
}

// ============================================================================
// radTNonlocalSIBC implementation
// ============================================================================

radTNonlocalSIBC::radTNonlocalSIBC()
    : crossSectionType_("rectangular")
    , width_(0.001)
    , height_(0.001)
    , conductivity_(1e7)
    , relPermeability_real_(1.0)
    , relPermeability_imag_(0.0)
    , useComplexPermeability_(false)
    , frequency_(0)
    , omega_(0)
    , skinDepth_(1e30)
    , n_panels_(0)
    , impedanceComputed_(false)
{
    fem2D_ = std::make_unique<radTCrossSection2DFEM>();
}

radTNonlocalSIBC::~radTNonlocalSIBC() {
}

void radTNonlocalSIBC::SetConductorGeometry(const std::string& crossSection,
                                             double width, double height) {
    crossSectionType_ = crossSection;
    width_ = width;
    height_ = (height > 0) ? height : width;
    impedanceComputed_ = false;

    // Generate mesh
    if (crossSection == "circular") {
        fem2D_->MeshCircle(width_ / 2, 8, 24);
    } else {
        fem2D_->MeshRectangle(width_, height_, 10, 10);
    }
}

void radTNonlocalSIBC::SetMaterial(double conductivity, double relPermeability) {
    conductivity_ = conductivity;
    relPermeability_real_ = relPermeability;
    relPermeability_imag_ = 0.0;
    useComplexPermeability_ = false;
    fem2D_->SetMaterial(conductivity, relPermeability);
    impedanceComputed_ = false;
    UpdateSkinDepth();
}

void radTNonlocalSIBC::SetComplexMaterial(double conductivity,
                                           double mu_prime_r,
                                           double mu_double_prime_r) {
    conductivity_ = conductivity;
    relPermeability_real_ = mu_prime_r;
    relPermeability_imag_ = mu_double_prime_r;
    useComplexPermeability_ = (mu_double_prime_r > 1e-10);
    fem2D_->SetComplexMaterial(conductivity, mu_prime_r, mu_double_prime_r);
    impedanceComputed_ = false;
    UpdateSkinDepth();
}

void radTNonlocalSIBC::SetFrequency(double frequency) {
    frequency_ = frequency;
    omega_ = 2.0 * RadConst::PI * frequency;
    fem2D_->SetFrequency(frequency);
    impedanceComputed_ = false;
    UpdateSkinDepth();
}

void radTNonlocalSIBC::UpdateSkinDepth() {
    if (conductivity_ > 0 && omega_ > 0) {
        // For complex μ = μ' - jμ", use |μ| for skin depth estimation
        double mu_abs = MU_0 * std::sqrt(relPermeability_real_ * relPermeability_real_ +
                                          relPermeability_imag_ * relPermeability_imag_);
        skinDepth_ = std::sqrt(2.0 / (omega_ * mu_abs * conductivity_));
    } else {
        skinDepth_ = 1e30;
    }
}

void radTNonlocalSIBC::ComputeImpedanceMatrix(int n_panels, std::vector<Complex>& Z_sibc) {
    n_panels_ = n_panels;

    // For nonlocal SIBC, we compute the full impedance matrix from 2D FEM
    fem2D_->AssembleSystem();
    fem2D_->ComputeImpedanceMatrix(Z_sibc_);

    // The Z_sibc matrix relates boundary H to boundary E
    // For surface panels along the wire, we need to map this to panel interactions

    // Simplified: for uniform cross-section along wire, the impedance matrix
    // relates currents at different axial positions

    // Map 2D boundary nodes to 3D surface panels
    int nBoundary = fem2D_->NumBoundaryNodes();

    // For now, average the impedance for each panel
    Z_sibc.resize(n_panels * n_panels);

    // Simplified model: diagonal dominance with coupling
    Complex Z_avg = GetLocalSurfaceImpedance();

    for (int i = 0; i < n_panels; ++i) {
        for (int j = 0; j < n_panels; ++j) {
            if (i == j) {
                // Diagonal: use local SIBC as approximation
                Z_sibc[i * n_panels + j] = Z_avg;
            } else {
                // Off-diagonal: coupling from nonlocal effects
                // Decays with axial distance
                int dist = std::abs(i - j);
                double decay = std::exp(-dist * skinDepth_ / (width_ / n_panels));
                Z_sibc[i * n_panels + j] = Z_avg * Complex(0.1 * decay, 0);
            }
        }
    }

    impedanceComputed_ = true;
}

void radTNonlocalSIBC::Apply(const std::vector<Complex>& K, std::vector<Complex>& Et) {
    if (!impedanceComputed_) {
        ComputeImpedanceMatrix(static_cast<int>(K.size()), Z_sibc_);
    }

    int n = static_cast<int>(K.size());
    Et.resize(n, Complex(0, 0));

    // Et = Z_sibc * K
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Et[i] += Z_sibc_[i * n + j] * K[j];
        }
    }
}

Complex radTNonlocalSIBC::GetLocalSurfaceImpedance() const {
    // For real permeability: Zs = (1+j) / (σδ) where δ = skin depth
    // For complex permeability μ = μ' - jμ":
    //   Zs = sqrt(jωμ / σ)
    //
    // Physical interpretation:
    //   - Real(Zs): Surface resistance, determines ohmic + magnetic loss
    //   - Imag(Zs): Surface reactance, determines energy storage

    if (conductivity_ <= 0 || omega_ <= 0) {
        return Complex(0, 0);
    }

    if (useComplexPermeability_) {
        // Complex permeability: Zs = sqrt(jω(μ' - jμ") / σ)
        Complex mu_complex = GetComplexRelPermeability() * MU_0;
        Complex Z_arg = Complex(0, omega_) * mu_complex / conductivity_;
        return std::sqrt(Z_arg);
    } else {
        // Real permeability: Zs = (1+j) / (σδ)
        if (skinDepth_ > 1e20) {
            return Complex(0, 0);
        }
        return Complex(1, 1) / (conductivity_ * skinDepth_);
    }
}

bool radTNonlocalSIBC::IsLocalSIBCSufficient() const {
    // Local SIBC is accurate when skin depth << conductor dimension
    double charSize = (crossSectionType_ == "circular") ? width_ : std::min(width_, height_);
    return skinDepth_ < 0.1 * charSize;
}

double radTNonlocalSIBC::GetDCResistance() const {
    // R_DC = 1 / (σ * A)
    double area;
    if (crossSectionType_ == "circular") {
        area = RadConst::PI * (width_ / 2) * (width_ / 2);
    } else {
        area = width_ * height_;
    }
    return 1.0 / (conductivity_ * area);
}

double radTNonlocalSIBC::GetInternalInductance() const {
    // L_internal = μ / (8π) for circular wire (high frequency)
    // For rectangular, approximately μ*width/(4*height)
    // Use real part of permeability for inductance (reactive energy storage)
    double mu = MU_0 * relPermeability_real_;

    if (crossSectionType_ == "circular") {
        return mu / (8 * RadConst::PI);
    } else {
        // Approximate for rectangular
        return mu * width_ / (4 * height_);
    }
}

// ============================================================================
// radTSIBCElement implementation
// ============================================================================

radTSIBCElement::radTSIBCElement()
    : area_(0)
{
    center_ = TVector3d(0, 0, 0);
    normal_ = TVector3d(0, 0, 1);
    tangent1_ = TVector3d(1, 0, 0);
    tangent2_ = TVector3d(0, 1, 0);
}

radTSIBCElement::~radTSIBCElement() {
}

void radTSIBCElement::SetGeometry(const TVector3d& center,
                                   const TVector3d& normal,
                                   double area,
                                   const TVector3d& tangent1,
                                   const TVector3d& tangent2) {
    center_ = center;
    normal_ = normal;
    area_ = area;
    tangent1_ = tangent1;
    tangent2_ = tangent2;
}

void radTSIBCElement::SetSIBC(std::shared_ptr<radTNonlocalSIBC> sibc) {
    sibc_ = sibc;
}

void radTSIBCElement::ComputeTangentialE(std::vector<Complex>& Et) const {
    if (sibc_) {
        sibc_->Apply(K_, Et);
    } else {
        Et.assign(K_.size(), Complex(0, 0));
    }
}

// ============================================================================
// radTSIBCSolver implementation
// ============================================================================

radTSIBCSolver::radTSIBCSolver()
    : frequency_(0)
    , conductivity_(1e7)
    , relPermeability_(1000)
    , compareWithLocal_(false)
    , localVsNonlocalError_(0)
{
    sibc_ = std::make_shared<radTNonlocalSIBC>();
}

radTSIBCSolver::~radTSIBCSolver() {
}

void radTSIBCSolver::AddElement(std::shared_ptr<radTSIBCElement> element) {
    elements_.push_back(element);
    element->SetSIBC(sibc_);
}

void radTSIBCSolver::SetFrequency(double frequency) {
    frequency_ = frequency;
    sibc_->SetFrequency(frequency);
}

void radTSIBCSolver::SetMaterial(double conductivity, double relPermeability) {
    conductivity_ = conductivity;
    relPermeability_ = relPermeability;
    sibc_->SetMaterial(conductivity, relPermeability);
}

void radTSIBCSolver::Solve() {
    BuildSystemMatrix();
    SolveLinearSystem();
    ExtractSolution();
}

Complex radTSIBCSolver::ComputeImpedance() {
    // Sum currents and voltages from elements
    Complex V_total(0, 0);
    Complex I_total(0, 0);

    for (const auto& elem : elements_) {
        const auto& K = elem->GetSurfaceCurrent();
        double area = elem->Area();

        // Current = K * perimeter (simplified)
        for (const auto& k : K) {
            I_total += k * area;
        }
    }

    // V = 1 (applied), Z = V / I
    if (std::abs(I_total) > 1e-15) {
        return Complex(1, 0) / I_total;
    }
    return Complex(0, 0);
}

std::vector<Complex> radTSIBCSolver::ImpedanceSweep(const std::vector<double>& frequencies) {
    std::vector<Complex> result(frequencies.size());

    for (size_t i = 0; i < frequencies.size(); ++i) {
        SetFrequency(frequencies[i]);
        Solve();
        result[i] = ComputeImpedance();
    }

    return result;
}

void radTSIBCSolver::BuildSystemMatrix() {
    int nElem = static_cast<int>(elements_.size());
    if (nElem == 0) return;

    // Simplified: 1 DOF per element (scalar current)
    int nDOF = nElem;

    systemMatrix_.resize(nDOF * nDOF, Complex(0, 0));
    rhs_.resize(nDOF, Complex(0, 0));
    solution_.resize(nDOF, Complex(0, 0));

    // Get SIBC impedance matrix
    std::vector<Complex> Z_sibc;
    sibc_->ComputeImpedanceMatrix(nElem, Z_sibc);

    // System matrix = Z_sibc (for simplified formulation)
    systemMatrix_ = Z_sibc;

    // RHS from excitation (1V applied)
    rhs_[0] = Complex(1, 0);
}

void radTSIBCSolver::SolveLinearSystem() {
    int n = static_cast<int>(rhs_.size());
    if (n == 0) return;

    solution_ = rhs_;
    std::vector<Complex> A = systemMatrix_;

    // Gaussian elimination
    for (int k = 0; k < n; ++k) {
        int maxRow = k;
        double maxVal = std::abs(A[k * n + k]);
        for (int i = k + 1; i < n; ++i) {
            double val = std::abs(A[i * n + k]);
            if (val > maxVal) {
                maxVal = val;
                maxRow = i;
            }
        }

        if (maxRow != k) {
            for (int j = 0; j < n; ++j) {
                std::swap(A[k * n + j], A[maxRow * n + j]);
            }
            std::swap(solution_[k], solution_[maxRow]);
        }

        for (int i = k + 1; i < n; ++i) {
            Complex factor = A[i * n + k] / A[k * n + k];
            for (int j = k; j < n; ++j) {
                A[i * n + j] -= factor * A[k * n + j];
            }
            solution_[i] -= factor * solution_[k];
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        for (int j = i + 1; j < n; ++j) {
            solution_[i] -= A[i * n + j] * solution_[j];
        }
        solution_[i] /= A[i * n + i];
    }
}

void radTSIBCSolver::ExtractSolution() {
    // Copy solution to elements
    for (size_t i = 0; i < elements_.size(); ++i) {
        auto& K = elements_[i]->GetSurfaceCurrent();
        K.resize(1);
        K[0] = solution_[i];
    }

    // Compare with local SIBC if requested
    if (compareWithLocal_) {
        Complex Z_local = sibc_->GetLocalSurfaceImpedance();
        Complex Z_nonlocal = ComputeImpedance();

        if (std::abs(Z_local) > 1e-15) {
            localVsNonlocalError_ = std::abs(Z_nonlocal - Z_local) / std::abs(Z_local);
        } else {
            localVsNonlocalError_ = 0;
        }
    }
}

} // namespace radia
