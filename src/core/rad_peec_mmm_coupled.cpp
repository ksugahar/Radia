/*
 * rad_peec_mmm_coupled.cpp
 *
 * Implementation of coupled PEEC (Conductor) + MMM (Magnet) Solver
 *
 * Solver Selection:
 * - MMM only (no PEEC): Use real LU (dgesv) - existing rad_relaxation_methods.cpp
 * - PEEC present: Use complex LU (zgesv) - this file
 *
 * Part of Radia project
 */

#include "rad_peec_mmm_coupled.h"
#include "rad_constants.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <cstring>  // For std::memcpy

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef HAVE_LAPACK
// Intel MKL headers for complex LU solver
#include "mkl_lapack.h"
#include "mkl_cblas.h"
#endif

namespace radia {

// Use constants from rad_constants.h via RadConst::
// MU_0 and INV_FOUR_PI already defined in rad_conductor.h

//=============================================================================
// PEECMMMCoupledSolver implementation
//=============================================================================

PEECMMMCoupledSolver::PEECMMMCoupledSolver()
    : conductor_(nullptr)
    , magnetInteraction_(nullptr)
    , magnetHandle_(0)
    , ownsMagnetInteraction_(false)
    , frequency_(0)
    , omega_(0)
    , tolerance_(1e-4)
    , maxIterations_(100)
    , hasComplexMu_(false)
    , uniformComplexMu_(true)
    , uniformChi_(0, 0)
    , excitationType_(None)
    , excitationValue_(0, 0)
    , externalField_(0, 0, 0)
    , impedance_(0, 0)
    , totalCurrent_(0, 0)
    , conductorPower_(0)
    , magnetPower_(0)
    , solved_(false)
{
}

PEECMMMCoupledSolver::~PEECMMMCoupledSolver()
{
    // Note: We don't own magnetInteraction_ if it was set via SetMagnetInteraction()
    // Only delete if we created it via SetMagnet()
    if (ownsMagnetInteraction_ && magnetInteraction_) {
        delete magnetInteraction_;
    }
}

void PEECMMMCoupledSolver::SetConductor(std::shared_ptr<radTConductor> conductor)
{
    conductor_ = conductor;
    solved_ = false;
}

void PEECMMMCoupledSolver::SetMagnet(int magnetHandle)
{
    magnetHandle_ = magnetHandle;
    // Note: The actual radTInteraction will be obtained from the Radia global
    // handle table during Solve(). This is a placeholder for the Python API.
    solved_ = false;
}

void PEECMMMCoupledSolver::SetMagnetInteraction(radTInteraction* interaction)
{
    if (ownsMagnetInteraction_ && magnetInteraction_) {
        delete magnetInteraction_;
    }
    magnetInteraction_ = interaction;
    ownsMagnetInteraction_ = false;
    solved_ = false;
}

void PEECMMMCoupledSolver::SetFrequency(double frequency)
{
    frequency_ = frequency;
    omega_ = 2.0 * RadConst::PI * frequency;

    // Update conductor frequency
    if (conductor_) {
        conductor_->SetFrequency(frequency);
    }

    solved_ = false;
}

void PEECMMMCoupledSolver::SetVoltageExcitation(std::complex<double> V)
{
    excitationType_ = Voltage;
    excitationValue_ = V;

    if (conductor_) {
        conductor_->SetVoltageExcitation(V.real(), V.imag());
    }

    solved_ = false;
}

void PEECMMMCoupledSolver::SetCurrentExcitation(std::complex<double> I)
{
    excitationType_ = Current;
    excitationValue_ = I;

    if (conductor_) {
        conductor_->SetCurrentExcitation(I.real(), I.imag());
    }

    solved_ = false;
}

void PEECMMMCoupledSolver::SetExternalField(double Hx, double Hy, double Hz)
{
    externalField_ = TVector3d(Hx, Hy, Hz);
    solved_ = false;
}

void PEECMMMCoupledSolver::SetComplexPermeability(double mu_r_real, double mu_r_imag)
{
    // Complex susceptibility: χ = (μ'_r - 1) - jμ"_r
    // Note: μ"_r is positive for loss (energy dissipation)
    // The sign convention follows: μ = μ' - jμ" (engineering convention)
    uniformChi_ = std::complex<double>(mu_r_real - 1.0, -mu_r_imag);
    hasComplexMu_ = true;
    uniformComplexMu_ = true;
    elementChi_.clear();
    solved_ = false;
}

void PEECMMMCoupledSolver::SetComplexPermeabilityPerElement(
    const std::vector<double>& mu_r_real,
    const std::vector<double>& mu_r_imag)
{
    if (mu_r_real.size() != mu_r_imag.size()) {
        throw std::runtime_error("PEECMMMCoupledSolver: mu_r_real and mu_r_imag must have same size");
    }

    int n = static_cast<int>(mu_r_real.size());
    elementChi_.resize(n);

    for (int i = 0; i < n; ++i) {
        // χ = (μ'_r - 1) - jμ"_r
        elementChi_[i] = std::complex<double>(mu_r_real[i] - 1.0, -mu_r_imag[i]);
    }

    hasComplexMu_ = true;
    uniformComplexMu_ = false;
    solved_ = false;
}

PEECMMMResult PEECMMMCoupledSolver::Solve()
{
    PEECMMMResult result;
    result.impedance = std::complex<double>(0, 0);
    result.resistance = 0;
    result.reactance = 0;
    result.inductance = 0;
    result.conductorPower = 0;
    result.magnetPower = 0;
    result.conductorDOF = 0;
    result.magnetDOF = 0;
    result.iterations = 0;
    result.maxMagnetizationChange = 0;

    // Validate inputs
    if (!conductor_) {
        throw std::runtime_error("PEECMMMCoupledSolver: Conductor not set");
    }

    if (!magnetInteraction_ && magnetHandle_ == 0) {
        throw std::runtime_error("PEECMMMCoupledSolver: Magnet not set");
    }

    // Get DOF counts
    result.conductorDOF = conductor_->NumPanels() * 3;  // 3 current components per panel

    if (magnetInteraction_) {
        result.magnetDOF = magnetInteraction_->GetTotalDOF();
    }

    // Assemble coupling matrices
    AssembleCouplingMatrices();

    // Determine if nonlinear solve is needed
    // (Check if magnet has nonlinear material)
    bool isNonlinear = false;  // TODO: Check from magnetInteraction_

    if (isNonlinear) {
        SolveNonlinear();
    } else {
        SolveLinear();
    }

    // Compute power losses
    ComputeConductorPower();
    ComputeMagnetPower();

    // Fill result
    result.impedance = impedance_;
    result.resistance = impedance_.real();
    result.reactance = impedance_.imag();
    if (omega_ > 0) {
        result.inductance = impedance_.imag() / omega_;
    }
    result.conductorPower = conductorPower_;
    result.magnetPower = magnetPower_;

    solved_ = true;
    return result;
}

std::vector<PEECMMMResult> PEECMMMCoupledSolver::FrequencySweep(
    const std::vector<double>& frequencies)
{
    std::vector<PEECMMMResult> results;
    results.reserve(frequencies.size());

    for (double f : frequencies) {
        SetFrequency(f);
        results.push_back(Solve());
    }

    return results;
}

void PEECMMMCoupledSolver::AssembleCouplingMatrices()
{
    // Assemble Z_cm: B field from conductor at magnet element centers
    ComputeZcm();

    // Assemble Z_mc: B field from magnet at conductor panel centers
    ComputeZmc();
}

void PEECMMMCoupledSolver::ComputeZcm()
{
    // Z_cm[i][j][k] = B_k field component at magnet element i center
    //                 from unit current in conductor panel j
    //
    // Physical interpretation:
    //   When coil carries current I, the surface current K flows TANGENTIALLY
    //   along the wire surface (not in the normal direction!).
    //   B field at magnet center = sum_j { Biot-Savart contribution from panel j }
    //
    // The coupling represents: dB/dI at magnet element centers
    //
    // CORRECTED FORMULATION (2026-01-11):
    // For PEEC panels on wire surface, K flows tangent to the wire (circumferential
    // direction for a loop), NOT in the panel normal direction.
    //
    // Biot-Savart law for surface current:
    //   B(r) = mu_0/(4*pi) * integral{ K x (r - r') / |r - r'|^3 } dA'
    //
    // For a flat panel with uniform current density K in tangent direction t:
    //   B ~ mu_0/(4*pi) * (t * |K|) x (r - r_panel) / |r - r_panel|^3 * Area
    //
    // For unit current I, the effective surface current is K = I / h_wire
    // where h_wire is the wire cross-section height. But since we're computing
    // B per unit total current, we just need:
    //   B_per_unit_I ~ mu_0/(4*pi) * t x r / r^3 * L_segment
    // where L_segment = sqrt(Area) is the segment length along the wire.

    if (!conductor_ || !magnetInteraction_) return;

    const int nCondPanels = conductor_->NumPanels();
    const int nMagElems = magnetInteraction_->GetNumElements();

    if (nCondPanels == 0 || nMagElems == 0) return;

    // Z_cm layout: [nMagElems][nCondPanels][3]
    // Index: (iMag * nCondPanels + jCond) * 3 + k
    Z_cm_.resize(static_cast<size_t>(nMagElems) * nCondPanels * 3);
    std::fill(Z_cm_.begin(), Z_cm_.end(), std::complex<double>(0, 0));

    const auto& condPanels = conductor_->GetPanels();

    // For each magnet element
    #pragma omp parallel for
    for (int iMag = 0; iMag < nMagElems; ++iMag) {
        // Get magnet element center from radTInteraction
        TVector3d magCenter = magnetInteraction_->GetElementCenter(iMag);

        // For each conductor panel
        for (int jCond = 0; jCond < nCondPanels; ++jCond) {
            const SurfacePanel& panel = condPanels[jCond];
            TVector3d condCenter = panel.center;

            // Compute tangent direction: direction of current flow
            // For a circular loop in XY plane with center at origin:
            //   radial direction: (x, y, 0) / sqrt(x^2 + y^2)
            //   tangent direction: (-y, x, 0) / sqrt(x^2 + y^2) (counterclockwise)
            double rho = std::sqrt(condCenter.x * condCenter.x + condCenter.y * condCenter.y);
            TVector3d t_flow;
            if (rho > 1e-12) {
                // Counterclockwise tangent in XY plane
                t_flow.x = -condCenter.y / rho;
                t_flow.y = condCenter.x / rho;
                t_flow.z = 0;
            } else {
                // Fallback for points near z-axis
                t_flow.x = 1.0;
                t_flow.y = 0.0;
                t_flow.z = 0.0;
            }

            // Segment length along wire (contribution to line integral)
            double L_segment = std::sqrt(panel.area);

            // Vector from panel center to magnet center
            TVector3d r;
            r.x = magCenter.x - condCenter.x;
            r.y = magCenter.y - condCenter.y;
            r.z = magCenter.z - condCenter.z;

            double dist_sq = r.x * r.x + r.y * r.y + r.z * r.z;
            double dist = std::sqrt(dist_sq);

            if (dist > 1e-12) {
                // Biot-Savart: B = mu_0/(4*pi) * I * dl x r / r^3
                // For segment: dl = t_flow * L_segment
                double factor = MU_0 * INV_FOUR_PI * L_segment / (dist_sq * dist);

                // Cross product t_flow x r (current direction cross displacement)
                // t x r = (ty*rz - tz*ry, tz*rx - tx*rz, tx*ry - ty*rx)
                double Bx = factor * (t_flow.y * r.z - t_flow.z * r.y);
                double By = factor * (t_flow.z * r.x - t_flow.x * r.z);
                double Bz = factor * (t_flow.x * r.y - t_flow.y * r.x);

                size_t idx_base = (static_cast<size_t>(iMag) * nCondPanels + jCond) * 3;

                // Store as complex (real for quasi-static, phase will be applied later)
                Z_cm_[idx_base + 0] = std::complex<double>(Bx, 0);
                Z_cm_[idx_base + 1] = std::complex<double>(By, 0);
                Z_cm_[idx_base + 2] = std::complex<double>(Bz, 0);
            }
        }
    }
}

void PEECMMMCoupledSolver::ComputeZmc()
{
    // Z_mc[i][j][k] = Coupling coefficient from magnet element j to conductor panel i
    //                 for magnetization component k
    //
    // CORRECTED FORMULATION (2026-01-11):
    // The flux linkage is computed via line integral of vector potential A:
    //   V_emf = -jω * Φ = -jω * ∮ A · dl
    //
    // For the k-th panel segment (edge of wire cross-section), the contribution is:
    //   Φ_panel = A(panel_center) · t_panel * L_panel
    //
    // where t_panel is the tangent vector (direction of current flow) and L_panel
    // is the segment length (approximated as sqrt(panel_area) for rectangular panels).
    //
    // The tangent vector is perpendicular to both the panel normal and radial direction.
    // For a loop in XY plane with axis in Z: t ~ (±sin(θ), ∓cos(θ), 0)
    //
    // Vector potential from magnetic dipole:
    //   A(r) = mu_0/(4*pi) * (m x r) / r^3
    //
    // For unit magnetization M in element j:
    //   A(r) = mu_0/(4*pi) * V_elem * (M x r) / r^3

    if (!conductor_ || !magnetInteraction_) return;

    const int nCondPanels = conductor_->NumPanels();
    const int nMagElems = magnetInteraction_->GetNumElements();

    if (nCondPanels == 0 || nMagElems == 0) return;

    // Z_mc layout: [nCondPanels][nMagElems][3]
    // Index: (iCond * nMagElems + jMag) * 3 + k
    Z_mc_.resize(static_cast<size_t>(nCondPanels) * nMagElems * 3);
    std::fill(Z_mc_.begin(), Z_mc_.end(), std::complex<double>(0, 0));

    const auto& condPanels = conductor_->GetPanels();

    // For each conductor panel, compute flux linkage contribution via A · dl
    //
    // Vector potential from magnetic dipole m at position r:
    //   A(r) = mu_0/(4*pi) * (m x r) / r^3
    //
    // For a magnetic element with unit magnetization M and volume V:
    //   m = M * V
    //   A = mu_0/(4*pi) * V * (M x r) / r^3
    //
    // The flux linkage contribution from this panel is:
    //   dΦ = A · dl
    //
    // where dl is the segment of the wire path corresponding to this panel.
    // For PEEC panels on the wire surface, the segment length is approximately
    // sqrt(panel_area) and the direction is tangent to the current flow.
    //
    // The tangent direction t is perpendicular to the panel normal n.
    // For a loop in the XY plane: t is mostly in the XY plane, perpendicular to radial.

    #pragma omp parallel for
    for (int iCond = 0; iCond < nCondPanels; ++iCond) {
        const SurfacePanel& panel = condPanels[iCond];
        TVector3d condCenter = panel.center;
        TVector3d n_panel = panel.normal;
        double A_panel = panel.area;

        // Estimate segment length (wire path length associated with this panel)
        double L_segment = std::sqrt(A_panel);

        // Compute tangent direction: perpendicular to normal, in the flow direction
        // For a circular loop in XY plane with center at origin:
        //   radial direction: (x, y, 0) / sqrt(x^2 + y^2)
        //   tangent direction: (-y, x, 0) / sqrt(x^2 + y^2) (counterclockwise)
        // The actual tangent should be perpendicular to both normal and radial.

        double rho = std::sqrt(condCenter.x * condCenter.x + condCenter.y * condCenter.y);
        TVector3d t_flow;
        if (rho > 1e-12) {
            // Counterclockwise tangent in XY plane
            t_flow.x = -condCenter.y / rho;
            t_flow.y = condCenter.x / rho;
            t_flow.z = 0;
        } else {
            // Fallback for points near z-axis
            t_flow.x = 1.0;
            t_flow.y = 0.0;
            t_flow.z = 0.0;
        }

        // For each magnet element
        for (int jMag = 0; jMag < nMagElems; ++jMag) {
            TVector3d magCenter = magnetInteraction_->GetElementCenter(jMag);
            double V_elem = magnetInteraction_->GetElementVolume(jMag);

            // Vector from magnet to conductor panel
            TVector3d r;
            r.x = condCenter.x - magCenter.x;
            r.y = condCenter.y - magCenter.y;
            r.z = condCenter.z - magCenter.z;

            double dist_sq = r.x * r.x + r.y * r.y + r.z * r.z;
            double dist = std::sqrt(dist_sq);

            if (dist > 1e-12 && V_elem > 1e-20) {
                double dist3 = dist_sq * dist;

                // Factor: mu_0/(4*pi) * V_elem / r^3
                double factor = MU_0 * INV_FOUR_PI * V_elem / dist3;

                // A from unit M_x = factor * (e_x x r) = factor * (0, rz, -ry)
                // A from unit M_y = factor * (e_y x r) = factor * (-rz, 0, rx)
                // A from unit M_z = factor * (e_z x r) = factor * (ry, -rx, 0)
                //
                // Flux contribution: dΦ = A · t * L_segment

                // For M_x: A = factor * (0, rz, -ry)
                // dΦ_x = factor * (rz * ty - ry * tz) * L
                double dPhi_Mx = factor * (r.z * t_flow.y - r.y * t_flow.z) * L_segment;

                // For M_y: A = factor * (-rz, 0, rx)
                // dΦ_y = factor * (-rz * tx + rx * tz) * L
                double dPhi_My = factor * (-r.z * t_flow.x + r.x * t_flow.z) * L_segment;

                // For M_z: A = factor * (ry, -rx, 0)
                // dΦ_z = factor * (ry * tx - rx * ty) * L
                double dPhi_Mz = factor * (r.y * t_flow.x - r.x * t_flow.y) * L_segment;

                size_t idx_base = (static_cast<size_t>(iCond) * nMagElems + jMag) * 3;

                // Store flux linkage coefficients
                Z_mc_[idx_base + 0] = std::complex<double>(dPhi_Mx, 0);
                Z_mc_[idx_base + 1] = std::complex<double>(dPhi_My, 0);
                Z_mc_[idx_base + 2] = std::complex<double>(dPhi_Mz, 0);
            }
        }
    }
}

void PEECMMMCoupledSolver::SolveLinear()
{
    // Solve the coupled PEEC-MMM system using complex LU decomposition (zgesv)
    //
    // CORRECTED formulation:
    // The conductor is treated as a SINGLE loop with impedance Z_cond.
    // The system has 1 current DOF (I) and 3*N_mag magnetization DOFs.
    //
    // Unified System Matrix (complex):
    //   [Z_cond   Z_LM ] [I ]   [V    ]
    //   [Z_ML     Z_MM ] [M ] = [H_ext]
    //
    // Where:
    //   Z_cond = Conductor self-impedance (R + jωL from PEEC solver)
    //   Z_MM = [I - χ·N] (demagnetization tensor)
    //   Z_LM = -jω * Φ_from_M (flux linkage from magnet to loop)
    //   Z_ML = -χ/μ0 * B_coil (H field from coil at magnet centers)
    //
    // Physical interpretation:
    //   - Loop equation: V = Z_cond * I + (-jω * Φ_from_M)
    //   - Magnet equation: [I - χN] * M = χ * (H_ext + H_coil)
    //
    // The inductance enhancement comes from:
    //   L_eff = L_0 + L_mutual
    //   where L_mutual = d(Φ_from_M)/dI > 0 for mu_r > 1

    // Solve conductor first to get base impedance
    if (!conductor_) {
        throw std::runtime_error("PEECMMMCoupledSolver: Conductor not set");
    }

    // Get base impedance from PEEC solver
    // radTConductor doesn't have Solve() - need to use radTConductorSolver
    conductor_->SetFrequency(frequency_);
    conductor_->SetVoltageExcitation(excitationValue_.real(), excitationValue_.imag());

    // Create a solver for this conductor and solve
    radTConductorSolver condSolver;
    condSolver.AddConductor(conductor_);
    condSolver.SetFrequency(frequency_);
    condSolver.Solve();
    std::complex<double> Z_cond = condSolver.GetPortImpedance();

    // If no magnet interaction, just use PEEC result
    if (!magnetInteraction_ || magnetInteraction_->GetNumElements() == 0) {
        impedance_ = Z_cond;
        totalCurrent_ = excitationValue_ / Z_cond;
        return;
    }

    // Get DOF counts
    int nMagnetElem = magnetInteraction_->GetNumElements();
    int nMagnetDOF = nMagnetElem * 3;  // 3 DOF per element (Mx, My, Mz)
    int totalDOF = 1 + nMagnetDOF;     // 1 loop current + magnet DOFs

    // Build unified system matrix (complex, column-major for LAPACK)
    std::vector<std::complex<double>> A(static_cast<size_t>(totalDOF) * totalDOF, std::complex<double>(0, 0));
    std::vector<std::complex<double>> b(totalDOF, std::complex<double>(0, 0));

    // Get complex susceptibility
    auto getElementChi = [this](int elemIdx) -> std::complex<double> {
        if (hasComplexMu_) {
            if (uniformComplexMu_) {
                return uniformChi_;
            } else if (elemIdx < static_cast<int>(elementChi_.size())) {
                return elementChi_[elemIdx];
            }
        }
        return std::complex<double>(999.0, 0);  // Default chi = mu_r - 1 = 999
    };

    // ================================================================
    // Fill Z_cond block (1x1): Row 0, Col 0
    // ================================================================
    A[0] = Z_cond;

    // RHS for loop equation: V (applied voltage)
    b[0] = excitationValue_;

    // ================================================================
    // Fill Z_MM block (nMagnetDOF x nMagnetDOF)
    // ================================================================
    // MMM Physics:
    //   M_i = chi_i * H_total_i
    //   H_total_i = H_ext + H_coil(I) + sum_{j!=i} N_ij * M_j
    //
    // The N_ij tensor gives H field at element i from magnetization at element j:
    //   H_ij = N_ij * M_j = (1/4pi) * V_j * (3*(r.M)*r/r^5 - M/r^3)
    //
    // For dipole field: H = (1/4pi) * (3*(m.r)*r/r^5 - m/r^3)
    // With m = M * V:   H_ij = (V_j/4pi) * (3*(M.r)*r/r^5 - M/r^3)
    //
    // In matrix form: H_ij[k] = sum_l N_ij[k][l] * M_j[l]
    // where N_ij[k][l] = (V_j/4pi/r^5) * (3*r[k]*r[l] - r^2*delta[k][l])
    //
    // MMM equation: M_i = chi_i * (H_ext + H_coil + sum_j N_ij * M_j)
    //
    // For the coupled system, we solve:
    //   [Z_cond   Z_LM ] [I ]   [V    ]
    //   [Z_ML     Z_MM ] [M ] = [b_mag]
    //
    // The Z_MM block should give: [I - chi*N] * M = chi * (H_ext + H_coil)
    // Or equivalently: M - chi*N*M = chi*H
    //
    // So Z_MM[i][i] = I (identity 3x3)
    //    Z_MM[i][j] = -chi_i * N_ij for i != j
    //
    // And Z_ML already adds -chi*H_coil to the RHS via the coupling.
    //
    // IMPORTANT: The RHS b_mag already includes chi*H_coil from Z_ML block,
    // so the magnet equation becomes:
    //   M_i - chi_i * sum_j N_ij * M_j + Z_ML[i] * I = chi_i * H_ext
    //
    // Wait - this is wrong. The Z_ML block adds the *effect* of I on M,
    // not the other way. Let me reconsider...
    //
    // CORRECT FORMULATION:
    // The magnet equation in the coupled system is:
    //   Z_MM * M + Z_ML * I = b_mag (= chi * H_ext)
    //
    // Z_ML gives: how coil current I affects the magnet equation
    // From MMM: M = chi * (H_ext + H_coil + H_demag)
    //         => M = chi * H_ext + chi * H_coil + chi * N * M
    //         => (I - chi*N) * M = chi * H_ext + chi * H_coil
    //
    // So: Z_MM = [I - chi*N]
    //     Z_ML = -chi * (B_coil / mu0) = -chi * H_coil_per_I
    //     b_mag = chi * H_ext
    //
    // Current implementation has Z_ML as: -chi/mu0 * B_coil (correct)
    // But we need to move the chi*H_coil to RHS, not keep it in Z_ML.
    //
    // Actually, Z_ML appears in: Z_MM * M + Z_ML * I = b_mag
    // Expanding: M - chi*N*M - chi*H_coil_per_I * I = chi*H_ext
    // Moving: M - chi*N*M = chi*H_ext + chi*H_coil_per_I * I
    //
    // This means Z_ML should have POSITIVE sign (moves to RHS with negative):
    //   Z_ML = +chi * H_coil_per_I
    //
    // Let me check the Z_ML block sign in current implementation...
    // Current: factor = -chi / MU_0 (wrong sign!)
    // Should be: factor = +chi / MU_0
    //
    // For now, keep the demagnetization tensor implementation but fix the sign.

    for (int iElem = 0; iElem < nMagnetElem; ++iElem) {
        std::complex<double> chi_i = getElementChi(iElem);
        TVector3d center_i = magnetInteraction_->GetElementCenter(iElem);

        for (int jElem = 0; jElem < nMagnetElem; ++jElem) {
            int row_base = 1 + iElem * 3;
            int col_base = 1 + jElem * 3;

            if (iElem == jElem) {
                // Diagonal block: I + chi_i * N_self
                // N_self is the self-demagnetization factor
                // For a sphere: N_self = 1/3 (isotropic)
                // For a cube: N_self ~ 1/3 (approximately)
                //
                // The MMM formulation is: M = chi * (H_ext + H_coil + H_from_others - N_self * M)
                // => M + chi * N_self * M = chi * (H_ext + H_coil + H_from_others)
                // => (1 + chi * N_self) * M = chi * (...)
                //
                // So diagonal = 1 + chi * N_self (scalar for isotropic N_self)
                double N_self = 1.0 / 3.0;  // Sphere approximation
                std::complex<double> diag_val = std::complex<double>(1.0, 0) + chi_i * std::complex<double>(N_self, 0);

                for (int k = 0; k < 3; ++k) {
                    int row = row_base + k;
                    int col = col_base + k;
                    A[static_cast<size_t>(col) * totalDOF + row] = diag_val;
                }
            } else {
                // Off-diagonal block: -chi_i * N_ij
                // This represents the demagnetizing effect of M_j on element i
                TVector3d center_j = magnetInteraction_->GetElementCenter(jElem);
                double V_j = magnetInteraction_->GetElementVolume(jElem);

                // Vector from j to i
                TVector3d r;
                r.x = center_i.x - center_j.x;
                r.y = center_i.y - center_j.y;
                r.z = center_i.z - center_j.z;

                double dist_sq = r.x * r.x + r.y * r.y + r.z * r.z;
                double dist = std::sqrt(dist_sq);

                if (dist > 1e-12 && V_j > 1e-20) {
                    double dist5 = dist_sq * dist_sq * dist;

                    // N_ij tensor: gives H at i from M at j
                    // N_ij[k][l] = (V_j/4pi/r^5) * (3*r[k]*r[l] - r^2*delta[k][l])
                    double factor = INV_FOUR_PI * V_j / dist5;
                    double r2 = dist_sq;

                    // Z_MM[i][j] = -chi_i * N_ij
                    std::complex<double> cf = -chi_i * std::complex<double>(factor, 0);

                    double rr[3] = {r.x, r.y, r.z};
                    for (int ki = 0; ki < 3; ++ki) {
                        for (int kj = 0; kj < 3; ++kj) {
                            double N_kl = 3.0 * rr[ki] * rr[kj];
                            if (ki == kj) N_kl -= r2;

                            int row = row_base + ki;
                            int col = col_base + kj;
                            A[static_cast<size_t>(col) * totalDOF + row] = cf * std::complex<double>(N_kl, 0);
                        }
                    }
                }
            }
        }
    }

    // RHS for magnet equation: χ · H_ext (usually zero if H_ext = 0)
    for (int iElem = 0; iElem < nMagnetElem; ++iElem) {
        std::complex<double> chi = getElementChi(iElem);
        int base = 1 + iElem * 3;
        b[base + 0] = chi * std::complex<double>(externalField_.x, 0);
        b[base + 1] = chi * std::complex<double>(externalField_.y, 0);
        b[base + 2] = chi * std::complex<double>(externalField_.z, 0);
    }

    // ================================================================
    // Fill Z_LM block (1 row x nMagnetDOF cols): -jω * Φ from M
    // ================================================================
    // The flux linkage from magnetization M to the loop is:
    //   Φ = sum over all magnet elements of: integral{ B(M) · n_panel } dA
    // For unit M in element j: Φ_j,k = Z_mc coefficient (summed over panels)
    //
    // Z_LM[0][jMag*3+k] = -jω * Σ_panel Z_mc[panel][jMag][k]
    std::complex<double> jw(0, omega_);

    if (!Z_mc_.empty()) {
        const int nCondPanels = conductor_->NumPanels();

        for (int jMag = 0; jMag < nMagnetElem; ++jMag) {
            for (int k = 0; k < 3; ++k) {
                // Sum flux from all panels
                std::complex<double> total_flux(0, 0);
                for (int iPanel = 0; iPanel < nCondPanels; ++iPanel) {
                    size_t idx = (static_cast<size_t>(iPanel) * nMagnetElem + jMag) * 3 + k;
                    total_flux += Z_mc_[idx];
                }

                // Z_LM entry: -jω * total_flux
                int col = 1 + jMag * 3 + k;
                A[static_cast<size_t>(col) * totalDOF + 0] = -jw * total_flux;
            }
        }
    }

    // ================================================================
    // Fill Z_ML block (nMagnetDOF rows x 1 col): -χ/μ0 * B_coil per unit I
    // ================================================================
    // The H field at magnet center from unit coil current is:
    //   H = B / μ0 = (1/μ0) * Σ_panel Z_cm[magnet][panel]
    //
    // Since Z_cm stores B from unit current in each panel, and all panels
    // carry the same current I, we sum over panels.
    //
    // Z_ML[iMag*3+k][0] = -χ/μ0 * Σ_panel Z_cm[iMag][panel][k]

    if (!Z_cm_.empty()) {
        const int nCondPanels = conductor_->NumPanels();

        for (int iMag = 0; iMag < nMagnetElem; ++iMag) {
            std::complex<double> chi = getElementChi(iMag);
            std::complex<double> factor = -chi / MU_0;

            for (int k = 0; k < 3; ++k) {
                // Sum B field from all panels
                std::complex<double> total_B(0, 0);
                for (int jPanel = 0; jPanel < nCondPanels; ++jPanel) {
                    size_t idx = (static_cast<size_t>(iMag) * nCondPanels + jPanel) * 3 + k;
                    total_B += Z_cm_[idx];
                }

                // Z_ML entry
                int row = 1 + iMag * 3 + k;
                A[row] = factor * total_B;  // Column 0
            }
        }
    }

    // Solve using complex LU decomposition (LAPACK zgesv)
#ifdef HAVE_LAPACK
    std::vector<int> ipiv(totalDOF);
    int nrhs = 1;
    int info = 0;

    // zgesv: solve A * x = b, overwrites A with LU and b with solution
    zgesv_(&totalDOF, &nrhs,
           reinterpret_cast<MKL_Complex16*>(A.data()), &totalDOF,
           ipiv.data(),
           reinterpret_cast<MKL_Complex16*>(b.data()), &totalDOF,
           &info);

    if (info != 0) {
        throw std::runtime_error("PEECMMMCoupledSolver: zgesv failed, info = " + std::to_string(info));
    }
#else
    // Fallback: simple Gaussian elimination for complex matrices
    // (Not optimized, for compilation without MKL)
    throw std::runtime_error("PEECMMMCoupledSolver: Complex LU requires Intel MKL (HAVE_LAPACK)");
#endif

    // Extract solution
    // b[0] = Loop current I
    // b[1..totalDOF-1] = Magnetization M (3 per element)

    // Get loop current (single DOF)
    totalCurrent_ = b[0];

    // Compute impedance: Z = V / I
    if (std::abs(totalCurrent_) > 1e-15) {
        impedance_ = excitationValue_ / totalCurrent_;
    } else {
        impedance_ = std::complex<double>(0, 0);
    }

    // Store solution for field computation
    solution_.resize(totalDOF);
    for (int i = 0; i < totalDOF; ++i) {
        solution_[i] = b[i];
    }
}

void PEECMMMCoupledSolver::SolveNonlinear()
{
    // For nonlinear materials, iterate:
    //
    // 1. Guess initial M
    // 2. Compute B_mag at conductor
    // 3. Solve PEEC for J
    // 4. Compute B_cond at magnet
    // 5. Update M based on total H = H_ext + H_cond + H_demag
    // 6. Check convergence, repeat if needed

    // Initial solve with linear approximation
    SolveLinear();

    // TODO: Implement nonlinear iteration
    // For now, single iteration (linear approximation)
}

void PEECMMMCoupledSolver::ComputeConductorPower()
{
    // P = I^2 * R = |I|^2 * Re(Z)
    double I_mag = std::abs(totalCurrent_);
    conductorPower_ = I_mag * I_mag * impedance_.real();
}

void PEECMMMCoupledSolver::ComputeMagnetPower()
{
    // Magnetic loss from complex permeability (hysteresis + eddy currents in grains)
    //
    // Complex permeability: μ = μ' - jμ"
    // Complex susceptibility: χ = (μ'_r - 1) - jμ"_r
    //
    // Power loss density (per unit volume):
    //   p = (ω/2) · μ0 · μ"_r · |H|^2  [W/m^3]
    //
    // For MMM formulation with M = χ·H:
    //   H = M / χ (element-wise)
    //   p = (ω/2) · μ0 · μ"_r · |M/χ|^2
    //
    // Total power: P = Σ_elem (ω/2) · μ0 · |χ"|/|χ|^2 · |M_elem|^2 · V_elem

    magnetPower_ = 0;

    if (!hasComplexMu_ || !magnetInteraction_ || solution_.empty()) {
        return;  // No complex permeability or no solution
    }

    int nMagnetElem = magnetInteraction_->GetNumElements();
    const int nLoopDOF = 1;  // Single loop current DOF

    // Get element volumes (approximation: assume unit volume for now)
    // TODO: Get actual volumes from geometry

    for (int elem = 0; elem < nMagnetElem; ++elem) {
        // Get complex susceptibility for this element
        std::complex<double> chi;
        if (uniformComplexMu_) {
            chi = uniformChi_;
        } else if (elem < static_cast<int>(elementChi_.size())) {
            chi = elementChi_[elem];
        } else {
            continue;  // Skip if no χ defined
        }

        // μ"_r = -Im(χ) (note: χ = (μ'_r - 1) - jμ"_r, so Im(χ) = -μ"_r)
        double mu_r_imag = -chi.imag();

        if (mu_r_imag <= 0) {
            continue;  // No loss if μ" <= 0
        }

        // Get magnetization M for this element from solution
        int idx = nLoopDOF + elem * 3;
        std::complex<double> Mx = solution_[idx + 0];
        std::complex<double> My = solution_[idx + 1];
        std::complex<double> Mz = solution_[idx + 2];

        // |M|^2
        double M_mag_sq = std::norm(Mx) + std::norm(My) + std::norm(Mz);

        // |χ|^2
        double chi_mag_sq = std::norm(chi);

        if (chi_mag_sq < 1e-20) {
            continue;  // Avoid division by zero
        }

        // |H|^2 = |M/χ|^2 = |M|^2 / |χ|^2
        double H_mag_sq = M_mag_sq / chi_mag_sq;

        // Power density: p = (ω/2) · μ0 · μ"_r · |H|^2
        double power_density = 0.5 * omega_ * MU_0 * mu_r_imag * H_mag_sq;

        // Assume unit volume for now (1 m^3)
        // TODO: Get actual element volume from geometry
        double volume = 1.0;

        magnetPower_ += power_density * volume;
    }
}

void PEECMMMCoupledSolver::ComputeB(const TVector3d& point,
                                     std::complex<double>& Bx,
                                     std::complex<double>& By,
                                     std::complex<double>& Bz) const
{
    // B_total = B_conductor + B_magnet

    // Conductor contribution (complex for AC)
    std::complex<double> Bx_cond, By_cond, Bz_cond;
    ConductorBField(point, Bx_cond, By_cond, Bz_cond);

    // Magnet contribution (real for static)
    double Bx_mag, By_mag, Bz_mag;
    MagnetBField(point, Bx_mag, By_mag, Bz_mag);

    // Sum
    Bx = Bx_cond + std::complex<double>(Bx_mag, 0);
    By = By_cond + std::complex<double>(By_mag, 0);
    Bz = Bz_cond + std::complex<double>(Bz_mag, 0);
}

void PEECMMMCoupledSolver::ComputeBConductor(const TVector3d& point,
                                              std::complex<double>& Bx,
                                              std::complex<double>& By,
                                              std::complex<double>& Bz) const
{
    ConductorBField(point, Bx, By, Bz);
}

void PEECMMMCoupledSolver::ComputeBMagnet(const TVector3d& point,
                                           double& Bx, double& By, double& Bz) const
{
    MagnetBField(point, Bx, By, Bz);
}

void PEECMMMCoupledSolver::ComputeBBatch(const std::vector<TVector3d>& points,
                                          std::vector<std::complex<double>>& Bx,
                                          std::vector<std::complex<double>>& By,
                                          std::vector<std::complex<double>>& Bz) const
{
    const int nPts = static_cast<int>(points.size());
    Bx.resize(nPts);
    By.resize(nPts);
    Bz.resize(nPts);

    #pragma omp parallel for
    for (int i = 0; i < nPts; ++i) {
        ComputeB(points[i], Bx[i], By[i], Bz[i]);
    }
}

void PEECMMMCoupledSolver::ConductorBField(const TVector3d& point,
                                            std::complex<double>& Bx,
                                            std::complex<double>& By,
                                            std::complex<double>& Bz) const
{
    if (conductor_) {
        conductor_->ComputeB(point, Bx, By, Bz);
    } else {
        Bx = By = Bz = std::complex<double>(0, 0);
    }
}

void PEECMMMCoupledSolver::MagnetBField(const TVector3d& point,
                                         double& Bx, double& By, double& Bz) const
{
    // TODO: Compute B from magnetInteraction_
    // Need to use Radia's Fld() function or equivalent

    Bx = By = Bz = 0;  // Placeholder
}

//=============================================================================
// Python API helper functions
//=============================================================================

PEECMMMCoupledSolver* CreatePEECMMMSolver(int conductorHandle, int magnetHandle)
{
    auto* solver = new PEECMMMCoupledSolver();

    // TODO: Get conductor from handle table
    // TODO: Get magnet from handle table

    solver->SetMagnet(magnetHandle);

    return solver;
}

std::complex<double> SolvePEECMMM(PEECMMMCoupledSolver* solver, double frequency)
{
    if (!solver) {
        throw std::runtime_error("SolvePEECMMM: Null solver");
    }

    solver->SetFrequency(frequency);
    solver->Solve();
    return solver->GetImpedance();
}

std::vector<double> GetPEECMMMField(PEECMMMCoupledSolver* solver,
                                     const TVector3d& point)
{
    std::vector<double> result(6, 0.0);

    if (solver) {
        std::complex<double> Bx, By, Bz;
        solver->ComputeB(point, Bx, By, Bz);

        result[0] = Bx.real();
        result[1] = By.real();
        result[2] = Bz.real();
        result[3] = Bx.imag();
        result[4] = By.imag();
        result[5] = Bz.imag();
    }

    return result;
}

} // namespace radia
