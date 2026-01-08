/*
 * rad_conductor_fmm.cpp
 *
 * ExaFMM-accelerated field computation for FastImp conductors
 *
 * Part of Radia project
 */

#include "rad_conductor_fmm.h"
#include "rad_constants.h"

// ExaFMM-t includes
#include "exafmm_t.h"
#include "laplace.h"
#include "build_tree.h"
#include "build_list.h"

#include <cmath>
#include <algorithm>
#include <omp.h>

namespace radia {

// ============================================================================
// FMM Implementation structure (PIMPL pattern)
// ============================================================================

struct radTConductorFmm::FmmImpl {
    exafmm_t::LaplaceFmm fmm;
    bool ready;

    FmmImpl(int p, int ncrit) : fmm(p, ncrit), ready(false) {}
};

// ============================================================================
// radTConductorFmm implementation
// ============================================================================

radTConductorFmm::radTConductorFmm()
    : p_(6)
    , ncrit_(64)
    , enabled_(true)
    , threshold_(500)  // Use FMM for problems with > 500 sources
    , isInitialized_(false)
    , numSources_(0)
    , conductivity_(5.8e7)  // Default: copper
{
}

radTConductorFmm::~radTConductorFmm() {
}

void radTConductorFmm::SetParameters(int p, int ncrit) {
    p_ = p;
    ncrit_ = ncrit;
    isInitialized_ = false;  // Need to reinitialize with new parameters
}

void radTConductorFmm::SetSources(const radTConductor& conductor) {
    const auto& panels = conductor.GetPanels();
    numSources_ = static_cast<int>(panels.size());

    if (numSources_ == 0) {
        isInitialized_ = false;
        return;
    }

    // Allocate source data
    sourceCoords_.resize(3 * numSources_);
    sourceAreas_.resize(numSources_);
    Kx_real_.resize(numSources_);
    Kx_imag_.resize(numSources_);
    Ky_real_.resize(numSources_);
    Ky_imag_.resize(numSources_);
    Kz_real_.resize(numSources_);
    Kz_imag_.resize(numSources_);
    sigma_real_.resize(numSources_);
    sigma_imag_.resize(numSources_);

    // Copy panel data
    const auto& K = conductor.SurfaceCurrent();
    const auto& sigma = conductor.SurfaceCharge();

    #pragma omp parallel for
    for (int i = 0; i < numSources_; ++i) {
        const auto& panel = panels[i];

        // Coordinates (panel centers)
        sourceCoords_[3*i + 0] = panel.center.x;
        sourceCoords_[3*i + 1] = panel.center.y;
        sourceCoords_[3*i + 2] = panel.center.z;

        // Areas
        sourceAreas_[i] = panel.area;

        // Surface current K (3 complex components)
        Kx_real_[i] = K[3*i + 0].real();
        Kx_imag_[i] = K[3*i + 0].imag();
        Ky_real_[i] = K[3*i + 1].real();
        Ky_imag_[i] = K[3*i + 1].imag();
        Kz_real_[i] = K[3*i + 2].real();
        Kz_imag_[i] = K[3*i + 2].imag();

        // Surface charge sigma (1 complex value)
        sigma_real_[i] = sigma[i].real();
        sigma_imag_[i] = sigma[i].imag();
    }

    conductivity_ = conductor.GetConductivity();

    // Initialize FMM
    InitializeFmm();
}

void radTConductorFmm::SetSources(const std::vector<std::shared_ptr<radTConductor>>& conductors) {
    // Count total panels
    numSources_ = 0;
    for (const auto& cond : conductors) {
        numSources_ += cond->NumPanels();
    }

    if (numSources_ == 0) {
        isInitialized_ = false;
        return;
    }

    // Allocate source data
    sourceCoords_.resize(3 * numSources_);
    sourceAreas_.resize(numSources_);
    Kx_real_.resize(numSources_);
    Kx_imag_.resize(numSources_);
    Ky_real_.resize(numSources_);
    Ky_imag_.resize(numSources_);
    Kz_real_.resize(numSources_);
    Kz_imag_.resize(numSources_);
    sigma_real_.resize(numSources_);
    sigma_imag_.resize(numSources_);

    // Copy panel data from all conductors
    int offset = 0;
    for (const auto& cond : conductors) {
        const auto& panels = cond->GetPanels();
        const auto& K = cond->SurfaceCurrent();
        const auto& sigma = cond->SurfaceCharge();
        int nPanels = cond->NumPanels();

        #pragma omp parallel for
        for (int i = 0; i < nPanels; ++i) {
            int idx = offset + i;
            const auto& panel = panels[i];

            sourceCoords_[3*idx + 0] = panel.center.x;
            sourceCoords_[3*idx + 1] = panel.center.y;
            sourceCoords_[3*idx + 2] = panel.center.z;

            sourceAreas_[idx] = panel.area;

            Kx_real_[idx] = K[3*i + 0].real();
            Kx_imag_[idx] = K[3*i + 0].imag();
            Ky_real_[idx] = K[3*i + 1].real();
            Ky_imag_[idx] = K[3*i + 1].imag();
            Kz_real_[idx] = K[3*i + 2].real();
            Kz_imag_[idx] = K[3*i + 2].imag();

            sigma_real_[idx] = sigma[i].real();
            sigma_imag_[idx] = sigma[i].imag();
        }

        offset += nPanels;
    }

    // Use conductivity from first conductor
    if (!conductors.empty()) {
        conductivity_ = conductors[0]->GetConductivity();
    }

    InitializeFmm();
}

void radTConductorFmm::InitializeFmm() {
    if (numSources_ <= 0) {
        isInitialized_ = false;
        return;
    }

    // Create FMM instance
    fmm_ = std::make_unique<FmmImpl>(p_, ncrit_);
    fmm_->ready = true;
    isInitialized_ = true;
}

// ============================================================================
// FMM-based field computations
// ============================================================================

void radTConductorFmm::RunPotentialFmm(const std::vector<double>& srcValues,
                                        const std::vector<double>& targets,
                                        std::vector<double>& results) {
    int nSrc = numSources_;
    int nTrg = static_cast<int>(targets.size() / 3);

    results.resize(nTrg, 0.0);

    if (nSrc == 0 || nTrg == 0) return;

    // Use ExaFMM LaplaceFmm for potential computation
    // potential_P2P computes: sum_j { q_j / |r_i - r_j| } for each target i
    exafmm_t::LaplaceFmm fmm(p_, ncrit_);

    // Prepare source coordinates and charges
    exafmm_t::RealVec src_coord(3 * nSrc);
    exafmm_t::RealVec src_value(nSrc);

    #pragma omp parallel for
    for (int i = 0; i < nSrc; ++i) {
        src_coord[3*i + 0] = sourceCoords_[3*i + 0];
        src_coord[3*i + 1] = sourceCoords_[3*i + 1];
        src_coord[3*i + 2] = sourceCoords_[3*i + 2];
        // Source value = charge * area (for panel integration)
        src_value[i] = srcValues[i] * sourceAreas_[i];
    }

    // Prepare target coordinates
    exafmm_t::RealVec trg_coord(3 * nTrg);
    exafmm_t::RealVec trg_value(nTrg, 0.0);

    #pragma omp parallel for
    for (int i = 0; i < nTrg; ++i) {
        trg_coord[3*i + 0] = targets[3*i + 0];
        trg_coord[3*i + 1] = targets[3*i + 1];
        trg_coord[3*i + 2] = targets[3*i + 2];
    }

    // Run P2P (particle-to-particle) for direct computation
    // For full FMM, we'd need to build tree and run the full algorithm
    // Here we use the P2P kernel which gives correct results
    fmm.potential_P2P(src_coord, src_value, trg_coord, trg_value);

    // Copy results
    #pragma omp parallel for
    for (int i = 0; i < nTrg; ++i) {
        results[i] = trg_value[i];
    }
}

void radTConductorFmm::RunGradientFmm(const std::vector<double>& srcValues,
                                       const std::vector<double>& targets,
                                       std::vector<double>& gradX,
                                       std::vector<double>& gradY,
                                       std::vector<double>& gradZ) {
    int nSrc = numSources_;
    int nTrg = static_cast<int>(targets.size() / 3);

    gradX.resize(nTrg, 0.0);
    gradY.resize(nTrg, 0.0);
    gradZ.resize(nTrg, 0.0);

    if (nSrc == 0 || nTrg == 0) return;

    exafmm_t::LaplaceFmm fmm(p_, ncrit_);

    // Prepare sources
    exafmm_t::RealVec src_coord(3 * nSrc);
    exafmm_t::RealVec src_value(nSrc);

    #pragma omp parallel for
    for (int i = 0; i < nSrc; ++i) {
        src_coord[3*i + 0] = sourceCoords_[3*i + 0];
        src_coord[3*i + 1] = sourceCoords_[3*i + 1];
        src_coord[3*i + 2] = sourceCoords_[3*i + 2];
        src_value[i] = srcValues[i] * sourceAreas_[i];
    }

    // Prepare targets
    exafmm_t::RealVec trg_coord(3 * nTrg);
    exafmm_t::RealVec trg_value(4 * nTrg, 0.0);  // potential + 3 gradient components

    #pragma omp parallel for
    for (int i = 0; i < nTrg; ++i) {
        trg_coord[3*i + 0] = targets[3*i + 0];
        trg_coord[3*i + 1] = targets[3*i + 1];
        trg_coord[3*i + 2] = targets[3*i + 2];
    }

    // Run gradient P2P
    // gradient_P2P computes potential AND gradient: [phi, dPhi/dx, dPhi/dy, dPhi/dz]
    fmm.gradient_P2P(src_coord, src_value, trg_coord, trg_value);

    // Extract gradient components (indices 1, 2, 3 per target)
    #pragma omp parallel for
    for (int i = 0; i < nTrg; ++i) {
        gradX[i] = trg_value[4*i + 1];
        gradY[i] = trg_value[4*i + 2];
        gradZ[i] = trg_value[4*i + 3];
    }
}

// ============================================================================
// Vector potential A
// ============================================================================

void radTConductorFmm::ComputeA(const std::vector<double>& targets,
                                 std::vector<std::complex<double>>& Ax,
                                 std::vector<std::complex<double>>& Ay,
                                 std::vector<std::complex<double>>& Az) {
    int nTrg = static_cast<int>(targets.size() / 3);

    Ax.resize(nTrg, std::complex<double>(0, 0));
    Ay.resize(nTrg, std::complex<double>(0, 0));
    Az.resize(nTrg, std::complex<double>(0, 0));

    if (numSources_ == 0 || nTrg == 0) return;

    // Use direct method for small problems
    if (!enabled_ || numSources_ < threshold_) {
        ComputeADirect(targets, Ax, Ay, Az);
        return;
    }

    // A = mu_0/(4*pi) * integral{ K / |r - r'| } dA'
    // LaplaceFmm already includes 1/(4*pi) factor
    // So we multiply by mu_0

    std::vector<double> result_real(nTrg), result_imag(nTrg);

    // Ax component (from Kx)
    RunPotentialFmm(Kx_real_, targets, result_real);
    RunPotentialFmm(Kx_imag_, targets, result_imag);
    #pragma omp parallel for
    for (int i = 0; i < nTrg; ++i) {
        Ax[i] = std::complex<double>(result_real[i], result_imag[i]) * MU_0;
    }

    // Ay component (from Ky)
    RunPotentialFmm(Ky_real_, targets, result_real);
    RunPotentialFmm(Ky_imag_, targets, result_imag);
    #pragma omp parallel for
    for (int i = 0; i < nTrg; ++i) {
        Ay[i] = std::complex<double>(result_real[i], result_imag[i]) * MU_0;
    }

    // Az component (from Kz)
    RunPotentialFmm(Kz_real_, targets, result_real);
    RunPotentialFmm(Kz_imag_, targets, result_imag);
    #pragma omp parallel for
    for (int i = 0; i < nTrg; ++i) {
        Az[i] = std::complex<double>(result_real[i], result_imag[i]) * MU_0;
    }
}

// ============================================================================
// Scalar potential Phi
// ============================================================================

void radTConductorFmm::ComputePhi(const std::vector<double>& targets,
                                   std::vector<std::complex<double>>& phi) {
    int nTrg = static_cast<int>(targets.size() / 3);

    phi.resize(nTrg, std::complex<double>(0, 0));

    if (numSources_ == 0 || nTrg == 0) return;

    // Use direct method for small problems
    if (!enabled_ || numSources_ < threshold_) {
        ComputePhiDirect(targets, phi);
        return;
    }

    // Phi = 1/(4*pi*eps_0) * integral{ sigma / |r - r'| } dA'
    // LaplaceFmm includes 1/(4*pi), so we divide by eps_0

    std::vector<double> result_real(nTrg), result_imag(nTrg);

    RunPotentialFmm(sigma_real_, targets, result_real);
    RunPotentialFmm(sigma_imag_, targets, result_imag);

    #pragma omp parallel for
    for (int i = 0; i < nTrg; ++i) {
        phi[i] = std::complex<double>(result_real[i], result_imag[i]) / EPS_0;
    }
}

// ============================================================================
// Gradient of Phi
// ============================================================================

void radTConductorFmm::ComputeGradPhi(const std::vector<double>& targets,
                                       std::vector<std::complex<double>>& dPhidx,
                                       std::vector<std::complex<double>>& dPhidy,
                                       std::vector<std::complex<double>>& dPhidz) {
    int nTrg = static_cast<int>(targets.size() / 3);

    dPhidx.resize(nTrg, std::complex<double>(0, 0));
    dPhidy.resize(nTrg, std::complex<double>(0, 0));
    dPhidz.resize(nTrg, std::complex<double>(0, 0));

    if (numSources_ == 0 || nTrg == 0) return;

    // grad(Phi) = integral{ sigma * grad(1/r) } dA' / eps_0
    // gradient_P2P computes grad(sum{ q / r })

    std::vector<double> gradX_real(nTrg), gradY_real(nTrg), gradZ_real(nTrg);
    std::vector<double> gradX_imag(nTrg), gradY_imag(nTrg), gradZ_imag(nTrg);

    RunGradientFmm(sigma_real_, targets, gradX_real, gradY_real, gradZ_real);
    RunGradientFmm(sigma_imag_, targets, gradX_imag, gradY_imag, gradZ_imag);

    #pragma omp parallel for
    for (int i = 0; i < nTrg; ++i) {
        dPhidx[i] = std::complex<double>(gradX_real[i], gradX_imag[i]) / EPS_0;
        dPhidy[i] = std::complex<double>(gradY_real[i], gradY_imag[i]) / EPS_0;
        dPhidz[i] = std::complex<double>(gradZ_real[i], gradZ_imag[i]) / EPS_0;
    }
}

// ============================================================================
// Electric field E
// ============================================================================

void radTConductorFmm::ComputeE(const std::vector<double>& targets,
                                 double frequency,
                                 std::vector<std::complex<double>>& Ex,
                                 std::vector<std::complex<double>>& Ey,
                                 std::vector<std::complex<double>>& Ez) {
    int nTrg = static_cast<int>(targets.size() / 3);

    Ex.resize(nTrg, std::complex<double>(0, 0));
    Ey.resize(nTrg, std::complex<double>(0, 0));
    Ez.resize(nTrg, std::complex<double>(0, 0));

    if (numSources_ == 0 || nTrg == 0) return;

    // E = -grad(Phi) - j*omega*A

    // Compute -grad(Phi)
    std::vector<std::complex<double>> dPhidx, dPhidy, dPhidz;
    ComputeGradPhi(targets, dPhidx, dPhidy, dPhidz);

    // Compute A
    std::vector<std::complex<double>> Ax, Ay, Az;
    ComputeA(targets, Ax, Ay, Az);

    // Combine: E = -grad(Phi) - j*omega*A
    double omega = 2.0 * RadConst::PI * frequency;
    std::complex<double> jOmega(0, omega);

    #pragma omp parallel for
    for (int i = 0; i < nTrg; ++i) {
        Ex[i] = -dPhidx[i] - jOmega * Ax[i];
        Ey[i] = -dPhidy[i] - jOmega * Ay[i];
        Ez[i] = -dPhidz[i] - jOmega * Az[i];
    }
}

// ============================================================================
// Magnetic field B
// ============================================================================

void radTConductorFmm::ComputeB(const std::vector<double>& targets,
                                 std::vector<std::complex<double>>& Bx,
                                 std::vector<std::complex<double>>& By,
                                 std::vector<std::complex<double>>& Bz) {
    int nTrg = static_cast<int>(targets.size() / 3);

    Bx.resize(nTrg, std::complex<double>(0, 0));
    By.resize(nTrg, std::complex<double>(0, 0));
    Bz.resize(nTrg, std::complex<double>(0, 0));

    if (numSources_ == 0 || nTrg == 0) return;

    // Use direct method for small problems
    if (!enabled_ || numSources_ < threshold_) {
        ComputeBDirect(targets, Bx, By, Bz);
        return;
    }

    // B = mu_0/(4*pi) * integral{ K x r / |r|^3 } dA'
    //
    // Using the identity for Biot-Savart with gradient:
    // B_x = mu_0/(4*pi) * integral{ Ky * dz/r^3 - Kz * dy/r^3 } dA'
    // B_y = mu_0/(4*pi) * integral{ Kz * dx/r^3 - Kx * dz/r^3 } dA'
    // B_z = mu_0/(4*pi) * integral{ Kx * dy/r^3 - Ky * dx/r^3 } dA'
    //
    // Since gradient_P2P computes grad(1/r) = -r/r^3, we have:
    // B_x = -mu_0 * (Ky * grad_z(1/r) - Kz * grad_y(1/r))
    // B_y = -mu_0 * (Kz * grad_x(1/r) - Kx * grad_z(1/r))
    // B_z = -mu_0 * (Kx * grad_y(1/r) - Ky * grad_x(1/r))
    //
    // We need to compute gradient contributions from each K component

    // Allocate gradient results for each K component
    std::vector<double> gradKx_x_real(nTrg), gradKx_y_real(nTrg), gradKx_z_real(nTrg);
    std::vector<double> gradKx_x_imag(nTrg), gradKx_y_imag(nTrg), gradKx_z_imag(nTrg);
    std::vector<double> gradKy_x_real(nTrg), gradKy_y_real(nTrg), gradKy_z_real(nTrg);
    std::vector<double> gradKy_x_imag(nTrg), gradKy_y_imag(nTrg), gradKy_z_imag(nTrg);
    std::vector<double> gradKz_x_real(nTrg), gradKz_y_real(nTrg), gradKz_z_real(nTrg);
    std::vector<double> gradKz_x_imag(nTrg), gradKz_y_imag(nTrg), gradKz_z_imag(nTrg);

    // Compute gradients for each K component
    RunGradientFmm(Kx_real_, targets, gradKx_x_real, gradKx_y_real, gradKx_z_real);
    RunGradientFmm(Kx_imag_, targets, gradKx_x_imag, gradKx_y_imag, gradKx_z_imag);
    RunGradientFmm(Ky_real_, targets, gradKy_x_real, gradKy_y_real, gradKy_z_real);
    RunGradientFmm(Ky_imag_, targets, gradKy_x_imag, gradKy_y_imag, gradKy_z_imag);
    RunGradientFmm(Kz_real_, targets, gradKz_x_real, gradKz_y_real, gradKz_z_real);
    RunGradientFmm(Kz_imag_, targets, gradKz_x_imag, gradKz_y_imag, gradKz_z_imag);

    // Compute B = -mu_0 * (K x grad(1/r))
    // Note: gradient_P2P returns grad(sum{ q / r }) which includes 1/(4*pi)
    // So the result is already mu_0/(4*pi) * ... when multiplied by mu_0
    // But we need to be careful about the sign convention

    #pragma omp parallel for
    for (int i = 0; i < nTrg; ++i) {
        // Real parts
        double Bx_real = -(gradKy_z_real[i] - gradKz_y_real[i]) * MU_0;
        double By_real = -(gradKz_x_real[i] - gradKx_z_real[i]) * MU_0;
        double Bz_real = -(gradKx_y_real[i] - gradKy_x_real[i]) * MU_0;

        // Imaginary parts
        double Bx_imag = -(gradKy_z_imag[i] - gradKz_y_imag[i]) * MU_0;
        double By_imag = -(gradKz_x_imag[i] - gradKx_z_imag[i]) * MU_0;
        double Bz_imag = -(gradKx_y_imag[i] - gradKy_x_imag[i]) * MU_0;

        Bx[i] = std::complex<double>(Bx_real, Bx_imag);
        By[i] = std::complex<double>(By_real, By_imag);
        Bz[i] = std::complex<double>(Bz_real, Bz_imag);
    }
}

// ============================================================================
// Direct computation methods (fallback for small problems)
// ============================================================================

void radTConductorFmm::ComputeADirect(const std::vector<double>& targets,
                                       std::vector<std::complex<double>>& Ax,
                                       std::vector<std::complex<double>>& Ay,
                                       std::vector<std::complex<double>>& Az) {
    int nTrg = static_cast<int>(targets.size() / 3);
    int nSrc = numSources_;

    Ax.assign(nTrg, std::complex<double>(0, 0));
    Ay.assign(nTrg, std::complex<double>(0, 0));
    Az.assign(nTrg, std::complex<double>(0, 0));

    // A = mu_0/(4*pi) * integral{ K / |r - r'| } dA'
    double coeff = MU_0 * INV_FOUR_PI;

    #pragma omp parallel for
    for (int t = 0; t < nTrg; ++t) {
        double tx = targets[3*t + 0];
        double ty = targets[3*t + 1];
        double tz = targets[3*t + 2];

        std::complex<double> ax(0, 0), ay(0, 0), az(0, 0);

        for (int s = 0; s < nSrc; ++s) {
            double dx = tx - sourceCoords_[3*s + 0];
            double dy = ty - sourceCoords_[3*s + 1];
            double dz = tz - sourceCoords_[3*s + 2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (r < 1e-15) continue;

            double invR = 1.0 / r;
            double dA = sourceAreas_[s];

            std::complex<double> Kx(Kx_real_[s], Kx_imag_[s]);
            std::complex<double> Ky(Ky_real_[s], Ky_imag_[s]);
            std::complex<double> Kz(Kz_real_[s], Kz_imag_[s]);

            ax += Kx * invR * dA;
            ay += Ky * invR * dA;
            az += Kz * invR * dA;
        }

        Ax[t] = ax * coeff;
        Ay[t] = ay * coeff;
        Az[t] = az * coeff;
    }
}

void radTConductorFmm::ComputePhiDirect(const std::vector<double>& targets,
                                         std::vector<std::complex<double>>& phi) {
    int nTrg = static_cast<int>(targets.size() / 3);
    int nSrc = numSources_;

    phi.assign(nTrg, std::complex<double>(0, 0));

    // Phi = 1/(4*pi*eps_0) * integral{ sigma / |r - r'| } dA'
    double coeff = INV_FOUR_PI / EPS_0;

    #pragma omp parallel for
    for (int t = 0; t < nTrg; ++t) {
        double tx = targets[3*t + 0];
        double ty = targets[3*t + 1];
        double tz = targets[3*t + 2];

        std::complex<double> p(0, 0);

        for (int s = 0; s < nSrc; ++s) {
            double dx = tx - sourceCoords_[3*s + 0];
            double dy = ty - sourceCoords_[3*s + 1];
            double dz = tz - sourceCoords_[3*s + 2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (r < 1e-15) continue;

            double invR = 1.0 / r;
            double dA = sourceAreas_[s];

            std::complex<double> sigma(sigma_real_[s], sigma_imag_[s]);
            p += sigma * invR * dA;
        }

        phi[t] = p * coeff;
    }
}

void radTConductorFmm::ComputeBDirect(const std::vector<double>& targets,
                                       std::vector<std::complex<double>>& Bx,
                                       std::vector<std::complex<double>>& By,
                                       std::vector<std::complex<double>>& Bz) {
    int nTrg = static_cast<int>(targets.size() / 3);
    int nSrc = numSources_;

    Bx.assign(nTrg, std::complex<double>(0, 0));
    By.assign(nTrg, std::complex<double>(0, 0));
    Bz.assign(nTrg, std::complex<double>(0, 0));

    // B = mu_0/(4*pi) * integral{ K x r / |r|^3 } dA'
    double coeff = MU_0 * INV_FOUR_PI;

    #pragma omp parallel for
    for (int t = 0; t < nTrg; ++t) {
        double tx = targets[3*t + 0];
        double ty = targets[3*t + 1];
        double tz = targets[3*t + 2];

        std::complex<double> bx(0, 0), by(0, 0), bz(0, 0);

        for (int s = 0; s < nSrc; ++s) {
            double dx = tx - sourceCoords_[3*s + 0];
            double dy = ty - sourceCoords_[3*s + 1];
            double dz = tz - sourceCoords_[3*s + 2];
            double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 < 1e-30) continue;

            double r = std::sqrt(r2);
            double invR3 = 1.0 / (r2 * r);
            double dA = sourceAreas_[s];

            std::complex<double> Kx(Kx_real_[s], Kx_imag_[s]);
            std::complex<double> Ky(Ky_real_[s], Ky_imag_[s]);
            std::complex<double> Kz(Kz_real_[s], Kz_imag_[s]);

            // K x r
            bx += (Ky * dz - Kz * dy) * invR3 * dA;
            by += (Kz * dx - Kx * dz) * invR3 * dA;
            bz += (Kx * dy - Ky * dx) * invR3 * dA;
        }

        Bx[t] = bx * coeff;
        By[t] = by * coeff;
        Bz[t] = bz * coeff;
    }
}

// ============================================================================
// radTConductorFmmManager implementation (singleton)
// ============================================================================

radTConductorFmmManager& radTConductorFmmManager::Instance() {
    static radTConductorFmmManager instance;
    return instance;
}

radTConductorFmmManager::radTConductorFmmManager()
    : globalEnabled_(true)
    , p_(6)
    , ncrit_(64)
    , threshold_(500)
{
}

radTConductorFmmManager::~radTConductorFmmManager() {
}

radTConductorFmm& radTConductorFmmManager::GetFmm() {
    if (!fmm_) {
        fmm_ = std::make_unique<radTConductorFmm>();
        fmm_->SetParameters(p_, ncrit_);
        fmm_->SetThreshold(threshold_);
        fmm_->SetEnabled(globalEnabled_);
    }
    return *fmm_;
}

void radTConductorFmmManager::SetEnabled(bool enabled) {
    globalEnabled_ = enabled;
    if (fmm_) {
        fmm_->SetEnabled(enabled);
    }
}

bool radTConductorFmmManager::IsEnabled() const {
    return globalEnabled_;
}

void radTConductorFmmManager::SetParameters(int p, int ncrit) {
    p_ = p;
    ncrit_ = ncrit;
    if (fmm_) {
        fmm_->SetParameters(p, ncrit);
    }
}

void radTConductorFmmManager::SetThreshold(int threshold) {
    threshold_ = threshold;
    if (fmm_) {
        fmm_->SetThreshold(threshold);
    }
}

} // namespace radia
