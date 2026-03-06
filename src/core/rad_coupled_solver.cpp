/*
 * rad_coupled_solver.cpp
 *
 * Coupled solver for conductor + magnetic material analysis
 *
 * Part of Radia project
 */

#include "rad_coupled_solver.h"
#include "rad_interaction.h"
#include "rad_constants.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace radia {

// ============================================================================
// radTCrossTerms implementation
// ============================================================================

radTCrossTerms::radTCrossTerms()
    : conductor_(nullptr)
    , interaction_(nullptr)
    , frequency_(0)
    , nCondDOF_(0)
    , nMagnDOF_(0)
{
}

radTCrossTerms::~radTCrossTerms() {
}

void radTCrossTerms::SetConductor(radTConductorSolver* conductor) {
    conductor_ = conductor;

    if (conductor_) {
        // Count total conductor DOFs
        nCondDOF_ = 0;
        // Each panel has K (scalar) + sigma (scalar) = 2 DOF
        // This matches BuildSystemMatrix in rad_conductor.cpp
    }
}

void radTCrossTerms::SetMagneticInteraction(radTInteraction* interaction) {
    interaction_ = interaction;

    if (interaction_) {
        nMagnDOF_ = interaction_->GetTotalDOF();
        CacheMagneticGeometry();
    }
}

void radTCrossTerms::SetFrequency(double frequency) {
    frequency_ = frequency;
    green_.SetFrequency(frequency);
}

void radTCrossTerms::CacheMagneticGeometry() {
    if (!interaction_) return;

    // Extract element centers and geometry from interaction object
    // This requires access to the underlying g3dRelax elements

    // Note: The actual implementation would access:
    // interaction_->g3dRelaxPtrVect to get element pointers
    // Each element has center coordinates and geometry

    // For now, this is a placeholder that will be filled when
    // we have proper access to the Radia element geometry
}

void radTCrossTerms::ComputeConductorToMagnetic(std::vector<Complex>& Z_cm) {
    if (!conductor_ || !interaction_) {
        Z_cm.clear();
        return;
    }

    // Z_cm computes how conductor currents K affect magnetic material H field
    //
    // Physical model:
    // H_from_conductor = mu_0/(4*pi) * integral{ K x r_hat / r^2 } dS_conductor
    //
    // This H field appears in the RHS of the magnetic material equation:
    // M = chi * (H_ext + H_interaction + H_from_conductor)

    int nCond = nCondDOF_;
    int nMagn = nMagnDOF_;

    Z_cm.resize(nMagn * nCond, Complex(0, 0));

    // Physical parameters
    double omega = 2.0 * RadConst::PI * frequency_;

    // For each magnetic element evaluation point...
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nMagn; ++i) {
        // Get magnetic element evaluation point
        // For 3DOF elements: center of element
        // For 6DOF MSC elements: face centers

        TVector3d evalPoint;
        // TODO: Get eval point from interaction_->m_hexaEvalPoints or element center

        // For each conductor panel...
        // Compute H contribution from panel current K

        // For now, placeholder computation
        // Real implementation needs conductor panel geometry
    }

    // Note: The matrix Z_cm has structure:
    // Z_cm[i * nCond + j] = dH_i / dK_j
    // where H_i is the H field at magnetic element i
    // and K_j is the current at conductor panel j
}

void radTCrossTerms::ComputeMagneticToConductor(std::vector<Complex>& Z_mc) {
    if (!conductor_ || !interaction_) {
        Z_mc.clear();
        return;
    }

    // Z_mc computes how magnetic material affects conductor E field
    //
    // Physical model for magnetic materials:
    // The magnetization M creates a scalar potential:
    // Phi_m = (1/4*pi) * integral{ sigma_m / r } dS_magnetic
    // where sigma_m = M . n_hat (magnetic surface charge)
    //
    // This contributes to the H field at conductor:
    // H_from_magnetic = -grad(Phi_m)
    //
    // And affects the conductor EFIE through B = mu_0 * (H + M)

    int nCond = nCondDOF_;
    int nMagn = nMagnDOF_;

    Z_mc.resize(nCond * nMagn, Complex(0, 0));

    // For each conductor panel evaluation point...
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nCond; ++i) {
        // Get conductor panel center
        // TODO: Get from conductor panels

        TVector3d evalPoint;

        // For each magnetic element...
        // Compute H contribution from magnetization

        // For now, placeholder computation
    }

    // Note: The matrix Z_mc has structure:
    // Z_mc[i * nMagn + j] = dE_i / dM_j (or d/d_sigma_m for MSC)
    // where E_i is the E field condition at conductor panel i
    // and M_j is the magnetization at magnetic element j
}

// ============================================================================
// radTCoupledSolver implementation
// ============================================================================

radTCoupledSolver::radTCoupledSolver()
    : conductor_(nullptr)
    , interaction_(nullptr)
    , frequency_(0)
    , formulation_(ConductorFormulation::MQS)
    , iterativeCoupling_(false)
    , couplingTol_(1e-6)
    , maxCouplingIter_(100)
{
}

radTCoupledSolver::~radTCoupledSolver() {
}

void radTCoupledSolver::SetConductorSolver(std::shared_ptr<radTConductorSolver> conductor) {
    conductor_ = conductor;
    crossTerms_.SetConductor(conductor.get());
}

void radTCoupledSolver::SetMagneticInteraction(radTInteraction* interaction) {
    interaction_ = interaction;
    crossTerms_.SetMagneticInteraction(interaction);
}

void radTCoupledSolver::SetFrequency(double frequency) {
    frequency_ = frequency;

    if (conductor_) {
        conductor_->SetFrequency(frequency);
    }
    crossTerms_.SetFrequency(frequency);
}

void radTCoupledSolver::SetFormulation(ConductorFormulation form) {
    formulation_ = form;

    if (conductor_) {
        conductor_->SetFormulation(form);
    }
}

void radTCoupledSolver::Solve() {
    BuildCoupledSystem();

    if (iterativeCoupling_) {
        SolveIterativeCoupling();
    } else {
        SolveDirectBlock();
    }

    ExtractConductorSolution();
    ExtractMagneticSolution();
}

radTCoupledSolver::Complex radTCoupledSolver::SolveImpedance(double frequency) {
    SetFrequency(frequency);
    Solve();

    // Compute total impedance from conductor solution
    if (conductor_) {
        return conductor_->SolveImpedance(frequency);
    }
    return Complex(0, 0);
}

std::vector<radTCoupledSolver::Complex> radTCoupledSolver::ImpedanceSweep(
    const std::vector<double>& frequencies) {

    std::vector<Complex> result(frequencies.size());

    for (size_t i = 0; i < frequencies.size(); ++i) {
        result[i] = SolveImpedance(frequencies[i]);
    }

    return result;
}

void radTCoupledSolver::BuildCoupledSystem() {
    // Get subsystem sizes
    int nCond = crossTerms_.NumConductorDOF();
    int nMagn = crossTerms_.NumMagneticDOF();
    int nTotal = nCond + nMagn;

    if (nTotal == 0) return;

    // Build cross-term matrices
    crossTerms_.ComputeConductorToMagnetic(Z_cm_);
    crossTerms_.ComputeMagneticToConductor(Z_mc_);

    // Build full system matrix for direct solve
    if (!iterativeCoupling_) {
        fullSystemMatrix_.resize(nTotal * nTotal, Complex(0, 0));
        fullRHS_.resize(nTotal, Complex(0, 0));
        fullSolution_.resize(nTotal, Complex(0, 0));

        // Block structure:
        // [Z_cc | Z_cm] [x_c]   [b_c]
        // [-----+-----] [---] = [---]
        // [Z_mc | Z_mm] [x_m]   [b_m]

        // Fill Z_cc block (top-left, nCond x nCond)
        // This is obtained from conductor_->systemMatrix_

        // Fill Z_mm block (bottom-right, nMagn x nMagn)
        // This is obtained from interaction_->m_flatInteractMatrix

        // Fill Z_cm block (top-right, nCond x nMagn)
        for (int i = 0; i < nCond; ++i) {
            for (int j = 0; j < nMagn; ++j) {
                fullSystemMatrix_[(j + nCond) * nTotal + i] = Z_cm_[i * nMagn + j];
            }
        }

        // Fill Z_mc block (bottom-left, nMagn x nCond)
        for (int i = 0; i < nMagn; ++i) {
            for (int j = 0; j < nCond; ++j) {
                fullSystemMatrix_[j * nTotal + (i + nCond)] = Z_mc_[i * nCond + j];
            }
        }

        // Fill RHS
        // b_c from conductor excitation
        // b_m from magnetic external field
    }
}

void radTCoupledSolver::SolveDirectBlock() {
    int nTotal = static_cast<int>(fullRHS_.size());
    if (nTotal == 0) return;

    // Copy for LAPACK
    std::vector<Complex> A = fullSystemMatrix_;
    fullSolution_ = fullRHS_;
    std::vector<int> ipiv(nTotal);

    // Solve using LAPACK zgesv
    #ifdef USE_MKL
    lapack_int info = LAPACKE_zgesv(LAPACK_COL_MAJOR, nTotal, 1,
                                     reinterpret_cast<lapack_complex_double*>(A.data()),
                                     nTotal,
                                     ipiv.data(),
                                     reinterpret_cast<lapack_complex_double*>(fullSolution_.data()),
                                     nTotal);
    if (info != 0) {
        throw std::runtime_error("Coupled solver: LAPACK zgesv failed");
    }
    #else
    // Fallback Gaussian elimination (same as rad_conductor.cpp)
    for (int k = 0; k < nTotal; ++k) {
        int maxRow = k;
        double maxVal = std::abs(A[k * nTotal + k]);
        for (int i = k + 1; i < nTotal; ++i) {
            double val = std::abs(A[i * nTotal + k]);
            if (val > maxVal) {
                maxVal = val;
                maxRow = i;
            }
        }

        if (maxRow != k) {
            for (int j = 0; j < nTotal; ++j) {
                std::swap(A[k * nTotal + j], A[maxRow * nTotal + j]);
            }
            std::swap(fullSolution_[k], fullSolution_[maxRow]);
        }

        for (int i = k + 1; i < nTotal; ++i) {
            Complex factor = A[i * nTotal + k] / A[k * nTotal + k];
            for (int j = k; j < nTotal; ++j) {
                A[i * nTotal + j] -= factor * A[k * nTotal + j];
            }
            fullSolution_[i] -= factor * fullSolution_[k];
        }
    }

    for (int i = nTotal - 1; i >= 0; --i) {
        for (int j = i + 1; j < nTotal; ++j) {
            fullSolution_[i] -= A[i * nTotal + j] * fullSolution_[j];
        }
        fullSolution_[i] /= A[i * nTotal + i];
    }
    #endif
}

void radTCoupledSolver::SolveIterativeCoupling() {
    // Iterative coupling: alternating between conductor and magnetic solves
    //
    // Algorithm:
    // 1. Initialize: solve conductor alone with b_c
    // 2. Compute H_from_conductor at magnetic elements
    // 3. Solve magnetic: M = chi * (H_ext + H_interact + H_from_conductor)
    // 4. Compute H_from_magnetic at conductor elements
    // 5. Update conductor RHS: b_c_new = b_c + Z_cm * M
    // 6. Solve conductor with updated RHS
    // 7. Check convergence, repeat 2-6 if needed

    int nCond = crossTerms_.NumConductorDOF();
    int nMagn = crossTerms_.NumMagneticDOF();

    // Initialize solution vectors
    solCond_.resize(nCond, Complex(0, 0));
    solMagn_.resize(nMagn, 0.0);

    std::vector<Complex> solCond_old(nCond);
    std::vector<double> solMagn_old(nMagn);

    for (int iter = 0; iter < maxCouplingIter_; ++iter) {
        // Store old solutions
        solCond_old = solCond_;
        solMagn_old = solMagn_;

        // Step 1: Solve conductor subsystem
        if (conductor_) {
            conductor_->Solve();
            // Extract solution from conductor
        }

        // Step 2: Compute coupling to magnetic
        // H_coupling = Z_cm^T * K (approximate)

        // Step 3: Solve magnetic subsystem
        // This would call the existing Radia relaxation solver
        // with modified external field

        // Step 4: Compute coupling to conductor
        // E_coupling = Z_mc * M

        // Step 5: Update conductor RHS

        // Step 6: Check convergence
        double maxChange = 0;
        for (int i = 0; i < nCond; ++i) {
            double change = std::abs(solCond_[i] - solCond_old[i]);
            maxChange = std::max(maxChange, change);
        }
        for (int i = 0; i < nMagn; ++i) {
            double change = std::abs(solMagn_[i] - solMagn_old[i]);
            maxChange = std::max(maxChange, change);
        }

        if (maxChange < couplingTol_) {
            break;
        }
    }
}

void radTCoupledSolver::ExtractConductorSolution() {
    int nCond = crossTerms_.NumConductorDOF();
    solCond_.resize(nCond);

    if (!iterativeCoupling_ && !fullSolution_.empty()) {
        // Extract from full solution
        for (int i = 0; i < nCond; ++i) {
            solCond_[i] = fullSolution_[i];
        }
    }
}

void radTCoupledSolver::ExtractMagneticSolution() {
    int nCond = crossTerms_.NumConductorDOF();
    int nMagn = crossTerms_.NumMagneticDOF();
    solMagn_.resize(nMagn);

    if (!iterativeCoupling_ && !fullSolution_.empty()) {
        // Extract from full solution (magnetic part is complex -> real)
        for (int i = 0; i < nMagn; ++i) {
            solMagn_[i] = fullSolution_[nCond + i].real();
        }
    }
}

// ============================================================================
// radTSolverFactory implementation
// ============================================================================

std::unique_ptr<radTCoupledSolver> radTSolverFactory::Create(SolverType type) {
    auto solver = std::make_unique<radTCoupledSolver>();

    switch (type) {
        case SolverType::ConductorOnly:
            // Configure for conductor-only
            break;

        case SolverType::MagneticOnly:
            // Configure for magnetic-only
            break;

        case SolverType::CoupledCM:
            // Default coupled configuration
            break;

        case SolverType::CoupledCMS:
            // With SIBC (Phase 3)
            break;
    }

    return solver;
}

radTSolverFactory::SolverType radTSolverFactory::DetectType(
    bool hasConductors, bool hasMagnetic, bool hasSIBC) {

    if (hasSIBC) {
        return SolverType::CoupledCMS;
    } else if (hasConductors && hasMagnetic) {
        return SolverType::CoupledCM;
    } else if (hasConductors) {
        return SolverType::ConductorOnly;
    } else {
        return SolverType::MagneticOnly;
    }
}

} // namespace radia
