/*
 * rad_peec_mmm_coupled.h
 *
 * Unified PEEC Loop-Star + MMM + MSC Coupled Solver
 *
 * This module provides coupled electromagnetic analysis combining:
 * - PEEC (Partial Element Equivalent Circuit) with Loop-Star decomposition
 * - MMM (Magnetic Moment Method) for magnetic materials
 * - MSC (Magnetic Surface Charge) for dielectric materials (future)
 * - ESIM (Effective Surface Impedance Method) for conductive materials
 *
 * Architecture (2026-01-10 redesign):
 * - Laplace kernel only (no Helmholtz) for quasi-static analysis
 * - Complete Loop-Star separation for low-frequency stability
 * - Frequency domain with complex material properties (mu'', epsilon'')
 *
 * Applications:
 * - WPT (Wireless Power Transfer) self-resonance analysis
 * - Coil impedance with magnetic core
 * - Future: Induction heating with ferromagnetic workpiece
 *
 * Unified System Matrix:
 *
 *   [Z_LL   Z_LS   Z_LM   0    ] [I_L ]   [V_L ]
 *   [Z_SL   Z_SS   0      Z_SE ] [I_S ] = [V_S ]
 *   [Z_ML   0      Z_MM   0    ] [M   ]   [H_ext]
 *   [0      Z_ES   0      Z_EE ] [P   ]   [D_ext]
 *
 *   I_L: Loop currents (solenoidal, div J = 0)
 *   I_S: Star currents (irrotational, charge-related)
 *   M:   Magnetization (MMM)
 *   P:   Polarization (MSC for dielectrics, future)
 *
 * Coupling Physics:
 *   Z_LL: Inductance via vector potential A (Laplace: mu0/(4*pi*r))
 *   Z_SS: Capacitance via scalar potential (Laplace: 1/(4*pi*epsilon*r))
 *   Z_LM: B field from M affects J (Loop-MMM coupling via -dA/dt)
 *   Z_ML: H field from J affects M (MMM-Loop coupling via Biot-Savart)
 *   Z_SE: E field from P affects sigma (Star-MSC for dielectric, future)
 *   Z_ES: E field from sigma affects P (MSC-Star for dielectric, future)
 *
 * ESIM Surface Impedance (Karl Hollaus formulation):
 *   Zs = (1+j) * Rs
 *   Rs = 1 / (sigma * delta)
 *   delta = sqrt(2 / (omega * mu0 * mu_r * sigma))
 *
 * Complex Material Properties:
 *   mu = mu0 * (mu_r' - j * mu_r'')  // Magnetic losses
 *   epsilon = epsilon0 * (epsilon_r' - j * epsilon_r'')  // Dielectric losses
 *
 * References:
 * [1] Z. Zhu et al., "Algorithms in FastImp", IEEE TCAD, 2005
 * [2] K. Hollaus, "Effective Surface Impedance...", 2024
 * [3] G. Vecchi, "Loop-Star decomposition", IEEE TAP, 1999
 *
 * Part of Radia project
 */

#ifndef RAD_PEEC_MMM_COUPLED_H
#define RAD_PEEC_MMM_COUPLED_H

#include "rad_conductor.h"
#include "rad_interaction.h"
#include "gmvect.h"
#include <complex>
#include <memory>
#include <vector>

namespace radia {

/**
 * @brief Result structure for coupled PEEC-MMM solver
 */
struct PEECMMMResult {
    std::complex<double> impedance;     // Port impedance [Ohm]
    double resistance;                   // Re(Z) [Ohm]
    double reactance;                    // Im(Z) [Ohm]
    double inductance;                   // Im(Z)/(2*pi*f) [H]
    double conductorPower;              // Ohmic loss in conductor [W]
    double magnetPower;                 // Eddy current loss in magnetic material [W] (if conductive)
    int conductorDOF;                   // Number of conductor DOF
    int magnetDOF;                      // Number of magnet DOF
    int iterations;                     // Number of nonlinear iterations (if applicable)
    double maxMagnetizationChange;      // Convergence indicator
};

/**
 * @brief Coupled PEEC-MMM solver
 *
 * Solves coupled electromagnetic system with:
 * - PEEC conductor model (surface currents J)
 * - MMM magnet model (magnetization M)
 */
class PEECMMMCoupledSolver {
public:
    PEECMMMCoupledSolver();
    ~PEECMMMCoupledSolver();

    // ========== Component setup ==========

    /**
     * @brief Set conductor element (PEEC)
     * @param conductor Pointer to conductor object
     */
    void SetConductor(std::shared_ptr<radTConductor> conductor);

    /**
     * @brief Set magnet element (MMM)
     * @param magnetHandle Radia handle to magnetic object
     *
     * The magnetic object should already have material applied via MatApl()
     */
    void SetMagnet(int magnetHandle);

    /**
     * @brief Set magnet element directly via Interaction object
     * @param interaction Pointer to radTInteraction object
     */
    void SetMagnetInteraction(radTInteraction* interaction);

    /**
     * @brief Get magnet handle
     */
    int GetMagnetHandle() const { return magnetHandle_; }

    /**
     * @brief Check if magnet interaction is set
     */
    bool HasMagnetInteraction() const { return magnetInteraction_ != nullptr; }

    /**
     * @brief Set ownership flag for magnet interaction
     */
    void SetOwnsMagnetInteraction(bool owns) { ownsMagnetInteraction_ = owns; }

    // ========== Analysis parameters ==========

    /**
     * @brief Set analysis frequency
     * @param frequency Frequency in Hz (0 for DC/static)
     */
    void SetFrequency(double frequency);
    double GetFrequency() const { return frequency_; }

    /**
     * @brief Set convergence tolerance for nonlinear iteration
     * @param tol Relative tolerance (default 1e-4)
     */
    void SetTolerance(double tol) { tolerance_ = tol; }

    /**
     * @brief Set maximum iterations for nonlinear solve
     * @param maxIter Maximum iterations (default 100)
     */
    void SetMaxIterations(int maxIter) { maxIterations_ = maxIter; }

    // ========== Matrix Symmetrization (for CLN model order reduction) ==========

    /**
     * @brief Enable/disable matrix symmetrization
     * @param useSymmetric If true, symmetrize the coupled matrix via variable scaling
     *
     * Symmetrization uses scaling: M' = sqrt(mu0 * V) * M
     * This makes Z_LM' = Z_ML'^T (symmetric) via reciprocity.
     *
     * Benefits:
     * - Enables CLN (Cauer Ladder Network) model order reduction
     * - Passivity guaranteed
     * - Physical interpretation as L-C ladder circuit
     *
     * The solution (I, M) is identical to non-symmetric formulation.
     */
    void SetUseSymmetric(bool useSymmetric) { useSymmetric_ = useSymmetric; }
    bool GetUseSymmetric() const { return useSymmetric_; }

    // ========== Complex Permeability (μ = μ' - jμ") ==========

    /**
     * @brief Set uniform complex relative permeability for all magnet elements
     * @param mu_r_real Real part of relative permeability (μ'_r)
     * @param mu_r_imag Imaginary part of relative permeability (μ"_r, positive for loss)
     *
     * The complex susceptibility is: χ = (μ'_r - 1) - jμ"_r
     * Magnetic loss power: P_mag = (ω/2) * μ0 * μ"_r * |H|^2 * Volume
     */
    void SetComplexPermeability(double mu_r_real, double mu_r_imag);

    /**
     * @brief Set complex relative permeability per element
     * @param mu_r_real Vector of real parts (μ'_r) for each element
     * @param mu_r_imag Vector of imaginary parts (μ"_r) for each element
     */
    void SetComplexPermeabilityPerElement(const std::vector<double>& mu_r_real,
                                          const std::vector<double>& mu_r_imag);

    /**
     * @brief Check if complex permeability is set
     */
    bool HasComplexPermeability() const { return hasComplexMu_; }

    // ========== Excitation ==========

    /**
     * @brief Set voltage excitation at conductor port
     * @param V Complex voltage [V]
     */
    void SetVoltageExcitation(std::complex<double> V);

    /**
     * @brief Set current excitation at conductor port
     * @param I Complex current [A]
     */
    void SetCurrentExcitation(std::complex<double> I);

    /**
     * @brief Set external magnetic field (background field)
     * @param Hx, Hy, Hz External H field components [A/m]
     */
    void SetExternalField(double Hx, double Hy, double Hz);

    // ========== Solve ==========

    /**
     * @brief Solve the coupled system
     * @return Result structure with impedance, power, etc.
     */
    PEECMMMResult Solve();

    /**
     * @brief Solve at multiple frequencies (frequency sweep)
     * @param frequencies Vector of frequencies [Hz]
     * @return Vector of results for each frequency
     */
    std::vector<PEECMMMResult> FrequencySweep(const std::vector<double>& frequencies);

    // ========== Results ==========

    /**
     * @brief Get port impedance after solve
     */
    std::complex<double> GetImpedance() const { return impedance_; }

    /**
     * @brief Get total current through conductor
     */
    std::complex<double> GetTotalCurrent() const { return totalCurrent_; }

    /**
     * @brief Get conductor power loss
     */
    double GetConductorPower() const { return conductorPower_; }

    /**
     * @brief Get magnet eddy current power loss (if conductive)
     */
    double GetMagnetPower() const { return magnetPower_; }

    // ========== Field computation ==========

    /**
     * @brief Compute total B field at a point (after solve)
     *
     * B_total = B_conductor + B_magnet
     *
     * @param point Evaluation point [m]
     * @param Bx, By, Bz Output B field components (complex for AC)
     */
    void ComputeB(const TVector3d& point,
                  std::complex<double>& Bx,
                  std::complex<double>& By,
                  std::complex<double>& Bz) const;

    /**
     * @brief Compute B field from conductor only
     */
    void ComputeBConductor(const TVector3d& point,
                           std::complex<double>& Bx,
                           std::complex<double>& By,
                           std::complex<double>& Bz) const;

    /**
     * @brief Compute B field from magnet only
     * @param point Evaluation point [m]
     * @param Bx, By, Bz Output B field components (real for static)
     */
    void ComputeBMagnet(const TVector3d& point,
                        double& Bx, double& By, double& Bz) const;

    /**
     * @brief Batch compute B field at multiple points (OpenMP parallelized)
     */
    void ComputeBBatch(const std::vector<TVector3d>& points,
                       std::vector<std::complex<double>>& Bx,
                       std::vector<std::complex<double>>& By,
                       std::vector<std::complex<double>>& Bz) const;

private:
    // Components
    std::shared_ptr<radTConductor> conductor_;
    radTInteraction* magnetInteraction_;
    int magnetHandle_;
    bool ownsMagnetInteraction_;  // True if we created the interaction

    // Analysis parameters
    double frequency_;
    double omega_;
    double tolerance_;
    int maxIterations_;
    bool useSymmetric_;  // Use symmetrized matrix formulation

    // Complex permeability: μ = μ' - jμ" → χ = (μ'_r - 1) - jμ"_r
    bool hasComplexMu_;                              // True if complex permeability is set
    bool uniformComplexMu_;                          // True if uniform (same for all elements)
    std::complex<double> uniformChi_;                // Uniform complex susceptibility
    std::vector<std::complex<double>> elementChi_;   // Per-element complex susceptibility

    // Excitation
    enum ExcitationType { None, Voltage, Current } excitationType_;
    std::complex<double> excitationValue_;
    TVector3d externalField_;

    // Coupling matrices
    // Z_cm[i][j] = B field from conductor panel j at magnet element i
    // Z_mc[i][j] = B field from magnet element j at conductor panel i
    std::vector<std::complex<double>> Z_cm_;  // Size: n_mag_elem * n_cond_panel
    std::vector<std::complex<double>> Z_mc_;  // Size: n_cond_panel * n_mag_elem

    // Solution
    std::complex<double> impedance_;
    std::complex<double> totalCurrent_;
    double conductorPower_;
    double magnetPower_;
    bool solved_;
    std::vector<std::complex<double>> solution_;  // Solution vector [I_L; M]

    // Internal methods
    void AssembleCouplingMatrices();
    void ComputeZcm();  // B_cond at magnet centers
    void ComputeZmc();  // B_mag at conductor panel centers

    void SolveLinear();     // Linear material case
    void SolveNonlinear();  // Nonlinear material case (iterate M)

    void ComputeConductorPower();
    void ComputeMagnetPower();

    // Helper: Compute B from conductor at a point
    void ConductorBField(const TVector3d& point,
                         std::complex<double>& Bx,
                         std::complex<double>& By,
                         std::complex<double>& Bz) const;

    // Helper: Compute B from magnet at a point (static)
    void MagnetBField(const TVector3d& point,
                      double& Bx, double& By, double& Bz) const;
};

// ========== Python API helper functions ==========

/**
 * @brief Create coupled solver from handles
 * @param conductorHandle Radia handle for conductor
 * @param magnetHandle Radia handle for magnet
 * @return Pointer to coupled solver
 */
PEECMMMCoupledSolver* CreatePEECMMMSolver(int conductorHandle, int magnetHandle);

/**
 * @brief Solve coupled system and return impedance
 * @param solver Pointer to solver
 * @param frequency Analysis frequency [Hz]
 * @return Complex impedance [Ohm]
 */
std::complex<double> SolvePEECMMM(PEECMMMCoupledSolver* solver, double frequency);

/**
 * @brief Get B field from coupled system
 * @param solver Pointer to solver
 * @param point Evaluation point
 * @return B field as [Bx_re, By_re, Bz_re, Bx_im, By_im, Bz_im]
 */
std::vector<double> GetPEECMMMField(PEECMMMCoupledSolver* solver,
                                     const TVector3d& point);

} // namespace radia

#endif // RAD_PEEC_MMM_COUPLED_H
