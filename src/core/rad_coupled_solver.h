/*
 * rad_coupled_solver.h
 *
 * Coupled solver for conductor + magnetic material analysis
 * Handles cross-terms between FastImp conductor formulation and Radia MSC
 *
 * References:
 * [1] Z. Zhu et al., "Algorithms in FastImp", IEEE TCAD, 2005
 * [2] S. Bilicz et al., "Nonlocal SIBC", ISEM 2023
 *
 * Architecture:
 * [Z_cc  Z_cm] [K    ]   [V_ext]
 * [Z_mc  Z_mm] [sigma] = [H_ext]
 *
 * Z_cc: Conductor-Conductor (FastImp)
 * Z_mm: Magnetic-Magnetic (Radia MSC)
 * Z_cm, Z_mc: Cross-coupling terms (this module)
 *
 * Part of Radia project
 */

#ifndef RAD_COUPLED_SOLVER_H
#define RAD_COUPLED_SOLVER_H

#include "rad_conductor.h"
#include "rad_green_fullwave.h"
#include "gmvect.h"
#include "rad_interaction.h"
#include <vector>
#include <complex>
#include <memory>
#include <functional>

namespace radia {

// Forward declaration
class radTg3dRelax;

/**
 * @brief Cross-term type for coupled system
 */
enum class CrossTermType {
    ConductorToMagnetic,  // Z_cm: Field from conductor current affects magnetization
    MagneticToConductor   // Z_mc: Field from magnetization affects conductor current
};

/**
 * @brief Cross-term matrix block between conductor and magnetic elements
 *
 * Computes interaction terms between conductor surface currents (K) and
 * magnetic material surface charges (sigma_m) or magnetization (M).
 */
class radTCrossTerms {
public:
    using Complex = std::complex<double>;

    radTCrossTerms();
    ~radTCrossTerms();

    /**
     * @brief Set conductor geometry
     * @param conductor Pointer to conductor solver
     */
    void SetConductor(radTConductorSolver* conductor);

    /**
     * @brief Set magnetic material interaction object
     * @param interaction Pointer to Radia interaction object
     */
    void SetMagneticInteraction(radTInteraction* interaction);

    /**
     * @brief Set operating frequency
     * @param frequency Frequency in Hz
     */
    void SetFrequency(double frequency);

    /**
     * @brief Compute Z_cm matrix block
     *
     * Field from conductor currents K at magnetic material evaluation points
     * H_cm = integral{ K x grad(G) } dS_conductor
     *
     * This field contributes to the RHS of the magnetic material equation.
     *
     * @param Z_cm Output matrix (size: n_magn x n_cond)
     */
    void ComputeConductorToMagnetic(std::vector<Complex>& Z_cm);

    /**
     * @brief Compute Z_mc matrix block
     *
     * Field from magnetization M at conductor evaluation points
     * H_mc = -grad(Phi_m) where Phi_m = (1/4pi) * integral{ sigma_m / r } dS_magn
     *
     * This field contributes to the conductor EFIE.
     *
     * @param Z_mc Output matrix (size: n_cond x n_magn)
     */
    void ComputeMagneticToConductor(std::vector<Complex>& Z_mc);

    /**
     * @brief Get number of conductor DOFs
     */
    int NumConductorDOF() const { return nCondDOF_; }

    /**
     * @brief Get number of magnetic material DOFs
     */
    int NumMagneticDOF() const { return nMagnDOF_; }

private:
    radTConductorSolver* conductor_;
    radTInteraction* interaction_;
    double frequency_;

    int nCondDOF_;      // Conductor DOFs (2 * n_panels for K + sigma)
    int nMagnDOF_;      // Magnetic DOFs (3 * n_elements for Mx, My, Mz or 6 for MSC)

    // Green's function calculator
    radTGreenFunction green_;

    // Cached magnetic element geometry
    std::vector<TVector3d> magnCenters_;   // Magnetic element centers
    std::vector<TVector3d> magnNormals_;   // Magnetic element normals (for MSC)
    std::vector<double> magnAreas_;         // Magnetic element areas

    void CacheMagneticGeometry();
};

/**
 * @brief Coupled solver for conductor + magnetic material systems
 *
 * Solves the coupled system:
 * [Z_cc  Z_cm] [x_c]   [b_c]
 * [Z_mc  Z_mm] [x_m] = [b_m]
 *
 * Uses block elimination or iterative coupling depending on problem size.
 */
class radTCoupledSolver {
public:
    using Complex = std::complex<double>;

    radTCoupledSolver();
    ~radTCoupledSolver();

    /**
     * @brief Add conductor solver
     */
    void SetConductorSolver(std::shared_ptr<radTConductorSolver> conductor);

    /**
     * @brief Set Radia magnetic interaction
     */
    void SetMagneticInteraction(radTInteraction* interaction);

    /**
     * @brief Set operating frequency
     */
    void SetFrequency(double frequency);

    /**
     * @brief Set formulation
     */
    void SetFormulation(ConductorFormulation form);

    /**
     * @brief Solve coupled system at current frequency
     */
    void Solve();

    /**
     * @brief Solve and return total impedance
     */
    Complex SolveImpedance(double frequency);

    /**
     * @brief Frequency sweep for impedance
     */
    std::vector<Complex> ImpedanceSweep(const std::vector<double>& frequencies);

    /**
     * @brief Enable/disable iterative coupling
     *
     * If disabled, uses direct block LU factorization.
     * If enabled, solves conductor and magnetic blocks iteratively.
     *
     * Default: false (direct solve)
     */
    void SetIterativeCoupling(bool enable) { iterativeCoupling_ = enable; }

    /**
     * @brief Set convergence tolerance for iterative coupling
     */
    void SetCouplingTolerance(double tol) { couplingTol_ = tol; }

    /**
     * @brief Set max iterations for iterative coupling
     */
    void SetMaxCouplingIterations(int maxIter) { maxCouplingIter_ = maxIter; }

    /**
     * @brief Get solution for conductor subsystem
     */
    const std::vector<Complex>& GetConductorSolution() const { return solCond_; }

    /**
     * @brief Get solution for magnetic subsystem
     */
    const std::vector<double>& GetMagneticSolution() const { return solMagn_; }

private:
    // Subsystem solvers
    std::shared_ptr<radTConductorSolver> conductor_;
    radTInteraction* interaction_;

    // Problem parameters
    double frequency_;
    ConductorFormulation formulation_;

    // Solution vectors
    std::vector<Complex> solCond_;    // Conductor solution (K, sigma)
    std::vector<double> solMagn_;      // Magnetic solution (M)

    // Cross-term calculator
    radTCrossTerms crossTerms_;

    // Coupling matrices
    std::vector<Complex> Z_cm_;       // Conductor -> Magnetic coupling
    std::vector<Complex> Z_mc_;       // Magnetic -> Conductor coupling

    // Full system matrix (if using direct solve)
    std::vector<Complex> fullSystemMatrix_;
    std::vector<Complex> fullRHS_;
    std::vector<Complex> fullSolution_;

    // Iterative coupling parameters
    bool iterativeCoupling_;
    double couplingTol_;
    int maxCouplingIter_;

    // Private methods
    void BuildCoupledSystem();
    void SolveDirectBlock();
    void SolveIterativeCoupling();

    void ExtractConductorSolution();
    void ExtractMagneticSolution();
};

/**
 * @brief Unified solver factory
 *
 * Creates appropriate solver based on problem type:
 * - Conductor only: FastImp solver
 * - Magnetic only: Radia MSC solver
 * - Coupled: Coupled solver
 */
class radTSolverFactory {
public:
    enum class SolverType {
        ConductorOnly,     // Pure conductor (FastImp)
        MagneticOnly,      // Pure magnetic (MSC)
        CoupledCM,         // Conductor + Magnetic (no SIBC)
        CoupledCMS         // Conductor + Magnetic + SIBC (Phase 3)
    };

    /**
     * @brief Create appropriate solver for problem
     * @param type Solver type
     * @return Unique pointer to solver
     */
    static std::unique_ptr<radTCoupledSolver> Create(SolverType type);

    /**
     * @brief Auto-detect solver type from problem contents
     * @param hasConductors True if problem has conductor elements
     * @param hasMagnetic True if problem has magnetic materials
     * @param hasSIBC True if problem has conductive magnetic materials
     * @return Detected solver type
     */
    static SolverType DetectType(bool hasConductors, bool hasMagnetic, bool hasSIBC);
};

} // namespace radia

#endif // RAD_COUPLED_SOLVER_H
