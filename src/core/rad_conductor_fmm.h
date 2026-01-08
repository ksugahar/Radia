/*
 * rad_conductor_fmm.h
 *
 * ExaFMM-accelerated field computation for FastImp conductors
 *
 * Uses LaplaceFmm kernel for all field types:
 * - A (vector potential): LaplaceFmm::potential_P2P
 * - Phi (scalar potential): LaplaceFmm::potential_P2P
 * - E (electric field): -grad(Phi) - j*omega*A (gradient_P2P + potential_P2P)
 * - B (magnetic field): K x grad(1/r) (gradient_P2P with cross product)
 *
 * References:
 * - ExaFMM-t: https://github.com/exafmm/exafmm-t
 * - FastImp: Z. Zhu et al., IEEE TCAD, 2005
 *
 * Part of Radia project
 */

#ifndef RAD_CONDUCTOR_FMM_H
#define RAD_CONDUCTOR_FMM_H

#include "rad_conductor.h"
#include <vector>
#include <complex>
#include <memory>

// Forward declarations for ExaFMM-t
namespace exafmm_t {
    class LaplaceFmm;
    template <typename T> struct Body;
    template <typename T> struct Node;
}

namespace radia {

/**
 * @brief FMM-accelerated field evaluator for conductors
 *
 * This class provides O(N) or O(N log N) field evaluation
 * for large-scale conductor problems using ExaFMM-t.
 */
class radTConductorFmm {
public:
    radTConductorFmm();
    ~radTConductorFmm();

    /**
     * @brief Set FMM parameters
     * @param p Expansion order (accuracy, default 6)
     * @param ncrit Number of particles per leaf (default 64)
     */
    void SetParameters(int p = 6, int ncrit = 64);

    /**
     * @brief Set source data from conductor panels
     * @param conductor Conductor with solved surface currents/charges
     */
    void SetSources(const radTConductor& conductor);

    /**
     * @brief Set multiple conductors as sources
     * @param conductors Vector of conductors
     */
    void SetSources(const std::vector<std::shared_ptr<radTConductor>>& conductors);

    /**
     * @brief Compute vector potential A at target points
     *
     * A = mu_0/(4*pi) * integral{ K / |r - r'| } dA'
     *
     * Uses LaplaceFmm::potential_P2P for each component of K
     *
     * @param targets Target point coordinates (3*nTargets values)
     * @param Ax, Ay, Az Output vector potential components (complex)
     */
    void ComputeA(const std::vector<double>& targets,
                  std::vector<std::complex<double>>& Ax,
                  std::vector<std::complex<double>>& Ay,
                  std::vector<std::complex<double>>& Az);

    /**
     * @brief Compute scalar potential Phi at target points
     *
     * Phi = 1/(4*pi*eps_0) * integral{ sigma / |r - r'| } dA'
     *
     * Uses LaplaceFmm::potential_P2P
     *
     * @param targets Target point coordinates
     * @param phi Output scalar potential (complex)
     */
    void ComputePhi(const std::vector<double>& targets,
                    std::vector<std::complex<double>>& phi);

    /**
     * @brief Compute gradient of scalar potential
     *
     * grad(Phi) = integral{ sigma * grad(1/r) } dA' / eps_0
     *
     * Uses LaplaceFmm::gradient_P2P
     *
     * @param targets Target point coordinates
     * @param dPhidx, dPhidy, dPhidz Output gradient components
     */
    void ComputeGradPhi(const std::vector<double>& targets,
                        std::vector<std::complex<double>>& dPhidx,
                        std::vector<std::complex<double>>& dPhidy,
                        std::vector<std::complex<double>>& dPhidz);

    /**
     * @brief Compute electric field E at target points
     *
     * E = -grad(Phi) - j*omega*A
     *
     * Uses gradient_P2P for -grad(Phi) and potential_P2P for A
     *
     * @param targets Target point coordinates
     * @param frequency Frequency in Hz (for j*omega*A term)
     * @param Ex, Ey, Ez Output electric field components (complex)
     */
    void ComputeE(const std::vector<double>& targets,
                  double frequency,
                  std::vector<std::complex<double>>& Ex,
                  std::vector<std::complex<double>>& Ey,
                  std::vector<std::complex<double>>& Ez);

    /**
     * @brief Compute magnetic field B at target points
     *
     * B = mu_0/(4*pi) * integral{ K x r / |r|^3 } dA'
     *   = mu_0 * curl(A)
     *
     * Using vector identity: K x grad(1/r) gives the curl contribution
     * B_i = epsilon_ijk * K_j * grad_k(1/r)
     *
     * Uses LaplaceFmm::gradient_P2P with cross product
     *
     * @param targets Target point coordinates
     * @param Bx, By, Bz Output magnetic field components (complex)
     */
    void ComputeB(const std::vector<double>& targets,
                  std::vector<std::complex<double>>& Bx,
                  std::vector<std::complex<double>>& By,
                  std::vector<std::complex<double>>& Bz);

    /**
     * @brief Check if FMM is available and initialized
     */
    bool IsReady() const { return isInitialized_; }

    /**
     * @brief Get number of source panels
     */
    int NumSources() const { return numSources_; }

    /**
     * @brief Enable/disable FMM (for comparison with direct method)
     */
    void SetEnabled(bool enabled) { enabled_ = enabled; }
    bool IsEnabled() const { return enabled_; }

    /**
     * @brief Get threshold for using FMM (use direct for small problems)
     */
    int GetThreshold() const { return threshold_; }
    void SetThreshold(int thresh) { threshold_ = thresh; }

private:
    // FMM parameters
    int p_;           // Expansion order
    int ncrit_;       // Particles per leaf
    bool enabled_;    // Enable FMM
    int threshold_;   // Minimum problem size for FMM

    // State
    bool isInitialized_;
    int numSources_;

    // Source data (panel centers and charges/currents)
    std::vector<double> sourceCoords_;   // 3*N coordinates
    std::vector<double> sourceAreas_;    // N panel areas

    // Surface current components (real and imaginary parts)
    std::vector<double> Kx_real_, Kx_imag_;
    std::vector<double> Ky_real_, Ky_imag_;
    std::vector<double> Kz_real_, Kz_imag_;

    // Surface charge (real and imaginary parts)
    std::vector<double> sigma_real_, sigma_imag_;

    // Physical parameters
    double conductivity_;

    // ExaFMM instance (opaque pointer to avoid header dependency)
    struct FmmImpl;
    std::unique_ptr<FmmImpl> fmm_;

    // Internal methods
    void InitializeFmm();
    void RunPotentialFmm(const std::vector<double>& srcValues,
                         const std::vector<double>& targets,
                         std::vector<double>& results);
    void RunGradientFmm(const std::vector<double>& srcValues,
                        const std::vector<double>& targets,
                        std::vector<double>& gradX,
                        std::vector<double>& gradY,
                        std::vector<double>& gradZ);

    // Direct computation fallback (for small problems or debugging)
    void ComputeADirect(const std::vector<double>& targets,
                        std::vector<std::complex<double>>& Ax,
                        std::vector<std::complex<double>>& Ay,
                        std::vector<std::complex<double>>& Az);
    void ComputePhiDirect(const std::vector<double>& targets,
                          std::vector<std::complex<double>>& phi);
    void ComputeBDirect(const std::vector<double>& targets,
                        std::vector<std::complex<double>>& Bx,
                        std::vector<std::complex<double>>& By,
                        std::vector<std::complex<double>>& Bz);
};

/**
 * @brief Global FMM accelerator instance for conductor field computation
 *
 * Provides singleton-like access for Python API integration
 */
class radTConductorFmmManager {
public:
    static radTConductorFmmManager& Instance();

    /**
     * @brief Get or create FMM accelerator for a conductor
     */
    radTConductorFmm& GetFmm();

    /**
     * @brief Enable/disable FMM acceleration globally
     */
    void SetEnabled(bool enabled);
    bool IsEnabled() const;

    /**
     * @brief Set FMM parameters
     */
    void SetParameters(int p, int ncrit);

    /**
     * @brief Set threshold for FMM usage
     */
    void SetThreshold(int threshold);

private:
    radTConductorFmmManager();
    ~radTConductorFmmManager();
    radTConductorFmmManager(const radTConductorFmmManager&) = delete;
    radTConductorFmmManager& operator=(const radTConductorFmmManager&) = delete;

    std::unique_ptr<radTConductorFmm> fmm_;
    bool globalEnabled_;
    int p_;
    int ncrit_;
    int threshold_;
};

} // namespace radia

#endif // RAD_CONDUCTOR_FMM_H
