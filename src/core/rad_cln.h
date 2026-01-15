/**
 * @file rad_cln.h
 * @brief CLN (Cauer Ladder Network) model order reduction
 *
 * Lanczos algorithm for reducing PEEC (R, L) matrices to CLN I form.
 *
 * CLN I form:
 *   - R_diag: diagonal resistance matrix
 *   - L_tridiag: tridiagonal inductance matrix
 *   - Z(s) = (R_diag + s*L_tridiag)^-1
 *
 * Usage:
 *   lanczos(K=R, N=L) -> R_diag (diagonal), L_tridiag (tridiagonal)
 */

#pragma once

#include <vector>
#include <complex>

namespace radia {
namespace cln {

/**
 * @brief Result of Lanczos transformation for CLN I model reduction
 */
struct LanczosResult {
    // CLN parameters
    std::vector<double> L;           // Series inductance (alpha), size n_output
    std::vector<double> R;           // Parallel resistance (beta), size n_output

    // Transformation matrix Q (row-major, n_input x n_output)
    // Q is K-orthogonal: Q^T * K * Q = R_diag (diagonal)
    // Note: Q^T * Q != I (not orthogonal in standard sense)
    std::vector<double> Q;

    // Transformed matrices
    std::vector<double> R_diag;      // Diagonalized R matrix (n_output x n_output)
    std::vector<double> L_tridiag;   // Tridiagonalized L matrix (n_output x n_output)

    // Dimensions
    int n_input;                     // Original dimension
    int n_output;                    // Reduced dimension

    // Status
    bool converged;                  // True if converged without breakdown
};

/**
 * @brief Lanczos algorithm for CLN I model reduction
 *
 * Transforms two Hermitian matrices (K, N) simultaneously:
 *   K -> R_diag (diagonal)
 *   N -> L_tridiag (tridiagonal)
 *
 * For DC PEEC: lanczos(K=R, N=L) produces CLN I form
 *
 * @param K Input matrix (resistance R), row-major, n x n
 * @param N Input matrix (inductance L), row-major, n x n
 * @param n Matrix dimension
 * @param n_iter Number of iterations (reduced dimension). -1 for full transform.
 * @param tol Convergence tolerance
 * @return LanczosResult with reduced CLN I model
 */
LanczosResult lanczos(
    const double* K,
    const double* N,
    int n,
    int n_iter = -1,
    double tol = 1e-30
);

/**
 * @brief Build tridiagonal matrix from diagonal values
 *
 * Constructs CLN I inductance matrix:
 *   L[i,i] = diag[i]
 *   L[i,i+1] = L[i+1,i] = -diag[i+1]
 *
 * @param diag Diagonal values, size n
 * @param n Dimension
 * @return Tridiagonal matrix, row-major, n x n
 */
std::vector<double> build_tridiagonal(const double* diag, int n);

/**
 * @brief Compute CLN I impedance at a single frequency
 *
 * Z(s) = 1 / I[0] where (R_diag + s*L_tridiag) * I = V, V[0]=1
 *
 * @param R_diag Diagonal resistance matrix, row-major, n x n
 * @param L_tridiag Tridiagonal inductance matrix, row-major, n x n
 * @param n Dimension
 * @param freq Frequency in Hz
 * @return Complex impedance
 */
std::complex<double> compute_cln_impedance(
    const double* R_diag,
    const double* L_tridiag,
    int n,
    double freq
);

/**
 * @brief Compute CLN I impedance over frequency sweep
 *
 * @param R_diag Diagonal resistance matrix, row-major, n x n
 * @param L_tridiag Tridiagonal inductance matrix, row-major, n x n
 * @param n Dimension
 * @param freqs Frequency array in Hz
 * @param n_freqs Number of frequencies
 * @param Z_out Output impedance array (pre-allocated, size n_freqs)
 */
void compute_cln_impedance_sweep(
    const double* R_diag,
    const double* L_tridiag,
    int n,
    const double* freqs,
    int n_freqs,
    std::complex<double>* Z_out
);

// ============================================================
// Loop-Star Coupling Transformation
// ============================================================

/**
 * @brief Transform Loop-Star coupling matrix with Lanczos Q matrix
 *
 * Computes: M_LS_reduced = Q^T * M_LS
 *
 * For Loop-Star PEEC, the coupling matrix M_LS (n_loop x n_star) relates
 * Loop currents to Star potentials. After Lanczos reduction:
 *   - Original: M_LS (n_loop x n_star)
 *   - Reduced:  M_LS' = Q^T * M_LS (n_reduced x n_star)
 *
 * @param Q Transformation matrix from Lanczos (n_loop x n_reduced), row-major
 * @param M_LS Loop-Star coupling matrix (n_loop x n_star), row-major
 * @param n_loop Number of loop elements (original dimension)
 * @param n_reduced Reduced dimension
 * @param n_star Number of star elements
 * @return M_LS_reduced (n_reduced x n_star), row-major
 */
std::vector<double> transform_coupling(
    const double* Q,
    const double* M_LS,
    int n_loop,
    int n_reduced,
    int n_star
);

/**
 * @brief Transform port vector with Lanczos Q matrix
 *
 * Computes: v_reduced = Q^T * v
 *
 * For port excitation vector v (n_loop), transforms to reduced space.
 *
 * @param Q Transformation matrix from Lanczos (n_loop x n_reduced), row-major
 * @param v Port vector (n_loop)
 * @param n_loop Number of loop elements (original dimension)
 * @param n_reduced Reduced dimension
 * @return v_reduced (n_reduced)
 */
std::vector<double> transform_port_vector(
    const double* Q,
    const double* v,
    int n_loop,
    int n_reduced
);

/**
 * @brief Full Loop-Star CLN system structure
 *
 * Represents the complete reduced PEEC system:
 *
 *   [R_diag + s*L_tridiag   s*M_LS_reduced] [I_L_reduced]   [V_L_reduced]
 *   [s*M_LS_reduced^T       P/s           ] [I_S        ] = [V_S        ]
 *
 * Where:
 *   - R_diag, L_tridiag: CLN I form (n_reduced x n_reduced)
 *   - M_LS_reduced: Transformed coupling (n_reduced x n_star)
 *   - P: Potential coefficient matrix (n_star x n_star, unchanged)
 */
struct LoopStarCLN {
    // CLN I form (Loop-Loop reduced)
    std::vector<double> R_diag;      // (n_reduced x n_reduced)
    std::vector<double> L_tridiag;   // (n_reduced x n_reduced)

    // Loop-Star coupling (reduced)
    std::vector<double> M_LS_reduced;  // (n_reduced x n_star)

    // Star-Star (unchanged, pointer to original)
    const double* P;                 // (n_star x n_star), NOT owned

    // Dimensions
    int n_reduced;
    int n_star;

    // Transformation matrix (for port vector transformation)
    std::vector<double> Q;           // (n_loop x n_reduced)
    int n_loop;
};

/**
 * @brief Compute full Loop-Star CLN impedance at a frequency
 *
 * Solves the complete PEEC system including capacitive effects:
 *
 *   Z_LL = R_diag + s*L_tridiag
 *   Z_LS = s*M_LS_reduced
 *   Z_SL = Z_LS^T (reciprocity)
 *   Z_SS = P/s
 *
 * @param cln LoopStarCLN structure
 * @param freq Frequency in Hz
 * @param port_loop Port vector in loop space (n_loop), will be transformed
 * @return Complex impedance
 */
std::complex<double> compute_loop_star_impedance(
    const LoopStarCLN& cln,
    double freq,
    const double* port_loop
);

// ============================================================
// ACA+ Low-Rank Approximation for P Matrix (HACApK)
// ============================================================

/**
 * @brief Result of ACA+ low-rank approximation
 *
 * P ≈ U * V^T where U: (n x k), V: (n x k)
 *
 * This enables Star-space reduction:
 * - Original: n_star charge panels -> n_star x n_star interactions
 * - Reduced: k charge modes -> k x k interactions
 */
struct ACAResult {
    std::vector<double> U;      // U matrix (n_star x k), row-major
    std::vector<double> V;      // V matrix (n_star x k), row-major
    int n;                      // Original dimension (n_star)
    int k;                      // Rank after compression
    double compression_ratio;   // k / n_star
    bool converged;             // True if ACA+ converged

    // Reconstruct full matrix P ≈ U * V^T (for verification)
    std::vector<double> reconstruct() const;
};

/**
 * @brief ACA+ low-rank approximation using HACApK
 *
 * Computes P ≈ U * V^T using ACA+ algorithm.
 *
 * @param P Input potential coefficient matrix (n x n), row-major
 * @param n Matrix dimension
 * @param eps ACA+ tolerance (lower = more accurate, higher = faster)
 * @param kmax Maximum rank (-1 for automatic: kmax = n)
 * @return ACAResult with low-rank factors U, V
 */
ACAResult aca_compress(
    const double* P,
    int n,
    double eps = 1e-4,
    int kmax = -1
);

/**
 * @brief Full Loop-Star CLN with ACA+ compressed Star
 *
 * Complete reduced PEEC system with Star-space compression:
 *
 *   [R_diag + s*L_tridiag   s*M_LS_mode] [I_L_reduced]   [V_L_reduced]
 *   [s*M_LS_mode^T          P_mode/s   ] [I_mode     ] = [V_mode     ]
 *
 * Where:
 *   - R_diag, L_tridiag: CLN I form (n_reduced x n_reduced)
 *   - M_LS_mode = M_LS_reduced * V: (n_reduced x k_star)
 *   - P_mode = U^T * P * V ≈ U^T * (U * V^T) * V = U^T * U * V^T * V
 *
 * This achieves double reduction:
 *   - Loop: n_loop -> n_reduced (Lanczos)
 *   - Star: n_star -> k_star (ACA+)
 */
struct LoopStarCLNCompressed {
    // CLN I form (Loop-Loop reduced)
    std::vector<double> R_diag;      // (n_reduced x n_reduced)
    std::vector<double> L_tridiag;   // (n_reduced x n_reduced)

    // Loop-Star coupling in mode space
    std::vector<double> M_LS_mode;   // (n_reduced x k_star)

    // Star-Star in mode space (P_mode = U^T * P * V)
    std::vector<double> P_mode;      // (k_star x k_star)

    // ACA+ factors for Star transformation
    std::vector<double> U;           // (n_star x k_star)
    std::vector<double> V;           // (n_star x k_star)

    // Dimensions
    int n_reduced;                   // Reduced loop dimension
    int n_star;                      // Original star dimension
    int k_star;                      // Compressed star dimension

    // Transformation matrices (for port vector)
    std::vector<double> Q;           // (n_loop x n_reduced)
    int n_loop;
};

/**
 * @brief Compute compressed Loop-Star CLN impedance
 *
 * @param cln LoopStarCLNCompressed structure
 * @param freq Frequency in Hz
 * @param port_loop Port vector in original loop space (n_loop)
 * @return Complex impedance
 */
std::complex<double> compute_compressed_impedance(
    const LoopStarCLNCompressed& cln,
    double freq,
    const double* port_loop
);

}  // namespace cln
}  // namespace radia
