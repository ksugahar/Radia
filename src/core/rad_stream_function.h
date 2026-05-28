/*-------------------------------------------------------------------------
 * rad_stream_function.h
 *
 * (ACA+)+TSVD accelerated least-norm solver -- generic, kernel-agnostic.
 *
 * Solves the underdetermined least-norm system  A phi = B  (M field points x
 * N basis sources, M < N) via a TSVD-regularized pseudo-inverse.  Naive TSVD
 * on the dense A is O(N M^2); we instead factor A ~= C D^T with ACA+ (rank
 * k_aca << M,N) then TSVD only the small factors -> ~ (M/k)^2 faster.
 *
 * KERNEL-AGNOSTIC BY DESIGN.  The matrix entry A(i,j) is supplied by the
 * caller as a callback `double entry(int i, int j)` (0-based), so the SAME
 * machinery serves any Radia source family using Radia's already-implemented
 * field computation:
 *   - coils      : Biot-Savart H/A from filaments (rad_biot_savart_filaments)
 *   - magnets /  : MMM / MSC field from magnetization (rad_field_unified,
 *     soft iron    rad_interaction)
 * No field kernel is embedded here -- this module does ONLY (ACA+)+TSVD.
 *
 * ACA+ itself is delegated to the in-repo HACApK C library (cHACApK_acaplus,
 * src/ext/HACApK) -- the single source of truth for ACA+ in Radia.  The TSVD
 * recompression (manuscript Method 2/3, IEEJ SA-25-020) is the only numerical
 * algorithm implemented in this file, because HACApK does not provide it.  The
 * validated Fortran reference coil_solver.f90 (method_aca_tsvd_1/2) is a
 * faithful port of the same HACApK ACA+, so this matches it to machine
 * precision for the coil kernel.
 *
 * Storage convention: ROW-MAJOR for the returned U / V (matches Radia's
 * CblasRowMajor / NumPy C-contiguous arrays).
 *-------------------------------------------------------------------------*/
#ifndef RAD_STREAM_FUNCTION_H
#define RAD_STREAM_FUNCTION_H

#include <functional>
#include <vector>

namespace radia {
namespace stream_function {

// Matrix-entry callback:  A(i,j), with i in [0,M) (field/observation index)
// and j in [0,N) (basis/source index).  Supplied by the caller from Radia's
// existing field computation (Biot-Savart, MMM/MSC, ...).
using EntryFn = std::function<double(int i, int j)>;

enum class Method {
    Method2,   // full re-SVD of both ACA factors (manuscript Method 2)
    Method3    // improved: SVD(C), E = diag(Sc) Vc^T D^T, SVD(E) (Method 3, default)
};

// Recompressed truncated SVD of A:  A ~= U diag(S) V^T, truncated to `modes`.
struct TSVDResult {
    int M = 0, N = 0, modes = 0;   // dimensions actually returned
    int k_aca = 0;                 // ACA+ rank found before truncation
    std::vector<double> U;         // M x modes, row-major
    std::vector<double> S;         // modes
    std::vector<double> V;         // N x modes, row-major
};

// ACA+ low-rank factorization  A ~= C D^T  (delegated to HACApK cHACApK_acaplus).
//   C: M x kmax col-major, D: N x kmax col-major (only first k_aca cols valid).
// `entry(i,j)` provides A(i,j) on demand.  Returns the rank k_aca.
int ACAPlus(int M, int N, const EntryFn& entry,
            int kmax, double aca_eps,
            std::vector<double>& C, std::vector<double>& D);

// (ACA+)+TSVD: factor A with ACA+, then TSVD the small factors.
// Returns A ~= U diag(S) V^T truncated to `modes` (<= k_aca).
TSVDResult ACATSVD(int M, int N, const EntryFn& entry,
                   int modes, int kmax, double aca_eps,
                   Method method = Method::Method3);

// Least-norm pseudo-inverse solve:  phi = V diag(1/S) U^T B, using the first
// k_mode (<= result.modes) singular triplets (TSVD regularization).
// B has length M; returns phi of length N.
std::vector<double> PseudoInverseSolve(const TSVDResult& result,
                                       const std::vector<double>& B,
                                       int k_mode);

}  // namespace stream_function
}  // namespace radia

#endif  // RAD_STREAM_FUNCTION_H
