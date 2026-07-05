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
 *   - magnets    : fixed-magnetization field kernels (rad_field_unified,
 *     soft iron    rad_interaction)
 * No field kernel is embedded here -- this module does ONLY (ACA+)+TSVD.
 *
 * ACA+ itself is delegated to the in-repo HACApK C library (cHACApK_acaplus,
 * src/ext/HACApK) -- the single source of truth for ACA+ in Radia.  The
 * recompression to a truncated SVD is the standard "SVD of a low-rank product"
 * (QR each tall-skinny ACA factor C, D, then ONE small SVD of the kt x kt core;
 * peer review JIAM-2026-36), the only numerical algorithm implemented in this
 * file because HACApK does not provide it.  It reproduces the dense TSVD of the
 * ACA approximation to machine precision.  (The legacy manuscript Method 2/3 --
 * two/three SVDs -- were removed; see memory/aca_tsvd_qr_recompression.md.)
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
// existing field computation (Biot-Savart, fixed-magnetization kernels, ...).
using EntryFn = std::function<double(int i, int j)>;

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

// (ACA+)+TSVD: factor A with ACA+, then recompress the small factors to a
// truncated SVD (standard QR-of-each-factor + one small SVD).
// Returns A ~= U diag(S) V^T truncated to `modes` (<= k_aca).
TSVDResult ACATSVD(int M, int N, const EntryFn& entry,
                   int modes, int kmax, double aca_eps);

// Least-norm pseudo-inverse solve:  phi = V diag(1/S) U^T B, using the first
// k_mode (<= result.modes) singular triplets (TSVD regularization).
// B has length M; returns phi of length N.
std::vector<double> PseudoInverseSolve(const TSVDResult& result,
                                       const std::vector<double>& B,
                                       int k_mode);

}  // namespace stream_function
}  // namespace radia

#endif  // RAD_STREAM_FUNCTION_H
