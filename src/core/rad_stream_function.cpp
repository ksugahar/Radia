/*-------------------------------------------------------------------------
 * rad_stream_function.cpp
 *
 * (ACA+)+TSVD accelerated least-norm solver -- generic, kernel-agnostic.
 * See rad_stream_function.h.
 *
 * This file implements ONLY:
 *   1. ACAPlus  -- a thin adapter that hands the caller's matrix-entry
 *      callback A(i,j) to HACApK's ACA+ (cHACApK_acaplus, the single source of
 *      truth for ACA+ in Radia).  No field kernel is embedded.
 *   2. ACATSVD  -- the standard SVD-of-a-low-rank-product recompression of the
 *      small ACA factors (QR of each tall-skinny factor + ONE small SVD; peer
 *      review JIAM-2026-36), which HACApK does not provide.
 *      LAPACKE_dgeqrf/dorgqr + dgesdd + cblas_dgemm.
 *   3. PseudoInverseSolve -- the least-norm back-substitution.
 *
 * Internally matrices are COLUMN-MAJOR (LAPACK / HACApK native: cHACApK_acaplus
 * writes zaa(ndl,kmax)=C, zab(ndt,kmax)=D column-major).  Only the final U / V
 * are converted to row-major for the NumPy interface.
 *-------------------------------------------------------------------------*/
#include "rad_stream_function.h"

#include <mkl.h>

// HACApK ACA+ (cHACApK_acaplus) + the per-call entry-function override
// (HACApK_set_entry_func / HACApK_clear_entry_func, defined in rad_hacapk.cpp).
// cHACApK_cpp.h wraps its includes in extern "C" for C++ use.
#include "../ext/HACApK/cHACApK_cpp.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

//-------------------------------------------------------------------------
// File-static handle to the caller's matrix-entry callback + a C-linkage
// trampoline that HACApK's ACA+ invokes (cHACApK_entry_ij -> g_entry_override
// -> here).  Set around the single cHACApK_acaplus call in ACAPlus
// (GIL-serialized; no concurrent MMM build), so a plain static pointer suffices.
//
// HACApK passes 1-based indices, row (field) index first, column (source)
// index second (cHACApK_calc_vec convention, identity lod); we convert to the
// 0-based (i,j) of the public EntryFn.
//-------------------------------------------------------------------------
namespace {
const radia::stream_function::EntryFn* g_stream_entry = nullptr;
}  // namespace

extern "C" double rad_stream_entry_func(int i, int j, int /*i_bemv*/) {
    return (*g_stream_entry)(i - 1, j - 1);
}

namespace radia {
namespace stream_function {

namespace {

// SVD (economy, jobz='S') of a column-major m x n matrix a (overwritten).
// On return: u is m x mn, s is mn, vt is mn x n, where mn = min(m,n).
// Caller sizes u (m*mn), s (mn), vt (mn*n).  Throws on LAPACK failure.
void SvdEcon(int m, int n, double* a, double* u, double* s, double* vt) {
    const int mn = std::min(m, n);
    const lapack_int info = LAPACKE_dgesdd(
        LAPACK_COL_MAJOR, 'S', m, n, a, m, s, u, m, vt, mn);
    if (info != 0)
        throw std::runtime_error("rad_stream_function: LAPACKE_dgesdd failed");
}

// Economy QR of a column-major m x n matrix `a` (m >= n), via LAPACKE_dgeqrf +
// LAPACKE_dorgqr.  On return: Q is m x n (orthonormal columns, col-major), R is
// n x n upper-triangular (col-major).  `a` is read-only (copied into Q first).
// This is the tall-skinny orthogonaliser of the standard "SVD of a low-rank
// product" recompression (cheaper than an SVD of the same factor).
void QrEcon(int m, int n, const double* a,
            std::vector<double>& Q, std::vector<double>& R) {
    Q.assign(static_cast<size_t>(m) * n, 0.0);
    std::copy(a, a + static_cast<size_t>(m) * n, Q.begin());   // dgeqrf is in-place
    std::vector<double> tau(std::min(m, n));
    lapack_int info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, m, n, Q.data(), m, tau.data());
    if (info != 0)
        throw std::runtime_error("rad_stream_function: LAPACKE_dgeqrf failed");
    // R = upper-triangular n x n block of the factored matrix.
    R.assign(static_cast<size_t>(n) * n, 0.0);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i <= j; ++i)
            R[i + static_cast<size_t>(j) * n] = Q[i + static_cast<size_t>(j) * m];
    // Overwrite Q with the explicit economy orthonormal factor (m x n).
    info = LAPACKE_dorgqr(LAPACK_COL_MAJOR, m, n, n, Q.data(), m, tau.data());
    if (info != 0)
        throw std::runtime_error("rad_stream_function: LAPACKE_dorgqr failed");
}

}  // namespace

//-------------------------------------------------------------------------
// ACA+ : delegated to HACApK (cHACApK_acaplus).  We only install the caller's
// entry callback as HACApK's matrix-entry source and translate parameters.
//   C = zaa (M x kmax col-major), D = zab (N x kmax col-major); only the
//   first k_aca columns are valid.  A(i,j) ~= sum_k C(i,k) D(j,k).
//
// Parameter mapping to reproduce the reference behaviour exactly:
//   param[61] = 1  -> ACA_EPS = pACA_EPS (absolute), apxnorm = first blknorm
//   param[64] = 1  -> minimum rank guard (matches the reference k>=1)
//   eps       = 1e-12 (relative convergence tolerance, = zeps in the reference)
//   pACA_EPS  = aca_eps (user absolute pivot/row/col threshold)
//   znrmmat   = 1.0 (unused when param[61]==1)
//-------------------------------------------------------------------------
int ACAPlus(int M, int N, const EntryFn& entry,
            int kmax, double aca_eps,
            std::vector<double>& C, std::vector<double>& D) {
    if (M <= 0 || N <= 0)
        throw std::invalid_argument("rad_stream_function: empty problem");
    if (!entry)
        throw std::invalid_argument("rad_stream_function: null entry callback");
    if (kmax < 1) kmax = std::min(M, N);
    kmax = std::min(kmax, std::min(M, N));

    // HACApK ACA+ output layout: zaa(ndl,kmax) col-major, zab(ndt,kmax) col-major.
    C.assign(static_cast<size_t>(M) * kmax, 0.0);
    D.assign(static_cast<size_t>(N) * kmax, 0.0);

    // Identity permutation, 1-based-indexable up to max(M,N) (cHACApK_calc_vec
    // reads lod[1..max(M,N)]).  lod[k]=k makes the raw index reach the callback.
    std::vector<int> lod(static_cast<size_t>(std::max(M, N)) + 1);
    for (size_t k = 0; k < lod.size(); ++k) lod[k] = static_cast<int>(k);

    std::vector<double> param(101, 0.0);
    param[61] = 1.0;  // ACA norm mode: absolute ACA_EPS, apxnorm = first blknorm
    param[64] = 1.0;  // minimum rank guard

    g_stream_entry = &entry;
    HACApK_set_entry_func(&rad_stream_entry_func);
    const int kt = cHACApK_acaplus(
        C.data(), D.data(), param.data(), M, N, /*nstrtl=*/1, /*nstrtt=*/1,
        lod.data(), /*i_bemv=*/0, kmax, /*eps=*/1.0e-12, /*znrmmat=*/1.0,
        /*pACA_EPS=*/aca_eps);
    HACApK_clear_entry_func();
    g_stream_entry = nullptr;
    return kt;
}

//-------------------------------------------------------------------------
TSVDResult ACATSVD(int M, int N, const EntryFn& entry,
                   int modes, int kmax, double aca_eps) {
    if (M <= 0 || N <= 0) throw std::invalid_argument("rad_stream_function: empty problem");

    std::vector<double> C, D;
    const int kt = ACAPlus(M, N, entry, kmax, aca_eps, C, D);
    if (kt < 1) throw std::runtime_error("rad_stream_function: ACA+ produced rank 0");

    const int m = std::min(modes, kt);  // cannot return more than kt triplets

    // Standard SVD-of-a-low-rank-product recompression (peer review JIAM-2026-36).
    // A ~= C D^T with C (M x kt), D (N x kt), both tall-skinny (kt <= min(M,N)).
    // Orthogonalise the factors by QR (cheaper than an SVD of the same factor):
    //   C = Qc Rc,  D = Qd Rd   ->   A = Qc (Rc Rd^T) Qd^T.
    // ONE small kt x kt SVD  Rc Rd^T = Um Sm VTm  then gives the exact SVD of A:
    //   A = (Qc Um) Sm (Qd VTm^T)^T,   U = Qc Um,  V = Qd VTm^T,  S = Sm.
    // (The legacy manuscript Method 2/3 -- two/three SVDs -- were removed here;
    //  the lesson lives in memory/aca_tsvd_qr_recompression.md.)
    std::vector<double> Qc, Rc, Qd, Rd;
    QrEcon(M, kt, C.data(), Qc, Rc);   // Qc (M x kt), Rc (kt x kt)
    QrEcon(N, kt, D.data(), Qd, Rd);   // Qd (N x kt), Rd (kt x kt)

    // Middle = Rc * Rd^T   (kt x kt)
    std::vector<double> Middle(static_cast<size_t>(kt) * kt);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, kt, kt, kt,
                1.0, Rc.data(), kt, Rd.data(), kt, 0.0, Middle.data(), kt);
    // ONE small SVD of the kt x kt Middle -> Um (kt x kt), Sm (kt), VTm (kt x kt)
    std::vector<double> Um(static_cast<size_t>(kt) * kt), Sm(kt), VTm(static_cast<size_t>(kt) * kt);
    SvdEcon(kt, kt, Middle.data(), Um.data(), Sm.data(), VTm.data());

    std::vector<double> Ufinal, Sfinal(m), Vfinal;  // col-major M x m, m, N x m
    for (int i = 0; i < m; ++i) Sfinal[i] = Sm[i];
    // U = Qc * Um(:, 0:m)   (M x m)
    Ufinal.assign(static_cast<size_t>(M) * m, 0.0);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, m, kt,
                1.0, Qc.data(), M, Um.data(), kt, 0.0, Ufinal.data(), M);
    // V = Qd * VTm^T(:, 0:m)   (N x m)
    Vfinal.assign(static_cast<size_t>(N) * m, 0.0);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, N, m, kt,
                1.0, Qd.data(), N, VTm.data(), kt, 0.0, Vfinal.data(), N);

    // Pack into row-major TSVDResult (NumPy C-contiguous).
    TSVDResult r;
    r.M = M; r.N = N; r.modes = m; r.k_aca = kt;
    r.S = Sfinal;
    r.U.assign(static_cast<size_t>(M) * m, 0.0);
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < m; ++j)
            r.U[static_cast<size_t>(i) * m + j] = Ufinal[i + static_cast<size_t>(j) * M];
    r.V.assign(static_cast<size_t>(N) * m, 0.0);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < m; ++j)
            r.V[static_cast<size_t>(i) * m + j] = Vfinal[i + static_cast<size_t>(j) * N];
    return r;
}

//-------------------------------------------------------------------------
std::vector<double> PseudoInverseSolve(const TSVDResult& t,
                                       const std::vector<double>& B, int k_mode) {
    if (static_cast<int>(B.size()) != t.M)
        throw std::invalid_argument("rad_stream_function: B length != M");
    const int km = std::max(1, std::min(k_mode, t.modes));
    // c[a] = (U^T B)[a] / S[a]   for a < km
    std::vector<double> c(km, 0.0);
    for (int a = 0; a < km; ++a) {
        double dot = 0.0;
        for (int i = 0; i < t.M; ++i) dot += t.U[static_cast<size_t>(i) * t.modes + a] * B[i];
        c[a] = (t.S[a] > 0.0) ? dot / t.S[a] : 0.0;
    }
    // phi = V(:, 0:km) * c
    std::vector<double> phi(t.N, 0.0);
    for (int i = 0; i < t.N; ++i) {
        double s = 0.0;
        for (int a = 0; a < km; ++a) s += t.V[static_cast<size_t>(i) * t.modes + a] * c[a];
        phi[i] = s;
    }
    return phi;
}

}  // namespace stream_function
}  // namespace radia
