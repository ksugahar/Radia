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
 *   2. ACATSVD  -- the TSVD recompression of the small ACA factors (manuscript
 *      Method 2/3), which HACApK does not provide.  LAPACKE_dgesdd + cblas_dgemm.
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
                   int modes, int kmax, double aca_eps, Method method) {
    if (M <= 0 || N <= 0) throw std::invalid_argument("rad_stream_function: empty problem");

    std::vector<double> C, D;
    const int kt = ACAPlus(M, N, entry, kmax, aca_eps, C, D);
    if (kt < 1) throw std::runtime_error("rad_stream_function: ACA+ produced rank 0");

    const int m = std::min(modes, kt);  // cannot return more than kt triplets

    // SVD of C (M x kt) -> Uc (M x kt), Sc (kt), VTc (kt x kt).
    std::vector<double> Uc(static_cast<size_t>(M) * kt), Sc(kt), VTc(static_cast<size_t>(kt) * kt);
    SvdEcon(M, kt, C.data(), Uc.data(), Sc.data(), VTc.data());  // uses first kt cols of C

    std::vector<double> Ufinal, Sfinal(m), Vfinal;  // col-major M x m, m, N x m

    if (method == Method::Method3) {
        // E = diag(Sc) * VTc * D^T   (kt x N)
        std::vector<double> ScVTc(static_cast<size_t>(kt) * kt);
        for (int j = 0; j < kt; ++j)
            for (int i = 0; i < kt; ++i)
                ScVTc[i + static_cast<size_t>(j) * kt] = Sc[i] * VTc[i + static_cast<size_t>(j) * kt];
        std::vector<double> E(static_cast<size_t>(kt) * N);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, kt, N, kt,
                    1.0, ScVTc.data(), kt, D.data(), N, 0.0, E.data(), kt);
        // SVD of E (kt x N) -> UE (kt x kt), SE (kt), VTE (kt x N).
        std::vector<double> UE(static_cast<size_t>(kt) * kt), SE(kt), VTE(static_cast<size_t>(kt) * N);
        SvdEcon(kt, N, E.data(), UE.data(), SE.data(), VTE.data());
        for (int i = 0; i < m; ++i) Sfinal[i] = SE[i];
        // U = Uc * UE(:, 0:m)   (M x m)
        Ufinal.assign(static_cast<size_t>(M) * m, 0.0);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, m, kt,
                    1.0, Uc.data(), M, UE.data(), kt, 0.0, Ufinal.data(), M);
        // V(i,j) = VTE(j,i)   (N x m)
        Vfinal.assign(static_cast<size_t>(N) * m, 0.0);
        for (int j = 0; j < m; ++j)
            for (int i = 0; i < N; ++i)
                Vfinal[i + static_cast<size_t>(j) * N] = VTE[j + static_cast<size_t>(i) * kt];
    } else {
        // Method2: also SVD D, combine via Middle.
        std::vector<double> Ud(static_cast<size_t>(N) * kt), Sd(kt), VTd(static_cast<size_t>(kt) * kt);
        SvdEcon(N, kt, D.data(), Ud.data(), Sd.data(), VTd.data());
        // VcTVd = VTc * VTd^T  (kt x kt)
        std::vector<double> VcTVd(static_cast<size_t>(kt) * kt);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, kt, kt, kt,
                    1.0, VTc.data(), kt, VTd.data(), kt, 0.0, VcTVd.data(), kt);
        // Middle(i,j) = Sc(i) * VcTVd(i,j) * Sd(j)
        std::vector<double> Middle(static_cast<size_t>(kt) * kt);
        for (int j = 0; j < kt; ++j)
            for (int i = 0; i < kt; ++i)
                Middle[i + static_cast<size_t>(j) * kt] = Sc[i] * VcTVd[i + static_cast<size_t>(j) * kt] * Sd[j];
        std::vector<double> Um(static_cast<size_t>(kt) * kt), Sm(kt), VTm(static_cast<size_t>(kt) * kt);
        SvdEcon(kt, kt, Middle.data(), Um.data(), Sm.data(), VTm.data());
        for (int i = 0; i < m; ++i) Sfinal[i] = Sm[i];
        // U = Uc * Um(:, 0:m)
        Ufinal.assign(static_cast<size_t>(M) * m, 0.0);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, m, kt,
                    1.0, Uc.data(), M, Um.data(), kt, 0.0, Ufinal.data(), M);
        // V = Ud * VTm^T(:, 0:m)
        Vfinal.assign(static_cast<size_t>(N) * m, 0.0);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, N, m, kt,
                    1.0, Ud.data(), N, VTm.data(), kt, 0.0, Vfinal.data(), N);
    }

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
