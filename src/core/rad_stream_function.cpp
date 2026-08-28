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
#include "rad_hacapk.h"

#include <mkl.h>

// HACApK ACA+ (cHACApK_acaplus) + the per-call entry-function override
// (HACApK_set_entry_func / HACApK_clear_entry_func, defined in rad_hacapk.cpp).
// cHACApK_cpp.h wraps its includes in extern "C" for C++ use.
#include "../ext/HACApK/cHACApK_cpp.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

//-------------------------------------------------------------------------
// File-static handle to the caller's matrix-entry callback + a C-linkage
// trampoline that HACApK's ACA+ invokes (cHACApK_entry_ij -> g_entry_override
// -> here).  Set around the single cHACApK_acaplus call in ACAPlus
// (GIL-serialized; no concurrent stream-function build), so a plain static pointer suffices.
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

TSVDResult DenseTSVD(int M, int N, const std::vector<double>& row_major,
                     int modes) {
    const int available = std::min(M, N);
    const int retained = std::min(modes, available);
    std::vector<double> matrix(static_cast<size_t>(M) * N);
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            matrix[i + static_cast<size_t>(j) * M] =
                row_major[static_cast<size_t>(i) * N + j];
    std::vector<double> u(static_cast<size_t>(M) * available);
    std::vector<double> s(available);
    std::vector<double> vt(static_cast<size_t>(available) * N);
    SvdEcon(M, N, matrix.data(), u.data(), s.data(), vt.data());
    TSVDResult result;
    result.M = M;
    result.N = N;
    result.modes = retained;
    result.k_aca = available;
    result.S.assign(s.begin(), s.begin() + retained);
    result.U.assign(static_cast<size_t>(M) * retained, 0.0);
    result.V.assign(static_cast<size_t>(N) * retained, 0.0);
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < retained; ++j)
            result.U[static_cast<size_t>(i) * retained + j] =
                u[i + static_cast<size_t>(j) * M];
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < retained; ++j)
            result.V[static_cast<size_t>(i) * retained + j] =
                vt[j + static_cast<size_t>(i) * available];
    return result;
}

std::vector<double> MatVec(int M, int N,
                           const std::vector<double>& matrix,
                           const std::vector<double>& value) {
    std::vector<double> result(M, 0.0);
    for (int i = 0; i < M; ++i) {
        double sum = 0.0;
        for (int j = 0; j < N; ++j)
            sum += matrix[static_cast<size_t>(i) * N + j] * value[j];
        result[i] = sum;
    }
    return result;
}

struct ResidualMetrics {
    double peak_to_peak = 0.0;
    double rms = 0.0;
    double max_abs = 0.0;
};

ResidualMetrics Metrics(const std::vector<double>& residual) {
    const auto range = std::minmax_element(residual.begin(), residual.end());
    double square_sum = 0.0;
    double max_abs = 0.0;
    for (double value : residual) {
        square_sum += value * value;
        max_abs = std::max(max_abs, std::abs(value));
    }
    return {*range.second - *range.first,
            std::sqrt(square_sum / static_cast<double>(residual.size())),
            max_abs};
}

bool TargetMet(const ResidualMetrics& metrics,
               const AbeBoundedOptions& options) {
    return (options.residual_peak_to_peak < 0.0 ||
            metrics.peak_to_peak <= options.residual_peak_to_peak) &&
           (options.residual_rms < 0.0 ||
            metrics.rms <= options.residual_rms);
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

    // ACA+ uses the same process-wide cHACApK_entry_ij override as the
    // matrix builders. Serialize the complete callback lifetime so a released-
    // GIL HDiv/PEEC build cannot replace the active entry source.
    std::lock_guard<std::mutex> operation_lock(
        RadHACApKCallback::OperationMutex());

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
    RadHACApKCallback::ClearCallbackException();
    HACApK_set_entry_func(&rad_stream_entry_func);
    const int kt = cHACApK_acaplus(
        C.data(), D.data(), param.data(), M, N, /*nstrtl=*/1, /*nstrtt=*/1,
        lod.data(), /*i_bemv=*/0, kmax, /*eps=*/1.0e-12, /*znrmmat=*/1.0,
        /*pACA_EPS=*/aca_eps);
    HACApK_clear_entry_func();
    g_stream_entry = nullptr;
    RadHACApKCallback::RethrowCallbackException();
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

//-------------------------------------------------------------------------
AbeBoundedResult SolveAbeBounded(
    int M, int N, const std::vector<double>& response,
    const std::vector<double>& target, const std::vector<double>& lower,
    const std::vector<double>& upper, const std::vector<double>& initial,
    const AbeBoundedOptions& options) {
    if (M <= 0 || N <= 0 ||
        response.size() != static_cast<size_t>(M) * N ||
        target.size() != static_cast<size_t>(M) ||
        lower.size() != static_cast<size_t>(N) ||
        upper.size() != static_cast<size_t>(N) ||
        initial.size() != static_cast<size_t>(N))
        throw std::invalid_argument(
            "rad_stream_function: invalid bounded Abe dimensions");
    const auto finite = [](double value) { return std::isfinite(value); };
    if (!std::all_of(response.begin(), response.end(), finite) ||
        !std::all_of(target.begin(), target.end(), finite) ||
        !std::all_of(initial.begin(), initial.end(), finite))
        throw std::invalid_argument(
            "rad_stream_function: bounded Abe arrays must be finite");
    for (int j = 0; j < N; ++j) {
        if (std::isnan(lower[j]) || std::isnan(upper[j]) ||
            lower[j] > upper[j] || initial[j] < lower[j] ||
            initial[j] > upper[j])
            throw std::invalid_argument(
                "rad_stream_function: invalid bounded Abe capacity");
    }
    if (options.max_iterations < 1 || !(options.relaxation > 0.0) ||
        options.relaxation > 1.0 || options.stagnation_tolerance < 0.0 ||
        options.relative_singular_threshold < 0.0 ||
        !(options.aca_eps > 0.0) ||
        (options.residual_peak_to_peak < 0.0 && options.residual_rms < 0.0))
        throw std::invalid_argument(
            "rad_stream_function: invalid bounded Abe options");

    const int max_rank = options.kmax <= 0
        ? std::min(M, N) : std::min(options.kmax, std::min(M, N));
    const int mode_count = options.modes <= 0
        ? max_rank : std::min(options.modes, max_rank);
    TSVDResult factor = options.dense_tsvd
        ? DenseTSVD(M, N, response, mode_count)
        : ACATSVD(M, N,
            [&response, N](int row, int col) {
                return response[static_cast<size_t>(row) * N + col];
            }, mode_count, max_rank, options.aca_eps);
    if (factor.modes < 1 || factor.S.empty() || !(factor.S.front() > 0.0))
        throw std::runtime_error(
            "rad_stream_function: bounded Abe factor has no usable mode");

    AbeBoundedResult result;
    result.factor = factor;
    std::vector<double> current = initial;
    const double singular_cutoff = options.relative_singular_threshold *
                                   factor.S.front();
    for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
        auto current_field = MatVec(M, N, response, current);
        std::vector<double> field_error(M);
        for (int i = 0; i < M; ++i)
            field_error[i] = target[i] - current_field[i];

        std::vector<double> correction(N, 0.0);
        std::vector<double> predicted_field = current_field;
        for (int mode = 0; mode < factor.modes; ++mode) {
            if (!(factor.S[mode] > singular_cutoff))
                continue;
            double strength = 0.0;
            for (int i = 0; i < M; ++i)
                strength += factor.U[static_cast<size_t>(i) *
                                     factor.modes + mode] * field_error[i];
            const double coefficient = strength / factor.S[mode];
            std::vector<double> mode_step(N);
            for (int j = 0; j < N; ++j) {
                mode_step[j] = factor.V[static_cast<size_t>(j) *
                                           factor.modes + mode] * coefficient;
                correction[j] += mode_step[j];
            }
            const auto mode_field = MatVec(M, N, response, mode_step);
            for (int i = 0; i < M; ++i)
                predicted_field[i] += mode_field[i];
            std::vector<double> predicted_residual(M);
            for (int i = 0; i < M; ++i)
                predicted_residual[i] = target[i] - predicted_field[i];
            if (TargetMet(Metrics(predicted_residual), options))
                break;
        }

        std::vector<double> bounded(N);
        int clipped = 0;
        double change = 0.0;
        for (int j = 0; j < N; ++j) {
            const double trial = current[j] +
                                 options.relaxation * correction[j];
            bounded[j] = std::clamp(trial, lower[j], upper[j]);
            if (bounded[j] != trial) ++clipped;
            change = std::max(change, std::abs(bounded[j] - current[j]));
        }
        const auto reconstructed = MatVec(M, N, response, bounded);
        std::vector<double> residual(M);
        for (int i = 0; i < M; ++i)
            residual[i] = target[i] - reconstructed[i];
        const auto metrics = Metrics(residual);
        result.clipped_history.push_back(clipped);
        result.potential_change_history.push_back(change);
        result.residual_peak_to_peak_history.push_back(metrics.peak_to_peak);
        result.residual_rms_history.push_back(metrics.rms);
        result.residual_max_abs_history.push_back(metrics.max_abs);
        current = std::move(bounded);
        result.reconstructed = reconstructed;
        result.residual = std::move(residual);
        result.residual_peak_to_peak = metrics.peak_to_peak;
        result.residual_rms = metrics.rms;
        result.residual_max_abs = metrics.max_abs;
        if (TargetMet(metrics, options)) {
            result.converged = true;
            result.stop_reason = "bounded_residual_target_met";
            break;
        }
        if (change <= options.stagnation_tolerance) {
            result.stop_reason = "bounded_stagnation";
            break;
        }
    }
    result.potential = std::move(current);
    result.iterations = static_cast<int>(result.clipped_history.size());
    if (result.stop_reason.empty())
        result.stop_reason = "bounded_max_iterations";
    return result;
}

}  // namespace stream_function
}  // namespace radia
