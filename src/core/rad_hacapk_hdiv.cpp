/* rad_hacapk_hdiv.cpp -- HACApK H-matrix for the symmetric HDiv-type VIM demag operator.
 * See rad_hacapk_hdiv.h.  The build/matvec/stats lifecycle is inherited from RadHACApKBase;
 * this file supplies only the HDiv kernel hooks (coordinates = face centroids, entry = the
 * charge-cluster Coulomb sum N[i][j] = sum_a sum_b B[a][i] G[a][b] B[b][j]). */
#include "rad_hacapk_hdiv.h"
#include "rad_parallel.h"
extern "C" {
#include "../ext/HACApK/cHACApK_base.h"
}
#include <cmath>
#include <utility>
#include <algorithm>
#include <core/taskmanager.hpp>   // ngcore::RegionTaskManager (parallel H-matvec under TaskManager)
#include <core/utils.hpp>         // ngcore::AtomicAdd
#include <unordered_map>
#include <map>
#include <memory>
#include <stdexcept>
#include <atomic>
#include <functional>
#include <cstddef>
#include <cstdint>
#include <chrono>
#include <cstdlib>
#include <limits>

extern "C" {
void HACApK_matvec_stats_reset(void);
void HACApK_matvec_stats_get(double *values, int n_values,
                             int64_t *counts, int n_counts);
}

#ifdef HAVE_LAPACK
#include "mkl_pardiso.h"          // PARDISO sparse-direct factor of the HDiv mass for the MASS RIESZ precond
#include "mkl_cblas.h"            // row-major cluster-space dense contractions
#include "mkl_lapacke.h"          // small cluster Rayleigh--Ritz eigensolve
namespace {
// PARDISO is entered from the Python-facing solve while NGSolve may already
// own a TaskManager pool.  MKL's process-wide thread setting is insufficient
// here: NGSolve or another MKL caller can replace it between factor and solve.
// Use MKL's per-calling-thread setting and restore the previous local value.
class PardisoMKLThreadGuard {
    int saved_;
public:
    explicit PardisoMKLThreadGuard(int nthreads)
        : saved_(mkl_set_num_threads_local(std::max(1, nthreads))) {}
    ~PardisoMKLThreadGuard() { mkl_set_num_threads_local(saved_); }
    PardisoMKLThreadGuard(const PardisoMKLThreadGuard&) = delete;
    PardisoMKLThreadGuard& operator=(const PardisoMKLThreadGuard&) = delete;
};

// RAII PARDISO SPD (mtype=2 real symmetric positive definite) factor of the HDiv mass M_mass,
// used as the MASS RIESZ preconditioner (z = M_mass^{-1} r) of the HDiv-VIM material CG / MINRES.  The
// mass is supplied as the FULL symmetric COO (mI,mJ,mV); only the UPPER triangle (j>=i) is assembled
// into the 0-based CSR PARDISO mtype=2 expects.  Follows the established sparse-direct PARDISO pattern in
// this repo.  Replaces the prior Python splu(M_mass) glue so
// the whole linear demag solve (H-matvec + mass solve + Krylov) runs in C++.
struct MassRieszPardiso {
    void* pt[64];
    MKL_INT iparm[64];
    MKL_INT n = 0, mtype = 2, maxfct = 1, mnum = 1, msglvl = 0;
    std::vector<MKL_INT> ia, ja;     // upper-triangular CSR, 0-based (iparm[34]=1); columns ascending
    std::vector<double>  a;
    bool factored = false;
    MassRieszPardiso() { for (int i = 0; i < 64; ++i) { pt[i] = nullptr; iparm[i] = 0; } }
    MassRieszPardiso(const MassRieszPardiso&) = delete;
    MassRieszPardiso& operator=(const MassRieszPardiso&) = delete;
    ~MassRieszPardiso() {
        if (factored) {
            // PARDISO is called from the HDiv Krylov loop while an NGSolve
            // TaskManager may be active.  Suspend its workers while PARDISO
            // owns the configured thread count; this avoids nested pools while
            // retaining parallel sparse factor/solve performance.
            ngcore::SuspendTaskManager stm;
            PardisoMKLThreadGuard mkl_guard(radia::GetMaxThreads());
            MKL_INT phase = -1, nrhs = 1, idum = 0, error = 0; double ddum = 0.0;
            pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia.data(), ja.data(),
                    &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
        }
    }
    // Assemble the upper-triangular CSR from the symmetric mass COO and factor (analyze phase 11 +
    // numeric phase 22).  Returns false on a PARDISO error (caller raises -- No-Fallbacks: a non-SPD
    // mass would be a setup bug, not a soft condition to paper over).
    bool Factor(const std::vector<int>& mI, const std::vector<int>& mJ,
                const std::vector<double>& mV, int n_face) {
        ngcore::SuspendTaskManager stm;
        PardisoMKLThreadGuard mkl_guard(radia::GetMaxThreads());
        n = n_face;
        std::vector<std::map<int, double>> row((size_t)n_face);   // std::map keeps columns ascending
        for (size_t k = 0; k < mV.size(); ++k) {
            int i = mI[k], j = mJ[k];
            if (i < 0 || i >= n_face || j < 0 || j >= n_face) continue;
            if (j < i) continue;                                 // upper triangle only (M_mass symmetric)
            row[(size_t)i][j] += mV[k];                           // merge any duplicate COO entries
        }
        ia.assign((size_t)n_face + 1, 0);
        for (int i = 0; i < n_face; ++i)
            ia[(size_t)i + 1] = ia[(size_t)i] + (MKL_INT)row[(size_t)i].size();
        MKL_INT nnz = ia[(size_t)n_face];
        ja.assign((size_t)nnz, 0); a.assign((size_t)nnz, 0.0);
        MKL_INT k = 0;
        for (int i = 0; i < n_face; ++i)
            for (const auto& kv : row[(size_t)i]) { ja[(size_t)k] = (MKL_INT)kv.first; a[(size_t)k] = kv.second; ++k; }
        pardisoinit(pt, &mtype, iparm);
        // PARDISO's own worker pool must not overlap the surrounding NGSolve
        // TaskManager.  iparm[2] is the C zero-based slot for the documented
        // number-of-processors control (Fortran iparm(3)).  The workers are
        // suspended above, so PARDISO may use the configured Radia count.
        iparm[2] = std::max<MKL_INT>(1, (MKL_INT)radia::GetMaxThreads());
        iparm[34] = 1;                                           // 0-based (C) indexing
        MKL_INT phase = 11, nrhs = 1, idum = 0, error = 0; double ddum = 0.0;
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a.data(), ia.data(), ja.data(),
                &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
        if (error == 0) {
            phase = 22;
            pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a.data(), ia.data(), ja.data(),
                    &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
        }
        if (error != 0) return false;
        factored = true;
        return true;
    }
    void Solve(const double* rhs, double* x) {                   // M_mass x = rhs (phase 33, single rhs)
        ngcore::SuspendTaskManager stm;
        PardisoMKLThreadGuard mkl_guard(radia::GetMaxThreads());
        MKL_INT phase = 33, nrhs = 1, idum = 0, error = 0;
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a.data(), ia.data(), ja.data(),
                &idum, &nrhs, iparm, &msglvl, const_cast<double*>(rhs), x, &error);
        if (error != 0)
            throw std::runtime_error("MassRieszPardiso: PARDISO solve phase failed");
    }
    void SolveMany(const double* rhs, double* x, int rhs_count) {
        // PARDISO stores dense right-hand sides column-major [n][nrhs].
        // Radia's public row-major [nrhs][n] buffer is byte-identical: each
        // right-hand side is one contiguous PARDISO column.
        if (rhs_count < 1)
            throw std::runtime_error(
                "MassRieszPardiso: rhs_count must be positive");
        ngcore::SuspendTaskManager stm;
        PardisoMKLThreadGuard mkl_guard(radia::GetMaxThreads());
        MKL_INT phase = 33, nrhs = static_cast<MKL_INT>(rhs_count);
        MKL_INT idum = 0, error = 0;
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a.data(),
                ia.data(), ja.data(), &idum, &nrhs, iparm, &msglvl,
                const_cast<double*>(rhs), x, &error);
        if (error != 0)
            throw std::runtime_error(
                "MassRieszPardiso: batched PARDISO solve phase failed");
    }
};
} // namespace

// Definition of the .h-forward-declared persistent factor holder (SINGLE entry on RadHACApKChargeGram):
// the PARDISO factor plus the exact COO arrays it was built from.  A solve call whose (n_face, mI, mJ, mV)
// compares EQUAL element-wise reuses the factor; any difference rebuilds and replaces the entry.  Only this
// translation unit sees the complete type (the members use the TU-local MassRieszPardiso above).
struct RadMassRieszCache {
    std::vector<int> keyI, keyJ;
    std::vector<double> keyV;
    int keyN = -1;
    MassRieszPardiso factor;
};

// The single get-or-build implementation used by SolveLinearMaterial (see the .h declaration for the
// pinning / single-resident contracts).  Hits on constant-mass chains: the Hantila
// hysteresis loop (W(nu0) fixed by construction) and the C++ scalar Picard (geometry-only M_mass; the
// scalar inv_chi lives outside the preconditioner).  Per-iteration TANGENT masses (the Python nu-secant /
// Newton W_tan) compare-miss and refactor exactly as the pre-cache code did -- with the old entry released
// FIRST, so their peak preconditioner memory stays at one factor.  Identical input -> identical factor ->
// bit-identical preconditioner: a cache hit changes timing only.
std::shared_ptr<RadMassRieszCache> RadHACApKChargeGram::EnsureMassRieszFactor(
    const std::vector<int>& mI, const std::vector<int>& mJ, const std::vector<double>& mV,
    int n_face, const char* caller, double* factor_s_accum, bool geometry_cache)
{
    std::shared_ptr<RadMassRieszCache>& slot =
        geometry_cache ? m_geometryMassRieszCache : m_massRieszCache;
    std::shared_ptr<RadMassRieszCache> keep = slot;   // pin the current entry
    const bool hit = keep && keep->keyN == n_face && keep->keyV == mV &&
                     keep->keyI == mI && keep->keyJ == mJ;
    if (hit) return keep;
    keep.reset();
    slot.reset();                      // release the OLD factor before building the new one
    const auto t0 = std::chrono::steady_clock::now();
    auto built = std::make_shared<RadMassRieszCache>();
    if (!built->factor.Factor(mI, mJ, mV, n_face))
        throw std::runtime_error(std::string(caller) +
            ": PARDISO SPD factor of the HDiv mass (mass Riesz preconditioner) failed");
    built->keyN = n_face; built->keyI = mI; built->keyJ = mJ; built->keyV = mV;
    if (factor_s_accum)
        *factor_s_accum += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    slot = built;
    return built;
}
#endif // HAVE_LAPACK
// Per-thread memo for the high-order QuadDot lives inside QuadDot (a function-local static thread_local map):
// PhiAtHO(src, m_qp[tgt][k]) (the expensive analytic base + inner subtraction loop) depends ONLY on
// (kind,host of tgt, src) -- IDENTICAL across the co-located monomials that share a host's outer points -- so
// the H-matrix fill otherwise recomputes it n_mono(host) times per source.  See QuadDot for the rationale.

//=========================================================================
// RadHACApKChargeGram -- charge-charge Coulomb Gram G as a HACApK H-matrix
//=========================================================================

static const double RAD_INV_FOUR_PI = 0.07957747154594766788;   // 1/(4 pi)

// --- high-order helpers: integer power + small matrix inverses (row-major) ---
static inline double rad_ipow(double b, int e) { double r = 1.0; for (int i = 0; i < e; ++i) r *= b; return r; }

static void rad_inv3x3(const double A[9], double Ai[9])   // inverse of a row-major 3x3 (A[r*3+c])
{
    const double det = A[0]*(A[4]*A[8]-A[5]*A[7]) - A[1]*(A[3]*A[8]-A[5]*A[6]) + A[2]*(A[3]*A[7]-A[4]*A[6]);
    const double iv = 1.0 / det;
    Ai[0] =  (A[4]*A[8]-A[5]*A[7])*iv; Ai[1] = -(A[1]*A[8]-A[2]*A[7])*iv; Ai[2] =  (A[1]*A[5]-A[2]*A[4])*iv;
    Ai[3] = -(A[3]*A[8]-A[5]*A[6])*iv; Ai[4] =  (A[0]*A[8]-A[2]*A[6])*iv; Ai[5] = -(A[0]*A[5]-A[2]*A[3])*iv;
    Ai[6] =  (A[3]*A[7]-A[4]*A[6])*iv; Ai[7] = -(A[0]*A[7]-A[1]*A[6])*iv; Ai[8] =  (A[0]*A[4]-A[1]*A[3])*iv;
}

static void rad_inv2x2(const double A[4], double Ai[4])    // inverse of a row-major 2x2
{
    const double iv = 1.0 / (A[0]*A[3] - A[1]*A[2]);
    Ai[0] =  A[3]*iv; Ai[1] = -A[1]*iv; Ai[2] = -A[2]*iv; Ai[3] =  A[0]*iv;
}

// Built-in 64-node Gauss-Duffy collapsed-cube tet rule (4 Gauss-Legendre pts/dim).  ref pts are
// barycentric (lam1,lam2,lam3) flat in `pts`, weights summing to 1/6 in `w` -> phys weight = w*|J|,
// |J| = 6*vol.  This is the SAME rule as radia.vim._core._gauss_duffy_tet(4) (so the C++ analytic
// charge-Gram matches the independent analytic reference).  Shared by the tet analytic ctor (outer quad on the
// tet itself) and the polytope ctor (outer quad on each centroid-fan sub-tet).
static void rad_gl4_duffy_tet(std::vector<double>& pts, std::vector<double>& w)
{
    static const double GL4x[4] = {0.06943184420297371, 0.33000947820757187,
                                   0.66999052179242813, 0.93056815579702629};   // 4-pt Gauss-Legendre on [0,1]
    static const double GL4w[4] = {0.17392742256872693, 0.32607257743127307,
                                   0.32607257743127307, 0.17392742256872693};
    pts.clear(); w.clear();
    pts.reserve(64 * 3); w.reserve(64);
    for (int ia = 0; ia < 4; ++ia)
        for (int ib = 0; ib < 4; ++ib)
            for (int ic = 0; ic < 4; ++ic) {
                const double aa = GL4x[ia], bb = GL4x[ib], cc = GL4x[ic];
                pts.push_back(aa);
                pts.push_back(bb * (1.0 - aa));
                pts.push_back(cc * (1.0 - aa) * (1.0 - bb));
                w.push_back(GL4w[ia] * GL4w[ib] * GL4w[ic] * (1.0 - aa) * (1.0 - aa) * (1.0 - bb));
            }
}

// Built-in Dunavant degree-5 symmetric triangle rule (7 nodes; bary (l1,l2,l3) + weight, weights sum to 1).
static const double RAD_DUN5[7][4] = {
    {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.225},
    {0.0597158717, 0.4701420641, 0.4701420641, 0.1323941527},
    {0.4701420641, 0.0597158717, 0.4701420641, 0.1323941527},
    {0.4701420641, 0.4701420641, 0.0597158717, 0.1323941527},
    {0.7974269853, 0.1012865073, 0.1012865073, 0.1259391805},
    {0.1012865073, 0.7974269853, 0.1012865073, 0.1259391805},
    {0.1012865073, 0.1012865073, 0.7974269853, 0.1259391805},
};

static void ValidateImageVectors(const std::vector<int>& image_masks,
                                 const std::vector<double>& image_signs);

RadHACApKChargeGram::RadHACApKChargeGram(std::vector<double> centroids,
                                         std::vector<double> measures,
                                         std::vector<double> self_energy)
    : m_cent(std::move(centroids)), m_meas(std::move(measures)), m_self(std::move(self_energy))
{
    m_n = (int)m_meas.size();
}

RadHACApKChargeGram::RadHACApKChargeGram(std::vector<double> points,
                                         std::vector<double> weights,
                                         double kernel_epsilon)
    : m_cent(std::move(points)), m_meas(std::move(weights)),
      m_sampledLaplace(true), m_sampledKernelEpsilon(kernel_epsilon)
{
    m_n = static_cast<int>(m_meas.size());
    if (m_n < 1 || m_cent.size() != static_cast<size_t>(3 * m_n))
        throw std::invalid_argument(
            "sampled Laplace HACApK requires points[n,3] and weights[n]");
    if (!std::isfinite(m_sampledKernelEpsilon) || m_sampledKernelEpsilon <= 0.0)
        throw std::invalid_argument("sampled Laplace kernel_epsilon must be positive");
    for (double value : m_cent)
        if (!std::isfinite(value))
            throw std::invalid_argument("sampled Laplace points must be finite");
    for (double weight : m_meas)
        if (!std::isfinite(weight) || weight <= 0.0)
            throw std::invalid_argument("sampled Laplace weights must be positive and finite");
}

RadHACApKChargeGram::RadHACApKChargeGram(std::vector<double> points,
                                         std::vector<double> weights,
                                         double kernel_epsilon,
                                         double reference_length)
    : m_cent(std::move(points)), m_meas(std::move(weights)),
      m_sampledPlanarLog(true), m_sampledKernelEpsilon(kernel_epsilon),
      m_sampledReferenceLength(reference_length)
{
    m_n = static_cast<int>(m_meas.size());
    if (m_n < 1 || m_cent.size() != static_cast<size_t>(3 * m_n))
        throw std::invalid_argument(
            "sampled planar-log HACApK requires points[n,3] and weights[n]");
    if (!std::isfinite(m_sampledKernelEpsilon) || m_sampledKernelEpsilon <= 0.0)
        throw std::invalid_argument("sampled planar-log kernel_epsilon must be positive");
    if (!std::isfinite(m_sampledReferenceLength) || m_sampledReferenceLength <= 0.0)
        throw std::invalid_argument("sampled planar-log reference_length must be positive");
    for (double value : m_cent)
        if (!std::isfinite(value))
            throw std::invalid_argument("sampled planar-log points must be finite");
    for (double weight : m_meas)
        if (!std::isfinite(weight) || weight <= 0.0)
            throw std::invalid_argument(
                "sampled planar-log weights must be positive and finite");
}

RadHACApKChargeGram::RadHACApKChargeGram(std::vector<double> cell_verts,
                                         std::vector<double> face_verts,
                                         int n_el, double near_factor,
                                         std::vector<int> image_masks, std::vector<double> image_signs,
                                         int far_quad)
    : m_n_el(n_el), m_analytic(true), m_near_factor(near_factor), m_far_quad(far_quad),
      m_cellV(std::move(cell_verts)), m_faceV(std::move(face_verts)),
      m_image_masks(std::move(image_masks)), m_image_signs(std::move(image_signs))
{
    ValidateImageVectors(m_image_masks, m_image_signs);
    const int n_bf = (int)(m_faceV.size() / 9);
    m_n = n_el + n_bf;
    m_cent.assign((size_t)m_n * 3, 0.0);
    m_meas.assign((size_t)m_n, 0.0);    // measure (cell vol / face area) -- for the near/far split monopole
    m_size.assign((size_t)m_n, 0.0);    // characteristic size (vol^1/3 / area^1/2) -- for the near criterion
    m_qp.resize(m_n);
    m_qw.resize(m_n);
    if (m_far_quad > 0) { m_qpf.resize(m_n); m_qwf.resize(m_n); }   // low-order FAR double-quad rule
    // degree-2 symmetric rules (weights sum to 1; scaled by measure below) -- the same rules the Python
    // Gauss point cloud / the validated prototype use: 4-pt tet (a,b barycentric), 3-pt tri (2/3,1/6).
    const double ta = 0.5854101966249685, tb = 0.1381966011250105;
    const double TETF[4][4] = {{ta,tb,tb,tb},{tb,ta,tb,tb},{tb,tb,ta,tb},{tb,tb,tb,ta}};
    const double TRIF[3][3] = {{2.0/3,1.0/6,1.0/6},{1.0/6,2.0/3,1.0/6},{1.0/6,1.0/6,2.0/3}};

    // Outer-quad rule on a CELL: a built-in 4-pt Gauss-Duffy collapsed-cube tet rule (4^3 = 64 nodes; ref-tet
    // barycentric (lam1,lam2,lam3) flat in ref_tet_pts, weights summing to 1/6 in ref_tet_w).  The order-0
    // charge is CONSTANT so the inner is the EXACT analytic PhiTet and the cell self-integral INT_T PhiTet dx
    // is smooth -- 4 pts/dim integrates it to ~1e-4.  (The old hardcoded equal-weight _bary_tet(3) rule
    // under-integrated the volume self-energy by ~6.5% -- invisible to every uniform-M demag golden because
    // div M = 0 there.  This is the same rule as radia.vim._vim._tet_ref(4).)
    std::vector<double> ref_tet_pts, ref_tet_w;
    rad_gl4_duffy_tet(ref_tet_pts, ref_tet_w);   // 64-node Gauss-Duffy tet rule (shared w/ the polytope ctor)
    const int nqt = (int)ref_tet_w.size();   // 64
    // Outer-quad rule on a FACE: Dunavant degree-5 symmetric triangle rule (RAD_DUN5; 7 nodes, sum to 1).
    const double (*DUN)[4] = RAD_DUN5;

    for (int a = 0; a < n_el; ++a) {
        const double* V = &m_cellV[(size_t)a * 12];   // 4 x 3
        double cx = 0, cy = 0, cz = 0;
        for (int i = 0; i < 4; ++i) { cx += V[3*i]; cy += V[3*i+1]; cz += V[3*i+2]; }
        m_cent[3*a] = cx / 4; m_cent[3*a+1] = cy / 4; m_cent[3*a+2] = cz / 4;
        double e1[3], e2[3], e3[3];
        for (int k = 0; k < 3; ++k) { e1[k] = V[3+k]-V[k]; e2[k] = V[6+k]-V[k]; e3[k] = V[9+k]-V[k]; }
        double cr[3] = {e2[1]*e3[2]-e2[2]*e3[1], e2[2]*e3[0]-e2[0]*e3[2], e2[0]*e3[1]-e2[1]*e3[0]};
        double vol = std::fabs(e1[0]*cr[0] + e1[1]*cr[1] + e1[2]*cr[2]) / 6.0;
        m_meas[a] = vol; m_size[a] = std::cbrt(vol);
        m_qp[a].resize(nqt);
        m_qw[a].resize(nqt);
        const double absJ = 6.0 * vol;     // |J| = 6*vol; ref_tet_w sums to 1/6 -> phys weights sum to vol
        for (int q = 0; q < nqt; ++q) {
            const double l1 = ref_tet_pts[3*q], l2 = ref_tet_pts[3*q+1], l3 = ref_tet_pts[3*q+2];
            rad_hdiv::Vec3 P;
            for (int k = 0; k < 3; ++k)
                P[k] = V[k] + l1*(V[3+k]-V[k]) + l2*(V[6+k]-V[k]) + l3*(V[9+k]-V[k]);
            m_qp[a][q] = P;
            m_qw[a][q] = ref_tet_w[q] * absJ;
        }
        if (m_far_quad > 0) {              // 4-pt degree-2 FAR rule (barycentric TETF; weights = vol/4)
            m_qpf[a].resize(4); m_qwf[a].resize(4);
            for (int q = 0; q < 4; ++q) {
                rad_hdiv::Vec3 P = {0, 0, 0};
                for (int i = 0; i < 4; ++i) for (int k = 0; k < 3; ++k) P[k] += TETF[q][i] * V[3*i+k];
                m_qpf[a][q] = P; m_qwf[a][q] = 0.25 * vol;
            }
        }
    }
    for (int b = 0; b < n_bf; ++b) {
        int a = n_el + b;
        const double* V = &m_faceV[(size_t)b * 9];    // 3 x 3
        double cx = 0, cy = 0, cz = 0;
        for (int i = 0; i < 3; ++i) { cx += V[3*i]; cy += V[3*i+1]; cz += V[3*i+2]; }
        m_cent[3*a] = cx / 3; m_cent[3*a+1] = cy / 3; m_cent[3*a+2] = cz / 3;
        double e1[3], e2[3];
        for (int k = 0; k < 3; ++k) { e1[k] = V[3+k]-V[k]; e2[k] = V[6+k]-V[k]; }
        double cr[3] = {e1[1]*e2[2]-e1[2]*e2[1], e1[2]*e2[0]-e1[0]*e2[2], e1[0]*e2[1]-e1[1]*e2[0]};
        double area = 0.5 * std::sqrt(cr[0]*cr[0] + cr[1]*cr[1] + cr[2]*cr[2]);
        m_meas[a] = area; m_size[a] = std::sqrt(area);
        m_qp[a].resize(7);
        m_qw[a].resize(7);
        for (int q = 0; q < 7; ++q) {
            rad_hdiv::Vec3 P = {0, 0, 0};
            for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) P[k] += DUN[q][i] * V[3*i+k];
            m_qp[a][q] = P;
            m_qw[a][q] = DUN[q][3] * area;
        }
        if (m_far_quad > 0) {              // 3-pt degree-2 FAR rule (barycentric TRIF; weights = area/3)
            m_qpf[a].resize(3); m_qwf[a].resize(3);
            for (int q = 0; q < 3; ++q) {
                rad_hdiv::Vec3 P = {0, 0, 0};
                for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) P[k] += TRIF[q][i] * V[3*i+k];
                m_qpf[a][q] = P; m_qwf[a][q] = area / 3.0;
            }
        }
    }
}

// POLYTOPE constructor (hex/wedge cells + quad faces): the triangulation is supplied from Python (cell
// hull tris / face sub-tris as flat triangle soups + CSR offsets).  Builds the SAME analytic charge Gram
// as the tet/triangle ctor, generalized to any flat-faced convex cell: cell outer quad = centroid-fan
// sub-tets (apex = cell_cent) each filled by the 64-node Gauss-Duffy rule; face outer quad = Dunavant-5
// per sub-triangle.  The source potential (PhiAt) is the divergence-theorem polytope potential (cell) /
// sum-of-sub-triangle Wilton potential (face), evaluated from m_srcTris.  Matches the independent analytic
// polytope reference entry-by-entry (same tris, same quad rules).
RadHACApKChargeGram::RadHACApKChargeGram(
    std::vector<double> cell_tris, std::vector<int> cell_troff,
    std::vector<double> cell_cent, std::vector<double> cell_meas,
    std::vector<double> face_tris, std::vector<int> face_troff,
    std::vector<double> face_cent, std::vector<double> face_meas,
    int n_el, double near_factor,
    std::vector<int> image_masks, std::vector<double> image_signs,
    int far_quad)
    : m_n_el(n_el), m_analytic(true), m_near_factor(near_factor), m_far_quad(far_quad), m_polytope(true),
      m_image_masks(std::move(image_masks)), m_image_signs(std::move(image_signs))
{
    ValidateImageVectors(m_image_masks, m_image_signs);
    const int n_cell = n_el;
    const int n_bf   = (int)face_meas.size();
    m_n = n_cell + n_bf;
    m_cent.assign((size_t)m_n * 3, 0.0);
    m_meas.assign((size_t)m_n, 0.0);
    m_size.assign((size_t)m_n, 0.0);
    m_qp.resize(m_n); m_qw.resize(m_n); m_srcTris.resize(m_n);
    if (m_far_quad > 0) { m_qpf.resize(m_n); m_qwf.resize(m_n); }   // low-order FAR rule on the sub-tets/sub-tris
    // degree-2 symmetric rules (same as the tet/tri ctor): 4-pt tet (barycentric), 3-pt tri.
    const double ta = 0.5854101966249685, tb = 0.1381966011250105;
    const double TETF[4][4] = {{ta,tb,tb,tb},{tb,ta,tb,tb},{tb,tb,ta,tb},{tb,tb,tb,ta}};
    const double TRIF[3][3] = {{2.0/3,1.0/6,1.0/6},{1.0/6,2.0/3,1.0/6},{1.0/6,1.0/6,2.0/3}};

    std::vector<double> ref_tet_pts, ref_tet_w;
    rad_gl4_duffy_tet(ref_tet_pts, ref_tet_w);          // 64-node Gauss-Duffy tet rule (shared)
    const int nqt = (int)ref_tet_w.size();

    auto get_tri = [](const std::vector<double>& soup, int t) {  // 9 doubles -> 3x Vec3
        std::array<rad_hdiv::Vec3, 3> T;
        for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) T[i][k] = soup[(size_t)t * 9 + 3 * i + k];
        return T;
    };

    // --- CELLS: centroid-fan outer quad (apex = cell_cent) + store hull tris for PhiAt ---
    for (int c = 0; c < n_cell; ++c) {
        const rad_hdiv::Vec3 cen = {cell_cent[3*c], cell_cent[3*c+1], cell_cent[3*c+2]};
        m_cent[3*c] = cen[0]; m_cent[3*c+1] = cen[1]; m_cent[3*c+2] = cen[2];
        m_meas[c] = cell_meas[c]; m_size[c] = std::cbrt(cell_meas[c]);
        const int t0 = cell_troff[c], t1 = cell_troff[c + 1];
        m_srcTris[c].reserve(t1 - t0);
        m_qp[c].reserve((size_t)(t1 - t0) * nqt);
        m_qw[c].reserve((size_t)(t1 - t0) * nqt);
        for (int t = t0; t < t1; ++t) {
            std::array<rad_hdiv::Vec3, 3> T = get_tri(cell_tris, t);
            m_srcTris[c].push_back(T);
            // sub-tet (cen, T0, T1, T2): tvol = |det([T0-cen, T1-cen, T2-cen])| / 6
            double a1[3], a2[3], a3[3];
            for (int k = 0; k < 3; ++k) { a1[k] = T[0][k]-cen[k]; a2[k] = T[1][k]-cen[k]; a3[k] = T[2][k]-cen[k]; }
            double cr[3] = {a2[1]*a3[2]-a2[2]*a3[1], a2[2]*a3[0]-a2[0]*a3[2], a2[0]*a3[1]-a2[1]*a3[0]};
            const double det6 = std::fabs(a1[0]*cr[0] + a1[1]*cr[1] + a1[2]*cr[2]);   // 6*tvol = |J|
            for (int q = 0; q < nqt; ++q) {
                const double l1 = ref_tet_pts[3*q], l2 = ref_tet_pts[3*q+1], l3 = ref_tet_pts[3*q+2];
                rad_hdiv::Vec3 P;
                for (int k = 0; k < 3; ++k) P[k] = cen[k] + l1*a1[k] + l2*a2[k] + l3*a3[k];
                m_qp[c].push_back(P);
                m_qw[c].push_back(ref_tet_w[q] * det6);     // phys weight = ref_w * |J|, sum over sub-tets = vol
            }
            if (m_far_quad > 0) {           // 4-pt degree-2 FAR rule on this sub-tet (cen,T0,T1,T2); w=tvol/4
                const double tvol = det6 / 6.0;
                for (int q = 0; q < 4; ++q) {
                    rad_hdiv::Vec3 P;
                    for (int k = 0; k < 3; ++k)
                        P[k] = TETF[q][0]*cen[k] + TETF[q][1]*T[0][k] + TETF[q][2]*T[1][k] + TETF[q][3]*T[2][k];
                    m_qpf[c].push_back(P); m_qwf[c].push_back(0.25 * tvol);
                }
            }
        }
    }
    // --- FACES: Dunavant-5 per sub-triangle + store sub-tris for PhiAt ---
    for (int b = 0; b < n_bf; ++b) {
        const int a = n_cell + b;
        m_cent[3*a] = face_cent[3*b]; m_cent[3*a+1] = face_cent[3*b+1]; m_cent[3*a+2] = face_cent[3*b+2];
        m_meas[a] = face_meas[b]; m_size[a] = std::sqrt(face_meas[b]);
        const int t0 = face_troff[b], t1 = face_troff[b + 1];
        m_srcTris[a].reserve(t1 - t0);
        m_qp[a].reserve((size_t)(t1 - t0) * 7);
        m_qw[a].reserve((size_t)(t1 - t0) * 7);
        for (int t = t0; t < t1; ++t) {
            std::array<rad_hdiv::Vec3, 3> T = get_tri(face_tris, t);
            m_srcTris[a].push_back(T);
            double e1[3], e2[3];
            for (int k = 0; k < 3; ++k) { e1[k] = T[1][k]-T[0][k]; e2[k] = T[2][k]-T[0][k]; }
            double cr[3] = {e1[1]*e2[2]-e1[2]*e2[1], e1[2]*e2[0]-e1[0]*e2[2], e1[0]*e2[1]-e1[1]*e2[0]};
            const double area = 0.5 * std::sqrt(cr[0]*cr[0] + cr[1]*cr[1] + cr[2]*cr[2]);
            for (int q = 0; q < 7; ++q) {
                rad_hdiv::Vec3 P = {0, 0, 0};
                for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) P[k] += RAD_DUN5[q][i] * T[i][k];
                m_qp[a].push_back(P);
                m_qw[a].push_back(RAD_DUN5[q][3] * area);
            }
            if (m_far_quad > 0) {           // 3-pt degree-2 FAR rule on this sub-triangle; w=area/3
                for (int q = 0; q < 3; ++q) {
                    rad_hdiv::Vec3 P = {0, 0, 0};
                    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) P[k] += TRIF[q][i] * T[i][k];
                    m_qpf[a].push_back(P); m_qwf[a].push_back(area / 3.0);
                }
            }
        }
    }
}

// CURVED POLYTOPE constructor (FULLY curved): curved CELL volume charge (sub-tets, CurvedTetMapMeasure outer
// quad + CurvedTetPotential in PhiAt) + curved FACE surface charge (sub-tris, CurvedTriMapMeasure +
// CurvedTriPotential).  The cell volume charge is DOMINANT (the lowest-order curved charge cannot represent uniform M exactly,
// div M != 0), so the cell MUST be curved.  cell_curved_nodes [n_cell_subtet*30] = 10 P2 nodes/sub-tet,
// cell_subtet_off [n_cell+1] CSR; ditto face_curved_nodes [n_bf_subtri*18] + face_subtri_off [n_bf+1].
RadHACApKChargeGram::RadHACApKChargeGram(
    std::vector<double> cell_curved_nodes, std::vector<int> cell_subtet_off,
    std::vector<double> cell_cent, std::vector<double> cell_meas,
    std::vector<double> face_curved_nodes, std::vector<int> face_subtri_off,
    std::vector<double> face_cent, std::vector<double> face_meas,
    std::vector<double> ref_tet_pts, std::vector<double> ref_tet_w,
    std::vector<double> ref_tri_pts, std::vector<double> ref_tri_w,
    std::vector<double> curve_gl, std::vector<double> curve_gw, int n_el)
    : m_n_el(n_el), m_analytic(true), m_near_factor(1e30), m_far_quad(0), m_polytope(true),
      m_gl(std::move(curve_gl)), m_gw(std::move(curve_gw))
{
    m_curved_face = true;
    const int n_cell = n_el;
    const int n_bf   = (int)face_meas.size();
    m_n = n_cell + n_bf;
    m_cent.assign((size_t)m_n * 3, 0.0);
    m_meas.assign((size_t)m_n, 0.0);
    m_size.assign((size_t)m_n, 0.0);
    m_qp.resize(m_n); m_qw.resize(m_n);
    m_srcCurvedTets.resize(m_n); m_srcCurvedTris.resize(m_n);
    const int nqt = (int)ref_tet_w.size();              // cell outer quad points per curved sub-tet
    const int nqr = (int)ref_tri_w.size();              // face outer quad points per curved sub-tri

    // --- CELLS: curved sub-tets (CurvedTetMapMeasure outer quad + CurvedTetPotential in PhiAt) ---
    for (int c = 0; c < n_cell; ++c) {
        m_cent[3*c] = cell_cent[3*c]; m_cent[3*c+1] = cell_cent[3*c+1]; m_cent[3*c+2] = cell_cent[3*c+2];
        m_meas[c] = cell_meas[c]; m_size[c] = std::cbrt(cell_meas[c]);
        const int t0 = cell_subtet_off[c], t1 = cell_subtet_off[c + 1];
        m_srcCurvedTets[c].reserve(t1 - t0);
        m_qp[c].reserve((size_t)(t1 - t0) * nqt);
        m_qw[c].reserve((size_t)(t1 - t0) * nqt);
        for (int t = t0; t < t1; ++t) {
            std::array<rad_hdiv::Vec3, 10> nd10;             // curved P2 sub-tet nodes [t*30] = 10 Vec3
            for (int i = 0; i < 10; ++i) for (int k = 0; k < 3; ++k)
                nd10[i][k] = cell_curved_nodes[(size_t)t * 30 + 3*i + k];
            m_srcCurvedTets[c].push_back(nd10);
            double nd[10][3];
            for (int i = 0; i < 10; ++i) for (int k = 0; k < 3; ++k) nd[i][k] = nd10[i][k];
            for (int q = 0; q < nqt; ++q) {                  // curved outer quad: CurvedTetMapMeasure at ref pts
                double X[3], dV;
                rad_hdiv::CurvedTetMapMeasure(nd, ref_tet_pts[3*q], ref_tet_pts[3*q+1], ref_tet_pts[3*q+2], X, dV);
                m_qp[c].push_back({ X[0], X[1], X[2] });
                m_qw[c].push_back(ref_tet_w[q] * dV);
            }
        }
    }
    // --- FACES: curved sub-tris (CurvedTriMapMeasure outer quad + CurvedTriPotential in PhiAt) ---
    for (int b = 0; b < n_bf; ++b) {
        const int a = n_cell + b;
        m_cent[3*a] = face_cent[3*b]; m_cent[3*a+1] = face_cent[3*b+1]; m_cent[3*a+2] = face_cent[3*b+2];
        m_meas[a] = face_meas[b]; m_size[a] = std::sqrt(face_meas[b]);
        const int t0 = face_subtri_off[b], t1 = face_subtri_off[b + 1];
        m_srcCurvedTris[a].reserve(t1 - t0);
        m_qp[a].reserve((size_t)(t1 - t0) * nqr);
        m_qw[a].reserve((size_t)(t1 - t0) * nqr);
        for (int t = t0; t < t1; ++t) {
            std::array<rad_hdiv::Vec3, 6> nd6;               // curved P2 sub-tri nodes [t*18] = 6 Vec3
            for (int i = 0; i < 6; ++i) for (int k = 0; k < 3; ++k)
                nd6[i][k] = face_curved_nodes[(size_t)t * 18 + 3*i + k];
            m_srcCurvedTris[a].push_back(nd6);
            double nd[6][3];
            for (int i = 0; i < 6; ++i) for (int k = 0; k < 3; ++k) nd[i][k] = nd6[i][k];
            for (int q = 0; q < nqr; ++q) {                  // curved outer quad: CurvedTriMapMeasure at ref pts
                double X[3], dA;
                rad_hdiv::CurvedTriMapMeasure(nd, ref_tri_pts[2*q], ref_tri_pts[2*q+1], X, dA);
                m_qp[a].push_back({ X[0], X[1], X[2] });
                m_qw[a].push_back(ref_tri_w[q] * dA);
            }
        }
    }
}

// GLOBAL monotonic build-id source for the QuadDot thread_local memo.  MUST be shared by EVERY
// RadHACApKChargeGram constructor: the memo (cache_owner) is a single function-local thread_local in
// QuadDot, so two constructors with INDEPENDENT counters would hand out colliding ids (each starting at 0)
// -> a high-order build (counter A, id 1) followed by a curved build (counter B, id 1) would NOT clear the
// thread_local cache (cache_owner == m_build_id) and the curved build would reuse the high-order build's
// stale PhiInner values -> nondeterministic per-element corruption under a warm (shared-TaskManager-region)
// threadpool.  A single global counter guarantees strictly-increasing, never-reused ids across all builds.
static long long NextChargeGramBuildId()
{
    static std::atomic<long long> s_id{0};
    return s_id.fetch_add(1) + 1;
}

static bool HexCacheStatsEnabledByEnv()
{
    const char* v = std::getenv("RADIA_HDIV_BLOCK_CACHE_STATS");
    if (!v || v[0] == '\0') v = std::getenv("RADIA_HDIV_HEX_CACHE_STATS");
    return v && v[0] != '\0' && v[0] != '0';
}

static inline void HexStatAdd(bool enabled, std::atomic<long long>& v)
{
    if (enabled) v.fetch_add(1, std::memory_order_relaxed);
}

static size_t HexBlockCacheLimit()
{
    static const size_t limit = []() -> size_t {
        const char* v = std::getenv("RADIA_HDIV_HEX_BLOCK_CACHE_LIMIT");
        if (!v || v[0] == '\0') return 200000u;
        const long long parsed = std::atoll(v);
        return parsed > 0 ? (size_t)parsed : 200000u;
    }();
    return limit;
}

static double HexFarOneSidedThreshold()
{
    static const double threshold = []() -> double {
        const char* v = std::getenv("RADIA_HDIV_HEX_FAR_ONESIDED");
        if (!v || v[0] == '\0') return 0.0;
        const double parsed = std::atof(v);
        return parsed >= 0.0 ? parsed : 0.0;
    }();
    return threshold;
}

static double WedgeFarOneSidedThreshold()
{
    static const double threshold = []() -> double {
        const char* v = std::getenv("RADIA_HDIV_WEDGE_FAR_ONESIDED");
        if (!v || v[0] == '\0') return 0.0;
        const double parsed = std::atof(v);
        return parsed >= 0.0 ? parsed : 0.0;
    }();
    return threshold;
}

// DISTORTED-pair far switch (C-4 fill speedup, 2026-08-09): a well-separated host pair whose geometry is
// NOT affine (Sculpt skin hexes, curved Q2 cells) still has a SMOOTH integrand, and the tensor-product
// far rule is geometry-map exact (Q2 point placement; Piola reference charge measure -- no Jacobian
// appears), so it applies verbatim.  Pairs farther than factor*(size_a+size_b) route to
// QuadBlockHexAffineFarProduct instead of the 6x6-sub graded machinery.  Default matches the accepted
// affine far gate (HEX_AFFINE_EXACT_NEAR_FACTOR = 1.0); <= 0 disables the switch (diagnostic A/B /
// regression-triage escape, same pattern as the one-sided thresholds above).
static double HexDistortedFarFactor()
{
    static const double factor = []() -> double {
        const char* v = std::getenv("RADIA_HDIV_HEX_DISTORTED_FAR_FACTOR");
        if (!v || v[0] == '\0') return 1.0;
        return std::atof(v);
    }();
    return factor;
}

// 0/off: disable wedge translation cache, 1: conservative cell-cell subset, 2/all/default: all translated hosts.
static int WedgeTransCacheScope()
{
    static const int scope = []() -> int {
        const char* v = std::getenv("RADIA_HDIV_WEDGE_TRANS_CACHE");
        if (!v || v[0] == '\0') return 2;
        if (v[0] == '0') return 0;
        if (v[0] == '1') return 1;
        if (v[0] == '2' || v[0] == 'a' || v[0] == 'A') return 2;
        return 2;
    }();
    return scope;
}

static bool HOFarOneSidedEnabled()
{
    static const bool enabled = []() -> bool {
        const char* v = std::getenv("RADIA_HDIV_HO_FAR_ONESIDED");
        return v && v[0] != '\0' && v[0] != '0';
    }();
    return enabled;
}

static bool HOAnalyticBlockEnabled()
{
    static const bool enabled = []() -> bool {
        const char* v = std::getenv("RADIA_HDIV_DISABLE_HO_ANALYTIC_BLOCK");
        return !v || v[0] == '\0' || v[0] == '0';
    }();
    return enabled;
}

static void ValidateImageVectors(const std::vector<int>& image_masks,
                                 const std::vector<double>& image_signs)
{
    if (image_masks.size() != image_signs.size()) {
        throw std::invalid_argument(
            "ChargeGram image_masks and image_signs must have the same length");
    }
}

void RadHACApKChargeGram::SetImageRotations(std::vector<double> angles)
{
    if (angles.empty()) { m_image_rot_angle.clear(); return; }
    if (angles.size() != m_image_masks.size())
        throw std::invalid_argument(
            "ChargeGram image rotation angles must have the same length as image_masks");
    for (size_t i = 0; i < angles.size(); ++i) {
        if (!std::isfinite(angles[i]))
            throw std::invalid_argument("ChargeGram image rotation angles must be finite");
        if (m_image_masks[i] == 0 && angles[i] == 0.0)
            throw std::invalid_argument(
                "ChargeGram: an image with mask 0 and rotation angle 0 is the IDENTITY -- it would "
                "double the direct term; give it a mirror mask or a non-zero rotation angle");
    }
    // INVERSE CLOSURE.  The Gram symmetrizes each image with its transpose, i.e. it evaluates
    // 0.5*(G_{T_i} + G_{T_i^-1}); the total only reproduces Sum_i s_i G_{T_i} when every image has its
    // INVERSE in the list carrying the SAME sign.  A full cyclic group satisfies this automatically
    // (theta_k pairs with theta_{N-k}, and the alternating pattern (-1)^k matches for even N).  A partial
    // or one-sided rotation list would be silently symmetrized into different physics -- reject it.
    const size_t n = angles.size();
    auto matrix = [&](size_t i, double m[9]) {
        const int mask = m_image_masks[i];
        const double mx = (mask & 1) ? -1.0 : 1.0;
        const double my = (mask & 2) ? -1.0 : 1.0;
        const double mz = (mask & 4) ? -1.0 : 1.0;
        const double c = std::cos(angles[i]), sn = std::sin(angles[i]);
        m[0] = mx*c; m[1] = -mx*sn; m[2] = 0.0;
        m[3] = my*sn; m[4] = my*c;  m[5] = 0.0;
        m[6] = 0.0;   m[7] = 0.0;   m[8] = mz;
    };
    for (size_t i = 0; i < n; ++i) {
        double a[9]; matrix(i, a);
        bool found = false;
        for (size_t j = 0; j < n && !found; ++j) {
            if (m_image_signs[j] != m_image_signs[i]) continue;
            double b[9]; matrix(j, b);
            double err = 0.0;
            for (int r = 0; r < 3; ++r)
                for (int c2 = 0; c2 < 3; ++c2)
                    err = std::max(err, std::fabs(b[3*c2 + r] - a[3*r + c2]));   // b == a^T ?
            found = err < 1e-12;
        }
        if (!found)
            throw std::invalid_argument(
                "ChargeGram image rotations must be closed under inversion with matching signs "
                "(image " + std::to_string(i) + " has no inverse partner in the list); pass the "
                "COMPLETE cyclic group, e.g. theta_k = 2*pi*k/N for k = 1..N-1");
    }
    m_image_rot_angle = std::move(angles);
}

// HIGH-ORDER constructor: polynomial charges (monomial basis per host).  See the header for the contract.
RadHACApKChargeGram::RadHACApKChargeGram(
    std::vector<double> cell_verts, std::vector<double> face_verts, int n_el,
    std::vector<int> charge_host, std::vector<int> charge_kind, std::vector<int> charge_expo,
    std::vector<double> ref_tet_pts, std::vector<double> ref_tet_w,
    std::vector<double> ref_tri_pts, std::vector<double> ref_tri_w,
    std::vector<double> ref_tet_pts_lo, std::vector<double> ref_tet_w_lo,
    std::vector<double> ref_tri_pts_lo, std::vector<double> ref_tri_w_lo,
    double ho_far_factor,
    std::vector<double> ref_tet_pts_in, std::vector<double> ref_tet_w_in,
    std::vector<double> ref_tri_pts_in, std::vector<double> ref_tri_w_in,
    std::vector<int> image_masks, std::vector<double> image_signs)
    : m_n_el(n_el), m_highorder(true), m_ho_far_factor(ho_far_factor),
      m_cellV(std::move(cell_verts)), m_faceV(std::move(face_verts)),
      m_image_masks(std::move(image_masks)), m_image_signs(std::move(image_signs)),
      m_host(std::move(charge_host)), m_kind(std::move(charge_kind)), m_expo(std::move(charge_expo))
{
    ValidateImageVectors(m_image_masks, m_image_signs);
    m_hexCacheStatsEnabled = HexCacheStatsEnabledByEnv();
    const int n_cell = n_el;
    const int n_bf   = (int)(m_faceV.size() / 9);
    m_n = (int)m_host.size();                       // number of polynomial CHARGES (the H-matrix dofs)
    m_build_id = NextChargeGramBuildId();           // GLOBAL unique id for the QuadDot memo (see NextChargeGramBuildId)
    // per-(kind,host) co-located charge count -> the QuadDot memo engages only where n_mono>1 (reuse exists);
    // skips e.g. p=1 volume (1 monomial/cell) so the cache never adds overhead where there is nothing to reuse.
    m_nmono.assign(m_n, 1);
    {
        std::unordered_map<long long, int> cnt;
        for (int a = 0; a < m_n; ++a) cnt[(long long)m_host[a]*2 + m_kind[a]]++;
        for (int a = 0; a < m_n; ++a) m_nmono[a] = cnt[(long long)m_host[a]*2 + m_kind[a]];
    }
    m_hoLocalOf.assign((size_t)m_n, 0);
    m_hoCellCharges.assign((size_t)n_cell, {});
    m_hoFaceCharges.assign((size_t)n_bf, {});
    for (int a = 0; a < m_n; ++a) {
        std::vector<int>& group = (m_kind[a] == 0) ? m_hoCellCharges[m_host[a]]
                                                   : m_hoFaceCharges[m_host[a]];
        m_hoLocalOf[a] = (int)group.size();
        group.push_back(a);
    }
    const int nqt = (int)ref_tet_w.size();
    const int nqr = (int)ref_tri_w.size();

    // per-CELL host: ref->phys affine inverse + mapped quadrature (outer & inner share this rule)
    m_cellInv.assign((size_t)n_cell * 9, 0.0);
    std::vector<std::vector<rad_hdiv::Vec3>> cellQP(n_cell);
    std::vector<std::vector<double>>          cellQW(n_cell);
    std::vector<rad_hdiv::Vec3> cellCent(n_cell);
    std::vector<double>         cellSize(n_cell);
    for (int c = 0; c < n_cell; ++c) {
        const double* V = &m_cellV[(size_t)c * 12];
        double E[9];                                // E[r*3+col] = e_{col}[r] = V[col+1][r]-V[0][r]
        for (int r = 0; r < 3; ++r) for (int col = 0; col < 3; ++col) E[r*3+col] = V[3*(col+1)+r] - V[r];
        rad_inv3x3(E, &m_cellInv[(size_t)c*9]);
        const double det = E[0]*(E[4]*E[8]-E[5]*E[7]) - E[1]*(E[3]*E[8]-E[5]*E[6]) + E[2]*(E[3]*E[7]-E[4]*E[6]);
        rad_hdiv::Vec3 cen = {0, 0, 0};
        for (int i = 0; i < 4; ++i) for (int k = 0; k < 3; ++k) cen[k] += V[3*i+k] / 4.0;
        cellCent[c] = cen;
        // FAR/NEAR size = bounding radius (max centroid->vertex distance), NOT cbrt(vol): the isotropic
        // cbrt(vol) UNDERESTIMATES the extent of high-aspect-ratio (needle/sliver) tets, so a TOUCHING pair
        // could satisfy r > ho_far_factor*(size_a+size_b) and be misclassified FAR -> routed to the
        // subtraction-free QuadDotFar on a near-SINGULAR integrand (wrong by ~1-5%, growing with aspect
        // ratio).  The bounding radius captures the long extent so touching pairs always stay NEAR.
        double rmax = 0.0;
        for (int i = 0; i < 4; ++i) {
            const double dvx = V[3*i] - cen[0], dvy = V[3*i+1] - cen[1], dvz = V[3*i+2] - cen[2];
            const double rr = std::sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
            if (rr > rmax) rmax = rr;
        }
        cellSize[c] = rmax;
        cellQP[c].resize(nqt); cellQW[c].resize(nqt);
        for (int q = 0; q < nqt; ++q) {
            const double a = ref_tet_pts[3*q], b = ref_tet_pts[3*q+1], cc = ref_tet_pts[3*q+2];
            rad_hdiv::Vec3 P;
            for (int k = 0; k < 3; ++k) P[k] = V[k] + a*(V[3+k]-V[k]) + b*(V[6+k]-V[k]) + cc*(V[9+k]-V[k]);
            cellQP[c][q] = P;
            cellQW[c][q] = ref_tet_w[q] * std::fabs(det);   // phys weight = ref_w * |J|, |J| = det = 6*vol
        }
    }
    // per-FACE host: 2x2 in-plane Gram inverse + mapped quadrature
    m_faceGinv.assign((size_t)n_bf * 4, 0.0);
    std::vector<std::vector<rad_hdiv::Vec3>> faceQP(n_bf);
    std::vector<std::vector<double>>          faceQW(n_bf);
    std::vector<rad_hdiv::Vec3> faceCent(n_bf);
    std::vector<double>         faceSize(n_bf);
    for (int f = 0; f < n_bf; ++f) {
        const double* V = &m_faceV[(size_t)f * 9];
        double a1[3], a2[3];
        for (int k = 0; k < 3; ++k) { a1[k] = V[3+k]-V[k]; a2[k] = V[6+k]-V[k]; }
        const double a1a2 = a1[0]*a2[0]+a1[1]*a2[1]+a1[2]*a2[2];
        double g[4] = { a1[0]*a1[0]+a1[1]*a1[1]+a1[2]*a1[2], a1a2,
                        a1a2,                                a2[0]*a2[0]+a2[1]*a2[1]+a2[2]*a2[2] };
        rad_inv2x2(g, &m_faceGinv[(size_t)f*4]);
        double cr[3] = {a1[1]*a2[2]-a1[2]*a2[1], a1[2]*a2[0]-a1[0]*a2[2], a1[0]*a2[1]-a1[1]*a2[0]};
        const double area = 0.5 * std::sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);
        rad_hdiv::Vec3 cen = {0, 0, 0};
        for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) cen[k] += V[3*i+k] / 3.0;
        faceCent[f] = cen;
        // FAR/NEAR size = bounding radius (max centroid->vertex distance), NOT sqrt(area) -- same reason as
        // the cell loop: a thin/elongated boundary face's sqrt(area) underestimates its extent.
        double rmax = 0.0;
        for (int i = 0; i < 3; ++i) {
            const double dvx = V[3*i] - cen[0], dvy = V[3*i+1] - cen[1], dvz = V[3*i+2] - cen[2];
            const double rr = std::sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
            if (rr > rmax) rmax = rr;
        }
        faceSize[f] = rmax;
        faceQP[f].resize(nqr); faceQW[f].resize(nqr);
        for (int q = 0; q < nqr; ++q) {
            const double u = ref_tri_pts[2*q], v = ref_tri_pts[2*q+1];
            rad_hdiv::Vec3 P;
            for (int k = 0; k < 3; ++k) P[k] = V[k] + u*a1[k] + v*a2[k];
            faceQP[f][q] = P;
            faceQW[f][q] = ref_tri_w[q] * (2.0 * area);     // phys weight = ref_w * |J|, |J| = 2*area
        }
    }
    InitHOPolynomialCoefficients();
    // INNER subtraction rule (B2 speedup): the subtraction remainder (m_src(y)-m_src(p)) is SMOOTH (the
    // singular part is carried EXACTLY by base = m_src(p)*PhiTet/TriPotential), so the inner sum tolerates a
    // COARSER rule than the outer (which must resolve the degree-p target monomial folded into m_qw).  When
    // the caller supplies ref_*_in, m_inP/m_inW use it; else they fall back to the outer rule (inner=outer).
    const int nqt_in = (int)ref_tet_w_in.size();
    const int nqr_in = (int)ref_tri_w_in.size();
    const bool use_inner = (nqt_in > 0 && nqr_in > 0);
    std::vector<std::vector<rad_hdiv::Vec3>> cellQP_in, faceQP_in;
    std::vector<std::vector<double>>          cellQW_in, faceQW_in;
    if (use_inner) {
        cellQP_in.resize(n_cell); cellQW_in.resize(n_cell);
        for (int c = 0; c < n_cell; ++c) {
            const double* V = &m_cellV[(size_t)c * 12];
            double E[9];
            for (int r = 0; r < 3; ++r) for (int col = 0; col < 3; ++col) E[r*3+col] = V[3*(col+1)+r] - V[r];
            const double det = E[0]*(E[4]*E[8]-E[5]*E[7]) - E[1]*(E[3]*E[8]-E[5]*E[6]) + E[2]*(E[3]*E[7]-E[4]*E[6]);
            cellQP_in[c].resize(nqt_in); cellQW_in[c].resize(nqt_in);
            for (int q = 0; q < nqt_in; ++q) {
                const double a = ref_tet_pts_in[3*q], b = ref_tet_pts_in[3*q+1], cc = ref_tet_pts_in[3*q+2];
                rad_hdiv::Vec3 P;
                for (int k = 0; k < 3; ++k) P[k] = V[k] + a*(V[3+k]-V[k]) + b*(V[6+k]-V[k]) + cc*(V[9+k]-V[k]);
                cellQP_in[c][q] = P;
                cellQW_in[c][q] = ref_tet_w_in[q] * std::fabs(det);
            }
        }
        faceQP_in.resize(n_bf); faceQW_in.resize(n_bf);
        for (int f = 0; f < n_bf; ++f) {
            const double* V = &m_faceV[(size_t)f * 9];
            double a1[3], a2[3];
            for (int k = 0; k < 3; ++k) { a1[k] = V[3+k]-V[k]; a2[k] = V[6+k]-V[k]; }
            double cr[3] = {a1[1]*a2[2]-a1[2]*a2[1], a1[2]*a2[0]-a1[0]*a2[2], a1[0]*a2[1]-a1[1]*a2[0]};
            const double area = 0.5 * std::sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);
            faceQP_in[f].resize(nqr_in); faceQW_in[f].resize(nqr_in);
            for (int q = 0; q < nqr_in; ++q) {
                const double u = ref_tri_pts_in[2*q], v = ref_tri_pts_in[2*q+1];
                rad_hdiv::Vec3 P;
                for (int k = 0; k < 3; ++k) P[k] = V[k] + u*a1[k] + v*a2[k];
                faceQP_in[f][q] = P;
                faceQW_in[f][q] = ref_tri_w_in[q] * (2.0 * area);
            }
        }
    }

    // per-CHARGE: host geometry + monomial-folded outer weights + the inner subtraction table
    m_cent.assign((size_t)m_n * 3, 0.0);
    m_size.assign((size_t)m_n, 0.0);
    m_qp.resize(m_n); m_qw.resize(m_n);
    m_inP.resize(m_n); m_inW.resize(m_n);
    for (int a = 0; a < m_n; ++a) {
        const int host = m_host[a];
        const std::vector<rad_hdiv::Vec3>& QP = (m_kind[a] == 0) ? cellQP[host] : faceQP[host];
        const std::vector<double>&         QW = (m_kind[a] == 0) ? cellQW[host] : faceQW[host];
        const rad_hdiv::Vec3& cen = (m_kind[a] == 0) ? cellCent[host] : faceCent[host];
        m_cent[3*a] = cen[0]; m_cent[3*a+1] = cen[1]; m_cent[3*a+2] = cen[2];
        m_size[a] = (m_kind[a] == 0) ? cellSize[host] : faceSize[host];
        m_qp[a] = QP;
        m_qw[a].resize(QP.size());
        for (size_t q = 0; q < QP.size(); ++q) {
            const double p[3] = {QP[q][0], QP[q][1], QP[q][2]};
            m_qw[a][q] = QW[q] * EvalMono(a, p);            // fold m_a(x_q) into the outer weight
        }
        if (use_inner) {          // B2: COARSER inner subtraction rule (smooth remainder, separate from outer)
            m_inP[a] = (m_kind[a] == 0) ? cellQP_in[host] : faceQP_in[host];
            m_inW[a] = (m_kind[a] == 0) ? cellQW_in[host] : faceQW_in[host];
        } else {
            m_inP[a] = QP;        // inner = outer (original behavior)
            m_inW[a] = QW;
        }
    }
    // precompute m_src(y_q) at the FIXED inner subtraction points -> bit-exact hoist of EvalMono out of the
    // hot PhiAtHO inner loop (the value depends only on (src,q), not on the outer point nor the tgt monomial,
    // yet was recomputed quad^3 times per entry AND for every entry / co-located monomial sharing the source).
    m_srcval.resize(m_n);
    for (int a = 0; a < m_n; ++a) {
        m_srcval[a].resize(m_inP[a].size());
        for (size_t q = 0; q < m_inP[a].size(); ++q) {
            const double y[3] = {m_inP[a][q][0], m_inP[a][q][1], m_inP[a][q][2]};
            m_srcval[a][q] = EvalMono(a, y);
        }
    }

    // ---- LOW-quad tables for the cheap FAR plain double-Gauss (near/far adaptive quadrature) ----
    // Built only when the caller supplies the LOW reference rules AND a finite far factor; otherwise the
    // far split is disabled and every pair uses the full high-quad subtraction (original behavior).
    const int nqt_lo = (int)ref_tet_w_lo.size();
    const int nqr_lo = (int)ref_tri_w_lo.size();
    if (m_ho_far_factor < 1e29 && nqt_lo > 0 && nqr_lo > 0) {
        std::vector<std::vector<rad_hdiv::Vec3>> cellQP_lo(n_cell), faceQP_lo(n_bf);
        std::vector<std::vector<double>>          cellQW_lo(n_cell), faceQW_lo(n_bf);
        for (int c = 0; c < n_cell; ++c) {
            const double* V = &m_cellV[(size_t)c * 12];
            double E[9];
            for (int r = 0; r < 3; ++r) for (int col = 0; col < 3; ++col) E[r*3+col] = V[3*(col+1)+r] - V[r];
            const double det = E[0]*(E[4]*E[8]-E[5]*E[7]) - E[1]*(E[3]*E[8]-E[5]*E[6]) + E[2]*(E[3]*E[7]-E[4]*E[6]);
            cellQP_lo[c].resize(nqt_lo); cellQW_lo[c].resize(nqt_lo);
            for (int q = 0; q < nqt_lo; ++q) {
                const double a = ref_tet_pts_lo[3*q], b = ref_tet_pts_lo[3*q+1], cc = ref_tet_pts_lo[3*q+2];
                rad_hdiv::Vec3 P;
                for (int k = 0; k < 3; ++k) P[k] = V[k] + a*(V[3+k]-V[k]) + b*(V[6+k]-V[k]) + cc*(V[9+k]-V[k]);
                cellQP_lo[c][q] = P;
                cellQW_lo[c][q] = ref_tet_w_lo[q] * std::fabs(det);
            }
        }
        for (int f = 0; f < n_bf; ++f) {
            const double* V = &m_faceV[(size_t)f * 9];
            double a1[3], a2[3];
            for (int k = 0; k < 3; ++k) { a1[k] = V[3+k]-V[k]; a2[k] = V[6+k]-V[k]; }
            double cr[3] = {a1[1]*a2[2]-a1[2]*a2[1], a1[2]*a2[0]-a1[0]*a2[2], a1[0]*a2[1]-a1[1]*a2[0]};
            const double area = 0.5 * std::sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);
            faceQP_lo[f].resize(nqr_lo); faceQW_lo[f].resize(nqr_lo);
            for (int q = 0; q < nqr_lo; ++q) {
                const double u = ref_tri_pts_lo[2*q], v = ref_tri_pts_lo[2*q+1];
                rad_hdiv::Vec3 P;
                for (int k = 0; k < 3; ++k) P[k] = V[k] + u*a1[k] + v*a2[k];
                faceQP_lo[f][q] = P;
                faceQW_lo[f][q] = ref_tri_w_lo[q] * (2.0 * area);
            }
        }
        m_qp_lo.resize(m_n); m_qw_lo.resize(m_n); m_inP_lo.resize(m_n); m_inW_lo.resize(m_n);
        for (int a = 0; a < m_n; ++a) {
            const int host = m_host[a];
            const std::vector<rad_hdiv::Vec3>& QPl = (m_kind[a] == 0) ? cellQP_lo[host] : faceQP_lo[host];
            const std::vector<double>&         QWl = (m_kind[a] == 0) ? cellQW_lo[host] : faceQW_lo[host];
            m_qp_lo[a] = QPl;
            m_qw_lo[a].resize(QPl.size());
            for (size_t q = 0; q < QPl.size(); ++q) {
                const double p[3] = {QPl[q][0], QPl[q][1], QPl[q][2]};
                m_qw_lo[a][q] = QWl[q] * EvalMono(a, p);   // fold m_a into the LOW outer weight
            }
            m_inP_lo[a] = QPl;       // LOW inner points (plain; m_b evaluated on the fly in QuadDotFar)
            m_inW_lo[a] = QWl;
        }
        m_srcval_lo.resize(m_n);     // precompute m_src(y_q) at the FIXED LOW inner points (bit-exact, for QuadDotFar)
        for (int a = 0; a < m_n; ++a) {
            m_srcval_lo[a].resize(m_inP_lo[a].size());
            for (size_t q = 0; q < m_inP_lo[a].size(); ++q) {
                const double y[3] = {m_inP_lo[a][q][0], m_inP_lo[a][q][1], m_inP_lo[a][q][2]};
                m_srcval_lo[a][q] = EvalMono(a, y);
            }
        }
    } else {
        m_ho_far_factor = 1e30;     // no LOW rule supplied -> disable the far split (every pair NEAR)
    }
}

RadHACApKChargeGram::RadHACApKChargeGram(
    std::vector<double> cell_verts, int n_el,
    std::vector<int> charge_host,
    std::vector<double> polynomial_coefficients,
    std::vector<int> polynomial_exponents,
    std::vector<double> ref_tet_pts,
    std::vector<double> ref_tet_w)
    : m_n_el(n_el), m_highorder(true), m_polyCombo(true),
      m_cellV(std::move(cell_verts)), m_host(std::move(charge_host)),
      m_comboCoeffs(std::move(polynomial_coefficients))
{
    if (n_el <= 0 || m_cellV.size() != static_cast<size_t>(n_el)*12)
        throw std::invalid_argument("HCurl polynomial Gram: invalid tetrahedron array");
    if (polynomial_exponents.empty() || polynomial_exponents.size()%3 != 0)
        throw std::invalid_argument("HCurl polynomial Gram: exponents must have shape (n_monomial,3)");
    m_comboNMono = static_cast<int>(polynomial_exponents.size()/3);
    m_comboExponents.resize(static_cast<size_t>(m_comboNMono));
    for (int m = 0; m < m_comboNMono; ++m) {
        auto& e = m_comboExponents[static_cast<size_t>(m)];
        e = {polynomial_exponents[3*m], polynomial_exponents[3*m+1],
             polynomial_exponents[3*m+2]};
        if (e[0] < 0 || e[1] < 0 || e[2] < 0 || e[0] + e[1] + e[2] > 18)
            throw std::invalid_argument("HCurl polynomial Gram: exponent degree must be in [0,18]");
    }
    m_n = static_cast<int>(m_host.size());
    if (m_n <= 0 || m_comboCoeffs.size() != static_cast<size_t>(m_n)*m_comboNMono)
        throw std::invalid_argument("HCurl polynomial Gram: coefficient array shape mismatch");
    if (ref_tet_pts.empty() || ref_tet_pts.size()%3 != 0 ||
        ref_tet_w.size() != ref_tet_pts.size()/3)
        throw std::invalid_argument("HCurl polynomial Gram: invalid outer tetrahedron rule");
    for (int host : m_host)
        if (host < 0 || host >= n_el)
            throw std::invalid_argument("HCurl polynomial Gram: charge host out of range");
    for (double value : m_comboCoeffs)
        if (!std::isfinite(value))
            throw std::invalid_argument("HCurl polynomial Gram: non-finite coefficient");

    m_build_id = NextChargeGramBuildId();
    m_kind.assign(static_cast<size_t>(m_n), 0);
    m_expo.assign(static_cast<size_t>(m_n)*3, 0);
    m_nmono.assign(static_cast<size_t>(m_n), 1);
    m_hoLocalOf.assign(static_cast<size_t>(m_n), 0);
    m_hoCellCharges.assign(static_cast<size_t>(n_el), {});
    m_hoFaceCharges.clear();
    for (int charge = 0; charge < m_n; ++charge) {
        auto& group = m_hoCellCharges[static_cast<size_t>(m_host[charge])];
        m_hoLocalOf[static_cast<size_t>(charge)] = static_cast<int>(group.size());
        group.push_back(charge);
    }
    for (int charge = 0; charge < m_n; ++charge)
        m_nmono[static_cast<size_t>(charge)] =
            static_cast<int>(m_hoCellCharges[static_cast<size_t>(m_host[charge])].size());

    const int nqt = static_cast<int>(ref_tet_w.size());
    m_cellInv.assign(static_cast<size_t>(n_el)*9, 0.0);
    std::vector<std::vector<rad_hdiv::Vec3>> cell_points(static_cast<size_t>(n_el));
    std::vector<std::vector<double>> cell_weights(static_cast<size_t>(n_el));
    std::vector<rad_hdiv::Vec3> cell_centers(static_cast<size_t>(n_el));
    std::vector<double> cell_sizes(static_cast<size_t>(n_el), 0.0);
    for (int cell = 0; cell < n_el; ++cell) {
        const double* V = &m_cellV[static_cast<size_t>(cell)*12];
        double E[9];
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                E[row*3+col] = V[3*(col+1)+row] - V[row];
        rad_inv3x3(E, &m_cellInv[static_cast<size_t>(cell)*9]);
        const double det = E[0]*(E[4]*E[8]-E[5]*E[7])
                         - E[1]*(E[3]*E[8]-E[5]*E[6])
                         + E[2]*(E[3]*E[7]-E[4]*E[6]);
        if (!std::isfinite(det) || std::fabs(det) < 1e-300)
            throw std::invalid_argument("HCurl polynomial Gram: degenerate tetrahedron");
        rad_hdiv::Vec3 center = {0.0, 0.0, 0.0};
        for (int vertex = 0; vertex < 4; ++vertex)
            for (int k = 0; k < 3; ++k) center[k] += V[3*vertex+k]/4.0;
        cell_centers[static_cast<size_t>(cell)] = center;
        double radius = 0.0;
        for (int vertex = 0; vertex < 4; ++vertex) {
            const double dx = V[3*vertex] - center[0];
            const double dy = V[3*vertex+1] - center[1];
            const double dz = V[3*vertex+2] - center[2];
            radius = std::max(radius, std::sqrt(dx*dx + dy*dy + dz*dz));
        }
        cell_sizes[static_cast<size_t>(cell)] = radius;
        auto& points = cell_points[static_cast<size_t>(cell)];
        auto& weights = cell_weights[static_cast<size_t>(cell)];
        points.resize(static_cast<size_t>(nqt));
        weights.resize(static_cast<size_t>(nqt));
        for (int q = 0; q < nqt; ++q) {
            const double xi0 = ref_tet_pts[3*q];
            const double xi1 = ref_tet_pts[3*q+1];
            const double xi2 = ref_tet_pts[3*q+2];
            for (int k = 0; k < 3; ++k)
                points[static_cast<size_t>(q)][k] = V[k]
                    + xi0*(V[3+k]-V[k]) + xi1*(V[6+k]-V[k])
                    + xi2*(V[9+k]-V[k]);
            weights[static_cast<size_t>(q)] = ref_tet_w[static_cast<size_t>(q)]*std::fabs(det);
        }
    }

    auto polynomial_at_reference = [&](int charge, const double xi[3]) {
        const double* coeff = &m_comboCoeffs[static_cast<size_t>(charge)*m_comboNMono];
        double value = 0.0;
        for (int m = 0; m < m_comboNMono; ++m) {
            const auto& e = m_comboExponents[static_cast<size_t>(m)];
            value += coeff[m] * rad_ipow(xi[0], e[0])
                              * rad_ipow(xi[1], e[1])
                              * rad_ipow(xi[2], e[2]);
        }
        return value;
    };

    m_cent.assign(static_cast<size_t>(m_n)*3, 0.0);
    m_size.assign(static_cast<size_t>(m_n), 0.0);
    m_qp.resize(static_cast<size_t>(m_n));
    m_qw.resize(static_cast<size_t>(m_n));
    for (int charge = 0; charge < m_n; ++charge) {
        const int host = m_host[charge];
        const auto& center = cell_centers[static_cast<size_t>(host)];
        m_cent[3*charge] = center[0];
        m_cent[3*charge+1] = center[1];
        m_cent[3*charge+2] = center[2];
        m_size[static_cast<size_t>(charge)] = cell_sizes[static_cast<size_t>(host)];
        m_qp[static_cast<size_t>(charge)] = cell_points[static_cast<size_t>(host)];
        auto& weighted = m_qw[static_cast<size_t>(charge)];
        weighted.resize(static_cast<size_t>(nqt));
        for (int q = 0; q < nqt; ++q) {
            const double xi[3] = {ref_tet_pts[3*q], ref_tet_pts[3*q+1], ref_tet_pts[3*q+2]};
            weighted[static_cast<size_t>(q)] = cell_weights[static_cast<size_t>(host)][static_cast<size_t>(q)]
                * polynomial_at_reference(charge, xi);
        }
    }
    m_hoAnalyticBlock = true;
    m_ho_far_factor = 1e30;
}

// ---- CURVED HIGH-ORDER (isoparametric P2) constructor: monomial-charge Gram on a mesh.Curve(2) geometry. ----
// Mirrors the flat HO build but uses the curved P2 map + curved measure for the OUTER quad (xi^expo folded at
// the REFERENCE point, no affine inverse) and the curved Duffy for the INNER potential (PhiInner -> PhiAtHO_
// Curved).  No analytic moments / inner-subtraction table.  Well-separated pairs may use a lower curved
// double-Gauss rule, matching the flat high-order near/far contract without crossing back into Python.
RadHACApKChargeGram::RadHACApKChargeGram(
    std::vector<double> cell_nodes, std::vector<double> face_nodes,
    std::vector<int> cell_vertices, std::vector<int> face_vertices,
    int n_el, int curve_order,
    std::vector<int> charge_host, std::vector<int> charge_kind, std::vector<int> charge_expo,
    std::vector<double> ref_tet_pts, std::vector<double> ref_tet_w,
    std::vector<double> ref_tri_pts, std::vector<double> ref_tri_w,
    std::vector<double> curve_gl, std::vector<double> curve_gw,
    std::vector<double> ref_tet_pts_lo, std::vector<double> ref_tet_w_lo,
    std::vector<double> ref_tri_pts_lo, std::vector<double> ref_tri_w_lo,
    double ho_far_factor,
    std::vector<int> image_masks, std::vector<double> image_signs,
    bool reference_density)
    : m_n_el(n_el), m_curved(true), m_curve_order(curve_order),
      m_curvedReferenceDensity(reference_density),
      m_cellNodes(std::move(cell_nodes)), m_faceNodes(std::move(face_nodes)),
      m_cellVertices(std::move(cell_vertices)), m_faceVertices(std::move(face_vertices)),
      m_gl(std::move(curve_gl)), m_gw(std::move(curve_gw)),
      m_highorder(true),
      m_image_masks(std::move(image_masks)), m_image_signs(std::move(image_signs)),
      m_host(std::move(charge_host)), m_kind(std::move(charge_kind)), m_expo(std::move(charge_expo))
{
    ValidateImageVectors(m_image_masks, m_image_signs);
    const int n_cell = n_el;
    const int n_bf   = (int)(m_faceNodes.size() / 18);
    if ((int)m_cellVertices.size() != 4*n_cell || (int)m_faceVertices.size() != 3*n_bf)
        throw std::invalid_argument("curved high-order ChargeGram: cell_vertices/face_vertices size mismatch");
    m_hexCacheStatsEnabled = HexCacheStatsEnabledByEnv();
    m_n = (int)m_host.size();
    m_build_id = NextChargeGramBuildId();           // GLOBAL unique id (shared with the high-order ctor)
    m_nmono.assign(m_n, 1);
    {
        std::unordered_map<long long, int> cnt;
        for (int a = 0; a < m_n; ++a) cnt[(long long)m_host[a]*2 + m_kind[a]]++;
        for (int a = 0; a < m_n; ++a) m_nmono[a] = cnt[(long long)m_host[a]*2 + m_kind[a]];
    }
    m_hoLocalOf.assign((size_t)m_n, 0);
    m_hoCellCharges.assign((size_t)n_cell, {});
    m_hoFaceCharges.assign((size_t)n_bf, {});
    for (int a = 0; a < m_n; ++a) {
        std::vector<int>& group = (m_kind[a] == 0) ? m_hoCellCharges[m_host[a]]
                                                   : m_hoFaceCharges[m_host[a]];
        m_hoLocalOf[a] = (int)group.size();
        group.push_back(a);
    }
    const int nqt = (int)ref_tet_w.size();
    const int nqr = (int)ref_tri_w.size();

    // per-HOST curved outer quad: physical points X(xi_q) + curved measure (ref_w * dV/dA, monomial folded
    // per-charge below); centroid + bounding radius from the P2 nodes (cluster-tree point / near-size).
    std::vector<std::vector<rad_hdiv::Vec3>> cellQP(n_cell), faceQP(n_bf);
    std::vector<std::vector<double>>          cellM(n_cell),  faceM(n_bf);
    std::vector<rad_hdiv::Vec3> cellCent(n_cell), faceCent(n_bf);
    std::vector<double>         cellSize(n_cell), faceSize(n_bf);
    for (int c = 0; c < n_cell; ++c) {
        const double (*nd)[3] = (const double(*)[3])&m_cellNodes[(size_t)c*30];
        cellQP[c].resize(nqt); cellM[c].resize(nqt);
        for (int q = 0; q < nqt; ++q) {
            double X[3], dV;
            rad_hdiv::CurvedTetMapMeasure(nd, ref_tet_pts[3*q], ref_tet_pts[3*q+1], ref_tet_pts[3*q+2], X, dV);
            cellQP[c][q] = { X[0], X[1], X[2] };
            cellM[c][q]  = ref_tet_w[q] * (m_curvedReferenceDensity ? 1.0 : dV);
        }
        rad_hdiv::Vec3 cen = {0, 0, 0};
        for (int i = 0; i < 4; ++i) for (int k = 0; k < 3; ++k) cen[k] += nd[i][k] / 4.0;
        cellCent[c] = cen;
        double rmax = 0.0;
        for (int i = 0; i < 10; ++i) {
            const double dx = nd[i][0]-cen[0], dy = nd[i][1]-cen[1], dz = nd[i][2]-cen[2];
            const double rr = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (rr > rmax) rmax = rr;
        }
        cellSize[c] = rmax;
    }
    for (int f = 0; f < n_bf; ++f) {
        const double (*nd)[3] = (const double(*)[3])&m_faceNodes[(size_t)f*18];
        faceQP[f].resize(nqr); faceM[f].resize(nqr);
        for (int q = 0; q < nqr; ++q) {
            double X[3], dA;
            rad_hdiv::CurvedTriMapMeasure(nd, ref_tri_pts[2*q], ref_tri_pts[2*q+1], X, dA);
            faceQP[f][q] = { X[0], X[1], X[2] };
            faceM[f][q]  = ref_tri_w[q] * (m_curvedReferenceDensity ? 1.0 : dA);
        }
        rad_hdiv::Vec3 cen = {0, 0, 0};
        for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) cen[k] += nd[i][k] / 3.0;
        faceCent[f] = cen;
        double rmax = 0.0;
        for (int i = 0; i < 6; ++i) {
            const double dx = nd[i][0]-cen[0], dy = nd[i][1]-cen[1], dz = nd[i][2]-cen[2];
            const double rr = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (rr > rmax) rmax = rr;
        }
        faceSize[f] = rmax;
    }

    // per-CHARGE: outer points = host's curved quad points; weight = host measure * monomial(xi_q) at the
    // REFERENCE point (curved uses the ref pt directly -- no affine inverse / EvalMono).
    m_cent.assign((size_t)m_n*3, 0.0);
    m_size.assign((size_t)m_n, 0.0);
    m_qp.resize(m_n); m_qw.resize(m_n);
    for (int a = 0; a < m_n; ++a) {
        const int host = m_host[a];
        const bool isCell = (m_kind[a] == 0);
        const std::vector<rad_hdiv::Vec3>& QP = isCell ? cellQP[host] : faceQP[host];
        const std::vector<double>&         QM = isCell ? cellM[host]  : faceM[host];
        const rad_hdiv::Vec3& cen = isCell ? cellCent[host] : faceCent[host];
        m_cent[3*a] = cen[0]; m_cent[3*a+1] = cen[1]; m_cent[3*a+2] = cen[2];
        m_size[a] = isCell ? cellSize[host] : faceSize[host];
        const int* e = &m_expo[(size_t)3*a];
        m_qp[a] = QP;
        m_qw[a].resize(QP.size());
        for (int q = 0; q < (int)QP.size(); ++q) {
            double mono;
            if (isCell) {
                mono = rad_ipow(ref_tet_pts[3*q], e[0]) * rad_ipow(ref_tet_pts[3*q+1], e[1])
                     * rad_ipow(ref_tet_pts[3*q+2], e[2]);
            } else {
                mono = rad_ipow(ref_tri_pts[2*q], e[0]) * rad_ipow(ref_tri_pts[2*q+1], e[1]);
            }
            m_qw[a][q] = QM[q] * mono;
        }
    }

    const bool any_low = !ref_tet_pts_lo.empty() || !ref_tet_w_lo.empty() ||
                         !ref_tri_pts_lo.empty() || !ref_tri_w_lo.empty();
    const bool valid_low = !ref_tet_w_lo.empty() && ref_tet_pts_lo.size() == 3*ref_tet_w_lo.size() &&
                           !ref_tri_w_lo.empty() && ref_tri_pts_lo.size() == 2*ref_tri_w_lo.size();
    if (any_low && !valid_low)
        throw std::invalid_argument("curved high-order ChargeGram: incomplete or inconsistent low quadrature");
    if (valid_low) {
        m_ho_far_factor = ho_far_factor;
        const int nqt_lo = (int)ref_tet_w_lo.size();
        const int nqr_lo = (int)ref_tri_w_lo.size();
        std::vector<std::vector<rad_hdiv::Vec3>> cellQP_lo(n_cell), faceQP_lo(n_bf);
        std::vector<std::vector<double>> cellM_lo(n_cell), faceM_lo(n_bf);
        for (int c = 0; c < n_cell; ++c) {
            const double (*nodes)[3] = (const double(*)[3])&m_cellNodes[(size_t)c*30];
            cellQP_lo[c].resize(nqt_lo); cellM_lo[c].resize(nqt_lo);
            for (int q = 0; q < nqt_lo; ++q) {
                double X[3], dV;
                rad_hdiv::CurvedTetMapMeasure(nodes, ref_tet_pts_lo[3*q], ref_tet_pts_lo[3*q+1],
                                              ref_tet_pts_lo[3*q+2], X, dV);
                cellQP_lo[c][q] = {X[0], X[1], X[2]};
                cellM_lo[c][q] = ref_tet_w_lo[q]*(m_curvedReferenceDensity ? 1.0 : dV);
            }
        }
        for (int f = 0; f < n_bf; ++f) {
            const double (*nodes)[3] = (const double(*)[3])&m_faceNodes[(size_t)f*18];
            faceQP_lo[f].resize(nqr_lo); faceM_lo[f].resize(nqr_lo);
            for (int q = 0; q < nqr_lo; ++q) {
                double X[3], dA;
                rad_hdiv::CurvedTriMapMeasure(nodes, ref_tri_pts_lo[2*q], ref_tri_pts_lo[2*q+1], X, dA);
                faceQP_lo[f][q] = {X[0], X[1], X[2]};
                faceM_lo[f][q] = ref_tri_w_lo[q]*(m_curvedReferenceDensity ? 1.0 : dA);
            }
        }
        m_qp_lo.resize(m_n); m_qw_lo.resize(m_n);
        m_inP_lo.resize(m_n); m_inW_lo.resize(m_n); m_srcval_lo.resize(m_n);
        for (int a = 0; a < m_n; ++a) {
            const int host = m_host[a];
            const bool is_cell = m_kind[a] == 0;
            const std::vector<rad_hdiv::Vec3>& points = is_cell ? cellQP_lo[host] : faceQP_lo[host];
            const std::vector<double>& measures = is_cell ? cellM_lo[host] : faceM_lo[host];
            const std::vector<double>& refs = is_cell ? ref_tet_pts_lo : ref_tri_pts_lo;
            const int ref_stride = is_cell ? 3 : 2;
            const int* e = &m_expo[(size_t)3*a];
            m_qp_lo[a] = points;
            m_qw_lo[a].resize(points.size());
            m_srcval_lo[a].resize(points.size());
            for (size_t q = 0; q < points.size(); ++q) {
                double mono = rad_ipow(refs[ref_stride*q], e[0]) * rad_ipow(refs[ref_stride*q+1], e[1]);
                if (is_cell) mono *= rad_ipow(refs[ref_stride*q+2], e[2]);
                m_qw_lo[a][q] = measures[q]*mono;
                m_srcval_lo[a][q] = mono;
            }
            m_inP_lo[a] = points;
            m_inW_lo[a] = measures;
        }
    } else {
        m_ho_far_factor = 1e30;
    }
}

static bool CurvedDirectEnabled()
{
    static const bool enabled = []() {
        const char* value = std::getenv("RADIA_HDIV_CURVED_DIRECT");
        return !value || value[0] == '\0' || value[0] != '0';
    }();
    return enabled;
}

// monomial m_charge at physical point p, via the host's REFERENCE barycentric coords (extrapolates for p
// outside the host -- the subtraction needs m_src(p) at the target's outer points)
double RadHACApKChargeGram::EvalMono(int charge, const double p[3]) const
{
    const int host = m_host[charge];
    const int* e = &m_expo[(size_t)3*charge];
    if (m_kind[charge] == 0) {                              // tet cell: lam1^i lam2^j lam3^k
        const double* V0 = &m_cellV[(size_t)host*12];
        const double* Inv = &m_cellInv[(size_t)host*9];
        const double d[3] = {p[0]-V0[0], p[1]-V0[1], p[2]-V0[2]};
        const double l0 = Inv[0]*d[0]+Inv[1]*d[1]+Inv[2]*d[2];
        const double l1 = Inv[3]*d[0]+Inv[4]*d[1]+Inv[5]*d[2];
        const double l2 = Inv[6]*d[0]+Inv[7]*d[1]+Inv[8]*d[2];
        if (m_polyCombo) {
            const double* coeff = &m_comboCoeffs[(size_t)charge*m_comboNMono];
            double value = 0.0;
            for (int m = 0; m < m_comboNMono; ++m) {
                const auto& ce = m_comboExponents[(size_t)m];
                value += coeff[m] * rad_ipow(l0, ce[0]) * rad_ipow(l1, ce[1])
                                  * rad_ipow(l2, ce[2]);
            }
            return value;
        }
        return rad_ipow(l0, e[0]) * rad_ipow(l1, e[1]) * rad_ipow(l2, e[2]);
    }
    const double* V = &m_faceV[(size_t)host*9];             // tri face: lam1^i lam2^j (in-plane ref coords)
    const double* Gi = &m_faceGinv[(size_t)host*4];
    const double d[3] = {p[0]-V[0], p[1]-V[1], p[2]-V[2]};
    const double a1d = (V[3]-V[0])*d[0]+(V[4]-V[1])*d[1]+(V[5]-V[2])*d[2];
    const double a2d = (V[6]-V[0])*d[0]+(V[7]-V[1])*d[1]+(V[8]-V[2])*d[2];
    const double l0 = Gi[0]*a1d + Gi[1]*a2d;
    const double l1 = Gi[2]*a1d + Gi[3]*a2d;
    return rad_ipow(l0, e[0]) * rad_ipow(l1, e[1]);
}

void RadHACApKChargeGram::InitHOPolynomialCoefficients()
{
    m_hoPolyDegree.assign((size_t)m_n, 0);
    m_hoPolyA.assign((size_t)m_n, 0.0);
    m_hoPolyB.assign((size_t)m_n * 3, 0.0);
    m_hoPolyC.assign((size_t)m_n * 9, 0.0);
    m_hoAnalyticBlock = true;
    bool has_quadratic_face_mode = false;
    for (int src = 0; src < m_n; ++src) {
        const int host = m_host[src];
        const int* e = &m_expo[(size_t)3*src];
        const int deg = e[0] + e[1] + e[2];
        m_hoPolyDegree[src] = deg;
        const int max_degree = (m_kind[src] == 0) ? 1 : 2;
        if (deg > max_degree) {
            m_hoAnalyticBlock = false;
            continue;
        }
        if (m_kind[src] == 1 && deg == 2) has_quadratic_face_mode = true;

        double beta[3][3] = {{0,0,0},{0,0,0},{0,0,0}}, V0[3];
        int ncoord;
        if (m_kind[src] == 0) {
            const double* V = &m_cellV[(size_t)host*12];
            const double* Inv = &m_cellInv[(size_t)host*9];
            for (int i = 0; i < 3; ++i) {
                beta[i][0]=Inv[3*i]; beta[i][1]=Inv[3*i+1]; beta[i][2]=Inv[3*i+2];
            }
            V0[0]=V[0]; V0[1]=V[1]; V0[2]=V[2]; ncoord = 3;
        } else {
            const double* V = &m_faceV[(size_t)host*9];
            const double* Gi = &m_faceGinv[(size_t)host*4];
            double a1[3], a2[3];
            for (int k=0;k<3;++k){ a1[k]=V[3+k]-V[k]; a2[k]=V[6+k]-V[k]; }
            for (int k=0;k<3;++k){
                beta[0][k]=Gi[0]*a1[k]+Gi[1]*a2[k];
                beta[1][k]=Gi[2]*a1[k]+Gi[3]*a2[k];
            }
            V0[0]=V[0]; V0[1]=V[1]; V0[2]=V[2]; ncoord = 2;
        }

        double facA[2], facB[2][3]; int nf = 0;
        for (int i = 0; i < ncoord; ++i) {
            for (int c = 0; c < e[i]; ++c) {
                facB[nf][0]=beta[i][0]; facB[nf][1]=beta[i][1]; facB[nf][2]=beta[i][2];
                facA[nf] = -(beta[i][0]*V0[0]+beta[i][1]*V0[1]+beta[i][2]*V0[2]);
                ++nf;
            }
        }
        double A = 1.0, B[3] = {0,0,0}, C[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        for (int f = 0; f < nf; ++f) {
            const double al = facA[f]; const double* be = facB[f];
            const double nA = A*al; double nB[3], nC[3][3];
            for (int k=0;k<3;++k) nB[k] = A*be[k] + al*B[k];
            for (int k=0;k<3;++k) for (int l=0;l<3;++l)
                nC[k][l] = al*C[k][l] + 0.5*(B[k]*be[l] + be[k]*B[l]);
            A = nA;
            for (int k=0;k<3;++k) {
                B[k]=nB[k];
                for (int l=0;l<3;++l) C[k][l]=nC[k][l];
            }
        }
        m_hoPolyA[src] = A;
        for (int k=0;k<3;++k) {
            m_hoPolyB[(size_t)3*src+k] = B[k];
            for (int l=0;l<3;++l) m_hoPolyC[(size_t)9*src+3*k+l] = C[k][l];
        }
    }
    // BDM1's one/three-charge groups are already cheap on the scalar memo path; whole-host blocks add more
    // cache bookkeeping than kernel work there.  The production win starts at BDM2's six quadratic face modes.
    m_hoAnalyticBlock = m_hoAnalyticBlock && has_quadratic_face_mode;
}

// EXACT analytic high-order inner potential INT_host(src) m_src(y)/|p-y| dy for FLAT panels, charge degree
// <= 2 (the hybrid's machine-precision branch -- replaces the point-subtraction PhiAtHO for order<=2, and is
// EXACT for self/adjacent/far alike, faster than the subtraction since there is NO inner quadrature loop).
// The affine-coord monomial m(y) = prod_i l_i(y)^e_i  (l_i = alpha_i + beta_i . y, beta_i the host
// barycentric gradient) is expanded as a PHYSICAL-coord polynomial A + B.y + y^T C y and contracted with the
// exact moment potentials  INT 1/R, INT y'/R, INT y'(x)y'/R  (rad_hdiv PhiTet/TetMoment1 for cells,
// TriPotential/TriMoment1/TriMoment2 for faces).  Validated to ~1e-14 vs an independent brute-force prototype.
// NOTE: a CELL (volume charge) only ever reaches degree p-1 <= 1 for order<=2, so TetMoment2 is not needed;
// CURVED panels OR tet degree>=2 (order>=3 volume) use the Duffy singular-quadrature path instead (validated
// in the independent Duffy prototypes; Python fail-loud guards order>2 until that path is ported).
double RadHACApKChargeGram::PhiAtHO_Analytic(int src, const double p[3]) const
{
    const int host = m_host[src];
    const int deg = m_hoPolyDegree[src];
    const double A = m_hoPolyA[src];
    const double* B = &m_hoPolyB[(size_t)3*src];
    const double* C = &m_hoPolyC[(size_t)9*src];
    if (m_kind[src] == 0) {                                  // cell: degree <= 1 for order<=2 (no TetMoment2 needed)
        double V[4][3]; const double* s=&m_cellV[(size_t)host*12];
        for (int i=0;i<4;++i) for (int k=0;k<3;++k) V[i][k]=s[3*i+k];
        const double I0 = rad_hdiv::PhiTet(V, p);
        if (deg == 0) return A * I0;
        double M1[3]; rad_hdiv::TetMoment1(V, p, M1);
        return A*I0 + B[0]*M1[0] + B[1]*M1[1] + B[2]*M1[2];
    }
    double V[3][3]; const double* s=&m_faceV[(size_t)host*9];
    for (int i=0;i<3;++i) for (int k=0;k<3;++k) V[i][k]=s[3*i+k];
    const double I0 = rad_hdiv::TriPotential(V, p);
    if (deg == 0) return A * I0;
    double M1[3]; rad_hdiv::TriMoment1(V, p, M1);
    double res = A*I0 + B[0]*M1[0] + B[1]*M1[1] + B[2]*M1[2];
    if (deg >= 2) {
        double M2[3][3]; rad_hdiv::TriMoment2(V, p, M2);
        for (int k=0;k<3;++k) for (int l=0;l<3;++l) res += C[3*k+l]*M2[k][l];
    }
    return res;
}

void RadHACApKChargeGram::PhiInnerHOHostVec(
    int kind, int host, const double p[3], const std::vector<int>& charges, double* values) const
{
    if (charges.empty()) return;
    if (m_polyCombo) {
        if (kind != 0) throw std::logic_error("HCurl polynomial Gram supports cell charges only");
        double V[4][3];
        const double* stored = &m_cellV[(size_t)host*12];
        for (int i = 0; i < 4; ++i)
            for (int k = 0; k < 3; ++k) V[i][k] = stored[3*i+k];
        thread_local std::vector<double> moments;
        moments.resize((size_t)m_comboNMono);
        rad_hdiv::TetReferencePotentialMoments(
            V, p, m_comboExponents, moments.data());
        for (size_t local = 0; local < charges.size(); ++local) {
            const int charge = charges[local];
            const double* coeff = &m_comboCoeffs[(size_t)charge*m_comboNMono];
            double value = 0.0;
            for (int m = 0; m < m_comboNMono; ++m)
                value += coeff[m]*moments[(size_t)m];
            values[local] = value;
        }
        return;
    }
    int max_degree = 0;
    for (int src : charges) max_degree = std::max(max_degree, m_hoPolyDegree[src]);

    double I0, M1[3] = {0,0,0}, M2[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    if (kind == 0) {
        double V[4][3]; const double* stored = &m_cellV[(size_t)host*12];
        for (int i=0;i<4;++i) for (int k=0;k<3;++k) V[i][k]=stored[3*i+k];
        I0 = rad_hdiv::PhiTet(V, p);
        if (max_degree >= 1) rad_hdiv::TetMoment1(V, p, M1);
    } else {
        double V[3][3]; const double* stored = &m_faceV[(size_t)host*9];
        for (int i=0;i<3;++i) for (int k=0;k<3;++k) V[i][k]=stored[3*i+k];
        I0 = rad_hdiv::TriPotential(V, p);
        if (max_degree >= 1) rad_hdiv::TriMoment1(V, p, M1);
        if (max_degree >= 2) rad_hdiv::TriMoment2(V, p, M2);
    }

    for (size_t local = 0; local < charges.size(); ++local) {
        const int src = charges[local];
        const double* B = &m_hoPolyB[(size_t)3*src];
        const double* C = &m_hoPolyC[(size_t)9*src];
        double value = m_hoPolyA[src]*I0;
        if (m_hoPolyDegree[src] >= 1)
            value += B[0]*M1[0] + B[1]*M1[1] + B[2]*M1[2];
        if (m_hoPolyDegree[src] >= 2)
            for (int k=0;k<3;++k) for (int l=0;l<3;++l) value += C[3*k+l]*M2[k][l];
        values[local] = value;
    }
}

// Curved P2 counterpart of PhiInnerHOHostVec.  A scalar PhiAtHO_Curved call repeats the closest-reference
// search, curved map/Jacobian, kernel distance, and Duffy loop for every co-located monomial.  NGSolve-style
// element kernels evaluate geometry once and all local shape functions together; do the same here and return
// the complete source-host vector in one pass.  This is algebraically the same signed reference Duffy rule as
// CurvedTetPotential/CurvedTriPotential, with only the monomial contractions vectorized across local modes.
void RadHACApKChargeGram::PhiInnerHOCurvedHostVec(
    int kind, int host, const double p[3], const std::vector<int>& charges, double* values) const
{
    std::fill(values, values + charges.size(), 0.0);
    if (charges.empty()) return;
    const double* gl = m_gl.data();
    const double* gw = m_gw.data();
    const int nq = (int)m_gl.size();

    if (kind == 0) {
        const double (*nodes)[3] = (const double(*)[3])&m_cellNodes[(size_t)host*30];
        double xi0[3];
        rad_hdiv::ClosestRefTet(nodes, p, xi0);
        static const double C[4][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
        static const int FC[4][3] = {{1,2,3},{0,3,2},{0,1,3},{2,1,0}};
        for (int f = 0; f < 4; ++f) {
            for (int lead = 0; lead < 3; ++lead) {
                const double* b1 = C[FC[f][lead]];
                const double* b2 = C[FC[f][(lead+1)%3]];
                const double* b3 = C[FC[f][(lead+2)%3]];
                double d1[3], d2[3], d3[3], e21[3], e32[3];
                for (int k = 0; k < 3; ++k) {
                    d1[k] = b1[k] - xi0[k]; d2[k] = b2[k] - xi0[k]; d3[k] = b3[k] - xi0[k];
                    e21[k] = b2[k] - b1[k]; e32[k] = b3[k] - b2[k];
                }
                const double cr[3] = {d2[1]*d3[2]-d2[2]*d3[1], d2[2]*d3[0]-d2[0]*d3[2],
                                      d2[0]*d3[1]-d2[1]*d3[0]};
                const double D = d1[0]*cr[0] + d1[1]*cr[1] + d1[2]*cr[2];
                if (std::fabs(D) < 1e-300) continue;
                for (int a = 0; a < nq; ++a) {
                    const double u = gl[a];
                    for (int b = 0; b < nq; ++b) {
                        const double v = gl[b];
                        for (int c = 0; c < nq; ++c) {
                            const double w = gl[c];
                            double z[3];
                            for (int k = 0; k < 3; ++k)
                                z[k] = xi0[k] + u*(d1[k] + v*(e21[k] + w*e32[k]));
                            double X[3], dV;
                            rad_hdiv::CurvedTetMapMeasure(nodes, z[0], z[1], z[2], X, dV);
                            const double dx = p[0]-X[0], dy = p[1]-X[1], dz = p[2]-X[2];
                            const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                            if (r < 1e-300) continue;
                            const double measure = m_curvedReferenceDensity ? 1.0 : dV;
                            const double common = (gw[a]*gw[b]*gw[c]/3.0)*(u*u*v*D)*measure/r;
                            for (size_t local = 0; local < charges.size(); ++local) {
                                const int* e = &m_expo[(size_t)3*charges[local]];
                                values[local] += common * rad_ipow(z[0], e[0]) * rad_ipow(z[1], e[1])
                                                        * rad_ipow(z[2], e[2]);
                            }
                        }
                    }
                }
            }
        }
        return;
    }

    const double (*nodes)[3] = (const double(*)[3])&m_faceNodes[(size_t)host*18];
    double xi0[2];
    rad_hdiv::ClosestRefTri(nodes, p, xi0);
    static const double C[3][2] = {{0,0},{1,0},{0,1}};
    for (int k = 0; k < 3; ++k) {
        const double* A = C[k];
        const double* B = C[(k+1)%3];
        const double e1x = A[0]-xi0[0], e1y = A[1]-xi0[1];
        const double e2x = B[0]-xi0[0], e2y = B[1]-xi0[1];
        const double sgn2 = e1x*e2y - e1y*e2x;
        for (int a = 0; a < nq; ++a) {
            const double u = gl[a];
            for (int b = 0; b < nq; ++b) {
                const double v = gl[b];
                const double xi = xi0[0] + u*e1x + u*v*(e2x-e1x);
                const double eta = xi0[1] + u*e1y + u*v*(e2y-e1y);
                double X[3], dA;
                rad_hdiv::CurvedTriMapMeasure(nodes, xi, eta, X, dA);
                const double dx = p[0]-X[0], dy = p[1]-X[1], dz = p[2]-X[2];
                const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (r < 1e-300) continue;
                const double measure = m_curvedReferenceDensity ? 1.0 : dA;
                const double common = gw[a]*gw[b]*(u*sgn2)*measure/r;
                for (size_t local = 0; local < charges.size(); ++local) {
                    const int* e = &m_expo[(size_t)3*charges[local]];
                    values[local] += common * rad_ipow(xi, e[0]) * rad_ipow(eta, e[1]);
                }
            }
        }
    }
}

bool RadHACApKChargeGram::CurvedHostsTouch(int kindA, int hostA, int kindB, int hostB) const
{
    const int* a = kindA == 0 ? &m_cellVertices[(size_t)4*hostA]
                              : &m_faceVertices[(size_t)3*hostA];
    const int* b = kindB == 0 ? &m_cellVertices[(size_t)4*hostB]
                              : &m_faceVertices[(size_t)3*hostB];
    const int na = kindA == 0 ? 4 : 3;
    const int nb = kindB == 0 ? 4 : 3;
    int shared = 0;
    for (int i = 0; i < na; ++i)
        for (int j = 0; j < nb; ++j)
            if (a[i] == b[j]) { ++shared; break; }
    // A cell-cell vertex-only intersection is a zero-dimensional singular set in the 3-D double integral
    // and the high-order direct product rule resolves it efficiently.  Keep Duffy whenever a boundary-face
    // charge participates: its lower-dimensional support makes a vertex singularity more visible in H.
    return shared >= 2 || (shared == 1 && (kindA == 1 || kindB == 1));
}

void RadHACApKChargeGram::PrecomputeCurvedTouchBlocks()
{
    const auto start = std::chrono::high_resolution_clock::now();
    const int n_cell = (int)(m_cellVertices.size()/4);
    const int n_face = (int)(m_faceVertices.size()/3);
    const int n_host = n_cell + n_face;
    m_curvedTouchBlockIndex.assign((size_t)n_host*n_host, -1);
    std::vector<std::pair<int,int>> pairs;
    for (int ga = 0; ga < n_host; ++ga) {
        const int kindA = ga < n_cell ? 0 : 1;
        const int hostA = ga < n_cell ? ga : ga-n_cell;
        for (int gb = ga; gb < n_host; ++gb) {
            const int kindB = gb < n_cell ? 0 : 1;
            const int hostB = gb < n_cell ? gb : gb-n_cell;
            if (!CurvedHostsTouch(kindA, hostA, kindB, hostB)) continue;
            m_curvedTouchBlockIndex[(size_t)ga*n_host + gb] = (int)pairs.size();
            pairs.emplace_back(ga, gb);
        }
    }
    m_curvedTouchBlocks.clear();
    m_curvedTouchBlocks.resize(pairs.size());
    ngcore::ParallelFor(ngcore::IntRange(pairs.size()), [&](size_t pair_index) {
        const int ga = pairs[pair_index].first, gb = pairs[pair_index].second;
        const int kindA = ga < n_cell ? 0 : 1, hostA = ga < n_cell ? ga : ga-n_cell;
        const int kindB = gb < n_cell ? 0 : 1, hostB = gb < n_cell ? gb : gb-n_cell;
        const int nA = kindA == 0 ? (int)m_hoCellCharges[hostA].size()
                                  : (int)m_hoFaceCharges[hostA].size();
        const int nB = kindB == 0 ? (int)m_hoCellCharges[hostB].size()
                                  : (int)m_hoFaceCharges[hostB].size();
        std::vector<double> ab = QuadBlockHOTet(kindA, hostA, kindB, hostB);
        std::vector<double> sym((size_t)nA*nB);
        if (ga == gb) {
            for (int la = 0; la < nA; ++la)
                for (int lb = 0; lb < nB; ++lb)
                    sym[(size_t)la*nB + lb] =
                        0.5*(ab[(size_t)la*nB + lb] + ab[(size_t)lb*nA + la]);
        } else {
            std::vector<double> ba = QuadBlockHOTet(kindB, hostB, kindA, hostA);
            for (int la = 0; la < nA; ++la)
                for (int lb = 0; lb < nB; ++lb)
                    sym[(size_t)la*nB + lb] =
                        0.5*(ab[(size_t)la*nB + lb] + ba[(size_t)lb*nA + la]);
        }
        m_curvedTouchBlocks[pair_index] = std::move(sym);
    });
    const auto stop = std::chrono::high_resolution_clock::now();
    m_curvedTouchBuildTime = std::chrono::duration<double>(stop-start).count();
}

bool RadHACApKChargeGram::CurvedTouchBlockValue(
    int kindA, int hostA, int localA, int kindB, int hostB, int localB, double& value) const
{
    if (m_curvedTouchBlockIndex.empty()) return false;
    const int n_cell = (int)(m_cellVertices.size()/4);
    const int n_host = n_cell + (int)(m_faceVertices.size()/3);
    const int ga = kindA == 0 ? hostA : n_cell+hostA;
    const int gb = kindB == 0 ? hostB : n_cell+hostB;
    const int lo = std::min(ga, gb), hi = std::max(ga, gb);
    const int slot = m_curvedTouchBlockIndex[(size_t)lo*n_host + hi];
    if (slot < 0) return false;
    const int kind_hi = hi < n_cell ? 0 : 1;
    const int host_hi = hi < n_cell ? hi : hi-n_cell;
    const int n_hi = kind_hi == 0 ? (int)m_hoCellCharges[host_hi].size()
                                  : (int)m_hoFaceCharges[host_hi].size();
    if (ga <= gb) value = m_curvedTouchBlocks[(size_t)slot][(size_t)localA*n_hi + localB];
    else value = m_curvedTouchBlocks[(size_t)slot][(size_t)localB*n_hi + localA];
    return true;
}

// Smooth non-touching curved host pair.  Geometry and curved measures were materialized once by the
// constructor; m_qw already folds each local monomial.  Evaluate the high-order tensor product directly,
// avoiding a closest-point search and a curved Duffy map for every target quadrature point.  Touching pairs
// stay on the singularity-resolving host-vector Duffy path.
std::vector<double> RadHACApKChargeGram::QuadBlockHOCurvedDirect(
    int kindT, int hostT, int kindS, int hostS) const
{
    const std::vector<int>& targets = kindT == 0 ? m_hoCellCharges[hostT] : m_hoFaceCharges[hostT];
    const std::vector<int>& sources = kindS == 0 ? m_hoCellCharges[hostS] : m_hoFaceCharges[hostS];
    const int nT = (int)targets.size(), nS = (int)sources.size();
    std::vector<double> block((size_t)nT*nS, 0.0), inner((size_t)nS, 0.0);
    if (nT == 0 || nS == 0) return block;
    const std::vector<rad_hdiv::Vec3>& target_points = m_qp[targets[0]];
    const std::vector<rad_hdiv::Vec3>& source_points = m_qp[sources[0]];
    for (size_t qt = 0; qt < target_points.size(); ++qt) {
        std::fill(inner.begin(), inner.end(), 0.0);
        const double x0 = target_points[qt][0], x1 = target_points[qt][1], x2 = target_points[qt][2];
        for (size_t qs = 0; qs < source_points.size(); ++qs) {
            const double dx = x0-source_points[qs][0], dy = x1-source_points[qs][1],
                         dz = x2-source_points[qs][2];
            const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (r < 1e-300) continue;
            const double invr = 1.0/r;
            for (int ls = 0; ls < nS; ++ls) inner[ls] += m_qw[sources[ls]][qs]*invr;
        }
        for (int lt = 0; lt < nT; ++lt) {
            const double wt = m_qw[targets[lt]][qt];
            double* row = &block[(size_t)lt*nS];
            for (int ls = 0; ls < nS; ++ls) row[ls] += wt*inner[ls];
        }
    }
    for (double& value : block) value *= RAD_INV_FOUR_PI;
    return block;
}

std::vector<double> RadHACApKChargeGram::QuadBlockHOTet(
    int kindT, int hostT, int kindS, int hostS) const
{
    const std::vector<int>& targets = (kindT == 0) ? m_hoCellCharges[hostT] : m_hoFaceCharges[hostT];
    const std::vector<int>& sources = (kindS == 0) ? m_hoCellCharges[hostS] : m_hoFaceCharges[hostS];
    const int nT = (int)targets.size(), nS = (int)sources.size();
    std::vector<double> block((size_t)nT*nS, 0.0), inner((size_t)nS, 0.0);
    if (nT == 0 || nS == 0) return block;
    if (m_curved && CurvedDirectEnabled() && !CurvedHostsTouch(kindT, hostT, kindS, hostS))
        return QuadBlockHOCurvedDirect(kindT, hostT, kindS, hostS);

    const std::vector<rad_hdiv::Vec3>& points = m_qp[targets[0]];
    for (size_t q = 0; q < points.size(); ++q) {
        const double p[3] = {points[q][0], points[q][1], points[q][2]};
        if (m_curved) PhiInnerHOCurvedHostVec(kindS, hostS, p, sources, inner.data());
        else PhiInnerHOHostVec(kindS, hostS, p, sources, inner.data());
        for (int lt = 0; lt < nT; ++lt) {
            const double weight = m_qw[targets[lt]][q];
            double* row = &block[(size_t)lt*nS];
            for (int ls = 0; ls < nS; ++ls) row[ls] += weight*inner[ls];
        }
    }
    for (double& value : block) value *= RAD_INV_FOUR_PI;
    return block;
}

namespace {
static void Inverse3Directional(const double A[9],const double dA[9],double I[9],double dI[9])
{
    rad_inv3x3(A,I);
    for(int i=0;i<3;++i)for(int j=0;j<3;++j){double s=0;for(int k=0;k<3;++k)for(int l=0;l<3;++l)s+=I[3*i+k]*dA[3*k+l]*I[3*l+j];dI[3*i+j]=-s;}
}
static double Det3Rate(const double A[9],const double dA[9])
{
    const double det=A[0]*(A[4]*A[8]-A[5]*A[7])-A[1]*(A[3]*A[8]-A[5]*A[6])+A[2]*(A[3]*A[7]-A[4]*A[6]);
    const double ddet=
        dA[0]*(A[4]*A[8]-A[5]*A[7])+A[0]*(dA[4]*A[8]+A[4]*dA[8]-dA[5]*A[7]-A[5]*dA[7])
       -dA[1]*(A[3]*A[8]-A[5]*A[6])-A[1]*(dA[3]*A[8]+A[3]*dA[8]-dA[5]*A[6]-A[5]*dA[6])
       +dA[2]*(A[3]*A[7]-A[4]*A[6])+A[2]*(dA[3]*A[7]+A[3]*dA[7]-dA[4]*A[6]-A[4]*dA[6]);
    return ddet/det;
}
static int MomentIndex3(int ax,int ay,int az)
{const int d=ax+ay+az;int n=0;for(int q=0;q<d;++q)n+=(q+1)*(q+2)/2;for(int x=0;x<ax;++x)n+=d-x+1;return n+ay;}
static void MulLinear3(double* p,double* dp,int& degree,const double f[4],const double df[4])
{
    double q[10]={},dq[10]={};
    for(int ax=0;ax<=degree;++ax)for(int ay=0;ay<=degree-ax;++ay){int az=degree-ax-ay; (void)az;}
    for(int total=0;total<=degree;++total)for(int ax=0;ax<=total;++ax)for(int ay=0;ay<=total-ax;++ay){int az=total-ax-ay,idx=MomentIndex3(ax,ay,az);for(int k=0;k<4;++k){int bx=ax+(k==1),by=ay+(k==2),bz=az+(k==3),j=MomentIndex3(bx,by,bz);q[j]+=p[idx]*f[k];dq[j]+=dp[idx]*f[k]+p[idx]*df[k];}}
    ++degree;std::copy(q,q+10,p);std::copy(dq,dq+10,dp);
}
}

std::vector<double> RadHACApKChargeGram::TetVolumeSelfBlockDirectionalDerivative(
    int host,const std::vector<double>& velocity) const
{
    if(!m_highorder||m_curved||m_hexmode||m_polyCombo)throw std::logic_error("TET volume derivative requires a flat polynomial TET charge Gram");
    if(host<0||host>=(int)m_hoCellCharges.size())throw std::out_of_range("TET host out of range");
    if(velocity.size()!=12)throw std::invalid_argument("vertex_velocity must have shape (4,3)");
    const auto& g=m_hoCellCharges[host];const int n=(int)g.size();std::vector<double> out((size_t)n*n),inner(n),innerv(n);
    const double* V0=&m_cellV[(size_t)host*12];double V[4][3],dV[4][3],E[9],dE[9],I[9],dI[9];
    for(int a=0;a<4;++a)for(int k=0;k<3;++k){V[a][k]=V0[3*a+k];dV[a][k]=velocity[3*a+k];}
    for(int k=0;k<3;++k)for(int j=0;j<3;++j){E[3*k+j]=V[j+1][k]-V[0][k];dE[3*k+j]=dV[j+1][k]-dV[0][k];}
    Inverse3Directional(E,dE,I,dI);const double rate=Det3Rate(E,dE);
    for(size_t q=0;q<m_qp[g[0]].size();++q){const auto& pp=m_qp[g[0]][q];double p[3]={pp[0],pp[1],pp[2]},xi[3]={};for(int i=0;i<3;++i)for(int k=0;k<3;++k)xi[i]+=I[3*i+k]*(p[k]-V[0][k]);double dpnt[3];for(int k=0;k<3;++k)dpnt[k]=dV[0][k]+xi[0]*(dV[1][k]-dV[0][k])+xi[1]*(dV[2][k]-dV[0][k])+xi[2]*(dV[3][k]-dV[0][k]);
        double mv[4],dm[4];rad_hdiv::TetPotentialMomentsDirectionalUpTo1(V,dV,p,dpnt,mv,dm);
        // The shared directional kernel returns the physical positive
        // Newtonian moments [1,x,y,z], matching PhiTet and m_qw.
        for(int ls=0;ls<n;++ls){const int* e=&m_expo[(size_t)3*g[ls]];const int deg=e[0]+e[1]+e[2];if(deg>1)throw std::logic_error("analytic TET volume derivative supports charge degree <= 1");double val=mv[0],der=dm[0];if(deg==1){int c=e[1]?1:(e[2]?2:0);double beta[3],dbeta[3],alpha=0,dalpha=0;for(int k=0;k<3;++k){beta[k]=I[3*c+k];dbeta[k]=dI[3*c+k];alpha-=beta[k]*V[0][k];dalpha-=dbeta[k]*V[0][k]+beta[k]*dV[0][k];}val=alpha*mv[0];der=dalpha*mv[0]+alpha*dm[0];for(int k=0;k<3;++k){val+=beta[k]*mv[k+1];der+=dbeta[k]*mv[k+1]+beta[k]*dm[k+1];}}innerv[ls]=val;inner[ls]=der;}
        for(int i=0;i<n;++i){const double w=m_qw[g[i]][q];for(int j=0;j<n;++j)out[(size_t)i*n+j]+=w*(inner[j]+rate*innerv[j]);}
    }
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){double x=.5*(out[(size_t)i*n+j]+out[(size_t)j*n+i]);out[(size_t)i*n+j]=out[(size_t)j*n+i]=x;}
    for(double& x:out)x*=RAD_INV_FOUR_PI;return out;
}

std::vector<double> RadHACApKChargeGram::TetFaceSelfBlockDirectionalDerivative(
    int host,const std::vector<double>& velocity) const
{
    if(!m_highorder||m_curved||m_hexmode||m_polyCombo)throw std::logic_error("TET face derivative requires a flat polynomial TET charge Gram");
    if(host<0||host>=(int)m_hoFaceCharges.size())throw std::out_of_range("TET face host out of range");
    if(velocity.size()!=9)throw std::invalid_argument("vertex_velocity must have shape (3,3)");
    const auto& g=m_hoFaceCharges[host];const int n=(int)g.size();std::vector<double> out((size_t)n*n),inner(n),innerv(n);
    const double* s=&m_faceV[(size_t)host*9];double V[3][3],dV[3][3],a[2][3],da[2][3];for(int i=0;i<3;++i)for(int k=0;k<3;++k){V[i][k]=s[3*i+k];dV[i][k]=velocity[3*i+k];}for(int c=0;c<2;++c)for(int k=0;k<3;++k){a[c][k]=V[c+1][k]-V[0][k];da[c][k]=dV[c+1][k]-dV[0][k];}
    double G[4],dG[4];for(int i=0;i<2;++i)for(int j=0;j<2;++j){G[2*i+j]=dG[2*i+j]=0;for(int k=0;k<3;++k){G[2*i+j]+=a[i][k]*a[j][k];dG[2*i+j]+=da[i][k]*a[j][k]+a[i][k]*da[j][k];}}
    const double det=G[0]*G[3]-G[1]*G[2];double GI[4]={G[3]/det,-G[1]/det,-G[2]/det,G[0]/det},dGI[4];for(int i=0;i<2;++i)for(int j=0;j<2;++j){double z=0;for(int k=0;k<2;++k)for(int l=0;l<2;++l)z+=GI[2*i+k]*dG[2*k+l]*GI[2*l+j];dGI[2*i+j]=-z;}
    double form[2][4]={},dform[2][4]={};for(int c=0;c<2;++c){for(int k=0;k<3;++k){form[c][k+1]=GI[2*c]*a[0][k]+GI[2*c+1]*a[1][k];dform[c][k+1]=dGI[2*c]*a[0][k]+GI[2*c]*da[0][k]+dGI[2*c+1]*a[1][k]+GI[2*c+1]*da[1][k];form[c][0]-=form[c][k+1]*V[0][k];dform[c][0]-=dform[c][k+1]*V[0][k]+form[c][k+1]*dV[0][k];}}
    double cr[3]={a[0][1]*a[1][2]-a[0][2]*a[1][1],a[0][2]*a[1][0]-a[0][0]*a[1][2],a[0][0]*a[1][1]-a[0][1]*a[1][0]},dcr[3]={da[0][1]*a[1][2]+a[0][1]*da[1][2]-da[0][2]*a[1][1]-a[0][2]*da[1][1],da[0][2]*a[1][0]+a[0][2]*da[1][0]-da[0][0]*a[1][2]-a[0][0]*da[1][2],da[0][0]*a[1][1]+a[0][0]*da[1][1]-da[0][1]*a[1][0]-a[0][1]*da[1][0]};double nn=0,cn=0;for(int k=0;k<3;++k){nn+=cr[k]*cr[k];cn+=cr[k]*dcr[k];}const double rate=cn/nn;
    for(size_t q=0;q<m_qp[g[0]].size();++q){const auto& pp=m_qp[g[0]][q];double p[3]={pp[0],pp[1],pp[2]},uv[2]={};for(int c=0;c<2;++c){uv[c]=form[c][0];for(int k=0;k<3;++k)uv[c]+=form[c][k+1]*p[k];}double dpnt[3];for(int k=0;k<3;++k)dpnt[k]=dV[0][k]+uv[0]*da[0][k]+uv[1]*da[1][k];double mv[10],dm[10];rad_hdiv::TriPotentialMomentsDirectionalUpTo2(V,dV,p,dpnt,mv,dm);
        for(int ls=0;ls<n;++ls){const int* e=&m_expo[(size_t)3*g[ls]];if(e[0]+e[1]>2||e[2])throw std::logic_error("analytic TET face derivative supports charge degree <= 2");double poly[10]={1},dpoly[10]={};int degree=0;for(int z=0;z<e[0];++z)MulLinear3(poly,dpoly,degree,form[0],dform[0]);for(int z=0;z<e[1];++z)MulLinear3(poly,dpoly,degree,form[1],dform[1]);double val=0,der=0;for(int total=0;total<=degree;++total)for(int ax=0;ax<=total;++ax)for(int ay=0;ay<=total-ax;++ay){int az=total-ax-ay,id=MomentIndex3(ax,ay,az);val+=poly[id]*mv[id];der+=dpoly[id]*mv[id]+poly[id]*dm[id];}innerv[ls]=val;inner[ls]=der;}
        for(int i=0;i<n;++i){double w=m_qw[g[i]][q];for(int j=0;j<n;++j)out[(size_t)i*n+j]+=w*(inner[j]+rate*innerv[j]);}}
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){double x=.5*(out[(size_t)i*n+j]+out[(size_t)j*n+i]);out[(size_t)i*n+j]=out[(size_t)j*n+i]=x;}for(double& x:out)x*=RAD_INV_FOUR_PI;return out;
}

std::vector<double> RadHACApKChargeGram::TetChargeGramDirectionalDerivative(
    const std::vector<double>& cell_velocity,const std::vector<double>& face_velocity) const
{
    return TetChargeGramDirectionalDerivativeImpl(cell_velocity,face_velocity,-1,-1);
}

std::vector<double> RadHACApKChargeGram::TetChargeGramDirectionalDerivativeImpl(
    const std::vector<double>& cell_velocity,const std::vector<double>& face_velocity,
    int selected_host_a,int selected_host_b) const
{
    if(!m_highorder||m_curved||m_hexmode||m_wedgemode)
        throw std::logic_error("TET ChargeGram derivative requires a flat polynomial TET charge Gram");
    const int nc=(int)m_hoCellCharges.size(),nf=(int)m_hoFaceCharges.size();
    if(cell_velocity.size()!=(size_t)nc*12)throw std::invalid_argument("cell_vertex_velocity must have shape (ncell,4,3)");
    if(face_velocity.size()!=(size_t)nf*9)throw std::invalid_argument("face_vertex_velocity must have shape (nface,3,3)");
    std::vector<double> out;

    // Target affine kinematics at a production outer point.  m_qw already
    // contains its physical measure, hence dw = rate*w.
    auto target_kinematics=[&](int kind,int host,const double p[3],double dp[3],double& rate){
        if(kind==0){
            const double*s=&m_cellV[(size_t)host*12],*ds=&cell_velocity[(size_t)host*12];double E[9],dE[9],I[9],dI[9];
            for(int k=0;k<3;++k)for(int j=0;j<3;++j){E[3*k+j]=s[3*(j+1)+k]-s[k];dE[3*k+j]=ds[3*(j+1)+k]-ds[k];}
            Inverse3Directional(E,dE,I,dI);double xi[3]={};for(int i=0;i<3;++i)for(int k=0;k<3;++k)xi[i]+=I[3*i+k]*(p[k]-s[k]);
            for(int k=0;k<3;++k)dp[k]=ds[k]+xi[0]*(ds[3+k]-ds[k])+xi[1]*(ds[6+k]-ds[k])+xi[2]*(ds[9+k]-ds[k]);
            rate=Det3Rate(E,dE);return;
        }
        const double*s=&m_faceV[(size_t)host*9],*ds=&face_velocity[(size_t)host*9];double a[2][3],da[2][3];
        for(int c=0;c<2;++c)for(int k=0;k<3;++k){a[c][k]=s[3*(c+1)+k]-s[k];da[c][k]=ds[3*(c+1)+k]-ds[k];}
        double G[4]={},dG[4]={};for(int i=0;i<2;++i)for(int j=0;j<2;++j)for(int k=0;k<3;++k){G[2*i+j]+=a[i][k]*a[j][k];dG[2*i+j]+=da[i][k]*a[j][k]+a[i][k]*da[j][k];}
        const double det=G[0]*G[3]-G[1]*G[2];double GI[4]={G[3]/det,-G[1]/det,-G[2]/det,G[0]/det};double uv[2]={};
        for(int c=0;c<2;++c)for(int k=0;k<3;++k)uv[c]+= (GI[2*c]*a[0][k]+GI[2*c+1]*a[1][k])*(p[k]-s[k]);
        for(int k=0;k<3;++k)dp[k]=ds[k]+uv[0]*da[0][k]+uv[1]*da[1][k];
        double cr[3]={a[0][1]*a[1][2]-a[0][2]*a[1][1],a[0][2]*a[1][0]-a[0][0]*a[1][2],a[0][0]*a[1][1]-a[0][1]*a[1][0]},dcr[3]={da[0][1]*a[1][2]+a[0][1]*da[1][2]-da[0][2]*a[1][1]-a[0][2]*da[1][1],da[0][2]*a[1][0]+a[0][2]*da[1][0]-da[0][0]*a[1][2]-a[0][0]*da[1][2],da[0][0]*a[1][1]+a[0][0]*da[1][1]-da[0][1]*a[1][0]-a[0][1]*da[1][0]};double nn=0,cn=0;for(int k=0;k<3;++k){nn+=cr[k]*cr[k];cn+=cr[k]*dcr[k];}rate=cn/nn;
    };

    // Analytic source-host potential and its directional derivative at a
    // moving target point.  This is the exact moment path used by
    // PhiInnerHOHostVec, augmented by forward directional algebra.
    auto source_inner=[&](int kind,int host,const double p[3],const double dp[3],std::vector<double>& val,std::vector<double>& der){
        const auto& g=kind==0?m_hoCellCharges[host]:m_hoFaceCharges[host];val.resize(g.size());der.resize(g.size());
        if(kind==0){
            const double*s=&m_cellV[(size_t)host*12],*ds=&cell_velocity[(size_t)host*12];double V[4][3],dV[4][3],E[9],dE[9],I[9],dI[9];for(int a=0;a<4;++a)for(int k=0;k<3;++k){V[a][k]=s[3*a+k];dV[a][k]=ds[3*a+k];}for(int k=0;k<3;++k)for(int j=0;j<3;++j){E[3*k+j]=V[j+1][k]-V[0][k];dE[3*k+j]=dV[j+1][k]-dV[0][k];}Inverse3Directional(E,dE,I,dI);
            if(m_polyCombo){
                std::vector<double> moments((size_t)m_comboNMono),dmoments((size_t)m_comboNMono);
                rad_hdiv::TetReferencePotentialMomentsDirectional(V,dV,p,dp,m_comboExponents,moments.data(),dmoments.data());
                for(size_t l=0;l<g.size();++l){const double*c=&m_comboCoeffs[(size_t)g[l]*m_comboNMono];val[l]=der[l]=0;for(int m=0;m<m_comboNMono;++m){val[l]+=c[m]*moments[m];der[l]+=c[m]*dmoments[m];}}
                return;
            }
            double mv[4],dm[4];rad_hdiv::TetPotentialMomentsDirectionalUpTo1(V,dV,p,dp,mv,dm);
            for(size_t l=0;l<g.size();++l){const int*e=&m_expo[(size_t)3*g[l]];const int deg=e[0]+e[1]+e[2];if(deg>1)throw std::logic_error("analytic TET volume derivative supports charge degree <= 1");val[l]=mv[0];der[l]=dm[0];if(deg==1){int c=e[1]?1:(e[2]?2:0);double alpha=0,dalpha=0;for(int k=0;k<3;++k){const double beta=I[3*c+k],dbeta=dI[3*c+k];alpha-=beta*V[0][k];dalpha-=dbeta*V[0][k]+beta*dV[0][k];val[l]+=beta*mv[k+1];der[l]+=dbeta*mv[k+1]+beta*dm[k+1];}val[l]+=(alpha-1.0)*mv[0];der[l]+=(dalpha)*mv[0]+(alpha-1.0)*dm[0];}}
            return;
        }
        const double*s=&m_faceV[(size_t)host*9],*ds=&face_velocity[(size_t)host*9];double V[3][3],dV[3][3],a[2][3],da[2][3];for(int i=0;i<3;++i)for(int k=0;k<3;++k){V[i][k]=s[3*i+k];dV[i][k]=ds[3*i+k];}for(int c=0;c<2;++c)for(int k=0;k<3;++k){a[c][k]=V[c+1][k]-V[0][k];da[c][k]=dV[c+1][k]-dV[0][k];}
        double G[4]={},dG[4]={};for(int i=0;i<2;++i)for(int j=0;j<2;++j)for(int k=0;k<3;++k){G[2*i+j]+=a[i][k]*a[j][k];dG[2*i+j]+=da[i][k]*a[j][k]+a[i][k]*da[j][k];}const double det=G[0]*G[3]-G[1]*G[2];double GI[4]={G[3]/det,-G[1]/det,-G[2]/det,G[0]/det},dGI[4];for(int i=0;i<2;++i)for(int j=0;j<2;++j){double z=0;for(int k=0;k<2;++k)for(int l=0;l<2;++l)z+=GI[2*i+k]*dG[2*k+l]*GI[2*l+j];dGI[2*i+j]=-z;}double form[2][4]={},dform[2][4]={};for(int c=0;c<2;++c)for(int k=0;k<3;++k){form[c][k+1]=GI[2*c]*a[0][k]+GI[2*c+1]*a[1][k];dform[c][k+1]=dGI[2*c]*a[0][k]+GI[2*c]*da[0][k]+dGI[2*c+1]*a[1][k]+GI[2*c+1]*da[1][k];form[c][0]-=form[c][k+1]*V[0][k];dform[c][0]-=dform[c][k+1]*V[0][k]+form[c][k+1]*dV[0][k];}
        double mv[10],dm[10];rad_hdiv::TriPotentialMomentsDirectionalUpTo2(V,dV,p,dp,mv,dm);for(size_t l=0;l<g.size();++l){const int*e=&m_expo[(size_t)3*g[l]];if(e[0]+e[1]>2||e[2])throw std::logic_error("analytic TET face derivative supports charge degree <= 2");double poly[10]={1},dpoly[10]={};int degree=0;for(int z=0;z<e[0];++z)MulLinear3(poly,dpoly,degree,form[0],dform[0]);for(int z=0;z<e[1];++z)MulLinear3(poly,dpoly,degree,form[1],dform[1]);val[l]=der[l]=0;for(int total=0;total<=degree;++total)for(int ax=0;ax<=total;++ax)for(int ay=0;ay<=total-ax;++ay){const int id=MomentIndex3(ax,ay,total-ax-ay);val[l]+=poly[id]*mv[id];der[l]+=dpoly[id]*mv[id]+poly[id]*dm[id];}}
    };
    auto directed=[&](int kt,int ht,int ks,int hs,int img){
        const auto&gt=kt==0?m_hoCellCharges[ht]:m_hoFaceCharges[ht];
        const auto&gs=ks==0?m_hoCellCharges[hs]:m_hoFaceCharges[hs];
        std::vector<double>b((size_t)gt.size()*gs.size()),v,d;
        for(size_t q=0;q<m_qp[gt[0]].size();++q){
            const auto&pp=m_qp[gt[0]][q];
            double p[3]={pp[0],pp[1],pp[2]},dp[3],rate;
            target_kinematics(kt,ht,p,dp,rate);
            // QuadDotRefl(tgt,src,img) evaluates the ordinary source
            // potential at T^-1(x).  The image transform is a FIXED linear
            // isometry, so the exact directional derivative is obtained by
            // mapping both the moving outer point and its velocity through
            // it; the measure rate is unchanged by the isometry.
            if(img!=0){
                double pm[3],dpm[3];
                ImageEvalPoint(img,p,pm);ImageEvalPoint(img,dp,dpm);
                for(int axis=0;axis<3;++axis){p[axis]=pm[axis];dp[axis]=dpm[axis];}
            }
            source_inner(ks,hs,p,dp,v,d);
            for(size_t i=0;i<gt.size();++i){
                const double w=m_qw[gt[i]][q];
                for(size_t j=0;j<gs.size();++j)
                    b[i*gs.size()+j]+=w*(d[j]+rate*v[j]);
            }
        }
        for(double&x:b)x*=RAD_INV_FOUR_PI;
        return b;
    };
    const int nh=nc+nf;
    auto pair_block=[&](int ga,int gb){
        const int ka=ga<nc?0:1,ha=ga<nc?ga:ga-nc;
        const int kb=gb<nc?0:1,hb=gb<nc?gb:gb-nc;
        const auto&A=ka==0?m_hoCellCharges[ha]:m_hoFaceCharges[ha];
        const auto&B=kb==0?m_hoCellCharges[hb]:m_hoFaceCharges[hb];
        std::vector<double> block;
        if(ga==gb&&!m_polyCombo)
            block=ka==0
                ?TetVolumeSelfBlockDirectionalDerivative(
                    ha,std::vector<double>(cell_velocity.begin()+(size_t)ha*12,
                                           cell_velocity.begin()+(size_t)(ha+1)*12))
                :TetFaceSelfBlockDirectionalDerivative(
                    ha,std::vector<double>(face_velocity.begin()+(size_t)ha*9,
                                           face_velocity.begin()+(size_t)(ha+1)*9));
        else{
            auto ab=directed(ka,ha,kb,hb,0);
            auto ba=directed(kb,hb,ka,ha,0);
            block.resize((size_t)A.size()*B.size());
            for(size_t i=0;i<A.size();++i)for(size_t j=0;j<B.size();++j)
                block[i*B.size()+j]=.5*(ab[i*B.size()+j]+ba[j*A.size()+i]);
        }
        // Differentiate the same symmetrized IMA fold used by the parent
        // ChargeGram.  Image self-hosts are ordinary reflected pairs, never
        // singular self-panel terms.
        for(size_t image=0;image<m_image_masks.size();++image){
            const auto ab=directed(ka,ha,kb,hb,(int)image+1);
            const auto ba=directed(kb,hb,ka,ha,(int)image+1);
            for(size_t i=0;i<A.size();++i)for(size_t j=0;j<B.size();++j)
                block[i*B.size()+j]+=m_image_signs[image]*.5*
                    (ab[i*B.size()+j]+ba[j*A.size()+i]);
        }
        return block;
    };
    if(selected_host_a>=0){if(selected_host_a>selected_host_b)std::swap(selected_host_a,selected_host_b);if(selected_host_b>=nh)throw std::out_of_range("TET derivative host out of range");return pair_block(selected_host_a,selected_host_b);}
    out.assign((size_t)m_n*m_n,0.0);
    for(int ga=0;ga<nh;++ga){const int ka=ga<nc?0:1,ha=ga<nc?ga:ga-nc;const auto&A=ka==0?m_hoCellCharges[ha]:m_hoFaceCharges[ha];for(int gb=ga;gb<nh;++gb){const int kb=gb<nc?0:1,hb=gb<nc?gb:gb-nc;const auto&B=kb==0?m_hoCellCharges[hb]:m_hoFaceCharges[hb];const auto block=pair_block(ga,gb);for(size_t i=0;i<A.size();++i)for(size_t j=0;j<B.size();++j){const double x=block[i*B.size()+j];out[(size_t)A[i]*m_n+B[j]]=x;out[(size_t)B[j]*m_n+A[i]]=x;}}}
    return out;
}

std::vector<double> RadHACApKChargeGram::TetChargeMapRowDirectionalRates(
    const std::vector<double>& cell_velocity,const std::vector<double>& face_velocity) const
{
    if(!m_highorder||m_curved||m_hexmode||m_wedgemode||m_polyCombo)
        throw std::logic_error("TET charge-map derivative requires a flat polynomial TET charge Gram");
    const int nc=(int)m_hoCellCharges.size(),nf=(int)m_hoFaceCharges.size();
    if(cell_velocity.size()!=(size_t)nc*12)throw std::invalid_argument("cell_vertex_velocity must have shape (ncell,4,3)");
    if(face_velocity.size()!=(size_t)nf*9)throw std::invalid_argument("face_vertex_velocity must have shape (nface,3,3)");
    std::vector<double> rates(m_n,0.0);
    for(int h=0;h<nc;++h){const double*s=&m_cellV[(size_t)h*12],*ds=&cell_velocity[(size_t)h*12];double E[9],dE[9];for(int k=0;k<3;++k)for(int j=0;j<3;++j){E[3*k+j]=s[3*(j+1)+k]-s[k];dE[3*k+j]=ds[3*(j+1)+k]-ds[k];}const double r=-Det3Rate(E,dE);for(int q:m_hoCellCharges[h])rates[q]=r;}
    for(int h=0;h<nf;++h){const double*s=&m_faceV[(size_t)h*9],*ds=&face_velocity[(size_t)h*9];double a[2][3],da[2][3];for(int c=0;c<2;++c)for(int k=0;k<3;++k){a[c][k]=s[3*(c+1)+k]-s[k];da[c][k]=ds[3*(c+1)+k]-ds[k];}double cr[3]={a[0][1]*a[1][2]-a[0][2]*a[1][1],a[0][2]*a[1][0]-a[0][0]*a[1][2],a[0][0]*a[1][1]-a[0][1]*a[1][0]},dcr[3]={da[0][1]*a[1][2]+a[0][1]*da[1][2]-da[0][2]*a[1][1]-a[0][2]*da[1][1],da[0][2]*a[1][0]+a[0][2]*da[1][0]-da[0][0]*a[1][2]-a[0][0]*da[1][2],da[0][0]*a[1][1]+a[0][0]*da[1][1]-da[0][1]*a[1][0]-a[0][1]*da[1][0]};double nn=0,cn=0;for(int k=0;k<3;++k){nn+=cr[k]*cr[k];cn+=cr[k]*dcr[k];}const double r=-cn/nn;for(int q:m_hoFaceCharges[h])rates[q]=r;}
    return rates;
}

// Duffy singular-quadrature inner potential INT_host(src) m_src(y)/|p-y| dy for the order>=3 / curved path
// (where the analytic moment kernels run out: a tet volume charge of degree>=2 needs TetMoment2; a surface
// charge of degree>=3 needs degree-3 moments).  6-pt Gauss-Legendre on signed radial sub-tets (cell) / signed
// sub-triangles (face) from x0 = closest point of the host to p; the Duffy Jacobian (u^2 for tet, u for tri)
// regularizes the 1/r, and the SIGNED sub-simplices telescope to INT_host for any x0 (inside / on / outside).
// Validated to ~1e-4 vs independent Duffy/brute-force prototypes.
double RadHACApKChargeGram::PhiAtHO_Duffy(int src, const double p[3]) const
{
    static const double GL[6] = {0.03376524289842399, 0.16939530676686777, 0.38069040695840156,
                                 0.61930959304159840, 0.83060469323313230, 0.96623475710157600};
    static const double GW[6] = {0.08566224618958520, 0.18038078652406930, 0.23395696728634550,
                                 0.23395696728634550, 0.18038078652406930, 0.08566224618958520};
    const int host = m_host[src];
    double acc = 0.0;
    if (m_kind[src] == 0) {                                   // ---- tet cell: 4 signed radial sub-tets ----
        double V[4][3]; const double* s = &m_cellV[(size_t)host*12];
        for (int i=0;i<4;++i) for (int k=0;k<3;++k) V[i][k] = s[3*i+k];
        // The signed sub-tets give the SIGNED-volume integral; the physical charge integral uses the ABSOLUTE
        // volume, so multiply by sign(host signed vol) (= -1 for a negatively-oriented mesh tet).
        double E0[3], E1[3], E2[3];
        for (int k=0;k<3;++k){ E0[k]=V[1][k]-V[0][k]; E1[k]=V[2][k]-V[0][k]; E2[k]=V[3][k]-V[0][k]; }
        const double hv = E0[0]*(E1[1]*E2[2]-E1[2]*E2[1]) - E0[1]*(E1[0]*E2[2]-E1[2]*E2[0])
                        + E0[2]*(E1[0]*E2[1]-E1[1]*E2[0]);
        const double sgn_host = (hv >= 0.0) ? 1.0 : -1.0;
        double x0[3]; rad_hdiv::ClosestPointTet(V, p, x0);
        static const int FC[4][3] = {{1,2,3},{0,3,2},{0,1,3},{2,1,0}};
        for (int f = 0; f < 4; ++f) {
            const double* b1 = V[FC[f][0]]; const double* b2 = V[FC[f][1]]; const double* b3 = V[FC[f][2]];
            double d1[3],d2[3],d3[3],e21[3],e32[3];
            for (int k=0;k<3;++k){ d1[k]=b1[k]-x0[k]; d2[k]=b2[k]-x0[k]; d3[k]=b3[k]-x0[k];
                                   e21[k]=b2[k]-b1[k]; e32[k]=b3[k]-b2[k]; }
            const double cr[3] = {d2[1]*d3[2]-d2[2]*d3[1], d2[2]*d3[0]-d2[0]*d3[2], d2[0]*d3[1]-d2[1]*d3[0]};
            const double D = d1[0]*cr[0]+d1[1]*cr[1]+d1[2]*cr[2];   // signed 6*vol(x0,b1,b2,b3)
            if (std::fabs(D) < 1e-300) continue;
            for (int a=0;a<6;++a){ const double u=GL[a];
                for (int b=0;b<6;++b){ const double v=GL[b];
                    for (int c=0;c<6;++c){ const double w=GL[c];
                        double y[3]; for (int k=0;k<3;++k) y[k]=x0[k]+u*(d1[k]+v*(e21[k]+w*e32[k]));
                        const double dx=p[0]-y[0], dy=p[1]-y[1], dz=p[2]-y[2];
                        const double r=std::sqrt(dx*dx+dy*dy+dz*dz);
                        if (r<1e-300) continue;
                        acc += GW[a]*GW[b]*GW[c]*(u*u*v*D)*EvalMono(src,y)/r;
                    }}}
        }
        return acc * sgn_host;
    }
    // ---- tri face: 3 signed sub-triangles from x0 = projection of p onto the face plane ----
    double V[3][3]; const double* s = &m_faceV[(size_t)host*9];
    for (int i=0;i<3;++i) for (int k=0;k<3;++k) V[i][k] = s[3*i+k];
    double e1[3], e2[3];
    for (int k=0;k<3;++k){ e1[k]=V[1][k]-V[0][k]; e2[k]=V[2][k]-V[0][k]; }
    double nrm[3] = {e1[1]*e2[2]-e1[2]*e2[1], e1[2]*e2[0]-e1[0]*e2[2], e1[0]*e2[1]-e1[1]*e2[0]};
    const double nl = std::sqrt(nrm[0]*nrm[0]+nrm[1]*nrm[1]+nrm[2]*nrm[2]);
    if (nl < 1e-300) return 0.0;
    nrm[0]/=nl; nrm[1]/=nl; nrm[2]/=nl;
    const double hh = (p[0]-V[0][0])*nrm[0]+(p[1]-V[0][1])*nrm[1]+(p[2]-V[0][2])*nrm[2];
    const double x0[3] = {p[0]-hh*nrm[0], p[1]-hh*nrm[1], p[2]-hh*nrm[2]};
    for (int kf = 0; kf < 3; ++kf) {
        const double* A = V[kf]; const double* B = V[(kf+1)%3];
        double ea[3], eb[3];
        for (int k=0;k<3;++k){ ea[k]=A[k]-x0[k]; eb[k]=B[k]-x0[k]; }
        const double cx[3] = {ea[1]*eb[2]-ea[2]*eb[1], ea[2]*eb[0]-ea[0]*eb[2], ea[0]*eb[1]-ea[1]*eb[0]};
        const double sgn2 = cx[0]*nrm[0]+cx[1]*nrm[1]+cx[2]*nrm[2];   // signed 2*area(x0,A,B)
        for (int a=0;a<6;++a){ const double u=GL[a];
            for (int b=0;b<6;++b){ const double v=GL[b];
                double y[3]; for (int k=0;k<3;++k) y[k]=x0[k]+u*ea[k]+u*v*(eb[k]-ea[k]);
                const double dx=p[0]-y[0], dy=p[1]-y[1], dz=p[2]-y[2];
                const double r=std::sqrt(dx*dx+dy*dy+dz*dz);
                if (r<1e-300) continue;
                acc += GW[a]*GW[b]*(u*sgn2)*EvalMono(src,y)/r;
            }}
    }
    return acc;
}

// CURVED (isoparametric P2) inner potential: the curved Duffy at the host's P2 nodes (the monomial is in the
// REFERENCE frame, so CurvedTet/TriPotential -- which fold xi^e and the curved measure -- is the full potential
// of source charge src at p).  No analytic moments exist on a curved element; this is the SOLE curved path.
double RadHACApKChargeGram::PhiAtHO_Curved(int src, const double p[3]) const
{
    const int host = m_host[src];
    const int* e = &m_expo[(size_t)3*src];
    const int nq = (int)m_gl.size();
    if (m_kind[src] == 0) {
        const double (*nd)[3] = (const double(*)[3])&m_cellNodes[(size_t)host*30];
        return rad_hdiv::CurvedTetPotential(
            nd, e[0], e[1], e[2], p, m_gl.data(), m_gw.data(), nq,
            !m_curvedReferenceDensity);
    }
    const double (*nd)[3] = (const double(*)[3])&m_faceNodes[(size_t)host*18];
    return rad_hdiv::CurvedTriPotential(
        nd, e[0], e[1], p, m_gl.data(), m_gw.data(), nq,
        !m_curvedReferenceDensity);
}

// Dispatch the high-order inner potential: CURVED -> the curved Duffy; else FLAT -> the EXACT analytic moment
// kernels where they suffice (charge degree<=2: a tet up to degree 1, a face up to degree 2), else the flat
// Duffy singular quadrature (order>=3).
double RadHACApKChargeGram::PhiInner(int src, const double p[3]) const
{
    if (m_curved) return PhiAtHO_Curved(src, p);
    const int* e = &m_expo[(size_t)3*src];
    const int deg = e[0] + e[1] + e[2];
    const bool analytic_ok = (m_kind[src] == 0) ? (deg <= 1) : (deg <= 2);
    return analytic_ok ? PhiAtHO_Analytic(src, p) : PhiAtHO_Duffy(src, p);
}

// polynomial-charge inner potential INT_host(src) m_src(y)/|p-y| dy by singularity SUBTRACTION reusing the
// exact constant-charge PhiTet/TriPotential: = m_src(p) Phi_host(p) + sum_q W_q (m_src(y_q) - m_src(p))/|p-y_q|.
double RadHACApKChargeGram::PhiAtHO(int src, const double p[3]) const
{
    const double msrc_p = EvalMono(src, p);
    const int host = m_host[src];
    double base;
    if (m_kind[src] == 0) {
        double V[4][3]; const double* s = &m_cellV[(size_t)host*12];
        for (int i = 0; i < 4; ++i) for (int k = 0; k < 3; ++k) V[i][k] = s[3*i+k];
        base = msrc_p * rad_hdiv::PhiTet(V, p);
    } else {
        double V[3][3]; const double* s = &m_faceV[(size_t)host*9];
        for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) V[i][k] = s[3*i+k];
        base = msrc_p * rad_hdiv::TriPotential(V, p);
    }
    const std::vector<rad_hdiv::Vec3>& Y = m_inP[src];
    const std::vector<double>&         W = m_inW[src];
    double rem = 0.0;
    for (size_t q = 0; q < Y.size(); ++q) {
        const double dx = p[0]-Y[q][0], dy = p[1]-Y[q][1], dz = p[2]-Y[q][2];
        const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 1e-300) continue;                          // p coincides with a node (self block) -> Phi term covers it
        rem += W[q] * (m_srcval[src][q] - msrc_p) / r;     // m_srcval[src][q] == EvalMono(src, Y[q]) (precomputed)
    }
    return base + rem;
}

double RadHACApKChargeGram::PhiAt(int src, const double p[3]) const
{
    if (m_polytope) {
        if (m_curved_face) {     // FULLY CURVED polytope: curved sub-tet (cell) / sub-tri (face) potentials
            const int nq = (int)m_gl.size();
            double tot = 0.0;
            if (src < m_n_el) {  // CELL volume charge: sum CurvedTetPotential over the cell's curved sub-tets
                for (const auto& nd10 : m_srcCurvedTets[src]) {
                    double nd[10][3];
                    for (int i = 0; i < 10; ++i) for (int k = 0; k < 3; ++k) nd[i][k] = nd10[i][k];
                    tot += rad_hdiv::CurvedTetPotential(nd, 0, 0, 0, p, m_gl.data(), m_gw.data(), nq);
                }
            } else {             // FACE surface charge: sum CurvedTriPotential over the face's curved sub-tris
                for (const auto& nd6 : m_srcCurvedTris[src]) {
                    double nd[6][3];
                    for (int i = 0; i < 6; ++i) for (int k = 0; k < 3; ++k) nd[i][k] = nd6[i][k];
                    tot += rad_hdiv::CurvedTriPotential(nd, 0, 0, p, m_gl.data(), m_gw.data(), nq);
                }
            }
            return tot;          // constant charge -> monomial exponent 0
        }
        const std::vector<std::array<rad_hdiv::Vec3, 3>>& tris = m_srcTris[src];
        if (src < m_n_el) {
            // CELL: divergence-theorem polytope potential = (1/2) sum_tri d_tri * TriPotential(tri,p),
            // d_tri = (T0 - p).n_out, n_out the OUTWARD unit normal (flipped via the cell centroid --
            // matches radia.vim._core._cell_hull_tris / _polytope_potential).
            const double cx = m_cent[3*src], cy = m_cent[3*src+1], cz = m_cent[3*src+2];
            double tot = 0.0;
            for (const auto& T : tris) {
                double e1[3], e2[3];
                for (int k = 0; k < 3; ++k) { e1[k] = T[1][k]-T[0][k]; e2[k] = T[2][k]-T[0][k]; }
                double n[3] = {e1[1]*e2[2]-e1[2]*e2[1], e1[2]*e2[0]-e1[0]*e2[2], e1[0]*e2[1]-e1[1]*e2[0]};
                const double nl = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
                if (nl < 1e-300) continue;
                n[0]/=nl; n[1]/=nl; n[2]/=nl;
                const double tcx = (T[0][0]+T[1][0]+T[2][0])/3.0, tcy = (T[0][1]+T[1][1]+T[2][1])/3.0,
                             tcz = (T[0][2]+T[1][2]+T[2][2])/3.0;
                if (n[0]*(tcx-cx) + n[1]*(tcy-cy) + n[2]*(tcz-cz) < 0.0) { n[0]=-n[0]; n[1]=-n[1]; n[2]=-n[2]; }
                const double d = (T[0][0]-p[0])*n[0] + (T[0][1]-p[1])*n[1] + (T[0][2]-p[2])*n[2];
                double V[3][3] = {{T[0][0],T[0][1],T[0][2]}, {T[1][0],T[1][1],T[1][2]}, {T[2][0],T[2][1],T[2][2]}};
                tot += d * rad_hdiv::TriPotential(V, p);
            }
            return 0.5 * tot;
        }
        double tot = 0.0;        // FACE: sum of flat sub-triangle Wilton potentials
        for (const auto& T : tris) {
            double V[3][3] = {{T[0][0],T[0][1],T[0][2]}, {T[1][0],T[1][1],T[1][2]}, {T[2][0],T[2][1],T[2][2]}};
            tot += rad_hdiv::TriPotential(V, p);
        }
        return tot;
    }
    if (src < m_n_el) {
        double V[4][3];
        const double* s = &m_cellV[(size_t)src * 12];
        for (int i = 0; i < 4; ++i) for (int k = 0; k < 3; ++k) V[i][k] = s[3*i+k];
        return rad_hdiv::PhiTet(V, p);
    }
    double V[3][3];
    const double* s = &m_faceV[(size_t)(src - m_n_el) * 9];
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) V[i][k] = s[3*i+k];
    return rad_hdiv::TriPotential(V, p);
}

double RadHACApKChargeGram::QuadDot(int tgt, int src) const
{
    const std::vector<rad_hdiv::Vec3>& P = m_qp[tgt];
    const std::vector<double>&         W = m_qw[tgt];
    if (m_highorder && m_nmono[tgt] > 1) {
        // CO-LOCATED MEMO (non-HACApK-path build speedup): the host carries m_nmono>1 monomials sharing these
        // outer points, so PhiAtHO(src, P[k]) is identical across them -> cache the Qout-length potential
        // vector per (kind,host,src) and reuse it BIT-EXACT for the host's other monomials.  Keyed on
        // (host,src) DIRECTLY (NOT cleared on host change): GetInteractionMatrixElement's symmetrization
        // 0.5*(QuadDot(a,b)+QuadDot(b,a)) alternates the tgt-host (host_a then host_b), so a clear-on-host
        // cache would thrash to zero hits.  With the (host,src) key both directions reuse: the row direction
        // across the co-located rows of a leaf, the col direction across the consecutive co-located cols.
        // Cap-based eviction bounds the per-thread working set; cleared on a new build (owner id).
        // CORRECTNESS DEPENDS on m_build_id being GLOBALLY unique across ALL constructors (NextChargeGramBuildId):
        // this thread_local memo outlives a single build (the TaskManager threadpool persists across builds in
        // one TM region), so a colliding id from a sibling constructor would skip the clear and reuse another
        // build's stale PhiInner values.  See NextChargeGramBuildId (2026-06-29 shared-TM corruption fix).
        static thread_local long long cache_owner = -1;
        static thread_local std::unordered_map<long long, std::vector<double>> cache;
        if (cache_owner != m_build_id) { cache.clear(); cache_owner = m_build_id; }
        const long long key = ((long long)(m_host[tgt]*2 + m_kind[tgt]) << 32) | (long long)(unsigned)src;
        auto it = cache.find(key);
        const std::vector<double>* phi;
        if (it != cache.end()) {
            phi = &it->second;
        } else {
            if (cache.size() > 32768u) cache.clear();   // memory cap (~16 MB/thread at Qout~64); rare flush
            std::vector<double> v(P.size());
            for (size_t k = 0; k < P.size(); ++k) {
                const double p[3] = {P[k][0], P[k][1], P[k][2]};
                v[k] = PhiInner(src, p);
            }
            phi = &cache.emplace(key, std::move(v)).first->second;
        }
        double s = 0.0;
        for (size_t k = 0; k < P.size(); ++k) s += W[k] * (*phi)[k];
        return s * RAD_INV_FOUR_PI;
    }
    double s = 0.0;
    for (size_t k = 0; k < P.size(); ++k) {
        const double p[3] = {P[k][0], P[k][1], P[k][2]};
        s += W[k] * (m_highorder ? PhiInner(src, p) : PhiAt(src, p));
    }
    return s * RAD_INV_FOUR_PI;
}

double RadHACApKChargeGram::QuadDotRefl(int tgt, int src, int img) const
{
    // IMA image entry G_img(tgt,src) = (1/4pi) INT_tgt Phi_{T(src)} = (1/4pi) INT_tgt Phi_src(T^-1(x))
    // (isometry |x - T(y)| = |T^-1(x) - y|), so we map tgt's outer points by the INVERSE image transform and
    // evaluate the UNmapped source potential there.  Always the full analytic PhiAt (the image of a
    // charge straddling the mirror is singular at the plane -> needs the exact through-singularity potential).
    const std::vector<rad_hdiv::Vec3>& P = m_qp[tgt];
    const std::vector<double>&         W = m_qw[tgt];
    double s = 0.0;
    for (size_t k = 0; k < P.size(); ++k) {
        const double p0[3] = {P[k][0], P[k][1], P[k][2]};
        double p[3]; ImageEvalPoint(img, p0, p);
        // m_highorder: the monomial-charge inner potential (host-agnostic potential-at-p, PhiAtHO_*);
        // m_analytic/polytope: the constant-charge PhiAt.  Image charge -> inverse-mapped eval point.
        s += W[k] * (m_highorder ? PhiInner(src, p) : PhiAt(src, p));
    }
    return s * RAD_INV_FOUR_PI;
}

double RadHACApKChargeGram::QuadDotFar(int tgt, int src) const
{
    // cheap FAR evaluation (near/far adaptive quadrature): plain LOW-quad double Gauss of
    //   (1/4pi) INT_tgt INT_src m_t(x) m_s(y) / |x-y|.
    // Only called for WELL-SEPARATED pairs, where 1/|x-y| is SMOOTH -> low-order Gauss is accurate and the
    // singularity-subtraction (PhiTet + inner subtraction) of the NEAR QuadDot is unnecessary.  m_t(x) is
    // folded into m_qw_lo (outer); m_s(y) is evaluated on the fly for the plain inner sum.
    const std::vector<rad_hdiv::Vec3>& Px = m_qp_lo[tgt];
    const std::vector<double>&         Wx = m_qw_lo[tgt];
    const std::vector<rad_hdiv::Vec3>& Py = m_inP_lo[src];
    const std::vector<double>&         Wy = m_inW_lo[src];
    double s = 0.0;
    for (size_t i = 0; i < Px.size(); ++i) {
        const double x0 = Px[i][0], x1 = Px[i][1], x2 = Px[i][2];
        double inner = 0.0;
        for (size_t j = 0; j < Py.size(); ++j) {
            const double dx = x0 - Py[j][0], dy = x1 - Py[j][1], dz = x2 - Py[j][2];
            const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (r < 1e-300) continue;                        // far pairs are well-separated; defensive only
            inner += Wy[j] * m_srcval_lo[src][j] / r;        // m_srcval_lo[src][j] == EvalMono(src, Py[j]) (precomputed)
        }
        s += Wx[i] * inner;                                  // Wx folds m_t(x)
    }
    return s * RAD_INV_FOUR_PI;
}

double RadHACApKChargeGram::QuadDotFarLow(int a, int b) const
{
    // Precision-preserving FAR evaluation (analytic mode, far_quad>0): a plain LOW-order double-quadrature of
    //   (1/4pi) INT_a INT_b 1/|x-y|
    // over the degree-2 rules (4-pt tet / 3-pt tri).  Far pairs are well-separated, so 1/|x-y| is smooth and
    // the degree-2 rule (exact through quadrupole moments) reproduces the all-analytic entry to O((size/r)^4)
    // -- vs the monopole's O((size/r)^2).  ~16 cheap evals/pair vs the NEAR QuadDot's ~1e3 transcendentals
    // (PhiTet/TriPotential).  Symmetric in (a,b) (1/r symmetric), so no 0.5*(.+.) needed.
    const std::vector<rad_hdiv::Vec3>& Pa = m_qpf[a];
    const std::vector<double>&         Wa = m_qwf[a];
    const std::vector<rad_hdiv::Vec3>& Pb = m_qpf[b];
    const std::vector<double>&         Wb = m_qwf[b];
    double s = 0.0;
    for (size_t i = 0; i < Pa.size(); ++i) {
        const double x0 = Pa[i][0], x1 = Pa[i][1], x2 = Pa[i][2];
        double inner = 0.0;
        for (size_t j = 0; j < Pb.size(); ++j) {
            const double dx = x0 - Pb[j][0], dy = x1 - Pb[j][1], dz = x2 - Pb[j][2];
            inner += Wb[j] / std::sqrt(dx*dx + dy*dy + dz*dz);
        }
        s += Wa[i] * inner;
    }
    return s * RAD_INV_FOUR_PI;
}

// ===================================================================== HEX BDM1/BDM2 mode
// Direct Q2 isoparametric geometry + the numpy-validated eig(M^-1 N)<=1 quadrature scheme (see the header
// ctor doc).  Reference tables: the unit hex [0,1]^3 with its Kuhn 6-sub-tet split (shared main diagonal
// 0-6) and the unit quad [0,1]^2 with its 2-sub-tri split -- the SAME tables as the Python contract.
// A future PYRAMID adds one more (ref sub-tet table + nodal map +
// monomial set) row here -- no structural change.
static const double HEXREF_V[8][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
static const int    HEXREF_TETS[6][4] = {{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}};
static const double QUADREF_V[4][2] = {{0,0},{1,0},{1,1},{0,1}};
static const int    QUADREF_TRIS[2][3] = {{0,1,2},{0,2,3}};

static inline void HexLag3(double t, double v[3], double d[3])
{
    v[0] = 2.0*(t-0.5)*(t-1.0); v[1] = 4.0*t*(1.0-t); v[2] = 2.0*t*(t-0.5);
    d[0] = 4.0*t-3.0;           d[1] = 4.0-8.0*t;     d[2] = 4.0*t-1.0;
}

void RadHACApKChargeGram::HexQ2Map(const double* nd27, const double xi[3], double X[3], double J[3][3])
{
    double vx[3], dx[3], vy[3], dy[3], vz[3], dz[3];
    HexLag3(xi[0], vx, dx); HexLag3(xi[1], vy, dy); HexLag3(xi[2], vz, dz);
    for (int k = 0; k < 3; ++k) { X[k] = 0.0; J[k][0] = J[k][1] = J[k][2] = 0.0; }
    for (int iz = 0; iz < 3; ++iz)
        for (int iy = 0; iy < 3; ++iy)
            for (int ix = 0; ix < 3; ++ix) {
                const double* nd = &nd27[3*(ix + 3*iy + 9*iz)];
                const double s  = vx[ix]*vy[iy]*vz[iz];
                const double gx = dx[ix]*vy[iy]*vz[iz];
                const double gy = vx[ix]*dy[iy]*vz[iz];
                const double gz = vx[ix]*vy[iy]*dz[iz];
                for (int k = 0; k < 3; ++k) {
                    X[k] += s*nd[k]; J[k][0] += gx*nd[k]; J[k][1] += gy*nd[k]; J[k][2] += gz*nd[k];
                }
            }
}

void RadHACApKChargeGram::QuadQ2Map(const double* nd9, const double uv[2], double X[3], double T[3][2])
{
    double vu[3], du[3], vv[3], dv[3];
    HexLag3(uv[0], vu, du); HexLag3(uv[1], vv, dv);
    for (int k = 0; k < 3; ++k) { X[k] = 0.0; T[k][0] = T[k][1] = 0.0; }
    for (int iv = 0; iv < 3; ++iv)
        for (int iu = 0; iu < 3; ++iu) {
            const double* nd = &nd9[3*(iu + 3*iv)];
            const double s = vu[iu]*vv[iv];
            for (int k = 0; k < 3; ++k) {
                X[k] += s*nd[k]; T[k][0] += du[iu]*vv[iv]*nd[k]; T[k][1] += vu[iu]*dv[iv]*nd[k];
            }
        }
}

// Values-only Q2 maps (no Jacobian): the Piola radial inner needs only X (REF measure -- no |det J|).
void RadHACApKChargeGram::HexQ2MapX(const double* nd27, const double xi[3], double X[3])
{
    double vx[3], dx[3], vy[3], dy[3], vz[3], dz[3];
    HexLag3(xi[0], vx, dx); HexLag3(xi[1], vy, dy); HexLag3(xi[2], vz, dz);
    X[0] = X[1] = X[2] = 0.0;
    for (int iz = 0; iz < 3; ++iz)
        for (int iy = 0; iy < 3; ++iy) {
            const double syz = vy[iy]*vz[iz];
            for (int ix = 0; ix < 3; ++ix) {
                const double* nd = &nd27[3*(ix + 3*iy + 9*iz)];
                const double s = vx[ix]*syz;
                X[0] += s*nd[0]; X[1] += s*nd[1]; X[2] += s*nd[2];
            }
        }
}

void RadHACApKChargeGram::QuadQ2MapX(const double* nd9, const double uv[2], double X[3])
{
    double vu[3], du[3], vv[3], dv[3];
    HexLag3(uv[0], vu, du); HexLag3(uv[1], vv, dv);
    X[0] = X[1] = X[2] = 0.0;
    for (int iv = 0; iv < 3; ++iv)
        for (int iu = 0; iu < 3; ++iu) {
            const double* nd = &nd9[3*(iu + 3*iv)];
            const double s = vu[iu]*vv[iv];
            X[0] += s*nd[0]; X[1] += s*nd[1]; X[2] += s*nd[2];
        }
}

// ============================================ WEDGE (PRISM) BDM1 geometry (2026-07-04) ===================
// Prism ref domain: (u,v) in the triangle {u>=0, v>=0, u+v<=1}, w in [0,1].  Corners 0-2 = bottom tri at
// z=0, 3-5 = top tri at z=1.  The 3-sub-tet split tiles the prism (each 6*vol_ref = 1; total 3*(1/6) =
// 1/2 = the prism ref volume).  A tri FACE ref = the same triangle (Tri6 corner order (1,0),(0,1),(0,0),
// matching TriSurfMap / D2_TRIREF); its single sub-tri IS the whole ref tri (2*area_ref = 1).
static const double WEDGEREF_V[6][3]  = {{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,0,1},{0,1,1}};
static const int    WEDGEREF_TETS[3][4] = {{0,1,2,5},{0,1,5,4},{0,4,5,3}};
static const double WTRIREF_V[3][2]   = {{1,0},{0,1},{0,0}};   // tri-face ref (Tri6 corner order)

// tri-P2 shape functions (barycentric quadratic; l0=u, l1=v, l2=1-u-v -- IDENTICAL to Tri6Map so the
// 18-node prism lattice node n = t + 6*iz uses the same (u,v) node layout as _TRI6_LAT on the Python side).
static inline void TriP2Shape(double u, double v, double s[6])
{
    const double l0 = u, l1 = v, l2 = 1.0 - u - v;
    s[0] = l0*(2*l0 - 1); s[1] = l1*(2*l1 - 1); s[2] = l2*(2*l2 - 1);
    s[3] = 4*l0*l1;       s[4] = 4*l1*l2;       s[5] = 4*l2*l0;
}

// 18-node prism map (tri-P2 (x) z-P2): node n = t + 6*iz, t = tri node 0..5, iz = z level 0..2.  Values-
// only (the Piola charge model never needs |det J|).
void RadHACApKChargeGram::WedgeQ2MapX(const double* nd18, const double xi[3], double X[3])
{
    double st[6]; TriP2Shape(xi[0], xi[1], st);
    double vz[3], dz[3]; HexLag3(xi[2], vz, dz);
    X[0] = X[1] = X[2] = 0.0;
    for (int iz = 0; iz < 3; ++iz)
        for (int t = 0; t < 6; ++t) {
            const double s = st[t]*vz[iz];
            const double* nd = &nd18[3*(t + 6*iz)];
            X[0] += s*nd[0]; X[1] += s*nd[1]; X[2] += s*nd[2];
        }
}

// 6-node quadratic surface-triangle map (a boundary tri cap lives in 3D): nd18 = 6 nodes x 3D.
void RadHACApKChargeGram::TriSurfMap(const double* nd18, const double uv[2], double X[3])
{
    double st[6]; TriP2Shape(uv[0], uv[1], st);
    X[0] = X[1] = X[2] = 0.0;
    for (int t = 0; t < 6; ++t) { X[0] += st[t]*nd18[3*t]; X[1] += st[t]*nd18[3*t+1]; X[2] += st[t]*nd18[3*t+2]; }
}

// 6*vol of the ref sub-tet s (WEDGEREF); 2*area of the whole tri-face ref (both = 1 for these splits, but
// computed generically for a future pyramid row).
static inline double WedgeSubSixVref(int s)
{
    const int* tv = WEDGEREF_TETS[s];
    double e[3][3];
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k) e[i][k] = WEDGEREF_V[tv[i+1]][k] - WEDGEREF_V[tv[0]][k];
    return std::fabs(e[0][0]*(e[1][1]*e[2][2]-e[1][2]*e[2][1]) - e[0][1]*(e[1][0]*e[2][2]-e[1][2]*e[2][0])
                     + e[0][2]*(e[1][0]*e[2][1]-e[1][1]*e[2][0]));
}
static inline double WTriSubTwoAref()
{
    const double a1u = WTRIREF_V[1][0]-WTRIREF_V[0][0], a1v = WTRIREF_V[1][1]-WTRIREF_V[0][1];
    const double a2u = WTRIREF_V[2][0]-WTRIREF_V[0][0], a2v = WTRIREF_V[2][1]-WTRIREF_V[0][1];
    return std::fabs(a1u*a2v - a1v*a2u);
}

// Forward decl: the wedge ctor (below) uses this file-static face-ref helper whose definition lives with
// the wedge compute block further down.
static void WFaceSubTriRef(int face_type, int s, double V[3][2]);

// Radial-cone face table of the ref sub-tet (vertex i's opposite face, oriented so the signed 6-vol D of
// (x0, b1, b2, b3) sums the tet exactly from any interior anchor) -- shared by the SELF radial and the
// static-SITE table generator.
static const int HEXTET_FC[4][3] = {{1,2,3},{0,3,2},{0,1,3},{2,1,0}};

static inline double HexDet3(const double J[3][3])
{
    return J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1]) - J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])
         + J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
}

static inline double HexSurfJ(const double T[3][2])
{
    const double cx = T[1][0]*T[2][1] - T[2][0]*T[1][1];
    const double cy = T[2][0]*T[0][1] - T[0][0]*T[2][1];
    const double cz = T[0][0]*T[1][1] - T[1][0]*T[0][1];
    return std::sqrt(cx*cx + cy*cy + cz*cz);
}

// Bary-rule scale factors: a rule whose weights sum to the UNIT-simplex measure (1/6 tet, 1/2 tri)
// integrates over a ref sub-simplex of measure V_sub as  sum (W_q * 6 V_sub) g  /  sum (W_q * 2 A_sub) g.
// For the Kuhn 6-tet / 2-tri splits both factors are exactly 1; computed generically (future pyramid rows).
static inline double HexSubSixVref(int s)
{
    const int* tv = HEXREF_TETS[s];
    double e[3][3];
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k) e[i][k] = HEXREF_V[tv[i+1]][k] - HEXREF_V[tv[0]][k];
    return std::fabs(e[0][0]*(e[1][1]*e[2][2]-e[1][2]*e[2][1]) - e[0][1]*(e[1][0]*e[2][2]-e[1][2]*e[2][0])
                     + e[0][2]*(e[1][0]*e[2][1]-e[1][1]*e[2][0]));
}

static inline double QuadSubTwoAref(int s)
{
    const int* tv = QUADREF_TRIS[s];
    const double a1u = QUADREF_V[tv[1]][0]-QUADREF_V[tv[0]][0], a1v = QUADREF_V[tv[1]][1]-QUADREF_V[tv[0]][1];
    const double a2u = QUADREF_V[tv[2]][0]-QUADREF_V[tv[0]][0], a2v = QUADREF_V[tv[2]][1]-QUADREF_V[tv[0]][1];
    return std::fabs(a1u*a2v - a1v*a2u);
}

double RadHACApKChargeGram::HexMonoEval(int charge, const double xi[3]) const
{
    const int* e = &m_expo[(size_t)3*charge];
    double v = 1.0;
    for (int d = 0; d < 3; ++d)
        for (int k = 0; k < e[d]; ++k) v *= xi[d];
    return v;
}

static constexpr int HEX_AFFINE_POLY_N = 84;
static constexpr int QUAD_AFFINE_POLY_N = 35;

static inline int HexPolyIdx(int ax, int ay, int az)
{
    const int deg = ax + ay + az;
    int idx = 0;
    for (int k = 0; k < deg; ++k) idx += (k + 1) * (k + 2) / 2;
    for (int a = 0; a < ax; ++a) idx += deg - a + 1;
    idx += ay;
    return idx;
}

static void HexPolyMulLinear(double* poly, int& deg, const double lin[4], int ncoeff)
{
    double tmp[HEX_AFFINE_POLY_N] = {};
    for (int total = 0; total <= deg; ++total) {
        for (int ax = 0; ax <= total; ++ax) {
            for (int ay = 0; ay <= total - ax; ++ay) {
                const int az = total - ax - ay;
                const double c = poly[HexPolyIdx(ax, ay, az)];
                if (c == 0.0) continue;
                tmp[HexPolyIdx(ax,     ay,     az    )] += c * lin[0];
                tmp[HexPolyIdx(ax + 1, ay,     az    )] += c * lin[1];
                tmp[HexPolyIdx(ax,     ay + 1, az    )] += c * lin[2];
                tmp[HexPolyIdx(ax,     ay,     az + 1)] += c * lin[3];
            }
        }
    }
    ++deg;
    for (int i = 0; i < ncoeff; ++i) poly[i] = tmp[i];
}

static void HexPolyMulLinearDirectional(double* poly,double* dpoly,int& deg,
                                        const double lin[4],const double dlin[4],int ncoeff)
{
    double tmp[HEX_AFFINE_POLY_N]={},dtmp[HEX_AFFINE_POLY_N]={};
    for(int total=0;total<=deg;++total)for(int ax=0;ax<=total;++ax)for(int ay=0;ay<=total-ax;++ay){
        const int az=total-ax-ay,idx=HexPolyIdx(ax,ay,az);const double c=poly[idx],dc=dpoly[idx];
        const int ids[4]={HexPolyIdx(ax,ay,az),HexPolyIdx(ax+1,ay,az),HexPolyIdx(ax,ay+1,az),HexPolyIdx(ax,ay,az+1)};
        for(int q=0;q<4;++q){tmp[ids[q]]+=c*lin[q];dtmp[ids[q]]+=dc*lin[q]+c*dlin[q];}
    }
    ++deg;for(int i=0;i<ncoeff;++i){poly[i]=tmp[i];dpoly[i]=dtmp[i];}
}

static bool HexAffineInverseForms(const double* nd27, double lin[3][4], double& inv_abs_det)
{
    const double* o = &nd27[0];
    const double* px = &nd27[3*2];
    const double* py = &nd27[3*6];
    const double* pz = &nd27[3*18];
    double A[3][3];
    for (int k = 0; k < 3; ++k) {
        A[k][0] = px[k] - o[k];
        A[k][1] = py[k] - o[k];
        A[k][2] = pz[k] - o[k];
    }
    const double det = HexDet3(A);
    double scale = 0.0;
    for (int j = 0; j < 3; ++j) {
        double n2 = 0.0;
        for (int k = 0; k < 3; ++k) n2 += A[k][j] * A[k][j];
        scale = std::max(scale, std::sqrt(n2));
    }
    if (std::fabs(det) < 1e-300 || scale < 1e-300) return false;
    const double tol = 1e-10 * scale + 1e-12;
    for (int iz = 0; iz < 3; ++iz)
        for (int iy = 0; iy < 3; ++iy)
            for (int ix = 0; ix < 3; ++ix) {
                const double xi = 0.5 * ix, eta = 0.5 * iy, zeta = 0.5 * iz;
                const double* p = &nd27[3*(ix + 3*iy + 9*iz)];
                double pred[3];
                for (int k = 0; k < 3; ++k) pred[k] = o[k] + xi*A[k][0] + eta*A[k][1] + zeta*A[k][2];
                const double err = std::sqrt((p[0]-pred[0])*(p[0]-pred[0]) +
                                             (p[1]-pred[1])*(p[1]-pred[1]) +
                                             (p[2]-pred[2])*(p[2]-pred[2]));
                if (err > tol) return false;
            }
    double B[3][3];
    const double id = 1.0 / det;
    B[0][0] =  (A[1][1]*A[2][2] - A[1][2]*A[2][1]) * id;
    B[0][1] = -(A[0][1]*A[2][2] - A[0][2]*A[2][1]) * id;
    B[0][2] =  (A[0][1]*A[1][2] - A[0][2]*A[1][1]) * id;
    B[1][0] = -(A[1][0]*A[2][2] - A[1][2]*A[2][0]) * id;
    B[1][1] =  (A[0][0]*A[2][2] - A[0][2]*A[2][0]) * id;
    B[1][2] = -(A[0][0]*A[1][2] - A[0][2]*A[1][0]) * id;
    B[2][0] =  (A[1][0]*A[2][1] - A[1][1]*A[2][0]) * id;
    B[2][1] = -(A[0][0]*A[2][1] - A[0][1]*A[2][0]) * id;
    B[2][2] =  (A[0][0]*A[1][1] - A[0][1]*A[1][0]) * id;
    for (int q = 0; q < 3; ++q) {
        lin[q][0] = -(B[q][0]*o[0] + B[q][1]*o[1] + B[q][2]*o[2]);
        lin[q][1] = B[q][0];
        lin[q][2] = B[q][1];
        lin[q][3] = B[q][2];
    }
    inv_abs_det = 1.0 / std::fabs(det);
    return true;
}

static bool QuadAffineInverseForms(const double* nd9, double lin[2][4], double& inv_surface_jac)
{
    const double* o = &nd9[0];
    const double* pu = &nd9[6];
    const double* pv = &nd9[18];
    double a[3], b[3];
    for (int k = 0; k < 3; ++k) { a[k] = pu[k] - o[k]; b[k] = pv[k] - o[k]; }
    const double aa = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    const double ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    const double bb = b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
    const double det = aa*bb - ab*ab;
    if (det <= 1e-300) return false;
    const double scale = std::max(std::sqrt(aa), std::sqrt(bb));
    const double tol = 1e-10*scale + 1e-12;
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i) {
            const double u = 0.5*i, v = 0.5*j;
            const double* q = &nd9[3*(i + 3*j)];
            const double dx = q[0] - (o[0] + u*a[0] + v*b[0]);
            const double dy = q[1] - (o[1] + u*a[1] + v*b[1]);
            const double dz = q[2] - (o[2] + u*a[2] + v*b[2]);
            if (dx*dx + dy*dy + dz*dz > tol*tol) return false;
        }
    const double i00 = bb/det, i01 = -ab/det, i11 = aa/det;
    for (int k = 0; k < 3; ++k) {
        lin[0][k+1] = i00*a[k] + i01*b[k];
        lin[1][k+1] = i01*a[k] + i11*b[k];
    }
    lin[0][0] = -(lin[0][1]*o[0] + lin[0][2]*o[1] + lin[0][3]*o[2]);
    lin[1][0] = -(lin[1][1]*o[0] + lin[1][2]*o[1] + lin[1][3]*o[2]);
    inv_surface_jac = 1.0/std::sqrt(det);
    return true;
}


static bool QuadAffineInverseFormsDirectional(const double* nd9,const double* dnd9,
    double lin[2][4],double dlin[2][4],double& inv_surface_jac,double& dinv_surface_jac)
{
    if(!QuadAffineInverseForms(nd9,lin,inv_surface_jac))return false;
    double a[3],b[3],da[3],db[3];for(int k=0;k<3;++k){a[k]=nd9[6+k]-nd9[k];b[k]=nd9[18+k]-nd9[k];da[k]=dnd9[6+k]-dnd9[k];db[k]=dnd9[18+k]-dnd9[k];}
    double aa=0,ab=0,bb=0,daa=0,dab=0,dbb=0;for(int k=0;k<3;++k){aa+=a[k]*a[k];ab+=a[k]*b[k];bb+=b[k]*b[k];daa+=2*a[k]*da[k];dab+=da[k]*b[k]+a[k]*db[k];dbb+=2*b[k]*db[k];}
    const double det=aa*bb-ab*ab,ddet=daa*bb+aa*dbb-2*ab*dab;
    dinv_surface_jac=-0.5*inv_surface_jac*ddet/det;
    const double Gi[2][2]={{bb/det,-ab/det},{-ab/det,aa/det}},dG[2][2]={{daa,dab},{dab,dbb}};
    double L[2][3],dC[2][3];for(int k=0;k<3;++k){L[0][k]=lin[0][k+1];L[1][k]=lin[1][k+1];dC[0][k]=da[k];dC[1][k]=db[k];}
    for(int q=0;q<2;++q)for(int k=0;k<3;++k){double v=0;for(int i=0;i<2;++i){v+=Gi[q][i]*dC[i][k];for(int j=0;j<2;++j)v-=Gi[q][i]*dG[i][j]*L[j][k];}dlin[q][k+1]=v;}
    for(int q=0;q<2;++q){dlin[q][0]=0;for(int k=0;k<3;++k)dlin[q][0]-=dlin[q][k+1]*nd9[k]+lin[q][k+1]*dnd9[k];}
    return true;
}

static bool HexAffineBasisChecked(const double* nd27, double origin[3], double A[3][3], double& det)
{
    const double* o = &nd27[0];
    const double* px = &nd27[3*2];
    const double* py = &nd27[3*6];
    const double* pz = &nd27[3*18];
    for (int k = 0; k < 3; ++k) {
        origin[k] = o[k];
        A[k][0] = px[k] - o[k];
        A[k][1] = py[k] - o[k];
        A[k][2] = pz[k] - o[k];
    }
    det = HexDet3(A);
    double scale = 0.0;
    for (int j = 0; j < 3; ++j) {
        double n2 = 0.0;
        for (int k = 0; k < 3; ++k) n2 += A[k][j] * A[k][j];
        scale = std::max(scale, std::sqrt(n2));
    }
    if (std::fabs(det) < 1e-300 || scale < 1e-300) return false;
    const double tol = 1e-10 * scale + 1e-12;
    for (int iz = 0; iz < 3; ++iz)
        for (int iy = 0; iy < 3; ++iy)
            for (int ix = 0; ix < 3; ++ix) {
                const double xi = 0.5 * ix, eta = 0.5 * iy, zeta = 0.5 * iz;
                const double* p = &nd27[3*(ix + 3*iy + 9*iz)];
                double pred[3];
                for (int k = 0; k < 3; ++k) pred[k] = origin[k] + xi*A[k][0] + eta*A[k][1] + zeta*A[k][2];
                const double err = std::sqrt((p[0]-pred[0])*(p[0]-pred[0]) +
                                             (p[1]-pred[1])*(p[1]-pred[1]) +
                                             (p[2]-pred[2])*(p[2]-pred[2]));
                if (err > tol) return false;
            }
    return true;
}

static bool HexInv3(const double A[3][3], double B[3][3])
{
    const double det = HexDet3(A);
    if (std::fabs(det) < 1e-300) return false;
    const double id = 1.0 / det;
    B[0][0] =  (A[1][1]*A[2][2] - A[1][2]*A[2][1]) * id;
    B[0][1] = -(A[0][1]*A[2][2] - A[0][2]*A[2][1]) * id;
    B[0][2] =  (A[0][1]*A[1][2] - A[0][2]*A[1][1]) * id;
    B[1][0] = -(A[1][0]*A[2][2] - A[1][2]*A[2][0]) * id;
    B[1][1] =  (A[0][0]*A[2][2] - A[0][2]*A[2][0]) * id;
    B[1][2] = -(A[0][0]*A[1][2] - A[0][2]*A[1][0]) * id;
    B[2][0] =  (A[1][0]*A[2][1] - A[1][1]*A[2][0]) * id;
    B[2][1] = -(A[0][0]*A[2][1] - A[0][1]*A[2][0]) * id;
    B[2][2] =  (A[0][0]*A[1][1] - A[0][1]*A[1][0]) * id;
    return true;
}

// Exact affine-cell inner is retained for self/near blocks, where the singularity and spectrum require it.
// Smooth far blocks use the complete-host tensor product below; applying the degree-six recurrence to every
// far outer point was the dominant BDM2 HEX H-matrix build cost.
static constexpr bool HEX_USE_AFFINE_EXACT_CELL_INNER = true;
static constexpr double HEX_AFFINE_EXACT_NEAR_FACTOR = 1.0;

// Build a Duffy-graded barycentric rule on a (dim+1)-vertex ref sub-simplex from the 1D rule (gl,gw),
// graded at LOCAL vertex `corner` (swap-permuted to Duffy vertex 0, matching the validated
// _ref_duffy_corner / _graded_outer_bary).  Appends (bary[nv], w_ref) pairs; w_ref sums to the ref
// simplex measure (1/6 tet, 1/2 tri).
static void HexDuffyBary(int dim, int corner, const std::vector<double>& gl, const std::vector<double>& gw,
                         std::vector<double>& bary_out, std::vector<double>& w_out)
{
    const int n = (int)gl.size();
    const int nv = dim + 1;
    // The Duffy APEX -- the node-accumulating, jac->0, singularity-RESOLVING vertex -- is barycentric
    // index 1 (bary = (L0, a, b, c) with a = u -> 1 at the apex; L0 = (1-u)(1-v)(1-w) -> 0 there), NOT
    // index 0.  Swap 1 <-> corner so the rule actually grades at the requested vertex.  (The old swap
    // 0 <-> corner graded at vertex 1 for every corner except corner==1 (then vertex 0) -- a latent
    // off-by-one inherited from the numpy prototype's mislabeled bary_std comment, MASKED by the
    // 1000-pt glin^3 inner density; exposed by the linearized-subtraction identity, whose a2-term
    // requires the remainder rule to resolve 1/|A dxi| at the graded corner -> face-self eig ~1e12.)
    int perm[4] = {0, 1, 2, 3};
    perm[1] = corner; perm[corner] = 1;              // swap 1 <-> corner (apex -> corner)
    if (dim == 3) {
        for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) for (int k = 0; k < n; ++k) {
            const double u = gl[i], v = gl[j], w = gl[k];
            const double a = u, b = v*(1.0-u), c = w*(1.0-u)*(1.0-v);
            const double jac = (1.0-u)*(1.0-u)*(1.0-v);
            double bstd[4] = {1.0-a-b-c, a, b, c};
            double b4[4];
            for (int t = 0; t < 4; ++t) b4[perm[t]] = bstd[t];
            for (int t = 0; t < 4; ++t) bary_out.push_back(b4[t]);
            w_out.push_back(gw[i]*gw[j]*gw[k]*jac);
        }
    } else {
        for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
            const double u = gl[i], v = gl[j];
            const double a = u, b = v*(1.0-u);
            const double jac = 1.0-u;
            double bstd[3] = {1.0-a-b, a, b};
            double b3[3];
            for (int t = 0; t < 3; ++t) b3[perm[t]] = bstd[t];
            for (int t = 0; t < 3; ++t) bary_out.push_back(b3[t]);
            w_out.push_back(gw[i]*gw[j]*jac);
        }
    }
    (void)nv;
}

RadHACApKChargeGram::RadHACApKChargeGram(
    std::vector<double> hex_cell_nodes, std::vector<double> quad_face_nodes,
    int n_el, int n_bf,
    std::vector<int> charge_host, std::vector<int> charge_kind, std::vector<int> charge_expo,
    std::vector<double> sym_tet_pts, std::vector<double> sym_tet_w,
    std::vector<double> sym_tri_pts, std::vector<double> sym_tri_w,
    std::vector<double> gl_out, std::vector<double> gw_out,
    std::vector<double> gl_in, std::vector<double> gw_in,
    std::vector<double> far_tet_pts, std::vector<double> far_tet_w,
    std::vector<double> far_tri_pts, std::vector<double> far_tri_w,
    double near_grade, double far_inner_factor,
    std::vector<int> image_masks, std::vector<double> image_signs)
    : m_n_el(n_el), m_hexmode(true), m_hex_n_bf(n_bf),
      m_hexNodes(std::move(hex_cell_nodes)), m_quadNodes(std::move(quad_face_nodes)),
      m_symTetP(std::move(sym_tet_pts)), m_symTetW(std::move(sym_tet_w)),
      m_symTriP(std::move(sym_tri_pts)), m_symTriW(std::move(sym_tri_w)),
      m_glOut(std::move(gl_out)), m_gwOut(std::move(gw_out)),
      m_glIn(std::move(gl_in)), m_gwIn(std::move(gw_in)),
      m_farTetP(std::move(far_tet_pts)), m_farTetW(std::move(far_tet_w)),
      m_farTriP(std::move(far_tri_pts)), m_farTriW(std::move(far_tri_w)),
      m_near_grade(near_grade), m_far_inner_factor(far_inner_factor),
      m_image_masks(std::move(image_masks)), m_image_signs(std::move(image_signs)),
      m_host(std::move(charge_host)), m_kind(std::move(charge_kind)), m_expo(std::move(charge_expo))
{
    ValidateImageVectors(m_image_masks, m_image_signs);
    m_n = (int)m_host.size();
    if (m_kind.size() != m_host.size() || m_expo.size() != 3*m_host.size())
        throw std::invalid_argument("HEX ChargeGram charge metadata sizes are inconsistent");
    for (int exponent : m_expo)
        if (exponent < 0 || exponent > 2)
            throw std::invalid_argument("HEX ChargeGram supports reference exponents in {0,1,2}");
    m_build_id = NextChargeGramBuildId();
    m_hexCacheStatsEnabled = HexCacheStatsEnabledByEnv();
    // ---- per-host sub-simplex physical geometry (corners, centroid, size) via the Q2 maps ----
    m_cellSubV.assign((size_t)n_el*6*4*3, 0.0); m_cellSubC.assign((size_t)n_el*6*3, 0.0);
    m_cellSubS.assign((size_t)n_el*6, 0.0);
    for (int c = 0; c < n_el; ++c) {
        const double* nd = &m_hexNodes[(size_t)c*81];
        for (int s = 0; s < 6; ++s) {
            double cen[3] = {0, 0, 0};
            for (int i = 0; i < 4; ++i) {
                const double* rv = HEXREF_V[HEXREF_TETS[s][i]];
                double X[3], J[3][3];
                const double xi[3] = {rv[0], rv[1], rv[2]};
                HexQ2Map(nd, xi, X, J);
                double* out = &m_cellSubV[(((size_t)c*6 + s)*4 + i)*3];
                out[0] = X[0]; out[1] = X[1]; out[2] = X[2];
                cen[0] += 0.25*X[0]; cen[1] += 0.25*X[1]; cen[2] += 0.25*X[2];
            }
            double* pc = &m_cellSubC[((size_t)c*6 + s)*3];
            pc[0] = cen[0]; pc[1] = cen[1]; pc[2] = cen[2];
            double sz = 0.0;
            for (int i = 0; i < 4; ++i) {
                const double* v = &m_cellSubV[(((size_t)c*6 + s)*4 + i)*3];
                const double dx = v[0]-cen[0], dy = v[1]-cen[1], dz = v[2]-cen[2];
                sz = std::max(sz, std::sqrt(dx*dx + dy*dy + dz*dz));
            }
            m_cellSubS[(size_t)c*6 + s] = sz;
        }
    }
    m_faceSubV.assign((size_t)n_bf*2*3*3, 0.0); m_faceSubC.assign((size_t)n_bf*2*3, 0.0);
    m_faceSubS.assign((size_t)n_bf*2, 0.0);
    for (int f = 0; f < n_bf; ++f) {
        const double* nd = &m_quadNodes[(size_t)f*27];
        for (int s = 0; s < 2; ++s) {
            double cen[3] = {0, 0, 0};
            for (int i = 0; i < 3; ++i) {
                const double* rv = QUADREF_V[QUADREF_TRIS[s][i]];
                double X[3], T[3][2];
                const double uv[2] = {rv[0], rv[1]};
                QuadQ2Map(nd, uv, X, T);
                double* out = &m_faceSubV[(((size_t)f*2 + s)*3 + i)*3];
                out[0] = X[0]; out[1] = X[1]; out[2] = X[2];
                cen[0] += X[0]/3.0; cen[1] += X[1]/3.0; cen[2] += X[2]/3.0;
            }
            double* pc = &m_faceSubC[((size_t)f*2 + s)*3];
            pc[0] = cen[0]; pc[1] = cen[1]; pc[2] = cen[2];
            double sz = 0.0;
            for (int i = 0; i < 3; ++i) {
                const double* v = &m_faceSubV[(((size_t)f*2 + s)*3 + i)*3];
                const double dx = v[0]-cen[0], dy = v[1]-cen[1], dz = v[2]-cen[2];
                sz = std::max(sz, std::sqrt(dx*dx + dy*dy + dz*dz));
            }
            m_faceSubS[(size_t)f*2 + s] = sz;
        }
    }
    // ---- per-charge host centroid/size (cluster-tree points + the near_hosts test) ----
    m_cent.assign((size_t)m_n*3, 0.0); m_size.assign((size_t)m_n, 0.0);
    for (int a = 0; a < m_n; ++a) {
        const int h = m_host[a];
        double cen[3] = {0, 0, 0};
        int ncorner;
        double corners[8][3];
        if (m_kind[a] == 0) {                       // hex corners of the 27-lattice: ix,iy,iz in {0,2}
            ncorner = 8;
            static const int cidx[8] = {0, 2, 6, 8, 18, 20, 24, 26};
            for (int i = 0; i < 8; ++i) {
                const double* nd = &m_hexNodes[(size_t)h*81 + 3*cidx[i]];
                for (int k = 0; k < 3; ++k) corners[i][k] = nd[k];
            }
        } else {                                    // quad corners of the 9-lattice
            ncorner = 4;
            static const int cidx[4] = {0, 2, 6, 8};
            for (int i = 0; i < 4; ++i) {
                const double* nd = &m_quadNodes[(size_t)h*27 + 3*cidx[i]];
                for (int k = 0; k < 3; ++k) corners[i][k] = nd[k];
            }
        }
        for (int i = 0; i < ncorner; ++i)
            for (int k = 0; k < 3; ++k) cen[k] += corners[i][k] / ncorner;
        double sz = 0.0;
        for (int i = 0; i < ncorner; ++i) {
            const double dx = corners[i][0]-cen[0], dy = corners[i][1]-cen[1], dz = corners[i][2]-cen[2];
            sz = std::max(sz, std::sqrt(dx*dx + dy*dy + dz*dz));
        }
        m_cent[3*a] = cen[0]; m_cent[3*a+1] = cen[1]; m_cent[3*a+2] = cen[2];
        m_size[a] = sz;
    }
    // ---- (kind,host)->local reverse maps for the block memo (co-located charges grouped per host) ----
    m_hexLocalOf.assign((size_t)m_n, 0);
    m_cellCharges.assign((size_t)n_el, {}); m_faceCharges.assign((size_t)n_bf, {});
    for (int a = 0; a < m_n; ++a) {
        std::vector<int>& grp = (m_kind[a] == 0) ? m_cellCharges[m_host[a]] : m_faceCharges[m_host[a]];
        m_hexLocalOf[a] = (int)grp.size();
        grp.push_back(a);
    }
    // ---- affine HEX source exact path, matching the flat-TET BDM2 design: reference Q2 charges become
    // physical-coordinate polynomials and the source inner is evaluated by analytic moments. ----
    m_hexAffineOrder = 1;
    for (int exponent : m_expo) m_hexAffineOrder = std::max(m_hexAffineOrder, exponent);
    const int affineAxisCount = m_hexAffineOrder + 1;
    m_hexAffineMonoCount = affineAxisCount*affineAxisCount*affineAxisCount;
    m_hexAffinePolyCount = (3*m_hexAffineOrder + 1)*(3*m_hexAffineOrder + 2)
                           *(3*m_hexAffineOrder + 3)/6;
    m_quadAffineMonoCount = affineAxisCount*affineAxisCount;
    m_quadAffinePolyCount = (2*m_hexAffineOrder + 1)*(2*m_hexAffineOrder + 2)
                            *(2*m_hexAffineOrder + 3)/6;
    m_hexAffineCell.assign((size_t)n_el, 0);
    m_hexAffineCoeff.assign((size_t)n_el * m_hexAffineMonoCount * m_hexAffinePolyCount, 0.0);
    if (HEX_USE_AFFINE_EXACT_CELL_INNER) {
        for (int c = 0; c < n_el; ++c) {
            double lin[3][4], inv_abs_det = 0.0;
            if (!HexAffineInverseForms(&m_hexNodes[(size_t)c*81], lin, inv_abs_det)) continue;
            m_hexAffineCell[c] = 1;
            for (int ez = 0; ez <= m_hexAffineOrder; ++ez)
              for (int ey = 0; ey <= m_hexAffineOrder; ++ey)
               for (int ex = 0; ex <= m_hexAffineOrder; ++ex) {
                const int mono = ex + affineAxisCount*ey + affineAxisCount*affineAxisCount*ez;
                double poly[HEX_AFFINE_POLY_N] = {};
                int deg = 0;
                poly[0] = 1.0;
                for (int repeat = 0; repeat < ex; ++repeat)
                    HexPolyMulLinear(poly, deg, lin[0], m_hexAffinePolyCount);
                for (int repeat = 0; repeat < ey; ++repeat)
                    HexPolyMulLinear(poly, deg, lin[1], m_hexAffinePolyCount);
                for (int repeat = 0; repeat < ez; ++repeat)
                    HexPolyMulLinear(poly, deg, lin[2], m_hexAffinePolyCount);
                double* dst = &m_hexAffineCoeff[
                    ((size_t)c*m_hexAffineMonoCount + mono)*m_hexAffinePolyCount];
                for (int i = 0; i < m_hexAffinePolyCount; ++i) dst[i] = inv_abs_det * poly[i];
            }
        }
    }
    m_quadAffineFace.assign((size_t)n_bf, 0);
    m_quadAffineCoeff.assign((size_t)n_bf * m_quadAffineMonoCount * m_quadAffinePolyCount, 0.0);
    for (int f = 0; f < n_bf; ++f) {
        double lin[2][4], inv_surface_jac = 0.0;
        if (!QuadAffineInverseForms(&m_quadNodes[(size_t)f*27], lin, inv_surface_jac)) continue;
        m_quadAffineFace[f] = 1;
        for (int ev = 0; ev <= m_hexAffineOrder; ++ev)
          for (int eu = 0; eu <= m_hexAffineOrder; ++eu) {
            const int mono = eu + affineAxisCount*ev;
            double poly[HEX_AFFINE_POLY_N] = {};
            int deg = 0;
            poly[0] = 1.0;
            for (int repeat = 0; repeat < eu; ++repeat)
                HexPolyMulLinear(poly, deg, lin[0], m_quadAffinePolyCount);
            for (int repeat = 0; repeat < ev; ++repeat)
                HexPolyMulLinear(poly, deg, lin[1], m_quadAffinePolyCount);
            double* dst = &m_quadAffineCoeff[
                ((size_t)f*m_quadAffineMonoCount + mono)*m_quadAffinePolyCount];
            for (int i = 0; i < m_quadAffinePolyCount; ++i) dst[i] = inv_surface_jac * poly[i];
        }
    }
    // Structured affine meshes (e.g. the cube scaling benchmark) repeat the same cell-cell block for
    // every equal lattice offset.  Detect that case once; GetHexBlock then caches cell-cell blocks by
    // integer offset instead of absolute host ids, preserving the exact quadrature while avoiding
    // translation-duplicate entry fills.
    m_hexUniformAffineCells = false;
    m_hexCellLattice.assign((size_t)n_el*3, 0);
    if (n_el > 0) {
        double o0[3], A0[3][3], det0 = 0.0, B0[3][3];
        bool ok = HexAffineBasisChecked(&m_hexNodes[0], o0, A0, det0) && HexInv3(A0, B0);
        double c0[3] = {o0[0] + 0.5*(A0[0][0] + A0[0][1] + A0[0][2]),
                        o0[1] + 0.5*(A0[1][0] + A0[1][1] + A0[1][2]),
                        o0[2] + 0.5*(A0[2][0] + A0[2][1] + A0[2][2])};
        double scale = 0.0;
        for (int j = 0; j < 3; ++j) {
            double n2 = 0.0;
            for (int k = 0; k < 3; ++k) n2 += A0[k][j] * A0[k][j];
            scale = std::max(scale, std::sqrt(n2));
        }
        const double tolA = 1e-10 * scale + 1e-12;
        for (int c = 0; ok && c < n_el; ++c) {
            double oc[3], Ac[3][3], detc = 0.0;
            ok = HexAffineBasisChecked(&m_hexNodes[(size_t)c*81], oc, Ac, detc);
            for (int i = 0; ok && i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    if (std::fabs(Ac[i][j] - A0[i][j]) > tolA) ok = false;
            if (!ok) break;
            double cc[3] = {oc[0] + 0.5*(Ac[0][0] + Ac[0][1] + Ac[0][2]),
                            oc[1] + 0.5*(Ac[1][0] + Ac[1][1] + Ac[1][2]),
                            oc[2] + 0.5*(Ac[2][0] + Ac[2][1] + Ac[2][2])};
            double q[3] = {0.0, 0.0, 0.0};
            for (int i = 0; i < 3; ++i)
                for (int k = 0; k < 3; ++k) q[i] += B0[i][k] * (cc[k] - c0[k]);
            for (int i = 0; i < 3; ++i) {
                const int qi = (int)std::llround(q[i]);
                if (std::fabs(q[i] - qi) > 1e-7) ok = false;
                m_hexCellLattice[(size_t)3*c + i] = qi;
            }
        }
        m_hexUniformAffineCells = ok;
    }
    // Extend the structured affine-cache idea from cell-cell blocks to every translated host pair
    // (cell-face and face-face too).  A host is cacheable when its Q2 lattice nodes, after subtracting the
    // host center, match one of a small set of templates and its center lies on the half-cell lattice of the
    // first affine cell.  The template id captures face orientation / parameterization, so blocks are reused
    // only for true translated copies.
    m_hexUniformTransHosts = false;
    m_hexHostTemplate.assign((size_t)n_el + (size_t)m_hex_n_bf, -1);
    m_hexHostLattice2.assign(((size_t)n_el + (size_t)m_hex_n_bf)*3, 0);
    if (m_hexUniformAffineCells && n_el > 0) {
        double o0[3], A0[3][3], det0 = 0.0, B0[3][3];
        bool ok = HexAffineBasisChecked(&m_hexNodes[0], o0, A0, det0) && HexInv3(A0, B0);
        double c0[3] = {o0[0] + 0.5*(A0[0][0] + A0[0][1] + A0[0][2]),
                        o0[1] + 0.5*(A0[1][0] + A0[1][1] + A0[1][2]),
                        o0[2] + 0.5*(A0[2][0] + A0[2][1] + A0[2][2])};
        double scale = 0.0;
        for (int j = 0; j < 3; ++j) {
            double n2 = 0.0;
            for (int k = 0; k < 3; ++k) n2 += A0[k][j] * A0[k][j];
            scale = std::max(scale, std::sqrt(n2));
        }
        const double tol = 1e-10 * scale + 1e-12;
        std::vector<std::vector<double>> cell_templates, face_templates;
        auto classify_host = [&](int kind, int h, const double* nd, int nnode, size_t host_index) {
            double ctr[3] = {0.0, 0.0, 0.0};
            for (int i = 0; i < nnode; ++i)
                for (int k = 0; k < 3; ++k) ctr[k] += nd[3*i + k];
            for (int k = 0; k < 3; ++k) ctr[k] /= (double)nnode;
            double q[3] = {0.0, 0.0, 0.0};
            for (int i = 0; i < 3; ++i)
                for (int k = 0; k < 3; ++k) q[i] += B0[i][k] * (ctr[k] - c0[k]);
            for (int i = 0; i < 3; ++i) {
                const int qi2 = (int)std::llround(2.0 * q[i]);
                if (std::fabs(2.0 * q[i] - qi2) > 2e-7) return false;
                m_hexHostLattice2[3*host_index + i] = qi2;
            }
            const std::vector<int>& grp = (kind == 0) ? m_cellCharges[h] : m_faceCharges[h];
            std::vector<double> sig;
            sig.reserve((size_t)nnode*3 + 1 + (size_t)grp.size()*3);
            for (int i = 0; i < nnode; ++i)
                for (int k = 0; k < 3; ++k) sig.push_back(nd[3*i + k] - ctr[k]);
            sig.push_back((double)grp.size());
            for (int g : grp)
                for (int k = 0; k < 3; ++k) sig.push_back((double)m_expo[(size_t)3*g + k]);
            auto& templ = (kind == 0) ? cell_templates : face_templates;
            int tid = -1;
            for (int t = 0; t < (int)templ.size(); ++t) {
                if (templ[t].size() != sig.size()) continue;
                bool same = true;
                for (size_t i = 0; i < sig.size(); ++i)
                    if (std::fabs(templ[t][i] - sig[i]) > tol) { same = false; break; }
                if (same) { tid = t; break; }
            }
            if (tid < 0) { tid = (int)templ.size(); templ.push_back(std::move(sig)); }
            m_hexHostTemplate[host_index] = tid;
            return true;
        };
        for (int c = 0; ok && c < n_el; ++c)
            ok = classify_host(0, c, &m_hexNodes[(size_t)c*81], 27, (size_t)c);
        for (int f = 0; ok && f < m_hex_n_bf; ++f)
            ok = classify_host(1, f, &m_quadNodes[(size_t)f*27], 9, (size_t)n_el + (size_t)f);
        m_hexUniformTransHosts = ok;
    }
    BuildHexSiteTables();   // static-site radial tables (non-self near inner) + mapped site positions
}

// WEDGE (PRISM) BDM1 ctor -- mirror of the hex ctor with 3-sub-tet prism cells + mixed tri/quad faces (see
// the header doc).  Reuses the hex-mode quadrature-table + block-serving members; fills only the wedge
// geometry.  Initializer list is in member DECLARATION order (m_n_el, the shared quad tables, the wedge
// nodes, then m_host/m_kind/m_expo) to avoid -Wreorder.
RadHACApKChargeGram::RadHACApKChargeGram(
    std::vector<double> wedge_cell_nodes, std::vector<double> face_nodes, std::vector<int> face_type,
    int n_el, int n_bf,
    std::vector<int> charge_host, std::vector<int> charge_kind, std::vector<int> charge_expo,
    std::vector<double> sym_tet_pts, std::vector<double> sym_tet_w,
    std::vector<double> sym_tri_pts, std::vector<double> sym_tri_w,
    std::vector<double> field_tri_pts, std::vector<double> field_tri_w,
    std::vector<double> gl_out, std::vector<double> gw_out,
    std::vector<double> gl_in, std::vector<double> gw_in,
    std::vector<double> far_tet_pts, std::vector<double> far_tet_w,
    std::vector<double> far_tri_pts, std::vector<double> far_tri_w,
    double near_grade, double far_inner_factor,
    std::vector<int> image_masks, std::vector<double> image_signs)
    : m_n_el(n_el),
      m_symTetP(std::move(sym_tet_pts)), m_symTetW(std::move(sym_tet_w)),
      m_symTriP(std::move(sym_tri_pts)), m_symTriW(std::move(sym_tri_w)),
      m_glOut(std::move(gl_out)), m_gwOut(std::move(gw_out)),
      m_glIn(std::move(gl_in)), m_gwIn(std::move(gw_in)),
      m_farTetP(std::move(far_tet_pts)), m_farTetW(std::move(far_tet_w)),
      m_farTriP(std::move(far_tri_pts)), m_farTriW(std::move(far_tri_w)),
      m_near_grade(near_grade), m_far_inner_factor(far_inner_factor),
      m_image_masks(std::move(image_masks)), m_image_signs(std::move(image_signs)),
      m_wedgemode(true), m_wedge_n_bf(n_bf),
      m_wCellNodes(std::move(wedge_cell_nodes)), m_wFaceNodes(std::move(face_nodes)),
      m_wFaceType(std::move(face_type)),
      m_wFieldTriP(std::move(field_tri_pts)), m_wFieldTriW(std::move(field_tri_w)),
      m_host(std::move(charge_host)), m_kind(std::move(charge_kind)), m_expo(std::move(charge_expo))
{
    ValidateImageVectors(m_image_masks, m_image_signs);
    m_n = (int)m_host.size();
    if (m_wFieldTriW.empty() || m_wFieldTriP.size() != 2*m_wFieldTriW.size())
        throw std::invalid_argument("WEDGE ChargeGram field triangle rule sizes are inconsistent");
    if (m_kind.size() != m_host.size() || m_expo.size() != 3*m_host.size())
        throw std::invalid_argument("WEDGE ChargeGram charge metadata sizes are inconsistent");
    for (int exponent : m_expo)
        if (exponent < 0 || exponent > 2)
            throw std::invalid_argument("WEDGE ChargeGram supports reference exponents in {0,1,2}");
    m_build_id = NextChargeGramBuildId();
    m_hexCacheStatsEnabled = HexCacheStatsEnabledByEnv();
    // ---- cell sub-tet physical geometry (3 sub-tets per prism) ----
    m_wCellSubV.assign((size_t)n_el*3*4*3, 0.0); m_wCellSubC.assign((size_t)n_el*3*3, 0.0);
    m_wCellSubS.assign((size_t)n_el*3, 0.0);
    for (int c = 0; c < n_el; ++c) {
        const double* nd = &m_wCellNodes[(size_t)c*54];
        for (int s = 0; s < 3; ++s) {
            double cen[3] = {0, 0, 0};
            for (int i = 0; i < 4; ++i) {
                const double* rv = WEDGEREF_V[WEDGEREF_TETS[s][i]];
                double X[3]; const double xi[3] = {rv[0], rv[1], rv[2]};
                WedgeQ2MapX(nd, xi, X);
                double* out = &m_wCellSubV[(((size_t)c*3 + s)*4 + i)*3];
                out[0] = X[0]; out[1] = X[1]; out[2] = X[2];
                cen[0] += 0.25*X[0]; cen[1] += 0.25*X[1]; cen[2] += 0.25*X[2];
            }
            double* pc = &m_wCellSubC[((size_t)c*3 + s)*3];
            pc[0] = cen[0]; pc[1] = cen[1]; pc[2] = cen[2];
            double sz = 0.0;
            for (int i = 0; i < 4; ++i) {
                const double* v = &m_wCellSubV[(((size_t)c*3 + s)*4 + i)*3];
                const double dx = v[0]-cen[0], dy = v[1]-cen[1], dz = v[2]-cen[2];
                sz = std::max(sz, std::sqrt(dx*dx + dy*dy + dz*dz));
            }
            m_wCellSubS[(size_t)c*3 + s] = sz;
        }
    }
    // ---- face sub-tri physical geometry (tri: 1 sub-tri; quad: 2 sub-tris) ----
    m_wFaceSubV.assign((size_t)n_bf*2*3*3, 0.0); m_wFaceSubC.assign((size_t)n_bf*2*3, 0.0);
    m_wFaceSubS.assign((size_t)n_bf*2, 0.0);
    for (int f = 0; f < n_bf; ++f) {
        const int ft = m_wFaceType[f];
        const int nsub = (ft == 0) ? 1 : 2;
        const double* nd = &m_wFaceNodes[(size_t)f*27];
        for (int s = 0; s < nsub; ++s) {
            double Vr[3][2]; WFaceSubTriRef(ft, s, Vr);
            double cen[3] = {0, 0, 0};
            for (int i = 0; i < 3; ++i) {
                double X[3];
                if (ft == 0) TriSurfMap(nd, Vr[i], X); else QuadQ2MapX(nd, Vr[i], X);
                double* out = &m_wFaceSubV[(((size_t)f*2 + s)*3 + i)*3];
                out[0] = X[0]; out[1] = X[1]; out[2] = X[2];
                cen[0] += X[0]/3.0; cen[1] += X[1]/3.0; cen[2] += X[2]/3.0;
            }
            double* pc = &m_wFaceSubC[((size_t)f*2 + s)*3];
            pc[0] = cen[0]; pc[1] = cen[1]; pc[2] = cen[2];
            double sz = 0.0;
            for (int i = 0; i < 3; ++i) {
                const double* v = &m_wFaceSubV[(((size_t)f*2 + s)*3 + i)*3];
                const double dx = v[0]-cen[0], dy = v[1]-cen[1], dz = v[2]-cen[2];
                sz = std::max(sz, std::sqrt(dx*dx + dy*dy + dz*dz));
            }
            m_wFaceSubS[(size_t)f*2 + s] = sz;
        }
    }
    // ---- per-charge host centroid/size (cluster-tree points + the near_hosts test) ----
    m_cent.assign((size_t)m_n*3, 0.0); m_size.assign((size_t)m_n, 0.0);
    for (int a = 0; a < m_n; ++a) {
        const int h = m_host[a];
        double cen[3] = {0, 0, 0};
        int ncorner; double corners[8][3];
        if (m_kind[a] == 0) {                       // 6 prism corners: tri nodes 0,1,2 at iz=0 (n=0,1,2) and iz=2 (n=12,13,14)
            ncorner = 6;
            static const int cidx[6] = {0, 1, 2, 12, 13, 14};
            for (int i = 0; i < 6; ++i) {
                const double* nd = &m_wCellNodes[(size_t)h*54 + 3*cidx[i]];
                for (int k = 0; k < 3; ++k) corners[i][k] = nd[k];
            }
        } else if (m_wFaceType[h] == 0) {           // tri face: 3 corners = tri nodes 0,1,2
            ncorner = 3;
            static const int cidx[3] = {0, 1, 2};
            for (int i = 0; i < 3; ++i) {
                const double* nd = &m_wFaceNodes[(size_t)h*27 + 3*cidx[i]];
                for (int k = 0; k < 3; ++k) corners[i][k] = nd[k];
            }
        } else {                                    // quad face: 4 corners of the 9-lattice
            ncorner = 4;
            static const int cidx[4] = {0, 2, 6, 8};
            for (int i = 0; i < 4; ++i) {
                const double* nd = &m_wFaceNodes[(size_t)h*27 + 3*cidx[i]];
                for (int k = 0; k < 3; ++k) corners[i][k] = nd[k];
            }
        }
        for (int i = 0; i < ncorner; ++i) for (int k = 0; k < 3; ++k) cen[k] += corners[i][k] / ncorner;
        double sz = 0.0;
        for (int i = 0; i < ncorner; ++i) {
            const double dx = corners[i][0]-cen[0], dy = corners[i][1]-cen[1], dz = corners[i][2]-cen[2];
            sz = std::max(sz, std::sqrt(dx*dx + dy*dy + dz*dz));
        }
        m_cent[3*a] = cen[0]; m_cent[3*a+1] = cen[1]; m_cent[3*a+2] = cen[2];
        m_size[a] = sz;
    }
    // ---- (kind,host)->local reverse maps ----
    m_hexLocalOf.assign((size_t)m_n, 0);
    m_cellCharges.assign((size_t)n_el, {}); m_faceCharges.assign((size_t)n_bf, {});
    for (int a = 0; a < m_n; ++a) {
        std::vector<int>& grp = (m_kind[a] == 0) ? m_cellCharges[m_host[a]] : m_faceCharges[m_host[a]];
        m_hexLocalOf[a] = (int)grp.size();
        grp.push_back(a);
    }
    // Structured prism meshes repeat the same wedge cell / boundary-face blocks under translation, but they
    // do not have the hex affine basis used above.  Detect only the conservative case: all host centers fit a
    // rectilinear lattice and each host's node cloud, after subtracting its center, matches one of a small set
    // of translated templates.  Generic/unstructured wedge meshes simply leave the ordinary host-id cache.
    m_hexUniformTransHosts = false;
    m_hexHostTemplate.assign((size_t)n_el + (size_t)n_bf, -1);
    m_hexHostLattice2.assign(((size_t)n_el + (size_t)n_bf)*3, 0);
    const int wedge_trans_scope = WedgeTransCacheScope();
    if (wedge_trans_scope > 0 && n_el > 0) {
        const size_t nhost = (size_t)n_el + (size_t)n_bf;
        std::vector<double> centers(nhost*3, 0.0);
        double lo[3] = {0.0, 0.0, 0.0}, hi[3] = {0.0, 0.0, 0.0};
        bool have_center = false;
        auto record_center = [&](size_t idx, const double* nd, int nnode) {
            double ctr[3] = {0.0, 0.0, 0.0};
            for (int i = 0; i < nnode; ++i)
                for (int k = 0; k < 3; ++k) ctr[k] += nd[3*i + k];
            for (int k = 0; k < 3; ++k) {
                ctr[k] /= (double)nnode;
                centers[3*idx + k] = ctr[k];
                if (!have_center) lo[k] = hi[k] = ctr[k];
                else { lo[k] = std::min(lo[k], ctr[k]); hi[k] = std::max(hi[k], ctr[k]); }
            }
            have_center = true;
        };
        for (int c = 0; c < n_el; ++c)
            record_center((size_t)c, &m_wCellNodes[(size_t)c*54], 18);
        for (int f = 0; f < n_bf; ++f) {
            const int nnode = (m_wFaceType[f] == 0) ? 6 : 9;
            record_center((size_t)n_el + (size_t)f, &m_wFaceNodes[(size_t)f*27], nnode);
        }
        double scale = 0.0;
        for (int k = 0; k < 3; ++k) scale = std::max(scale, hi[k] - lo[k]);
        const double tol = 1e-10 * scale + 1e-12;
        double step[3] = {1.0, 1.0, 1.0};
        bool ok = have_center && scale > 0.0;
        for (int k = 0; ok && k < 3; ++k) {
            std::vector<double> coord;
            coord.reserve(nhost);
            for (size_t i = 0; i < nhost; ++i) coord.push_back(centers[3*i + k]);
            std::sort(coord.begin(), coord.end());
            double min_step = 0.0;
            for (size_t i = 1; i < coord.size(); ++i) {
                const double d = coord[i] - coord[i-1];
                if (d > tol && (min_step == 0.0 || d < min_step)) min_step = d;
            }
            step[k] = (min_step > 0.0) ? min_step : 1.0;
        }
        std::vector<std::vector<double>> cell_templates, face_templates;
        const size_t max_templates = 96;
        auto classify_host = [&](int kind, int h, const double* nd, int nnode, size_t host_index) {
            const double* ctr = &centers[3*host_index];
            for (int k = 0; k < 3; ++k) {
                const double q = 2.0 * (ctr[k] - lo[k]) / step[k];
                if (q < -2.0e9 || q > 2.0e9) return false;
                const int qi2 = (int)std::llround(q);
                if (std::fabs(q - (double)qi2) > 2e-6) return false;
                m_hexHostLattice2[3*host_index + k] = qi2;
            }
            const std::vector<int>& grp = (kind == 0) ? m_cellCharges[h] : m_faceCharges[h];
            std::vector<double> sig;
            sig.reserve((size_t)nnode*3 + 1 + (size_t)grp.size()*3);
            for (int i = 0; i < nnode; ++i)
                for (int k = 0; k < 3; ++k) sig.push_back(nd[3*i + k] - ctr[k]);
            sig.push_back((double)grp.size());
            for (int g : grp)
                for (int k = 0; k < 3; ++k) sig.push_back((double)m_expo[(size_t)3*g + k]);
            auto& templ = (kind == 0) ? cell_templates : face_templates;
            int tid = -1;
            for (int t = 0; t < (int)templ.size(); ++t) {
                if (templ[t].size() != sig.size()) continue;
                bool same = true;
                for (size_t i = 0; i < sig.size(); ++i)
                    if (std::fabs(templ[t][i] - sig[i]) > tol) { same = false; break; }
                if (same) { tid = t; break; }
            }
            if (tid < 0) {
                if (templ.size() >= max_templates) return false;
                tid = (int)templ.size();
                templ.push_back(std::move(sig));
            }
            m_hexHostTemplate[host_index] = tid;
            return true;
        };
        for (int c = 0; ok && c < n_el; ++c)
            ok = classify_host(0, c, &m_wCellNodes[(size_t)c*54], 18, (size_t)c);
        if (wedge_trans_scope >= 2) {
            for (int f = 0; ok && f < n_bf; ++f) {
                const int nnode = (m_wFaceType[f] == 0) ? 6 : 9;
                ok = classify_host(1, f, &m_wFaceNodes[(size_t)f*27], nnode, (size_t)n_el + (size_t)f);
            }
        }
        m_hexUniformTransHosts = ok;
    }
    BuildWedgeSiteTables();
}

// Ref coords of anchor site k of cell sub-tet s (hex-ref frame): 0-3 corners, 4-9 edge midpoints
// ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3)), 10-13 face centers (HEXTET_FC order), 14 centroid.
static void HexSiteRef(int s, int k, double x0[3])
{
    const int* tv = HEXREF_TETS[s];
    double V[4][3];
    for (int i = 0; i < 4; ++i) for (int d = 0; d < 3; ++d) V[i][d] = HEXREF_V[tv[i]][d];
    static const int E[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
    if (k < 4)       for (int d = 0; d < 3; ++d) x0[d] = V[k][d];
    else if (k < 10) for (int d = 0; d < 3; ++d) x0[d] = 0.5*(V[E[k-4][0]][d] + V[E[k-4][1]][d]);
    else if (k < 14) {
        const int* f = HEXTET_FC[k-10];
        for (int d = 0; d < 3; ++d) x0[d] = (V[f[0]][d] + V[f[1]][d] + V[f[2]][d])/3.0;
    } else            for (int d = 0; d < 3; ++d) x0[d] = 0.25*(V[0][d]+V[1][d]+V[2][d]+V[3][d]);
}

// Ref uv coords of anchor site k of face sub-tri s: 0-2 corners, 3-5 edge midpoints ((0,1),(1,2),(2,0)),
// 6 centroid.
static void QuadSiteRef(int s, int k, double u0[2])
{
    const int* tv = QUADREF_TRIS[s];
    double V[3][2];
    for (int i = 0; i < 3; ++i) for (int d = 0; d < 2; ++d) V[i][d] = QUADREF_V[tv[i]][d];
    if (k < 3)      for (int d = 0; d < 2; ++d) u0[d] = V[k][d];
    else if (k < 6) for (int d = 0; d < 2; ++d) u0[d] = 0.5*(V[k-3][d] + V[(k-2)%3][d]);
    else            for (int d = 0; d < 2; ++d) u0[d] = (V[0][d]+V[1][d]+V[2][d])/3.0;
}

// Build the host-INDEPENDENT static-site radial tables (see the header doc): for each (cell sub, site) /
// (face sub, site), the radial-cone nodes from the site are FIXED ref points, so the Q2 shape values S,
// the Q1 monomial values M, and the signed Piola weights w are precomputed once; a call is then one
// nq x 27|9 "GEMV" (X = S @ nodes) + nq kernel evals.  Cones whose base face contains the site have D = 0
// and are skipped (corner sites keep 1 of 4 cones, edge mids 2, face centers 3, centroid 4).  Also fills
// the per-host MAPPED site positions used by the nearest-site pick.
void RadHACApKChargeGram::BuildHexSiteTables()
{
    const int nR = (int)m_glIn.size();
    const double* GL = m_glIn.data();
    const double* GW = m_gwIn.data();
    m_cellSiteRad.assign(6*15, HexSiteRad());
    for (int s = 0; s < 6; ++s) {
        const int* tv = HEXREF_TETS[s];
        double V[4][3];
        for (int i = 0; i < 4; ++i) for (int d = 0; d < 3; ++d) V[i][d] = HEXREF_V[tv[i]][d];
        double E0[3], E1[3], E2[3];
        for (int d = 0; d < 3; ++d) { E0[d] = V[1][d]-V[0][d]; E1[d] = V[2][d]-V[0][d]; E2[d] = V[3][d]-V[0][d]; }
        const double hv = E0[0]*(E1[1]*E2[2]-E1[2]*E2[1]) - E0[1]*(E1[0]*E2[2]-E1[2]*E2[0])
                        + E0[2]*(E1[0]*E2[1]-E1[1]*E2[0]);
        const double sgnT = (hv >= 0.0) ? 1.0 : -1.0;
        for (int k = 0; k < 15; ++k) {
            HexSiteRad& R = m_cellSiteRad[(size_t)s*15 + k];
            double x0[3];
            HexSiteRef(s, k, x0);
            for (int f = 0; f < 4; ++f) {
                const double* b1 = V[HEXTET_FC[f][0]];
                const double* b2 = V[HEXTET_FC[f][1]];
                const double* b3 = V[HEXTET_FC[f][2]];
                double d1[3], d2[3], d3[3], e21[3], e32[3];
                for (int d = 0; d < 3; ++d) {
                    d1[d] = b1[d]-x0[d]; d2[d] = b2[d]-x0[d]; d3[d] = b3[d]-x0[d];
                    e21[d] = b2[d]-b1[d]; e32[d] = b3[d]-b2[d];
                }
                const double cr[3] = {d2[1]*d3[2]-d2[2]*d3[1], d2[2]*d3[0]-d2[0]*d3[2], d2[0]*d3[1]-d2[1]*d3[0]};
                const double D = d1[0]*cr[0] + d1[1]*cr[1] + d1[2]*cr[2];
                if (std::fabs(D) < 1e-12) continue;              // degenerate cone: site lies on this face
                for (int a = 0; a < nR; ++a) { const double u = GL[a];
                    for (int b = 0; b < nR; ++b) { const double v = GL[b];
                        for (int c = 0; c < nR; ++c) { const double w = GL[c];
                            double y[3];
                            for (int d = 0; d < 3; ++d) y[d] = x0[d] + u*(d1[d] + v*(e21[d] + w*e32[d]));
                            R.w.push_back(sgnT*GW[a]*GW[b]*GW[c]*(u*u*v*D));
                            double vx[3], dxu[3], vy[3], dyu[3], vz[3], dzu[3];
                            HexLag3(y[0], vx, dxu); HexLag3(y[1], vy, dyu); HexLag3(y[2], vz, dzu);
                            for (int iz = 0; iz < 3; ++iz)
                                for (int iy = 0; iy < 3; ++iy)
                                    for (int ix = 0; ix < 3; ++ix) R.S.push_back(vx[ix]*vy[iy]*vz[iz]);
                            const double px[3] = {1.0, y[0], y[0]*y[0]};
                            const double py[3] = {1.0, y[1], y[1]*y[1]};
                            const double pz[3] = {1.0, y[2], y[2]*y[2]};
                            for (int ez = 0; ez < 3; ++ez)
                                for (int ey = 0; ey < 3; ++ey)
                                    for (int ex = 0; ex < 3; ++ex)
                                        R.M.push_back(px[ex]*py[ey]*pz[ez]);
                        }
                    }
                }
            }
            R.nq = (int)R.w.size();
        }
    }
    m_faceSiteRad.assign(2*7, HexSiteRad());
    for (int s = 0; s < 2; ++s) {
        const int* tv = QUADREF_TRIS[s];
        double V[3][2];
        for (int i = 0; i < 3; ++i) for (int d = 0; d < 2; ++d) V[i][d] = QUADREF_V[tv[i]][d];
        for (int k = 0; k < 7; ++k) {
            HexSiteRad& R = m_faceSiteRad[(size_t)s*7 + k];
            double u0[2];
            QuadSiteRef(s, k, u0);
            for (int kf = 0; kf < 3; ++kf) {
                const double* A = V[kf]; const double* B = V[(kf+1)%3];
                const double ea[2] = {A[0]-u0[0], A[1]-u0[1]};
                const double eb[2] = {B[0]-u0[0], B[1]-u0[1]};
                const double s2 = ea[0]*eb[1] - ea[1]*eb[0];
                if (std::fabs(s2) < 1e-12) continue;             // degenerate cone: site lies on this edge
                for (int a = 0; a < nR; ++a) { const double u = GL[a];
                    for (int b = 0; b < nR; ++b) { const double v = GL[b];
                        const double yu = u0[0] + u*(ea[0] + v*(eb[0]-ea[0]));
                        const double yv = u0[1] + u*(ea[1] + v*(eb[1]-ea[1]));
                        R.w.push_back(GW[a]*GW[b]*(u*s2));       // QUADREF_TRIS are CCW: signed s2 sums to +
                        double vu[3], duu[3], vv[3], dvu[3];
                        HexLag3(yu, vu, duu); HexLag3(yv, vv, dvu);
                        for (int iv = 0; iv < 3; ++iv)
                            for (int iu = 0; iu < 3; ++iu) R.S.push_back(vu[iu]*vv[iv]);
                        const double pu[3] = {1.0, yu, yu*yu};
                        const double pv[3] = {1.0, yv, yv*yv};
                        for (int ev = 0; ev < 3; ++ev)
                            for (int eu = 0; eu < 3; ++eu) R.M.push_back(pu[eu]*pv[ev]);
                    }
                }
            }
            R.nq = (int)R.w.size();
        }
    }
    // ---- mapped site positions per host (the nearest-site pick is a physical distance test) ----
    m_cellSiteX.assign((size_t)m_n_el*6*15*3, 0.0);
    for (int c = 0; c < m_n_el; ++c) {
        const double* nd = &m_hexNodes[(size_t)c*81];
        for (int s = 0; s < 6; ++s)
            for (int k = 0; k < 15; ++k) {
                double x0[3], X[3];
                HexSiteRef(s, k, x0);
                HexQ2MapX(nd, x0, X);
                double* out = &m_cellSiteX[(((size_t)c*6 + s)*15 + k)*3];
                out[0] = X[0]; out[1] = X[1]; out[2] = X[2];
            }
    }
    m_faceSiteX.assign((size_t)m_hex_n_bf*2*7*3, 0.0);
    for (int f = 0; f < m_hex_n_bf; ++f) {
        const double* nd = &m_quadNodes[(size_t)f*27];
        for (int s = 0; s < 2; ++s)
            for (int k = 0; k < 7; ++k) {
                double u0[2], X[3];
                QuadSiteRef(s, k, u0);
                QuadQ2MapX(nd, u0, X);
                double* out = &m_faceSiteX[(((size_t)f*2 + s)*7 + k)*3];
                out[0] = X[0]; out[1] = X[1]; out[2] = X[2];
            }
    }
    m_hex_state_sum = HexStateChecksum();   // heap-stomp canary: everything a block compute reads
}

// Checksum of every hex-mode member array the block computation reads (heap-stomp canary; see header).
double RadHACApKChargeGram::HexStateChecksum() const
{
    double s = 0.0;
    for (const auto& kv : HexStateBreakdown()) s += kv.second;
    return s;
}

// Per-array checksum breakdown (flake forensics: which array differs between two instances).
std::vector<std::pair<std::string, double>> RadHACApKChargeGram::HexStateBreakdown() const
{
    std::vector<std::pair<std::string, double>> out;
    auto add = [&out](const char* name, const std::vector<double>& v) {
        double s = 0.0;
        for (double x : v) s += x;
        out.emplace_back(name, s);
    };
    auto addi = [&out](const char* name, const std::vector<int>& v) {
        double s = 0.0;
        for (int x : v) s += (double)x;
        out.emplace_back(name, s);
    };
    add("hexNodes", m_hexNodes); add("quadNodes", m_quadNodes);
    add("symTetP", m_symTetP); add("symTetW", m_symTetW);
    add("symTriP", m_symTriP); add("symTriW", m_symTriW);
    add("glOut", m_glOut); add("gwOut", m_gwOut); add("glIn", m_glIn); add("gwIn", m_gwIn);
    add("farTetP", m_farTetP); add("farTetW", m_farTetW);
    add("farTriP", m_farTriP); add("farTriW", m_farTriW);
    add("cellSubC", m_cellSubC); add("cellSubS", m_cellSubS); add("cellSubV", m_cellSubV);
    add("faceSubC", m_faceSubC); add("faceSubS", m_faceSubS); add("faceSubV", m_faceSubV);
    {
        double s = 0.0;
        for (unsigned char x : m_hexAffineCell) s += (double)x;
        out.emplace_back("hexAffineCell", s);
    }
    out.emplace_back("hexAffineOrder", (double)m_hexAffineOrder);
    out.emplace_back("hexAffineMonoCount", (double)m_hexAffineMonoCount);
    out.emplace_back("hexAffinePolyCount", (double)m_hexAffinePolyCount);
    add("hexAffineCoeff", m_hexAffineCoeff);
    {
        double s = 0.0;
        for (unsigned char x : m_quadAffineFace) s += (double)x;
        out.emplace_back("quadAffineFace", s);
    }
    out.emplace_back("quadAffineMonoCount", (double)m_quadAffineMonoCount);
    out.emplace_back("quadAffinePolyCount", (double)m_quadAffinePolyCount);
    add("quadAffineCoeff", m_quadAffineCoeff);
    out.emplace_back("hexUniformAffineCells", m_hexUniformAffineCells ? 1.0 : 0.0);
    addi("hexCellLattice", m_hexCellLattice);
    out.emplace_back("hexUniformTransHosts", m_hexUniformTransHosts ? 1.0 : 0.0);
    addi("hexHostTemplate", m_hexHostTemplate);
    addi("hexHostLattice2", m_hexHostLattice2);
    add("cent", m_cent); add("size", m_size);
    addi("host", m_host); addi("kind", m_kind); addi("expo", m_expo); addi("hexLocalOf", m_hexLocalOf);
    {
        double s = 0.0;
        for (const HexSiteRad& R : m_cellSiteRad) { s += R.nq; for (double x : R.S) s += x; for (double x : R.M) s += x; for (double x : R.w) s += x; }
        out.emplace_back("cellSiteRad", s);
    }
    {
        double s = 0.0;
        for (const HexSiteRad& R : m_faceSiteRad) { s += R.nq; for (double x : R.S) s += x; for (double x : R.M) s += x; for (double x : R.w) s += x; }
        out.emplace_back("faceSiteRad", s);
    }
    add("cellSiteX", m_cellSiteX); add("faceSiteX", m_faceSiteX);
    // 2D planar mode arrays (empty in the hex mode and vice versa)
    add("d2CellMap", m_d2CellMap); add("d2EdgeMap", m_d2EdgeMap);
    out.emplace_back("d2GeometryOrder", static_cast<double>(m_d2GeometryOrder));
    addi("d2CellType", m_d2CellType);
    add("d2SymTriP", m_d2SymTriP); add("d2SymTriW", m_d2SymTriW);
    add("d2GlE", m_d2GlE); add("d2GwE", m_d2GwE);
    add("d2FarTriP", m_d2FarTriP); add("d2FarTriW", m_d2FarTriW);
    add("d2CellSubC", m_d2CellSubC); add("d2CellSubS", m_d2CellSubS);
    add("d2EdgeC", m_d2EdgeC); add("d2EdgeS", m_d2EdgeS);
    add("d2CellSiteX", m_d2CellSiteX); add("d2EdgeSiteX", m_d2EdgeSiteX);
    // WEDGE (PRISM) mode arrays (empty in the hex/2D modes and vice versa -> the hex checksum is unchanged)
    add("wCellNodes", m_wCellNodes); add("wFaceNodes", m_wFaceNodes); addi("wFaceType", m_wFaceType);
    add("wCellSubC", m_wCellSubC); add("wCellSubS", m_wCellSubS); add("wCellSubV", m_wCellSubV);
    add("wFaceSubC", m_wFaceSubC); add("wFaceSubS", m_wFaceSubS); add("wFaceSubV", m_wFaceSubV);
    {
        double s = 0.0;
        for (const HexSiteRad& R : m_wCellSiteRad)     { s += R.nq; for (double x : R.S) s += x; for (double x : R.M) s += x; for (double x : R.w) s += x; }
        for (const HexSiteRad& R : m_wFaceSiteRadTri)  { s += R.nq; for (double x : R.S) s += x; for (double x : R.M) s += x; for (double x : R.w) s += x; }
        for (const HexSiteRad& R : m_wFaceSiteRadQuad) { s += R.nq; for (double x : R.S) s += x; for (double x : R.M) s += x; for (double x : R.w) s += x; }
        out.emplace_back("wSiteRad", s);
    }
    add("wCellSiteX", m_wCellSiteX); add("wFaceSiteX", m_wFaceSiteX);
    return out;
}

void RadHACApKChargeGram::ResetHexCacheStats()
{
    m_hexBlockLookups.store(0, std::memory_order_relaxed);
    m_hexBlockHits.store(0, std::memory_order_relaxed);
    m_hexBlockMisses.store(0, std::memory_order_relaxed);
    m_hexBlockClears.store(0, std::memory_order_relaxed);
    m_hexTransBlockLookups.store(0, std::memory_order_relaxed);
    m_hexTransBlockHits.store(0, std::memory_order_relaxed);
    m_hexTransBlockMisses.store(0, std::memory_order_relaxed);
    m_hexTransBlockClears.store(0, std::memory_order_relaxed);
    m_hexSymBlockLookups.store(0, std::memory_order_relaxed);
    m_hexSymBlockHits.store(0, std::memory_order_relaxed);
    m_hexSymBlockMisses.store(0, std::memory_order_relaxed);
    m_hexSymBlockClears.store(0, std::memory_order_relaxed);
    m_hexSymTransBlockLookups.store(0, std::memory_order_relaxed);
    m_hexSymTransBlockHits.store(0, std::memory_order_relaxed);
    m_hexSymTransBlockMisses.store(0, std::memory_order_relaxed);
    m_hexSymTransBlockClears.store(0, std::memory_order_relaxed);
    m_hoSymBlockLookups.store(0, std::memory_order_relaxed);
    m_hoSymBlockHits.store(0, std::memory_order_relaxed);
    m_hoSymBlockMisses.store(0, std::memory_order_relaxed);
    m_hoSymBlockClears.store(0, std::memory_order_relaxed);
    m_hexBlkAffineNear.store(0, std::memory_order_relaxed);
    m_hexBlkAffineFar.store(0, std::memory_order_relaxed);
    m_hexBlkDistortedFar.store(0, std::memory_order_relaxed);
    m_hexBlkGeneralNear.store(0, std::memory_order_relaxed);
    m_hexBlkGeneralFar.store(0, std::memory_order_relaxed);
    m_hexNsAffineNear.store(0, std::memory_order_relaxed);
    m_hexNsAffineFar.store(0, std::memory_order_relaxed);
    m_hexNsDistortedFar.store(0, std::memory_order_relaxed);
    m_hexNsGeneralNear.store(0, std::memory_order_relaxed);
    m_hexNsGeneralFar.store(0, std::memory_order_relaxed);
    m_hexGeneralSharedLookups.store(0, std::memory_order_relaxed);
    m_hexGeneralSharedHits.store(0, std::memory_order_relaxed);
    m_hexGeneralSharedMisses.store(0, std::memory_order_relaxed);
}

std::vector<std::pair<std::string, double>> RadHACApKChargeGram::HexCacheStats() const
{
    auto ld = [](const std::atomic<long long>& v) { return (double)v.load(std::memory_order_relaxed); };
    std::vector<std::pair<std::string, double>> out;
    out.emplace_back("hex_cache_stats_enabled", m_hexCacheStatsEnabled ? 1.0 : 0.0);
    out.emplace_back("hex_far_one_sided_threshold", HexFarOneSidedThreshold());
    out.emplace_back("hex_affine_exact_near_factor", HEX_AFFINE_EXACT_NEAR_FACTOR);
    out.emplace_back("hex_distorted_far_factor", HexDistortedFarFactor());
    // QuadBlockHex dispatch profile: computed blocks + accumulated wall seconds per branch (thread-summed
    // across the ParallelFor workers, so the seconds compare branch-to-branch, not to the build wall clock).
    out.emplace_back("hex_blk_affine_near", ld(m_hexBlkAffineNear));
    out.emplace_back("hex_blk_affine_far", ld(m_hexBlkAffineFar));
    out.emplace_back("hex_blk_distorted_far", ld(m_hexBlkDistortedFar));
    out.emplace_back("hex_blk_general_near", ld(m_hexBlkGeneralNear));
    out.emplace_back("hex_blk_general_far", ld(m_hexBlkGeneralFar));
    out.emplace_back("hex_s_affine_near", ld(m_hexNsAffineNear) * 1e-9);
    out.emplace_back("hex_s_affine_far", ld(m_hexNsAffineFar) * 1e-9);
    out.emplace_back("hex_s_distorted_far", ld(m_hexNsDistortedFar) * 1e-9);
    out.emplace_back("hex_s_general_near", ld(m_hexNsGeneralNear) * 1e-9);
    out.emplace_back("hex_s_general_far", ld(m_hexNsGeneralFar) * 1e-9);
    out.emplace_back("hex_general_shared_lookups", ld(m_hexGeneralSharedLookups));
    out.emplace_back("hex_general_shared_hits", ld(m_hexGeneralSharedHits));
    out.emplace_back("hex_general_shared_misses", ld(m_hexGeneralSharedMisses));
    out.emplace_back("wedge_far_one_sided_threshold", WedgeFarOneSidedThreshold());
    out.emplace_back("wedge_trans_cache_scope", (double)WedgeTransCacheScope());
    out.emplace_back("wedge_trans_cache_enabled", WedgeTransCacheScope() > 0 ? 1.0 : 0.0);
    out.emplace_back("ho_far_one_sided_enabled", HOFarOneSidedEnabled() ? 1.0 : 0.0);
    out.emplace_back("ho_analytic_block_available", m_hoAnalyticBlock ? 1.0 : 0.0);
    out.emplace_back("ho_analytic_block_enabled",
                     m_hoAnalyticBlock && HOAnalyticBlockEnabled() ? 1.0 : 0.0);
    out.emplace_back("hex_block_lookups", ld(m_hexBlockLookups));
    out.emplace_back("hex_block_hits", ld(m_hexBlockHits));
    out.emplace_back("hex_block_misses", ld(m_hexBlockMisses));
    out.emplace_back("hex_block_clears", ld(m_hexBlockClears));
    out.emplace_back("hex_trans_block_lookups", ld(m_hexTransBlockLookups));
    out.emplace_back("hex_trans_block_hits", ld(m_hexTransBlockHits));
    out.emplace_back("hex_trans_block_misses", ld(m_hexTransBlockMisses));
    out.emplace_back("hex_trans_block_clears", ld(m_hexTransBlockClears));
    out.emplace_back("hex_sym_block_lookups", ld(m_hexSymBlockLookups));
    out.emplace_back("hex_sym_block_hits", ld(m_hexSymBlockHits));
    out.emplace_back("hex_sym_block_misses", ld(m_hexSymBlockMisses));
    out.emplace_back("hex_sym_block_clears", ld(m_hexSymBlockClears));
    out.emplace_back("hex_sym_trans_block_lookups", ld(m_hexSymTransBlockLookups));
    out.emplace_back("hex_sym_trans_block_hits", ld(m_hexSymTransBlockHits));
    out.emplace_back("hex_sym_trans_block_misses", ld(m_hexSymTransBlockMisses));
    out.emplace_back("hex_sym_trans_block_clears", ld(m_hexSymTransBlockClears));
    out.emplace_back("ho_sym_block_lookups", ld(m_hoSymBlockLookups));
    out.emplace_back("ho_sym_block_hits", ld(m_hoSymBlockHits));
    out.emplace_back("ho_sym_block_misses", ld(m_hoSymBlockMisses));
    out.emplace_back("ho_sym_block_clears", ld(m_hoSymBlockClears));
    const double btot = ld(m_hexBlockLookups);
    const double ttot = ld(m_hexTransBlockLookups);
    const double sbtot = ld(m_hexSymBlockLookups);
    const double sttot = ld(m_hexSymTransBlockLookups);
    const double hotot = ld(m_hoSymBlockLookups);
    out.emplace_back("hex_block_hit_rate", btot > 0.0 ? ld(m_hexBlockHits) / btot : 0.0);
    out.emplace_back("hex_trans_block_hit_rate", ttot > 0.0 ? ld(m_hexTransBlockHits) / ttot : 0.0);
    out.emplace_back("hex_sym_block_hit_rate", sbtot > 0.0 ? ld(m_hexSymBlockHits) / sbtot : 0.0);
    out.emplace_back("hex_sym_trans_block_hit_rate", sttot > 0.0 ? ld(m_hexSymTransBlockHits) / sttot : 0.0);
    out.emplace_back("ho_sym_block_hit_rate", hotot > 0.0 ? ld(m_hoSymBlockHits) / hotot : 0.0);
    out.emplace_back("curved_touch_blocks", (double)m_curvedTouchBlocks.size());
    out.emplace_back("curved_touch_build_time", m_curvedTouchBuildTime);
    return out;
}

// A materialized quadrature cloud on one sub-simplex: physical points, geometry weights (rule weight x
// scale x |det J| -- everything EXCEPT the charge monomial), and the hex/quad REF coords (for the
// per-charge monomial).  Cached per (kind, host, sub, corner/rule): the cloud depends only on geometry,
// so it is reused across ALL outer points selecting the same grading corner AND all co-located charges
// (the numpy-validated src_cache pattern; ~2 orders of magnitude fewer Q2-map evals on near pairs).
struct HexQuadCloud { std::vector<double> pts, wgeo, xi; };

// Materialize the cloud for sub-simplex `sub` of the host with nodes `nd` from a bary rule.  full_bary:
// the rule stores nv coords/point (graded Duffy); else nv-1 lam coords (the fixed far/sym tables).
static void HexBuildCloud(const double* nd, bool cell, int sub, const double* baryP, const double* baryW,
                          int nq, bool full_bary, HexQuadCloud& out)
{
    const int* tv = cell ? HEXREF_TETS[sub] : QUADREF_TRIS[sub];
    const int nv = cell ? 4 : 3;
    const double scale = cell ? HexSubSixVref(sub) : QuadSubTwoAref(sub);
    out.pts.resize((size_t)nq*3); out.wgeo.resize(nq); out.xi.resize((size_t)nq*3);
    for (int q = 0; q < nq; ++q) {
        double bary[4];
        if (full_bary) {
            for (int t = 0; t < nv; ++t) bary[t] = baryP[(size_t)nv*q + t];
        } else {
            double lsum = 0.0;
            for (int t = 1; t < nv; ++t) { bary[t] = baryP[(size_t)(nv-1)*q + (t-1)]; lsum += bary[t]; }
            bary[0] = 1.0 - lsum;
        }
        if (cell) {
            double xi[3] = {0, 0, 0};
            for (int t = 0; t < 4; ++t)
                for (int k = 0; k < 3; ++k) xi[k] += bary[t]*HEXREF_V[tv[t]][k];
            double X[3], J[3][3];
            RadHACApKChargeGram::HexQ2Map(nd, xi, X, J);
            for (int k = 0; k < 3; ++k) { out.pts[(size_t)3*q+k] = X[k]; out.xi[(size_t)3*q+k] = xi[k]; }
            out.wgeo[q] = baryW[q]*scale;              // REF measure (Piola-exact charge: no |det J|)
        } else {
            double uv[2] = {0, 0};
            for (int t = 0; t < 3; ++t)
                for (int k = 0; k < 2; ++k) uv[k] += bary[t]*QUADREF_V[tv[t]][k];
            double X[3], T[3][2];
            RadHACApKChargeGram::QuadQ2Map(nd, uv, X, T);
            out.pts[(size_t)3*q] = X[0]; out.pts[(size_t)3*q+1] = X[1]; out.pts[(size_t)3*q+2] = X[2];
            out.xi[(size_t)3*q] = uv[0]; out.xi[(size_t)3*q+1] = uv[1]; out.xi[(size_t)3*q+2] = 0.0;
            out.wgeo[q] = baryW[q]*scale;              // REF measure (Piola-exact charge: no surf J)
        }
    }
}

// thread_local cloud-cache key (build_id-guarded like the QuadDot memo; see NextChargeGramBuildId).
// kind(1b at 62) | outer(1b at 61) | graded(1b at 60) | host(<<8) | sub(<<2) | corner (far rule: corner=3
// non-graded; graded inner/outer set bit 60, corner in 0..3).
static inline long long HexCloudKey(int kind, bool outer, bool graded, int host, int sub, int corner)
{
    return ((long long)kind << 62) | ((long long)(outer ? 1 : 0) << 61) | ((long long)(graded ? 1 : 0) << 60)
         | ((long long)host << 8) | ((long long)sub << 2) | (long long)corner;
}

// SHARED_PTR values (2026-07-03 crash fix): the capacity clear below fires on ~20k-charge meshes (a
// 1000-hex cube wants ~43k outer clouds > the 32768 cap; <=8^3 stays under -- which is why the bug slept
// through every gate).  QuadBlockHex HOLDS its outer cloud across inner calls that fetch far clouds, so a
// by-value cache whose clear() destroys storage turned that hold into a use-after-free (0xC0000005 at
// n=10, reproduced 2/2 on the committed binary).  shared_ptr makes the clear safe: in-flight holders keep
// their cloud alive; the cache only drops its refs.
static thread_local long long s_hex_cloud_owner = -1;
static thread_local std::unordered_map<long long, std::shared_ptr<const HexQuadCloud>> s_hex_cloud_cache;

static std::shared_ptr<const HexQuadCloud> HexGetCloud(long long build_id, long long key,
                                                       const std::function<void(HexQuadCloud&)>& make)
{
    if (s_hex_cloud_owner != build_id) { s_hex_cloud_cache.clear(); s_hex_cloud_owner = build_id; }
    auto it = s_hex_cloud_cache.find(key);
    if (it == s_hex_cloud_cache.end()) {
        if (s_hex_cloud_cache.size() > 32768u) s_hex_cloud_cache.clear();   // safe: holders own shared_ptrs
        auto c = std::make_shared<HexQuadCloud>();
        make(*c);
        it = s_hex_cloud_cache.emplace(key, std::move(c)).first;
    }
    return it->second;
}

// Vectorized inner: INT over sub `subB` of src host (kindS,hS) of mono_b(y)/|p-y| dy for ALL source local
// charges srcG[], accumulated into inn[ls].  FAR field point -> the cheap cached far cloud (smooth 1/r);
// NEAR -> the static-SITE radial (PhiInnerHexSiteVec: precomputed ref tables anchored at the nearest
// site) -- the same exact radial cone tiling as the self path, served at shape-"GEMV" cost.
// NON-SELF near inner: the static-SITE radial (see the header doc).  Nearest mapped site anchors the
// precomputed ref-space radial tables; the per-call work is one nq x 27|9 shape "GEMV" + nq kernel evals.
void RadHACApKChargeGram::PhiInnerHexSiteVec(int kindS, int hS, int subB, const double p[3],
                                             const std::vector<int>& srcG, double* inn) const
{
    const bool cell = (kindS == 0);
    const int nsite = cell ? 15 : 7;
    const double* sx = cell ? &m_cellSiteX[(((size_t)hS*6 + subB)*15)*3]
                            : &m_faceSiteX[(((size_t)hS*2 + subB)*7)*3];
    int best = 0; double bd = 1e300;
    for (int k = 0; k < nsite; ++k) {
        const double dx = p[0]-sx[3*k], dy = p[1]-sx[3*k+1], dz = p[2]-sx[3*k+2];
        const double d = dx*dx + dy*dy + dz*dz;
        if (d < bd) { bd = d; best = k; }
    }
    const HexSiteRad& R = cell ? m_cellSiteRad[(size_t)subB*15 + best] : m_faceSiteRad[(size_t)subB*7 + best];
    const double* nd = cell ? &m_hexNodes[(size_t)hS*81] : &m_quadNodes[(size_t)hS*27];
    const int nn = cell ? 27 : 9;
    const int nm = cell ? 27 : 9;
    const int nS = (int)srcG.size();
    std::vector<int> col((size_t)nS);
    for (int ls = 0; ls < nS; ++ls) {
        const int* e = &m_expo[(size_t)3*srcG[ls]];
        col[ls] = e[0] + 3*e[1] + (cell ? 9*e[2] : 0);
    }
    for (int q = 0; q < R.nq; ++q) {
        const double* Sq = &R.S[(size_t)q*nn];
        double X0 = 0.0, X1 = 0.0, X2 = 0.0;
        for (int n2 = 0; n2 < nn; ++n2) {
            const double s = Sq[n2]; const double* v = &nd[3*n2];
            X0 += s*v[0]; X1 += s*v[1]; X2 += s*v[2];
        }
        const double dx = p[0]-X0, dy = p[1]-X1, dz = p[2]-X2;
        const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 1e-300) continue;
        const double g = R.w[q]/r;
        const double* Mq = &R.M[(size_t)q*nm];
        for (int ls = 0; ls < nS; ++ls) inn[ls] += g*Mq[col[ls]];
    }
}

void RadHACApKChargeGram::PhiInnerHexAffineCellSubVec(int hS, int subB, const double p[3],
                                                      const std::vector<int>& srcG, double* inn) const
{
    const size_t sid = (size_t)hS*6 + subB;
    const double* sv = &m_cellSubV[sid*4*3];
    double V[4][3];
    for (int i = 0; i < 4; ++i)
        for (int k = 0; k < 3; ++k) V[i][k] = sv[3*i + k];
    double moments[HEX_AFFINE_POLY_N];
    if (m_hexAffineOrder == 1) rad_hdiv::TetPotentialMomentsUpTo3(V, p, moments);
    else                       rad_hdiv::TetPotentialMomentsUpTo6(V, p, moments);
    for (int ls = 0; ls < (int)srcG.size(); ++ls) {
        const int* e = &m_expo[(size_t)3*srcG[ls]];
        const int axisCount = m_hexAffineOrder + 1;
        const int mono = e[0] + axisCount*e[1] + axisCount*axisCount*e[2];
        const double* coeff = &m_hexAffineCoeff[
            ((size_t)hS*m_hexAffineMonoCount + mono)*m_hexAffinePolyCount];
        double s = 0.0;
        for (int i = 0; i < m_hexAffinePolyCount; ++i) s += coeff[i] * moments[i];
        inn[ls] += s;
    }
}
void RadHACApKChargeGram::PhiInnerHexAffineFaceSubVec(int hS, int subB, const double p[3],
                                                      const std::vector<int>& srcG, double* inn) const
{
    // A flat Q2 face is affine on the whole quad.  Convert all Q2 reference monomials to physical
    // polynomials (degree <= 4) once in the constructor and apply the same analytic surface-moment
    // strategy used by the flat-TET BDM2 path.
    const int* tv = QUADREF_TRIS[subB];
    double V[3][3];
    for (int i = 0; i < 3; ++i) {
        const double* v = &m_faceSubV[(((size_t)hS * 2 + subB) * 3 + i) * 3];
        for (int k = 0; k < 3; ++k) V[i][k] = v[k];
    }
    double moments[QUAD_AFFINE_POLY_N];
    if (m_hexAffineOrder == 1) rad_hdiv::TriPotentialMomentsUpTo2(V, p, moments);
    else                       rad_hdiv::TriPotentialMomentsUpTo4(V, p, moments);
    for (int ls = 0; ls < (int)srcG.size(); ++ls) {
        const int* e = &m_expo[(size_t)3*srcG[ls]];
        const int axisCount = m_hexAffineOrder + 1;
        const int mono = e[0] + axisCount*e[1];
        const double* coeff = &m_quadAffineCoeff[
            ((size_t)hS*m_quadAffineMonoCount + mono)*m_quadAffinePolyCount];
        double value = 0.0;
        for (int i = 0; i < m_quadAffinePolyCount; ++i) value += coeff[i]*moments[i];
        inn[ls] += value;
    }
}

void RadHACApKChargeGram::PhiInnerHexAffineCellVec(int hS, const double p[3],
                                                   const std::vector<int>& srcG, double* inn) const
{
    for (int subB = 0; subB < 6; ++subB)
        PhiInnerHexAffineCellSubVec(hS, subB, p, srcG, inn);
}

void RadHACApKChargeGram::PhiInnerHexAffineFaceVec(int hS, const double p[3],
                                                   const std::vector<int>& srcG, double* inn) const
{
    for (int subB = 0; subB < 2; ++subB)
        PhiInnerHexAffineFaceSubVec(hS, subB, p, srcG, inn);
}

// Smooth affine host pairs do not need degree-six analytic potential recurrences at every outer point.
// Integrate both complete reference hosts with the same tensor Gauss rule instead.  The cube/quad rule is
// invariant under axis reflections (unlike a fixed sub-tet diagonal), and all local charge modes share the
// kernel evaluation.  Near/self pairs stay on QuadBlockHexAffineProduct's analytic source potential.
// NOT affine-only despite the name: points are placed by the full Q2 map (HexQ2MapX/QuadQ2MapX) and the
// Piola reference charge measure makes the Jacobian drop out, so the same rule serves well-separated
// DISTORTED/curved pairs too (the HexDistortedFarFactor dispatch in QuadBlockHex, C-4 fill speedup).
const RadHACApKChargeGram::HexFarRule& RadHACApKChargeGram::GetHexFarRule(int kind, int host) const
{
    const unsigned long long key =
        (unsigned long long)(kind != 0) | ((unsigned long long)(unsigned)host << 1);
    {
        std::shared_lock<std::shared_mutex> rl(m_hexFarRuleMutex);
        auto it = m_hexFarRuleCache.find(key);
        if (it != m_hexFarRuleCache.end()) return it->second;
    }
    const std::vector<int>& charges = (kind == 0) ? m_cellCharges[host] : m_faceCharges[host];
    HexFarRule rule;
    const bool cell = (kind == 0);
    const int n1 = (int)m_glOut.size();
    const int np = cell ? n1*n1*n1 : n1*n1;
    const double* nodes = cell ? &m_hexNodes[(size_t)host*81]
                               : &m_quadNodes[(size_t)host*27];
    rule.np = np;
    rule.n_local = (int)charges.size();
    rule.x.resize((size_t)np*3);
    rule.w.resize(np);
    rule.values.resize((size_t)np*charges.size());
    int q = 0;
    for (int iz = 0; iz < (cell ? n1 : 1); ++iz)
        for (int iy = 0; iy < n1; ++iy)
            for (int ix = 0; ix < n1; ++ix, ++q) {
                const double xi[3] = {
                    m_glOut[ix], m_glOut[iy], cell ? m_glOut[iz] : 0.0};
                if (cell) HexQ2MapX(nodes, xi, &rule.x[(size_t)3*q]);
                else {
                    const double uv[2] = {xi[0], xi[1]};
                    QuadQ2MapX(nodes, uv, &rule.x[(size_t)3*q]);
                }
                rule.w[q] = m_gwOut[ix]*m_gwOut[iy]
                          * (cell ? m_gwOut[iz] : 1.0);
                for (int local = 0; local < (int)charges.size(); ++local)
                    rule.values[(size_t)q*charges.size() + local] =
                        HexMonoEval(charges[local], xi);
            }
    std::unique_lock<std::shared_mutex> wl(m_hexFarRuleMutex);
    auto it = m_hexFarRuleCache.emplace(key, std::move(rule)).first;   // racing first insert wins
    return it->second;
}

std::vector<double> RadHACApKChargeGram::QuadBlockHexAffineFarProduct(
    int kindT, int hT, int kindS, int hS, int img) const
{
    const std::vector<int>& tgtG = (kindT == 0) ? m_cellCharges[hT] : m_faceCharges[hT];
    const std::vector<int>& srcG = (kindS == 0) ? m_cellCharges[hS] : m_faceCharges[hS];
    const int nT = (int)tgtG.size(), nS = (int)srcG.size();
    std::vector<double> block((size_t)nT*nS, 0.0);
    if (nT == 0 || nS == 0) return block;

    const HexFarRule& target = GetHexFarRule(kindT, hT);
    const HexFarRule& source = GetHexFarRule(kindS, hS);
    std::vector<double> inner((size_t)nS, 0.0);
    const int nqT = target.np, nqS = source.np;
    for (int qt = 0; qt < nqT; ++qt) {
        double reflected[3];
        ImageEvalPoint(img, &target.x[(size_t)3*qt], reflected);
        std::fill(inner.begin(), inner.end(), 0.0);
        for (int qs = 0; qs < nqS; ++qs) {
            const double dx = reflected[0] - source.x[(size_t)3*qs];
            const double dy = reflected[1] - source.x[(size_t)3*qs + 1];
            const double dz = reflected[2] - source.x[(size_t)3*qs + 2];
            const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (r < 1e-300) continue;
            const double kernel_weight = source.w[qs]/r;
            const double* source_values = &source.values[(size_t)qs*nS];
            for (int ls = 0; ls < nS; ++ls)
                inner[ls] += kernel_weight*source_values[ls];
        }
        const double* target_values = &target.values[(size_t)qt*nT];
        for (int lt = 0; lt < nT; ++lt) {
            const double outer_weight = target.w[qt]*target_values[lt];
            double* row = &block[(size_t)lt*nS];
            for (int ls = 0; ls < nS; ++ls) row[ls] += outer_weight*inner[ls];
        }
    }
    for (double& value : block) value *= RAD_INV_FOUR_PI;
    return block;
}

void RadHACApKChargeGram::PhiInnerHexSubVec(int kindS, int hS, int subB, const double p[3],
                                            const std::vector<int>& srcG, double* inn) const
{
    const bool cell = (kindS == 0);
    const double* nd = cell ? &m_hexNodes[(size_t)hS*81] : &m_quadNodes[(size_t)hS*27];
    if (HEX_USE_AFFINE_EXACT_CELL_INNER && cell
            && hS >= 0 && hS < (int)m_hexAffineCell.size() && m_hexAffineCell[hS]) {
        PhiInnerHexAffineCellSubVec(hS, subB, p, srcG, inn);
        return;
    }
    if (!cell) {
        bool affine = true;
        const double* p0 = &nd[0];
        const double* px = &nd[6];
        const double* py = &nd[18];
        const double tol = 1e-10 * std::max({std::sqrt((px[0]-p0[0])*(px[0]-p0[0])
                                                       + (px[1]-p0[1])*(px[1]-p0[1])
                                                       + (px[2]-p0[2])*(px[2]-p0[2])),
                                             std::sqrt((py[0]-p0[0])*(py[0]-p0[0])
                                                       + (py[1]-p0[1])*(py[1]-p0[1])
                                                       + (py[2]-p0[2])*(py[2]-p0[2]))}) + 1e-12;
        for (int j = 0; affine && j < 3; ++j)
            for (int i = 0; i < 3; ++i) {
                const double xi = 0.5*i, eta = 0.5*j;
                const double* q = &nd[3*(i + 3*j)];
                const double dx = q[0] - (p0[0] + xi*(px[0]-p0[0]) + eta*(py[0]-p0[0]));
                const double dy = q[1] - (p0[1] + xi*(px[1]-p0[1]) + eta*(py[1]-p0[1]));
                const double dz = q[2] - (p0[2] + xi*(px[2]-p0[2]) + eta*(py[2]-p0[2]));
                affine = dx*dx + dy*dy + dz*dz <= tol*tol;
            }
        if (affine
                && hS >= 0 && hS < (int)m_quadAffineFace.size() && m_quadAffineFace[hS]) {
            PhiInnerHexAffineFaceSubVec(hS, subB, p, srcG, inn);
            return;
        }
    }
    const size_t sid = cell ? ((size_t)hS*6 + subB) : ((size_t)hS*2 + subB);
    const double* cs = cell ? &m_cellSubC[sid*3] : &m_faceSubC[sid*3];
    const double  sz = cell ? m_cellSubS[sid] : m_faceSubS[sid];
    const double dxc = p[0]-cs[0], dyc = p[1]-cs[1], dzc = p[2]-cs[2];
    const bool far_pt = std::sqrt(dxc*dxc + dyc*dyc + dzc*dzc) > m_far_inner_factor*sz;
    if (!far_pt) {
        PhiInnerHexSiteVec(kindS, hS, subB, p, srcG, inn);
        return;
    }
    const std::shared_ptr<const HexQuadCloud> cl =
        HexGetCloud(m_build_id, HexCloudKey(cell ? 0 : 1, false, false, hS, subB, 3),
        [&](HexQuadCloud& c) {
            if (cell) HexBuildCloud(nd, true, subB, m_farTetP.data(), m_farTetW.data(),
                                    (int)m_farTetW.size(), false, c);
            else      HexBuildCloud(nd, false, subB, m_farTriP.data(), m_farTriW.data(),
                                    (int)m_farTriW.size(), false, c);
        });
    const int nq = (int)cl->wgeo.size();
    const int nS = (int)srcG.size();
    for (int q = 0; q < nq; ++q) {
        const double dx = p[0]-cl->pts[3*q], dy = p[1]-cl->pts[3*q+1], dz = p[2]-cl->pts[3*q+2];
        const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 1e-300) continue;
        const double gr = cl->wgeo[q]/r;
        const double* xi = &cl->xi[3*q];
        for (int ls = 0; ls < nS; ++ls) inn[ls] += gr*HexMonoEval(srcG[ls], xi);
    }
}

// Affine self/near product rule for the hex charge block.  The legacy path partitions the target cell/face
// along fixed diagonals before applying a barycentric rule; those diagonals are not invariant under mirror
// reflection.  For flat Q2 hosts the geometry map is affine, so integrating over the complete reference
// cube/quad with a symmetric tensor rule removes that artificial choice.  The source potential is summed
// over all reference sub-tets/sub-triangles using the exact polynomial moment kernels above.  Smooth affine
// far blocks use QuadBlockHexAffineFarProduct; curved hosts continue through QuadBlockHex's graded path.
std::vector<double> RadHACApKChargeGram::QuadBlockHexAffineProduct(int kindT, int hT, int kindS, int hS, int img) const
{
    const std::vector<int>& tgtG = (kindT == 0) ? m_cellCharges[hT] : m_faceCharges[hT];
    const std::vector<int>& srcG = (kindS == 0) ? m_cellCharges[hS] : m_faceCharges[hS];
    const int nT = (int)tgtG.size(), nS = (int)srcG.size();
    std::vector<double> blk((size_t)nT*nS, 0.0);
    if (nT == 0 || nS == 0) return blk;

    auto reflpt = [this, img](const double* v, double* o) { ImageEvalPoint(img, v, o); };
    const double* ndT = (kindT == 0) ? &m_hexNodes[(size_t)hT*81] : &m_quadNodes[(size_t)hT*27];
    std::vector<double> inn((size_t)nS, 0.0);
    std::vector<double> xi(3, 0.0);
    const int nq = (int)m_glOut.size();
    for (int iz = 0; iz < (kindT == 0 ? nq : 1); ++iz) {
        xi[2] = (kindT == 0) ? m_glOut[iz] : 0.0;
        for (int iy = 0; iy < nq; ++iy) {
            xi[1] = m_glOut[iy];
            for (int ix = 0; ix < nq; ++ix) {
                xi[0] = m_glOut[ix];
                double p[3];
                if (kindT == 0) HexQ2MapX(ndT, xi.data(), p);
                else {
                    const double uv[2] = {xi[0], xi[1]};
                    QuadQ2MapX(ndT, uv, p);
                }
                double peval[3]; reflpt(p, peval);
                std::fill(inn.begin(), inn.end(), 0.0);
                if (kindS == 0) PhiInnerHexAffineCellVec(hS, peval, srcG, inn.data());
                else            PhiInnerHexAffineFaceVec(hS, peval, srcG, inn.data());
                const double wg = m_gwOut[ix] * m_gwOut[iy] * ((kindT == 0) ? m_gwOut[iz] : 1.0);
                for (int lt = 0; lt < nT; ++lt) {
                    const double wl = wg * HexMonoEval(tgtG[lt], xi.data());
                    double* row = &blk[(size_t)lt*nS];
                    for (int ls = 0; ls < nS; ++ls) row[ls] += wl * inn[ls];
                }
            }
        }
    }
    for (double& v : blk) v *= RAD_INV_FOUR_PI;
    return blk;
}

// 2D closest point on a (ref-space) triangle -- the clamp for the radial anchor on faces.
static void ClosestPointTri2D(const double V[3][2], const double p[2], double out[2])
{
    const double e1u = V[1][0]-V[0][0], e1v = V[1][1]-V[0][1];
    const double e2u = V[2][0]-V[0][0], e2v = V[2][1]-V[0][1];
    const double det = e1u*e2v - e1v*e2u;
    if (std::fabs(det) > 1e-300) {
        const double pu = p[0]-V[0][0], pv = p[1]-V[0][1];
        const double l1 = ( pu*e2v - pv*e2u)/det;
        const double l2 = (-pu*e1v + pv*e1u)/det;
        if (l1 >= 0.0 && l2 >= 0.0 && l1 + l2 <= 1.0) { out[0] = p[0]; out[1] = p[1]; return; }
    }
    double best = 1e300;
    for (int e = 0; e < 3; ++e) {
        const double* A = V[e]; const double* B = V[(e+1)%3];
        const double du = B[0]-A[0], dv = B[1]-A[1];
        const double L2 = du*du + dv*dv;
        double t = (L2 > 1e-300) ? ((p[0]-A[0])*du + (p[1]-A[1])*dv)/L2 : 0.0;
        t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
        const double qu = A[0]+t*du, qv = A[1]+t*dv;
        const double d = (p[0]-qu)*(p[0]-qu) + (p[1]-qv)*(p[1]-qv);
        if (d < best) { best = d; out[0] = qu; out[1] = qv; }
    }
}

// SELF inner: the tet path's PhiAtHO_Duffy RADIAL signed decomposition ported to the REF frame (see the
// header doc).  Anchor x0 = xiT, the outer point's own ref coords (the pulled-back kernel 1/|p-X(xi)|
// peaks there -- exact, no inverse), clamped into the ref sub-simplex; 4 signed radial sub-tets (3 signed
// sub-tris on faces) from x0 with the Duffy apex AT x0: the u^2 (u) volume element kills the 1/r peak
// exactly, and the map's warp enters only as a SMOOTH per-point factor -- robust on strongly distorted
// and curved hexes alike (the corner-graded-cloud / linearized-subtraction schemes oscillated +-3%,
// eig 1.02-1.11 > 1, on the real Cubit cylinder).  SELF-ONLY since 2026-07-03: non-self near calls take
// PhiInnerHexSiteVec (static-site radial; the per-outer-point Newton-anchor branch was removed with
// them).  m_glIn/m_gwIn = the radial 1D Gauss rule.
void RadHACApKChargeGram::PhiInnerHexRadialVec(int kindS, int hS, int subB, const double p[3],
                                               const double* xiT, const std::vector<int>& srcG,
                                               double* inn) const
{
    if (!xiT)
        throw std::logic_error("PhiInnerHexRadialVec: xiT required (SELF-only; non-self near uses the site radial)");
    const bool cell = (kindS == 0);
    if (HEX_USE_AFFINE_EXACT_CELL_INNER && cell
            && hS >= 0 && hS < (int)m_hexAffineCell.size() && m_hexAffineCell[hS]) {
        PhiInnerHexAffineCellSubVec(hS, subB, p, srcG, inn);
        return;
    }
    if (!cell
            && hS >= 0 && hS < (int)m_quadAffineFace.size() && m_quadAffineFace[hS]) {
        PhiInnerHexAffineFaceSubVec(hS, subB, p, srcG, inn);
        return;
    }
    const double* nd = cell ? &m_hexNodes[(size_t)hS*81] : &m_quadNodes[(size_t)hS*27];
    const int nR = (int)m_glIn.size();
    const double* GL = m_glIn.data();
    const double* GW = m_gwIn.data();
    const int nS = (int)srcG.size();
    std::vector<double> acc((size_t)nS, 0.0);

    if (cell) {
        const int* tv = HEXREF_TETS[subB];
        double V[4][3];
        for (int i = 0; i < 4; ++i) for (int k = 0; k < 3; ++k) V[i][k] = HEXREF_V[tv[i]][k];
        // ---- anchor: xiT (the outer point's own ref coords -- the self kernel peaks there), clamped ----
        const double xr[3] = {xiT[0], xiT[1], xiT[2]};
        double x0[3];
        rad_hdiv::ClosestPointTet(V, xr, x0);                    // clamp into the ref sub-tet
        // ---- orientation of the ref sub-tet (computed, not assumed) ----
        double E0[3], E1[3], E2[3];
        for (int k = 0; k < 3; ++k) { E0[k] = V[1][k]-V[0][k]; E1[k] = V[2][k]-V[0][k]; E2[k] = V[3][k]-V[0][k]; }
        const double hv = E0[0]*(E1[1]*E2[2]-E1[2]*E2[1]) - E0[1]*(E1[0]*E2[2]-E1[2]*E2[0])
                        + E0[2]*(E1[0]*E2[1]-E1[1]*E2[0]);
        const double sgnT = (hv >= 0.0) ? 1.0 : -1.0;
        // ---- 4 signed radial sub-tets from x0 (the PhiAtHO_Duffy pattern, in REF space) ----
        for (int f = 0; f < 4; ++f) {
            const double* b1 = V[HEXTET_FC[f][0]]; const double* b2 = V[HEXTET_FC[f][1]];
            const double* b3 = V[HEXTET_FC[f][2]];
            double d1[3], d2[3], d3[3], e21[3], e32[3];
            for (int k = 0; k < 3; ++k) {
                d1[k] = b1[k]-x0[k]; d2[k] = b2[k]-x0[k]; d3[k] = b3[k]-x0[k];
                e21[k] = b2[k]-b1[k]; e32[k] = b3[k]-b2[k];
            }
            const double cr[3] = {d2[1]*d3[2]-d2[2]*d3[1], d2[2]*d3[0]-d2[0]*d3[2], d2[0]*d3[1]-d2[1]*d3[0]};
            const double D = d1[0]*cr[0] + d1[1]*cr[1] + d1[2]*cr[2];   // signed 6*vol(x0,b1,b2,b3), REF
            if (std::fabs(D) < 1e-300) continue;
            for (int a = 0; a < nR; ++a) { const double u = GL[a];
                for (int b = 0; b < nR; ++b) { const double v = GL[b];
                    for (int c = 0; c < nR; ++c) { const double w = GL[c];
                        double y[3];
                        for (int k = 0; k < 3; ++k) y[k] = x0[k] + u*(d1[k] + v*(e21[k] + w*e32[k]));
                        double X[3];
                        HexQ2MapX(nd, y, X);                     // Piola: values-only, no Jacobian
                        const double dx = p[0]-X[0], dy = p[1]-X[1], dz = p[2]-X[2];
                        const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                        if (r < 1e-300) continue;
                        const double wq = GW[a]*GW[b]*GW[c]*(u*u*v*D)/r;   // REF measure (Piola)
                        for (int ls = 0; ls < nS; ++ls) acc[ls] += wq*HexMonoEval(srcG[ls], y);
                    }
                }
            }
        }
        for (int ls = 0; ls < nS; ++ls) inn[ls] += sgnT*acc[ls];
    } else {
        const int* tvq = QUADREF_TRIS[subB];
        double V2[3][2];
        for (int i = 0; i < 3; ++i) for (int k = 0; k < 2; ++k) V2[i][k] = QUADREF_V[tvq[i]][k];
        // ---- anchor: xiT (the outer point's own ref uv coords), clamped ----
        const double ur[2] = {xiT[0], xiT[1]};
        double u0[2];
        ClosestPointTri2D(V2, ur, u0);                           // clamp into the ref sub-tri
        // ---- 3 signed radial sub-tris from u0 (PhiAtHO_Duffy face pattern, in REF uv space) ----
        for (int kf = 0; kf < 3; ++kf) {
            const double* A = V2[kf]; const double* B = V2[(kf+1)%3];
            const double ea[2] = {A[0]-u0[0], A[1]-u0[1]};
            const double eb[2] = {B[0]-u0[0], B[1]-u0[1]};
            const double s2 = ea[0]*eb[1] - ea[1]*eb[0];         // signed 2*area(u0, A, B), REF uv
            if (std::fabs(s2) < 1e-300) continue;
            for (int a = 0; a < nR; ++a) { const double u = GL[a];
                for (int b = 0; b < nR; ++b) { const double v = GL[b];
                    const double yuv[2] = {u0[0] + u*(ea[0] + v*(eb[0]-ea[0])),
                                           u0[1] + u*(ea[1] + v*(eb[1]-ea[1]))};
                    double X[3];
                    QuadQ2MapX(nd, yuv, X);                      // Piola: values-only, no Jacobian
                    const double dx = p[0]-X[0], dy = p[1]-X[1], dz = p[2]-X[2];
                    const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                    if (r < 1e-300) continue;
                    const double wq = GW[a]*GW[b]*(u*s2)/r;               // REF measure (Piola)
                    const double y3[3] = {yuv[0], yuv[1], 0.0};
                    for (int ls = 0; ls < nS; ++ls) acc[ls] += wq*HexMonoEval(srcG[ls], y3);
                }
            }
        }
        // QUADREF_TRIS are CCW (+) in the ref uv frame; the signed s2 pieces sum to the + integral.
        for (int ls = 0; ls < nS; ++ls) inn[ls] += acc[ls];
    }
}

void RadHACApKChargeGram::DPhiInnerHexRadialCellVec(
    int hS, int subB, const double p[3], const double dp[3], const double* xiT,
    const double* node_velocity, const std::vector<int>& srcG, double* dinn,
    const std::vector<double>* radial_points,
    const std::vector<double>* radial_weights) const
{
    if (!xiT) throw std::logic_error("DPhiInnerHexRadialCellVec: xiT required");
    if ((radial_points == nullptr) != (radial_weights == nullptr))
        throw std::invalid_argument("radial points and weights must be provided together");
    const std::vector<double>& gl = radial_points ? *radial_points : m_glIn;
    const std::vector<double>& gw = radial_weights ? *radial_weights : m_gwIn;
    if (gl.empty() || gl.size() != gw.size())
        throw std::invalid_argument("radial points and weights must have equal nonzero length");
    const double* nd = &m_hexNodes[(size_t)hS*81];
    const int* tv = HEXREF_TETS[subB];
    double V[4][3];
    for (int i=0;i<4;++i) for(int k=0;k<3;++k) V[i][k]=HEXREF_V[tv[i]][k];
    double x0[3]; rad_hdiv::ClosestPointTet(V, xiT, x0);
    double E0[3],E1[3],E2[3];
    for(int k=0;k<3;++k){E0[k]=V[1][k]-V[0][k];E1[k]=V[2][k]-V[0][k];E2[k]=V[3][k]-V[0][k];}
    const double hv=E0[0]*(E1[1]*E2[2]-E1[2]*E2[1])-E0[1]*(E1[0]*E2[2]-E1[2]*E2[0])+E0[2]*(E1[0]*E2[1]-E1[1]*E2[0]);
    const double sgnT=hv>=0.0?1.0:-1.0;
    const int nS=(int)srcG.size(), nR=(int)gl.size();
    std::vector<double> acc((size_t)nS,0.0);
    for(int f=0;f<4;++f){
        const double* b1=V[HEXTET_FC[f][0]], *b2=V[HEXTET_FC[f][1]], *b3=V[HEXTET_FC[f][2]];
        double d1[3],d2[3],d3[3],e21[3],e32[3];
        for(int k=0;k<3;++k){d1[k]=b1[k]-x0[k];d2[k]=b2[k]-x0[k];d3[k]=b3[k]-x0[k];e21[k]=b2[k]-b1[k];e32[k]=b3[k]-b2[k];}
        const double cr[3]={d2[1]*d3[2]-d2[2]*d3[1],d2[2]*d3[0]-d2[0]*d3[2],d2[0]*d3[1]-d2[1]*d3[0]};
        const double D=d1[0]*cr[0]+d1[1]*cr[1]+d1[2]*cr[2];
        if(std::fabs(D)<1e-300) continue;
        for(int a=0;a<nR;++a){const double u=gl[a];for(int b=0;b<nR;++b){const double v=gl[b];for(int c=0;c<nR;++c){const double w=gl[c];
            double y[3];for(int k=0;k<3;++k)y[k]=x0[k]+u*(d1[k]+v*(e21[k]+w*e32[k]));
            double X[3],dX[3]; HexQ2MapX(nd,y,X); HexQ2MapX(node_velocity,y,dX);
            const double R[3]={p[0]-X[0],p[1]-X[1],p[2]-X[2]};
            const double dR[3]={dp[0]-dX[0],dp[1]-dX[1],dp[2]-dX[2]};
            const double r2=R[0]*R[0]+R[1]*R[1]+R[2]*R[2]; if(r2<1e-300)continue;
            const double dk=-(R[0]*dR[0]+R[1]*dR[1]+R[2]*dR[2])/(r2*std::sqrt(r2));
            const double wq=gw[a]*gw[b]*gw[c]*(u*u*v*D)*dk;
            for(int ls=0;ls<nS;++ls)acc[ls]+=wq*HexMonoEval(srcG[ls],y);
        }}}
    }
    for(int ls=0;ls<nS;++ls)dinn[ls]+=sgnT*acc[ls];
}

void RadHACApKChargeGram::DPhiInnerHexRadialFaceVec(
    int hS, int subB, const double p[3], const double dp[3], const double* xiT,
    const double* node_velocity, const std::vector<int>& srcG, double* dinn) const
{
    if (!xiT) throw std::logic_error("DPhiInnerHexRadialFaceVec: xiT required");
    const double* nd=&m_quadNodes[(size_t)hS*27];
    const int* tv=QUADREF_TRIS[subB]; double V[3][2];
    for(int i=0;i<3;++i)for(int k=0;k<2;++k)V[i][k]=QUADREF_V[tv[i]][k];
    const double ur[2]={xiT[0],xiT[1]}; double x0[2]; ClosestPointTri2D(V,ur,x0);
    const int n=(int)srcG.size(),nr=(int)m_glIn.size(); std::vector<double> acc((size_t)n,0.0);
    for(int f=0;f<3;++f){
        const double* A=V[f],*B=V[(f+1)%3];
        const double ea[2]={A[0]-x0[0],A[1]-x0[1]},eb[2]={B[0]-x0[0],B[1]-x0[1]};
        const double s2=ea[0]*eb[1]-ea[1]*eb[0]; if(std::fabs(s2)<1e-300)continue;
        for(int a=0;a<nr;++a){const double u=m_glIn[a];for(int b=0;b<nr;++b){const double v=m_glIn[b];
            const double y[2]={x0[0]+u*(ea[0]+v*(eb[0]-ea[0])),x0[1]+u*(ea[1]+v*(eb[1]-ea[1]))};
            double X[3],dX[3]; QuadQ2MapX(nd,y,X); QuadQ2MapX(node_velocity,y,dX);
            const double R[3]={p[0]-X[0],p[1]-X[1],p[2]-X[2]},dR[3]={dp[0]-dX[0],dp[1]-dX[1],dp[2]-dX[2]};
            const double r2=R[0]*R[0]+R[1]*R[1]+R[2]*R[2]; if(r2<1e-300)continue;
            const double dk=-(R[0]*dR[0]+R[1]*dR[1]+R[2]*dR[2])/(r2*std::sqrt(r2));
            const double wq=m_gwIn[a]*m_gwIn[b]*(u*s2)*dk; const double y3[3]={y[0],y[1],0.0};
            for(int j=0;j<n;++j)acc[j]+=wq*HexMonoEval(srcG[j],y3);
        }}
    }
    for(int j=0;j<n;++j)dinn[j]+=acc[j];
}

void RadHACApKChargeGram::DPhiInnerHexSubVec(int kindS,int hS,int subB,
    const double p[3],const double dp[3],const double* velocity,
    const std::vector<int>& srcG,double* dinn) const
{
    const bool cell=kindS==0;const size_t sid=cell?((size_t)hS*6+subB):((size_t)hS*2+subB);
    if(!cell&&hS<(int)m_quadAffineFace.size()&&m_quadAffineFace[hS]){
        const double* nd=&m_quadNodes[(size_t)hS*27];double lin[2][4],dlin[2][4],invj=0,dinvj=0;
        if(!QuadAffineInverseFormsDirectional(nd,velocity,lin,dlin,invj,dinvj))throw std::logic_error("affine HEX face derivative has a singular geometry map");
        const int np=m_quadAffinePolyCount,axis=m_hexAffineOrder+1;const int* tv=QUADREF_TRIS[subB];double V[3][3],dV[3][3];
        for(int a=0;a<3;++a){const double uv[2]={QUADREF_V[tv[a]][0],QUADREF_V[tv[a]][1]};QuadQ2MapX(nd,uv,V[a]);QuadQ2MapX(velocity,uv,dV[a]);}
        double mv[QUAD_AFFINE_POLY_N]={},dm[QUAD_AFFINE_POLY_N]={};if(m_hexAffineOrder==1)rad_hdiv::TriPotentialMomentsDirectionalUpTo2(V,dV,p,dp,mv,dm);else rad_hdiv::TriPotentialMomentsDirectionalUpTo4(V,dV,p,dp,mv,dm);
        for(int ls=0;ls<(int)srcG.size();++ls){const int* e=&m_expo[(size_t)3*srcG[ls]];double poly[HEX_AFFINE_POLY_N]={},dpoly[HEX_AFFINE_POLY_N]={};int deg=0;poly[0]=1;for(int q=0;q<e[0];++q)HexPolyMulLinearDirectional(poly,dpoly,deg,lin[0],dlin[0],np);for(int q=0;q<e[1];++q)HexPolyMulLinearDirectional(poly,dpoly,deg,lin[1],dlin[1],np);const int mono=e[0]+axis*e[1];const double* coeff=&m_quadAffineCoeff[((size_t)hS*m_quadAffineMonoCount+mono)*np];for(int k=0;k<np;++k){const double dc=dinvj*poly[k]+invj*dpoly[k];dinn[ls]+=dc*mv[k]+coeff[k]*dm[k];}}
        return;
    }
    const double* cs=cell?&m_cellSubC[sid*3]:&m_faceSubC[sid*3];const double sz=cell?m_cellSubS[sid]:m_faceSubS[sid];
    const double dx=p[0]-cs[0],dy=p[1]-cs[1],dz=p[2]-cs[2];const bool far=std::sqrt(dx*dx+dy*dy+dz*dz)>m_far_inner_factor*sz;
    const double* nd=cell?&m_hexNodes[(size_t)hS*81]:&m_quadNodes[(size_t)hS*27];const int n=(int)srcG.size();
    if(far){const auto cl=HexGetCloud(m_build_id,HexCloudKey(cell?0:1,false,false,hS,subB,3),[&](HexQuadCloud& c){if(cell)HexBuildCloud(nd,true,subB,m_farTetP.data(),m_farTetW.data(),(int)m_farTetW.size(),false,c);else HexBuildCloud(nd,false,subB,m_farTriP.data(),m_farTriW.data(),(int)m_farTriW.size(),false,c);});
        for(int q=0;q<(int)cl->wgeo.size();++q){const double* xi=&cl->xi[3*q];double dX[3];if(cell)HexQ2MapX(velocity,xi,dX);else{const double uv[2]={xi[0],xi[1]};QuadQ2MapX(velocity,uv,dX);}const double R[3]={p[0]-cl->pts[3*q],p[1]-cl->pts[3*q+1],p[2]-cl->pts[3*q+2]},dR[3]={dp[0]-dX[0],dp[1]-dX[1],dp[2]-dX[2]};const double r2=R[0]*R[0]+R[1]*R[1]+R[2]*R[2];if(r2<1e-300)continue;const double w=-cl->wgeo[q]*(R[0]*dR[0]+R[1]*dR[1]+R[2]*dR[2])/(r2*std::sqrt(r2));for(int j=0;j<n;++j)dinn[j]+=w*HexMonoEval(srcG[j],xi);}return;}
    const int nsite=cell?15:7;const double* sx=cell?&m_cellSiteX[(((size_t)hS*6+subB)*15)*3]:&m_faceSiteX[(((size_t)hS*2+subB)*7)*3];int best=0;double bd=1e300;for(int k=0;k<nsite;++k){const double x=p[0]-sx[3*k],y=p[1]-sx[3*k+1],z=p[2]-sx[3*k+2],d=x*x+y*y+z*z;if(d<bd){bd=d;best=k;}}
    const HexSiteRad& Q=cell?m_cellSiteRad[(size_t)subB*15+best]:m_faceSiteRad[(size_t)subB*7+best];const int nn=cell?27:9,nm=cell?27:9;std::vector<int> col(n);for(int j=0;j<n;++j){const int* e=&m_expo[(size_t)3*srcG[j]];col[j]=e[0]+3*e[1]+(cell?9*e[2]:0);}for(int q=0;q<Q.nq;++q){const double* S=&Q.S[(size_t)q*nn];double X[3]={},dX[3]={};for(int a=0;a<nn;++a)for(int k=0;k<3;++k){X[k]+=S[a]*nd[3*a+k];dX[k]+=S[a]*velocity[3*a+k];}const double R[3]={p[0]-X[0],p[1]-X[1],p[2]-X[2]},dR[3]={dp[0]-dX[0],dp[1]-dX[1],dp[2]-dX[2]};const double r2=R[0]*R[0]+R[1]*R[1]+R[2]*R[2];if(r2<1e-300)continue;const double w=-Q.w[q]*(R[0]*dR[0]+R[1]*dR[1]+R[2]*dR[2])/(r2*std::sqrt(r2));const double* M=&Q.M[(size_t)q*nm];for(int j=0;j<n;++j)dinn[j]+=w*M[col[j]];}
}

std::vector<double> RadHACApKChargeGram::QuadBlockHexDirectionalDerivative(
    int kindT,int hT,int kindS,int hS,const double* velocityT,const double* velocityS,int mask) const
{
    if(mask!=0)throw std::invalid_argument("shape derivative currently requires direct image 0");
    const auto& tg=kindT==0?m_cellCharges[hT]:m_faceCharges[hT];const auto& sg=kindS==0?m_cellCharges[hS]:m_faceCharges[hS];const int nt=(int)tg.size(),ns=(int)sg.size();std::vector<double> out((size_t)nt*ns,0),inn(ns);
    if(kindT==kindS&&hT==hS)return kindT==0?HexVolumeSelfBlockDirectionalDerivative(hT,std::vector<double>(velocityT,velocityT+81)):HexFaceSelfBlockDirectionalDerivative(hT,std::vector<double>(velocityT,velocityT+27));
    if(kindT==1&&kindS==1&&m_quadAffineFace[hT]&&m_quadAffineFace[hS]){
        const double* nd=&m_quadNodes[(size_t)hT*27];const int nq=(int)m_glOut.size();double xi[3]={0,0,0};
        for(int iy=0;iy<nq;++iy){xi[1]=m_glOut[iy];for(int ix=0;ix<nq;++ix){xi[0]=m_glOut[ix];const double uv[2]={xi[0],xi[1]};double p[3],dp[3];QuadQ2MapX(nd,uv,p);QuadQ2MapX(velocityT,uv,dp);std::fill(inn.begin(),inn.end(),0.0);for(int sub=0;sub<2;++sub)DPhiInnerHexSubVec(1,hS,sub,p,dp,velocityS,sg,inn.data());const double wg=m_gwOut[ix]*m_gwOut[iy];for(int i=0;i<nt;++i){const double w=wg*HexMonoEval(tg[i],xi);for(int j=0;j<ns;++j)out[(size_t)i*ns+j]+=w*inn[j];}}}
        for(double& x:out)x*=RAD_INV_FOUR_PI;return out;
    }
    const bool ct=kindT==0,cs=kindS==0;const int nst=ct?6:2,nss=cs?6:2,nv=ct?4:3;const double* ndT=ct?&m_hexNodes[(size_t)hT*81]:&m_quadNodes[(size_t)hT*27];const int rt=tg[0],rs=sg[0];const double dh[3]={m_cent[3*rt]-m_cent[3*rs],m_cent[3*rt+1]-m_cent[3*rs+1],m_cent[3*rt+2]-m_cent[3*rs+2]};const bool nearHosts=std::sqrt(dh[0]*dh[0]+dh[1]*dh[1]+dh[2]*dh[2])<=m_near_grade*(m_size[rt]+m_size[rs]);
    for(int a=0;a<nst;++a){const size_t ia=ct?(size_t)hT*6+a:(size_t)hT*2+a;const double* va=ct?&m_cellSubV[ia*12]:&m_faceSubV[ia*9];const double* ca=ct?&m_cellSubC[ia*3]:&m_faceSubC[ia*3];const double sa=ct?m_cellSubS[ia]:m_faceSubS[ia];for(int b=0;b<nss;++b){const size_t ib=cs?(size_t)hS*6+b:(size_t)hS*2+b;const double* cb=cs?&m_cellSubC[ib*3]:&m_faceSubC[ib*3];const double sb=cs?m_cellSubS[ib]:m_faceSubS[ib];const double dx=ca[0]-cb[0],dy=ca[1]-cb[1],dz=ca[2]-cb[2];const bool near=nearHosts&&std::sqrt(dx*dx+dy*dy+dz*dz)<=m_near_grade*(sa+sb);std::shared_ptr<const HexQuadCloud> oc;
        if(!near)oc=HexGetCloud(m_build_id,HexCloudKey(ct?0:1,true,false,hT,a,3),[&](HexQuadCloud& c){if(ct)HexBuildCloud(ndT,true,a,m_symTetP.data(),m_symTetW.data(),(int)m_symTetW.size(),false,c);else HexBuildCloud(ndT,false,a,m_symTriP.data(),m_symTriW.data(),(int)m_symTriW.size(),false,c);});else{int corner=0;double best=1e300;for(int i=0;i<nv;++i){const double x=va[3*i]-cb[0],y=va[3*i+1]-cb[1],z=va[3*i+2]-cb[2],d=x*x+y*y+z*z;if(d<best){best=d;corner=i;}}oc=HexGetCloud(m_build_id,HexCloudKey(ct?0:1,true,true,hT,a,corner),[&](HexQuadCloud& c){std::vector<double> gb,gw;HexDuffyBary(ct?3:2,corner,m_glOut,m_gwOut,gb,gw);HexBuildCloud(ndT,ct,a,gb.data(),gw.data(),(int)gw.size(),true,c);});}
        for(int q=0;q<(int)oc->wgeo.size();++q){const double* xi=&oc->xi[3*q];double dp[3];if(ct)HexQ2MapX(velocityT,xi,dp);else{const double uv[2]={xi[0],xi[1]};QuadQ2MapX(velocityT,uv,dp);}std::fill(inn.begin(),inn.end(),0);const double p[3]={oc->pts[3*q],oc->pts[3*q+1],oc->pts[3*q+2]};DPhiInnerHexSubVec(kindS,hS,b,p,dp,velocityS,sg,inn.data());for(int i=0;i<nt;++i){const double w=oc->wgeo[q]*HexMonoEval(tg[i],xi);for(int j=0;j<ns;++j)out[(size_t)i*ns+j]+=w*inn[j];}}
    }}for(double& x:out)x*=RAD_INV_FOUR_PI;return out;
}

std::vector<double> RadHACApKChargeGram::HexChargeGramDirectionalDerivative(
    const std::vector<double>& cellVelocity,const std::vector<double>& faceVelocity) const
{
    if(!m_hexmode||m_wedgemode||m_d2)throw std::logic_error("HEX ChargeGram derivative requires a 3D HEX Gram");
    if(cellVelocity.size()!=m_cellCharges.size()*81)throw std::invalid_argument("cell_node_velocity must have shape (ncell,27,3)");
    if(faceVelocity.size()!=m_faceCharges.size()*27)throw std::invalid_argument("face_node_velocity must have shape (nface,9,3)");
    std::vector<double> dense((size_t)m_n*m_n,0.0);const int nh=(int)m_cellCharges.size()+(int)m_faceCharges.size();
    auto kind=[&](int h){return h<(int)m_cellCharges.size()?0:1;};auto host=[&](int h){return kind(h)==0?h:h-(int)m_cellCharges.size();};auto vel=[&](int k,int h){return k==0?&cellVelocity[(size_t)h*81]:&faceVelocity[(size_t)h*27];};auto grp=[&](int k,int h)->const std::vector<int>&{return k==0?m_cellCharges[h]:m_faceCharges[h];};
    for(int A=0;A<nh;++A){const int ka=kind(A),ha=host(A);for(int B=A;B<nh;++B){const int kb=kind(B),hb=host(B);const auto ab=QuadBlockHexDirectionalDerivative(ka,ha,kb,hb,vel(ka,ha),vel(kb,hb));const auto& ga=grp(ka,ha);const auto& gb=grp(kb,hb);if(A==B){for(int i=0;i<(int)ga.size();++i)for(int j=0;j<(int)ga.size();++j)dense[(size_t)ga[i]*m_n+ga[j]]=ab[(size_t)i*ga.size()+j];continue;}const auto ba=QuadBlockHexDirectionalDerivative(kb,hb,ka,ha,vel(kb,hb),vel(ka,ha));for(int i=0;i<(int)ga.size();++i)for(int j=0;j<(int)gb.size();++j){const double x=.5*(ab[(size_t)i*gb.size()+j]+ba[(size_t)j*ga.size()+i]);dense[(size_t)ga[i]*m_n+gb[j]]=x;dense[(size_t)gb[j]*m_n+ga[i]]=x;}}}
    return dense;
}

RadHACApKChargeGramDerivative::RadHACApKChargeGramDerivative(
    const RadHACApKChargeGram& parent,ChargeDerivativeFamily family,
    std::vector<double> cellVelocity,std::vector<double> faceVelocity)
    :m_parent(parent),m_family(family),m_cell_velocity(std::move(cellVelocity)),
     m_face_velocity(std::move(faceVelocity))
{
    if(m_family==ChargeDerivativeFamily::Hex){
        if(m_cell_velocity.size()!=m_parent.m_cellCharges.size()*81||m_face_velocity.size()!=m_parent.m_faceCharges.size()*27)
            throw std::invalid_argument("HEX derivative velocities have wrong host dimensions");
    }else if(m_family==ChargeDerivativeFamily::Tet){
        if(m_cell_velocity.size()!=m_parent.m_hoCellCharges.size()*12||m_face_velocity.size()!=m_parent.m_hoFaceCharges.size()*9)
            throw std::invalid_argument("TET derivative velocities have wrong host dimensions");
    }else if(m_cell_velocity.size()!=m_parent.m_cellCharges.size()*54||m_face_velocity.size()!=m_parent.m_faceCharges.size()*27)
        throw std::invalid_argument("WEDGE derivative velocities have wrong host dimensions");
    static std::atomic<unsigned long long> next{1};m_cache_token=next.fetch_add(1,std::memory_order_relaxed);
}

void RadHACApKChargeGramDerivative::ExtractCoordinates()
{
    m_n_elem=m_parent.m_n;m_ndof=m_parent.m_n;m_coordinates=m_parent.m_cent;
}

double RadHACApKChargeGramDerivative::GetInteractionMatrixElement(int i,int j) const
{
    if(i<0||j<0||i>=m_parent.m_n||j>=m_parent.m_n)throw std::out_of_range("ChargeGram derivative index out of range");
    if(m_family==ChargeDerivativeFamily::Wedge)return m_parent.WedgeChargeGramDirectionalDerivativeElement(i,j,m_cell_velocity,m_face_velocity,m_cache_token);
    if(m_family==ChargeDerivativeFamily::Tet){
        const int nc=(int)m_parent.m_hoCellCharges.size();
        int oa=m_parent.m_kind[i]==0?m_parent.m_host[i]:nc+m_parent.m_host[i];
        int ob=m_parent.m_kind[j]==0?m_parent.m_host[j]:nc+m_parent.m_host[j];
        int li=m_parent.m_hoLocalOf[i],lj=m_parent.m_hoLocalOf[j];
        if(oa>ob){std::swap(oa,ob);std::swap(li,lj);}
        const long long key=((long long)oa<<32)|(unsigned)ob;
        {
            std::lock_guard<std::mutex> lock(m_cache_mutex);
            auto it=m_block_cache.find(key);
            if(it!=m_block_cache.end()){
                const int kb=ob<nc?0:1,hb=ob<nc?ob:ob-nc;
                const int nb=kb==0?(int)m_parent.m_hoCellCharges[hb].size():(int)m_parent.m_hoFaceCharges[hb].size();
                return it->second[(size_t)li*nb+lj];
            }
        }
        auto computed=m_parent.TetChargeGramDirectionalDerivativeImpl(
            m_cell_velocity,m_face_velocity,oa,ob);
        std::lock_guard<std::mutex> lock(m_cache_mutex);
        if(m_block_cache.size()>HexBlockCacheLimit())m_block_cache.clear();
        auto it=m_block_cache.emplace(key,std::move(computed)).first;
        const int kb=ob<nc?0:1,hb=ob<nc?ob:ob-nc;
        const int nb=kb==0?(int)m_parent.m_hoCellCharges[hb].size():(int)m_parent.m_hoFaceCharges[hb].size();
        return it->second[(size_t)li*nb+lj];
    }
    const int ki=m_parent.m_kind[i],kj=m_parent.m_kind[j],hi=m_parent.m_host[i],hj=m_parent.m_host[j];
    int ga=ki==0?hi:(int)m_parent.m_cellCharges.size()+hi,gb=kj==0?hj:(int)m_parent.m_cellCharges.size()+hj;
    const bool transpose=ga>gb;if(transpose){std::swap(ga,gb);}
    const int ka=ga<(int)m_parent.m_cellCharges.size()?0:1,ha=ka==0?ga:ga-(int)m_parent.m_cellCharges.size();
    const int kb=gb<(int)m_parent.m_cellCharges.size()?0:1,hb=kb==0?gb:gb-(int)m_parent.m_cellCharges.size();
    const long long key=((long long)ga<<32)|(unsigned)gb;
    {std::lock_guard<std::mutex> lock(m_cache_mutex);auto it=m_block_cache.find(key);if(it!=m_block_cache.end()){
        const auto& B=kb==0?m_parent.m_cellCharges[hb]:m_parent.m_faceCharges[hb];const int li=m_parent.m_hexLocalOf[i],lj=m_parent.m_hexLocalOf[j];return transpose?it->second[(size_t)lj*B.size()+li]:it->second[(size_t)li*B.size()+lj];}}
    std::vector<double> computed;
    {
        const double* va=ka==0?&m_cell_velocity[(size_t)ha*81]:&m_face_velocity[(size_t)ha*27];
        const double* vb=kb==0?&m_cell_velocity[(size_t)hb*81]:&m_face_velocity[(size_t)hb*27];
        auto ab=m_parent.QuadBlockHexDirectionalDerivative(ka,ha,kb,hb,va,vb);
        if(ga!=gb){auto ba=m_parent.QuadBlockHexDirectionalDerivative(kb,hb,ka,ha,vb,va);const auto& A=ka==0?m_parent.m_cellCharges[ha]:m_parent.m_faceCharges[ha];const auto& B=kb==0?m_parent.m_cellCharges[hb]:m_parent.m_faceCharges[hb];for(int a=0;a<(int)A.size();++a)for(int b=0;b<(int)B.size();++b)ab[(size_t)a*B.size()+b]=.5*(ab[(size_t)a*B.size()+b]+ba[(size_t)b*A.size()+a]);}
        computed=std::move(ab);
    }
    std::lock_guard<std::mutex> lock(m_cache_mutex);auto it=m_block_cache.emplace(key,std::move(computed)).first;
    const auto& A=ka==0?m_parent.m_cellCharges[ha]:m_parent.m_faceCharges[ha];const auto& B=kb==0?m_parent.m_cellCharges[hb]:m_parent.m_faceCharges[hb];
    const int li=m_parent.m_hexLocalOf[i],lj=m_parent.m_hexLocalOf[j];return transpose?it->second[(size_t)lj*B.size()+li]:it->second[(size_t)li*B.size()+lj];
}

void RadHACApKChargeGramDerivative::ClearEntryCache()
{
    std::lock_guard<std::mutex> lock(m_cache_mutex);
    m_block_cache.clear();
    m_block_cache.rehash(0);
}

std::unique_ptr<RadHACApKChargeGramDerivative> RadHACApKChargeGram::BuildDirectionalDerivativeOperator(
    ChargeDerivativeFamily family,const std::vector<double>& cellVelocity,
    const std::vector<double>& faceVelocity,const RadHACApKParams& params) const
{
    auto result=std::make_unique<RadHACApKChargeGramDerivative>(*this,family,cellVelocity,faceVelocity);
    if(!result->BuildHMatrix(params))throw std::runtime_error("failed to build ChargeGram derivative H-matrix");
    // ACA fill may visit many host pairs.  Their exact blocks are only a
    // construction cache; retaining them would silently recreate dense-like
    // storage beside the completed H-matrix.
    result->ClearEntryCache();
    return result;
}

std::vector<double> RadHACApKChargeGram::DirectionalDerivativeContractions(
    ChargeDerivativeFamily family,int nDirections,
    const std::vector<double>& cellVelocity,const std::vector<double>& faceVelocity,
    const std::vector<double>& left,const std::vector<double>& right) const
{
    if(nDirections<1)throw std::invalid_argument("n_directions must be positive");
    if(left.size()!=(size_t)m_n||right.size()!=(size_t)m_n)
        throw std::invalid_argument("left/right vectors must match ChargeGram ndof");
    size_t cellStride=0,faceStride=0,cellNodeStride=0,faceNodeStride=0;
    if(family==ChargeDerivativeFamily::Hex){cellNodeStride=81;faceNodeStride=27;cellStride=m_cellCharges.size()*81;faceStride=m_faceCharges.size()*27;}
    else if(family==ChargeDerivativeFamily::Tet){cellNodeStride=12;faceNodeStride=9;cellStride=m_hoCellCharges.size()*12;faceStride=m_hoFaceCharges.size()*9;}
    else {cellNodeStride=54;faceNodeStride=27;cellStride=m_cellCharges.size()*54;faceStride=m_faceCharges.size()*27;}
    if(cellVelocity.size()!=(size_t)nDirections*cellStride||faceVelocity.size()!=(size_t)nDirections*faceStride)
        throw std::invalid_argument("batched derivative velocity shape mismatch");
    const auto* leaves=static_cast<const st_cHACApK_leafmtxp_t*>(m_leafmtxp);
    const auto* control=static_cast<const st_cHACApK_lcontrol_t*>(m_control);
    const bool useTree=IsValid()&&leaves&&control&&control->lod;
    std::vector<double> result((size_t)nDirections,0.0);
    // Stream one direction at a time.  Exact-entry kernels read through the
    // same parent ChargeGram quadrature/cache state, which is not reentrant.
    for(int kk=0;kk<nDirections;++kk){
        std::vector<double> cv(cellVelocity.begin()+kk*cellStride,cellVelocity.begin()+(kk+1)*cellStride);
        std::vector<double> fv(faceVelocity.begin()+kk*faceStride,faceVelocity.begin()+(kk+1)*faceStride);
        auto hostActive=[](const std::vector<double>& velocity,size_t offset,size_t count){
            for(size_t p=offset;p<offset+count;++p)if(velocity[p]!=0.0)return true;return false;};
        std::vector<unsigned char> active((size_t)m_n,0);
        for(int charge=0;charge<m_n;++charge){const int kind=m_kind[charge],host=m_host[charge];active[charge]=
            kind==0?hostActive(cv,(size_t)host*cellNodeStride,cellNodeStride)
                   :hostActive(fv,(size_t)host*faceNodeStride,faceNodeStride);}
        RadHACApKChargeGramDerivative derivative(*this,family,std::move(cv),std::move(fv));
        long double sum=0.0L;
        if(!useTree){
            for(int i=0;i<m_n;++i){sum+=(long double)left[i]*derivative.GetInteractionMatrixElement(i,i)*right[i];for(int j=i+1;j<m_n;++j){const double a=derivative.GetInteractionMatrixElement(i,j);sum+=(long double)a*((long double)left[i]*right[j]+(long double)left[j]*right[i]);}}
            result[kk]=(double)sum;continue;
        }
        for(int ip=1;ip<=leaves->nlf;++ip){
            const st_cHACApK_leafmtx_t* leaf=leaves->st_lf[ip];
            if(!leaf||leaf->nstrtl>leaf->nstrtt)continue;
            const int nr=leaf->ndl,nc=leaf->ndt,r0=leaf->nstrtl,c0=leaf->nstrtt;
            bool blockActive=false;
            for(int i=0;i<nr&&!blockActive;++i)blockActive=active[(size_t)control->lod[r0+i]-1]!=0;
            for(int j=0;j<nc&&!blockActive;++j)blockActive=active[(size_t)control->lod[c0+j]-1]!=0;
            if(!blockActive)continue;
            const bool upper=r0<c0;
            auto rowDof=[&](int i){return control->lod[r0+i]-1;};
            auto colDof=[&](int j){return control->lod[c0+j]-1;};
            auto addRankOne=[&](const std::vector<double>&u,const std::vector<double>&v){
                long double lu=0,ur=0,lv=0,vr=0;
                for(int i=0;i<nr;++i){const int g=rowDof(i);lu+=(long double)left[g]*u[i];ur+=(long double)u[i]*right[g];}
                for(int j=0;j<nc;++j){const int g=colDof(j);lv+=(long double)left[g]*v[j];vr+=(long double)v[j]*right[g];}
                sum+=lu*vr;if(upper)sum+=lv*ur;
            };
            if(leaf->ltmtx!=1){
                for(int i=0;i<nr;++i){const int gi=rowDof(i);for(int j=0;j<nc;++j){const int gj=colDof(j);const double a=derivative.GetInteractionMatrixElement(gi,gj);sum+=(long double)left[gi]*a*right[gj];if(upper)sum+=(long double)left[gj]*a*right[gi];}}
                continue;
            }
            // ACA on this admissible derivative leaf.  The rank-one factors
            // are contracted immediately and discarded: no dG H-matrix is
            // materialised, while entry values remain the analytic kernel.
            const int rankLimit=std::min({nr,nc,m_derivativeMaxRank});
            std::vector<std::vector<double>> us,vs;us.reserve(rankLimit);vs.reserve(rankLimit);
            std::vector<unsigned char> used((size_t)nr,0);int pivotRow=0;
            long double approximationNorm2=0;
            for(int rank=0;rank<rankLimit;++rank){
                std::vector<double> v((size_t)nc),u((size_t)nr);
                int pivotCol=-1;double pivotAbs=0;
                // A zero residual row need not imply a zero block.  Search
                // unused cluster rows until a stable pivot is found.
                for(int attempt=0;attempt<nr&&pivotCol<0;++attempt){
                    while(pivotRow<nr&&used[pivotRow])++pivotRow;
                    if(pivotRow>=nr){pivotRow=0;while(pivotRow<nr&&used[pivotRow])++pivotRow;}
                    if(pivotRow>=nr)break;
                    for(int j=0;j<nc;++j){double value=derivative.GetInteractionMatrixElement(rowDof(pivotRow),colDof(j));for(size_t s=0;s<us.size();++s)value-=us[s][pivotRow]*vs[s][j];v[j]=value;if(std::abs(value)>pivotAbs){pivotAbs=std::abs(value);pivotCol=j;}}
                    if(pivotAbs<=1e-30){used[pivotRow]=1;pivotCol=-1;pivotAbs=0;++pivotRow;}
                }
                if(pivotCol<0)break;
                const double pivot=v[pivotCol];
                for(int i=0;i<nr;++i){double value=derivative.GetInteractionMatrixElement(rowDof(i),colDof(pivotCol));for(size_t s=0;s<us.size();++s)value-=us[s][i]*vs[s][pivotCol];u[i]=value/pivot;}
                used[pivotRow]=1;
                addRankOne(u,v);
                long double un2=0,vn2=0,cross=0;for(double x:u)un2+=(long double)x*x;for(double x:v)vn2+=(long double)x*x;
                for(size_t s=0;s<us.size();++s){long double uu=0,vv=0;for(int i=0;i<nr;++i)uu+=(long double)us[s][i]*u[i];for(int j=0;j<nc;++j)vv+=(long double)vs[s][j]*v[j];cross+=uu*vv;}
                const long double termNorm2=un2*vn2;approximationNorm2+=termNorm2+2*cross;
                us.push_back(std::move(u));vs.push_back(std::move(v));
                int next=-1;double nextAbs=0;for(int i=0;i<nr;++i)if(!used[i]&&std::abs(us.back()[i])>nextAbs){nextAbs=std::abs(us.back()[i]);next=i;}pivotRow=next<0?0:next;
                if(rank>0&&termNorm2<=m_derivativeAcaEps*m_derivativeAcaEps*std::max((long double)0,approximationNorm2))break;
            }
        }
        result[kk]=(double)sum;
    }
    return result;
}

std::vector<double> RadHACApKChargeGram::DirectionalDerivativeContractionsMany(
    ChargeDerivativeFamily family,int nDirections,int nLeft,
    const std::vector<double>& cellVelocity,const std::vector<double>& faceVelocity,
    const std::vector<double>& left,const std::vector<double>& right) const
{
    if(nDirections<1||nLeft<1)
        throw std::invalid_argument("n_directions and n_left must be positive");
    if(left.size()!=(size_t)nLeft*m_n||right.size()!=(size_t)m_n)
        throw std::invalid_argument("left matrix/right vector must match ChargeGram ndof");
    size_t cellStride=0,faceStride=0,cellNodeStride=0,faceNodeStride=0;
    if(family==ChargeDerivativeFamily::Hex){cellNodeStride=81;faceNodeStride=27;cellStride=m_cellCharges.size()*81;faceStride=m_faceCharges.size()*27;}
    else if(family==ChargeDerivativeFamily::Tet){cellNodeStride=12;faceNodeStride=9;cellStride=m_hoCellCharges.size()*12;faceStride=m_hoFaceCharges.size()*9;}
    else {cellNodeStride=54;faceNodeStride=27;cellStride=m_cellCharges.size()*54;faceStride=m_faceCharges.size()*27;}
    if(cellVelocity.size()!=(size_t)nDirections*cellStride||faceVelocity.size()!=(size_t)nDirections*faceStride)
        throw std::invalid_argument("batched derivative velocity shape mismatch");
    const auto* leaves=static_cast<const st_cHACApK_leafmtxp_t*>(m_leafmtxp);
    const auto* control=static_cast<const st_cHACApK_lcontrol_t*>(m_control);
    const bool useTree=IsValid()&&leaves&&control&&control->lod;
    std::vector<double> result((size_t)nLeft*nDirections,0.0);
    auto hostActive=[](const std::vector<double>& velocity,size_t offset,size_t count){
        for(size_t p=offset;p<offset+count;++p)if(velocity[p]!=0.0)return true;
        return false;
    };
    std::vector<std::unique_ptr<RadHACApKChargeGramDerivative>> derivatives((size_t)nDirections);
    std::vector<std::vector<unsigned char>> activeByDirection((size_t)nDirections);
    for(int kk=0;kk<nDirections;++kk){
        std::vector<double> cv(cellVelocity.begin()+(size_t)kk*cellStride,
                               cellVelocity.begin()+(size_t)(kk+1)*cellStride);
        std::vector<double> fv(faceVelocity.begin()+(size_t)kk*faceStride,
                               faceVelocity.begin()+(size_t)(kk+1)*faceStride);
        auto& active=activeByDirection[(size_t)kk];active.assign((size_t)m_n,0);
        for(int charge=0;charge<m_n;++charge){const int kind=m_kind[charge],host=m_host[charge];active[charge]=
            kind==0?hostActive(cv,(size_t)host*cellNodeStride,cellNodeStride)
                   :hostActive(fv,(size_t)host*faceNodeStride,faceNodeStride);}
        derivatives[(size_t)kk]=std::make_unique<RadHACApKChargeGramDerivative>(
            *this,family,std::move(cv),std::move(fv));
    }
    if(!useTree){
        for(int kk=0;kk<nDirections;++kk){
            auto& derivative=*derivatives[kk];
            std::vector<long double> sums((size_t)nLeft,0.0L);
            for(int i=0;i<m_n;++i){
                const double diagonal=derivative.GetInteractionMatrixElement(i,i);
                for(int ll=0;ll<nLeft;++ll)sums[ll]+=(long double)left[(size_t)ll*m_n+i]*diagonal*right[i];
                for(int j=i+1;j<m_n;++j){
                    const double a=derivative.GetInteractionMatrixElement(i,j);
                    for(int ll=0;ll<nLeft;++ll){const double* l=&left[(size_t)ll*m_n];sums[ll]+=(long double)a*((long double)l[i]*right[j]+(long double)l[j]*right[i]);}
                }
            }
            for(int ll=0;ll<nLeft;++ll)result[(size_t)ll*nDirections+kk]=(double)sums[ll];
        }
        return result;
    }

    // Keep the multi-left contraction batched.  Flat-TET exact-entry kernels
    // read only immutable parent geometry, while every directional derivative
    // owns a mutex-protected host-pair cache.  Their leaf chunks can therefore
    // run concurrently without materialising a derivative H-matrix.  Retain
    // the historical serial traversal for HEX/WEDGE: those kernels still use
    // shared parent caches whose re-entrancy is not part of this contract.
    const int nLeaves=leaves->nlf;
    const int targetTasks=std::max(1,4*radia::GetMaxThreads());
    const int chunksPerDirection=std::min(nLeaves,
        std::max(1,(targetTasks+nDirections-1)/nDirections));
    const int nTasks=nDirections*chunksPerDirection;
    std::vector<long double> partial((size_t)nTasks*nLeft,0.0L);
    auto contractTask=[&](size_t task){
        const int kk=(int)task/chunksPerDirection;
        const int chunk=(int)task%chunksPerDirection;
        const int first=1+(int)((long long)chunk*nLeaves/chunksPerDirection);
        const int last=1+(int)((long long)(chunk+1)*nLeaves/chunksPerDirection);
        auto& derivative=*derivatives[(size_t)kk];
        const auto& active=activeByDirection[(size_t)kk];
        std::vector<long double> sums((size_t)nLeft,0.0L);
        for(int ip=first;ip<last;++ip){
            const st_cHACApK_leafmtx_t* leaf=leaves->st_lf[ip];
            if(!leaf||leaf->nstrtl>leaf->nstrtt)continue;
            const int nr=leaf->ndl,nc=leaf->ndt,r0=leaf->nstrtl,c0=leaf->nstrtt;
            bool blockActive=false;
            for(int i=0;i<nr&&!blockActive;++i)blockActive=active[(size_t)control->lod[r0+i]-1]!=0;
            for(int j=0;j<nc&&!blockActive;++j)blockActive=active[(size_t)control->lod[c0+j]-1]!=0;
            if(!blockActive)continue;
            const bool upper=r0<c0;
            auto rowDof=[&](int i){return control->lod[r0+i]-1;};
            auto colDof=[&](int j){return control->lod[c0+j]-1;};
            auto addRankOne=[&](const std::vector<double>&u,const std::vector<double>&v){
                long double ur=0,vr=0;
                for(int i=0;i<nr;++i)ur+=(long double)u[i]*right[rowDof(i)];
                for(int j=0;j<nc;++j)vr+=(long double)v[j]*right[colDof(j)];
                for(int ll=0;ll<nLeft;++ll){
                    const double* l=&left[(size_t)ll*m_n];long double lu=0,lv=0;
                    for(int i=0;i<nr;++i)lu+=(long double)l[rowDof(i)]*u[i];
                    for(int j=0;j<nc;++j)lv+=(long double)l[colDof(j)]*v[j];
                    sums[ll]+=lu*vr;if(upper)sums[ll]+=lv*ur;
                }
            };
            if(leaf->ltmtx!=1){
                for(int i=0;i<nr;++i){const int gi=rowDof(i);for(int j=0;j<nc;++j){const int gj=colDof(j);const double a=derivative.GetInteractionMatrixElement(gi,gj);for(int ll=0;ll<nLeft;++ll){const double* l=&left[(size_t)ll*m_n];sums[ll]+=(long double)l[gi]*a*right[gj];if(upper)sums[ll]+=(long double)l[gj]*a*right[gi];}}}
                continue;
            }
            const int rankLimit=std::min({nr,nc,m_derivativeMaxRank});
            std::vector<std::vector<double>> us,vs;us.reserve(rankLimit);vs.reserve(rankLimit);
            std::vector<unsigned char> used((size_t)nr,0);int pivotRow=0;
            long double approximationNorm2=0;
            for(int rank=0;rank<rankLimit;++rank){
                std::vector<double> v((size_t)nc),u((size_t)nr);
                int pivotCol=-1;double pivotAbs=0;
                for(int attempt=0;attempt<nr&&pivotCol<0;++attempt){
                    while(pivotRow<nr&&used[pivotRow])++pivotRow;
                    if(pivotRow>=nr){pivotRow=0;while(pivotRow<nr&&used[pivotRow])++pivotRow;}
                    if(pivotRow>=nr)break;
                    for(int j=0;j<nc;++j){double value=derivative.GetInteractionMatrixElement(rowDof(pivotRow),colDof(j));for(size_t s=0;s<us.size();++s)value-=us[s][pivotRow]*vs[s][j];v[j]=value;if(std::abs(value)>pivotAbs){pivotAbs=std::abs(value);pivotCol=j;}}
                    if(pivotAbs<=1e-30){used[pivotRow]=1;pivotCol=-1;pivotAbs=0;++pivotRow;}
                }
                if(pivotCol<0)break;
                const double pivot=v[pivotCol];
                for(int i=0;i<nr;++i){double value=derivative.GetInteractionMatrixElement(rowDof(i),colDof(pivotCol));for(size_t s=0;s<us.size();++s)value-=us[s][i]*vs[s][pivotCol];u[i]=value/pivot;}
                used[pivotRow]=1;addRankOne(u,v);
                long double un2=0,vn2=0,cross=0;for(double x:u)un2+=(long double)x*x;for(double x:v)vn2+=(long double)x*x;
                for(size_t s=0;s<us.size();++s){long double uu=0,vv=0;for(int i=0;i<nr;++i)uu+=(long double)us[s][i]*u[i];for(int j=0;j<nc;++j)vv+=(long double)vs[s][j]*v[j];cross+=uu*vv;}
                const long double termNorm2=un2*vn2;approximationNorm2+=termNorm2+2*cross;
                us.push_back(std::move(u));vs.push_back(std::move(v));
                int next=-1;double nextAbs=0;for(int i=0;i<nr;++i)if(!used[i]&&std::abs(us.back()[i])>nextAbs){nextAbs=std::abs(us.back()[i]);next=i;}pivotRow=next<0?0:next;
                if(rank>0&&termNorm2<=m_derivativeAcaEps*m_derivativeAcaEps*std::max((long double)0,approximationNorm2))break;
            }
        }
        for(int ll=0;ll<nLeft;++ll)partial[task*(size_t)nLeft+ll]=sums[ll];
    };
    if(family==ChargeDerivativeFamily::Tet&&nTasks>1){
        ngcore::RegionTaskManager task_manager;
        ngcore::ParallelFor(ngcore::IntRange(nTasks),[&](int task){
            contractTask((size_t)task);
        });
    }else{
        for(size_t task=0;task<(size_t)nTasks;++task)contractTask(task);
    }
    for(int kk=0;kk<nDirections;++kk)for(int chunk=0;chunk<chunksPerDirection;++chunk){
        const size_t task=(size_t)kk*chunksPerDirection+chunk;
        for(int ll=0;ll<nLeft;++ll)
            result[(size_t)ll*nDirections+kk]+=(double)partial[task*(size_t)nLeft+ll];
    }
    return result;
}

// The whole DIRECTED host-pair block (target host (kindT,hT) outer x source host (kindS,hS) inner) for every
// local charge pair, computed in ONE pass.  All near/far/grading decisions are host+sub geometric (identical
// across the block), so the per-entry value is served from ONE block computation -- the expensive kernel
// work on each (outer pt, inner pt) is shared across all nT*nS monomial combos.  Returns [nT*nS]
// row-major, INV4PI folded.
std::vector<double> RadHACApKChargeGram::QuadBlockHex(int kindT, int hT, int kindS, int hS, int img) const
{
    const std::vector<int>& tgtG = (kindT == 0) ? m_cellCharges[hT] : m_faceCharges[hT];
    const std::vector<int>& srcG = (kindS == 0) ? m_cellCharges[hS] : m_faceCharges[hS];
    const int nT = (int)tgtG.size(), nS = (int)srcG.size();
    std::vector<double> blk((size_t)nT*nS, 0.0);
    if (nT == 0 || nS == 0) return blk;
    auto face_affine = [this](int h) {
        const double* nd = &m_quadNodes[(size_t)h*27];
        const double* p0 = &nd[0];
        const double* px = &nd[6];
        const double* py = &nd[18];
        const double tol = 1e-10 * std::max({
            std::sqrt((px[0]-p0[0])*(px[0]-p0[0]) + (px[1]-p0[1])*(px[1]-p0[1]) + (px[2]-p0[2])*(px[2]-p0[2])),
            std::sqrt((py[0]-p0[0])*(py[0]-p0[0]) + (py[1]-p0[1])*(py[1]-p0[1]) + (py[2]-p0[2])*(py[2]-p0[2]))}) + 1e-12;
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i) {
                const double u = 0.5*i, v = 0.5*j;
                const double* q = &nd[3*(i + 3*j)];
                const double dx = q[0] - (p0[0] + u*(px[0]-p0[0]) + v*(py[0]-p0[0]));
                const double dy = q[1] - (p0[1] + u*(px[1]-p0[1]) + v*(py[1]-p0[1]));
                const double dz = q[2] - (p0[2] + u*(px[2]-p0[2]) + v*(py[2]-p0[2]));
                if (dx*dx + dy*dy + dz*dz > tol*tol) return false;
            }
        return true;
    };
    const bool affineT = kindT == 0 && hT >= 0 && hT < (int)m_hexAffineCell.size() && m_hexAffineCell[hT]
                       || kindT != 0 && face_affine(hT);
    const bool affineS = kindS == 0 && hS >= 0 && hS < (int)m_hexAffineCell.size() && m_hexAffineCell[hS]
                       || kindS != 0 && face_affine(hS);
    // Host-pair separation on the image-mapped source centroid: drives the affine near/far product
    // switch, the distorted-pair far switch, and the general path's near_hosts grading below.
    const int repA = tgtG[0], repB = srcG[0];
    double repBc[3];
    ImageEvalPoint(img, &m_cent[(size_t)3*repB], repBc);
    const double sep = std::sqrt(
        (m_cent[(size_t)3*repA]     - repBc[0])*(m_cent[(size_t)3*repA]     - repBc[0])
      + (m_cent[(size_t)3*repA + 1] - repBc[1])*(m_cent[(size_t)3*repA + 1] - repBc[1])
      + (m_cent[(size_t)3*repA + 2] - repBc[2])*(m_cent[(size_t)3*repA + 2] - repBc[2]));
    auto timed = [](std::atomic<long long>& blocks, std::atomic<long long>& ns_acc,
                    auto&& fn) -> std::vector<double> {
        const auto t0 = std::chrono::steady_clock::now();
        std::vector<double> out = fn();
        blocks.fetch_add(1, std::memory_order_relaxed);
        ns_acc.fetch_add((long long)std::chrono::duration_cast<std::chrono::nanoseconds>(
                             std::chrono::steady_clock::now() - t0).count(),
                         std::memory_order_relaxed);
        return out;
    };
    if (affineT && affineS) {
        const bool exact_near = sep <= HEX_AFFINE_EXACT_NEAR_FACTOR*(m_size[repA] + m_size[repB]);
        return exact_near
            ? timed(m_hexBlkAffineNear, m_hexNsAffineNear,
                    [&]{ return QuadBlockHexAffineProduct(kindT, hT, kindS, hS, img); })
            : timed(m_hexBlkAffineFar, m_hexNsAffineFar,
                    [&]{ return QuadBlockHexAffineFarProduct(kindT, hT, kindS, hS, img); });
    }
    // DISTORTED-pair far switch (see HexDistortedFarFactor): the tensor far product is geometry-map
    // exact (Q2 point placement, Piola reference charge measure), so a well-separated pair with a
    // distorted/curved host needs no 6x6-sub graded machinery either.  Self pairs have sep == 0 and
    // always stay on the graded/radial path.
    {
        const double fac = HexDistortedFarFactor();
        if (fac > 0.0 && sep > fac*(m_size[repA] + m_size[repB]))
            return timed(m_hexBlkDistortedFar, m_hexNsDistortedFar,
                         [&]{ return QuadBlockHexAffineFarProduct(kindT, hT, kindS, hS, img); });
    }
    const auto t_general0 = std::chrono::steady_clock::now();
    const bool cellT = (kindT == 0), cellS = (kindS == 0);
    const int nsubT = cellT ? 6 : 2, nsubS = cellS ? 6 : 2;
    // IMA (img>0): couple the TARGET host with the source host mapped by the image transform T.  By the
    // isometry |T(x)-y| = |x-T^-1(y)|, the image block = INT_target m_t(x) Phi_source(T^-1(x)) dx:
    // map the source-side geometry driving the near/far grading, map the outer eval point before the
    // source-potential eval, and never treat the pair as SELF.  img==0 => reflpt identity => byte-identical.
    auto reflpt = [this, img](const double* v, double* o){ ImageEvalPoint(img, v, o); };
    const int rt = tgtG[0], rs = srcG[0];      // representative charges (host-level cent/size)
    const bool near_hosts = (img == 0 && kindT == kindS && hT == hS)
                            || sep <= m_near_grade*(m_size[rt] + m_size[rs]);   // sep == reflected r_h above
    // SELF host pair: the inner takes the RADIAL decomposition with the EXACT anchor xiT (the outer
    // point's own ref coords -- no Newton).  The OUTER grading below is UNCHANGED -- it is required by the
    // Q1 charge degree regardless of how the inner is computed (exact inner + regular outer -> eig 1.088).
    // self_pair: the (reflected) source host coincides with the target host -> use the EXACT self-radial.
    // img==0 is always self for hT==hS.  img>0 is self ONLY when the host is INVARIANT under the transform
    // (lies ON the mirror plane, T(host)==host -- e.g. a z=0 cut FACE reflected across z=0 is itself; under a
    // cyclic ROTATION only an on-axis host qualifies).  Then
    // the reflected self-term EXACTLY equals the direct self-term (same exact radial quadrature), so the LARGE
    // on-plane cut-face charge (sigma = M.n ~ |M| when M is perpendicular to the plane) CANCELS exactly for
    // sign -1, instead of leaving a quadrature-mismatch residual (the ~1.5% hex / ~29% wedge antisymmetric-
    // plane error, 2026-07-05: direct self used the exact radial, reflected self used the ~1e-3 site-radial).
    // An OFF-plane self host (a z>0 cell) reflects to a genuine image elsewhere -> NOT self.
    bool self_pair = (kindT == kindS && hT == hS);
    if (self_pair && img != 0) {
        double rc_[3]; reflpt(&m_cent[3*rt], rc_);
        const double d_ = std::abs(rc_[0]-m_cent[3*rt]) + std::abs(rc_[1]-m_cent[3*rt+1]) + std::abs(rc_[2]-m_cent[3*rt+2]);
        self_pair = (d_ < 1e-6 * m_size[rt] + 1e-12);   // T(host)==host <=> host is invariant
    }
    const int nqreg = cellT ? (int)m_symTetW.size() : (int)m_symTriW.size();
    const double* ndT = cellT ? &m_hexNodes[(size_t)hT*81] : &m_quadNodes[(size_t)hT*27];
    const int nvT = cellT ? 4 : 3;
    std::vector<double> inn(nS), owt(nT);
    for (int sA = 0; sA < nsubT; ++sA) {
        const size_t sidA = cellT ? ((size_t)hT*6 + sA) : ((size_t)hT*2 + sA);
        const double szA = cellT ? m_cellSubS[sidA] : m_faceSubS[sidA];
        const double* subVA = cellT ? &m_cellSubV[sidA*4*3] : &m_faceSubV[sidA*3*3];
        for (int sB = 0; sB < nsubS; ++sB) {
            const size_t sidB = cellS ? ((size_t)hS*6 + sB) : ((size_t)hS*2 + sB);
            const double* cB0 = cellS ? &m_cellSubC[sidB*3] : &m_faceSubC[sidB*3];
            double cB[3]; reflpt(cB0, cB);   // mapped source sub-centroid drives near/far + Duffy corner (img>0)
            const double szB = cellS ? m_cellSubS[sidB] : m_faceSubS[sidB];
            const double* cA = cellT ? &m_cellSubC[sidA*3] : &m_faceSubC[sidA*3];
            const double dx = cA[0]-cB[0], dy = cA[1]-cB[1], dz = cA[2]-cB[2];
            const bool near_sub = near_hosts &&
                std::sqrt(dx*dx + dy*dy + dz*dz) <= m_near_grade*(szA + szB);
            // OUTER geometry cloud on target sub sA (monomial-FREE): regular symmetric or graded toward
            // cB.  HELD as a shared_ptr: the inner calls below fetch far clouds from the same cache, and
            // its capacity clear must not invalidate this hold (the n=10 0xC0000005 use-after-free).
            std::shared_ptr<const HexQuadCloud> oc;
            if (!near_sub) {
                oc = HexGetCloud(m_build_id, HexCloudKey(cellT ? 0 : 1, true, false, hT, sA, 3),
                    [&](HexQuadCloud& c) {
                        if (cellT) HexBuildCloud(ndT, true, sA, m_symTetP.data(), m_symTetW.data(), nqreg, false, c);
                        else       HexBuildCloud(ndT, false, sA, m_symTriP.data(), m_symTriW.data(), nqreg, false, c);
                    });
            } else {
                int corner = 0; double best = 1e300;
                for (int i = 0; i < nvT; ++i) {
                    const double ddx = subVA[3*i]-cB[0], ddy = subVA[3*i+1]-cB[1], ddz = subVA[3*i+2]-cB[2];
                    const double d = ddx*ddx + ddy*ddy + ddz*ddz;
                    if (d < best) { best = d; corner = i; }
                }
                oc = HexGetCloud(m_build_id, HexCloudKey(cellT ? 0 : 1, true, true, hT, sA, corner),
                    [&](HexQuadCloud& c) {
                        std::vector<double> gb, gw;
                        HexDuffyBary(cellT ? 3 : 2, corner, m_glOut, m_gwOut, gb, gw);
                        HexBuildCloud(ndT, cellT, sA, gb.data(), gw.data(), (int)gw.size(), true, c);
                    });
            }
            const int nqo = (int)oc->wgeo.size();
            for (int q = 0; q < nqo; ++q) {
                const double pq[3] = {oc->pts[3*q], oc->pts[3*q+1], oc->pts[3*q+2]};
                const double* xiT = &oc->xi[3*q];
                for (int ls = 0; ls < nS; ++ls) inn[ls] = 0.0;
                double peval[3]; reflpt(pq, peval);                                             // inverse-map for the (image) source eval
                if (self_pair) PhiInnerHexRadialVec(kindS, hS, sB, pq, xiT, srcG, inn.data());  // radial, exact anchor
                else           PhiInnerHexSubVec(kindS, hS, sB, peval, srcG, inn.data());       // far cloud / radial (peval==pq if img==0)

                const double wg = oc->wgeo[q];
                for (int lt = 0; lt < nT; ++lt) owt[lt] = wg*HexMonoEval(tgtG[lt], xiT);
                for (int lt = 0; lt < nT; ++lt) {
                    const double wl = owt[lt];
                    double* row = &blk[(size_t)lt*nS];
                    for (int ls = 0; ls < nS; ++ls) row[ls] += wl*inn[ls];
                }
            }
        }
    }
    for (double& v : blk) v *= RAD_INV_FOUR_PI;
    (near_hosts ? m_hexBlkGeneralNear : m_hexBlkGeneralFar)
        .fetch_add(1, std::memory_order_relaxed);
    (near_hosts ? m_hexNsGeneralNear : m_hexNsGeneralFar)
        .fetch_add((long long)std::chrono::duration_cast<std::chrono::nanoseconds>(
                       std::chrono::steady_clock::now() - t_general0).count(),
                   std::memory_order_relaxed);
    return blk;
}

std::vector<double> RadHACApKChargeGram::HexVolumeSelfBlockDirectionalDerivative(
    int host, const std::vector<double>& node_velocity) const
{
    if (!m_hexmode || m_wedgemode) throw std::logic_error(
        "HexVolumeSelfBlockDirectionalDerivative requires a 3D HEX charge Gram");
    if (host < 0 || host >= (int)m_cellCharges.size()) throw std::out_of_range("HEX host out of range");
    if (node_velocity.size() != 81) throw std::invalid_argument("node_velocity must have shape (27,3)");
    const std::vector<int>& g=m_cellCharges[host]; const int n=(int)g.size();
    std::vector<double> out((size_t)n*n,0.0), dinn(n), owt(n);
    if(n==0)return out;
    const double* nd=&m_hexNodes[(size_t)host*81];
    const bool affine = host < (int)m_hexAffineCell.size() && m_hexAffineCell[host];
    if (affine) {
        // Match QuadBlockHexAffineProduct's complete-cube outer rule.  Pure
        // dilation follows exactly from the degree -1 Laplace homogeneity;
        // other velocities differentiate the radial source integral directly.
        double origin[3], geometry_gradient[3][3], geometry_det;
        double inverse_gradient[3][3];
        if (!HexAffineBasisChecked(nd, origin, geometry_gradient, geometry_det)
                || !HexInv3(geometry_gradient, inverse_gradient))
            throw std::logic_error("affine HEX derivative has a singular geometry map");
        double velocity_gradient[3][3];
        const int axis_node[3] = {2, 6, 18};
        for (int physical = 0; physical < 3; ++physical)
            for (int reference = 0; reference < 3; ++reference)
                velocity_gradient[physical][reference] =
                    node_velocity[3*axis_node[reference] + physical]
                    - node_velocity[physical];
        double trace = 0.0;
        for (int physical = 0; physical < 3; ++physical)
            for (int reference = 0; reference < 3; ++reference)
                trace += velocity_gradient[physical][reference]
                       * inverse_gradient[reference][physical];
        const double isotropic_rate = trace/3.0;
        std::vector<double> residual_velocity(81);
        for (int i = 0; i < 81; ++i)
            residual_velocity[i] = node_velocity[i] - isotropic_rate*nd[i];
        double velocity_scale = 1.0;
        for (double value : node_velocity)
            velocity_scale = std::max(velocity_scale, std::abs(value));
        const double rigid_tolerance = 64.0*std::numeric_limits<double>::epsilon()
                                     * velocity_scale;
        bool isotropic_plus_translation = true;
        for (int node = 1; node < 27 && isotropic_plus_translation; ++node)
            for (int physical = 0; physical < 3; ++physical)
                if (std::abs(residual_velocity[3*node + physical]
                           - residual_velocity[physical]) > rigid_tolerance) {
                    isotropic_plus_translation = false;
                    break;
                }
        if (isotropic_plus_translation) {
            const auto exact_block = QuadBlockHexAffineProduct(0, host, 0, host, 0);
            for (std::size_t i = 0; i < out.size(); ++i)
                out[i] = -isotropic_rate*exact_block[i];
            return out;
        }

        struct GaussRule {
            std::vector<double> points;
            std::vector<double> weights;
        };
        static const GaussRule radial = [] {
            // This derivative is an offline shape operation.  A fixed rule
            // keeps its accuracy independent of the solve-time inner rule.
            constexpr int order = 24;
            GaussRule rule;
            rule.points.resize(order);
            rule.weights.resize(order);
            const double pi = std::acos(-1.0);
            for (int i = 0; i < (order + 1)/2; ++i) {
                double z = std::cos(pi*(i + 0.75)/(order + 0.5));
                double derivative = 0.0;
                for (;;) {
                    double p1 = 1.0, p2 = 0.0;
                    for (int degree = 1; degree <= order; ++degree) {
                        const double p3 = p2;
                        p2 = p1;
                        p1 = ((2.0*degree - 1.0)*z*p2 - (degree - 1.0)*p3)/degree;
                    }
                    derivative = order*(z*p1 - p2)/(z*z - 1.0);
                    const double next = z - p1/derivative;
                    if (std::abs(next - z) <= 4.0*std::numeric_limits<double>::epsilon()) {
                        z = next;
                        break;
                    }
                    z = next;
                }
                const double weight = 1.0/((1.0 - z*z)*derivative*derivative);
                rule.points[i] = 0.5*(1.0 - z);
                rule.points[order - 1 - i] = 0.5*(1.0 + z);
                rule.weights[i] = weight;
                rule.weights[order - 1 - i] = weight;
            }
            return rule;
        }();
        const int nq = (int)m_glOut.size();
        double xi[3];
        for (int iz = 0; iz < nq; ++iz) {
            xi[2] = m_glOut[iz];
            for (int iy = 0; iy < nq; ++iy) {
                xi[1] = m_glOut[iy];
                for (int ix = 0; ix < nq; ++ix) {
                    xi[0] = m_glOut[ix];
                    double p[3], dp[3];
                    HexQ2MapX(nd, xi, p);
                    HexQ2MapX(node_velocity.data(), xi, dp);
                    std::fill(dinn.begin(), dinn.end(), 0.0);
                    for (int subB = 0; subB < 6; ++subB)
                        DPhiInnerHexRadialCellVec(
                            host, subB, p, dp, xi, node_velocity.data(),
                            g, dinn.data(), &radial.points, &radial.weights);
                    const double wg = m_gwOut[ix]*m_gwOut[iy]*m_gwOut[iz];
                    for (int lt = 0; lt < n; ++lt)
                        owt[lt] = wg*HexMonoEval(g[lt], xi);
                    for (int lt = 0; lt < n; ++lt)
                        for (int ls = 0; ls < n; ++ls)
                            out[(size_t)lt*n+ls] += owt[lt]*dinn[ls];
                }
            }
        }
    } else {
    for(int sA=0;sA<6;++sA){
        const size_t sidA=(size_t)host*6+sA; const double* subVA=&m_cellSubV[sidA*12]; const double* cA=&m_cellSubC[sidA*3];
        for(int sB=0;sB<6;++sB){
            const size_t sidB=(size_t)host*6+sB; const double* cB=&m_cellSubC[sidB*3];
            const double dcx=cA[0]-cB[0],dcy=cA[1]-cB[1],dcz=cA[2]-cB[2];
            const bool near_sub=std::sqrt(dcx*dcx+dcy*dcy+dcz*dcz)<=m_near_grade*(m_cellSubS[sidA]+m_cellSubS[sidB]);
            int corner=0; double best=1e300;
            for(int i=0;i<4;++i){const double dx=subVA[3*i]-cB[0],dy=subVA[3*i+1]-cB[1],dz=subVA[3*i+2]-cB[2];const double d=dx*dx+dy*dy+dz*dz;if(d<best){best=d;corner=i;}}
            std::shared_ptr<const HexQuadCloud> oc;
            if(near_sub) oc=HexGetCloud(m_build_id,HexCloudKey(0,true,true,host,sA,corner),[&](HexQuadCloud& c){std::vector<double> gb,gw;HexDuffyBary(3,corner,m_glOut,m_gwOut,gb,gw);HexBuildCloud(nd,true,sA,gb.data(),gw.data(),(int)gw.size(),true,c);});
            else oc=HexGetCloud(m_build_id,HexCloudKey(0,true,false,host,sA,3),[&](HexQuadCloud& c){HexBuildCloud(nd,true,sA,m_symTetP.data(),m_symTetW.data(),(int)m_symTetW.size(),false,c);});
            for(int q=0;q<(int)oc->wgeo.size();++q){
                const double* xi=&oc->xi[3*q]; const double p[3]={oc->pts[3*q],oc->pts[3*q+1],oc->pts[3*q+2]};
                double dp[3]; HexQ2MapX(node_velocity.data(),xi,dp);
                std::fill(dinn.begin(),dinn.end(),0.0);
                DPhiInnerHexRadialCellVec(host,sB,p,dp,xi,node_velocity.data(),g,dinn.data());
                for(int lt=0;lt<n;++lt)owt[lt]=oc->wgeo[q]*HexMonoEval(g[lt],xi);
                for(int lt=0;lt<n;++lt)for(int ls=0;ls<n;++ls)out[(size_t)lt*n+ls]+=owt[lt]*dinn[ls];
            }
        }
    }
    }
    for(int i=0;i<n;++i) for(int j=i+1;j<n;++j) {
        const double value=0.5*(out[(size_t)i*n+j]+out[(size_t)j*n+i]);
        out[(size_t)i*n+j]=value; out[(size_t)j*n+i]=value;
    }
    for(double& v:out)v*=RAD_INV_FOUR_PI;
    return out;
}

std::vector<double> RadHACApKChargeGram::HexFaceSelfBlockDirectionalDerivative(
    int host, const std::vector<double>& node_velocity) const
{
    if (!m_hexmode || m_wedgemode) throw std::logic_error(
        "HexFaceSelfBlockDirectionalDerivative requires a 3D HEX charge Gram");
    if(host<0||host>=(int)m_faceCharges.size())throw std::out_of_range("HEX face host out of range");
    if(node_velocity.size()!=27)throw std::invalid_argument("node_velocity must have shape (9,3)");
    const std::vector<int>& g=m_faceCharges[host]; const int n=(int)g.size();
    std::vector<double> out((size_t)n*n,0.0),dinn(n); if(n==0)return out;
    const double* nd=&m_quadNodes[(size_t)host*27];
    const bool affine=host<(int)m_quadAffineFace.size()&&m_quadAffineFace[host];
    if(affine){
        double lin[2][4],dlin[2][4],invj=0,dinvj=0;
        if(!QuadAffineInverseFormsDirectional(nd,node_velocity.data(),lin,dlin,invj,dinvj))
            throw std::logic_error("affine HEX face derivative has a singular geometry map");
        const int axisCount=m_hexAffineOrder+1,np=m_quadAffinePolyCount;
        std::vector<double> dcoeff((size_t)n*np,0.0);
        for(int ls=0;ls<n;++ls){const int* e=&m_expo[(size_t)3*g[ls]];double poly[HEX_AFFINE_POLY_N]={},dp[HEX_AFFINE_POLY_N]={};int deg=0;poly[0]=1;
            for(int q=0;q<e[0];++q)HexPolyMulLinearDirectional(poly,dp,deg,lin[0],dlin[0],np);
            for(int q=0;q<e[1];++q)HexPolyMulLinearDirectional(poly,dp,deg,lin[1],dlin[1],np);
            for(int k=0;k<np;++k)dcoeff[(size_t)ls*np+k]=dinvj*poly[k]+invj*dp[k];
        }
        const int nq=(int)m_glOut.size();double xi[3]={0,0,0};
        for(int iy=0;iy<nq;++iy){xi[1]=m_glOut[iy];for(int ix=0;ix<nq;++ix){xi[0]=m_glOut[ix];
            const double uv[2]={xi[0],xi[1]};double p[3],dpnt[3];QuadQ2MapX(nd,uv,p);QuadQ2MapX(node_velocity.data(),uv,dpnt);
            std::fill(dinn.begin(),dinn.end(),0.0);
            for(int sub=0;sub<2;++sub){const int* tv=QUADREF_TRIS[sub];double V[3][3],dV[3][3];
                for(int a=0;a<3;++a){const double ruv[2]={QUADREF_V[tv[a]][0],QUADREF_V[tv[a]][1]};QuadQ2MapX(nd,ruv,V[a]);QuadQ2MapX(node_velocity.data(),ruv,dV[a]);}
                double mv[QUAD_AFFINE_POLY_N]={},dm[QUAD_AFFINE_POLY_N]={};
                if(m_hexAffineOrder==1)rad_hdiv::TriPotentialMomentsDirectionalUpTo2(V,dV,p,dpnt,mv,dm);
                else rad_hdiv::TriPotentialMomentsDirectionalUpTo4(V,dV,p,dpnt,mv,dm);
                for(int ls=0;ls<n;++ls){const int* e=&m_expo[(size_t)3*g[ls]];const int mono=e[0]+axisCount*e[1];const double* c=&m_quadAffineCoeff[((size_t)host*m_quadAffineMonoCount+mono)*np];
                    for(int k=0;k<np;++k)dinn[ls]+=dcoeff[(size_t)ls*np+k]*mv[k]+c[k]*dm[k];}
            }
            const double wg=m_gwOut[ix]*m_gwOut[iy];for(int i=0;i<n;++i){const double w=wg*HexMonoEval(g[i],xi);for(int j=0;j<n;++j)out[(size_t)i*n+j]+=w*dinn[j];}
        }}
    } else {
    for(int sA=0;sA<2;++sA){
        const size_t sidA=(size_t)host*2+sA; const double* subVA=&m_faceSubV[sidA*9],*cA=&m_faceSubC[sidA*3];
        for(int sB=0;sB<2;++sB){
            const size_t sidB=(size_t)host*2+sB; const double* cB=&m_faceSubC[sidB*3];
            const double dx=cA[0]-cB[0],dy=cA[1]-cB[1],dz=cA[2]-cB[2];
            const bool near=std::sqrt(dx*dx+dy*dy+dz*dz)<=m_near_grade*(m_faceSubS[sidA]+m_faceSubS[sidB]);
            int corner=0;double best=1e300;for(int i=0;i<3;++i){const double x=subVA[3*i]-cB[0],y=subVA[3*i+1]-cB[1],z=subVA[3*i+2]-cB[2],d=x*x+y*y+z*z;if(d<best){best=d;corner=i;}}
            std::shared_ptr<const HexQuadCloud> oc;
            if(near)oc=HexGetCloud(m_build_id,HexCloudKey(1,true,true,host,sA,corner),[&](HexQuadCloud& c){std::vector<double> gb,gw;HexDuffyBary(2,corner,m_glOut,m_gwOut,gb,gw);HexBuildCloud(nd,false,sA,gb.data(),gw.data(),(int)gw.size(),true,c);});
            else oc=HexGetCloud(m_build_id,HexCloudKey(1,true,false,host,sA,3),[&](HexQuadCloud& c){HexBuildCloud(nd,false,sA,m_symTriP.data(),m_symTriW.data(),(int)m_symTriW.size(),false,c);});
            for(int q=0;q<(int)oc->wgeo.size();++q){const double* xi=&oc->xi[3*q];const double p[3]={oc->pts[3*q],oc->pts[3*q+1],oc->pts[3*q+2]};
                const double uv[2]={xi[0],xi[1]};double dp[3];QuadQ2MapX(node_velocity.data(),uv,dp);std::fill(dinn.begin(),dinn.end(),0.0);
                DPhiInnerHexRadialFaceVec(host,sB,p,dp,xi,node_velocity.data(),g,dinn.data());
                for(int i=0;i<n;++i){const double w=oc->wgeo[q]*HexMonoEval(g[i],xi);for(int j=0;j<n;++j)out[(size_t)i*n+j]+=w*dinn[j];}
            }
        }
    }
    }
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){const double v=.5*(out[(size_t)i*n+j]+out[(size_t)j*n+i]);out[(size_t)i*n+j]=out[(size_t)j*n+i]=v;}
    for(double& v:out)v*=RAD_INV_FOUR_PI; return out;
}

// thread_local block cache (build_id-guarded, same discipline as the cloud cache).  Keyed by the directed
// (kindT,hT,kindS,hS); a HACApK dense leaf touches all nT*nS entries of a host pair -> computed once, reused.
struct HexBlockKey {
    int kindT;
    int hT;
    int kindS;
    int hS;
    int img;                                           // IMA image index (0 = direct block)
    bool operator==(const HexBlockKey& o) const
    {
        return kindT == o.kindT && hT == o.hT && kindS == o.kindS && hS == o.hS && img == o.img;
    }
};

struct HexBlockKeyHash {
    std::size_t operator()(const HexBlockKey& k) const
    {
        std::size_t h = 1469598103934665603ull;
        auto mix = [&](int v) {
            h ^= static_cast<std::size_t>(static_cast<unsigned int>(v));
            h *= 1099511628211ull;
        };
        mix(k.kindT); mix(k.hT); mix(k.kindS); mix(k.hS); mix(k.img);
        return h;
    }
};

struct HexTransBlockKey {
    int kindT;
    int tmplT;
    int kindS;
    int tmplS;
    int dx2;
    int dy2;
    int dz2;
    int sameHost;
    bool operator==(const HexTransBlockKey& o) const {
        return kindT == o.kindT && tmplT == o.tmplT && kindS == o.kindS && tmplS == o.tmplS &&
               dx2 == o.dx2 && dy2 == o.dy2 && dz2 == o.dz2 && sameHost == o.sameHost;
    }
};

struct HexTransBlockKeyHash {
    std::size_t operator()(const HexTransBlockKey& k) const
    {
        std::size_t h = 1469598103934665603ull;
        auto mix = [&](int v) {
            h ^= static_cast<std::size_t>(static_cast<unsigned int>(v));
            h *= 1099511628211ull;
        };
        mix(k.kindT); mix(k.tmplT); mix(k.kindS); mix(k.tmplS);
        mix(k.dx2); mix(k.dy2); mix(k.dz2); mix(k.sameHost);
        return h;
    }
};

static thread_local long long s_hex_block_owner = -1;
static thread_local std::unordered_map<HexBlockKey, std::vector<double>, HexBlockKeyHash> s_hex_block_cache;
static thread_local long long s_hex_trans_block_owner = -1;
static thread_local std::unordered_map<HexTransBlockKey, std::vector<double>, HexTransBlockKeyHash> s_hex_trans_block_cache;
static thread_local long long s_hex_sym_block_owner = -1;
static thread_local std::unordered_map<HexBlockKey, std::vector<double>, HexBlockKeyHash> s_hex_sym_block_cache;
static thread_local long long s_hex_sym_trans_block_owner = -1;
static thread_local std::unordered_map<HexTransBlockKey, std::vector<double>, HexTransBlockKeyHash> s_hex_sym_trans_block_cache;
static thread_local long long s_ho_tet_sym_block_owner = -1;
static thread_local std::unordered_map<HexBlockKey, std::vector<double>, HexBlockKeyHash> s_ho_tet_sym_block_cache;

// True when the pair would take QuadBlockHex's general path (graded near ~30 ms/block, site-radial
// mid band ~0.5 ms/block): neither the affine-affine product nor the distorted-pair far product
// applies.  Used ONLY to route the block into the instance-shared general cache -- a borderline
// disagreement with QuadBlockHex's own dispatch merely caches a block in the other tier (see the
// header member doc), so this test may use the cheap ctor-precomputed affinity flags
// (m_hexAffineCell / m_quadAffineFace) instead of QuadBlockHex's per-call face_affine lattice
// re-check.
bool RadHACApKChargeGram::HexPairTakesGeneralPath(int kindT, int hT, int kindS, int hS, int img) const
{
    const std::vector<int>& tgtG = (kindT == 0) ? m_cellCharges[hT] : m_faceCharges[hT];
    const std::vector<int>& srcG = (kindS == 0) ? m_cellCharges[hS] : m_faceCharges[hS];
    if (tgtG.empty() || srcG.empty()) return false;
    const bool affT = (kindT == 0)
        ? (hT >= 0 && hT < (int)m_hexAffineCell.size() && m_hexAffineCell[hT])
        : (hT >= 0 && hT < (int)m_quadAffineFace.size() && m_quadAffineFace[hT]);
    const bool affS = (kindS == 0)
        ? (hS >= 0 && hS < (int)m_hexAffineCell.size() && m_hexAffineCell[hS])
        : (hS >= 0 && hS < (int)m_quadAffineFace.size() && m_quadAffineFace[hS]);
    if (affT && affS) return false;
    const double fac = HexDistortedFarFactor();
    if (fac <= 0.0) return true;
    const int repA = tgtG[0], repB = srcG[0];
    double repBc[3];
    ImageEvalPoint(img, &m_cent[(size_t)3*repB], repBc);
    const double sep = std::sqrt(
        (m_cent[(size_t)3*repA]     - repBc[0])*(m_cent[(size_t)3*repA]     - repBc[0])
      + (m_cent[(size_t)3*repA + 1] - repBc[1])*(m_cent[(size_t)3*repA + 1] - repBc[1])
      + (m_cent[(size_t)3*repA + 2] - repBc[2])*(m_cent[(size_t)3*repA + 2] - repBc[2]));
    return sep <= fac*(m_size[repA] + m_size[repB]);
}

const std::vector<double>& RadHACApKChargeGram::GetHexBlock(int kindT, int hT, int kindS, int hS, int img) const
{
    const int wedge_scope = m_wedgemode ? WedgeTransCacheScope() : 2;
    const bool use_trans_cache = std::getenv("RADIA_HDIV_DISABLE_TRANS_CACHE") == nullptr &&
                                 !m_d2 && m_hexUniformTransHosts && img == 0 &&
                                 (!m_wedgemode || wedge_scope >= 2 || (kindT == 0 && kindS == 0));
    if (!use_trans_cache && !m_d2 && !m_wedgemode
            && HexPairTakesGeneralPath(kindT, hT, kindS, hS, img)) {
        // Expensive general-path block (graded near or site-radial mid band): serve from the
        // instance-shared cache so the fill workers do not each recompute it (the per-thread caches
        // below duplicate ~2x).  Key packs (kinds, image index, hosts) into 64 bits; hosts < 2^28 and
        // the image index < 64 (a 63-fold cyclic group is far past any practical machine sector count).
        if (img < 0 || img > 63)
            throw std::invalid_argument("ChargeGram: image index out of range (max 63 images)");
        const unsigned long long key =
            (unsigned long long)(kindT != 0) | ((unsigned long long)(kindS != 0) << 1)
          | ((unsigned long long)(unsigned)img << 2)
          | ((unsigned long long)(unsigned)hT << 8) | ((unsigned long long)(unsigned)hS << 36);
        m_hexGeneralSharedLookups.fetch_add(1, std::memory_order_relaxed);
        {
            std::shared_lock<std::shared_mutex> rl(m_hexGeneralSharedMutex);
            auto it = m_hexGeneralSharedCache.find(key);
            if (it != m_hexGeneralSharedCache.end()) {
                m_hexGeneralSharedHits.fetch_add(1, std::memory_order_relaxed);
                return it->second;
            }
        }
        std::vector<double> blk = QuadBlockHex(kindT, hT, kindS, hS, img);
        m_hexGeneralSharedMisses.fetch_add(1, std::memory_order_relaxed);
        std::unique_lock<std::shared_mutex> wl(m_hexGeneralSharedMutex);
        auto it = m_hexGeneralSharedCache.emplace(key, std::move(blk)).first;   // racing first insert wins
        return it->second;
    }
    if (use_trans_cache) {
        HexStatAdd(m_hexCacheStatsEnabled, m_hexTransBlockLookups);
        if (s_hex_trans_block_owner != m_build_id) {
            s_hex_trans_block_cache.clear();
            s_hex_trans_block_owner = m_build_id;
            HexStatAdd(m_hexCacheStatsEnabled, m_hexTransBlockClears);
        }
        const size_t itHost = (kindT == 0) ? (size_t)hT : (size_t)m_n_el + (size_t)hT;
        const size_t isHost = (kindS == 0) ? (size_t)hS : (size_t)m_n_el + (size_t)hS;
        const int* lt = &m_hexHostLattice2[3*itHost];
        const int* ls = &m_hexHostLattice2[3*isHost];
        const HexTransBlockKey tkey{kindT, m_hexHostTemplate[itHost], kindS, m_hexHostTemplate[isHost],
                                    lt[0] - ls[0], lt[1] - ls[1], lt[2] - ls[2],
                                    (kindT == kindS && hT == hS) ? 1 : 0};
        auto it = s_hex_trans_block_cache.find(tkey);
        if (it == s_hex_trans_block_cache.end()) {
            HexStatAdd(m_hexCacheStatsEnabled, m_hexTransBlockMisses);
            if (s_hex_trans_block_cache.size() > HexBlockCacheLimit()) {
                s_hex_trans_block_cache.clear();
                HexStatAdd(m_hexCacheStatsEnabled, m_hexTransBlockClears);
            }
            it = s_hex_trans_block_cache.emplace(tkey,
                    m_wedgemode ? QuadBlockWedge(kindT, hT, kindS, hS, img)
                                : QuadBlockHex(kindT, hT, kindS, hS, img)).first;
        } else HexStatAdd(m_hexCacheStatsEnabled, m_hexTransBlockHits);
        return it->second;
    }
    HexStatAdd(m_hexCacheStatsEnabled, m_hexBlockLookups);
    if (s_hex_block_owner != m_build_id) {
        s_hex_block_cache.clear();
        s_hex_block_owner = m_build_id;
        HexStatAdd(m_hexCacheStatsEnabled, m_hexBlockClears);
    }
    const HexBlockKey key{kindT, hT, kindS, hS, img};
    auto it = s_hex_block_cache.find(key);
    if (it == s_hex_block_cache.end()) {
        HexStatAdd(m_hexCacheStatsEnabled, m_hexBlockMisses);
        if (s_hex_block_cache.size() > HexBlockCacheLimit()) {
            s_hex_block_cache.clear();
            HexStatAdd(m_hexCacheStatsEnabled, m_hexBlockClears);
        }
        it = s_hex_block_cache.emplace(key,
                   m_d2        ? QuadBlock2D(kindT, hT, kindS, hS, img)
                 : m_wedgemode ? QuadBlockWedge(kindT, hT, kindS, hS, img)
                 :               QuadBlockHex(kindT, hT, kindS, hS, img)).first;
    } else HexStatAdd(m_hexCacheStatsEnabled, m_hexBlockHits);
    return it->second;
}

// GetHexSymBlock(A,B,img) = 0.5*(G_T(A,B) + G_T(B,A)^T) = 0.5*(G_T(A,B) + G_{T^-1}(A,B)), because the
// directed BA block with transform T is the AB block with T^-1 (isometry).  For a MIRROR T == T^-1, so this is
// a pure quadrature symmetrization.  For a cyclic ROTATION the two differ, and the total Sum_i s_i * symblock_i
// equals the intended Sum_i s_i * G_{T_i} only when the image set is closed under inversion with matching
// signs -- which SetImageRotations enforces (a full cyclic group is closed automatically).
const std::vector<double>& RadHACApKChargeGram::GetHexSymBlock(int kindA, int hA, int kindB, int hB, int img) const
{
    auto nlocal = [&](int kind, int host) {
        return (kind == 0) ? (int)m_cellCharges[host].size() : (int)m_faceCharges[host].size();
    };
    const int nA0 = nlocal(kindA, hA);
    const int nB0 = nlocal(kindB, hB);
    auto make_sym = [&]() {
        std::vector<double> ab = GetHexBlock(kindA, hA, kindB, hB, img);
        std::vector<double> ba = GetHexBlock(kindB, hB, kindA, hA, img);
        std::vector<double> sym((size_t)nA0 * (size_t)nB0);
        for (int la = 0; la < nA0; ++la)
            for (int lb = 0; lb < nB0; ++lb)
                sym[(size_t)la*nB0 + lb] =
                    0.5 * (ab[(size_t)la*nB0 + lb] + ba[(size_t)lb*nA0 + la]);
        return sym;
    };
    auto far_one_sided = [&]() {
        // HEX/WEDGE host-pair FAR blocks are smooth low-order double integrals in physical space.  The
        // A directed AB block is the transpose of BA only in exact arithmetic.  Finite product quadrature is
        // orientation-dependent, so production averages both directions to preserve explicit-reflection
        // invariance.  A positive environment threshold is retained only for diagnostic timing experiments.
        const double threshold = m_wedgemode ? WedgeFarOneSidedThreshold() : HexFarOneSidedThreshold();
        if (threshold <= 0.0 || img != 0 || m_d2) return false;
        const std::vector<int>& aG = (kindA == 0) ? m_cellCharges[hA] : m_faceCharges[hA];
        const std::vector<int>& bG = (kindB == 0) ? m_cellCharges[hB] : m_faceCharges[hB];
        if (aG.empty() || bG.empty()) return false;
        const int a = aG[0], b = bG[0];
        const double dx = m_cent[3*a] - m_cent[3*b];
        const double dy = m_cent[3*a + 1] - m_cent[3*b + 1];
        const double dz = m_cent[3*a + 2] - m_cent[3*b + 2];
        return std::sqrt(dx*dx + dy*dy + dz*dz) > threshold * (m_size[a] + m_size[b]);
    };
    auto make_one_sided = [&]() {
        return GetHexBlock(kindA, hA, kindB, hB, img);
    };

    const int wedge_scope = m_wedgemode ? WedgeTransCacheScope() : 2;
    const bool use_trans_cache = std::getenv("RADIA_HDIV_DISABLE_TRANS_CACHE") == nullptr &&
                                 !m_d2 && m_hexUniformTransHosts && img == 0 &&
                                 (!m_wedgemode || wedge_scope >= 2 || (kindA == 0 && kindB == 0));
    if (use_trans_cache) {
        HexStatAdd(m_hexCacheStatsEnabled, m_hexSymTransBlockLookups);
        if (s_hex_sym_trans_block_owner != m_build_id) {
            s_hex_sym_trans_block_cache.clear();
            s_hex_sym_trans_block_owner = m_build_id;
            HexStatAdd(m_hexCacheStatsEnabled, m_hexSymTransBlockClears);
        }
        const size_t iaHost = (kindA == 0) ? (size_t)hA : (size_t)m_n_el + (size_t)hA;
        const size_t ibHost = (kindB == 0) ? (size_t)hB : (size_t)m_n_el + (size_t)hB;
        const int* la = &m_hexHostLattice2[3*iaHost];
        const int* lb = &m_hexHostLattice2[3*ibHost];
        const HexTransBlockKey key{kindA, m_hexHostTemplate[iaHost], kindB, m_hexHostTemplate[ibHost],
                                   la[0] - lb[0], la[1] - lb[1], la[2] - lb[2],
                                   (kindA == kindB && hA == hB) ? 1 : 0};
        auto it = s_hex_sym_trans_block_cache.find(key);
        if (it == s_hex_sym_trans_block_cache.end()) {
            HexStatAdd(m_hexCacheStatsEnabled, m_hexSymTransBlockMisses);
            if (s_hex_sym_trans_block_cache.size() > HexBlockCacheLimit()) {
                s_hex_sym_trans_block_cache.clear();
                HexStatAdd(m_hexCacheStatsEnabled, m_hexSymTransBlockClears);
            }
            it = s_hex_sym_trans_block_cache.emplace(key, far_one_sided() ? make_one_sided() : make_sym()).first;
        } else {
            HexStatAdd(m_hexCacheStatsEnabled, m_hexSymTransBlockHits);
        }
        return it->second;
    }

    HexStatAdd(m_hexCacheStatsEnabled, m_hexSymBlockLookups);
    if (s_hex_sym_block_owner != m_build_id) {
        s_hex_sym_block_cache.clear();
        s_hex_sym_block_owner = m_build_id;
        HexStatAdd(m_hexCacheStatsEnabled, m_hexSymBlockClears);
    }
    const HexBlockKey key{kindA, hA, kindB, hB, img};
    auto it = s_hex_sym_block_cache.find(key);
    if (it == s_hex_sym_block_cache.end()) {
        HexStatAdd(m_hexCacheStatsEnabled, m_hexSymBlockMisses);
        if (s_hex_sym_block_cache.size() > HexBlockCacheLimit()) {
            s_hex_sym_block_cache.clear();
            HexStatAdd(m_hexCacheStatsEnabled, m_hexSymBlockClears);
        }
        it = s_hex_sym_block_cache.emplace(key, far_one_sided() ? make_one_sided() : make_sym()).first;
    } else {
        HexStatAdd(m_hexCacheStatsEnabled, m_hexSymBlockHits);
    }
    return it->second;
}

const std::vector<double>& RadHACApKChargeGram::GetHOTetSymBlock(
    int kindA, int hostA, int kindB, int hostB) const
{
    HexStatAdd(m_hexCacheStatsEnabled, m_hoSymBlockLookups);
    if (s_ho_tet_sym_block_owner != m_build_id) {
        s_ho_tet_sym_block_cache.clear();
        s_ho_tet_sym_block_owner = m_build_id;
        HexStatAdd(m_hexCacheStatsEnabled, m_hoSymBlockClears);
    }
    const HexBlockKey key{kindA, hostA, kindB, hostB, 0};
    auto it = s_ho_tet_sym_block_cache.find(key);
    if (it == s_ho_tet_sym_block_cache.end()) {
        HexStatAdd(m_hexCacheStatsEnabled, m_hoSymBlockMisses);
        if (s_ho_tet_sym_block_cache.size() > HexBlockCacheLimit()) {
            s_ho_tet_sym_block_cache.clear();
            HexStatAdd(m_hexCacheStatsEnabled, m_hoSymBlockClears);
        }
        const int nA = (kindA == 0) ? (int)m_hoCellCharges[hostA].size()
                                    : (int)m_hoFaceCharges[hostA].size();
        const int nB = (kindB == 0) ? (int)m_hoCellCharges[hostB].size()
                                    : (int)m_hoFaceCharges[hostB].size();
        std::vector<double> ab = QuadBlockHOTet(kindA, hostA, kindB, hostB);
        const bool direct_curved = m_curved && CurvedDirectEnabled() &&
                                   !CurvedHostsTouch(kindA, hostA, kindB, hostB);
        const bool same_host = kindA == kindB && hostA == hostB;
        std::vector<double> ba;
        if (!direct_curved && !same_host)
            ba = QuadBlockHOTet(kindB, hostB, kindA, hostA);
        std::vector<double> sym((size_t)nA*nB);
        for (int localA = 0; localA < nA; ++localA)
            for (int localB = 0; localB < nB; ++localB)
                sym[(size_t)localA*nB + localB] = direct_curved
                    ? ab[(size_t)localA*nB + localB]
                    : 0.5*(ab[(size_t)localA*nB + localB] +
                          (same_host ? ab[(size_t)localB*nA + localA]
                                     : ba[(size_t)localB*nA + localA]));
        it = s_ho_tet_sym_block_cache.emplace(key, std::move(sym)).first;
        if (!same_host) {
            // The symmetrized reverse block is exactly the transpose.  Store
            // it now so a later H-matrix leaf with the opposite cluster
            // ordering does not repeat both expensive directed integrations.
            const HexBlockKey reverse_key{kindB, hostB, kindA, hostA, 0};
            if (s_ho_tet_sym_block_cache.find(reverse_key) ==
                s_ho_tet_sym_block_cache.end()) {
                std::vector<double> reverse((size_t)nB*nA);
                const auto& forward = it->second;
                for (int localA = 0; localA < nA; ++localA)
                    for (int localB = 0; localB < nB; ++localB)
                        reverse[(size_t)localB*nA + localA] =
                            forward[(size_t)localA*nB + localB];
                s_ho_tet_sym_block_cache.emplace(reverse_key, std::move(reverse));
                it = s_ho_tet_sym_block_cache.find(key);
            }
        }
    } else HexStatAdd(m_hexCacheStatsEnabled, m_hoSymBlockHits);
    return it->second;
}

// ===================================================================== 2D PLANAR mode (2026-07-03)
// Motor cross-section layer (memory hdiv-vim-tri-quad-motor): kernel -ln(r)/(2pi), charges = -div M on
// tri/quad cells + M.n on boundary edges, Piola-exact REF measures, regular symmetric outer (the log
// kernel's single-layer potentials are continuous -- numpy-validated that NO graded outer is needed),
// radial-cone inner for near/self, cheap far cloud otherwise.  See the header ctor doc.
static const double D2_TRIREF_V[3][2] = {{1, 0}, {0, 1}, {0, 0}};   // NGSolve trig reference

static inline double D2Pow(double x, int exponent)
{
    if (exponent == 0) return 1.0;
    if (exponent == 1) return x;
    if (exponent == 2) return x*x;
    return x*x*x;
}

static inline double D2MonoCell(const int* e, const double xi[2])
{
    return D2Pow(xi[0], e[0]) * D2Pow(xi[1], e[1]);
}

static inline double D2MonoEdge(const int* e, double t)
{
    return D2Pow(t, e[0]);
}

void RadHACApKChargeGram::D2CellMap(
    int cell_type, const double* coeff, const double xi[2], double X[2]) const
{
    X[0] = X[1] = 0.0;
    int k = 0;
    if (cell_type == 0) {
        for (int total = 0; total <= m_d2GeometryOrder; ++total)
            for (int j = 0; j <= total; ++j, ++k) {
                const int i = total - j;
                const double s = D2Pow(xi[0], i)*D2Pow(xi[1], j);
                X[0] += s*coeff[2*k]; X[1] += s*coeff[2*k + 1];
            }
    } else {
        for (int j = 0; j <= m_d2GeometryOrder; ++j)
            for (int i = 0; i <= m_d2GeometryOrder; ++i, ++k) {
                const double s = D2Pow(xi[0], i)*D2Pow(xi[1], j);
                X[0] += s*coeff[2*k]; X[1] += s*coeff[2*k + 1];
            }
    }
}

void RadHACApKChargeGram::D2EdgeMap(const double* coeff, double t, double X[2]) const
{
    X[0] = X[1] = 0.0;
    double p = 1.0;
    for (int k = 0; k <= m_d2GeometryOrder; ++k, p *= t) {
        X[0] += p*coeff[2*k]; X[1] += p*coeff[2*k + 1];
    }
}

void RadHACApKChargeGram::D2EdgeTangent(const double* coeff, double t, double T[2]) const
{
    T[0] = T[1] = 0.0;
    double p = 1.0;
    for (int k = 1; k <= m_d2GeometryOrder; ++k, p *= t) {
        T[0] += k*p*coeff[2*k]; T[1] += k*p*coeff[2*k + 1];
    }
}

static void ValidateD2GeometryOrder(int geometry_order)
{
    if (geometry_order < 1 || geometry_order > 3)
        throw std::invalid_argument("2D ChargeGram: geometry_order must be in {1,2,3}");
}

static int D2CellNSub(int cell_type)
{
    return (cell_type == 1) ? 4 : 1;
}

// Sub-triangle reference vertices.  A quad uses the four-triangle centre fan,
// not a fixed diagonal: the set is invariant under every square reflection and
// therefore preserves IMA parity in the assembled Gram operator.
static void D2SubTri(int cell_type, int s, double V[3][2])
{
    if (cell_type == 0) {
        for (int i = 0; i < 3; ++i) { V[i][0] = D2_TRIREF_V[i][0]; V[i][1] = D2_TRIREF_V[i][1]; }
    } else {
        const int a = s & 3, b = (s + 1) & 3;
        V[0][0] = 0.5; V[0][1] = 0.5;
        V[1][0] = QUADREF_V[a][0]; V[1][1] = QUADREF_V[a][1];
        V[2][0] = QUADREF_V[b][0]; V[2][1] = QUADREF_V[b][1];
    }
}

// anchor site k (0-6) of a ref sub-tri: corners, edge mids ((0,1),(1,2),(2,0)), centroid
static void D2SiteRef(const double V[3][2], int k, double x0[2])
{
    if (k < 3)      { x0[0] = V[k][0]; x0[1] = V[k][1]; }
    else if (k < 6) {
        const int a = k - 3, b = (k - 2) % 3;
        x0[0] = 0.5*(V[a][0] + V[b][0]); x0[1] = 0.5*(V[a][1] + V[b][1]);
    } else          { x0[0] = (V[0][0]+V[1][0]+V[2][0])/3.0; x0[1] = (V[0][1]+V[1][1]+V[2][1])/3.0; }
}

RadHACApKChargeGram::RadHACApKChargeGram(int /*dim2_tag*/, int geometry_order,
    std::vector<double> cell_map, std::vector<int> cell_type, std::vector<double> edge_map,
    int n_el, int n_be,
    std::vector<int> charge_host, std::vector<int> charge_kind, std::vector<int> charge_expo,
    std::vector<double> sym_tri_pts, std::vector<double> sym_tri_w,
    std::vector<double> gl_quad, std::vector<double> gw_quad,
    std::vector<double> gl_edge, std::vector<double> gw_edge,
    std::vector<double> gl_in, std::vector<double> gw_in,
    std::vector<double> far_tri_pts, std::vector<double> far_tri_w,
    double near_grade, double far_inner_factor,
    std::vector<int> image_masks, std::vector<double> image_signs)
    : m_n_el(n_el),
      m_glIn(std::move(gl_in)), m_gwIn(std::move(gw_in)),
      m_near_grade(near_grade), m_far_inner_factor(far_inner_factor),
      m_image_masks(std::move(image_masks)), m_image_signs(std::move(image_signs)),
      m_host(std::move(charge_host)), m_kind(std::move(charge_kind)), m_expo(std::move(charge_expo))
{
    ValidateImageVectors(m_image_masks, m_image_signs);
    ValidateD2GeometryOrder(geometry_order);
    m_d2 = true;
    m_d2_n_be = n_be;
    m_d2GeometryOrder = geometry_order;
    m_d2CellMapStride = 2*(geometry_order + 1)*(geometry_order + 1);
    m_d2EdgeMapStride = 2*(geometry_order + 1);
    m_d2CellMap = std::move(cell_map);
    m_d2CellType  = std::move(cell_type);
    m_d2EdgeMap = std::move(edge_map);
    if (m_d2CellType.size() != static_cast<size_t>(n_el)
            || m_d2CellMap.size() != static_cast<size_t>(n_el*m_d2CellMapStride)
            || m_d2EdgeMap.size() != static_cast<size_t>(n_be*m_d2EdgeMapStride))
        throw std::invalid_argument("2D ChargeGram: geometry map array size mismatch");
    for (int ct : m_d2CellType)
        if (ct < 0 || ct > 1)
            throw std::invalid_argument("2D ChargeGram: cell_type must be 0 (tri) or 1 (quad)");
    m_d2SymTriP = std::move(sym_tri_pts); m_d2SymTriW = std::move(sym_tri_w);
    m_d2GlQ = std::move(gl_quad); m_d2GwQ = std::move(gw_quad);
    m_d2GlE = std::move(gl_edge); m_d2GwE = std::move(gw_edge);
    m_d2FarTriP = std::move(far_tri_pts); m_d2FarTriW = std::move(far_tri_w);
    if (m_d2GlQ.empty() || m_d2GlQ.size() != m_d2GwQ.size())
        throw std::invalid_argument("2D ChargeGram: quad outer Gauss nodes/weights must be non-empty and equal-sized");
    m_n = (int)m_host.size();
    m_build_id = NextChargeGramBuildId();
    m_hexCacheStatsEnabled = HexCacheStatsEnabledByEnv();
    // ---- per-sub geometry: centroid/size (near test) + mapped anchor sites ----
    m_d2CellSubC.assign((size_t)n_el*4*2, 0.0); m_d2CellSubS.assign((size_t)n_el*4, 0.0);
    m_d2CellSiteX.assign((size_t)n_el*4*7*2, 0.0);
    for (int c = 0; c < n_el; ++c) {
        const double* map = &m_d2CellMap[(size_t)c*m_d2CellMapStride];
        const int ct = m_d2CellType[c];
        const int nsub = D2CellNSub(ct);
        for (int s = 0; s < nsub; ++s) {
            double V[3][2];
            D2SubTri(ct, s, V);
            double cen[2] = {0, 0}, P[3][2];
            for (int i = 0; i < 3; ++i) {
                D2CellMap(ct, map, V[i], P[i]);
                cen[0] += P[i][0]/3.0; cen[1] += P[i][1]/3.0;
            }
            double* pc = &m_d2CellSubC[((size_t)c*4 + s)*2];
            pc[0] = cen[0]; pc[1] = cen[1];
            double sz = 0.0;
            for (int i = 0; i < 3; ++i) {
                const double dx = P[i][0]-cen[0], dy = P[i][1]-cen[1];
                sz = std::max(sz, std::sqrt(dx*dx + dy*dy));
            }
            m_d2CellSubS[(size_t)c*4 + s] = sz;
            for (int k = 0; k < 7; ++k) {
                double x0[2], X[2];
                D2SiteRef(V, k, x0);
                D2CellMap(ct, map, x0, X);
                double* out = &m_d2CellSiteX[(((size_t)c*4 + s)*7 + k)*2];
                out[0] = X[0]; out[1] = X[1];
            }
        }
    }
    m_d2EdgeC.assign((size_t)n_be*2, 0.0); m_d2EdgeS.assign((size_t)n_be, 0.0);
    m_d2EdgeSiteX.assign((size_t)n_be*3*2, 0.0);
    for (int f = 0; f < n_be; ++f) {
        const double* map = &m_d2EdgeMap[(size_t)f*m_d2EdgeMapStride];
        double P0[2], P1[2], Pm[2];
        D2EdgeMap(map, 0.0, P0); D2EdgeMap(map, 1.0, P1); D2EdgeMap(map, 0.5, Pm);
        m_d2EdgeC[(size_t)f*2] = Pm[0]; m_d2EdgeC[(size_t)f*2 + 1] = Pm[1];
        const double dx = P1[0]-P0[0], dy = P1[1]-P0[1];
        m_d2EdgeS[f] = 0.5*std::sqrt(dx*dx + dy*dy);
        double* sx = &m_d2EdgeSiteX[(size_t)f*3*2];
        sx[0] = P0[0]; sx[1] = P0[1]; sx[2] = P1[0]; sx[3] = P1[1]; sx[4] = Pm[0]; sx[5] = Pm[1];
    }
    // ---- per-charge centroid/size (cluster-tree points; z = 0) + (kind,host) reverse maps ----
    m_cent.assign((size_t)m_n*3, 0.0); m_size.assign((size_t)m_n, 0.0);
    for (int a = 0; a < m_n; ++a) {
        const int h = m_host[a];
        if (m_kind[a] == 0) {
            const int ct = m_d2CellType[h];
            const int nsub = D2CellNSub(ct);
            double cen[2] = {0, 0}, sz = 0.0;
            for (int s = 0; s < nsub; ++s) {
                cen[0] += m_d2CellSubC[((size_t)h*4 + s)*2] / nsub;
                cen[1] += m_d2CellSubC[((size_t)h*4 + s)*2 + 1] / nsub;
            }
            for (int s = 0; s < nsub; ++s) {
                const double dx = m_d2CellSubC[((size_t)h*4 + s)*2] - cen[0];
                const double dy = m_d2CellSubC[((size_t)h*4 + s)*2 + 1] - cen[1];
                sz = std::max(sz, m_d2CellSubS[(size_t)h*4 + s] + std::sqrt(dx*dx + dy*dy));
            }
            m_cent[3*a] = cen[0]; m_cent[3*a + 1] = cen[1]; m_size[a] = sz;
        } else {
            m_cent[3*a] = m_d2EdgeC[(size_t)h*2]; m_cent[3*a + 1] = m_d2EdgeC[(size_t)h*2 + 1];
            m_size[a] = m_d2EdgeS[h];
        }
    }
    m_hexLocalOf.assign((size_t)m_n, 0);
    m_cellCharges.assign((size_t)n_el, {}); m_faceCharges.assign((size_t)n_be, {});
    for (int a = 0; a < m_n; ++a) {
        std::vector<int>& grp = (m_kind[a] == 0) ? m_cellCharges[m_host[a]] : m_faceCharges[m_host[a]];
        m_hexLocalOf[a] = (int)grp.size();
        grp.push_back(a);
    }
    m_hex_state_sum = HexStateChecksum();   // instance-integrity canary (shared with the hex mode)
}

// 2D inner: INT over sub subB of source (kindS,hS) of m_b(eta)*(-ln|p-X(eta)|) d(ref eta).
void RadHACApKChargeGram::PhiInner2DVec(int kindS, int hS, int subB, const double p[2],
                                        const double* xiT, const std::vector<int>& srcG,
                                        double* inn) const
{
    const int nS = (int)srcG.size();
    const int nR = (int)m_glIn.size();
    const double* GL = m_glIn.data();
    const double* GW = m_gwIn.data();
    if (kindS == 0) {
        const int ct = m_d2CellType[hS];
        const double* map = &m_d2CellMap[(size_t)hS*m_d2CellMapStride];
        double V[3][2];
        D2SubTri(ct, subB, V);
        const double* cs = &m_d2CellSubC[((size_t)hS*4 + subB)*2];
        const double sz = m_d2CellSubS[(size_t)hS*4 + subB];
        const double dxc = p[0]-cs[0], dyc = p[1]-cs[1];
        if (std::sqrt(dxc*dxc + dyc*dyc) > m_far_inner_factor*sz) {
            // FAR: smooth -ln(r), the fixed bary tri rule mapped on the fly (2D is cheap)
            const int nq = (int)m_d2FarTriW.size();
            for (int q = 0; q < nq; ++q) {
                const double l1 = m_d2FarTriP[2*q], l2 = m_d2FarTriP[2*q + 1];
                const double xi[2] = {V[0][0] + l1*(V[1][0]-V[0][0]) + l2*(V[2][0]-V[0][0]),
                                      V[0][1] + l1*(V[1][1]-V[0][1]) + l2*(V[2][1]-V[0][1])};
                double X[2];
                D2CellMap(ct, map, xi, X);
                const double dx = p[0]-X[0], dy = p[1]-X[1];
                const double r = std::sqrt(dx*dx + dy*dy);
                if (r < 1e-300) continue;
                const double e1u = V[1][0]-V[0][0], e1v = V[1][1]-V[0][1];
                const double e2u = V[2][0]-V[0][0], e2v = V[2][1]-V[0][1];
                const double sc = std::fabs(e1u*e2v - e1v*e2u);          // 2*A_sub(ref); rule W sums 1/2
                const double g = 2.0*m_d2FarTriW[q]*0.5*sc*(-std::log(r));
                for (int ls = 0; ls < nS; ++ls) {
                    const int* e = &m_expo[(size_t)3*srcG[ls]];
                    inn[ls] += g*D2MonoCell(e, xi);
                }
            }
            return;
        }
        // NEAR/SELF: signed radial cones from the anchor (xiT on the self host, nearest site else)
        double x0[2];
        if (xiT) {
            const double xr[2] = {xiT[0], xiT[1]};
            ClosestPointTri2D(V, xr, x0);
        } else {
            const double* sx = &m_d2CellSiteX[(((size_t)hS*4 + subB)*7)*2];
            int best = 0; double bd = 1e300;
            for (int k = 0; k < 7; ++k) {
                const double dx = p[0]-sx[2*k], dy = p[1]-sx[2*k + 1];
                const double d = dx*dx + dy*dy;
                if (d < bd) { bd = d; best = k; }
            }
            D2SiteRef(V, best, x0);
        }
        for (int kf = 0; kf < 3; ++kf) {
            const double* A = V[kf]; const double* B = V[(kf + 1) % 3];
            const double ea[2] = {A[0]-x0[0], A[1]-x0[1]};
            const double eb[2] = {B[0]-x0[0], B[1]-x0[1]};
            const double s2 = ea[0]*eb[1] - ea[1]*eb[0];
            if (std::fabs(s2) < 1e-14) continue;
            for (int a2 = 0; a2 < nR; ++a2) { const double u = GL[a2];
                for (int b2 = 0; b2 < nR; ++b2) { const double v = GL[b2];
                    const double xi[2] = {x0[0] + u*(ea[0] + v*(eb[0]-ea[0])),
                                          x0[1] + u*(ea[1] + v*(eb[1]-ea[1]))};
                    double X[2];
                    D2CellMap(ct, map, xi, X);
                    const double dx = p[0]-X[0], dy = p[1]-X[1];
                    const double r = std::sqrt(dx*dx + dy*dy);
                    if (r < 1e-300) continue;
                    const double wq = GW[a2]*GW[b2]*(u*s2)*(-std::log(r));
                    for (int ls = 0; ls < nS; ++ls) {
                        const int* e = &m_expo[(size_t)3*srcG[ls]];
                        inn[ls] += wq*D2MonoCell(e, xi);
                    }
                }
            }
        }
        return;
    }
    // EDGE source: INT_0^1 t^e * (-ln|p-X(t)|) dt.  SELF (xiT set): the log singularity sits at the
    // OUTER point's own parameter t* -- split [0,t*] + [t*,1] and grade each piece INTO t* (s = t* -/+
    // len*g^2 turns the integrand into the smooth u*ln(u) class).  Near non-self: grade toward the
    // nearest endpoint (the projection of p); far: plain Gauss.
    const double* map = &m_d2EdgeMap[(size_t)hS*m_d2EdgeMapStride];
    const double* ec = &m_d2EdgeC[(size_t)hS*2];
    const double es = m_d2EdgeS[hS];
    const double dxc = p[0]-ec[0], dyc = p[1]-ec[1];
    const bool far_pt = !xiT && std::sqrt(dxc*dxc + dyc*dyc) > m_far_inner_factor*es;
    const int nq = (int)m_d2GwE.size();
    auto accum = [&](double t, double w) {
        double X[2];
        D2EdgeMap(map, t, X);
        const double dx = p[0]-X[0], dy = p[1]-X[1];
        const double r = std::sqrt(dx*dx + dy*dy);
        if (r < 1e-300) return;
        const double g2 = w*(-std::log(r));
        for (int ls = 0; ls < nS; ++ls) {
            const int* e = &m_expo[(size_t)3*srcG[ls]];
            inn[ls] += g2*D2MonoEdge(e, t);
        }
    };
    if (far_pt) {
        for (int q = 0; q < nq; ++q) accum(m_d2GlE[q], m_d2GwE[q]);
        return;
    }
    // NEAR/SELF: split-grade around the kernel peak's parameter ts -- xiT[0] on the self edge (exact),
    // else the PROJECTION of p onto the (quadratic) edge: chord initial guess + a short Newton on
    // (X(t)-p).X'(t) = 0.  Grading toward an ENDPOINT instead (the first implementation) mis-resolves
    // every near pair whose peak is interior (cell outer points facing their own boundary edge) -- that
    // overestimate is exactly what leaked eig > 1 on the structured quad mesh.
    double ts;
    if (xiT) {
        ts = std::min(1.0, std::max(0.0, xiT[0]));
    } else {
        double P0[2], P1[2];
        D2EdgeMap(map, 0.0, P0); D2EdgeMap(map, 1.0, P1);
        const double du = P1[0]-P0[0], dv = P1[1]-P0[1];
        const double L2 = du*du + dv*dv;
        ts = (L2 > 1e-300) ? ((p[0]-P0[0])*du + (p[1]-P0[1])*dv)/L2 : 0.5;
        ts = std::min(1.0, std::max(0.0, ts));
        for (int it = 0; it < 4; ++it) {
            double X[2], T[2];
            D2EdgeMap(map, ts, X); D2EdgeTangent(map, ts, T);
            const double X0 = X[0], X1 = X[1], T0 = T[0], T1 = T[1];
            const double g1 = (X0-p[0])*T0 + (X1-p[1])*T1;  // d/dt |X-p|^2 / 2
            const double h2 = T0*T0 + T1*T1;                // + curvature term dropped (small, quadratic map)
            if (h2 < 1e-300) break;
            ts = std::min(1.0, std::max(0.0, ts - g1/h2));
        }
    }
    for (int side = 0; side < 2; ++side) {
        const double len = side ? (1.0 - ts) : ts;
        if (len < 1e-14) continue;
        for (int q = 0; q < nq; ++q) {
            const double g = m_d2GlE[q];
            const double t = side ? (ts + len*g*g) : (ts - len*g*g);
            accum(t, 2.0*g*m_d2GwE[q]*len);
        }
    }
}

// Whole DIRECTED 2D host-pair block (target outer x source inner), 1/(2pi) folded.  Regular symmetric
// outer everywhere (numpy-validated: the log kernel needs no graded outer); the SELF host pair passes the
// outer point's own ref coords as the inner anchor.
std::vector<double> RadHACApKChargeGram::QuadBlock2D(
    int kindT, int hT, int kindS, int hS, int img) const
{
    const std::vector<int>& tgtG = (kindT == 0) ? m_cellCharges[hT] : m_faceCharges[hT];
    const std::vector<int>& srcG = (kindS == 0) ? m_cellCharges[hS] : m_faceCharges[hS];
    const int nT = (int)tgtG.size(), nS = (int)srcG.size();
    std::vector<double> blk((size_t)nT*nS, 0.0);
    if (nT == 0 || nS == 0) return blk;
    const bool self_pair = (img == 0 && kindT == kindS && hT == hS);
    const int nsubS = (kindS == 0) ? D2CellNSub(m_d2CellType[hS]) : 1;
    std::vector<double> inn(nS);
    auto accumulate = [&](const double xiA[2], double wg, const double Xp[2]) {
        // In-plane image transform: T^-1 applied to the eval point (the 2D rotation about +z is the
        // in-plane rotation, so ImageEvalPoint's z component is simply unused here).
        double Pin[3] = {Xp[0], Xp[1], 0.0}, Pout[3];
        ImageEvalPoint(img, Pin, Pout);
        const double Peval[2] = {Pout[0], Pout[1]};
        for (int sB = 0; sB < nsubS; ++sB) {
            for (int ls = 0; ls < nS; ++ls) inn[ls] = 0.0;
            PhiInner2DVec(kindS, hS, sB, Peval, self_pair ? xiA : nullptr, srcG, inn.data());
            for (int lt = 0; lt < nT; ++lt) {
                const int* e = &m_expo[(size_t)3*tgtG[lt]];
                const double ma = (kindT == 0) ? D2MonoCell(e, xiA) : D2MonoEdge(e, xiA[0]);
                double* row = &blk[(size_t)lt*nS];
                for (int ls = 0; ls < nS; ++ls) row[ls] += wg*ma*inn[ls];
            }
        }
    };
    if (kindT == 0) {
        const int ct = m_d2CellType[hT];
        const double* map = &m_d2CellMap[(size_t)hT*m_d2CellMapStride];
        if (ct == 1) {
            // Tensor Gauss is invariant under xi -> 1-xi / eta -> 1-eta and avoids both the
            // fixed-diagonal defect and the vertex-order dependence of a Duffy triangle rule.
            const int nq = (int)m_d2GlQ.size();
            for (int j = 0; j < nq; ++j) for (int i = 0; i < nq; ++i) {
                const double xiA[2] = {m_d2GlQ[i], m_d2GlQ[j]};
                double Xp[2];
                D2CellMap(ct, map, xiA, Xp);
                accumulate(xiA, m_d2GwQ[i]*m_d2GwQ[j], Xp);
            }
        } else {
            const int nq = (int)m_d2SymTriW.size();
            double V[3][2];
            D2SubTri(ct, 0, V);
            const double e1u = V[1][0]-V[0][0], e1v = V[1][1]-V[0][1];
            const double e2u = V[2][0]-V[0][0], e2v = V[2][1]-V[0][1];
            const double sc = std::fabs(e1u*e2v - e1v*e2u);
            for (int q = 0; q < nq; ++q) {
                const double l1 = m_d2SymTriP[2*q], l2 = m_d2SymTriP[2*q + 1];
                const double xiA[2] = {V[0][0] + l1*e1u + l2*e2u, V[0][1] + l1*e1v + l2*e2v};
                double Xp[2];
                D2CellMap(ct, map, xiA, Xp);
                accumulate(xiA, m_d2SymTriW[q]*sc, Xp);
            }
        }
    } else {
        const double* map = &m_d2EdgeMap[(size_t)hT*m_d2EdgeMapStride];
        const int nq = (int)m_d2GwE.size();
        for (int q = 0; q < nq; ++q) {
            const double t = m_d2GlE[q];
            const double xiA[2] = {t, 0.0};
            double Xp[2];
            D2EdgeMap(map, t, Xp);
            accumulate(xiA, m_d2GwE[q], Xp);
        }
    }
    const double INV2PI = 1.0/(2.0*3.14159265358979323846);
    for (double& v : blk) v *= INV2PI;
    return blk;
}

void RadHACApKChargeGram::ExtractCoordinates()
{
    m_n_elem = m_n;
    m_ndof   = m_n;
    m_coordinates = m_cent;   // [n*3] charge centroids (the cluster-tree points)
}

extern "C" void cHACApK_set_sym_fill(int flag);   // cHACApK_base.c (skip strictly-lower leaves at fill)

bool RadHACApKChargeGram::BuildHMatrix(const RadHACApKParams& params)
{
    // Symmetric fill: the Gram's applies all route through MatVecSym (see the header doc), so skip the
    // strictly-lower leaves at fill time -- ~2x build, identical upper leaves (MatVecSym bit-identical).
    // Set/reset around this ONE build; base BuildHMatrix returns bool (no exceptions cross the C fill).
    ResetHexCacheStats();
    m_derivativeAcaEps=params.aca_eps;
    m_derivativeMaxRank=params.max_rank;
    // Charge-basis normalization (see the header doc at MatVecSym): the
    // sigma pre-pass + the m_fillNormalized toggle live in OnBeforeBuild --
    // the base build fixes m_n via ExtractCoordinates FIRST, then calls
    // OnBeforeBuild, then fills; only the reset after the fill lives here.
    cHACApK_set_sym_fill(1);
    const bool ok = RadHACApKBase::BuildHMatrix(params);
    cHACApK_set_sym_fill(0);
    m_fillNormalized = false;
    return ok;
}

void RadHACApKChargeGram::ComputeChargeSigma()
{
    m_chargeSigma.assign((size_t)m_n, 1.0);
    m_chargeSigmaInv.assign((size_t)m_n, 1.0);
    {
        // C++ HACApK self-wrap policy: this pre-pass runs BEFORE the base
        // build's own region, so it stands up its own.
        ngcore::RegionTaskManager rtm(
            std::max(1, (int)ngcore::TaskManager::GetMaxThreads()));
        ngcore::ParallelFor(ngcore::IntRange(m_n), [&](size_t p) {
            const double d = GetInteractionMatrixElementRaw((int)p, (int)p);
            if (d > 0.0 && std::isfinite(d)) {
                const double s = std::sqrt(d);
                m_chargeSigma[p] = s;
                m_chargeSigmaInv[p] = 1.0 / s;
            }
        });
    }
    m_sigmaActive = true;
}

void RadHACApKChargeGram::MatVecSym(const std::vector<double>& x,
                                    std::vector<double>& y)
{
    if (!m_sigmaActive) { RadHACApKBase::MatVecSym(x, y); return; }
    // stored leaves hold Ghat = S^-1 G S^-1; the physical apply is
    // G x = S (Ghat (S x)).
    std::vector<double> xs((size_t)m_n);
    for (int p = 0; p < m_n; ++p)
        xs[(size_t)p] = x[(size_t)p] * m_chargeSigma[(size_t)p];
    RadHACApKBase::MatVecSym(xs, y);
    for (int p = 0; p < m_n; ++p)
        y[(size_t)p] *= m_chargeSigma[(size_t)p];
}

void RadHACApKChargeGram::MatVecSymMany(const std::vector<double>& x,
                                        int nrhs, std::vector<double>& y)
{
    if (!m_sigmaActive) { RadHACApKBase::MatVecSymMany(x, nrhs, y); return; }
    std::vector<double> xs(x.size());
    for (int r = 0; r < nrhs; ++r)
        for (int p = 0; p < m_n; ++p)
            xs[(size_t)r * m_n + p] =
                x[(size_t)r * m_n + p] * m_chargeSigma[(size_t)p];
    RadHACApKBase::MatVecSymMany(xs, nrhs, y);
    for (int r = 0; r < nrhs; ++r)
        for (int p = 0; p < m_n; ++p)
            y[(size_t)r * m_n + p] *= m_chargeSigma[(size_t)p];
}

void RadHACApKChargeGram::OnBeforeBuild()
{
    if (m_curved) PrecomputeCurvedTouchBlocks();
    // Charge-basis normalization: runs after the base build's
    // ExtractCoordinates (m_n is final here) and before the fill, so the
    // fill's entry oracle serves Ghat = S^-1 G S^-1 with O(1) dynamic range.
    ComputeChargeSigma();
    m_fillNormalized = true;
}

// ============================================ WEDGE (PRISM) BDM1 compute (2026-07-04) ===================
// A faithful mirror of the hex-mode compute path (BuildHexSiteTables / PhiInnerHex{Site,Sub,Radial}Vec /
// QuadBlockHex) with two structural changes: (1) the CELL is a prism -> 3 sub-tets (WEDGEREF_TETS), 18-node
// map (WedgeQ2MapX); (2) the boundary FACE is MIXED -> a per-face type (m_wFaceType) selects tri (1 sub-tri,
// 6-node TriSurfMap) vs quad (2 sub-tris, 9-node QuadQ2MapX -- reused from hex).  The block memo, cloud
// cache (HexQuadCloud / HexGetCloud), leaf helpers (HexMonoEval, HexDuffyBary, ClosestPointTet /
// ClosestPointTri2D, the HexSiteRad struct), and the whole solver surface are shared verbatim, so the
// golden hex path is byte-for-byte untouched.

// ref sub-tri vertices of face host of type `face_type`, sub `s` (tri: the whole ref tri; quad: 2 sub-tris)
static void WFaceSubTriRef(int face_type, int s, double V[3][2])
{
    if (face_type == 0) {
        for (int i = 0; i < 3; ++i) { V[i][0] = WTRIREF_V[i][0]; V[i][1] = WTRIREF_V[i][1]; }
    } else {
        const int* tv = QUADREF_TRIS[s];
        for (int i = 0; i < 3; ++i) { V[i][0] = QUADREF_V[tv[i]][0]; V[i][1] = QUADREF_V[tv[i]][1]; }
    }
}

// Materialize a quadrature cloud on sub-simplex `sub` of a wedge host (mirror of HexBuildCloud): cell ->
// WEDGEREF_TETS + WedgeQ2MapX (18-node); face -> WFaceSubTriRef + (tri: TriSurfMap 6-node / quad: QuadQ2MapX
// 9-node).  wgeo = ruleW * ref-measure scale (Piola: no |det J|).
static void WedgeBuildCloud(const double* nd, int kind, int face_type, int sub,
                            const double* baryP, const double* baryW, int nq, bool full_bary,
                            HexQuadCloud& out)
{
    const bool cell = (kind == 0);
    out.pts.resize((size_t)nq*3); out.wgeo.resize(nq); out.xi.resize((size_t)nq*3);
    if (cell) {
        const int* tv = WEDGEREF_TETS[sub];
        const double scale = WedgeSubSixVref(sub);
        for (int q = 0; q < nq; ++q) {
            double bary[4];
            if (full_bary) { for (int t = 0; t < 4; ++t) bary[t] = baryP[(size_t)4*q + t]; }
            else { double ls = 0.0; for (int t = 1; t < 4; ++t) { bary[t] = baryP[(size_t)3*q + (t-1)]; ls += bary[t]; } bary[0] = 1.0 - ls; }
            double xi[3] = {0, 0, 0};
            for (int t = 0; t < 4; ++t) for (int k = 0; k < 3; ++k) xi[k] += bary[t]*WEDGEREF_V[tv[t]][k];
            double X[3]; RadHACApKChargeGram::WedgeQ2MapX(nd, xi, X);
            for (int k = 0; k < 3; ++k) { out.pts[(size_t)3*q+k] = X[k]; out.xi[(size_t)3*q+k] = xi[k]; }
            out.wgeo[q] = baryW[q]*scale;
        }
    } else {
        double V[3][2]; WFaceSubTriRef(face_type, sub, V);
        const double scale = (face_type == 0) ? WTriSubTwoAref() : QuadSubTwoAref(sub);
        for (int q = 0; q < nq; ++q) {
            double bary[3];
            if (full_bary) { for (int t = 0; t < 3; ++t) bary[t] = baryP[(size_t)3*q + t]; }
            else { double ls = 0.0; for (int t = 1; t < 3; ++t) { bary[t] = baryP[(size_t)2*q + (t-1)]; ls += bary[t]; } bary[0] = 1.0 - ls; }
            double uv[2] = {0, 0};
            for (int t = 0; t < 3; ++t) for (int k = 0; k < 2; ++k) uv[k] += bary[t]*V[t][k];
            double X[3];
            if (face_type == 0) RadHACApKChargeGram::TriSurfMap(nd, uv, X);
            else                RadHACApKChargeGram::QuadQ2MapX(nd, uv, X);
            out.pts[(size_t)3*q] = X[0]; out.pts[(size_t)3*q+1] = X[1]; out.pts[(size_t)3*q+2] = X[2];
            out.xi[(size_t)3*q] = uv[0]; out.xi[(size_t)3*q+1] = uv[1]; out.xi[(size_t)3*q+2] = 0.0;
            out.wgeo[q] = baryW[q]*scale;
        }
    }
}

// Ref coords of anchor site k of cell sub-tet s (WEDGE ref frame): 0-3 corners, 4-9 edge mids, 10-13 face
// centers (HEXTET_FC order), 14 centroid -- identical layout to HexSiteRef.
static void WedgeCellSiteRef(int s, int k, double x0[3])
{
    const int* tv = WEDGEREF_TETS[s];
    double V[4][3];
    for (int i = 0; i < 4; ++i) for (int d = 0; d < 3; ++d) V[i][d] = WEDGEREF_V[tv[i]][d];
    static const int E[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
    if (k < 4)       for (int d = 0; d < 3; ++d) x0[d] = V[k][d];
    else if (k < 10) for (int d = 0; d < 3; ++d) x0[d] = 0.5*(V[E[k-4][0]][d] + V[E[k-4][1]][d]);
    else if (k < 14) { const int* f = HEXTET_FC[k-10]; for (int d = 0; d < 3; ++d) x0[d] = (V[f[0]][d]+V[f[1]][d]+V[f[2]][d])/3.0; }
    else             for (int d = 0; d < 3; ++d) x0[d] = 0.25*(V[0][d]+V[1][d]+V[2][d]+V[3][d]);
}

// Ref uv coords of anchor site k of a face sub-tri with explicit verts V: 0-2 corners, 3-5 edge mids, 6 centroid.
static void WTriSiteRef(const double V[3][2], int k, double u0[2])
{
    if (k < 3)      for (int d = 0; d < 2; ++d) u0[d] = V[k][d];
    else if (k < 6) for (int d = 0; d < 2; ++d) u0[d] = 0.5*(V[k-3][d] + V[(k-2)%3][d]);
    else            for (int d = 0; d < 2; ++d) u0[d] = (V[0][d]+V[1][d]+V[2][d])/3.0;
}

// Build the host-INDEPENDENT static-site radial tables for the wedge (mirror of BuildHexSiteTables): cell
// 3 sub-tets x 15 sites (18-wide shape S, 8-wide Q1 monomial M); tri-face 1 sub-tri x 7 sites (6-wide S,
// 4-wide M); quad-face 2 sub-tris x 7 sites (9-wide S, 4-wide M).  Plus the per-host mapped site positions.
void RadHACApKChargeGram::BuildWedgeSiteTables()
{
    const int nR = (int)m_glIn.size();
    const double* GL = m_glIn.data();
    const double* GW = m_gwIn.data();
    // ---- cell site tables (3 sub-tets) ----
    m_wCellSiteRad.assign(3*15, HexSiteRad());
    for (int s = 0; s < 3; ++s) {
        const int* tv = WEDGEREF_TETS[s];
        double V[4][3];
        for (int i = 0; i < 4; ++i) for (int d = 0; d < 3; ++d) V[i][d] = WEDGEREF_V[tv[i]][d];
        double E0[3], E1[3], E2[3];
        for (int d = 0; d < 3; ++d) { E0[d] = V[1][d]-V[0][d]; E1[d] = V[2][d]-V[0][d]; E2[d] = V[3][d]-V[0][d]; }
        const double hv = E0[0]*(E1[1]*E2[2]-E1[2]*E2[1]) - E0[1]*(E1[0]*E2[2]-E1[2]*E2[0])
                        + E0[2]*(E1[0]*E2[1]-E1[1]*E2[0]);
        const double sgnT = (hv >= 0.0) ? 1.0 : -1.0;
        for (int k = 0; k < 15; ++k) {
            HexSiteRad& R = m_wCellSiteRad[(size_t)s*15 + k];
            double x0[3]; WedgeCellSiteRef(s, k, x0);
            for (int f = 0; f < 4; ++f) {
                const double* b1 = V[HEXTET_FC[f][0]]; const double* b2 = V[HEXTET_FC[f][1]]; const double* b3 = V[HEXTET_FC[f][2]];
                double d1[3], d2[3], d3[3], e21[3], e32[3];
                for (int d = 0; d < 3; ++d) { d1[d] = b1[d]-x0[d]; d2[d] = b2[d]-x0[d]; d3[d] = b3[d]-x0[d]; e21[d] = b2[d]-b1[d]; e32[d] = b3[d]-b2[d]; }
                const double cr[3] = {d2[1]*d3[2]-d2[2]*d3[1], d2[2]*d3[0]-d2[0]*d3[2], d2[0]*d3[1]-d2[1]*d3[0]};
                const double D = d1[0]*cr[0] + d1[1]*cr[1] + d1[2]*cr[2];
                if (std::fabs(D) < 1e-12) continue;
                for (int a = 0; a < nR; ++a) { const double u = GL[a];
                    for (int b = 0; b < nR; ++b) { const double v = GL[b];
                        for (int c = 0; c < nR; ++c) { const double w = GL[c];
                            double y[3];
                            for (int d = 0; d < 3; ++d) y[d] = x0[d] + u*(d1[d] + v*(e21[d] + w*e32[d]));
                            R.w.push_back(sgnT*GW[a]*GW[b]*GW[c]*(u*u*v*D));
                            double st[6], vz[3], dz[3]; TriP2Shape(y[0], y[1], st); HexLag3(y[2], vz, dz);
                            for (int iz = 0; iz < 3; ++iz) for (int t = 0; t < 6; ++t) R.S.push_back(st[t]*vz[iz]);
                            const double px[3] = {1.0, y[0], y[0]*y[0]};
                            const double py[3] = {1.0, y[1], y[1]*y[1]};
                            const double pz[3] = {1.0, y[2], y[2]*y[2]};
                            for (int ez = 0; ez < 3; ++ez)
                                for (int ey = 0; ey < 3; ++ey)
                                    for (int ex = 0; ex < 3; ++ex)
                                        R.M.push_back(px[ex]*py[ey]*pz[ez]);
                        }
                    }
                }
            }
            R.nq = (int)R.w.size();
        }
    }
    // ---- face site tables: tri (1 sub-tri, 6-wide S) + quad (2 sub-tris, 9-wide S) ----
    auto build_face = [&](std::vector<HexSiteRad>& tab, int nsub, int face_type) {
        tab.assign((size_t)nsub*7, HexSiteRad());
        for (int s = 0; s < nsub; ++s) {
            double V[3][2]; WFaceSubTriRef(face_type, s, V);
            for (int k = 0; k < 7; ++k) {
                HexSiteRad& R = tab[(size_t)s*7 + k];
                double u0[2]; WTriSiteRef(V, k, u0);
                for (int kf = 0; kf < 3; ++kf) {
                    const double* A = V[kf]; const double* B = V[(kf+1)%3];
                    const double ea[2] = {A[0]-u0[0], A[1]-u0[1]};
                    const double eb[2] = {B[0]-u0[0], B[1]-u0[1]};
                    const double s2 = ea[0]*eb[1] - ea[1]*eb[0];
                    if (std::fabs(s2) < 1e-12) continue;
                    for (int a = 0; a < nR; ++a) { const double u = GL[a];
                        for (int b = 0; b < nR; ++b) { const double v = GL[b];
                            const double yu = u0[0] + u*(ea[0] + v*(eb[0]-ea[0]));
                            const double yv = u0[1] + u*(ea[1] + v*(eb[1]-ea[1]));
                            R.w.push_back(GW[a]*GW[b]*(u*s2));
                            if (face_type == 0) { double st[6]; TriP2Shape(yu, yv, st); for (int t = 0; t < 6; ++t) R.S.push_back(st[t]); }
                            else { double vu[3], duu[3], vv[3], dvu[3]; HexLag3(yu, vu, duu); HexLag3(yv, vv, dvu);
                                   for (int iv = 0; iv < 3; ++iv) for (int iu = 0; iu < 3; ++iu) R.S.push_back(vu[iu]*vv[iv]); }
                            const double pu[3] = {1.0, yu, yu*yu};
                            const double pv[3] = {1.0, yv, yv*yv};
                            for (int ev = 0; ev < 3; ++ev)
                                for (int eu = 0; eu < 3; ++eu) R.M.push_back(pu[eu]*pv[ev]);
                        }
                    }
                }
                R.nq = (int)R.w.size();
            }
        }
    };
    build_face(m_wFaceSiteRadTri, 1, 0);
    build_face(m_wFaceSiteRadQuad, 2, 1);
    // ---- mapped site positions per host (nearest-site pick is a physical distance test) ----
    m_wCellSiteX.assign((size_t)m_n_el*3*15*3, 0.0);
    for (int c = 0; c < m_n_el; ++c) {
        const double* nd = &m_wCellNodes[(size_t)c*54];
        for (int s = 0; s < 3; ++s)
            for (int k = 0; k < 15; ++k) {
                double x0[3], X[3]; WedgeCellSiteRef(s, k, x0); WedgeQ2MapX(nd, x0, X);
                double* out = &m_wCellSiteX[(((size_t)c*3 + s)*15 + k)*3];
                out[0] = X[0]; out[1] = X[1]; out[2] = X[2];
            }
    }
    m_wFaceSiteX.assign((size_t)m_wedge_n_bf*2*7*3, 0.0);
    for (int f = 0; f < m_wedge_n_bf; ++f) {
        const int ft = m_wFaceType[f];
        const int nsub = (ft == 0) ? 1 : 2;
        const double* nd = &m_wFaceNodes[(size_t)f*27];
        for (int s = 0; s < nsub; ++s) {
            double V[3][2]; WFaceSubTriRef(ft, s, V);
            for (int k = 0; k < 7; ++k) {
                double u0[2], X[3]; WTriSiteRef(V, k, u0);
                if (ft == 0) TriSurfMap(nd, u0, X); else QuadQ2MapX(nd, u0, X);
                double* out = &m_wFaceSiteX[(((size_t)f*2 + s)*7 + k)*3];
                out[0] = X[0]; out[1] = X[1]; out[2] = X[2];
            }
        }
    }
    m_hex_state_sum = HexStateChecksum();
}

// NON-SELF near inner (static-SITE radial): mirror of PhiInnerHexSiteVec with the mixed-face S/M widths.
void RadHACApKChargeGram::PhiInnerWedgeSiteVec(int kindS, int hS, int subB, const double p[3],
                                               const std::vector<int>& srcG, double* inn) const
{
    const bool cell = (kindS == 0);
    const int ft = cell ? -1 : m_wFaceType[hS];
    const double* sx = cell ? &m_wCellSiteX[(((size_t)hS*3 + subB)*15)*3]
                            : &m_wFaceSiteX[(((size_t)hS*2 + subB)*7)*3];
    const int nsite = cell ? 15 : 7;
    int best = 0; double bd = 1e300;
    const double site_tie_tol = srcG.empty() ? 0.0
        : 1e-13 * (m_size[srcG[0]] + 1.0) * (m_size[srcG[0]] + 1.0);
    for (int k = 0; k < nsite; ++k) {
        const double dx = p[0]-sx[3*k], dy = p[1]-sx[3*k+1], dz = p[2]-sx[3*k+2];
        const double d = dx*dx + dy*dy + dz*dz;
        if (d < bd - site_tie_tol) { bd = d; best = k; }
    }
    const HexSiteRad& R = cell ? m_wCellSiteRad[(size_t)subB*15 + best]
                               : (ft == 0 ? m_wFaceSiteRadTri[(size_t)best]
                                          : m_wFaceSiteRadQuad[(size_t)subB*7 + best]);
    const double* nd = cell ? &m_wCellNodes[(size_t)hS*54] : &m_wFaceNodes[(size_t)hS*27];
    const int nn = cell ? 18 : (ft == 0 ? 6 : 9);
    const int nm = cell ? 27 : 9;
    const int nS = (int)srcG.size();
    std::vector<int> col((size_t)nS);
    for (int ls = 0; ls < nS; ++ls) {
        const int* e = &m_expo[(size_t)3*srcG[ls]];
        col[ls] = e[0] + 3*e[1] + (cell ? 9*e[2] : 0);
    }
    for (int q = 0; q < R.nq; ++q) {
        const double* Sq = &R.S[(size_t)q*nn];
        double X0 = 0.0, X1 = 0.0, X2 = 0.0;
        for (int n2 = 0; n2 < nn; ++n2) { const double s = Sq[n2]; const double* v = &nd[3*n2]; X0 += s*v[0]; X1 += s*v[1]; X2 += s*v[2]; }
        const double dx = p[0]-X0, dy = p[1]-X1, dz = p[2]-X2;
        const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 1e-300) continue;
        const double g = R.w[q]/r;
        const double* Mq = &R.M[(size_t)q*nm];
        for (int ls = 0; ls < nS; ++ls) inn[ls] += g*Mq[col[ls]];
    }
}

static bool WedgeQuadAffineCoefficients(const double* nodes,const double* velocity,
                                        int eu,int ev,int order,double* coeff,double* derivative);

// Far field point -> cheap cached far cloud; else -> the static-site radial.  Mirror of PhiInnerHexSubVec.
void RadHACApKChargeGram::PhiInnerWedgeSubVec(int kindS, int hS, int subB, const double p[3],
                                              const std::vector<int>& srcG, double* inn) const
{
    const bool cell = (kindS == 0);
    const int ft = cell ? -1 : m_wFaceType[hS];
    // An affine quad has an exact physical-polynomial source potential.  Use
    // the same moment path for non-self pairs as for the self radial entry;
    // this removes an otherwise artificial site/radial accuracy boundary.
    if (!cell && ft == 1) {
        const double* qnd=&m_wFaceNodes[(size_t)hS*27];
        double probe[QUAD_AFFINE_POLY_N]={};
        if (WedgeQuadAffineCoefficients(qnd,nullptr,0,0,1,probe,nullptr)) {
            int order=1;for(int charge:srcG){const int*e=&m_expo[(size_t)3*charge];order=std::max({order,e[0],e[1]});}
            const int np=(2*order+1)*(2*order+2)*(2*order+3)/6;
            const int*tv=QUADREF_TRIS[subB];double V[3][3];
            for(int a=0;a<3;++a){const double uv[2]={QUADREF_V[tv[a]][0],QUADREF_V[tv[a]][1]};QuadQ2MapX(qnd,uv,V[a]);}
            double moments[QUAD_AFFINE_POLY_N]={};
            if(order==1)rad_hdiv::TriPotentialMomentsUpTo2(V,p,moments);else rad_hdiv::TriPotentialMomentsUpTo4(V,p,moments);
            for(int j=0;j<(int)srcG.size();++j){const int*e=&m_expo[(size_t)3*srcG[j]];double coeff[QUAD_AFFINE_POLY_N]={};WedgeQuadAffineCoefficients(qnd,nullptr,e[0],e[1],order,coeff,nullptr);for(int k=0;k<np;++k)inn[j]+=coeff[k]*moments[k];}
            return;
        }
    }
    const size_t sid = cell ? ((size_t)hS*3 + subB) : ((size_t)hS*2 + subB);
    const double* cs = cell ? &m_wCellSubC[sid*3] : &m_wFaceSubC[sid*3];
    const double  sz = cell ? m_wCellSubS[sid] : m_wFaceSubS[sid];
    const double dxc = p[0]-cs[0], dyc = p[1]-cs[1], dzc = p[2]-cs[2];
    const bool far_pt = std::sqrt(dxc*dxc + dyc*dyc + dzc*dzc) > m_far_inner_factor*sz;
    if (!far_pt) { PhiInnerWedgeSiteVec(kindS, hS, subB, p, srcG, inn); return; }
    const double* nd = cell ? &m_wCellNodes[(size_t)hS*54] : &m_wFaceNodes[(size_t)hS*27];
    const std::shared_ptr<const HexQuadCloud> cl =
        HexGetCloud(m_build_id, HexCloudKey(cell ? 0 : 1, false, false, hS, subB, 3),
        [&](HexQuadCloud& c) {
            if (cell) WedgeBuildCloud(nd, 0, -1, subB, m_farTetP.data(), m_farTetW.data(), (int)m_farTetW.size(), false, c);
            else      WedgeBuildCloud(nd, 1, ft, subB, m_farTriP.data(), m_farTriW.data(), (int)m_farTriW.size(), false, c);
        });
    const int nq = (int)cl->wgeo.size();
    const int nS = (int)srcG.size();
    for (int q = 0; q < nq; ++q) {
        const double dx = p[0]-cl->pts[3*q], dy = p[1]-cl->pts[3*q+1], dz = p[2]-cl->pts[3*q+2];
        const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 1e-300) continue;
        const double gr = cl->wgeo[q]/r;
        const double* xi = &cl->xi[3*q];
        for (int ls = 0; ls < nS; ++ls) inn[ls] += gr*HexMonoEval(srcG[ls], xi);
    }
}

static bool WedgeQuadAffineCoefficients(const double* nodes,const double* velocity,
                                        int eu,int ev,int order,double* coeff,double* derivative)
{
    double lin[2][4],invj=0;const int np=(2*order+1)*(2*order+2)*(2*order+3)/6;
    if(!velocity){if(!QuadAffineInverseForms(nodes,lin,invj))return false;double poly[HEX_AFFINE_POLY_N]={};int deg=0;poly[0]=1;for(int q=0;q<eu;++q)HexPolyMulLinear(poly,deg,lin[0],np);for(int q=0;q<ev;++q)HexPolyMulLinear(poly,deg,lin[1],np);for(int k=0;k<np;++k)coeff[k]=invj*poly[k];return true;}
    double dlin[2][4],dinvj=0;if(!QuadAffineInverseFormsDirectional(nodes,velocity,lin,dlin,invj,dinvj))return false;double poly[HEX_AFFINE_POLY_N]={},dpoly[HEX_AFFINE_POLY_N]={};int deg=0;poly[0]=1;for(int q=0;q<eu;++q)HexPolyMulLinearDirectional(poly,dpoly,deg,lin[0],dlin[0],np);for(int q=0;q<ev;++q)HexPolyMulLinearDirectional(poly,dpoly,deg,lin[1],dlin[1],np);for(int k=0;k<np;++k){coeff[k]=invj*poly[k];derivative[k]=dinvj*poly[k]+invj*dpoly[k];}return true;
}

void RadHACApKChargeGram::DPhiInnerWedgeSubVec(int kindS,int hS,int subB,const double p[3],const double dp[3],
                                               const double* velocity,const std::vector<int>& srcG,double* dinn) const
{
    const bool cell=kindS==0;const int ft=cell?-1:m_wFaceType[hS];const size_t sid=cell?((size_t)hS*3+subB):((size_t)hS*2+subB);
    if(!cell&&ft==1){
        const double*qnd=&m_wFaceNodes[(size_t)hS*27];double probe[QUAD_AFFINE_POLY_N]={};
        if(WedgeQuadAffineCoefficients(qnd,nullptr,0,0,1,probe,nullptr)){
            int order=1;for(int charge:srcG){const int*e=&m_expo[(size_t)3*charge];order=std::max({order,e[0],e[1]});}const int np=(2*order+1)*(2*order+2)*(2*order+3)/6;
            const int*tv=QUADREF_TRIS[subB];double V[3][3],dV[3][3];for(int a=0;a<3;++a){const double uv[2]={QUADREF_V[tv[a]][0],QUADREF_V[tv[a]][1]};QuadQ2MapX(qnd,uv,V[a]);QuadQ2MapX(velocity,uv,dV[a]);}
            double mv[QUAD_AFFINE_POLY_N]={},dm[QUAD_AFFINE_POLY_N]={};if(order==1)rad_hdiv::TriPotentialMomentsDirectionalUpTo2(V,dV,p,dp,mv,dm);else rad_hdiv::TriPotentialMomentsDirectionalUpTo4(V,dV,p,dp,mv,dm);
            for(int j=0;j<(int)srcG.size();++j){const int*e=&m_expo[(size_t)3*srcG[j]];double c[QUAD_AFFINE_POLY_N]={},dc[QUAD_AFFINE_POLY_N]={};WedgeQuadAffineCoefficients(qnd,velocity,e[0],e[1],order,c,dc);for(int k=0;k<np;++k)dinn[j]+=dc[k]*mv[k]+c[k]*dm[k];}return;
        }
    }
    const double*cs=cell?&m_wCellSubC[sid*3]:&m_wFaceSubC[sid*3];const double sz=cell?m_wCellSubS[sid]:m_wFaceSubS[sid];const double dx=p[0]-cs[0],dy=p[1]-cs[1],dz=p[2]-cs[2];
    const bool far=std::sqrt(dx*dx+dy*dy+dz*dz)>m_far_inner_factor*sz;const double*nd=cell?&m_wCellNodes[(size_t)hS*54]:&m_wFaceNodes[(size_t)hS*27];
    if(far){HexQuadCloud c;if(cell)WedgeBuildCloud(nd,0,-1,subB,m_farTetP.data(),m_farTetW.data(),(int)m_farTetW.size(),false,c);else WedgeBuildCloud(nd,1,ft,subB,m_farTriP.data(),m_farTriW.data(),(int)m_farTriW.size(),false,c);for(int q=0;q<(int)c.wgeo.size();++q){const double*xi=&c.xi[3*q];double dX[3];if(cell)WedgeQ2MapX(velocity,xi,dX);else{const double uv[2]={xi[0],xi[1]};if(ft==0)TriSurfMap(velocity,uv,dX);else QuadQ2MapX(velocity,uv,dX);}const double R[3]={p[0]-c.pts[3*q],p[1]-c.pts[3*q+1],p[2]-c.pts[3*q+2]},dR[3]={dp[0]-dX[0],dp[1]-dX[1],dp[2]-dX[2]},r2=R[0]*R[0]+R[1]*R[1]+R[2]*R[2];if(r2<1e-300)continue;const double dk=-(R[0]*dR[0]+R[1]*dR[1]+R[2]*dR[2])/(r2*std::sqrt(r2));for(int j=0;j<(int)srcG.size();++j)dinn[j]+=c.wgeo[q]*dk*HexMonoEval(srcG[j],xi);}return;}
    const double*sx=cell?&m_wCellSiteX[(((size_t)hS*3+subB)*15)*3]:&m_wFaceSiteX[(((size_t)hS*2+subB)*7)*3];const int nsite=cell?15:7;int best=0;double bd=1e300;for(int k=0;k<nsite;++k){const double x=p[0]-sx[3*k],y=p[1]-sx[3*k+1],z=p[2]-sx[3*k+2],d=x*x+y*y+z*z;if(d<bd){bd=d;best=k;}}
    const HexSiteRad&R=cell?m_wCellSiteRad[(size_t)subB*15+best]:(ft==0?m_wFaceSiteRadTri[(size_t)best]:m_wFaceSiteRadQuad[(size_t)subB*7+best]);const int nn=cell?18:(ft==0?6:9),nm=cell?27:9;std::vector<int>col(srcG.size());for(int j=0;j<(int)srcG.size();++j){const int*e=&m_expo[(size_t)3*srcG[j]];col[j]=e[0]+3*e[1]+(cell?9*e[2]:0);}for(int q=0;q<R.nq;++q){const double*S=&R.S[(size_t)q*nn];double X[3]={},dX[3]={};for(int i=0;i<nn;++i)for(int k=0;k<3;++k){X[k]+=S[i]*nd[3*i+k];dX[k]+=S[i]*velocity[3*i+k];}const double rv[3]={p[0]-X[0],p[1]-X[1],p[2]-X[2]},dv[3]={dp[0]-dX[0],dp[1]-dX[1],dp[2]-dX[2]},r2=rv[0]*rv[0]+rv[1]*rv[1]+rv[2]*rv[2];if(r2<1e-300)continue;const double dk=-(rv[0]*dv[0]+rv[1]*dv[1]+rv[2]*dv[2])/(r2*std::sqrt(r2));const double*M=&R.M[(size_t)q*nm];for(int j=0;j<(int)srcG.size();++j)dinn[j]+=R.w[q]*dk*M[col[j]];}
}

// SELF inner: the exact-anchor (xiT) REF-frame radial decomposition.  Mirror of PhiInnerHexRadialVec.
void RadHACApKChargeGram::PhiInnerWedgeRadialVec(int kindS, int hS, int subB, const double p[3],
                                                 const double* xiT, const std::vector<int>& srcG,
                                                 double* inn) const
{
    if (!xiT) throw std::logic_error("PhiInnerWedgeRadialVec: xiT required (SELF-only)");
    const bool cell = (kindS == 0);
    const int ft = cell ? -1 : m_wFaceType[hS];
    if(!cell && ft==1){
        const double* qnd=&m_wFaceNodes[(size_t)hS*27];double probe[QUAD_AFFINE_POLY_N]={};
        if(WedgeQuadAffineCoefficients(qnd,nullptr,0,0,1,probe,nullptr)){
            int order=1;for(int charge:srcG){const int*e=&m_expo[(size_t)3*charge];order=std::max({order,e[0],e[1]});}
            const int np=(2*order+1)*(2*order+2)*(2*order+3)/6;
            const int*tv=QUADREF_TRIS[subB];double V[3][3];for(int a=0;a<3;++a){const double uv[2]={QUADREF_V[tv[a]][0],QUADREF_V[tv[a]][1]};QuadQ2MapX(qnd,uv,V[a]);}
            double moments[QUAD_AFFINE_POLY_N]={};if(order==1)rad_hdiv::TriPotentialMomentsUpTo2(V,p,moments);else rad_hdiv::TriPotentialMomentsUpTo4(V,p,moments);
            for(int j=0;j<(int)srcG.size();++j){const int*e=&m_expo[(size_t)3*srcG[j]];double coeff[QUAD_AFFINE_POLY_N]={};WedgeQuadAffineCoefficients(qnd,nullptr,e[0],e[1],order,coeff,nullptr);for(int k=0;k<np;++k)inn[j]+=coeff[k]*moments[k];}
            return;
        }
    }
    const double* nd = cell ? &m_wCellNodes[(size_t)hS*54] : &m_wFaceNodes[(size_t)hS*27];
    const int nR = (int)m_glIn.size();
    const double* GL = m_glIn.data();
    const double* GW = m_gwIn.data();
    const int nS = (int)srcG.size();
    std::vector<double> acc((size_t)nS, 0.0);
    if (cell) {
        const int* tv = WEDGEREF_TETS[subB];
        double V[4][3];
        for (int i = 0; i < 4; ++i) for (int k = 0; k < 3; ++k) V[i][k] = WEDGEREF_V[tv[i]][k];
        const double xr[3] = {xiT[0], xiT[1], xiT[2]};
        double x0[3]; rad_hdiv::ClosestPointTet(V, xr, x0);
        double E0[3], E1[3], E2[3];
        for (int k = 0; k < 3; ++k) { E0[k] = V[1][k]-V[0][k]; E1[k] = V[2][k]-V[0][k]; E2[k] = V[3][k]-V[0][k]; }
        const double hv = E0[0]*(E1[1]*E2[2]-E1[2]*E2[1]) - E0[1]*(E1[0]*E2[2]-E1[2]*E2[0])
                        + E0[2]*(E1[0]*E2[1]-E1[1]*E2[0]);
        const double sgnT = (hv >= 0.0) ? 1.0 : -1.0;
        for (int f = 0; f < 4; ++f) {
            const double* b1 = V[HEXTET_FC[f][0]]; const double* b2 = V[HEXTET_FC[f][1]]; const double* b3 = V[HEXTET_FC[f][2]];
            double d1[3], d2[3], d3[3], e21[3], e32[3];
            for (int k = 0; k < 3; ++k) { d1[k] = b1[k]-x0[k]; d2[k] = b2[k]-x0[k]; d3[k] = b3[k]-x0[k]; e21[k] = b2[k]-b1[k]; e32[k] = b3[k]-b2[k]; }
            const double cr[3] = {d2[1]*d3[2]-d2[2]*d3[1], d2[2]*d3[0]-d2[0]*d3[2], d2[0]*d3[1]-d2[1]*d3[0]};
            const double D = d1[0]*cr[0] + d1[1]*cr[1] + d1[2]*cr[2];
            if (std::fabs(D) < 1e-300) continue;
            for (int a = 0; a < nR; ++a) { const double u = GL[a];
                for (int b = 0; b < nR; ++b) { const double v = GL[b];
                    for (int c = 0; c < nR; ++c) { const double w = GL[c];
                        double y[3];
                        for (int k = 0; k < 3; ++k) y[k] = x0[k] + u*(d1[k] + v*(e21[k] + w*e32[k]));
                        double X[3]; WedgeQ2MapX(nd, y, X);
                        const double dx = p[0]-X[0], dy = p[1]-X[1], dz = p[2]-X[2];
                        const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                        if (r < 1e-300) continue;
                        const double wq = GW[a]*GW[b]*GW[c]*(u*u*v*D)/r;
                        for (int ls = 0; ls < nS; ++ls) acc[ls] += wq*HexMonoEval(srcG[ls], y);
                    }
                }
            }
        }
        for (int ls = 0; ls < nS; ++ls) inn[ls] += sgnT*acc[ls];
    } else {
        double V2[3][2]; WFaceSubTriRef(ft, subB, V2);
        const double ur[2] = {xiT[0], xiT[1]};
        double u0[2]; ClosestPointTri2D(V2, ur, u0);
        for (int kf = 0; kf < 3; ++kf) {
            const double* A = V2[kf]; const double* B = V2[(kf+1)%3];
            const double ea[2] = {A[0]-u0[0], A[1]-u0[1]};
            const double eb[2] = {B[0]-u0[0], B[1]-u0[1]};
            const double s2 = ea[0]*eb[1] - ea[1]*eb[0];
            if (std::fabs(s2) < 1e-300) continue;
            for (int a = 0; a < nR; ++a) { const double u = GL[a];
                for (int b = 0; b < nR; ++b) { const double v = GL[b];
                    const double yuv[2] = {u0[0] + u*(ea[0] + v*(eb[0]-ea[0])), u0[1] + u*(ea[1] + v*(eb[1]-ea[1]))};
                    double X[3];
                    if (ft == 0) TriSurfMap(nd, yuv, X); else QuadQ2MapX(nd, yuv, X);
                    const double dx = p[0]-X[0], dy = p[1]-X[1], dz = p[2]-X[2];
                    const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                    if (r < 1e-300) continue;
                    const double wq = GW[a]*GW[b]*(u*s2)/r;
                    const double y3[3] = {yuv[0], yuv[1], 0.0};
                    for (int ls = 0; ls < nS; ++ls) acc[ls] += wq*HexMonoEval(srcG[ls], y3);
                }
            }
        }
        for (int ls = 0; ls < nS; ++ls) inn[ls] += acc[ls];
    }
}

void RadHACApKChargeGram::DPhiInnerWedgeRadialVec(
    int kind, int host, int sub, const double p[3], const double dp[3], const double* xiT,
    const double* velocity, const std::vector<int>& srcG, double* dinn) const
{
    if(!xiT) throw std::logic_error("DPhiInnerWedgeRadialVec: xiT required");
    const bool cell=kind==0; const int ft=cell?-1:m_wFaceType[host];
    const double* nd=cell?&m_wCellNodes[(size_t)host*54]:&m_wFaceNodes[(size_t)host*27];
    const int nr=(int)m_glIn.size(),ns=(int)srcG.size(); std::vector<double> acc(ns);
    auto add=[&](const double X[3],const double dX[3],double w,const double y[3]){
        const double R[3]={p[0]-X[0],p[1]-X[1],p[2]-X[2]},dR[3]={dp[0]-dX[0],dp[1]-dX[1],dp[2]-dX[2]};
        const double r2=R[0]*R[0]+R[1]*R[1]+R[2]*R[2];if(r2<1e-300)return;
        const double dk=-(R[0]*dR[0]+R[1]*dR[1]+R[2]*dR[2])/(r2*std::sqrt(r2));
        for(int j=0;j<ns;++j)acc[j]+=w*dk*HexMonoEval(srcG[j],y);
    };
    if(cell){const int*tv=WEDGEREF_TETS[sub];double V[4][3];for(int i=0;i<4;++i)for(int k=0;k<3;++k)V[i][k]=WEDGEREF_V[tv[i]][k];double x0[3];rad_hdiv::ClosestPointTet(V,xiT,x0);
        double E0[3],E1[3],E2[3];for(int k=0;k<3;++k){E0[k]=V[1][k]-V[0][k];E1[k]=V[2][k]-V[0][k];E2[k]=V[3][k]-V[0][k];}
        const double hv=E0[0]*(E1[1]*E2[2]-E1[2]*E2[1])-E0[1]*(E1[0]*E2[2]-E1[2]*E2[0])+E0[2]*(E1[0]*E2[1]-E1[1]*E2[0]),sgn=hv>=0?1.0:-1.0;
        for(int f=0;f<4;++f){const double*b1=V[HEXTET_FC[f][0]],*b2=V[HEXTET_FC[f][1]],*b3=V[HEXTET_FC[f][2]];double d1[3],d2[3],d3[3],e21[3],e32[3];for(int k=0;k<3;++k){d1[k]=b1[k]-x0[k];d2[k]=b2[k]-x0[k];d3[k]=b3[k]-x0[k];e21[k]=b2[k]-b1[k];e32[k]=b3[k]-b2[k];}const double cr[3]={d2[1]*d3[2]-d2[2]*d3[1],d2[2]*d3[0]-d2[0]*d3[2],d2[0]*d3[1]-d2[1]*d3[0]},D=d1[0]*cr[0]+d1[1]*cr[1]+d1[2]*cr[2];if(std::abs(D)<1e-300)continue;
            for(int a=0;a<nr;++a){const double u=m_glIn[a];for(int b=0;b<nr;++b){const double v=m_glIn[b];for(int c=0;c<nr;++c){const double w=m_glIn[c];double y[3];for(int k=0;k<3;++k)y[k]=x0[k]+u*(d1[k]+v*(e21[k]+w*e32[k]));double X[3],dX[3];WedgeQ2MapX(nd,y,X);WedgeQ2MapX(velocity,y,dX);add(X,dX,sgn*m_gwIn[a]*m_gwIn[b]*m_gwIn[c]*u*u*v*D,y);}}}
        }
    }else{double V[3][2];WFaceSubTriRef(ft,sub,V);const double ur[2]={xiT[0],xiT[1]};double x0[2];ClosestPointTri2D(V,ur,x0);for(int f=0;f<3;++f){const double*A=V[f],*B=V[(f+1)%3],ea[2]={A[0]-x0[0],A[1]-x0[1]},eb[2]={B[0]-x0[0],B[1]-x0[1]},s2=ea[0]*eb[1]-ea[1]*eb[0];if(std::abs(s2)<1e-300)continue;for(int a=0;a<nr;++a){const double u=m_glIn[a];for(int b=0;b<nr;++b){const double v=m_glIn[b],y2[2]={x0[0]+u*(ea[0]+v*(eb[0]-ea[0])),x0[1]+u*(ea[1]+v*(eb[1]-ea[1]))};double X[3],dX[3];if(ft==0){TriSurfMap(nd,y2,X);TriSurfMap(velocity,y2,dX);}else{QuadQ2MapX(nd,y2,X);QuadQ2MapX(velocity,y2,dX);}const double y[3]={y2[0],y2[1],0};add(X,dX,m_gwIn[a]*m_gwIn[b]*u*s2,y);}}}}
    for(int j=0;j<ns;++j)dinn[j]+=acc[j];
}

// Directed host-pair block (mirror of QuadBlockHex) with mixed-face sub counts / node strides.
std::vector<double> RadHACApKChargeGram::QuadBlockWedge(int kindT, int hT, int kindS, int hS, int img) const
{
    const std::vector<int>& tgtG = (kindT == 0) ? m_cellCharges[hT] : m_faceCharges[hT];
    const std::vector<int>& srcG = (kindS == 0) ? m_cellCharges[hS] : m_faceCharges[hS];
    const int nT = (int)tgtG.size(), nS = (int)srcG.size();
    std::vector<double> blk((size_t)nT*nS, 0.0);
    if (nT == 0 || nS == 0) return blk;
    const bool cellT = (kindT == 0), cellS = (kindS == 0);
    const int ftT = cellT ? -1 : m_wFaceType[hT];
    const int ftS = cellS ? -1 : m_wFaceType[hS];
    const int nsubT = cellT ? 3 : (ftT == 0 ? 1 : 2);
    const int nsubS = cellS ? 3 : (ftS == 0 ? 1 : 2);
    // IMA (img>0): see QuadBlockHex -- couple the TARGET host with the image-mapped source host (map the
    // grading geometry + the outer eval point, never SELF).  img==0 => identity => byte-identical direct block.
    auto reflpt = [this, img](const double* v, double* o){ ImageEvalPoint(img, v, o); };
    const int rt = tgtG[0], rs = srcG[0];
    double rsc[3]; reflpt(&m_cent[3*rs], rsc);
    const double dxh = m_cent[3*rt]-rsc[0], dyh = m_cent[3*rt+1]-rsc[1], dzh = m_cent[3*rt+2]-rsc[2];
    const double r_h = std::sqrt(dxh*dxh + dyh*dyh + dzh*dzh);
    // self_pair: the (reflected) source host coincides with the target host -> use the EXACT self-radial.
    // img==0 is always self for hT==hS.  img>0 is self ONLY when the host is INVARIANT under the transform
    // (lies ON the mirror plane, T(host)==host -- e.g. a z=0 cut FACE reflected across z=0 is itself; under a
    // cyclic ROTATION only an on-axis host qualifies).  Then
    // the reflected self-term EXACTLY equals the direct self-term (same exact radial quadrature), so the LARGE
    // on-plane cut-face charge (sigma = M.n ~ |M| when M is perpendicular to the plane) CANCELS exactly for
    // sign -1, instead of leaving a quadrature-mismatch residual (the ~1.5% hex / ~29% wedge antisymmetric-
    // plane error, 2026-07-05: direct self used the exact radial, reflected self used the ~1e-3 site-radial).
    // An OFF-plane self host (a z>0 cell) reflects to a genuine image elsewhere -> NOT self.
    bool self_pair = (kindT == kindS && hT == hS);
    if (self_pair && img != 0) {
        double rc_[3]; reflpt(&m_cent[3*rt], rc_);
        const double d_ = std::abs(rc_[0]-m_cent[3*rt]) + std::abs(rc_[1]-m_cent[3*rt+1]) + std::abs(rc_[2]-m_cent[3*rt+2]);
        self_pair = (d_ < 1e-6 * m_size[rt] + 1e-12);   // T(host)==host <=> host is invariant
    }
    const bool near_hosts = self_pair || r_h <= m_near_grade*(m_size[rt] + m_size[rs]);
    const double* ndT = cellT ? &m_wCellNodes[(size_t)hT*54] : &m_wFaceNodes[(size_t)hT*27];
    const int nvT = cellT ? 4 : 3;
    std::vector<double> inn(nS), owt(nT);
    for (int sA = 0; sA < nsubT; ++sA) {
        const size_t sidA = cellT ? ((size_t)hT*3 + sA) : ((size_t)hT*2 + sA);
        const double szA = cellT ? m_wCellSubS[sidA] : m_wFaceSubS[sidA];
        const double* subVA = cellT ? &m_wCellSubV[sidA*4*3] : &m_wFaceSubV[sidA*3*3];
        const double* cA = cellT ? &m_wCellSubC[sidA*3] : &m_wFaceSubC[sidA*3];
        for (int sB = 0; sB < nsubS; ++sB) {
            const size_t sidB = cellS ? ((size_t)hS*3 + sB) : ((size_t)hS*2 + sB);
            const double* cB0 = cellS ? &m_wCellSubC[sidB*3] : &m_wFaceSubC[sidB*3];
            double cB[3]; reflpt(cB0, cB);   // mapped source sub-centroid (img>0)
            const double szB = cellS ? m_wCellSubS[sidB] : m_wFaceSubS[sidB];
            const double dx = cA[0]-cB[0], dy = cA[1]-cB[1], dz = cA[2]-cB[2];
            const bool near_sub = near_hosts && std::sqrt(dx*dx + dy*dy + dz*dz) <= m_near_grade*(szA + szB);
            std::shared_ptr<const HexQuadCloud> oc;
            if (!near_sub) {
                const int nqreg = cellT ? (int)m_symTetW.size() : (int)m_symTriW.size();
                oc = HexGetCloud(m_build_id, HexCloudKey(cellT ? 0 : 1, true, false, hT, sA, 3),
                    [&](HexQuadCloud& c) {
                        if (cellT) WedgeBuildCloud(ndT, 0, -1, sA, m_symTetP.data(), m_symTetW.data(), nqreg, false, c);
                        else       WedgeBuildCloud(ndT, 1, ftT, sA, m_symTriP.data(), m_symTriW.data(), nqreg, false, c);
                    });
            } else {
                int corner = 0; double best = 1e300;
                const double corner_tie_tol = 1e-13 * (szA + szB + 1.0) * (szA + szB + 1.0);
                for (int i = 0; i < nvT; ++i) {
                    const double ddx = subVA[3*i]-cB[0], ddy = subVA[3*i+1]-cB[1], ddz = subVA[3*i+2]-cB[2];
                    const double d = ddx*ddx + ddy*ddy + ddz*ddz;
                    if (d < best - corner_tie_tol) { best = d; corner = i; }
                }
                oc = HexGetCloud(m_build_id, HexCloudKey(cellT ? 0 : 1, true, true, hT, sA, corner),
                    [&](HexQuadCloud& c) {
                        std::vector<double> gb, gw;
                        HexDuffyBary(cellT ? 3 : 2, corner, m_glOut, m_gwOut, gb, gw);
                        WedgeBuildCloud(ndT, cellT ? 0 : 1, ftT, sA, gb.data(), gw.data(), (int)gw.size(), true, c);
                    });
            }
            const int nqo = (int)oc->wgeo.size();
            for (int q = 0; q < nqo; ++q) {
                const double pq[3] = {oc->pts[3*q], oc->pts[3*q+1], oc->pts[3*q+2]};
                const double* xiT = &oc->xi[3*q];
                for (int ls = 0; ls < nS; ++ls) inn[ls] = 0.0;
                double peval[3]; reflpt(pq, peval);                                             // reflect for the (image) source eval
                if (self_pair) PhiInnerWedgeRadialVec(kindS, hS, sB, pq, xiT, srcG, inn.data());
                else           PhiInnerWedgeSubVec(kindS, hS, sB, peval, srcG, inn.data());      // peval==pq if mask==0
                const double wg = oc->wgeo[q];
                for (int lt = 0; lt < nT; ++lt) owt[lt] = wg*HexMonoEval(tgtG[lt], xiT);
                for (int lt = 0; lt < nT; ++lt) {
                    const double wl = owt[lt];
                    double* row = &blk[(size_t)lt*nS];
                    for (int ls = 0; ls < nS; ++ls) row[ls] += wl*inn[ls];
                }
            }
        }
    }
    for (double& v : blk) v *= RAD_INV_FOUR_PI;
    return blk;
}

static void FinishWedgeSelfDerivative(std::vector<double>& a,int n){for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){const double v=.5*(a[(size_t)i*n+j]+a[(size_t)j*n+i]);a[(size_t)i*n+j]=a[(size_t)j*n+i]=v;}for(double&v:a)v*=RAD_INV_FOUR_PI;}

std::vector<double> RadHACApKChargeGram::WedgeVolumeSelfBlockDirectionalDerivative(int host,const std::vector<double>& vel) const
{
    if(!m_wedgemode)throw std::logic_error("WedgeVolumeSelfBlockDirectionalDerivative requires a WEDGE charge Gram");if(host<0||host>=(int)m_cellCharges.size())throw std::out_of_range("WEDGE host out of range");if(vel.size()!=54)throw std::invalid_argument("node_velocity must have shape (18,3)");
    const auto&g=m_cellCharges[host];const int n=(int)g.size();std::vector<double>out((size_t)n*n),di(n);const double*nd=&m_wCellNodes[(size_t)host*54];
    for(int sa=0;sa<3;++sa){const size_t ia=(size_t)host*3+sa;const double*va=&m_wCellSubV[ia*12],*ca=&m_wCellSubC[ia*3];for(int sb=0;sb<3;++sb){const size_t ib=(size_t)host*3+sb;const double*cb=&m_wCellSubC[ib*3];int corner=0;double best=1e300;for(int i=0;i<4;++i){const double x=va[3*i]-cb[0],y=va[3*i+1]-cb[1],z=va[3*i+2]-cb[2],d=x*x+y*y+z*z;if(d<best){best=d;corner=i;}}const double x=ca[0]-cb[0],y=ca[1]-cb[1],z=ca[2]-cb[2];const bool near=std::sqrt(x*x+y*y+z*z)<=m_near_grade*(m_wCellSubS[ia]+m_wCellSubS[ib]);std::shared_ptr<const HexQuadCloud>oc;
        if(near)oc=HexGetCloud(m_build_id,HexCloudKey(0,true,true,host,sa,corner),[&](HexQuadCloud&c){std::vector<double>b,w;HexDuffyBary(3,corner,m_glOut,m_gwOut,b,w);WedgeBuildCloud(nd,0,-1,sa,b.data(),w.data(),(int)w.size(),true,c);});else oc=HexGetCloud(m_build_id,HexCloudKey(0,true,false,host,sa,3),[&](HexQuadCloud&c){WedgeBuildCloud(nd,0,-1,sa,m_symTetP.data(),m_symTetW.data(),(int)m_symTetW.size(),false,c);});
        for(int q=0;q<(int)oc->wgeo.size();++q){const double*xi=&oc->xi[3*q],p[3]={oc->pts[3*q],oc->pts[3*q+1],oc->pts[3*q+2]};double dp[3];WedgeQ2MapX(vel.data(),xi,dp);std::fill(di.begin(),di.end(),0);DPhiInnerWedgeRadialVec(0,host,sb,p,dp,xi,vel.data(),g,di.data());for(int i=0;i<n;++i){const double w=oc->wgeo[q]*HexMonoEval(g[i],xi);for(int j=0;j<n;++j)out[(size_t)i*n+j]+=w*di[j];}}
    }}FinishWedgeSelfDerivative(out,n);return out;
}

std::vector<double> RadHACApKChargeGram::WedgeFaceSelfBlockDirectionalDerivative(int host,const std::vector<double>&vel) const
{
    if(!m_wedgemode)throw std::logic_error("WedgeFaceSelfBlockDirectionalDerivative requires a WEDGE charge Gram");if(host<0||host>=(int)m_faceCharges.size())throw std::out_of_range("WEDGE face host out of range");const int ft=m_wFaceType[host],nn=ft==0?6:9;if(vel.size()!=(size_t)3*nn)throw std::invalid_argument(ft==0?"node_velocity must have shape (6,3)":"node_velocity must have shape (9,3)");
    const auto&g=m_faceCharges[host];const int n=(int)g.size(),ns=ft==0?1:2;std::vector<double>out((size_t)n*n),di(n);const double*nd=&m_wFaceNodes[(size_t)host*27];
    double affine_probe[QUAD_AFFINE_POLY_N]={};
    const bool affine_quad=ft==1&&WedgeQuadAffineCoefficients(nd,nullptr,0,0,1,affine_probe,nullptr);
    if(affine_quad){
        int order=1;for(int charge:g){const int*e=&m_expo[(size_t)3*charge];order=std::max({order,e[0],e[1]});}const int np=(2*order+1)*(2*order+2)*(2*order+3)/6;std::vector<double>coeff((size_t)n*np),dc((size_t)n*np),inner(n);
        for(int j=0;j<n;++j){const int*e=&m_expo[(size_t)3*g[j]];if(!WedgeQuadAffineCoefficients(nd,vel.data(),e[0],e[1],order,&coeff[(size_t)j*np],&dc[(size_t)j*np]))throw std::logic_error("affine WEDGE quad face derivative has a singular geometry map");}
        for(int sa=0;sa<2;++sa){const size_t ia=(size_t)host*2+sa;const double*va=&m_wFaceSubV[ia*9],*ca=&m_wFaceSubC[ia*3];for(int sb=0;sb<2;++sb){const size_t ib=(size_t)host*2+sb;const double*cb=&m_wFaceSubC[ib*3];const double tie=1e-13*(m_wFaceSubS[ia]+m_wFaceSubS[ib]+1.0)*(m_wFaceSubS[ia]+m_wFaceSubS[ib]+1.0);int corner=0;double best=1e300;for(int i=0;i<3;++i){const double x=va[3*i]-cb[0],y=va[3*i+1]-cb[1],z=va[3*i+2]-cb[2],d=x*x+y*y+z*z;if(d<best-tie){best=d;corner=i;}}const double x=ca[0]-cb[0],y=ca[1]-cb[1],z=ca[2]-cb[2];const bool near=std::sqrt(x*x+y*y+z*z)<=m_near_grade*(m_wFaceSubS[ia]+m_wFaceSubS[ib]);std::shared_ptr<const HexQuadCloud>oc;
            if(near)oc=HexGetCloud(m_build_id,HexCloudKey(1,true,true,host,sa,corner),[&](HexQuadCloud&c){std::vector<double>b,w;HexDuffyBary(2,corner,m_glOut,m_gwOut,b,w);WedgeBuildCloud(nd,1,1,sa,b.data(),w.data(),(int)w.size(),true,c);});else oc=HexGetCloud(m_build_id,HexCloudKey(1,true,false,host,sa,3),[&](HexQuadCloud&c){WedgeBuildCloud(nd,1,1,sa,m_symTriP.data(),m_symTriW.data(),(int)m_symTriW.size(),false,c);});
            const int*tv=QUADREF_TRIS[sb];double V[3][3],dV[3][3];for(int a=0;a<3;++a){const double r[2]={QUADREF_V[tv[a]][0],QUADREF_V[tv[a]][1]};QuadQ2MapX(nd,r,V[a]);QuadQ2MapX(vel.data(),r,dV[a]);}
            for(int q=0;q<(int)oc->wgeo.size();++q){const double*xi=&oc->xi[3*q],p[3]={oc->pts[3*q],oc->pts[3*q+1],oc->pts[3*q+2]},uv[2]={xi[0],xi[1]};double dp[3];QuadQ2MapX(vel.data(),uv,dp);double mv[QUAD_AFFINE_POLY_N]={},dm[QUAD_AFFINE_POLY_N]={};if(order==1)rad_hdiv::TriPotentialMomentsDirectionalUpTo2(V,dV,p,dp,mv,dm);else rad_hdiv::TriPotentialMomentsDirectionalUpTo4(V,dV,p,dp,mv,dm);std::fill(inner.begin(),inner.end(),0);for(int j=0;j<n;++j)for(int k=0;k<np;++k)inner[j]+=dc[(size_t)j*np+k]*mv[k]+coeff[(size_t)j*np+k]*dm[k];for(int i=0;i<n;++i){const double w=oc->wgeo[q]*HexMonoEval(g[i],xi);for(int j=0;j<n;++j)out[(size_t)i*n+j]+=w*inner[j];}}
        }}FinishWedgeSelfDerivative(out,n);return out;
    }
    for(int sa=0;sa<ns;++sa){const size_t ia=(size_t)host*2+sa;const double*va=&m_wFaceSubV[ia*9],*ca=&m_wFaceSubC[ia*3];for(int sb=0;sb<ns;++sb){const size_t ib=(size_t)host*2+sb;const double*cb=&m_wFaceSubC[ib*3];int corner=0;double best=1e300;for(int i=0;i<3;++i){const double x=va[3*i]-cb[0],y=va[3*i+1]-cb[1],z=va[3*i+2]-cb[2],d=x*x+y*y+z*z;if(d<best){best=d;corner=i;}}const double x=ca[0]-cb[0],y=ca[1]-cb[1],z=ca[2]-cb[2];const bool near=std::sqrt(x*x+y*y+z*z)<=m_near_grade*(m_wFaceSubS[ia]+m_wFaceSubS[ib]);std::shared_ptr<const HexQuadCloud>oc;
        if(near)oc=HexGetCloud(m_build_id,HexCloudKey(1,true,true,host,sa,corner),[&](HexQuadCloud&c){std::vector<double>b,w;HexDuffyBary(2,corner,m_glOut,m_gwOut,b,w);WedgeBuildCloud(nd,1,ft,sa,b.data(),w.data(),(int)w.size(),true,c);});else oc=HexGetCloud(m_build_id,HexCloudKey(1,true,false,host,sa,3),[&](HexQuadCloud&c){WedgeBuildCloud(nd,1,ft,sa,m_symTriP.data(),m_symTriW.data(),(int)m_symTriW.size(),false,c);});
        for(int q=0;q<(int)oc->wgeo.size();++q){const double*xi=&oc->xi[3*q],p[3]={oc->pts[3*q],oc->pts[3*q+1],oc->pts[3*q+2]},uv[2]={xi[0],xi[1]};double dp[3];if(ft==0)TriSurfMap(vel.data(),uv,dp);else QuadQ2MapX(vel.data(),uv,dp);std::fill(di.begin(),di.end(),0);DPhiInnerWedgeRadialVec(1,host,sb,p,dp,xi,vel.data(),g,di.data());for(int i=0;i<n;++i){const double w=oc->wgeo[q]*HexMonoEval(g[i],xi);for(int j=0;j<n;++j)out[(size_t)i*n+j]+=w*di[j];}}
    }}FinishWedgeSelfDerivative(out,n);return out;
}

std::vector<double> RadHACApKChargeGram::QuadBlockWedgeDirectionalDerivative(
    int kindT,int hT,int kindS,int hS,const double* velocityT,const double* velocityS,int mask) const
{
    if(mask!=0)throw std::invalid_argument("WEDGE shape derivative currently requires direct mask 0");
    const auto&tg=kindT==0?m_cellCharges[hT]:m_faceCharges[hT];
    const auto&sg=kindS==0?m_cellCharges[hS]:m_faceCharges[hS];
    const int nt=(int)tg.size(),ns=(int)sg.size();
    if(kindT==kindS&&hT==hS){
        if(kindT==0)return WedgeVolumeSelfBlockDirectionalDerivative(hT,std::vector<double>(velocityT,velocityT+54));
        const int nn=m_wFaceType[hT]==0?18:27;
        return WedgeFaceSelfBlockDirectionalDerivative(hT,std::vector<double>(velocityT,velocityT+nn));
    }
    std::vector<double>out((size_t)nt*ns,0.0),inn((size_t)ns);
    const bool cellT=kindT==0,cellS=kindS==0;
    const int ftT=cellT?-1:m_wFaceType[hT],ftS=cellS?-1:m_wFaceType[hS];
    const int nsubT=cellT?3:(ftT==0?1:2),nsubS=cellS?3:(ftS==0?1:2);
    const double*ndT=cellT?&m_wCellNodes[(size_t)hT*54]:&m_wFaceNodes[(size_t)hT*27];
    const int nvT=cellT?4:3;
    const int rt=tg[0],rs=sg[0];
    const double dh[3]={m_cent[3*rt]-m_cent[3*rs],m_cent[3*rt+1]-m_cent[3*rs+1],m_cent[3*rt+2]-m_cent[3*rs+2]};
    const bool nearHosts=std::sqrt(dh[0]*dh[0]+dh[1]*dh[1]+dh[2]*dh[2])<=m_near_grade*(m_size[rt]+m_size[rs]);
    for(int sa=0;sa<nsubT;++sa){
        const size_t ia=cellT?(size_t)hT*3+sa:(size_t)hT*2+sa;
        const double*va=cellT?&m_wCellSubV[ia*12]:&m_wFaceSubV[ia*9];
        const double*ca=cellT?&m_wCellSubC[ia*3]:&m_wFaceSubC[ia*3];
        const double sza=cellT?m_wCellSubS[ia]:m_wFaceSubS[ia];
        for(int sb=0;sb<nsubS;++sb){
            const size_t ib=cellS?(size_t)hS*3+sb:(size_t)hS*2+sb;
            const double*cb=cellS?&m_wCellSubC[ib*3]:&m_wFaceSubC[ib*3];
            const double szb=cellS?m_wCellSubS[ib]:m_wFaceSubS[ib];
            const double dx=ca[0]-cb[0],dy=ca[1]-cb[1],dz=ca[2]-cb[2];
            const bool near=nearHosts&&std::sqrt(dx*dx+dy*dy+dz*dz)<=m_near_grade*(sza+szb);
            std::shared_ptr<const HexQuadCloud>oc;
            if(!near){
                oc=HexGetCloud(m_build_id,HexCloudKey(cellT?0:1,true,false,hT,sa,3),[&](HexQuadCloud&c){
                    if(cellT)WedgeBuildCloud(ndT,0,-1,sa,m_symTetP.data(),m_symTetW.data(),(int)m_symTetW.size(),false,c);
                    else WedgeBuildCloud(ndT,1,ftT,sa,m_symTriP.data(),m_symTriW.data(),(int)m_symTriW.size(),false,c);
                });
            }else{
                int corner=0;double best=1e300;
                const double tie=1e-13*(sza+szb+1.0)*(sza+szb+1.0);
                for(int i=0;i<nvT;++i){const double x=va[3*i]-cb[0],y=va[3*i+1]-cb[1],z=va[3*i+2]-cb[2],d=x*x+y*y+z*z;if(d<best-tie){best=d;corner=i;}}
                oc=HexGetCloud(m_build_id,HexCloudKey(cellT?0:1,true,true,hT,sa,corner),[&](HexQuadCloud&c){
                    std::vector<double>b,w;HexDuffyBary(cellT?3:2,corner,m_glOut,m_gwOut,b,w);
                    WedgeBuildCloud(ndT,cellT?0:1,ftT,sa,b.data(),w.data(),(int)w.size(),true,c);
                });
            }
            for(int q=0;q<(int)oc->wgeo.size();++q){
                const double*xi=&oc->xi[3*q];double dp[3];
                if(cellT)WedgeQ2MapX(velocityT,xi,dp);else{const double uv[2]={xi[0],xi[1]};if(ftT==0)TriSurfMap(velocityT,uv,dp);else QuadQ2MapX(velocityT,uv,dp);}
                std::fill(inn.begin(),inn.end(),0.0);
                const double p[3]={oc->pts[3*q],oc->pts[3*q+1],oc->pts[3*q+2]};
                DPhiInnerWedgeSubVec(kindS,hS,sb,p,dp,velocityS,sg,inn.data());
                for(int i=0;i<nt;++i){const double w=oc->wgeo[q]*HexMonoEval(tg[i],xi);for(int j=0;j<ns;++j)out[(size_t)i*ns+j]+=w*inn[j];}
            }
        }
    }
    for(double&x:out)x*=RAD_INV_FOUR_PI;
    return out;
}

std::vector<double> RadHACApKChargeGram::WedgeChargeGramDirectionalDerivative(
    const std::vector<double>&cellVelocity,const std::vector<double>&faceVelocity) const
{
    if(!m_wedgemode||m_d2)throw std::logic_error("WEDGE ChargeGram derivative requires a 3D WEDGE Gram");
    if(cellVelocity.size()!=m_cellCharges.size()*54)throw std::invalid_argument("cell_node_velocity must have shape (ncell,18,3)");
    if(faceVelocity.size()!=m_faceCharges.size()*27)throw std::invalid_argument("face_node_velocity must have shape (nface,9,3), with triangular faces padded to nine nodes");
    std::vector<double>dense((size_t)m_n*m_n,0.0);
    const int nc=(int)m_cellCharges.size(),nf=(int)m_faceCharges.size(),nh=nc+nf;
    auto kind=[&](int h){return h<nc?0:1;};auto host=[&](int h){return h<nc?h:h-nc;};
    auto vel=[&](int k,int h){return k==0?&cellVelocity[(size_t)h*54]:&faceVelocity[(size_t)h*27];};
    auto grp=[&](int k,int h)->const std::vector<int>&{return k==0?m_cellCharges[h]:m_faceCharges[h];};
    for(int A=0;A<nh;++A){const int ka=kind(A),ha=host(A);for(int B=A;B<nh;++B){const int kb=kind(B),hb=host(B);const auto&ga=grp(ka,ha);const auto&gb=grp(kb,hb);
        const auto ab=QuadBlockWedgeDirectionalDerivative(ka,ha,kb,hb,vel(ka,ha),vel(kb,hb));
        if(A==B){for(int i=0;i<(int)ga.size();++i)for(int j=0;j<(int)ga.size();++j)dense[(size_t)ga[i]*m_n+ga[j]]=ab[(size_t)i*ga.size()+j];continue;}
        const auto ba=QuadBlockWedgeDirectionalDerivative(kb,hb,ka,ha,vel(kb,hb),vel(ka,ha));
        for(int i=0;i<(int)ga.size();++i)for(int j=0;j<(int)gb.size();++j){const double x=.5*(ab[(size_t)i*gb.size()+j]+ba[(size_t)j*ga.size()+i]);dense[(size_t)ga[i]*m_n+gb[j]]=x;dense[(size_t)gb[j]*m_n+ga[i]]=x;}
    }}
    return dense;
}

namespace {
struct WedgeDerivativeCacheOwner {
    long long build=-1;
    unsigned long long token=0;
    const double* cells=nullptr;
    const double* faces=nullptr;
    bool operator==(const WedgeDerivativeCacheOwner&o)const{return build==o.build&&token==o.token&&cells==o.cells&&faces==o.faces;}
};
static thread_local WedgeDerivativeCacheOwner s_wedge_derivative_owner;
static thread_local std::unordered_map<HexBlockKey,std::vector<double>,HexBlockKeyHash> s_wedge_derivative_blocks;
static thread_local WedgeDerivativeCacheOwner s_tet_derivative_owner;
static thread_local std::unordered_map<HexBlockKey,std::vector<double>,HexBlockKeyHash> s_tet_derivative_blocks;
}

const std::vector<double>& RadHACApKChargeGram::GetWedgeDirectionalSymBlock(
    int kindA,int hostA,int kindB,int hostB,
    const std::vector<double>&cellVelocity,const std::vector<double>&faceVelocity,
    unsigned long long cacheToken) const
{
    const WedgeDerivativeCacheOwner owner{m_build_id,cacheToken,cellVelocity.data(),faceVelocity.data()};
    if(!(s_wedge_derivative_owner==owner)){s_wedge_derivative_blocks.clear();s_wedge_derivative_owner=owner;}
    const HexBlockKey key{kindA,hostA,kindB,hostB,0};
    auto it=s_wedge_derivative_blocks.find(key);
    if(it!=s_wedge_derivative_blocks.end())return it->second;
    if(s_wedge_derivative_blocks.size()>HexBlockCacheLimit())s_wedge_derivative_blocks.clear();
    const auto&ga=kindA==0?m_cellCharges[hostA]:m_faceCharges[hostA];
    const auto&gb=kindB==0?m_cellCharges[hostB]:m_faceCharges[hostB];
    auto vel=[&](int kind,int host){return kind==0?&cellVelocity[(size_t)host*54]:&faceVelocity[(size_t)host*27];};
    auto ab=QuadBlockWedgeDirectionalDerivative(kindA,hostA,kindB,hostB,vel(kindA,hostA),vel(kindB,hostB));
    std::vector<double>sym((size_t)ga.size()*gb.size());
    if(kindA==kindB&&hostA==hostB) sym=std::move(ab);
    else{
        auto ba=QuadBlockWedgeDirectionalDerivative(kindB,hostB,kindA,hostA,vel(kindB,hostB),vel(kindA,hostA));
        for(size_t i=0;i<ga.size();++i)for(size_t j=0;j<gb.size();++j)sym[i*gb.size()+j]=.5*(ab[i*gb.size()+j]+ba[j*ga.size()+i]);
    }
    return s_wedge_derivative_blocks.emplace(key,std::move(sym)).first->second;
}

double RadHACApKChargeGram::WedgeChargeGramDirectionalDerivativeElement(
    int row,int col,const std::vector<double>&cellVelocity,const std::vector<double>&faceVelocity,
    unsigned long long cacheToken) const
{
    if(!m_wedgemode||m_d2)throw std::logic_error("WEDGE derivative entry requires a 3D WEDGE Gram");
    if(row<0||row>=m_n||col<0||col>=m_n)throw std::out_of_range("WEDGE derivative entry index out of range");
    if(cellVelocity.size()!=m_cellCharges.size()*54)throw std::invalid_argument("cell_node_velocity must have shape (ncell,18,3)");
    if(faceVelocity.size()!=m_faceCharges.size()*27)throw std::invalid_argument("face_node_velocity must have shape (nface,9,3)");
    int ka=m_kind[row],ha=m_host[row],kb=m_kind[col],hb=m_host[col];
    int la=m_hexLocalOf[row],lb=m_hexLocalOf[col];
    const int oa=ka==0?ha:(int)m_cellCharges.size()+ha;
    const int ob=kb==0?hb:(int)m_cellCharges.size()+hb;
    if(oa>ob){std::swap(ka,kb);std::swap(ha,hb);std::swap(la,lb);}
    const int nb=kb==0?(int)m_cellCharges[hb].size():(int)m_faceCharges[hb].size();
    const auto&block=GetWedgeDirectionalSymBlock(ka,ha,kb,hb,cellVelocity,faceVelocity,cacheToken);
    return block[(size_t)la*nb+lb];
}

double RadHACApKChargeGram::TetChargeGramDirectionalDerivativeElement(
    int row,int col,const std::vector<double>&cellVelocity,const std::vector<double>&faceVelocity,
    unsigned long long cacheToken) const
{
    if(!m_highorder||m_curved||m_hexmode||m_wedgemode)
        throw std::logic_error("TET derivative entry requires a flat polynomial TET Gram");
    if(row<0||row>=m_n||col<0||col>=m_n)throw std::out_of_range("TET derivative entry index out of range");
    if(cellVelocity.size()!=m_hoCellCharges.size()*12)throw std::invalid_argument("cell_vertex_velocity must have shape (ncell,4,3)");
    if(faceVelocity.size()!=m_hoFaceCharges.size()*9)throw std::invalid_argument("face_vertex_velocity must have shape (nface,3,3)");
    const int nc=(int)m_hoCellCharges.size();
    int oa=m_kind[row]==0?m_host[row]:nc+m_host[row],ob=m_kind[col]==0?m_host[col]:nc+m_host[col];
    int la=m_hoLocalOf[row],lb=m_hoLocalOf[col];
    if(oa>ob){std::swap(oa,ob);std::swap(la,lb);}
    const WedgeDerivativeCacheOwner owner{m_build_id,cacheToken,cellVelocity.data(),faceVelocity.data()};
    if(!(s_tet_derivative_owner==owner)){s_tet_derivative_blocks.clear();s_tet_derivative_owner=owner;}
    const HexBlockKey key{0,oa,0,ob,0};auto it=s_tet_derivative_blocks.find(key);
    if(it==s_tet_derivative_blocks.end()){
        if(s_tet_derivative_blocks.size()>HexBlockCacheLimit())s_tet_derivative_blocks.clear();
        it=s_tet_derivative_blocks.emplace(key,TetChargeGramDirectionalDerivativeImpl(cellVelocity,faceVelocity,oa,ob)).first;
    }
    const int kb=ob<nc?0:1,hb=ob<nc?ob:ob-nc;
    const int nb=kb==0?(int)m_hoCellCharges[hb].size():(int)m_hoFaceCharges[hb].size();
    return it->second[(size_t)la*nb+lb];
}

double RadHACApKChargeGram::GetInteractionMatrixElement(int a, int b) const
{
    // Physical entry, EXCEPT while the H-matrix fill is running: the fill's
    // oracle serves the normalized Ghat = sigma_a^-1 sigma_b^-1 * raw so the
    // stored leaves/ACA carry O(1) dynamic range (see the MatVecSym header
    // doc for the roundoff-amplification incident this prevents).  Every
    // caller outside the fill -- Jacobi-diagonal builders, diagnostics,
    // reciprocity gates -- keeps the raw physical Gram.
    const double raw = GetInteractionMatrixElementRaw(a, b);
    if (m_fillNormalized)
        return raw * (m_chargeSigmaInv[(size_t)a] * m_chargeSigmaInv[(size_t)b]);
    return raw;
}

double RadHACApKChargeGram::GetInteractionMatrixElementRaw(int a, int b) const
{
    // Fail-loud bounds guard (2026-07-03 flake hunt): a HACApK-side index bug (1-based lod handling /
    // buffer overrun) would otherwise read garbage hosts and produce plausible-but-wrong blocks.
    if (a < 0 || a >= m_n || b < 0 || b >= m_n)
        throw std::out_of_range("ChargeGram entry index out of range: a=" + std::to_string(a)
                                + " b=" + std::to_string(b) + " n=" + std::to_string(m_n));
    if (m_sampledLaplace) {
        const double dx = m_cent[3 * a] - m_cent[3 * b];
        const double dy = m_cent[3 * a + 1] - m_cent[3 * b + 1];
        const double dz = m_cent[3 * a + 2] - m_cent[3 * b + 2];
        const double eps2 = m_sampledKernelEpsilon * m_sampledKernelEpsilon;
        return m_meas[a] * m_meas[b] * RAD_INV_FOUR_PI /
               std::sqrt(dx * dx + dy * dy + dz * dz + eps2);
    }
    if (m_sampledPlanarLog) {
        const double dx = m_cent[3 * a] - m_cent[3 * b];
        const double dy = m_cent[3 * a + 1] - m_cent[3 * b + 1];
        const double eps2 = m_sampledKernelEpsilon * m_sampledKernelEpsilon;
        const double distance = std::sqrt(dx * dx + dy * dy + eps2);
        return -2.0 * RAD_INV_FOUR_PI * m_meas[a] * m_meas[b]
             * std::log(distance / m_sampledReferenceLength);
    }
    if (m_d2) {
        // 2D planar mode: served block-wise like the hex mode, symmetrized 0.5*(AB + BA).  Each scalar
        // is read BEFORE the next GetHexBlock fetch -- the memo's capacity clear would otherwise leave a
        // dangling reference (the same use-after-free family as the cloud-cache n=10 crash).
        const int kA = m_kind[a], hA = m_host[a], kB = m_kind[b], hB = m_host[b];
        const int la = m_hexLocalOf[a], lb = m_hexLocalOf[b];
        const int nB = (kB == 0) ? (int)m_cellCharges[hB].size() : (int)m_faceCharges[hB].size();
        double base = GetHexSymBlock(kA, hA, kB, hB)[(size_t)la*nB + lb];
        for (size_t i = 0; i < m_image_masks.size(); ++i)
            base += m_image_signs[i]
                  * GetHexSymBlock(kA, hA, kB, hB, (int)i + 1)[(size_t)la*nB + lb];
        return base;
    }
    if (m_hexmode || m_wedgemode) {
        // HEX / WEDGE BDM1: the pair-graded scheme (near subs -> both-domains-graded Duffy outer; far -> the
        // regular symmetric outer; inner always graded/far-dispatched), symmetrized like the other modes.
        // The wedge mode shares this block-serving path verbatim (GetHexBlock -> QuadBlockWedge dispatch).
        // Served from the whole-host-pair block memo (the 64x co-location win) -- bit-identical to
        // the symmetrized 0.5*(block_AB + block_BA) per-entry value, kernel work shared per block.
        // Each scalar is read BEFORE the next GetHexBlock fetch: the memo's capacity clear (fires on
        // ~20k-charge meshes) would otherwise leave a dangling reference.
        const int kA = m_kind[a], hA = m_host[a], kB = m_kind[b], hB = m_host[b];
        const int la = m_hexLocalOf[a], lb = m_hexLocalOf[b];
        const int nB = (kB == 0) ? (int)m_cellCharges[hB].size() : (int)m_faceCharges[hB].size();
        double base = GetHexSymBlock(kA, hA, kB, hB)[(size_t)la*nB + lb];
        // IMA: fold in the mirror-image blocks (the source host REFLECTED on each image mask) so a reduced
        // (1/2,1/4,1/8) symmetry model reproduces the full model -- G_IMA = G + sum_i sign_i*0.5*(refl(a,b)+refl(b,a)).
        // Read each block scalar into a double BEFORE the next GetHexBlock fetch (the memo's capacity clear would
        // otherwise dangle the returned reference -- the same use-after-free family as the direct vAB/vBA above).
        for (size_t i = 0; i < m_image_masks.size(); ++i)
            base += m_image_signs[i] * GetHexSymBlock(kA, hA, kB, hB, (int)i + 1)[(size_t)la*nB + lb];
        return base;
    }
    if (m_highorder) {
        // polynomial charges, symmetrized; the HACApK ACA compresses the well-separated low-rank blocks.
        // NEAR/FAR adaptive quadrature: a well-separated pair uses the cheap LOW-quad plain double-Gauss
        // (QuadDotFar) -- the kernel is smooth there so the expensive HIGH-quad singularity-subtraction is
        // unnecessary; NEAR/self pairs keep the full QuadDot.  This is NOT a monopole far (zero-mean modes
        // have zero monopole) -- it is just a lower quadrature order where the integrand is smooth.
        // m_ho_far_factor = 1e30 (no LOW rule supplied) => every pair NEAR => original all-high-quad path.
        const bool far_pair = a != b && m_ho_far_factor < 1e29 &&
            [&]{ const double dx = m_cent[3*a]-m_cent[3*b], dy = m_cent[3*a+1]-m_cent[3*b+1],
                              dz = m_cent[3*a+2]-m_cent[3*b+2];
                 return std::sqrt(dx*dx + dy*dy + dz*dz) >
                        m_ho_far_factor * (m_size[a] + m_size[b]); }();
        double base;
        const bool use_host_block = m_polyCombo || m_curved ||
            (m_hoAnalyticBlock && HOAnalyticBlockEnabled());
        if (use_host_block && !far_pair) {
            // Flat order<=2 and curved P2 NEAR/self: all monomial charges on a host share the geometry and
            // kernel work at every outer point.  Build and memoize the complete symmetrized host-pair block;
            // the curved path evaluates one Duffy rule for all local modes instead of one rule per entry.
            const int kindA = m_kind[a], hostA = m_host[a], kindB = m_kind[b], hostB = m_host[b];
            const int localA = m_hoLocalOf[a], localB = m_hoLocalOf[b];
            const int nB = (kindB == 0) ? (int)m_hoCellCharges[hostB].size()
                                       : (int)m_hoFaceCharges[hostB].size();
            if (!m_curved ||
                !CurvedTouchBlockValue(kindA, hostA, localA, kindB, hostB, localB, base))
                base = GetHOTetSymBlock(kindA, hostA, kindB, hostB)[(size_t)localA*nB + localB];
        } else if (a == b) {
            base = QuadDot(a, a);                                // curved self fallback
        } else if (far_pair) {
            // FAR: cheap low-quad plain double-Gauss.  Production evaluates both directed rules because a
            // one-sided finite rule is not invariant under explicit mesh reflection.  Keep the environment
            // switch only as a diagnostic timing probe.
            base = HOFarOneSidedEnabled()
                ? QuadDotFar(a, b)
                : 0.5 * (QuadDotFar(a, b) + QuadDotFar(b, a));
        } else {
            base = 0.5 * (QuadDot(a, b) + QuadDot(b, a));         // NEAR: full high-quad subtraction
        }
        // IMA: fold in the mirror-image charge interactions (QuadDotRefl uses PhiInner in this mode) so a
        // reduced (1/2, 1/4, 1/8) symmetry model reproduces the full model -- G_IMA = G + sum_i sign_i *
        // 0.5*(refl(a,b)+refl(b,a)).  Empty image => plain highorder.
        for (size_t i = 0; i < m_image_masks.size(); ++i)
            base += m_image_signs[i] * 0.5 *
                    (QuadDotRefl(a, b, (int)i + 1) + QuadDotRefl(b, a, (int)i + 1));
        return base;
    }
    if (m_analytic) {
        // Diagonal = the analytic self (the Wilton/phi_tet potential is exact through the 1/r singularity).
        double base;
        if (a == b) {
            base = QuadDot(a, a);
        } else {
            // NEAR/FAR split (build speedup): the analytic entry 0.5*(outer-quad_a . Phi_b + outer-quad_b .
            // Phi_a) is EXPENSIVE (PhiTet/TriPotential per outer point) and only matters for NEAR pairs
            // (the non-uniform-M / div M != 0 interaction).  FAR pairs use either a cheap centroid-MONOPOLE
            // (far_quad=0, O((size/r)^2) -- breaks symmetry slightly) or a low-order DOUBLE-QUADRATURE of 1/r
            // (far_quad>0, O((size/r)^4) -- reproduces the all-analytic Gram, the precision-preserving speedup).
            // near_factor = 1e30 (default) => all pairs NEAR => all-analytic (matches the dense
            // analytic reference golden); near_factor ~ 2 gives the fast split.
            const double dx = m_cent[3*a]     - m_cent[3*b];
            const double dy = m_cent[3*a + 1] - m_cent[3*b + 1];
            const double dz = m_cent[3*a + 2] - m_cent[3*b + 2];
            const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (r <= m_near_factor * (m_size[a] + m_size[b]))
                base = 0.5 * (QuadDot(a, b) + QuadDot(b, a));    // NEAR: exact analytic
            else if (m_far_quad > 0)
                base = QuadDotFarLow(a, b);                      // FAR: precision-preserving low-order double-quad
            else
                base = m_meas[a] * m_meas[b] * RAD_INV_FOUR_PI / r;  // FAR: cheap centroid-monopole
        }
        // IMA: fold in the mirror-image charge interactions (always full analytic) so the reduced
        // (1/2,1/4,1/8) model reproduces the full model: G_IMA = G + sum_i sign_i*0.5*(refl(a,b)+refl(b,a)).
        for (size_t i = 0; i < m_image_masks.size(); ++i)
            base += m_image_signs[i] * 0.5 *
                    (QuadDotRefl(a, b, (int)i + 1) + QuadDotRefl(b, a, (int)i + 1));
        return base;
    }
    if (a == b) return m_self[a];
    double dx = m_cent[3 * a + 0] - m_cent[3 * b + 0];
    double dy = m_cent[3 * a + 1] - m_cent[3 * b + 1];
    double dz = m_cent[3 * a + 2] - m_cent[3 * b + 2];
    return m_meas[a] * m_meas[b] * RAD_INV_FOUR_PI / std::sqrt(dx * dx + dy * dy + dz * dz);
}

std::vector<double> RadHACApKChargeGram::SolveLinearMaterial(
    const std::vector<int>& B_indptr, const std::vector<int>& B_indices,
    const std::vector<double>& B_data, int n_face,
    const std::vector<int>& mI, const std::vector<int>& mJ, const std::vector<double>& mV,
    double inv_chi, const std::vector<double>& prec, const std::vector<double>& rhs,
    double tol, int maxit, int& iters_out, bool mass_riesz, bool symmetric,
    const std::vector<double>* x0,
    const std::vector<double>* coarse_basis,
    const std::vector<double>* coarse_applied,
    const std::vector<double>* coarse_factor,
    int coarse_dim)
{
    using Clock = std::chrono::steady_clock;
    auto elapsed = [](Clock::time_point a, Clock::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    };
    m_lastSolveTiming = SolveTiming();
    HACApK_matvec_stats_reset();
    const auto t_total0 = Clock::now();
    const int n_charge = (int)B_indptr.size() - 1;     // B is n_charge x n_face (CSR over charges)
    if ((int)rhs.size() != n_face)
        throw std::runtime_error("SolveLinearMaterial: rhs size mismatch");
    const bool balanced_cluster_prec = coarse_dim > 0;
    if (coarse_dim < 0 ||
        (balanced_cluster_prec &&
         (!coarse_basis || !coarse_applied || !coarse_factor ||
          coarse_basis->size() != static_cast<size_t>(coarse_dim) * n_face ||
          coarse_applied->size() != static_cast<size_t>(coarse_dim) * n_face ||
          coarse_factor->size() != static_cast<size_t>(coarse_dim) * coarse_dim)))
        throw std::runtime_error(
            "SolveLinearMaterial: invalid cluster coarse preconditioner");
    const bool constrained = m_operatorConstrained.size() == (size_t)n_face &&
        std::any_of(m_operatorConstrained.begin(), m_operatorConstrained.end(),
                    [](unsigned char value) { return value != 0; });
    const bool configured_charge_map = m_operatorChargeConfigured &&
        &B_indptr == &m_operatorBIndptr && &B_indices == &m_operatorBIndices &&
        &B_data == &m_operatorBData && m_operatorBTIndptr.size() == (size_t)n_face + 1;
    auto project = [&](std::vector<double>& value) {
        if (!constrained) return;
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t i) {
            if (m_operatorConstrained[i]) value[i] = 0.0;
        });
    };
    // TaskManager self-wrap (AGENTS.md "Parallelization: NGSolve TaskManager"): keep the pool up across
    // the whole CG loop so the Gram H-matvec is parallel without a caller `with TaskManager()`.
    ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
#ifdef HAVE_LAPACK
    // The solve is entered from an active NGSolve TaskManager.  Keep MKL's
    // process setting serial for the non-PARDISO part of this solve; each
    // PARDISO phase temporarily installs its configured local thread count
    // while SuspendTaskManager keeps the NGSolve workers asleep.  HACApK's
    // outer matvec therefore remains TaskManager-parallel without nested
    // worker pools, while the sparse factor/solve can still scale.
    radia::MKLThreadGuard solve_mkl_guard(1);
#endif
    // MASS RIESZ preconditioner (the default 'auto' path): z = M_mass^{-1} r via a single PARDISO SPD
    // factor of the HDiv mass (built once, applied per iteration).  ~3-5x fewer iters than the diagonal
    // Jacobi (the diagonal under-resolves the HDiv mass off-diagonal coupling) and nearly mu_r-flat.  When
    // mass_riesz is false the legacy diagonal Jacobi z = r/prec is used (linear_solver="cpp-cg").
#ifdef HAVE_LAPACK
    // PERSISTENT factor (2026-07-10): get-or-build via EnsureMassRieszFactor (exact-COO key; hit =
    // reuse, factor_s += 0; miss = release-then-refactor).  mrKeep pins the entry for the CG loop.
    std::shared_ptr<RadMassRieszCache> mrKeep;
    MassRieszPardiso* mr = nullptr;
    if (mass_riesz) {
        const std::vector<int>* factorI = &mI;
        const std::vector<int>* factorJ = &mJ;
        const std::vector<double>* factorV = &mV;
        std::vector<int> projectedI, projectedJ;
        std::vector<double> projectedV;
        if (constrained) {
            projectedI.reserve(mI.size() + (size_t)n_face);
            projectedJ.reserve(mJ.size() + (size_t)n_face);
            projectedV.reserve(mV.size() + (size_t)n_face);
            for (size_t k = 0; k < mV.size(); ++k) {
                const int i = mI[k], j = mJ[k];
                if (m_operatorConstrained[(size_t)i] || m_operatorConstrained[(size_t)j]) continue;
                projectedI.push_back(i); projectedJ.push_back(j); projectedV.push_back(mV[k]);
            }
            for (int i = 0; i < n_face; ++i) if (m_operatorConstrained[(size_t)i]) {
                projectedI.push_back(i); projectedJ.push_back(i); projectedV.push_back(1.0);
            }
            factorI = &projectedI; factorJ = &projectedJ; factorV = &projectedV;
        }
        mrKeep = EnsureMassRieszFactor(*factorI, *factorJ, *factorV, n_face, "SolveLinearMaterial",
                                       &m_lastSolveTiming.factor_s);
        mr = &mrKeep->factor;
    }
#else
    if (mass_riesz)
        throw std::runtime_error("SolveLinearMaterial: mass Riesz preconditioner requires MKL PARDISO "
                                 "(HAVE_LAPACK)");
#endif
    const bool configured_mass = m_operatorMassConfigured &&
        &mI == &m_operatorMassI && &mJ == &m_operatorMassJ && &mV == &m_operatorMassV &&
        m_operatorMassIndptr.size() == (size_t)n_face + 1;
    std::vector<int> local_mass_indptr, local_mass_col;
    std::vector<double> local_mass_val;
    if (!configured_mass) {
        local_mass_indptr.assign((size_t)n_face + 1, 0);
        for (size_t k = 0; k < mV.size(); ++k) {
            const int i = mI[k];
            if (i >= 0 && i < n_face && mJ[k] >= 0 && mJ[k] < n_face)
                ++local_mass_indptr[(size_t)i + 1];
        }
        for (int i = 0; i < n_face; ++i)
            local_mass_indptr[(size_t)i + 1] += local_mass_indptr[(size_t)i];
        local_mass_col.resize((size_t)local_mass_indptr[(size_t)n_face]);
        local_mass_val.resize((size_t)local_mass_indptr[(size_t)n_face]);
        std::vector<int> mass_cur = local_mass_indptr;
        for (size_t k = 0; k < mV.size(); ++k) {
            const int i = mI[k], j = mJ[k];
            if (i < 0 || i >= n_face || j < 0 || j >= n_face) continue;
            const int p = mass_cur[(size_t)i]++;
            local_mass_col[(size_t)p] = j;
            local_mass_val[(size_t)p] = mV[k];
        }
    }
    const std::vector<int>& mass_indptr = configured_mass
        ? m_operatorMassIndptr : local_mass_indptr;
    const std::vector<int>& mass_col = configured_mass
        ? m_operatorMassIndices : local_mass_col;
    const std::vector<double>& mass_val = configured_mass
        ? m_operatorMassData : local_mass_val;
    auto coarse_solve = [&](std::vector<double>& value) {
        for (int i = 0; i < coarse_dim; ++i) {
            for (int j = 0; j < i; ++j)
                value[static_cast<size_t>(i)] -=
                    (*coarse_factor)[static_cast<size_t>(i)*coarse_dim+j] *
                    value[static_cast<size_t>(j)];
            value[static_cast<size_t>(i)] /=
                (*coarse_factor)[static_cast<size_t>(i)*coarse_dim+i];
        }
        for (int ii = coarse_dim; ii-- > 0;) {
            for (int j = ii + 1; j < coarse_dim; ++j)
                value[static_cast<size_t>(ii)] -=
                    (*coarse_factor)[static_cast<size_t>(j)*coarse_dim+ii] *
                    value[static_cast<size_t>(j)];
            value[static_cast<size_t>(ii)] /=
                (*coarse_factor)[static_cast<size_t>(ii)*coarse_dim+ii];
        }
    };
    std::vector<double> coarse_r_work(static_cast<size_t>(coarse_dim), 0.0);
    std::vector<double> coarse_d_work(static_cast<size_t>(coarse_dim), 0.0);
    auto applyPrec = [&](const std::vector<double>& rr, std::vector<double>& zz) {
        const auto t0 = Clock::now();
#ifdef HAVE_LAPACK
        if (mass_riesz) {
            mr->Solve(rr.data(), zz.data());
            project(zz);
            m_lastSolveTiming.prec_s += elapsed(t0, Clock::now());
            ++m_lastSolveTiming.prec_count;
            return;
        }
#endif
        if (!balanced_cluster_prec) {
            ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) {
                zz[f] = rr[f] / prec[f];
            });
        }
        else {
            // Balanced two-level preconditioner
            //   P^T D^-1 P + Z (Z^T A Z)^-1 Z^T,
            // P = I - A Z (Z^T A Z)^-1 Z^T.
            // Z and A*Z come directly from the H-matrix cluster tree and are
            // precomputed once per RHS batch.  This is true deflation during
            // every PCG iteration, not merely a one-off coarse initial guess.
            std::fill(coarse_r_work.begin(), coarse_r_work.end(), 0.0);
#ifdef HAVE_LAPACK
            cblas_dgemv(CblasRowMajor, CblasNoTrans, coarse_dim, n_face,
                        1.0, coarse_basis->data(), n_face, rr.data(), 1,
                        0.0, coarse_r_work.data(), 1);
#else
            for (int c = 0; c < coarse_dim; ++c) {
                const double* zc = coarse_basis->data() + static_cast<size_t>(c)*n_face;
                for (int f = 0; f < n_face; ++f)
                    coarse_r_work[static_cast<size_t>(c)] += zc[f] * rr[static_cast<size_t>(f)];
            }
#endif
            coarse_solve(coarse_r_work);
#ifdef HAVE_LAPACK
            std::copy(rr.begin(), rr.end(), zz.begin());
            cblas_dgemv(CblasRowMajor, CblasTrans, coarse_dim, n_face,
                        -1.0, coarse_applied->data(), n_face,
                        coarse_r_work.data(), 1, 1.0, zz.data(), 1);
            ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) {
                zz[f] /= prec[f];
            });
            std::fill(coarse_d_work.begin(), coarse_d_work.end(), 0.0);
            cblas_dgemv(CblasRowMajor, CblasNoTrans, coarse_dim, n_face,
                        1.0, coarse_applied->data(), n_face, zz.data(), 1,
                        0.0, coarse_d_work.data(), 1);
#else
            ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) {
                double projected = rr[f];
                for (int c = 0; c < coarse_dim; ++c)
                    projected -= (*coarse_applied)[static_cast<size_t>(c)*n_face+f] *
                                 coarse_r_work[static_cast<size_t>(c)];
                zz[f] = projected / prec[f];
            });
            std::fill(coarse_d_work.begin(), coarse_d_work.end(), 0.0);
            for (int c = 0; c < coarse_dim; ++c) {
                const double* azc = coarse_applied->data() + static_cast<size_t>(c)*n_face;
                for (int f = 0; f < n_face; ++f)
                    coarse_d_work[static_cast<size_t>(c)] += azc[f] * zz[static_cast<size_t>(f)];
            }
#endif
            coarse_solve(coarse_d_work);
#ifdef HAVE_LAPACK
            for (int c = 0; c < coarse_dim; ++c)
                coarse_d_work[static_cast<size_t>(c)] =
                    coarse_r_work[static_cast<size_t>(c)] -
                    coarse_d_work[static_cast<size_t>(c)];
            cblas_dgemv(CblasRowMajor, CblasTrans, coarse_dim, n_face,
                        1.0, coarse_basis->data(), n_face,
                        coarse_d_work.data(), 1, 1.0, zz.data(), 1);
#else
            ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) {
                for (int c = 0; c < coarse_dim; ++c)
                    zz[f] += (*coarse_basis)[static_cast<size_t>(c)*n_face+f] *
                             (coarse_r_work[static_cast<size_t>(c)] -
                              coarse_d_work[static_cast<size_t>(c)]);
            });
#endif
        }
        m_lastSolveTiming.prec_s += elapsed(t0, Clock::now());
        ++m_lastSolveTiming.prec_count;
    };
    // A x = inv_chi*(M_mass x) + B^T (G (B x)), with G applied as the charge-Gram H-matvec.
    std::vector<double> q((size_t)n_charge), Gq((size_t)n_charge);
    auto applyA = [&](const std::vector<double>& x, std::vector<double>& y) {
        const auto ta0 = Clock::now();
        std::fill(q.begin(), q.end(), 0.0);
        const auto tb0 = Clock::now();
        ngcore::ParallelFor(ngcore::IntRange(n_charge), [&](size_t a) {
            double s = 0.0;
            for (int k = B_indptr[a]; k < B_indptr[a + 1]; ++k) s += B_data[k] * x[B_indices[k]];
            q[a] = s;
        });
        const auto tb1 = Clock::now();
        std::fill(Gq.begin(), Gq.end(), 0.0);
        const auto tg0 = Clock::now();
        if (symmetric) MatVecSym(q, Gq);               // EXACTLY symmetric -> CG-valid Gram apply
        else           MatVec(q, Gq);                  // shadowed: also MatVecSym (sym-fill leaves lower empty)
        const auto tg1 = Clock::now();
        y.assign((size_t)n_face, 0.0);
        const auto tt0 = Clock::now();
        if (configured_charge_map) {
            ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) {
                double sum = 0.0;
                for (int k = m_operatorBTIndptr[f]; k < m_operatorBTIndptr[f + 1]; ++k)
                    sum += m_operatorBTData[(size_t)k] * Gq[(size_t)m_operatorBTIndices[(size_t)k]];
                y[f] = sum;
            });
        }
        else {
            ngcore::ParallelFor(ngcore::IntRange(n_charge), [&](size_t a) {
                const double ga = Gq[a];
                for (int k = B_indptr[a]; k < B_indptr[a + 1]; ++k)
                    ngcore::AtomicAdd(y[(size_t)B_indices[(size_t)k]], B_data[(size_t)k] * ga);
            });
        }
        const auto tt1 = Clock::now();
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t i) {
            double s = 0.0;
            for (int k = mass_indptr[i]; k < mass_indptr[i + 1]; ++k)
                s += mass_val[(size_t)k] * x[mass_col[(size_t)k]];
            y[i] += inv_chi * s;
        });
        project(y);
        const auto ta1 = Clock::now();
        const double bx = elapsed(tb0, tb1);
        const double gm = elapsed(tg0, tg1);
        const double bt = elapsed(tt0, tt1);
        const double ma = elapsed(tt1, ta1);
        const double total = elapsed(ta0, ta1);
        m_lastSolveTiming.bx_s += bx;
        m_lastSolveTiming.gmatvec_s += gm;
        m_lastSolveTiming.btx_s += bt;
        m_lastSolveTiming.mass_s += ma;
        m_lastSolveTiming.ax_total_s += total;
        m_lastSolveTiming.ax_other_s += total - bx - gm - bt - ma;
        ++m_lastSolveTiming.apply_count;
    };
    constexpr int dot_block_size = 4096;
    const int dot_blocks = (n_face + dot_block_size - 1) / dot_block_size;
    std::vector<double> dot_partial((size_t)dot_blocks, 0.0);
    auto dot = [&](const std::vector<double>& a, const std::vector<double>& b) {
        const auto t0 = Clock::now();
        auto sum_block = [&](size_t block) {
            const int begin = static_cast<int>(block) * dot_block_size;
            const int end = std::min(n_face, begin + dot_block_size);
            double sum = 0.0, correction = 0.0;
            for (int f = begin; f < end; ++f) {
                const double value = a[(size_t)f] * b[(size_t)f];
                const double next = sum + value;
                correction += std::fabs(sum) >= std::fabs(value)
                    ? (sum - next) + value : (value - next) + sum;
                sum = next;
            }
            dot_partial[block] = sum + correction;
        };
        if (dot_blocks == 1) sum_block(0);
        else ngcore::ParallelFor(ngcore::IntRange(dot_blocks), sum_block);
        double s = 0.0, correction = 0.0;
        for (double value : dot_partial) {
            const double next = s + value;
            correction += std::fabs(s) >= std::fabs(value)
                ? (s - next) + value : (value - next) + s;
            s = next;
        }
        m_lastSolveTiming.dot_s += elapsed(t0, Clock::now());
        ++m_lastSolveTiming.dot_count;
        return s + correction;
    };
    // Preconditioned conjugate gradients (SPD system; M^{-1} = mass Riesz or 1/prec diagonal Jacobi).
    std::vector<double> rhs_projected = rhs;
    project(rhs_projected);
    std::vector<double> x((size_t)n_face, 0.0), r = rhs_projected, z((size_t)n_face), p((size_t)n_face), Ap;
    if (x0) {
        if ((int)x0->size() != n_face)
            throw std::runtime_error("SolveLinearMaterial: x0 size mismatch");
        x = *x0;
        project(x);
        applyA(x, Ap);
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { r[f] = rhs_projected[f] - Ap[f]; });
    }
    applyPrec(r, z);
    p = z;
    double rz = dot(r, z);
    double bnorm = dot(rhs_projected, rhs_projected);
    bnorm = std::sqrt(bnorm); if (bnorm == 0.0) bnorm = 1.0;
    constexpr int residual_refresh_period = 1000;
    auto recomputeResidual = [&]() {
        applyA(x, Ap);
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) {
            r[f] = rhs_projected[f] - Ap[f];
        });
        project(r);
    };
    int it = 0;
    for (; it < maxit; ++it) {
        double rnorm = dot(r, r);
        if (std::sqrt(rnorm) <= tol * bnorm) {
            // Never accept convergence from the recursive residual alone.
            // Long Jacobi-PCG runs on compressed charge-Gram operators have
            // exhibited a three-decade drift after ~2000 iterations.  Check
            // the true residual and restart from it if the recurrence was
            // optimistic.  This keeps finite differences out of production
            // while making the linear tolerance an actual solver contract.
            recomputeResidual();
            rnorm = dot(r, r);
            if (std::sqrt(rnorm) <= tol * bnorm) break;
            applyPrec(r, z);
            p = z;
            rz = dot(r, z);
            continue;
        }
        applyA(p, Ap);
        double pAp = dot(p, Ap);
        double alpha = rz / pAp;
        const auto tu0 = Clock::now();
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { x[f] += alpha * p[f]; r[f] -= alpha * Ap[f]; });
        const bool refresh = ((it + 1) % residual_refresh_period) == 0;
        if (refresh) recomputeResidual();
        applyPrec(r, z);
        double rz_new = dot(r, z);
        if (refresh) {
            p = z;
        }
        else {
            double beta = rz_new / rz;
            ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { p[f] = z[f] + beta * p[f]; });
        }
        project(x); project(r); project(p);
        m_lastSolveTiming.pcg_update_s += elapsed(tu0, Clock::now());
        rz = rz_new;
    }
    iters_out = it;
    m_lastSolveTiming.total_s = elapsed(t_total0, Clock::now());
    {
        double mv[8] = {0.0};
        int64_t mc[8] = {0};
        HACApK_matvec_stats_get(mv, 8, mc, 8);
        m_lastSolveTiming.hmatvec_total_s = mv[0];
        m_lastSolveTiming.hmatvec_zero_s = mv[1];
        m_lastSolveTiming.hmatvec_permute_s = mv[2];
        m_lastSolveTiming.hmatvec_leaf_s = mv[3];
        m_lastSolveTiming.hmatvec_reduce_s = mv[4];
        m_lastSolveTiming.hmatvec_meta_s = mv[5];
        m_lastSolveTiming.hmatvec_lowrank_flop_est = mv[6];
        m_lastSolveTiming.hmatvec_dense_flop_est = mv[7];
        m_lastSolveTiming.hmatvec_calls = (double)mc[0];
        m_lastSolveTiming.hmatvec_lowrank_leaves = (double)mc[1];
        m_lastSolveTiming.hmatvec_dense_leaves = (double)mc[2];
        m_lastSolveTiming.hmatvec_mirrored_upper_leaves = (double)mc[3];
        m_lastSolveTiming.hmatvec_diagonal_leaves = (double)mc[4];
        m_lastSolveTiming.hmatvec_skipped_lower_leaves = (double)mc[5];
        m_lastSolveTiming.hmatvec_last_nd = (double)mc[6];
        m_lastSolveTiming.hmatvec_last_nthr = (double)mc[7];
    }
    return x;
}

std::vector<double> RadHACApKChargeGram::ApplyDemagOperator(
    const std::vector<int>& B_indptr, const std::vector<int>& B_indices,
    const std::vector<double>& B_data, int n_face,
    const std::vector<double>& x, bool symmetric)
{
    const int n_charge = static_cast<int>(B_indptr.size()) - 1;
    if (n_charge != m_ndof)
        throw std::runtime_error("ApplyDemagOperator: B row count must equal charge-Gram ndof");
    if (n_face < 0 || static_cast<int>(x.size()) != n_face)
        throw std::runtime_error("ApplyDemagOperator: x size mismatch");
    if (B_indices.size() != B_data.size() || B_indptr.empty() ||
        B_indptr.front() != 0 || B_indptr.back() != static_cast<int>(B_data.size()))
        throw std::runtime_error("ApplyDemagOperator: invalid B CSR arrays");
    for (int col : B_indices)
        if (col < 0 || col >= n_face)
            throw std::runtime_error("ApplyDemagOperator: B column index out of range");

    ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
    std::vector<double> q(static_cast<size_t>(n_charge), 0.0);
    std::vector<double> Gq(static_cast<size_t>(n_charge), 0.0);
    std::vector<double> y(static_cast<size_t>(n_face), 0.0);
    ngcore::ParallelFor(ngcore::IntRange(n_charge), [&](size_t a) {
        double sum = 0.0;
        for (int k = B_indptr[a]; k < B_indptr[a + 1]; ++k)
            sum += B_data[static_cast<size_t>(k)] * x[static_cast<size_t>(B_indices[static_cast<size_t>(k)])];
        q[a] = sum;
    });
    if (symmetric) MatVecSym(q, Gq);
    else MatVec(q, Gq);
    ngcore::ParallelFor(ngcore::IntRange(n_charge), [&](size_t a) {
        const double ga = Gq[a];
        for (int k = B_indptr[a]; k < B_indptr[a + 1]; ++k)
            ngcore::AtomicAdd(y[static_cast<size_t>(B_indices[static_cast<size_t>(k)])],
                              B_data[static_cast<size_t>(k)] * ga);
    });
    return y;
}

std::vector<double> RadHACApKChargeGram::ApplyMassRiesz(
    const std::vector<int>& mI, const std::vector<int>& mJ,
    const std::vector<double>& mV, int n_face,
    const std::vector<double>& rhs)
{
    if (n_face < 0 || static_cast<int>(rhs.size()) != n_face)
        throw std::runtime_error("ApplyMassRiesz: rhs size mismatch");
    if (mI.size() != mJ.size() || mI.size() != mV.size())
        throw std::runtime_error("ApplyMassRiesz: mass COO array size mismatch");
#ifdef HAVE_LAPACK
    ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
    radia::MKLThreadGuard solve_mkl_guard(1);
    const bool constrained = m_operatorConstrained.size() == (size_t)n_face &&
        std::any_of(m_operatorConstrained.begin(), m_operatorConstrained.end(),
                    [](unsigned char value) { return value != 0; });
    const std::vector<int>* factorI = &mI;
    const std::vector<int>* factorJ = &mJ;
    const std::vector<double>* factorV = &mV;
    std::vector<int> projectedI, projectedJ;
    std::vector<double> projectedV;
    std::vector<double> projectedRhs = rhs;
    if (constrained) {
        projectedI.reserve(mI.size() + (size_t)n_face);
        projectedJ.reserve(mJ.size() + (size_t)n_face);
        projectedV.reserve(mV.size() + (size_t)n_face);
        for (size_t k = 0; k < mV.size(); ++k) {
            const int i = mI[k], j = mJ[k];
            if (m_operatorConstrained[(size_t)i] || m_operatorConstrained[(size_t)j]) continue;
            projectedI.push_back(i); projectedJ.push_back(j); projectedV.push_back(mV[k]);
        }
        for (int i = 0; i < n_face; ++i) if (m_operatorConstrained[(size_t)i]) {
            projectedI.push_back(i); projectedJ.push_back(i); projectedV.push_back(1.0);
            projectedRhs[(size_t)i] = 0.0;
        }
        factorI = &projectedI; factorJ = &projectedJ; factorV = &projectedV;
    }
    auto keep = EnsureMassRieszFactor(*factorI, *factorJ, *factorV, n_face,
                                      "ApplyMassRiesz", nullptr,
                                      /*geometry_cache=*/true);
    std::vector<double> x(static_cast<size_t>(n_face), 0.0);
    keep->factor.Solve(projectedRhs.data(), x.data());
    if (constrained)
        for (int i = 0; i < n_face; ++i) if (m_operatorConstrained[(size_t)i]) x[(size_t)i] = 0.0;
    return x;
#else
    throw std::runtime_error("ApplyMassRiesz requires MKL PARDISO (HAVE_LAPACK)");
#endif
}

void RadHACApKChargeGram::ConfigureChargeMap(
    std::vector<int> B_indptr, std::vector<int> B_indices,
    std::vector<double> B_data, int n_face)
{
    ConfigureVectorChargeMap(
        std::move(B_indptr), std::move(B_indices), std::move(B_data), n_face, 1);
    if (!m_image_masks.empty() && m_kind.size() == (size_t)m_ndof) {
        const int dimension = m_d2 ? 2 : 3;
        bool positive_axis[3] = {false, false, false};
        for (size_t image = 0; image < m_image_masks.size(); ++image)
            for (int axis = 0; axis < dimension; ++axis)
                if (m_image_masks[image] == (1 << axis) && m_image_signs[image] > 0.0)
                    positive_axis[axis] = true;
        const std::vector<double>* face_nodes = nullptr;
        int stride = 0;
        if (m_d2) { face_nodes = &m_d2EdgeMap; stride = m_d2EdgeMapStride; }
        else if (m_hexmode) { face_nodes = &m_quadNodes; stride = 27; }
        else if (m_wedgemode) { face_nodes = &m_wFaceNodes; stride = 27; }
        else if (m_curved) { face_nodes = &m_faceNodes; stride = 18; }
        else if (m_highorder) { face_nodes = &m_faceV; stride = 9; }
        double scale = 1.0;
        if (face_nodes)
            for (double value : *face_nodes) scale = std::max(scale, std::fabs(value));
        const double plane_tol = 128.0 * std::numeric_limits<double>::epsilon() * scale;
        for (int a = 0; face_nodes && a < m_ndof; ++a) {
            if (m_kind[(size_t)a] != 1) continue;
            const int host = m_host[(size_t)a];
            const double* nodes = &(*face_nodes)[(size_t)host*(size_t)stride];
            int node_count = 0;
            if (m_d2) node_count = m_d2GeometryOrder + 1;
            else if (m_hexmode) node_count = 9;
            else if (m_wedgemode) node_count = m_wFaceType[(size_t)host] == 0 ? 6 : 9;
            else if (m_curved) node_count = 6;
            else node_count = 3;
            bool on_positive_plane = false;
            for (int axis = 0; axis < dimension; ++axis) if (positive_axis[axis]) {
                bool on_plane = true;
                for (int node = 0; node < node_count; ++node)
                    on_plane = on_plane && std::fabs(nodes[dimension*node + axis]) <= plane_tol;
                on_positive_plane = on_positive_plane || on_plane;
            }
            if (!on_positive_plane) continue;
            for (int k = m_operatorBIndptr[(size_t)a];
                 k < m_operatorBIndptr[(size_t)a + 1]; ++k)
                if (std::fabs(m_operatorBData[(size_t)k]) > 0.0)
                    m_operatorConstrained[(size_t)m_operatorBIndices[(size_t)k]] = 1;
        }
    }
}

void RadHACApKChargeGram::ConfigureVectorChargeMap(
    std::vector<int> B_indptr, std::vector<int> B_indices,
    std::vector<double> B_data, int n_face, int n_components)
{
    if (n_face < 0 || n_components < 1 ||
        static_cast<int>(B_indptr.size()) != n_components * m_ndof + 1)
        throw std::runtime_error(
            "ConfigureVectorChargeMap: B row count must equal components times charge-Gram ndof");
    if (B_indices.size() != B_data.size() || B_indptr.empty() || B_indptr.front() != 0 ||
        B_indptr.back() != static_cast<int>(B_data.size()))
        throw std::runtime_error("ConfigureVectorChargeMap: invalid B CSR arrays");
    for (size_t row = 0; row + 1 < B_indptr.size(); ++row) {
        if (B_indptr[row] > B_indptr[row + 1])
            throw std::runtime_error("ConfigureVectorChargeMap: B_indptr must be nondecreasing");
    }
    for (int col : B_indices) {
        if (col < 0 || col >= n_face)
            throw std::runtime_error("ConfigureVectorChargeMap: B column index out of range");
    }
    if (m_operatorMassConfigured && m_operatorNFace != n_face)
        throw std::runtime_error(
            "ConfigureVectorChargeMap: n_face differs from configured mass matrix");
    m_operatorBTIndptr.assign((size_t)n_components * (size_t)n_face + 1, 0);
    for (int component = 0; component < n_components; ++component)
        for (int row = 0; row < m_ndof; ++row)
            for (int k = B_indptr[(size_t)component * (size_t)m_ndof + (size_t)row];
                 k < B_indptr[(size_t)component * (size_t)m_ndof + (size_t)row + 1]; ++k)
                ++m_operatorBTIndptr[(size_t)component * (size_t)n_face +
                                     (size_t)B_indices[(size_t)k] + 1];
    for (size_t row = 0; row + 1 < m_operatorBTIndptr.size(); ++row)
        m_operatorBTIndptr[row + 1] += m_operatorBTIndptr[row];
    m_operatorBTIndices.resize(B_indices.size());
    m_operatorBTData.resize(B_data.size());
    std::vector<int> bt_cursor = m_operatorBTIndptr;
    for (int component = 0; component < n_components; ++component)
        for (int row = 0; row < m_ndof; ++row) {
            const size_t flat_row = (size_t)component * (size_t)m_ndof + (size_t)row;
            for (int k = B_indptr[flat_row]; k < B_indptr[flat_row + 1]; ++k) {
                const int col = B_indices[(size_t)k];
                const int dst = bt_cursor[(size_t)component * (size_t)n_face +
                                          (size_t)col]++;
                m_operatorBTIndices[(size_t)dst] = row;
                m_operatorBTData[(size_t)dst] = B_data[(size_t)k];
            }
        }
    m_operatorConstrained.assign((size_t)n_face, 0);
    m_operatorBIndptr = std::move(B_indptr);
    m_operatorBIndices = std::move(B_indices);
    m_operatorBData = std::move(B_data);
    m_operatorChargeComponents = n_components;
    m_operatorNFace = n_face;
    m_operatorChargeConfigured = true;
}

int RadHACApKChargeGram::ConfiguredConstraintCount() const
{
    return (int)std::count_if(m_operatorConstrained.begin(), m_operatorConstrained.end(),
                              [](unsigned char value) { return value != 0; });
}

void RadHACApKChargeGram::SetConfiguredConstraints(
    const std::vector<int>& dofs, bool preserve_existing)
{
    if (!m_operatorChargeConfigured || m_operatorNFace < 0)
        throw std::runtime_error(
            "SetConfiguredConstraints: charge map must be configured first");
    if (m_operatorConstrained.size() != static_cast<size_t>(m_operatorNFace))
        m_operatorConstrained.assign(static_cast<size_t>(m_operatorNFace), 0);
    else if (!preserve_existing)
        std::fill(m_operatorConstrained.begin(), m_operatorConstrained.end(), 0);
    for (int dof : dofs) {
        if (dof < 0 || dof >= m_operatorNFace)
            throw std::out_of_range(
                "SetConfiguredConstraints: DOF index out of range");
        m_operatorConstrained[static_cast<size_t>(dof)] = 1;
    }
}

void RadHACApKChargeGram::ConfigureMassMatrix(
    std::vector<int> mI, std::vector<int> mJ,
    std::vector<double> mV, int n_face)
{
    if (n_face < 0 || mI.size() != mJ.size() || mI.size() != mV.size())
        throw std::runtime_error("ConfigureMassMatrix: invalid COO arrays");
    if (m_operatorChargeConfigured && m_operatorNFace != n_face)
        throw std::runtime_error("ConfigureMassMatrix: n_face differs from configured charge map");
    for (size_t k = 0; k < mV.size(); ++k) {
        if (mI[k] < 0 || mI[k] >= n_face || mJ[k] < 0 || mJ[k] >= n_face)
            throw std::runtime_error("ConfigureMassMatrix: COO index out of range");
        if (!std::isfinite(mV[k]))
            throw std::runtime_error("ConfigureMassMatrix: non-finite COO value");
    }
    m_operatorMassI = std::move(mI);
    m_operatorMassJ = std::move(mJ);
    m_operatorMassV = std::move(mV);
    m_operatorMassIndptr.assign((size_t)n_face + 1, 0);
    for (int row : m_operatorMassI)
        ++m_operatorMassIndptr[(size_t)row + 1];
    for (int row = 0; row < n_face; ++row)
        m_operatorMassIndptr[(size_t)row + 1] += m_operatorMassIndptr[(size_t)row];
    m_operatorMassIndices.resize(m_operatorMassJ.size());
    m_operatorMassData.resize(m_operatorMassV.size());
    std::vector<int> mass_cursor = m_operatorMassIndptr;
    for (size_t k = 0; k < m_operatorMassV.size(); ++k) {
        const int dst = mass_cursor[(size_t)m_operatorMassI[k]]++;
        m_operatorMassIndices[(size_t)dst] = m_operatorMassJ[k];
        m_operatorMassData[(size_t)dst] = m_operatorMassV[k];
    }
    m_operatorNFace = n_face;
    if (m_operatorConstrained.size() != (size_t)n_face)
        m_operatorConstrained.assign((size_t)n_face, 0);
    m_operatorMassConfigured = true;
    m_operatorMassIsGeometry = false;
}

void RadHACApKChargeGram::ConfigureGeometryMassMatrix(
    std::vector<int> mI, std::vector<int> mJ,
    std::vector<double> mV, int n_face)
{
    if (n_face < 0 || mI.size() != mJ.size() || mI.size() != mV.size())
        throw std::runtime_error("ConfigureGeometryMassMatrix: invalid COO arrays");
    if (m_operatorChargeConfigured && m_operatorNFace != n_face)
        throw std::runtime_error("ConfigureGeometryMassMatrix: n_face differs from configured charge map");
    for (size_t k = 0; k < mV.size(); ++k) {
        if (mI[k] < 0 || mI[k] >= n_face || mJ[k] < 0 || mJ[k] >= n_face)
            throw std::runtime_error("ConfigureGeometryMassMatrix: COO index out of range");
        if (!std::isfinite(mV[k]))
            throw std::runtime_error("ConfigureGeometryMassMatrix: non-finite COO value");
    }
    m_operatorGeometryMassI = std::move(mI);
    m_operatorGeometryMassJ = std::move(mJ);
    m_operatorGeometryMassV = std::move(mV);
    m_operatorGeometryMassIndptr.assign((size_t)n_face + 1, 0);
    for (int row : m_operatorGeometryMassI)
        ++m_operatorGeometryMassIndptr[(size_t)row + 1];
    for (int row = 0; row < n_face; ++row)
        m_operatorGeometryMassIndptr[(size_t)row + 1] += m_operatorGeometryMassIndptr[(size_t)row];
    m_operatorGeometryMassIndices.resize(m_operatorGeometryMassJ.size());
    m_operatorGeometryMassData.resize(m_operatorGeometryMassV.size());
    std::vector<int> geometry_mass_cursor = m_operatorGeometryMassIndptr;
    for (size_t k = 0; k < m_operatorGeometryMassV.size(); ++k) {
        const int dst = geometry_mass_cursor[(size_t)m_operatorGeometryMassI[k]]++;
        m_operatorGeometryMassIndices[(size_t)dst] = m_operatorGeometryMassJ[k];
        m_operatorGeometryMassData[(size_t)dst] = m_operatorGeometryMassV[k];
    }
    m_operatorNFace = n_face;
    if (m_operatorConstrained.size() != (size_t)n_face)
        m_operatorConstrained.assign((size_t)n_face, 0);
    m_operatorGeometryMassConfigured = true;
    m_operatorMassIsGeometry = (
        m_operatorMassConfigured &&
        m_operatorMassI == m_operatorGeometryMassI &&
        m_operatorMassJ == m_operatorGeometryMassJ &&
        m_operatorMassV == m_operatorGeometryMassV);
}

bool RadHACApKChargeGram::RestoreGeometryMassMatrix()
{
    if (!m_operatorGeometryMassConfigured)
        throw std::runtime_error(
            "RestoreGeometryMassMatrix: geometry mass matrix is not configured");
    if (m_operatorMassIsGeometry)
        return false;
    ConfigureMassMatrix(
        m_operatorGeometryMassI, m_operatorGeometryMassJ,
        m_operatorGeometryMassV, m_operatorNFace);
    m_operatorMassIsGeometry = true;
    return true;
}

std::vector<double> RadHACApKChargeGram::ApplyConfiguredDemag(
    const std::vector<double>& x, bool symmetric)
{
    if (!m_operatorChargeConfigured)
        throw std::runtime_error("ApplyConfiguredDemag: charge map is not configured");
    if ((int)x.size() != m_operatorNFace)
        throw std::runtime_error("ApplyConfiguredDemag: x size mismatch");
    std::vector<double> y((size_t)m_operatorNFace, 0.0);
    ApplyConfiguredDemag(x.data(), y.data(), symmetric);
    return y;
}

void RadHACApKChargeGram::ApplyConfiguredDemag(
    const double* x, double* y, bool symmetric)
{
    ApplyConfiguredDemagImpl(x, y, 1.0, false, symmetric);
}

void RadHACApKChargeGram::ApplyConfiguredDemagAdd(
    double scale, const double* x, double* y, bool symmetric)
{
    ApplyConfiguredDemagImpl(x, y, scale, true, symmetric);
}

void RadHACApKChargeGram::ApplyConfiguredDemagImpl(
    const double* x, double* y, double scale, bool add, bool symmetric,
    bool respect_constraints)
{
    if (!m_operatorChargeConfigured)
        throw std::runtime_error("ApplyConfiguredDemag: charge map is not configured");
    if (!x || !y)
        throw std::runtime_error("ApplyConfiguredDemag: null vector data");

    ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
    // BaseMatrix matvecs are called repeatedly from NGSolve Krylov solvers.
    // Keep one workspace per calling thread so repeated applications allocate
    // neither charge vectors nor output temporaries and remain re-entrant.
    static thread_local std::vector<double> q;
    static thread_local std::vector<double> Gq;
    q.resize((size_t)m_ndof);
    Gq.resize((size_t)m_ndof);
    double* const q_data = q.data();
    if (!add) std::fill(y, y + m_operatorNFace, 0.0);
    for (int component = 0; component < m_operatorChargeComponents; ++component) {
        const size_t row_offset = (size_t)component * (size_t)m_ndof;
        ngcore::ParallelFor(ngcore::IntRange(m_ndof), [&](size_t a) {
            double sum = 0.0;
            const size_t row = row_offset + a;
            for (int k = m_operatorBIndptr[row]; k < m_operatorBIndptr[row + 1]; ++k)
                sum += m_operatorBData[(size_t)k] *
                       x[(size_t)m_operatorBIndices[(size_t)k]];
            q_data[a] = sum;
        });
        std::fill(Gq.begin(), Gq.end(), 0.0);
        if (symmetric) MatVecSym(q, Gq);
        else MatVec(q, Gq);
        const double* const Gq_data = Gq.data();
        const size_t bt_offset = (size_t)component * (size_t)m_operatorNFace;
        ngcore::ParallelFor(ngcore::IntRange(m_operatorNFace), [&](size_t f) {
            double sum = 0.0;
            const size_t row = bt_offset + f;
            for (int k = m_operatorBTIndptr[row]; k < m_operatorBTIndptr[row + 1]; ++k)
                sum += m_operatorBTData[(size_t)k] *
                       Gq_data[(size_t)m_operatorBTIndices[(size_t)k]];
            if (!respect_constraints || !m_operatorConstrained[f])
                y[f] += scale * sum;
        });
    }
}

std::vector<double> RadHACApKChargeGram::ApplyConfiguredGeometryMass(
    const std::vector<double>& x)
{
    if ((int)x.size() != m_operatorNFace)
        throw std::runtime_error("ApplyConfiguredGeometryMass: x size mismatch");
    std::vector<double> y((size_t)m_operatorNFace, 0.0);
    ApplyConfiguredGeometryMass(x.data(), y.data());
    return y;
}

void RadHACApKChargeGram::ApplyConfiguredGeometryMass(const double* x, double* y)
{
    if (!m_operatorGeometryMassConfigured)
        throw std::runtime_error("ApplyConfiguredGeometryMass: geometry mass matrix is not configured");
    if (!x || !y)
        throw std::runtime_error("ApplyConfiguredGeometryMass: null vector data");
    ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
    ngcore::ParallelFor(ngcore::IntRange(m_operatorNFace), [&](size_t row) {
        double value = 0.0;
        for (int k = m_operatorGeometryMassIndptr[row];
             k < m_operatorGeometryMassIndptr[row + 1]; ++k)
            value += m_operatorGeometryMassData[(size_t)k] *
                     x[(size_t)m_operatorGeometryMassIndices[(size_t)k]];
        y[row] = m_operatorConstrained[row] ? 0.0 : value;
    });
}

std::vector<double> RadHACApKChargeGram::ApplyConfiguredMassRiesz(
    const std::vector<double>& rhs)
{
    if (!m_operatorGeometryMassConfigured)
        throw std::runtime_error("ApplyConfiguredMassRiesz: geometry mass matrix is not configured");
    return ApplyMassRiesz(m_operatorGeometryMassI, m_operatorGeometryMassJ,
                          m_operatorGeometryMassV,
                          m_operatorNFace, rhs);
}

std::vector<double> RadHACApKChargeGram::SolveConfiguredLinearMaterial(
    double inv_chi, const std::vector<double>& rhs, double tol, int maxit,
    int& iters_out, bool mass_riesz, bool symmetric, const std::vector<double>* x0)
{
    if (!m_operatorChargeConfigured || !m_operatorMassConfigured)
        throw std::runtime_error("SolveConfiguredLinearMaterial: charge map and mass matrix must be configured");
    std::vector<double> no_prec;
    return SolveLinearMaterial(m_operatorBIndptr, m_operatorBIndices, m_operatorBData,
                               m_operatorNFace, m_operatorMassI, m_operatorMassJ,
                               m_operatorMassV, inv_chi, no_prec, rhs, tol, maxit,
                               iters_out, mass_riesz, symmetric, x0);
}

std::vector<double> RadHACApKChargeGram::SolveConfiguredLinearMaterialAutoPrec(
    double inv_chi, const std::vector<double>& rhs, double tol, int maxit,
    int& iters_out, double& prec_min, double& prec_max,
    const std::vector<double>* x0)
{
    if (!m_operatorChargeConfigured || !m_operatorMassConfigured)
        throw std::runtime_error(
            "SolveConfiguredLinearMaterialAutoPrec: charge map and mass matrix must be configured");
    const int n_face = m_operatorNFace;
    if (static_cast<int>(rhs.size()) != n_face)
        throw std::runtime_error("SolveConfiguredLinearMaterialAutoPrec: rhs size mismatch");

    std::vector<double> mass_diag(static_cast<size_t>(n_face), 0.0);
    for (size_t k = 0; k < m_operatorMassV.size(); ++k) {
        if (m_operatorMassI[k] == m_operatorMassJ[k])
            mass_diag[static_cast<size_t>(m_operatorMassI[k])] += m_operatorMassV[k];
    }
    std::vector<std::vector<int>> support_id(static_cast<size_t>(n_face));
    std::vector<std::vector<double>> support_value(static_cast<size_t>(n_face));
    for (int a = 0; a < m_ndof; ++a) {
        for (int k = m_operatorBIndptr[static_cast<size_t>(a)];
             k < m_operatorBIndptr[static_cast<size_t>(a) + 1]; ++k) {
            const int f = m_operatorBIndices[static_cast<size_t>(k)];
            support_id[static_cast<size_t>(f)].push_back(a);
            support_value[static_cast<size_t>(f)].push_back(m_operatorBData[static_cast<size_t>(k)]);
        }
    }
    std::vector<double> prec(static_cast<size_t>(n_face), 0.0);
    {
        ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) {
            double ndiag = 0.0;
            const auto& ids = support_id[f];
            const auto& vals = support_value[f];
            for (size_t p = 0; p < ids.size(); ++p)
                for (size_t q = 0; q < ids.size(); ++q)
                    ndiag += vals[p] * vals[q] * GetInteractionMatrixElement(ids[p], ids[q]);
            double value = inv_chi * mass_diag[f] + ndiag;
            if (!(value > 0.0) || !std::isfinite(value)) value = 1.0;
            prec[f] = value;
        });
    }
    prec_min = n_face ? prec[0] : 0.0;
    prec_max = prec_min;
    for (double value : prec) {
        prec_min = std::min(prec_min, value);
        prec_max = std::max(prec_max, value);
    }
    return SolveLinearMaterial(m_operatorBIndptr, m_operatorBIndices, m_operatorBData,
                               n_face, m_operatorMassI, m_operatorMassJ, m_operatorMassV,
                               inv_chi, prec, rhs, tol, maxit, iters_out,
                                /*mass_riesz=*/false, /*symmetric=*/true, x0);
}

std::vector<double> RadHACApKChargeGram::ApplyConfiguredLinearMaterialOperator(
    double inv_chi, const std::vector<double>& x, bool respect_constraints)
{
    if (!m_operatorChargeConfigured || !m_operatorMassConfigured)
        throw std::runtime_error(
            "ApplyConfiguredLinearMaterialOperator: charge map and mass matrix must be configured");
    if (static_cast<int>(x.size()) != m_operatorNFace)
        throw std::runtime_error(
            "ApplyConfiguredLinearMaterialOperator: x size mismatch");
    if (!std::isfinite(inv_chi) || inv_chi < 0.0)
        throw std::runtime_error(
            "ApplyConfiguredLinearMaterialOperator: inv_chi must be finite and nonnegative");
    std::vector<double> y(static_cast<size_t>(m_operatorNFace), 0.0);
    ApplyConfiguredDemagImpl(x.data(), y.data(), 1.0, false, true,
                             respect_constraints);
    ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
    ngcore::ParallelFor(ngcore::IntRange(m_operatorNFace), [&](size_t row) {
        if (respect_constraints && m_operatorConstrained[row]) {
            y[row] = 0.0;
            return;
        }
        double value = 0.0;
        for (int k = m_operatorMassIndptr[row];
             k < m_operatorMassIndptr[row + 1]; ++k)
            value += m_operatorMassData[static_cast<size_t>(k)] *
                     x[static_cast<size_t>(m_operatorMassIndices[static_cast<size_t>(k)])];
        y[row] += inv_chi * value;
    });
    return y;
}

std::vector<double> RadHACApKChargeGram::ApplyConfiguredLinearMaterialOperatorMany(
    double inv_chi, const std::vector<double>& x, int nrhs,
    bool respect_constraints)
{
    if (!m_operatorChargeConfigured || !m_operatorMassConfigured)
        throw std::runtime_error(
            "ApplyConfiguredLinearMaterialOperatorMany: charge map and mass matrix must be configured");
    const int n_face = m_operatorNFace;
    if (nrhs < 1 || static_cast<int64_t>(x.size()) !=
            static_cast<int64_t>(nrhs)*n_face)
        throw std::runtime_error(
            "ApplyConfiguredLinearMaterialOperatorMany: x must be row-major [nrhs][n_face]");
    if (!std::isfinite(inv_chi) || inv_chi < 0.0)
        throw std::runtime_error(
            "ApplyConfiguredLinearMaterialOperatorMany: inv_chi must be finite and nonnegative");
    const size_t face_total = static_cast<size_t>(nrhs)*n_face;
    std::vector<double> y(face_total, 0.0);
    for (int component = 0; component < m_operatorChargeComponents; ++component) {
        std::vector<double> charge(static_cast<size_t>(nrhs)*m_ndof, 0.0);
        const size_t row_offset = static_cast<size_t>(component)*m_ndof;
        {
            ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
            ngcore::ParallelFor(
                ngcore::IntRange(static_cast<size_t>(nrhs)*m_ndof),
                [&](size_t index) {
                    const int rhs = static_cast<int>(index/m_ndof);
                    const int row = static_cast<int>(index%m_ndof);
                    double value = 0.0;
                    const size_t mapped = row_offset+static_cast<size_t>(row);
                    for (int k = m_operatorBIndptr[mapped];
                         k < m_operatorBIndptr[mapped+1]; ++k)
                        value += m_operatorBData[static_cast<size_t>(k)] *
                            x[static_cast<size_t>(rhs)*n_face+
                              m_operatorBIndices[static_cast<size_t>(k)]];
                    charge[index] = value;
                });
        }
        std::vector<double> gcharge;
        MatVecSymMany(charge, nrhs, gcharge);
        const size_t bt_offset = static_cast<size_t>(component)*n_face;
        {
            ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
            ngcore::ParallelFor(ngcore::IntRange(face_total), [&](size_t index) {
                const int rhs = static_cast<int>(index/n_face);
                const int face = static_cast<int>(index%n_face);
                if (respect_constraints &&
                    m_operatorConstrained[static_cast<size_t>(face)]) return;
                double value = 0.0;
                const size_t mapped = bt_offset+static_cast<size_t>(face);
                for (int k = m_operatorBTIndptr[mapped];
                     k < m_operatorBTIndptr[mapped+1]; ++k)
                    value += m_operatorBTData[static_cast<size_t>(k)] *
                        gcharge[static_cast<size_t>(rhs)*m_ndof+
                                m_operatorBTIndices[static_cast<size_t>(k)]];
                y[index] += value;
            });
        }
    }
    {
        ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
        ngcore::ParallelFor(ngcore::IntRange(face_total), [&](size_t index) {
            const int rhs = static_cast<int>(index/n_face);
            const int row = static_cast<int>(index%n_face);
            if (respect_constraints &&
                m_operatorConstrained[static_cast<size_t>(row)]) {
                y[index] = 0.0;
                return;
            }
            double value = 0.0;
            for (int k = m_operatorMassIndptr[static_cast<size_t>(row)];
                 k < m_operatorMassIndptr[static_cast<size_t>(row)+1]; ++k)
                value += m_operatorMassData[static_cast<size_t>(k)] *
                    x[static_cast<size_t>(rhs)*n_face+
                      m_operatorMassIndices[static_cast<size_t>(k)]];
            y[index] += inv_chi*value;
        });
    }
    return y;
}

std::vector<double>
RadHACApKChargeGram::ConfiguredLinearMaterialElementBlocks(
    double inv_chi, const std::vector<int>& candidate_dofs,
    const std::vector<int>& block_offsets)
{
    if (!m_operatorChargeConfigured || !m_operatorMassConfigured)
        throw std::runtime_error(
            "ConfiguredLinearMaterialElementBlocks: charge map and mass matrix must be configured");
    if (!std::isfinite(inv_chi) || inv_chi < 0.0)
        throw std::runtime_error(
            "ConfiguredLinearMaterialElementBlocks: inv_chi must be finite and nonnegative");
    if (block_offsets.size() < 2 || block_offsets.front() != 0 ||
        block_offsets.back() != static_cast<int>(candidate_dofs.size()))
        throw std::runtime_error(
            "ConfiguredLinearMaterialElementBlocks: offsets must start at zero and end at candidate size");
    const int n_face = m_operatorNFace;
    std::vector<unsigned char> seen(static_cast<size_t>(n_face), 0);
    for (int dof : candidate_dofs) {
        if (dof < 0 || dof >= n_face || seen[static_cast<size_t>(dof)])
            throw std::runtime_error(
                "ConfiguredLinearMaterialElementBlocks: candidate DOFs must be unique and in range");
        seen[static_cast<size_t>(dof)] = 1;
    }
    const int n_block = static_cast<int>(block_offsets.size()) - 1;
    std::vector<size_t> value_offsets(static_cast<size_t>(n_block) + 1, 0);
    for (int block = 0; block < n_block; ++block) {
        const int begin = block_offsets[static_cast<size_t>(block)];
        const int end = block_offsets[static_cast<size_t>(block) + 1];
        if (begin < 0 || end <= begin ||
            end > static_cast<int>(candidate_dofs.size()))
            throw std::runtime_error(
                "ConfiguredLinearMaterialElementBlocks: offsets must define nonempty increasing blocks");
        const size_t width = static_cast<size_t>(end - begin);
        value_offsets[static_cast<size_t>(block) + 1] =
            value_offsets[static_cast<size_t>(block)] + width * width;
    }
    std::vector<double> values(value_offsets.back(), 0.0);
    {
        ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
        ngcore::ParallelFor(ngcore::IntRange(n_block), [&](size_t block_index) {
            const int block = static_cast<int>(block_index);
            const int begin = block_offsets[block_index];
            const int end = block_offsets[block_index + 1];
            const int width = end - begin;
            const size_t output_begin = value_offsets[block_index];
            double* output = values.data() + output_begin;

            std::unordered_map<int, int> local_face;
            local_face.reserve(static_cast<size_t>(2 * width));
            for (int local = 0; local < width; ++local)
                local_face.emplace(
                    candidate_dofs[static_cast<size_t>(begin + local)], local);

            // Exact element-local material mass from the configured CSR.
            for (int local_row = 0; local_row < width; ++local_row) {
                const int row = candidate_dofs[static_cast<size_t>(
                    begin + local_row)];
                for (int k = m_operatorMassIndptr[static_cast<size_t>(row)];
                     k < m_operatorMassIndptr[static_cast<size_t>(row) + 1]; ++k) {
                    const auto found = local_face.find(
                        m_operatorMassIndices[static_cast<size_t>(k)]);
                    if (found != local_face.end())
                        output[static_cast<size_t>(local_row) * width +
                               found->second] +=
                            inv_chi * m_operatorMassData[static_cast<size_t>(k)];
                }
            }

            // For each vector-charge component, gather only the scalar-charge
            // modes touched by this element, form its exact small G block, and
            // contract B_e^T G_e B_e.  No global H-matrix matvec is performed.
            for (int component = 0; component < m_operatorChargeComponents;
                 ++component) {
                std::vector<int> charges;
                std::unordered_map<int, int> local_charge;
                for (int local = 0; local < width; ++local) {
                    const int face = candidate_dofs[static_cast<size_t>(
                        begin + local)];
                    const size_t mapped = static_cast<size_t>(component) *
                        n_face + static_cast<size_t>(face);
                    for (int k = m_operatorBTIndptr[mapped];
                         k < m_operatorBTIndptr[mapped + 1]; ++k) {
                        const int charge =
                            m_operatorBTIndices[static_cast<size_t>(k)];
                        if (local_charge.emplace(
                                charge, static_cast<int>(charges.size())).second)
                            charges.push_back(charge);
                    }
                }
                const int n_charge = static_cast<int>(charges.size());
                if (n_charge == 0) continue;
                std::vector<double> local_b(
                    static_cast<size_t>(n_charge) * width, 0.0);
                for (int local = 0; local < width; ++local) {
                    const int face = candidate_dofs[static_cast<size_t>(
                        begin + local)];
                    const size_t mapped = static_cast<size_t>(component) *
                        n_face + static_cast<size_t>(face);
                    for (int k = m_operatorBTIndptr[mapped];
                         k < m_operatorBTIndptr[mapped + 1]; ++k) {
                        const int charge =
                            m_operatorBTIndices[static_cast<size_t>(k)];
                        const int q = local_charge.find(charge)->second;
                        local_b[static_cast<size_t>(q) * width + local] +=
                            m_operatorBTData[static_cast<size_t>(k)];
                    }
                }
                std::vector<double> gb(
                    static_cast<size_t>(n_charge) * width, 0.0);
                for (int q = 0; q < n_charge; ++q)
                    for (int r = 0; r < n_charge; ++r) {
                        const double g = GetInteractionMatrixElement(
                            charges[static_cast<size_t>(q)],
                            charges[static_cast<size_t>(r)]);
                        const double* b_row = local_b.data() +
                            static_cast<size_t>(r) * width;
                        double* gb_row = gb.data() +
                            static_cast<size_t>(q) * width;
                        for (int column = 0; column < width; ++column)
                            gb_row[column] += g * b_row[column];
                    }
                for (int row = 0; row < width; ++row)
                    for (int column = 0; column < width; ++column) {
                        double value = 0.0;
                        for (int q = 0; q < n_charge; ++q)
                            value += local_b[static_cast<size_t>(q) * width + row] *
                                     gb[static_cast<size_t>(q) * width + column];
                        output[static_cast<size_t>(row) * width + column] += value;
                    }
            }
            for (int row = 0; row < width; ++row)
                for (int column = row + 1; column < width; ++column) {
                    const size_t rc = static_cast<size_t>(row) * width + column;
                    const size_t cr = static_cast<size_t>(column) * width + row;
                    const double symmetric = 0.5 * (output[rc] + output[cr]);
                    output[rc] = symmetric;
                    output[cr] = symmetric;
                }
        });
    }
    if (!std::all_of(values.begin(), values.end(),
                     [](double value) { return std::isfinite(value); }))
        throw std::runtime_error(
            "ConfiguredLinearMaterialElementBlocks: non-finite local block");
    return values;
}

std::vector<int>
RadHACApKChargeGram::ConfiguredLinearMaterialCandidateClusters(
    const std::vector<int>& candidate_dofs,
    const std::vector<int>& block_offsets, int requested_clusters,
    int& n_cluster_out)
{
    if (!m_operatorChargeConfigured)
        throw std::runtime_error(
            "ConfiguredLinearMaterialCandidateClusters: charge map must be configured");
    if (requested_clusters < 1)
        throw std::runtime_error(
            "ConfiguredLinearMaterialCandidateClusters: requested_clusters must be positive");
    if (block_offsets.size() < 2 || block_offsets.front() != 0 ||
        block_offsets.back() != static_cast<int>(candidate_dofs.size()))
        throw std::runtime_error(
            "ConfiguredLinearMaterialCandidateClusters: offsets must start at zero and end at candidate size");
    const int n_face = m_operatorNFace;
    std::vector<unsigned char> seen(static_cast<size_t>(n_face), 0);
    for (int dof : candidate_dofs) {
        if (dof < 0 || dof >= n_face || seen[static_cast<size_t>(dof)])
            throw std::runtime_error(
                "ConfiguredLinearMaterialCandidateClusters: candidate DOFs must be unique and in range");
        seen[static_cast<size_t>(dof)] = 1;
    }
    const int n_block = static_cast<int>(block_offsets.size()) - 1;
    for (int block = 0; block < n_block; ++block)
        if (block_offsets[static_cast<size_t>(block)] < 0 ||
            block_offsets[static_cast<size_t>(block) + 1] <=
                block_offsets[static_cast<size_t>(block)] ||
            block_offsets[static_cast<size_t>(block) + 1] >
                static_cast<int>(candidate_dofs.size()))
            throw std::runtime_error(
                "ConfiguredLinearMaterialCandidateClusters: offsets must define nonempty increasing blocks");

    auto* leaf = static_cast<st_cHACApK_leafmtxp>(m_leafmtxp);
    auto* ctl = static_cast<st_cHACApK_lcontrol>(m_control);
    if (!leaf || !leaf->st_clt_root || !ctl || !ctl->lod)
        throw std::runtime_error(
            "ConfiguredLinearMaterialCandidateClusters: H-matrix cluster tree is unavailable");

    // Use actual tree nodes, not a geometric reconstruction in Python.  The
    // same largest-node splitting rule is used by the cluster Ritz solver, so
    // screening and H-matvec share one spatial hierarchy and permutation.
    std::vector<st_cHACApK_cluster> clusters;
    clusters.push_back(leaf->st_clt_root);
    while (static_cast<int>(clusters.size()) < requested_clusters) {
        int split = -1, largest = -1;
        for (int i = 0; i < static_cast<int>(clusters.size()); ++i) {
            const auto node = clusters[static_cast<size_t>(i)];
            if (node && node->nnson > 0 &&
                static_cast<int>(clusters.size()) + node->nnson - 1 <=
                    requested_clusters &&
                node->nsize > largest) {
                split = i;
                largest = node->nsize;
            }
        }
        if (split < 0) break;
        const auto node = clusters[static_cast<size_t>(split)];
        clusters.erase(clusters.begin() + split);
        for (int child = 1; child <= node->nnson; ++child)
            clusters.push_back(node->pc_sons[child]);
    }
    n_cluster_out = static_cast<int>(clusters.size());
    std::vector<int> charge_cluster(static_cast<size_t>(m_ndof), -1);
    for (int cluster = 0; cluster < n_cluster_out; ++cluster) {
        const auto node = clusters[static_cast<size_t>(cluster)];
        for (int pos = node->nstrt; pos < node->nstrt + node->nsize; ++pos) {
            if (pos < 1 || pos > m_ndof)
                throw std::runtime_error(
                    "ConfiguredLinearMaterialCandidateClusters: invalid cluster range");
            const int charge = ctl->lod[pos] - 1;
            if (charge < 0 || charge >= m_ndof)
                throw std::runtime_error(
                    "ConfiguredLinearMaterialCandidateClusters: invalid cluster permutation");
            charge_cluster[static_cast<size_t>(charge)] = cluster;
        }
    }

    std::vector<int> labels(static_cast<size_t>(n_block), 0);
    ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
    ngcore::ParallelFor(ngcore::IntRange(n_block), [&](size_t block_index) {
        std::vector<double> weights(static_cast<size_t>(n_cluster_out), 0.0);
        const int begin = block_offsets[block_index];
        const int end = block_offsets[block_index + 1];
        for (int local = begin; local < end; ++local) {
            const int face = candidate_dofs[static_cast<size_t>(local)];
            for (int component = 0; component < m_operatorChargeComponents;
                 ++component) {
                const size_t mapped = static_cast<size_t>(component) * n_face +
                    static_cast<size_t>(face);
                for (int entry = m_operatorBTIndptr[mapped];
                     entry < m_operatorBTIndptr[mapped + 1]; ++entry) {
                    const int charge =
                        m_operatorBTIndices[static_cast<size_t>(entry)];
                    const int cluster =
                        charge_cluster[static_cast<size_t>(charge)];
                    if (cluster >= 0)
                        weights[static_cast<size_t>(cluster)] += std::fabs(
                            m_operatorBTData[static_cast<size_t>(entry)]);
                }
            }
        }
        const auto best = std::max_element(weights.begin(), weights.end());
        if (best == weights.end() || !(*best > 0.0))
            throw std::runtime_error(
                "ConfiguredLinearMaterialCandidateClusters: candidate block has no charge support");
        labels[block_index] = static_cast<int>(best - weights.begin());
    });
    return labels;
}

RadHACApKChargeGram::CandidateSchurReduction
RadHACApKChargeGram::ReduceConfiguredCandidateSchur(
    double inv_chi, const std::vector<int>& candidate_dofs,
    const std::vector<double>& rhs, const std::vector<double>& state,
    const std::vector<double>& response_matrix,
    const std::vector<double>& adjoints, int n_response,
    double tol, int maxit, int solve_batch_size, bool mass_riesz)
{
    using Clock = std::chrono::steady_clock;
    if (!m_operatorChargeConfigured || !m_operatorMassConfigured)
        throw std::runtime_error(
            "ReduceConfiguredCandidateSchur: charge map and mass matrix must be configured");
    const int n_face = m_operatorNFace;
    const int nc = static_cast<int>(candidate_dofs.size());
    if (!std::isfinite(inv_chi) || inv_chi < 0.0 ||
        !std::isfinite(tol) || tol <= 0.0 || tol >= 1.0 || maxit < 1)
        throw std::runtime_error(
            "ReduceConfiguredCandidateSchur: invalid material or solver parameters");
    if (nc < 1 || n_response < 1 || solve_batch_size < 1)
        throw std::runtime_error(
            "ReduceConfiguredCandidateSchur: candidate, response, and batch dimensions must be positive");
    if (static_cast<int>(rhs.size()) != n_face ||
        static_cast<int>(state.size()) != n_face ||
        static_cast<int64_t>(response_matrix.size()) !=
            static_cast<int64_t>(n_response)*n_face ||
        response_matrix.size() != adjoints.size())
        throw std::runtime_error(
            "ReduceConfiguredCandidateSchur: state/response array shape mismatch");
    std::vector<unsigned char> seen(static_cast<size_t>(n_face), 0);
    for (int dof : candidate_dofs) {
        if (dof < 0 || dof >= n_face || seen[static_cast<size_t>(dof)])
            throw std::runtime_error(
                "ReduceConfiguredCandidateSchur: candidate DOFs must be unique and in range");
        if (!m_operatorConstrained[static_cast<size_t>(dof)])
            throw std::runtime_error(
                "ReduceConfiguredCandidateSchur: candidate DOFs must be outside the active set");
        seen[static_cast<size_t>(dof)] = 1;
    }

    CandidateSchurReduction result;
    result.n_candidate = nc;
    result.n_response = n_response;
    const auto operator_started = Clock::now();
    std::vector<double> basis(static_cast<size_t>(nc)*n_face, 0.0);
    for (int column = 0; column < nc; ++column)
        basis[static_cast<size_t>(column)*n_face+
              candidate_dofs[static_cast<size_t>(column)]] = 1.0;
    std::vector<double> columns = ApplyConfiguredLinearMaterialOperatorMany(
        inv_chi, basis, nc, /*respect_constraints=*/false);
    basis.clear(); basis.shrink_to_fit();
    result.operator_s = std::chrono::duration<double>(
        Clock::now()-operator_started).count();

    // A_aE in row-major [candidate][face], with inactive rows projected out.
    std::vector<double> coupling = columns;
    {
        ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
        ngcore::ParallelFor(
            ngcore::IntRange(static_cast<size_t>(nc)*n_face),
            [&](size_t index) {
                if (m_operatorConstrained[index%static_cast<size_t>(n_face)])
                    coupling[index] = 0.0;
            });
    }
    const auto solve_started = Clock::now();
    // A broken-element mass matrix has no active/candidate cross block, hence
    // A_aE is entirely the charge interaction B_a^T G B_E.  BDM1 contains a
    // substantial divergence-free local subspace, so its 36 coefficient
    // columns are often a much smaller-rank collection of active RHS vectors.
    // Compute a thin native TSVD once, solve only its orthonormal right-singular
    // rows, then reconstruct A_aa^-1 A_aE.  This is algebraically exact for the
    // retained row space and keeps the final whole-element full solve as the
    // physical acceptance gate.
    std::vector<double> solve_rhs = coupling;
    std::vector<double> left_singular;
    std::vector<double> singular_values;
    int coupling_rank = nc;
#ifdef HAVE_LAPACK
    {
        double coupling_square = 0.0;
        for (double value : coupling) coupling_square += value*value;
        if (coupling_square == 0.0) {
            coupling_rank = 0;
            solve_rhs.clear();
        }
        else {
            std::vector<double> factor_input = coupling;
            singular_values.assign(static_cast<size_t>(nc), 0.0);
            left_singular.assign(static_cast<size_t>(nc)*nc, 0.0);
            std::vector<double> right_singular(
                static_cast<size_t>(nc)*n_face, 0.0);
            const int info = LAPACKE_dgesdd(
                LAPACK_ROW_MAJOR, 'S', nc, n_face, factor_input.data(), n_face,
                singular_values.data(), left_singular.data(), nc,
                right_singular.data(), n_face);
            if (info != 0)
                throw std::runtime_error(
                    "ReduceConfiguredCandidateSchur: coupling TSVD failed");
            const double scale = singular_values.front();
            const double relative_cutoff = std::max(1.0e-13, 1.0e-2*tol);
            coupling_rank = 0;
            double total_square = 0.0, discarded_square = 0.0;
            for (double value : singular_values) total_square += value*value;
            for (int mode = 0; mode < nc; ++mode) {
                const double value = singular_values[static_cast<size_t>(mode)];
                if (value > relative_cutoff*scale)
                    ++coupling_rank;
                else
                    discarded_square += value*value;
            }
            solve_rhs.assign(
                right_singular.begin(),right_singular.begin()+
                static_cast<size_t>(coupling_rank)*n_face);
            result.coupling_relative_truncation_error =
                std::sqrt(discarded_square/total_square);
        }
    }
#endif
    result.coupling_rank = coupling_rank;
    std::vector<double> solved_modes(
        static_cast<size_t>(coupling_rank)*n_face, 0.0);
    result.coupling_mode_iterations.reserve(
        static_cast<size_t>(coupling_rank));
    for (int begin = 0; begin < coupling_rank; begin += solve_batch_size) {
        const int count = std::min(solve_batch_size, coupling_rank-begin);
        const auto first = solve_rhs.begin()+static_cast<size_t>(begin)*n_face;
        std::vector<double> chunk(first, first+static_cast<size_t>(count)*n_face);
        std::vector<int> iterations;
        double pmin = 0.0, pmax = 0.0, coarse_setup_s = 0.0, projection_s = 0.0;
        int coarse_dim = 0, recycle_dim = 0;
        std::vector<double> solved = SolveConfiguredLinearMaterialAutoPrecMany(
            inv_chi, chunk, count, tol, maxit,
            /*cluster_coarse_size=*/0, /*cluster_deflation_size=*/0,
            /*recycle_size=*/0, iterations, pmin, pmax,
            coarse_dim, recycle_dim, coarse_setup_s, projection_s,
            mass_riesz, nullptr);
        std::copy(solved.begin(), solved.end(),
                  solved_modes.begin()+static_cast<size_t>(begin)*n_face);
        result.coupling_mode_iterations.insert(
            result.coupling_mode_iterations.end(),
            iterations.begin(), iterations.end());
    }
    // ``iters`` predates coupling compression and is candidate-length public
    // output.  Keep its shape stable as a conservative upper-bound summary;
    // exact retained-mode counts are exposed separately.
    const int iteration_upper_bound = result.coupling_mode_iterations.empty()
        ? 0 : *std::max_element(result.coupling_mode_iterations.begin(),
                                result.coupling_mode_iterations.end());
    result.iterations.assign(static_cast<size_t>(nc), iteration_upper_bound);
    std::vector<double> active_solutions(static_cast<size_t>(nc)*n_face, 0.0);
#ifdef HAVE_LAPACK
    if (!left_singular.empty() && coupling_rank > 0) {
        std::vector<double> scaled_left(
            static_cast<size_t>(nc)*coupling_rank, 0.0);
        for (int row = 0; row < nc; ++row)
            for (int mode = 0; mode < coupling_rank; ++mode)
                scaled_left[static_cast<size_t>(row)*coupling_rank+mode] =
                    left_singular[static_cast<size_t>(row)*nc+mode]*
                    singular_values[static_cast<size_t>(mode)];
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            nc, n_face, coupling_rank, 1.0,
            scaled_left.data(), coupling_rank,
            solved_modes.data(), n_face, 0.0,
            active_solutions.data(), n_face);
    }
    else
#endif
        active_solutions.swap(solved_modes);
    result.solve_s = std::chrono::duration<double>(
        Clock::now()-solve_started).count();

    const auto contraction_started = Clock::now();
    result.schur.assign(static_cast<size_t>(nc)*nc, 0.0);
    for (int p = 0; p < nc; ++p)
        for (int q = 0; q < nc; ++q)
            result.schur[static_cast<size_t>(p)*nc+q] =
                columns[static_cast<size_t>(q)*n_face+
                        candidate_dofs[static_cast<size_t>(p)]];
    result.rhs.resize(static_cast<size_t>(nc));
    result.response.assign(static_cast<size_t>(n_response)*nc, 0.0);
    for (int p = 0; p < nc; ++p) {
        result.rhs[static_cast<size_t>(p)] =
            rhs[static_cast<size_t>(candidate_dofs[static_cast<size_t>(p)])];
        for (int output = 0; output < n_response; ++output)
            result.response[static_cast<size_t>(output)*nc+p] =
                response_matrix[static_cast<size_t>(output)*n_face+
                                candidate_dofs[static_cast<size_t>(p)]];
    }
#ifdef HAVE_LAPACK
    if (coupling_rank > 0) {
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
            nc, nc, n_face, -1.0, coupling.data(), n_face,
            active_solutions.data(), n_face, 1.0, result.schur.data(), nc);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, nc, n_face,
            -1.0, coupling.data(), n_face, state.data(), 1,
            1.0, result.rhs.data(), 1);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
            n_response, nc, n_face, -1.0, adjoints.data(), n_face,
            coupling.data(), n_face, 1.0, result.response.data(), nc);
    }
#else
    for (int p = 0; p < nc; ++p) {
        for (int face = 0; face < n_face; ++face)
            result.rhs[static_cast<size_t>(p)] -=
                coupling[static_cast<size_t>(p)*n_face+face]*
                state[static_cast<size_t>(face)];
        for (int q = 0; q < nc; ++q)
            for (int face = 0; face < n_face; ++face)
                result.schur[static_cast<size_t>(p)*nc+q] -=
                    coupling[static_cast<size_t>(p)*n_face+face]*
                    active_solutions[static_cast<size_t>(q)*n_face+face];
        for (int output = 0; output < n_response; ++output)
            for (int face = 0; face < n_face; ++face)
                result.response[static_cast<size_t>(output)*nc+p] -=
                    adjoints[static_cast<size_t>(output)*n_face+face]*
                    coupling[static_cast<size_t>(p)*n_face+face];
    }
#endif
    for (int p = 0; p < nc; ++p)
        for (int q = p+1; q < nc; ++q) {
            const double value = 0.5*(
                result.schur[static_cast<size_t>(p)*nc+q]+
                result.schur[static_cast<size_t>(q)*nc+p]);
            result.schur[static_cast<size_t>(p)*nc+q] = value;
            result.schur[static_cast<size_t>(q)*nc+p] = value;
        }
    result.contraction_s = std::chrono::duration<double>(
        Clock::now()-contraction_started).count();
    return result;
}

std::vector<double> RadHACApKChargeGram::SolveConfiguredLinearMaterialAutoPrecMany(
    double inv_chi, const std::vector<double>& rhs, int nrhs,
    double tol, int maxit, int cluster_coarse_size,
    int cluster_deflation_size, int recycle_size,
    std::vector<int>& iters_out, double& prec_min, double& prec_max,
    int& coarse_dim_out, int& recycle_dim_out,
    double& coarse_setup_s, double& projection_s,
    bool mass_riesz,
    const std::vector<double>* x0)
{
    using Clock = std::chrono::steady_clock;
    if (!m_operatorChargeConfigured || !m_operatorMassConfigured)
        throw std::runtime_error(
            "SolveConfiguredLinearMaterialAutoPrecMany: charge map and mass matrix must be configured");
    const int n_face = m_operatorNFace;
    if (nrhs < 1 || static_cast<int64_t>(rhs.size()) !=
            static_cast<int64_t>(nrhs) * n_face)
        throw std::runtime_error(
            "SolveConfiguredLinearMaterialAutoPrecMany: rhs must be row-major [nrhs][n_face]");
    if (x0 && x0->size() != rhs.size())
        throw std::runtime_error(
            "SolveConfiguredLinearMaterialAutoPrecMany: x0 shape mismatch");
    if (cluster_coarse_size < 0 || cluster_deflation_size < 0 || recycle_size < 0)
        throw std::runtime_error(
            "SolveConfiguredLinearMaterialAutoPrecMany: coarse/recycle sizes must be non-negative");
    if (cluster_coarse_size > 0 && cluster_deflation_size < 1)
        throw std::runtime_error(
            "SolveConfiguredLinearMaterialAutoPrecMany: positive cluster size requires deflation size");
    if (mass_riesz && cluster_coarse_size > 0)
        throw std::runtime_error(
            "SolveConfiguredLinearMaterialAutoPrecMany: mass-Riesz block PCG "
            "currently requires cluster_coarse_size=0");

    // Build the exact system diagonal once for the complete RHS batch.  The
    // former Python loop rebuilt this support map and diagonal for every
    // state/adjoint even though all columns share A.
    std::vector<double> mass_diag(static_cast<size_t>(n_face), 0.0);
    for (size_t k = 0; k < m_operatorMassV.size(); ++k)
        if (m_operatorMassI[k] == m_operatorMassJ[k])
            mass_diag[static_cast<size_t>(m_operatorMassI[k])] += m_operatorMassV[k];
    std::vector<std::vector<int>> support_id(static_cast<size_t>(n_face));
    std::vector<std::vector<double>> support_value(static_cast<size_t>(n_face));
    for (int a = 0; a < m_ndof; ++a)
        for (int k = m_operatorBIndptr[static_cast<size_t>(a)];
             k < m_operatorBIndptr[static_cast<size_t>(a) + 1]; ++k) {
            const int f = m_operatorBIndices[static_cast<size_t>(k)];
            support_id[static_cast<size_t>(f)].push_back(a);
            support_value[static_cast<size_t>(f)].push_back(
                m_operatorBData[static_cast<size_t>(k)]);
        }
    std::vector<double> prec(static_cast<size_t>(n_face), 0.0);
    {
        ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) {
            double ndiag = 0.0;
            const auto& ids = support_id[f];
            const auto& vals = support_value[f];
            for (size_t p = 0; p < ids.size(); ++p)
                for (size_t q = 0; q < ids.size(); ++q)
                    ndiag += vals[p] * vals[q] *
                        GetInteractionMatrixElement(ids[p], ids[q]);
            double value = inv_chi * mass_diag[f] + ndiag;
            if (!(value > 0.0) || !std::isfinite(value)) value = 1.0;
            prec[f] = value;
        });
    }
    prec_min = n_face ? prec[0] : 0.0;
    prec_max = prec_min;
    for (double value : prec) {
        prec_min = std::min(prec_min, value);
        prec_max = std::max(prec_max, value);
    }

#ifdef HAVE_LAPACK
    std::shared_ptr<RadMassRieszCache> block_mr_keep;
    MassRieszPardiso* block_mr = nullptr;
    if (mass_riesz) {
        block_mr_keep = EnsureMassRieszFactor(
            m_operatorMassI, m_operatorMassJ, m_operatorMassV, n_face,
            "SolveConfiguredLinearMaterialAutoPrecMany", nullptr);
        block_mr = &block_mr_keep->factor;
    }
#else
    if (mass_riesz)
        throw std::runtime_error(
            "SolveConfiguredLinearMaterialAutoPrecMany: mass Riesz requires "
            "MKL PARDISO (HAVE_LAPACK)");
#endif

    auto dot = [&](const double* a, const double* b) {
        double sum = 0.0, correction = 0.0;
        for (int i = 0; i < n_face; ++i) {
            const double value = a[i] * b[i];
            const double next = sum + value;
            correction += std::fabs(sum) >= std::fabs(value)
                ? (sum - next) + value : (value - next) + sum;
            sum = next;
        }
        return sum + correction;
    };
    auto apply_system = [&](const double* input, std::vector<double>& output) {
        output.assign(static_cast<size_t>(n_face), 0.0);
        ApplyConfiguredDemag(input, output.data(), true);
        ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t row) {
            if (m_operatorConstrained[row]) { output[row] = 0.0; return; }
            double value = 0.0;
            for (int k = m_operatorMassIndptr[row];
                 k < m_operatorMassIndptr[row + 1]; ++k)
                value += m_operatorMassData[static_cast<size_t>(k)] *
                         input[static_cast<size_t>(m_operatorMassIndices[static_cast<size_t>(k)])];
            output[row] += inv_chi * value;
        });
    };
    auto cholesky = [](std::vector<double>& matrix, int n, const char* who) {
        double max_diag = 0.0;
        for (int i = 0; i < n; ++i)
            max_diag = std::max(max_diag, std::fabs(matrix[static_cast<size_t>(i)*n+i]));
        const double floor = std::max(1.0e-30, 1.0e-13 * max_diag);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                double value = 0.5 * (matrix[static_cast<size_t>(i)*n+j] +
                                      matrix[static_cast<size_t>(j)*n+i]);
                for (int k = 0; k < j; ++k)
                    value -= matrix[static_cast<size_t>(i)*n+k] *
                             matrix[static_cast<size_t>(j)*n+k];
                if (i == j) {
                    if (!(value > floor) || !std::isfinite(value))
                        throw std::runtime_error(std::string(who) +
                            ": Galerkin matrix is not numerically SPD");
                    matrix[static_cast<size_t>(i)*n+j] = std::sqrt(value);
                }
                else
                    matrix[static_cast<size_t>(i)*n+j] =
                        value / matrix[static_cast<size_t>(j)*n+j];
            }
            for (int j = i + 1; j < n; ++j)
                matrix[static_cast<size_t>(i)*n+j] = 0.0;
        }
    };
    auto chol_solve = [](const std::vector<double>& factor, int n,
                         std::vector<double>& value) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < i; ++j)
                value[static_cast<size_t>(i)] -=
                    factor[static_cast<size_t>(i)*n+j] * value[static_cast<size_t>(j)];
            value[static_cast<size_t>(i)] /= factor[static_cast<size_t>(i)*n+i];
        }
        for (int ii = n; ii-- > 0;) {
            for (int j = ii + 1; j < n; ++j)
                value[static_cast<size_t>(ii)] -=
                    factor[static_cast<size_t>(j)*n+ii] * value[static_cast<size_t>(j)];
            value[static_cast<size_t>(ii)] /= factor[static_cast<size_t>(ii)*n+ii];
        }
    };

    // Select actual nodes from the preserved binary cluster tree.  Repeatedly
    // split the largest current node until the requested coarse dimension is
    // reached; ranges and lod are the same permutation used by H-matvec.
    const auto t_coarse0 = Clock::now();
    std::vector<st_cHACApK_cluster> clusters;
    if (cluster_coarse_size > 0) {
        auto* leaf = static_cast<st_cHACApK_leafmtxp>(m_leafmtxp);
        auto* ctl = static_cast<st_cHACApK_lcontrol>(m_control);
        if (!leaf || !leaf->st_clt_root || !ctl || !ctl->lod)
            throw std::runtime_error(
                "SolveConfiguredLinearMaterialAutoPrecMany: H-matrix cluster tree is unavailable");
        clusters.push_back(leaf->st_clt_root);
        while (static_cast<int>(clusters.size()) < cluster_coarse_size) {
            int split = -1, largest = -1;
            for (int i = 0; i < static_cast<int>(clusters.size()); ++i) {
                const auto node = clusters[static_cast<size_t>(i)];
                if (node && node->nnson > 0 &&
                    static_cast<int>(clusters.size()) + node->nnson - 1 <= cluster_coarse_size &&
                    node->nsize > largest) {
                    split = i; largest = node->nsize;
                }
            }
            if (split < 0) break;
            const auto node = clusters[static_cast<size_t>(split)];
            clusters.erase(clusters.begin() + split);
            for (int child = 1; child <= node->nnson; ++child)
                clusters.push_back(node->pc_sons[child]);
        }
        std::vector<int> charge_cluster(static_cast<size_t>(m_ndof), -1);
        for (int c = 0; c < static_cast<int>(clusters.size()); ++c) {
            const auto node = clusters[static_cast<size_t>(c)];
            for (int pos = node->nstrt; pos < node->nstrt + node->nsize; ++pos) {
                if (pos < 1 || pos > m_ndof)
                    throw std::runtime_error(
                        "SolveConfiguredLinearMaterialAutoPrecMany: invalid cluster range");
                const int charge = ctl->lod[pos] - 1;
                if (charge < 0 || charge >= m_ndof)
                    throw std::runtime_error(
                        "SolveConfiguredLinearMaterialAutoPrecMany: invalid cluster permutation");
                charge_cluster[static_cast<size_t>(charge)] = c;
            }
        }
        // Lift each charge-tree aggregate through D^-1 B^T.  A piecewise
        // constant FACE vector is not an HDiv aggregate and was measured to
        // worsen the 40k-DOF condition number.  B^T 1_cluster is instead the
        // aggregate boundary flux of that charge cluster (interior faces
        // cancel); diagonal scaling makes it the natural Jacobi coarse mode.
        std::vector<double> raw_z(clusters.size() * static_cast<size_t>(n_face), 0.0);
        for (int charge = 0; charge < m_ndof; ++charge) {
            const int c = charge_cluster[static_cast<size_t>(charge)];
            if (c < 0) continue;
            for (int k = m_operatorBIndptr[static_cast<size_t>(charge)];
                 k < m_operatorBIndptr[static_cast<size_t>(charge) + 1]; ++k) {
                const int face = m_operatorBIndices[static_cast<size_t>(k)];
                if (!m_operatorConstrained[static_cast<size_t>(face)])
                    raw_z[static_cast<size_t>(c)*n_face+face] +=
                        m_operatorBData[static_cast<size_t>(k)];
            }
        }
        std::vector<double> z;
        z.reserve(clusters.size() * static_cast<size_t>(n_face));
        int dim = 0;
        std::vector<double> candidate(static_cast<size_t>(n_face), 0.0);
        for (size_t c = 0; c < clusters.size(); ++c) {
            for (int face = 0; face < n_face; ++face)
                candidate[static_cast<size_t>(face)] =
                    raw_z[c*static_cast<size_t>(n_face)+face] / prec[static_cast<size_t>(face)];
            const double raw_norm = std::sqrt(std::max(
                0.0, dot(candidate.data(), candidate.data())));
            for (int pass = 0; pass < 2; ++pass)
                for (int previous = 0; previous < dim; ++previous) {
                    const double* zp = z.data() + static_cast<size_t>(previous)*n_face;
                    const double coefficient = dot(zp, candidate.data());
                    for (int face = 0; face < n_face; ++face)
                        candidate[static_cast<size_t>(face)] -= coefficient * zp[face];
                }
            const double norm = std::sqrt(std::max(
                0.0, dot(candidate.data(), candidate.data())));
            if (raw_norm > 1.0e-30 && norm > 1.0e-10 * raw_norm) {
                for (double& value : candidate) value /= norm;
                z.insert(z.end(), candidate.begin(), candidate.end());
                ++dim;
            }
        }
        clusters.clear();
        clusters.shrink_to_fit();
        raw_z.clear();
        raw_z.shrink_to_fit();
        std::vector<double> az(static_cast<size_t>(dim) * n_face, 0.0);
        for (int c = 0; c < dim; ++c) {
            std::vector<double> applied;
            apply_system(z.data() + static_cast<size_t>(c)*n_face, applied);
            std::copy(applied.begin(), applied.end(),
                      az.begin() + static_cast<size_t>(c)*n_face);
        }
        std::vector<double> coarse_factor(static_cast<size_t>(dim)*dim, 0.0);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                coarse_factor[static_cast<size_t>(i)*dim+j] = dot(
                    z.data() + static_cast<size_t>(i)*n_face,
                    az.data() + static_cast<size_t>(j)*n_face);

        // The aggregate span contains both smooth and high-energy boundary
        // fluxes.  Deflating all of it was measured slower and increased the
        // 40k-DOF iteration count.  Contract the actual cluster-tree Galerkin
        // matrix and retain only its lowest Rayleigh--Ritz modes.
        const int keep_requested = std::min(dim, cluster_deflation_size);
        if (keep_requested < dim) {
#ifdef HAVE_LAPACK
            std::vector<double> eigenvectors = coarse_factor;
            std::vector<double> eigenvalues(static_cast<size_t>(dim), 0.0);
            const int info = LAPACKE_dsyev(
                LAPACK_ROW_MAJOR, 'V', 'U', dim, eigenvectors.data(), dim,
                eigenvalues.data());
            if (info != 0)
                throw std::runtime_error(
                    "SolveConfiguredLinearMaterialAutoPrecMany: cluster Rayleigh-Ritz eigensolve failed");
            // Compressed H-matrices can leave roundoff-scale negative Ritz
            // values even while all solved physical directions are positive.
            // Never inject those into an SPD PCG preconditioner: retain the
            // lowest strictly positive modes above a relative spectral floor.
            double spectral_scale = 0.0;
            for (double value : eigenvalues)
                spectral_scale = std::max(spectral_scale, std::fabs(value));
            const double spectral_floor = std::max(1.0e-30, 1.0e-11*spectral_scale);
            std::vector<int> selected;
            for (int i = 0; i < dim && static_cast<int>(selected.size()) < keep_requested; ++i)
                if (eigenvalues[static_cast<size_t>(i)] > spectral_floor)
                    selected.push_back(i);
            if (selected.empty())
                throw std::runtime_error(
                    "SolveConfiguredLinearMaterialAutoPrecMany: no positive cluster Ritz mode");
            const int keep = static_cast<int>(selected.size());
            std::vector<double> selected_z(static_cast<size_t>(keep)*n_face, 0.0);
            std::vector<double> selected_az(static_cast<size_t>(keep)*n_face, 0.0);
            for (int mode = 0; mode < keep; ++mode)
                for (int aggregate = 0; aggregate < dim; ++aggregate) {
                    // LAPACKE row-major returns the mathematical eigenvector
                    // matrix: eigenvectors are columns.
                    const double coefficient =
                        eigenvectors[static_cast<size_t>(aggregate)*dim+
                                     selected[static_cast<size_t>(mode)]];
                    const double* za = z.data() + static_cast<size_t>(aggregate)*n_face;
                    const double* aza = az.data() + static_cast<size_t>(aggregate)*n_face;
                    double* zs = selected_z.data() + static_cast<size_t>(mode)*n_face;
                    double* azs = selected_az.data() + static_cast<size_t>(mode)*n_face;
                    for (int face = 0; face < n_face; ++face) {
                        zs[face] += coefficient * za[face];
                        azs[face] += coefficient * aza[face];
                    }
                }
            z.swap(selected_z);
            az.swap(selected_az);
            dim = keep;
            coarse_factor.assign(static_cast<size_t>(dim)*dim, 0.0);
            for (int i = 0; i < dim; ++i)
                for (int j = 0; j < dim; ++j)
                    coarse_factor[static_cast<size_t>(i)*dim+j] = dot(
                        z.data() + static_cast<size_t>(i)*n_face,
                        az.data() + static_cast<size_t>(j)*n_face);
#else
            throw std::runtime_error(
                "SolveConfiguredLinearMaterialAutoPrecMany: cluster Rayleigh-Ritz requires LAPACK");
#endif
        }
        coarse_dim_out = dim;
        if (dim > 0) cholesky(coarse_factor, dim, "cluster coarse correction");
        coarse_setup_s = std::chrono::duration<double>(Clock::now() - t_coarse0).count();

        std::vector<double> solutions(static_cast<size_t>(nrhs)*n_face, 0.0);
        std::vector<double> recycle_u, recycle_au;
        int recycle_dim = 0;
        iters_out.assign(static_cast<size_t>(nrhs), 0);
        projection_s = 0.0;
        for (int column = 0; column < nrhs; ++column) {
            const double* b = rhs.data() + static_cast<size_t>(column)*n_face;
            std::vector<double> seed(static_cast<size_t>(n_face), 0.0);
            if (x0)
                std::copy(x0->begin() + static_cast<size_t>(column)*n_face,
                          x0->begin() + static_cast<size_t>(column+1)*n_face,
                          seed.begin());
            const auto tp0 = Clock::now();
            std::vector<double> applied;
            apply_system(seed.data(), applied);
            std::vector<double> residual(static_cast<size_t>(n_face));
            for (int face = 0; face < n_face; ++face)
                residual[static_cast<size_t>(face)] = b[face] - applied[static_cast<size_t>(face)];
            if (dim > 0) {
                std::vector<double> coeff(static_cast<size_t>(dim), 0.0);
                for (int c = 0; c < dim; ++c)
                    coeff[static_cast<size_t>(c)] = dot(
                        z.data() + static_cast<size_t>(c)*n_face, residual.data());
                chol_solve(coarse_factor, dim, coeff);
                for (int c = 0; c < dim; ++c) {
                    const double value = coeff[static_cast<size_t>(c)];
                    const double* zc = z.data() + static_cast<size_t>(c)*n_face;
                    const double* azc = az.data() + static_cast<size_t>(c)*n_face;
                    for (int face = 0; face < n_face; ++face) {
                        seed[static_cast<size_t>(face)] += value * zc[face];
                        residual[static_cast<size_t>(face)] -= value * azc[face];
                    }
                }
            }
            if (recycle_dim > 0) {
                std::vector<double> factor(static_cast<size_t>(recycle_dim)*recycle_dim, 0.0);
                std::vector<double> coeff(static_cast<size_t>(recycle_dim), 0.0);
                for (int i = 0; i < recycle_dim; ++i) {
                    const double* ui = recycle_u.data() + static_cast<size_t>(i)*n_face;
                    coeff[static_cast<size_t>(i)] = dot(ui, residual.data());
                    for (int j = 0; j < recycle_dim; ++j)
                        factor[static_cast<size_t>(i)*recycle_dim+j] = dot(
                            ui, recycle_au.data() + static_cast<size_t>(j)*n_face);
                }
                cholesky(factor, recycle_dim, "Ritz recycle correction");
                chol_solve(factor, recycle_dim, coeff);
                for (int c = 0; c < recycle_dim; ++c) {
                    const double value = coeff[static_cast<size_t>(c)];
                    const double* uc = recycle_u.data() + static_cast<size_t>(c)*n_face;
                    for (int face = 0; face < n_face; ++face)
                        seed[static_cast<size_t>(face)] += value * uc[face];
                }
            }
            projection_s += std::chrono::duration<double>(Clock::now() - tp0).count();
            std::vector<double> one_rhs(b, b + n_face);
            int iterations = 0;
            std::vector<double> solution = SolveLinearMaterial(
                m_operatorBIndptr, m_operatorBIndices, m_operatorBData,
                n_face, m_operatorMassI, m_operatorMassJ, m_operatorMassV,
                inv_chi, prec, one_rhs, tol, maxit, iterations,
                /*mass_riesz=*/false, /*symmetric=*/true, &seed,
                &z, &az, &coarse_factor, dim);
            iters_out[static_cast<size_t>(column)] = iterations;
            std::copy(solution.begin(), solution.end(),
                      solutions.begin() + static_cast<size_t>(column)*n_face);

            if (recycle_dim < recycle_size) {
                std::vector<double> candidate = solution;
                for (int c = 0; c < recycle_dim; ++c) {
                    const double* uc = recycle_u.data() + static_cast<size_t>(c)*n_face;
                    const double coeff = dot(uc, candidate.data());
                    for (int face = 0; face < n_face; ++face)
                        candidate[static_cast<size_t>(face)] -= coeff * uc[face];
                }
                const double norm = std::sqrt(std::max(0.0, dot(candidate.data(), candidate.data())));
                const double solution_norm = std::sqrt(std::max(0.0, dot(solution.data(), solution.data())));
                if (norm > 1.0e-10 * std::max(1.0, solution_norm)) {
                    for (double& value : candidate) value /= norm;
                    std::vector<double> a_candidate;
                    apply_system(candidate.data(), a_candidate);
                    recycle_u.insert(recycle_u.end(), candidate.begin(), candidate.end());
                    recycle_au.insert(recycle_au.end(), a_candidate.begin(), a_candidate.end());
                    ++recycle_dim;
                }
            }
        }
        recycle_dim_out = recycle_dim;
        return solutions;
    }

    // Explicitly disabled coarse space: block PCG.  State and adjoint columns
    // share one Krylov space, so a single BLAS-3 H-matrix traversal advances
    // every right-hand side.  The public storage remains row-major
    // [nrhs][n_face]; the small dense matrices below use the mathematical
    // column convention internally.
    coarse_dim_out = 0;
    coarse_setup_s = std::chrono::duration<double>(Clock::now() - t_coarse0).count();
    projection_s = 0.0;
    recycle_dim_out = 0;
    const size_t total = static_cast<size_t>(nrhs)*n_face;
    auto project_many = [&](std::vector<double>& value) {
        ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
        ngcore::ParallelFor(ngcore::IntRange(total), [&](size_t index) {
            if (m_operatorConstrained[index % static_cast<size_t>(n_face)])
                value[index] = 0.0;
        });
    };
    auto apply_preconditioner_many = [&](const std::vector<double>& input,
                                         std::vector<double>& output,
                                         const std::vector<int>& rows) {
        std::fill(output.begin(), output.end(), 0.0);
#ifdef HAVE_LAPACK
        if (mass_riesz) {
            block_mr->SolveMany(input.data(), output.data(), nrhs);
            project_many(output);
            return;
        }
#endif
        for (int row : rows)
            for (int face = 0; face < n_face; ++face)
                output[static_cast<size_t>(row)*n_face+face] =
                    input[static_cast<size_t>(row)*n_face+face] /
                    prec[static_cast<size_t>(face)];
    };
    auto apply_system_many = [&](const std::vector<double>& input,
                                 std::vector<double>& output) {
        std::vector<double> charge(static_cast<size_t>(nrhs)*m_ndof, 0.0);
        {
            ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
            ngcore::ParallelFor(ngcore::IntRange(static_cast<size_t>(nrhs)*m_ndof),
                [&](size_t index) {
                    const int column = static_cast<int>(index / m_ndof);
                    const int row = static_cast<int>(index % m_ndof);
                    double value = 0.0;
                    for (int k = m_operatorBIndptr[static_cast<size_t>(row)];
                         k < m_operatorBIndptr[static_cast<size_t>(row)+1]; ++k)
                        value += m_operatorBData[static_cast<size_t>(k)] *
                            input[static_cast<size_t>(column)*n_face+
                                  m_operatorBIndices[static_cast<size_t>(k)]];
                    charge[index] = value;
                });
        }
        std::vector<double> gcharge;
        MatVecSymMany(charge, nrhs, gcharge);
        output.assign(total, 0.0);
        {
            ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
            ngcore::ParallelFor(ngcore::IntRange(total), [&](size_t index) {
                const int column = static_cast<int>(index / n_face);
                const int face = static_cast<int>(index % n_face);
                if (m_operatorConstrained[static_cast<size_t>(face)]) return;
                double value = 0.0;
                for (int k = m_operatorBTIndptr[static_cast<size_t>(face)];
                     k < m_operatorBTIndptr[static_cast<size_t>(face)+1]; ++k)
                    value += m_operatorBTData[static_cast<size_t>(k)] *
                        gcharge[static_cast<size_t>(column)*m_ndof+
                                m_operatorBTIndices[static_cast<size_t>(k)]];
                double mass_value = 0.0;
                for (int k = m_operatorMassIndptr[static_cast<size_t>(face)];
                     k < m_operatorMassIndptr[static_cast<size_t>(face)+1]; ++k)
                    mass_value += m_operatorMassData[static_cast<size_t>(k)] *
                        input[static_cast<size_t>(column)*n_face+
                              m_operatorMassIndices[static_cast<size_t>(k)]];
                output[index] = value + inv_chi*mass_value;
            });
        }
    };
    auto row_dot = [&](const std::vector<double>& a,
                       const std::vector<double>& b, int row_a, int row_b) {
        return dot(a.data()+static_cast<size_t>(row_a)*n_face,
                   b.data()+static_cast<size_t>(row_b)*n_face);
    };
    auto block_gram = [&](const std::vector<double>& left,
                          const std::vector<double>& right,
                          const std::vector<int>& rows) {
        const int dim = static_cast<int>(rows.size());
        std::vector<double> result(static_cast<size_t>(dim)*dim, 0.0);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                result[static_cast<size_t>(i)*dim+j] =
                    row_dot(left, right, rows[static_cast<size_t>(i)],
                            rows[static_cast<size_t>(j)]);
        return result;
    };
    auto symmetric_pseudoinverse_apply = [&](std::vector<double> matrix,
                                              const std::vector<double>& right,
                                              int dim, const char* who) {
#ifdef HAVE_LAPACK
        for (int i = 0; i < dim; ++i)
            for (int j = i+1; j < dim; ++j) {
                const double value = 0.5*(
                    matrix[static_cast<size_t>(i)*dim+j] +
                    matrix[static_cast<size_t>(j)*dim+i]);
                matrix[static_cast<size_t>(i)*dim+j] = value;
                matrix[static_cast<size_t>(j)*dim+i] = value;
            }
        std::vector<double> eigenvalues(static_cast<size_t>(dim), 0.0);
        const int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', dim,
                                       matrix.data(), dim, eigenvalues.data());
        if (info != 0)
            throw std::runtime_error(std::string(who)+": eigensolve failed");
        const double scale = std::max(1.0e-300,
            *std::max_element(eigenvalues.begin(), eigenvalues.end()));
        const double floor = 1.0e-13*scale;
        std::vector<double> transformed(static_cast<size_t>(dim)*dim, 0.0);
        std::vector<double> result(static_cast<size_t>(dim)*dim, 0.0);
        int rank = 0;
        for (int mode = 0; mode < dim; ++mode) {
            const double lambda = eigenvalues[static_cast<size_t>(mode)];
            if (!(lambda > floor)) continue;
            ++rank;
            for (int column = 0; column < dim; ++column)
                for (int row = 0; row < dim; ++row)
                    transformed[static_cast<size_t>(mode)*dim+column] +=
                        matrix[static_cast<size_t>(row)*dim+mode] *
                        right[static_cast<size_t>(row)*dim+column] / lambda;
        }
        if (rank == 0)
            throw std::runtime_error(std::string(who)+": numerical rank is zero");
        for (int row = 0; row < dim; ++row)
            for (int column = 0; column < dim; ++column)
                for (int mode = 0; mode < dim; ++mode)
                    result[static_cast<size_t>(row)*dim+column] +=
                        matrix[static_cast<size_t>(row)*dim+mode] *
                        transformed[static_cast<size_t>(mode)*dim+column];
        return result;
#else
        cholesky(matrix, dim, who);
        std::vector<double> result(static_cast<size_t>(dim)*dim, 0.0);
        for (int column = 0; column < dim; ++column) {
            std::vector<double> value(static_cast<size_t>(dim), 0.0);
            for (int row = 0; row < dim; ++row)
                value[static_cast<size_t>(row)] =
                    right[static_cast<size_t>(row)*dim+column];
            chol_solve(matrix, dim, value);
            for (int row = 0; row < dim; ++row)
                result[static_cast<size_t>(row)*dim+column] = value[static_cast<size_t>(row)];
        }
        return result;
#endif
    };
    auto add_block_product = [&](std::vector<double>& target,
                                 const std::vector<double>& source,
                                 const std::vector<double>& coefficients,
                                 const std::vector<int>& rows, double scale) {
        const int dim = static_cast<int>(rows.size());
        ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
        ngcore::ParallelFor(ngcore::IntRange(static_cast<size_t>(dim)*n_face),
            [&](size_t index) {
                const int local_target = static_cast<int>(index/n_face);
                const int face = static_cast<int>(index%n_face);
                double value = 0.0;
                for (int local_source = 0; local_source < dim; ++local_source)
                    value += source[static_cast<size_t>(
                        rows[static_cast<size_t>(local_source)])*n_face+face] *
                        coefficients[static_cast<size_t>(local_source)*dim+
                                     local_target];
                target[static_cast<size_t>(rows[static_cast<size_t>(local_target)])*
                       n_face+face] += scale*value;
            });
    };

    // Normalize columns before the block recurrence.  The physical state and
    // point-field adjoints differ by many orders of magnitude; without this
    // scaling their 7x7 Gram matrices lose the smaller columns numerically.
    std::vector<double> projected_rhs = rhs;
    project_many(projected_rhs);
    std::vector<double> rhs_scale(static_cast<size_t>(nrhs), 1.0);
    std::vector<int> active;
    iters_out.assign(static_cast<size_t>(nrhs), maxit);
    for (int row = 0; row < nrhs; ++row) {
        rhs_scale[static_cast<size_t>(row)] = std::sqrt(std::max(
            0.0, row_dot(projected_rhs, projected_rhs, row, row)));
        if (rhs_scale[static_cast<size_t>(row)] > 0.0) {
            const double inverse = 1.0/rhs_scale[static_cast<size_t>(row)];
            for (int face = 0; face < n_face; ++face)
                projected_rhs[static_cast<size_t>(row)*n_face+face] *= inverse;
            active.push_back(row);
        }
        else iters_out[static_cast<size_t>(row)] = 0;
    }
    std::vector<double> solutions(total, 0.0);
    if (x0)
        for (int row : active)
            for (int face = 0; face < n_face; ++face)
                solutions[static_cast<size_t>(row)*n_face+face] =
                    (*x0)[static_cast<size_t>(row)*n_face+face] /
                    rhs_scale[static_cast<size_t>(row)];
    project_many(solutions);
    std::vector<double> applied, residual(total), zvec(total), pvec(total), apvec;
    apply_system_many(solutions, applied);
    for (size_t i = 0; i < total; ++i)
        residual[i] = projected_rhs[i]-applied[i];
    project_many(residual);

    auto restart_block = [&]() {
        apply_preconditioner_many(residual, zvec, active);
        pvec = zvec;
        return block_gram(residual, zvec, active);
    };
    std::vector<double> gamma = restart_block();
    constexpr int refresh_period = 250;
    // A block space naturally loses rank as correlated state/adjoint columns
    // share converged spectral components.  Use it as a bounded common-Krylov
    // startup, then finish each physical column with the mature true-residual
    // scalar PCG.  This keeps rank deflation from turning into a hard failure
    // and caps the cost if BLAS-3 is not profitable on a particular CPU.
    constexpr int shared_krylov_limit = 500;
    const auto shared_t0 = Clock::now();
    int block_iterations = 0;
    bool block_breakdown = false;
    for (int iteration = 0;
         iteration < std::min(maxit, shared_krylov_limit) && !active.empty();
         ++iteration) {
        block_iterations = iteration+1;
        apply_system_many(pvec, apvec);
        const int dim = static_cast<int>(active.size());
        const std::vector<double> delta = block_gram(pvec, apvec, active);
        std::vector<double> alpha;
        try {
            alpha = symmetric_pseudoinverse_apply(
                delta, gamma, dim, "block PCG search matrix");
        }
        catch (const std::runtime_error& error) {
            if (std::string(error.what()).find("numerical rank is zero") ==
                    std::string::npos)
                throw;
            block_breakdown = true;
            break;
        }
        add_block_product(solutions, pvec, alpha, active, 1.0);
        add_block_product(residual, apvec, alpha, active, -1.0);
        project_many(solutions); project_many(residual);

        bool candidate = ((iteration+1)%refresh_period)==0;
        for (int row : active)
            candidate = candidate || std::sqrt(std::max(
                0.0, row_dot(residual, residual, row, row))) <= tol;
        if (candidate) {
            apply_system_many(solutions, applied);
            for (size_t i = 0; i < total; ++i)
                residual[i] = projected_rhs[i]-applied[i];
            project_many(residual);
            std::vector<int> remaining;
            for (int row : active) {
                const double norm = std::sqrt(std::max(
                    0.0, row_dot(residual, residual, row, row)));
                if (norm <= tol)
                    iters_out[static_cast<size_t>(row)] = iteration+1;
                else remaining.push_back(row);
            }
            active.swap(remaining);
            if (active.empty()) break;
            gamma = restart_block();
            continue;
        }

        apply_preconditioner_many(residual, zvec, active);
        const std::vector<double> gamma_new = block_gram(residual, zvec, active);
        std::vector<double> beta;
        try {
            beta = symmetric_pseudoinverse_apply(
                gamma, gamma_new, dim, "block PCG residual matrix");
        }
        catch (const std::runtime_error& error) {
            if (std::string(error.what()).find("numerical rank is zero") ==
                    std::string::npos)
                throw;
            block_breakdown = true;
            break;
        }
        std::vector<double> next_p = zvec;
        add_block_product(next_p, pvec, beta, active, 1.0);
        pvec.swap(next_p);
        project_many(pvec);
        gamma = gamma_new;
    }
    projection_s = std::chrono::duration<double>(Clock::now()-shared_t0).count();
    if (!active.empty()) {
        // Finish from the shared block iterate.  The exact same A and Jacobi
        // diagonal are reused, and SolveLinearMaterial retains its periodic
        // true-residual replacement and final convergence gate.
        const int remaining_iterations = std::max(1, maxit-block_iterations);
        for (int row : active) {
            const auto offset = static_cast<size_t>(row)*n_face;
            std::vector<double> one_rhs(
                projected_rhs.begin()+offset,
                projected_rhs.begin()+offset+n_face);
            std::vector<double> seed(
                solutions.begin()+offset, solutions.begin()+offset+n_face);
            int scalar_iterations = 0;
            std::vector<double> solution = SolveLinearMaterial(
                m_operatorBIndptr, m_operatorBIndices, m_operatorBData,
                n_face, m_operatorMassI, m_operatorMassJ, m_operatorMassV,
                inv_chi, prec, one_rhs, tol, remaining_iterations,
                scalar_iterations, /*mass_riesz=*/mass_riesz,
                /*symmetric=*/true, &seed, nullptr, nullptr, nullptr, 0);
            iters_out[static_cast<size_t>(row)] =
                block_iterations+scalar_iterations;
            std::copy(solution.begin(), solution.end(),
                      solutions.begin()+offset);
        }
    }
    (void)block_breakdown; // reported indirectly by per-column iteration counts
    for (int row = 0; row < nrhs; ++row)
        for (int face = 0; face < n_face; ++face)
            solutions[static_cast<size_t>(row)*n_face+face] *=
                rhs_scale[static_cast<size_t>(row)];
    return solutions;
}

std::vector<double> RadHACApKChargeGram::ConfiguredFieldFunctionalRows(
    const std::vector<double>& observations,
    const std::vector<double>& weights, int n_rows) const
{
    if (!m_operatorChargeConfigured)
        throw std::runtime_error(
            "ConfiguredFieldFunctionalRows: charge map is not configured");
    if (m_d2)
        throw std::runtime_error(
            "ConfiguredFieldFunctionalRows: planar geometry is not supported");
    if (observations.empty() || observations.size()%3 != 0)
        throw std::invalid_argument(
            "ConfiguredFieldFunctionalRows: observations must have shape (n,3)");
    const int n_observations = static_cast<int>(observations.size()/3);
    if (n_rows < 1 || weights.size() !=
            static_cast<size_t>(n_rows)*observations.size())
        throw std::invalid_argument(
            "ConfiguredFieldFunctionalRows: weights must have shape (n_rows,n,3)");
    for (double value : observations)
        if (!std::isfinite(value))
            throw std::invalid_argument(
                "ConfiguredFieldFunctionalRows: observations must be finite");
    for (double value : weights)
        if (!std::isfinite(value))
            throw std::invalid_argument(
                "ConfiguredFieldFunctionalRows: weights must be finite");
    if (m_hexmode) {
        // An affine HEX BDM1 charge is a physical cubic polynomial in the
        // cell and a physical quadratic polynomial on each quad facet.  Split
        // the cell into the canonical six Kuhn tetrahedra and every facet into
        // two triangles, then reuse the exact analytic TET/TRI field kernels.
        // This is essential near a long pole face: ordinary volume quadrature
        // of the reciprocal 1/r^3 dipole field is not convergent at practical
        // orders when the orbit lies only a fraction of an element away.
        for (int a = 0; a < m_ndof; ++a) {
            const int* e = &m_expo[static_cast<size_t>(3*a)];
            const int degree = e[0]+e[1]+e[2];
            if ((m_kind[static_cast<size_t>(a)] == 0 && degree > 3) ||
                (m_kind[static_cast<size_t>(a)] == 1 &&
                 (degree > 2 || e[2] != 0)))
                throw std::runtime_error(
                    "ConfiguredFieldFunctionalRows: affine HEX supports "
                    "volume degree <= 3 and facet degree <= 2");
        }

        std::vector<double> cell_forms(static_cast<size_t>(m_n_el)*12);
        std::vector<double> cell_inv_jac(static_cast<size_t>(m_n_el));
        std::vector<double> cell_corners(static_cast<size_t>(m_n_el)*24);
        for (int host = 0; host < m_n_el; ++host) {
            const double* nodes = &m_hexNodes[static_cast<size_t>(host)*81];
            double forms[3][4], inv_jac = 0.0;
            if (!HexAffineInverseForms(nodes, forms, inv_jac))
                throw std::runtime_error(
                    "ConfiguredFieldFunctionalRows: exact HEX rows require "
                    "flat affine Q1 geometry");
            std::copy_n(&forms[0][0], 12,
                        &cell_forms[static_cast<size_t>(host)*12]);
            cell_inv_jac[static_cast<size_t>(host)] = inv_jac;
            for (int vertex = 0; vertex < 8; ++vertex) {
                double X[3], J[3][3];
                HexQ2Map(nodes, HEXREF_V[vertex], X, J);
                std::copy_n(X, 3,
                    &cell_corners[(static_cast<size_t>(host)*8+vertex)*3]);
            }
        }

        std::vector<double> face_forms(static_cast<size_t>(m_hex_n_bf)*8);
        std::vector<double> face_inv_jac(static_cast<size_t>(m_hex_n_bf));
        std::vector<double> face_corners(static_cast<size_t>(m_hex_n_bf)*12);
        for (int host = 0; host < m_hex_n_bf; ++host) {
            const double* nodes = &m_quadNodes[static_cast<size_t>(host)*27];
            double forms[2][4], inv_jac = 0.0;
            if (!QuadAffineInverseForms(nodes, forms, inv_jac))
                throw std::runtime_error(
                    "ConfiguredFieldFunctionalRows: exact HEX rows require "
                    "flat affine quad facets");
            std::copy_n(&forms[0][0], 8,
                        &face_forms[static_cast<size_t>(host)*8]);
            face_inv_jac[static_cast<size_t>(host)] = inv_jac;
            for (int vertex = 0; vertex < 4; ++vertex) {
                double X[3], T[3][2];
                QuadQ2Map(nodes, QUADREF_V[vertex], X, T);
                std::copy_n(X, 3,
                    &face_corners[(static_cast<size_t>(host)*4+vertex)*3]);
            }
        }

        std::vector<std::vector<std::pair<int, double>>> charge_by_face(
            static_cast<size_t>(m_operatorNFace));
        for (int a = 0; a < m_ndof; ++a) {
            for (int entry = m_operatorBIndptr[static_cast<size_t>(a)];
                 entry < m_operatorBIndptr[static_cast<size_t>(a)+1]; ++entry) {
                const int face = m_operatorBIndices[static_cast<size_t>(entry)];
                if (face < 0 || face >= m_operatorNFace)
                    throw std::runtime_error(
                        "ConfiguredFieldFunctionalRows: charge-map column is out of range");
                const double value = m_operatorBData[static_cast<size_t>(entry)];
                if (value != 0.0)
                    charge_by_face[static_cast<size_t>(face)].emplace_back(a, value);
            }
        }

        struct WeightedRow {
            int row;
            double value[3];
        };
        std::vector<std::vector<WeightedRow>> rows_by_observation(
            static_cast<size_t>(n_observations));
        for (int row = 0; row < n_rows; ++row)
            for (int observation = 0; observation < n_observations; ++observation) {
                const size_t offset =
                    (static_cast<size_t>(row)*n_observations+observation)*3;
                if (weights[offset] != 0.0 || weights[offset+1] != 0.0 ||
                        weights[offset+2] != 0.0)
                    rows_by_observation[static_cast<size_t>(observation)].push_back(
                        {row,{weights[offset],weights[offset+1],weights[offset+2]}});
            }

        // Convert every reference tensor-product charge mode to the shared
        // physical monomial ordering once.  Observation evaluation below is
        // then grouped by geometry host, so all local modes reuse the same
        // analytic triangle/tetrahedron moments.
        std::vector<double> mode_polynomial(static_cast<size_t>(m_ndof)*20,0.0);
        for(int a=0;a<m_ndof;++a){
            const int host=m_host[static_cast<size_t>(a)];
            const bool volume=m_kind[static_cast<size_t>(a)]==0;
            const double* forms=volume
                ?&cell_forms[static_cast<size_t>(host)*12]
                :&face_forms[static_cast<size_t>(host)*8];
            double* polynomial=&mode_polynomial[static_cast<size_t>(a)*20];
            polynomial[0]=volume
                ?cell_inv_jac[static_cast<size_t>(host)]
                :face_inv_jac[static_cast<size_t>(host)];
            int degree=0;
            const int axes=volume?3:2;
            const int* e=&m_expo[static_cast<size_t>(3*a)];
            for(int axis=0;axis<axes;++axis)
                for(int power=0;power<e[axis];++power)
                    HexPolyMulLinear(polynomial,degree,&forms[4*axis],20);
        }

        constexpr double inv_four_pi =
            0.079577471545947667884441881686257181;
        std::vector<double> charge_rows(
            static_cast<size_t>(n_rows)*m_ndof,0.0);
        ngcore::RegionTaskManager task_manager;
        const int n_hosts=m_n_el+m_hex_n_bf;
        ngcore::ParallelFor(ngcore::IntRange(n_hosts), [&](int combined_host) {
            const bool volume=combined_host<m_n_el;
            const int host=volume?combined_host:combined_host-m_n_el;
            const auto& host_charges=volume
                ?m_cellCharges[static_cast<size_t>(host)]
                :m_faceCharges[static_cast<size_t>(host)];
            if(host_charges.empty())return;
            const int n_basis=volume?20:10;
            auto evaluate_base=[&](const double target[3],double out[20][3]){
                double correction[20][3]={};
                for(int index=0;index<20;++index)
                    for(int k=0;k<3;++k)out[index][k]=0.0;
                const int count=volume?6:2;
                for(int sub=0;sub<count;++sub){
                    double term[20][3]={};
                    if(volume){
                        double vertices[4][3];
                        const double* corners=&cell_corners[static_cast<size_t>(host)*24];
                        for(int local=0;local<4;++local)for(int k=0;k<3;++k)
                            vertices[local][k]=corners[3*HEXREF_TETS[sub][local]+k];
                        rad_hdiv::TetVolFieldCubicBasis(vertices,target,term);
                    }else{
                        double vertices[3][3];
                        const double* corners=&face_corners[static_cast<size_t>(host)*12];
                        for(int local=0;local<3;++local)for(int k=0;k<3;++k)
                            vertices[local][k]=corners[3*QUADREF_TRIS[sub][local]+k];
                        rad_hdiv::QuadTriFieldBasis(vertices,target,term);
                    }
                    for(int index=0;index<n_basis;++index)
                        for(int k=0;k<3;++k){
                        const double next=out[index][k]+term[index][k];
                        correction[index][k]+=
                            std::fabs(out[index][k])>=std::fabs(term[index][k])
                            ?(out[index][k]-next)+term[index][k]
                            :(term[index][k]-next)+out[index][k];
                        out[index][k]=next;
                    }
                }
                for(int index=0;index<n_basis;++index)
                    for(int k=0;k<3;++k)out[index][k]+=correction[index][k];
            };
            std::vector<double> sums(
                static_cast<size_t>(host_charges.size())*n_rows,0.0);
            std::vector<double> corrections(sums.size(),0.0);
            for(int observation=0;observation<n_observations;++observation){
                const double* target=&observations[static_cast<size_t>(observation)*3];
                double total[20][3];evaluate_base(target,total);
                for(size_t image=0;image<m_image_masks.size();++image){
                    // E_{T(sigma)}(x) = T E_sigma(T^-1 x): inverse-map the eval point, forward-map the vector.
                    const int img=(int)image+1;
                    double reflected[3];ImageEvalPoint(img,target,reflected);
                    double term[20][3];
                    evaluate_base(reflected,term);
                    for(int index=0;index<n_basis;++index){
                        double mapped[3];ImageApplyVector(img,term[index],mapped);
                        for(int k=0;k<3;++k)
                            total[index][k]+=m_image_signs[image]*mapped[k];
                    }
                }
                for(const auto& weighted:
                        rows_by_observation[static_cast<size_t>(observation)]){
                    double projected[20];
                    for(int index=0;index<n_basis;++index)
                        projected[index]=weighted.value[0]*total[index][0]
                            +weighted.value[1]*total[index][1]
                            +weighted.value[2]*total[index][2];
                    for(size_t local=0;local<host_charges.size();++local){
                        const int a=host_charges[local];
                        const double* polynomial=
                            &mode_polynomial[static_cast<size_t>(a)*20];
                        double term=0.0;
                        for(int index=0;index<n_basis;++index)
                            term+=polynomial[index]*projected[index];
                        const size_t slot=local*static_cast<size_t>(n_rows)
                            +weighted.row;
                        double& sum=sums[slot];
                        double& correction=corrections[slot];
                        const double next=sum+term;
                        correction+=std::fabs(sum)>=std::fabs(term)
                            ?(sum-next)+term:(term-next)+sum;
                        sum=next;
                    }
                }
            }
            for(size_t local=0;local<host_charges.size();++local){
                const int a=host_charges[local];
                for(int row=0;row<n_rows;++row){
                    const size_t slot=local*static_cast<size_t>(n_rows)+row;
                    charge_rows[static_cast<size_t>(row)*m_ndof+a]=inv_four_pi*
                        (sums[slot]+corrections[slot]);
                }
            }
        });

        std::vector<double> output(
            static_cast<size_t>(n_rows)*m_operatorNFace,0.0);
        ngcore::ParallelFor(
            ngcore::IntRange(n_rows*m_operatorNFace),[&](int linear){
                const int face=linear%m_operatorNFace;
                const int row=linear/m_operatorNFace;
                double sum=0.0,correction=0.0;
                for(const auto& item:charge_by_face[static_cast<size_t>(face)]){
                    const double term=item.second*charge_rows[
                        static_cast<size_t>(row)*m_ndof+item.first];
                    const double next=sum+term;
                    correction+=std::fabs(sum)>=std::fabs(term)
                        ?(sum-next)+term:(term-next)+sum;
                    sum=next;
                }
                output[static_cast<size_t>(linear)]=sum+correction;
            });
        return output;
    }
    if (!m_highorder || m_curved || m_wedgemode)
        throw std::runtime_error(
            "ConfiguredFieldFunctionalRows: exact sparse observation rows currently require flat affine TET geometry");
    for (int a = 0; a < m_ndof; ++a) {
        const int* e = &m_expo[static_cast<size_t>(3*a)];
        const int degree = e[0]+e[1]+e[2];
        if ((m_kind[static_cast<size_t>(a)] == 0 && degree > 1) ||
            (m_kind[static_cast<size_t>(a)] == 1 && degree > 2))
            throw std::runtime_error(
                "ConfiguredFieldFunctionalRows: flat TET charge degree is unsupported");
    }

    // Transpose the configured sparse charge map once.  A broken RT/BDM TET
    // coefficient touches only its cell divergence mode and four face modes;
    // retaining that locality avoids the old O(n_face*n_charge) zero scan.
    std::vector<std::vector<std::pair<int, double>>> charge_by_face(
        static_cast<size_t>(m_operatorNFace));
    for (int a = 0; a < m_ndof; ++a) {
        for (int entry = m_operatorBIndptr[static_cast<size_t>(a)];
             entry < m_operatorBIndptr[static_cast<size_t>(a)+1]; ++entry) {
            const int face = m_operatorBIndices[static_cast<size_t>(entry)];
            if (face < 0 || face >= m_operatorNFace)
                throw std::runtime_error(
                    "ConfiguredFieldFunctionalRows: charge-map column is out of range");
            const double value = m_operatorBData[static_cast<size_t>(entry)];
            if (value != 0.0)
                charge_by_face[static_cast<size_t>(face)].emplace_back(a, value);
        }
    }

    std::vector<double> output(
        static_cast<size_t>(n_rows)*m_operatorNFace, 0.0);
    constexpr double inv_four_pi =
        0.079577471545947667884441881686257181;
    ngcore::RegionTaskManager task_manager;
    ngcore::ParallelFor(ngcore::IntRange(m_operatorNFace), [&](int face) {
        const auto& column = charge_by_face[static_cast<size_t>(face)];
        if (column.empty()) return;
        std::vector<double> volume;
        std::vector<double> surface;
        volume.reserve(column.size()*16);
        surface.reserve(column.size()*22);
        for (const auto& item : column) {
            const int a = item.first;
            const double coefficient = item.second;
            const int host = m_host[static_cast<size_t>(a)];
            const int* e = &m_expo[static_cast<size_t>(3*a)];
            if (m_kind[static_cast<size_t>(a)] == 0) {
                const size_t offset = volume.size();
                volume.resize(offset+16, 0.0);
                double* block = volume.data()+offset;
                std::copy_n(&m_cellV[static_cast<size_t>(host)*12], 12, block);
                const int degree = e[0]+e[1]+e[2];
                if (degree == 0) {
                    block[12] = coefficient;
                } else {
                    const int axis = e[0] ? 0 : (e[1] ? 1 : 2);
                    const double* inv = &m_cellInv[static_cast<size_t>(host)*9];
                    double gradient[3] = {
                        coefficient*inv[3*axis],
                        coefficient*inv[3*axis+1],
                        coefficient*inv[3*axis+2]
                    };
                    block[12] = -(gradient[0]*block[0] +
                                  gradient[1]*block[1] +
                                  gradient[2]*block[2]);
                    for (int k = 0; k < 3; ++k) block[13+k] = gradient[k];
                }
                continue;
            }

            const size_t offset = surface.size();
            surface.resize(offset+22, 0.0);
            double* block = surface.data()+offset;
            const double* vertices =
                &m_faceV[static_cast<size_t>(host)*9];
            std::copy_n(vertices, 9, block);
            const double e1[3] = {
                vertices[3]-vertices[0], vertices[4]-vertices[1],
                vertices[5]-vertices[2]};
            const double e2[3] = {
                vertices[6]-vertices[0], vertices[7]-vertices[1],
                vertices[8]-vertices[2]};
            const double* gi =
                &m_faceGinv[static_cast<size_t>(host)*4];
            double L[2][3];
            for (int k = 0; k < 3; ++k) {
                L[0][k] = gi[0]*e1[k]+gi[1]*e2[k];
                L[1][k] = gi[2]*e1[k]+gi[3]*e2[k];
            }
            const double b[2] = {
                -(L[0][0]*vertices[0]+L[0][1]*vertices[1]+L[0][2]*vertices[2]),
                -(L[1][0]*vertices[0]+L[1][1]*vertices[1]+L[1][2]*vertices[2])
            };
            const int i = e[0], j = e[1];
            if (i == 0 && j == 0) {
                block[9] = coefficient;
            } else if (i+j == 1) {
                const int axis = i ? 0 : 1;
                block[9] = coefficient*b[axis];
                for (int k = 0; k < 3; ++k)
                    block[10+k] = coefficient*L[axis][k];
            } else {
                const int first = i == 2 ? 0 : (j == 2 ? 1 : 0);
                const int second = i == 2 ? 0 : (j == 2 ? 1 : 1);
                block[9] = coefficient*b[first]*b[second];
                for (int k = 0; k < 3; ++k)
                    block[10+k] = coefficient*(
                        b[first]*L[second][k]+b[second]*L[first][k]);
                for (int r = 0; r < 3; ++r)
                    for (int c = 0; c < 3; ++c)
                        block[13+3*r+c] = first == second
                            ? coefficient*L[first][r]*L[first][c]
                            : 0.5*coefficient*(
                                L[first][r]*L[second][c]+L[second][r]*L[first][c]);
            }
        }
        auto evaluator = rad_hdiv::HDivFieldEvaluator::FromTet(
            std::move(volume), std::move(surface),
            m_image_masks, m_image_signs);
        evaluator->SetImageRotations(m_image_rot_angle);
        std::vector<double> field(observations.size());
        evaluator->EvaluateSerial(
            observations.data(), static_cast<size_t>(n_observations),
            field.data(), rad_hdiv::HDivFieldEvaluator::Algorithm::Direct);
        for (int row = 0; row < n_rows; ++row) {
            double value = 0.0;
            double correction = 0.0;
            const size_t offset =
                static_cast<size_t>(row)*observations.size();
            for (size_t index = 0; index < observations.size(); ++index) {
                const double term = weights[offset+index]*field[index];
                const double next = value+term;
                correction += std::fabs(value) >= std::fabs(term)
                    ? (value-next)+term : (term-next)+value;
                value = next;
            }
            output[static_cast<size_t>(row)*m_operatorNFace+face] =
                inv_four_pi*(value+correction);
        }
    });
    return output;
}

std::vector<double>
RadHACApKChargeGram::ConfiguredFieldFunctionalRowsDirectionalDerivative(
    const std::vector<double>& observations,
    const std::vector<double>& weights, int n_rows, int n_modes,
    const std::vector<double>& cell_velocity,
    const std::vector<double>& face_velocity) const
{
    if (!m_operatorChargeConfigured)
        throw std::runtime_error(
            "ConfiguredFieldFunctionalRowsDirectionalDerivative: charge map is not configured");
    if (m_d2 || !m_highorder || m_curved || m_hexmode || m_wedgemode || m_polyCombo)
        throw std::runtime_error(
            "ConfiguredFieldFunctionalRowsDirectionalDerivative: exact analytic rows require flat affine TET geometry");
    if (observations.empty() || observations.size()%3 != 0)
        throw std::invalid_argument(
            "ConfiguredFieldFunctionalRowsDirectionalDerivative: observations must have shape (n,3)");
    const int n_observations=static_cast<int>(observations.size()/3);
    if (n_rows<1 || weights.size()!=static_cast<size_t>(n_rows)*observations.size())
        throw std::invalid_argument(
            "ConfiguredFieldFunctionalRowsDirectionalDerivative: weights must have shape (n_rows,n,3)");
    if (n_modes<1)
        throw std::invalid_argument(
            "ConfiguredFieldFunctionalRowsDirectionalDerivative: n_modes must be positive");
    const int n_cells=static_cast<int>(m_hoCellCharges.size());
    const int n_faces=static_cast<int>(m_hoFaceCharges.size());
    if (cell_velocity.size()!=static_cast<size_t>(n_modes)*n_cells*12)
        throw std::invalid_argument(
            "ConfiguredFieldFunctionalRowsDirectionalDerivative: cell_vertex_velocity must have shape (n_modes,ncell,4,3)");
    if (face_velocity.size()!=static_cast<size_t>(n_modes)*n_faces*9)
        throw std::invalid_argument(
            "ConfiguredFieldFunctionalRowsDirectionalDerivative: face_vertex_velocity must have shape (n_modes,nface,3,3)");
    for(double value:observations)if(!std::isfinite(value))
        throw std::invalid_argument(
            "ConfiguredFieldFunctionalRowsDirectionalDerivative: observations must be finite");
    for(double value:weights)if(!std::isfinite(value))
        throw std::invalid_argument(
            "ConfiguredFieldFunctionalRowsDirectionalDerivative: weights must be finite");
    for(double value:cell_velocity)if(!std::isfinite(value))
        throw std::invalid_argument(
            "ConfiguredFieldFunctionalRowsDirectionalDerivative: cell velocities must be finite");
    for(double value:face_velocity)if(!std::isfinite(value))
        throw std::invalid_argument(
            "ConfiguredFieldFunctionalRowsDirectionalDerivative: face velocities must be finite");
    for(int a=0;a<m_ndof;++a){
        const int* e=&m_expo[static_cast<size_t>(3*a)];
        const int degree=e[0]+e[1]+e[2];
        if((m_kind[static_cast<size_t>(a)]==0&&degree>1)
           ||(m_kind[static_cast<size_t>(a)]==1&&degree>2))
            throw std::runtime_error(
                "ConfiguredFieldFunctionalRowsDirectionalDerivative: flat TET charge degree is unsupported");
    }

    // The Piola B row and its host measure have opposite logarithmic rates.
    // Building the physical polynomial with coefficient (1,rate) therefore
    // differentiates the complete charge-row response, including dB, while
    // retaining the fixed reference monomial exactly.
    std::vector<double> rates(static_cast<size_t>(n_modes)*m_ndof);
    for(int mode=0;mode<n_modes;++mode){
        const auto cell_begin=cell_velocity.begin()+static_cast<size_t>(mode)*n_cells*12;
        const auto face_begin=face_velocity.begin()+static_cast<size_t>(mode)*n_faces*9;
        const auto local=TetChargeMapRowDirectionalRates(
            std::vector<double>(cell_begin,cell_begin+static_cast<size_t>(n_cells)*12),
            std::vector<double>(face_begin,face_begin+static_cast<size_t>(n_faces)*9));
        std::copy(local.begin(),local.end(),rates.begin()+static_cast<size_t>(mode)*m_ndof);
    }

    std::vector<double> charge_rows(
        static_cast<size_t>(n_modes)*n_rows*m_ndof,0.0);
    constexpr double inv_four_pi=0.079577471545947667884441881686257181;
    const double zero_direction[3]={0.0,0.0,0.0};
    ngcore::RegionTaskManager task_manager;
    ngcore::ParallelFor(ngcore::IntRange(n_modes*m_ndof),[&](int linear){
        const int mode=linear/m_ndof;
        const int a=linear-mode*m_ndof;
        const int host=m_host[static_cast<size_t>(a)];
        const int* exponent=&m_expo[static_cast<size_t>(3*a)];
        const double rate=rates[static_cast<size_t>(mode)*m_ndof+a];
        double V4[4][3]{},dV4[4][3]{},rho0=0.0,drho0=0.0,g[3]{},dg[3]{};
        double V3[3][3]{},dV3[3][3]{},sigma0=0.0,dsigma0=0.0;
        double slope[3]{},dslope[3]{},hessian[3][3]{},dhessian[3][3]{};
        const bool volume=m_kind[static_cast<size_t>(a)]==0;
        if(volume){
            const double* vertices=&m_cellV[static_cast<size_t>(host)*12];
            const double* velocity=&cell_velocity[
                (static_cast<size_t>(mode)*n_cells+host)*12];
            for(int i=0;i<4;++i)for(int k=0;k<3;++k){
                V4[i][k]=vertices[3*i+k];dV4[i][k]=velocity[3*i+k];
            }
            const int degree=exponent[0]+exponent[1]+exponent[2];
            if(degree==0){rho0=1.0;drho0=rate;}
            else{
                const int axis=exponent[0]?0:(exponent[1]?1:2);
                const double* inverse=&m_cellInv[static_cast<size_t>(host)*9];
                double dE[3][3]{};
                for(int physical=0;physical<3;++physical)
                    for(int reference=0;reference<3;++reference)
                        dE[physical][reference]=dV4[reference+1][physical]-dV4[0][physical];
                double dinverse[3][3]{};
                for(int i=0;i<3;++i)for(int j=0;j<3;++j)
                    for(int k=0;k<3;++k)for(int l=0;l<3;++l)
                        dinverse[i][j]-=inverse[3*i+k]*dE[k][l]*inverse[3*l+j];
                for(int k=0;k<3;++k){
                    g[k]=inverse[3*axis+k];
                    dg[k]=rate*g[k]+dinverse[axis][k];
                    rho0-=g[k]*V4[0][k];
                    drho0-=dg[k]*V4[0][k]+g[k]*dV4[0][k];
                }
            }
        }else{
            const double* vertices=&m_faceV[static_cast<size_t>(host)*9];
            const double* velocity=&face_velocity[
                (static_cast<size_t>(mode)*n_faces+host)*9];
            for(int i=0;i<3;++i)for(int k=0;k<3;++k){
                V3[i][k]=vertices[3*i+k];dV3[i][k]=velocity[3*i+k];
            }
            double edge[2][3],dedge[2][3];
            for(int i=0;i<2;++i)for(int k=0;k<3;++k){
                edge[i][k]=V3[i+1][k]-V3[0][k];
                dedge[i][k]=dV3[i+1][k]-dV3[0][k];
            }
            double dgram[2][2]{};
            for(int i=0;i<2;++i)for(int j=0;j<2;++j)
                for(int k=0;k<3;++k)
                    dgram[i][j]+=dedge[i][k]*edge[j][k]+edge[i][k]*dedge[j][k];
            const double* inverse=&m_faceGinv[static_cast<size_t>(host)*4];
            double dinverse[2][2]{};
            for(int i=0;i<2;++i)for(int j=0;j<2;++j)
                for(int k=0;k<2;++k)for(int l=0;l<2;++l)
                    dinverse[i][j]-=inverse[2*i+k]*dgram[k][l]*inverse[2*l+j];
            double L[2][3]{},dL[2][3]{},b[2]{},db[2]{};
            for(int i=0;i<2;++i)for(int k=0;k<3;++k){
                for(int j=0;j<2;++j){
                    L[i][k]+=inverse[2*i+j]*edge[j][k];
                    dL[i][k]+=dinverse[i][j]*edge[j][k]
                              +inverse[2*i+j]*dedge[j][k];
                }
                b[i]-=L[i][k]*V3[0][k];
                db[i]-=dL[i][k]*V3[0][k]+L[i][k]*dV3[0][k];
            }
            const int i=exponent[0],j=exponent[1];
            if(i==0&&j==0){sigma0=1.0;dsigma0=rate;}
            else if(i+j==1){
                const int axis=i?0:1;
                sigma0=b[axis];dsigma0=rate*b[axis]+db[axis];
                for(int k=0;k<3;++k){
                    slope[k]=L[axis][k];
                    dslope[k]=rate*L[axis][k]+dL[axis][k];
                }
            }else{
                const int first=i==2?0:(j==2?1:0);
                const int second=i==2?0:(j==2?1:1);
                sigma0=b[first]*b[second];
                dsigma0=rate*sigma0+db[first]*b[second]+b[first]*db[second];
                for(int k=0;k<3;++k){
                    slope[k]=b[first]*L[second][k]+b[second]*L[first][k];
                    dslope[k]=rate*slope[k]+db[first]*L[second][k]
                        +b[first]*dL[second][k]+db[second]*L[first][k]
                        +b[second]*dL[first][k];
                }
                for(int r0=0;r0<3;++r0)for(int c0=0;c0<3;++c0){
                    if(first==second){
                        hessian[r0][c0]=L[first][r0]*L[first][c0];
                        dhessian[r0][c0]=rate*hessian[r0][c0]
                            +dL[first][r0]*L[first][c0]
                            +L[first][r0]*dL[first][c0];
                    }else{
                        hessian[r0][c0]=0.5*(L[first][r0]*L[second][c0]
                                            +L[second][r0]*L[first][c0]);
                        dhessian[r0][c0]=rate*hessian[r0][c0]+0.5*(
                            dL[first][r0]*L[second][c0]
                            +L[first][r0]*dL[second][c0]
                            +dL[second][r0]*L[first][c0]
                            +L[second][r0]*dL[first][c0]);
                    }
                }
            }
        }

        std::vector<double> sums(static_cast<size_t>(n_rows),0.0);
        std::vector<double> corrections(static_cast<size_t>(n_rows),0.0);
        auto evaluate=[&](const double target[3],double derivative[3]){
            double field[3];
            if(volume)rad_hdiv::TetVolFieldLinearDirectional(
                V4,dV4,target,zero_direction,rho0,drho0,g,dg,field,derivative);
            else rad_hdiv::QuadTriFieldDirectional(
                V3,dV3,target,zero_direction,sigma0,dsigma0,slope,dslope,
                hessian,dhessian,field,derivative);
        };
        for(int observation=0;observation<n_observations;++observation){
            const double* target=&observations[static_cast<size_t>(observation)*3];
            double total[3];evaluate(target,total);
            for(size_t image=0;image<m_image_masks.size();++image){
                // E_{T(sigma)}(x) = T E_sigma(T^-1 x): inverse-map the eval point, forward-map the vector.
                const int img=(int)image+1;
                double reflected[3],term[3],mapped[3];
                ImageEvalPoint(img,target,reflected);
                evaluate(reflected,term);
                ImageApplyVector(img,term,mapped);
                for(int k=0;k<3;++k)
                    total[k]+=m_image_signs[image]*mapped[k];
            }
            for(int row=0;row<n_rows;++row){
                const size_t offset=(static_cast<size_t>(row)*n_observations+observation)*3;
                const double term=weights[offset]*total[0]
                    +weights[offset+1]*total[1]+weights[offset+2]*total[2];
                const double next=sums[static_cast<size_t>(row)]+term;
                corrections[static_cast<size_t>(row)]+=std::fabs(sums[static_cast<size_t>(row)])>=std::fabs(term)
                    ?(sums[static_cast<size_t>(row)]-next)+term
                    :(term-next)+sums[static_cast<size_t>(row)];
                sums[static_cast<size_t>(row)]=next;
            }
        }
        for(int row=0;row<n_rows;++row)
            charge_rows[(static_cast<size_t>(mode)*n_rows+row)*m_ndof+a]
                =inv_four_pi*(sums[static_cast<size_t>(row)]+corrections[static_cast<size_t>(row)]);
    });

    std::vector<std::vector<std::pair<int,double>>> charge_by_face(
        static_cast<size_t>(m_operatorNFace));
    for(int a=0;a<m_ndof;++a)
        for(int entry=m_operatorBIndptr[static_cast<size_t>(a)];
            entry<m_operatorBIndptr[static_cast<size_t>(a)+1];++entry){
            const int face=m_operatorBIndices[static_cast<size_t>(entry)];
            if(face<0||face>=m_operatorNFace)
                throw std::runtime_error(
                    "ConfiguredFieldFunctionalRowsDirectionalDerivative: charge-map column is out of range");
            const double coefficient=m_operatorBData[static_cast<size_t>(entry)];
            if(coefficient!=0.0)
                charge_by_face[static_cast<size_t>(face)].emplace_back(a,coefficient);
        }
    std::vector<double> output(
        static_cast<size_t>(n_modes)*n_rows*m_operatorNFace,0.0);
    ngcore::ParallelFor(ngcore::IntRange(n_modes*n_rows*m_operatorNFace),[&](int linear){
        const int face=linear%m_operatorNFace;
        const int outer=linear/m_operatorNFace;
        const int row=outer%n_rows;
        const int mode=outer/n_rows;
        double sum=0.0,correction=0.0;
        for(const auto& item:charge_by_face[static_cast<size_t>(face)]){
            const double term=item.second*charge_rows[
                (static_cast<size_t>(mode)*n_rows+row)*m_ndof+item.first];
            const double next=sum+term;
            correction+=std::fabs(sum)>=std::fabs(term)
                ?(sum-next)+term:(term-next)+sum;
            sum=next;
        }
        output[static_cast<size_t>(linear)]=sum+correction;
    });
    return output;
}

std::vector<double>
RadHACApKChargeGram::ConfiguredFieldValuesShapeDerivative(
    const std::vector<double>& observations,
    const std::vector<double>& magnetization,
    const std::vector<double>& magnetization_jacobian,
    int n_modes,
    const std::vector<double>& cell_velocity,
    const std::vector<double>& face_velocity) const
{
    if (!m_operatorChargeConfigured)
        throw std::runtime_error(
            "ConfiguredFieldValuesShapeDerivative: charge map is not configured");
    if (m_d2 || !m_highorder || m_curved || m_hexmode || m_wedgemode
            || m_polyCombo)
        throw std::runtime_error(
            "ConfiguredFieldValuesShapeDerivative: exact analytic values "
            "require flat affine TET geometry");
    if (observations.empty() || observations.size()%3 != 0)
        throw std::invalid_argument(
            "ConfiguredFieldValuesShapeDerivative: observations must have shape (n,3)");
    const int n_observations = static_cast<int>(observations.size()/3);
    if (n_modes < 1 || magnetization.size()
            != static_cast<size_t>(m_operatorNFace)
            || magnetization_jacobian.size()
            != static_cast<size_t>(n_modes)*m_operatorNFace)
        throw std::invalid_argument(
            "ConfiguredFieldValuesShapeDerivative: state arrays must have "
            "shape (n_face) and (n_modes,n_face)");
    const int n_cells = static_cast<int>(m_hoCellCharges.size());
    const int n_faces = static_cast<int>(m_hoFaceCharges.size());
    if (cell_velocity.size() != static_cast<size_t>(n_modes)*n_cells*12
            || face_velocity.size()
                != static_cast<size_t>(n_modes)*n_faces*9)
        throw std::invalid_argument(
            "ConfiguredFieldValuesShapeDerivative: velocity array shape mismatch");
    for (double value : observations)
        if (!std::isfinite(value))
            throw std::invalid_argument(
                "ConfiguredFieldValuesShapeDerivative: observations must be finite");
    for (double value : magnetization)
        if (!std::isfinite(value))
            throw std::invalid_argument(
                "ConfiguredFieldValuesShapeDerivative: state must be finite");
    for (double value : magnetization_jacobian)
        if (!std::isfinite(value))
            throw std::invalid_argument(
                "ConfiguredFieldValuesShapeDerivative: state Jacobian must be finite");

    std::vector<double> charge(static_cast<size_t>(m_ndof), 0.0);
    std::vector<double> dcharge(
        static_cast<size_t>(n_modes)*m_ndof, 0.0);
    for (int a = 0; a < m_ndof; ++a) {
        for (int entry = m_operatorBIndptr[static_cast<size_t>(a)];
             entry < m_operatorBIndptr[static_cast<size_t>(a)+1]; ++entry) {
            const int face = m_operatorBIndices[static_cast<size_t>(entry)];
            const double coefficient =
                m_operatorBData[static_cast<size_t>(entry)];
            charge[static_cast<size_t>(a)] += coefficient
                * magnetization[static_cast<size_t>(face)];
            for (int mode = 0; mode < n_modes; ++mode)
                dcharge[static_cast<size_t>(mode)*m_ndof+a] += coefficient
                    * magnetization_jacobian[
                        static_cast<size_t>(mode)*m_operatorNFace+face];
        }
    }
    for (int mode = 0; mode < n_modes; ++mode) {
        const auto cell_begin = cell_velocity.begin()
            + static_cast<size_t>(mode)*n_cells*12;
        const auto face_begin = face_velocity.begin()
            + static_cast<size_t>(mode)*n_faces*9;
        const auto rates = TetChargeMapRowDirectionalRates(
            std::vector<double>(
                cell_begin, cell_begin+static_cast<size_t>(n_cells)*12),
            std::vector<double>(
                face_begin, face_begin+static_cast<size_t>(n_faces)*9));
        for (int a = 0; a < m_ndof; ++a)
            dcharge[static_cast<size_t>(mode)*m_ndof+a]
                += rates[static_cast<size_t>(a)]*charge[static_cast<size_t>(a)];
    }

    struct VolumeSource {
        double V[4][3]{}, dV[4][3]{};
        double rho0 = 0.0, drho0 = 0.0;
        double g[3]{}, dg[3]{};
    };
    struct SurfaceSource {
        double V[3][3]{}, dV[3][3]{};
        double sigma0 = 0.0, dsigma0 = 0.0;
        double slope[3]{}, dslope[3]{};
        double hessian[3][3]{}, dhessian[3][3]{};
    };
    std::vector<VolumeSource> volumes(
        static_cast<size_t>(n_modes)*m_ndof);
    std::vector<SurfaceSource> surfaces(
        static_cast<size_t>(n_modes)*m_ndof);
    std::vector<unsigned char> is_volume(static_cast<size_t>(m_ndof));

    for (int a = 0; a < m_ndof; ++a) {
        const int* exponent = &m_expo[static_cast<size_t>(3*a)];
        const int degree = exponent[0]+exponent[1]+exponent[2];
        const bool volume = m_kind[static_cast<size_t>(a)] == 0;
        is_volume[static_cast<size_t>(a)] = volume ? 1 : 0;
        if ((volume && degree > 1) || (!volume && degree > 2))
            throw std::runtime_error(
                "ConfiguredFieldValuesShapeDerivative: unsupported charge degree");
        const int host = m_host[static_cast<size_t>(a)];
        for (int mode = 0; mode < n_modes; ++mode) {
            const size_t linear = static_cast<size_t>(mode)*m_ndof+a;
            const double coefficient = charge[static_cast<size_t>(a)];
            const double dcoefficient = dcharge[linear];
            if (volume) {
                VolumeSource& source = volumes[linear];
                const double* vertices = &m_cellV[static_cast<size_t>(host)*12];
                const double* velocity = &cell_velocity[
                    (static_cast<size_t>(mode)*n_cells+host)*12];
                for (int i = 0; i < 4; ++i)
                    for (int k = 0; k < 3; ++k) {
                        source.V[i][k] = vertices[3*i+k];
                        source.dV[i][k] = velocity[3*i+k];
                    }
                if (degree == 0) {
                    source.rho0 = coefficient;
                    source.drho0 = dcoefficient;
                } else {
                    const int axis = exponent[0] ? 0 : (exponent[1] ? 1 : 2);
                    const double* inverse = &m_cellInv[static_cast<size_t>(host)*9];
                    double dE[3][3]{}, dinverse[3][3]{};
                    for (int physical = 0; physical < 3; ++physical)
                        for (int reference = 0; reference < 3; ++reference)
                            dE[physical][reference] =
                                source.dV[reference+1][physical]
                                - source.dV[0][physical];
                    for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j)
                            for (int k = 0; k < 3; ++k)
                                for (int l = 0; l < 3; ++l)
                                    dinverse[i][j] -= inverse[3*i+k]
                                        * dE[k][l]*inverse[3*l+j];
                    for (int k = 0; k < 3; ++k) {
                        source.g[k] = coefficient*inverse[3*axis+k];
                        source.dg[k] = dcoefficient*inverse[3*axis+k]
                            + coefficient*dinverse[axis][k];
                        source.rho0 -= source.g[k]*source.V[0][k];
                        source.drho0 -= source.dg[k]*source.V[0][k]
                            + source.g[k]*source.dV[0][k];
                    }
                }
            } else {
                SurfaceSource& source = surfaces[linear];
                const double* vertices = &m_faceV[static_cast<size_t>(host)*9];
                const double* velocity = &face_velocity[
                    (static_cast<size_t>(mode)*n_faces+host)*9];
                for (int i = 0; i < 3; ++i)
                    for (int k = 0; k < 3; ++k) {
                        source.V[i][k] = vertices[3*i+k];
                        source.dV[i][k] = velocity[3*i+k];
                    }
                double edge[2][3], dedge[2][3], dgram[2][2]{};
                for (int i = 0; i < 2; ++i)
                    for (int k = 0; k < 3; ++k) {
                        edge[i][k] = source.V[i+1][k]-source.V[0][k];
                        dedge[i][k] = source.dV[i+1][k]-source.dV[0][k];
                    }
                for (int i = 0; i < 2; ++i)
                    for (int j = 0; j < 2; ++j)
                        for (int k = 0; k < 3; ++k)
                            dgram[i][j] += dedge[i][k]*edge[j][k]
                                + edge[i][k]*dedge[j][k];
                const double* inverse = &m_faceGinv[static_cast<size_t>(host)*4];
                double dinverse[2][2]{}, L[2][3]{}, dL[2][3]{}, b[2]{}, db[2]{};
                for (int i = 0; i < 2; ++i)
                    for (int j = 0; j < 2; ++j)
                        for (int k = 0; k < 2; ++k)
                            for (int l = 0; l < 2; ++l)
                                dinverse[i][j] -= inverse[2*i+k]
                                    * dgram[k][l]*inverse[2*l+j];
                for (int i = 0; i < 2; ++i)
                    for (int k = 0; k < 3; ++k) {
                        for (int j = 0; j < 2; ++j) {
                            L[i][k] += inverse[2*i+j]*edge[j][k];
                            dL[i][k] += dinverse[i][j]*edge[j][k]
                                + inverse[2*i+j]*dedge[j][k];
                        }
                        b[i] -= L[i][k]*source.V[0][k];
                        db[i] -= dL[i][k]*source.V[0][k]
                            + L[i][k]*source.dV[0][k];
                    }
                const int i = exponent[0], j = exponent[1];
                if (i == 0 && j == 0) {
                    source.sigma0 = coefficient;
                    source.dsigma0 = dcoefficient;
                } else if (i+j == 1) {
                    const int axis = i ? 0 : 1;
                    source.sigma0 = coefficient*b[axis];
                    source.dsigma0 = dcoefficient*b[axis]
                        + coefficient*db[axis];
                    for (int k = 0; k < 3; ++k) {
                        source.slope[k] = coefficient*L[axis][k];
                        source.dslope[k] = dcoefficient*L[axis][k]
                            + coefficient*dL[axis][k];
                    }
                } else {
                    const int first = i == 2 ? 0 : (j == 2 ? 1 : 0);
                    const int second = i == 2 ? 0 : (j == 2 ? 1 : 1);
                    source.sigma0 = coefficient*b[first]*b[second];
                    source.dsigma0 = dcoefficient*b[first]*b[second]
                        + coefficient*(db[first]*b[second]
                            + b[first]*db[second]);
                    for (int k = 0; k < 3; ++k) {
                        const double value = b[first]*L[second][k]
                            + b[second]*L[first][k];
                        const double derivative = db[first]*L[second][k]
                            + b[first]*dL[second][k]
                            + db[second]*L[first][k]
                            + b[second]*dL[first][k];
                        source.slope[k] = coefficient*value;
                        source.dslope[k] = dcoefficient*value
                            + coefficient*derivative;
                    }
                    for (int r = 0; r < 3; ++r)
                        for (int c = 0; c < 3; ++c) {
                            const double value = first == second
                                ? L[first][r]*L[first][c]
                                : 0.5*(L[first][r]*L[second][c]
                                    + L[second][r]*L[first][c]);
                            const double derivative = first == second
                                ? dL[first][r]*L[first][c]
                                    + L[first][r]*dL[first][c]
                                : 0.5*(dL[first][r]*L[second][c]
                                    + L[first][r]*dL[second][c]
                                    + dL[second][r]*L[first][c]
                                    + L[second][r]*dL[first][c]);
                            source.hessian[r][c] = coefficient*value;
                            source.dhessian[r][c] = dcoefficient*value
                                + coefficient*derivative;
                        }
                }
            }
        }
    }

    constexpr double inv_four_pi =
        0.079577471545947667884441881686257181;
    const double zero_direction[3] = {0.0, 0.0, 0.0};
    std::vector<double> output(
        static_cast<size_t>(n_modes)*observations.size(), 0.0);
    ngcore::RegionTaskManager task_manager;
    ngcore::ParallelFor(
        ngcore::IntRange(n_modes*n_observations), [&](int linear) {
            const int mode = linear/n_observations;
            const int observation = linear-mode*n_observations;
            const double* target = &observations[
                static_cast<size_t>(observation)*3];
            double sum[3]{}, correction[3]{};
            auto add = [&](const double value[3], double scale) {
                for (int k = 0; k < 3; ++k) {
                    const double term = scale*value[k];
                    const double next = sum[k]+term;
                    correction[k] += std::fabs(sum[k]) >= std::fabs(term)
                        ? (sum[k]-next)+term : (term-next)+sum[k];
                    sum[k] = next;
                }
            };
            for (int a = 0; a < m_ndof; ++a) {
                const size_t source_index =
                    static_cast<size_t>(mode)*m_ndof+a;
                auto evaluate = [&](const double point[3], double value[3]) {
                    double base_value[3];
                    if (is_volume[static_cast<size_t>(a)]) {
                        const VolumeSource& source = volumes[source_index];
                        rad_hdiv::TetVolFieldLinearDirectional(
                            source.V, source.dV, point, zero_direction,
                            source.rho0, source.drho0,
                            source.g, source.dg, base_value, value);
                    } else {
                        const SurfaceSource& source = surfaces[source_index];
                        rad_hdiv::QuadTriFieldDirectional(
                            source.V, source.dV, point, zero_direction,
                            source.sigma0, source.dsigma0,
                            source.slope, source.dslope,
                            source.hessian, source.dhessian,
                            base_value, value);
                    }
                };
                double value[3];
                evaluate(target, value);
                add(value, 1.0);
                for (size_t image = 0; image < m_image_masks.size(); ++image) {
                    const int img = static_cast<int>(image)+1;
                    double reflected[3], term[3], mapped[3];
                    ImageEvalPoint(img, target, reflected);
                    evaluate(reflected, term);
                    ImageApplyVector(img, term, mapped);
                    add(mapped, m_image_signs[image]);
                }
            }
            const size_t offset =
                (static_cast<size_t>(mode)*n_observations+observation)*3;
            for (int k = 0; k < 3; ++k)
                output[offset+k] = inv_four_pi*(sum[k]+correction[k]);
        });
    return output;
}

std::shared_ptr<rad_hdiv::HDivFieldEvaluator>
RadHACApKChargeGram::CreateConfiguredFieldEvaluator(
    const std::vector<double>& magnetization,
    const rad_hdiv::FieldEvaluatorOptions& options) const
{
    if (!m_operatorChargeConfigured)
        throw std::runtime_error("CreateConfiguredFieldEvaluator: charge map is not configured");
    if (static_cast<int>(magnetization.size()) != m_operatorNFace)
        throw std::runtime_error("CreateConfiguredFieldEvaluator: magnetization size mismatch");
    if (m_d2)
        throw std::runtime_error("CreateConfiguredFieldEvaluator: use the planar 2D field evaluator");

    std::vector<double> charge(static_cast<size_t>(m_ndof), 0.0);
    for (int a = 0; a < m_ndof; ++a) {
        double value = 0.0;
        for (int k = m_operatorBIndptr[static_cast<size_t>(a)];
             k < m_operatorBIndptr[static_cast<size_t>(a) + 1]; ++k)
            value += m_operatorBData[static_cast<size_t>(k)] *
                     magnetization[static_cast<size_t>(m_operatorBIndices[static_cast<size_t>(k)])];
        charge[static_cast<size_t>(a)] = value;
    }

    // Flat TET BDM1/BDM2: convert reference monomials to physical polynomials once and retain the exact
    // analytic volume/triangle field kernels at every target, including near-surface targets.
    bool analytic_tet = m_highorder && !m_curved && !m_hexmode && !m_wedgemode;
    if (analytic_tet) {
        for (int a = 0; a < m_ndof; ++a) {
            const int* e = &m_expo[static_cast<size_t>(3*a)];
            const int degree = e[0] + e[1] + e[2];
            if ((m_kind[a] == 0 && degree > 1) || (m_kind[a] == 1 && degree > 2)) {
                analytic_tet = false;
                break;
            }
        }
    }
    if (analytic_tet) {
        const int n_cells = m_n_el;
        const int n_faces = static_cast<int>(m_faceV.size()/9);
        std::vector<double> volume(static_cast<size_t>(n_cells)*16, 0.0);
        std::vector<double> surface(static_cast<size_t>(n_faces)*22, 0.0);
        for (int c = 0; c < n_cells; ++c)
            std::copy_n(&m_cellV[static_cast<size_t>(c)*12], 12,
                        &volume[static_cast<size_t>(c)*16]);
        for (int f = 0; f < n_faces; ++f)
            std::copy_n(&m_faceV[static_cast<size_t>(f)*9], 9,
                        &surface[static_cast<size_t>(f)*22]);

        for (int a = 0; a < m_ndof; ++a) {
            const double coefficient = charge[static_cast<size_t>(a)];
            if (coefficient == 0.0) continue;
            const int host = m_host[static_cast<size_t>(a)];
            const int* e = &m_expo[static_cast<size_t>(3*a)];
            if (m_kind[static_cast<size_t>(a)] == 0) {
                double* out = &volume[static_cast<size_t>(host)*16];
                const int degree = e[0] + e[1] + e[2];
                if (degree == 0) {
                    out[12] += coefficient;
                } else {
                    const int axis = e[0] ? 0 : (e[1] ? 1 : 2);
                    const double* inv = &m_cellInv[static_cast<size_t>(host)*9];
                    double gradient[3] = {
                        coefficient*inv[3*axis],
                        coefficient*inv[3*axis + 1],
                        coefficient*inv[3*axis + 2]
                    };
                    out[12] -= gradient[0]*out[0] + gradient[1]*out[1] + gradient[2]*out[2];
                    for (int k = 0; k < 3; ++k) out[13+k] += gradient[k];
                }
                continue;
            }

            double* out = &surface[static_cast<size_t>(host)*22];
            const double* vertices = &m_faceV[static_cast<size_t>(host)*9];
            const double e1[3] = {vertices[3]-vertices[0], vertices[4]-vertices[1], vertices[5]-vertices[2]};
            const double e2[3] = {vertices[6]-vertices[0], vertices[7]-vertices[1], vertices[8]-vertices[2]};
            const double* gi = &m_faceGinv[static_cast<size_t>(host)*4];
            double L[2][3];
            for (int k = 0; k < 3; ++k) {
                L[0][k] = gi[0]*e1[k] + gi[1]*e2[k];
                L[1][k] = gi[2]*e1[k] + gi[3]*e2[k];
            }
            const double b[2] = {
                -(L[0][0]*vertices[0] + L[0][1]*vertices[1] + L[0][2]*vertices[2]),
                -(L[1][0]*vertices[0] + L[1][1]*vertices[1] + L[1][2]*vertices[2])
            };
            const int i = e[0], j = e[1];
            if (i == 0 && j == 0) {
                out[9] += coefficient;
            } else if (i + j == 1) {
                const int axis = i ? 0 : 1;
                out[9] += coefficient*b[axis];
                for (int k = 0; k < 3; ++k) out[10+k] += coefficient*L[axis][k];
            } else {
                const int first = (i == 2) ? 0 : (j == 2 ? 1 : 0);
                const int second = (i == 2) ? 0 : (j == 2 ? 1 : 1);
                out[9] += coefficient*b[first]*b[second];
                for (int k = 0; k < 3; ++k)
                    out[10+k] += coefficient*(b[first]*L[second][k] + b[second]*L[first][k]);
                for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) {
                    const double h = (first == second)
                        ? coefficient*L[first][r]*L[first][c]
                        : 0.5*coefficient*(L[first][r]*L[second][c] + L[second][r]*L[first][c]);
                    out[13 + 3*r + c] += h;
                }
            }
        }
        auto evaluator = rad_hdiv::HDivFieldEvaluator::FromTet(
            std::move(volume), std::move(surface), m_image_masks, m_image_signs, options);
        evaluator->SetImageRotations(m_image_rot_angle);
        return evaluator;
    }

    // Straight affine HEX BDM1: retain an exact near-field representation.
    // The reference volume charge has total degree <= 3 and each reference
    // facet charge has degree <= 2.  Under an affine map these remain physical
    // linear/quadratic polynomials, so the canonical 6-TET / 2-TRI split can
    // reuse the analytic kernels.  The old fixed tensor cloud is adequate for
    // smooth Gram far blocks but is not a convergent near-pole field evaluator.
    bool analytic_hex = m_hexmode;
    bool configured_hex_rt0 = m_hexmode;
    for (int a = 0; analytic_hex && a < m_ndof; ++a) {
        const int* e = &m_expo[static_cast<size_t>(3*a)];
        const int degree = e[0]+e[1]+e[2];
        configured_hex_rt0 = configured_hex_rt0 && degree == 0;
        analytic_hex = m_kind[static_cast<size_t>(a)] == 0
            ? degree <= 3
            : (degree <= 2 && e[2] == 0);
    }
    if (analytic_hex) {
        std::vector<double> volume;
        std::vector<double> surface;
        volume.reserve(static_cast<size_t>(m_n_el)*6*32);
        surface.reserve(static_cast<size_t>(m_hex_n_bf)*2*22);
        for (int host = 0; analytic_hex && host < m_n_el; ++host) {
            const double* nodes = &m_hexNodes[static_cast<size_t>(host)*81];
            double forms[3][4], inv_jac = 0.0;
            if (!HexAffineInverseForms(nodes, forms, inv_jac)) {
                analytic_hex = false;
                break;
            }
            double corners[8][3];
            for (int vertex = 0; vertex < 8; ++vertex) {
                double J[3][3];
                HexQ2Map(nodes, HEXREF_V[vertex], corners[vertex], J);
            }
            double polynomial[20] = {};
            for (int charge_id : m_cellCharges[static_cast<size_t>(host)]) {
                const double factor = charge[static_cast<size_t>(charge_id)]*inv_jac;
                const int* e = &m_expo[static_cast<size_t>(3*charge_id)];
                double mode[20] = {};
                mode[0] = factor;
                int degree = 0;
                for (int axis = 0; axis < 3; ++axis)
                    for (int power = 0; power < e[axis]; ++power)
                        HexPolyMulLinear(mode,degree,forms[axis],20);
                for (int index = 0; index < 20; ++index)
                    polynomial[index] += mode[index];
            }
            bool empty = true;
            for (double value : polynomial) empty = empty && value == 0.0;
            if (empty) continue;
            for (int sub = 0; sub < 6; ++sub) {
                const size_t offset = volume.size();
                volume.resize(offset+32, 0.0);
                for (int local = 0; local < 4; ++local)
                    for (int axis = 0; axis < 3; ++axis)
                        volume[offset+3*local+axis] =
                            corners[HEXREF_TETS[sub][local]][axis];
                std::copy_n(polynomial,20,volume.data()+offset+12);
            }
        }

        for (int host = 0; analytic_hex && host < m_hex_n_bf; ++host) {
            const double* nodes = &m_quadNodes[static_cast<size_t>(host)*27];
            double forms[2][4], inv_jac = 0.0;
            if (!QuadAffineInverseForms(nodes, forms, inv_jac)) {
                analytic_hex = false;
                break;
            }
            double corners[4][3];
            for (int vertex = 0; vertex < 4; ++vertex) {
                double T[3][2];
                QuadQ2Map(nodes, QUADREF_V[vertex], corners[vertex], T);
            }
            double sigma0 = 0.0, slope[3] = {0.0, 0.0, 0.0};
            double hessian[3][3] = {};
            for (int charge_id : m_faceCharges[static_cast<size_t>(host)]) {
                const double factor = charge[static_cast<size_t>(charge_id)]*inv_jac;
                const int* e = &m_expo[static_cast<size_t>(3*charge_id)];
                const int i = e[0], j = e[1];
                if (i == 0 && j == 0) {
                    sigma0 += factor;
                } else if (i+j == 1) {
                    const int axis = i ? 0 : 1;
                    sigma0 += factor*forms[axis][0];
                    for (int k = 0; k < 3; ++k)
                        slope[k] += factor*forms[axis][k+1];
                } else {
                    const int first = i == 2 ? 0 : (j == 2 ? 1 : 0);
                    const int second = i == 2 ? 0 : (j == 2 ? 1 : 1);
                    const double* f = forms[first];
                    const double* s = forms[second];
                    sigma0 += factor*f[0]*s[0];
                    for (int k = 0; k < 3; ++k)
                        slope[k] += factor*(f[0]*s[k+1]+s[0]*f[k+1]);
                    for (int r = 0; r < 3; ++r)
                        for (int c = 0; c < 3; ++c)
                            hessian[r][c] += first == second
                                ? factor*f[r+1]*f[c+1]
                                : 0.5*factor*(f[r+1]*s[c+1]+s[r+1]*f[c+1]);
                }
            }
            bool empty = sigma0 == 0.0;
            for (int k = 0; k < 3; ++k) empty = empty && slope[k] == 0.0;
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                    empty = empty && hessian[r][c] == 0.0;
            if (empty) continue;
            for (int sub = 0; sub < 2; ++sub) {
                const size_t offset = surface.size();
                surface.resize(offset+22, 0.0);
                for (int local = 0; local < 3; ++local)
                    for (int axis = 0; axis < 3; ++axis)
                        surface[offset+3*local+axis] =
                            corners[QUADREF_TRIS[sub][local]][axis];
                surface[offset+9] = sigma0;
                for (int k = 0; k < 3; ++k)
                    surface[offset+10+k] = slope[k];
                for (int r = 0; r < 3; ++r)
                    for (int c = 0; c < 3; ++c)
                        surface[offset+13+3*r+c] = hessian[r][c];
            }
        }
        if (analytic_hex) {
            if (!volume.empty() || !surface.empty()) {
                auto evaluator = rad_hdiv::HDivFieldEvaluator::FromPolynomialTet(
                    std::move(volume), std::move(surface),
                    m_image_masks, m_image_signs, options);
                evaluator->SetImageRotations(m_image_rot_angle);
                return evaluator;
            }
            auto empty_evaluator = rad_hdiv::HDivFieldEvaluator::FromCloud(
                {0.0, 0.0, 0.0}, {0.0},
                m_image_masks, m_image_signs, options);
            empty_evaluator->SetImageRotations(m_image_rot_angle);
            return empty_evaluator;
        }
    }
    if (configured_hex_rt0)
        throw std::runtime_error(
            "CreateConfiguredFieldEvaluator: HEX RT0 direct field requires "
            "an affine Q1 geometry; non-affine/warped HEX must be refined or "
            "promoted to the body-fitted BDM field stage");

    // Curved HEX/WEDGE use their shared host quadrature directly.  Co-located charge modes are combined at
    // each physical point, so the source count scales with elements, not with element modes.
    if (m_hexmode || m_wedgemode) {
        std::vector<double> xyz;
        std::vector<double> strength;
        auto append_cloud = [&](const HexQuadCloud& cloud, const std::vector<int>& group) {
            for (size_t q = 0; q < cloud.wgeo.size(); ++q) {
                const double* xi = &cloud.xi[3*q];
                double density = 0.0;
                for (int charge_id : group)
                    density += charge[static_cast<size_t>(charge_id)]*HexMonoEval(charge_id, xi);
                xyz.push_back(cloud.pts[3*q]);
                xyz.push_back(cloud.pts[3*q + 1]);
                xyz.push_back(cloud.pts[3*q + 2]);
                strength.push_back(cloud.wgeo[q]*density);
            }
        };
        if (m_hexmode) {
            // A tensor-product rule is invariant under every signed permutation of the HEX reference
            // axes.  That is essential for the IMA field contract: the explicit lower-half element and
            // the reflected upper-half element must materialize the same physical point sources.  The
            // Kuhn 6-tet / 2-tri split used by the near Gram integration is deliberately not reused here,
            // because its diagonal selects an orientation and produced O(1e-4) full-vs-image field drift.
            for (int host = 0; host < m_n_el; ++host) {
                const double* nodes = &m_hexNodes[static_cast<size_t>(host)*81];
                HexQuadCloud cloud;
                const size_t ng = m_glOut.size();
                cloud.pts.reserve(3*ng*ng*ng);
                cloud.xi.reserve(3*ng*ng*ng);
                cloud.wgeo.reserve(ng*ng*ng);
                for (size_t iz = 0; iz < ng; ++iz)
                    for (size_t iy = 0; iy < ng; ++iy)
                        for (size_t ix = 0; ix < ng; ++ix) {
                            const double xi[3] = {m_glOut[ix], m_glOut[iy], m_glOut[iz]};
                            double X[3], J[3][3];
                            HexQ2Map(nodes, xi, X, J);
                            cloud.pts.insert(cloud.pts.end(), X, X+3);
                            cloud.xi.insert(cloud.xi.end(), xi, xi+3);
                            cloud.wgeo.push_back(m_gwOut[ix]*m_gwOut[iy]*m_gwOut[iz]);
                        }
                append_cloud(cloud, m_cellCharges[static_cast<size_t>(host)]);
            }
            for (int host = 0; host < m_hex_n_bf; ++host) {
                const double* nodes = &m_quadNodes[static_cast<size_t>(host)*27];
                HexQuadCloud cloud;
                const size_t ng = m_glOut.size();
                cloud.pts.reserve(3*ng*ng);
                cloud.xi.reserve(3*ng*ng);
                cloud.wgeo.reserve(ng*ng);
                for (size_t iv = 0; iv < ng; ++iv)
                    for (size_t iu = 0; iu < ng; ++iu) {
                        const double uv[2] = {m_glOut[iu], m_glOut[iv]};
                        double X[3], T[3][2];
                        QuadQ2Map(nodes, uv, X, T);
                        cloud.pts.insert(cloud.pts.end(), X, X+3);
                        cloud.xi.push_back(uv[0]);
                        cloud.xi.push_back(uv[1]);
                        cloud.xi.push_back(0.0);
                        cloud.wgeo.push_back(m_gwOut[iu]*m_gwOut[iv]);
                    }
                append_cloud(cloud, m_faceCharges[static_cast<size_t>(host)]);
            }
        } else {
            // Prism product quadrature is invariant under triangle-vertex permutations and axial
            // reflection.  Keep the 3-tet decomposition for singular Gram integration only.
            for (int host = 0; host < m_n_el; ++host) {
                const double* nodes = &m_wCellNodes[static_cast<size_t>(host)*54];
                HexQuadCloud cloud;
                const size_t nt = m_wFieldTriW.size(), ng = m_glOut.size();
                cloud.pts.reserve(3*nt*ng);
                cloud.xi.reserve(3*nt*ng);
                cloud.wgeo.reserve(nt*ng);
                for (size_t iw = 0; iw < ng; ++iw)
                    for (size_t iq = 0; iq < nt; ++iq) {
                        const double xi[3] = {
                            m_wFieldTriP[2*iq], m_wFieldTriP[2*iq + 1], m_glOut[iw]
                        };
                        double X[3];
                        WedgeQ2MapX(nodes, xi, X);
                        cloud.pts.insert(cloud.pts.end(), X, X+3);
                        cloud.xi.insert(cloud.xi.end(), xi, xi+3);
                        cloud.wgeo.push_back(m_wFieldTriW[iq]*m_gwOut[iw]);
                    }
                append_cloud(cloud, m_cellCharges[static_cast<size_t>(host)]);
            }
            for (int host = 0; host < m_wedge_n_bf; ++host) {
                const int face_type = m_wFaceType[static_cast<size_t>(host)];
                const double* nodes = &m_wFaceNodes[static_cast<size_t>(host)*27];
                HexQuadCloud cloud;
                if (face_type == 0) {
                    const size_t nt = m_wFieldTriW.size();
                    cloud.pts.reserve(3*nt);
                    cloud.xi.reserve(3*nt);
                    cloud.wgeo.reserve(nt);
                    for (size_t iq = 0; iq < nt; ++iq) {
                        const double uv[2] = {m_wFieldTriP[2*iq], m_wFieldTriP[2*iq + 1]};
                        double X[3];
                        TriSurfMap(nodes, uv, X);
                        cloud.pts.insert(cloud.pts.end(), X, X+3);
                        cloud.xi.push_back(uv[0]);
                        cloud.xi.push_back(uv[1]);
                        cloud.xi.push_back(0.0);
                        cloud.wgeo.push_back(m_wFieldTriW[iq]);
                    }
                } else {
                    const size_t ng = m_glOut.size();
                    cloud.pts.reserve(3*ng*ng);
                    cloud.xi.reserve(3*ng*ng);
                    cloud.wgeo.reserve(ng*ng);
                    for (size_t iv = 0; iv < ng; ++iv)
                        for (size_t iu = 0; iu < ng; ++iu) {
                            const double uv[2] = {m_glOut[iu], m_glOut[iv]};
                            double X[3];
                            QuadQ2MapX(nodes, uv, X);
                            cloud.pts.insert(cloud.pts.end(), X, X+3);
                            cloud.xi.push_back(uv[0]);
                            cloud.xi.push_back(uv[1]);
                            cloud.xi.push_back(0.0);
                            cloud.wgeo.push_back(m_gwOut[iu]*m_gwOut[iv]);
                        }
                }
                append_cloud(cloud, m_faceCharges[static_cast<size_t>(host)]);
            }
        }
        if (strength.empty()) {
            xyz = {0.0, 0.0, 0.0};
            strength = {0.0};
        }
        auto evaluator = rad_hdiv::HDivFieldEvaluator::FromCloud(
            std::move(xyz), std::move(strength), m_image_masks, m_image_signs, options);
        evaluator->SetImageRotations(m_image_rot_angle);
        return evaluator;
    }

    // Curved TET: retain P2 geometry and the combined BDM1/BDM2 reference
    // polynomial per host.  The persistent evaluator integrates exact element
    // leaves at observation time and uses prebuilt moments only for accepted
    // tree nodes.  This avoids both Python source packing and the near-field
    // error of freezing a low-order quadrature cloud at solve time.
    const int n_cells = static_cast<int>(m_cellNodes.size()/30);
    const int n_faces = static_cast<int>(m_faceNodes.size()/18);
    std::vector<double> volume(static_cast<size_t>(n_cells)*34, 0.0);
    std::vector<double> surface(static_cast<size_t>(n_faces)*24, 0.0);
    for (int host = 0; host < n_cells; ++host)
        std::copy_n(&m_cellNodes[static_cast<size_t>(host)*30], 30,
                    &volume[static_cast<size_t>(host)*34]);
    for (int host = 0; host < n_faces; ++host)
        std::copy_n(&m_faceNodes[static_cast<size_t>(host)*18], 18,
                    &surface[static_cast<size_t>(host)*24]);
    for (int a = 0; a < m_ndof; ++a) {
        const double coefficient = charge[static_cast<size_t>(a)];
        if (coefficient == 0.0) continue;
        const int host = m_host[static_cast<size_t>(a)];
        const int* e = &m_expo[static_cast<size_t>(3*a)];
        if (m_kind[static_cast<size_t>(a)] == 0) {
            int local = 0;
            if (e[0] == 1 && e[1] == 0 && e[2] == 0) local = 1;
            else if (e[0] == 0 && e[1] == 1 && e[2] == 0) local = 2;
            else if (e[0] == 0 && e[1] == 0 && e[2] == 1) local = 3;
            else if (e[0] != 0 || e[1] != 0 || e[2] != 0)
                throw std::runtime_error("CreateConfiguredFieldEvaluator: curved volume charge degree > 1");
            volume[static_cast<size_t>(host)*34 + 30 + local] += coefficient;
        } else {
            int local = -1;
            if (e[0] == 0 && e[1] == 0) local = 0;
            else if (e[0] == 0 && e[1] == 1) local = 1;
            else if (e[0] == 0 && e[1] == 2) local = 2;
            else if (e[0] == 1 && e[1] == 0) local = 3;
            else if (e[0] == 1 && e[1] == 1) local = 4;
            else if (e[0] == 2 && e[1] == 0) local = 5;
            if (local < 0)
                throw std::runtime_error("CreateConfiguredFieldEvaluator: curved surface charge degree > 2");
            surface[static_cast<size_t>(host)*24 + 18 + local] += coefficient;
        }
    }
    auto evaluator = rad_hdiv::HDivFieldEvaluator::FromCurvedTet(
        std::move(volume), std::move(surface), m_gl, m_gw,
        m_image_masks, m_image_signs, options);
    evaluator->SetImageRotations(m_image_rot_angle);
    return evaluator;
}

std::shared_ptr<rad_planar_charges::PlanarFieldEvaluator>
RadHACApKChargeGram::CreateConfiguredPlanarFieldEvaluator(
    const std::vector<double>& magnetization) const
{
    if (!m_operatorChargeConfigured)
        throw std::runtime_error("CreateConfiguredPlanarFieldEvaluator: charge map is not configured");
    if (!m_d2)
        throw std::runtime_error("CreateConfiguredPlanarFieldEvaluator: Gram is not planar");
    if (static_cast<int>(magnetization.size()) != m_operatorNFace)
        throw std::runtime_error("CreateConfiguredPlanarFieldEvaluator: magnetization size mismatch");

    std::vector<double> charge(static_cast<size_t>(m_ndof), 0.0);
    for (int a = 0; a < m_ndof; ++a)
        for (int k = m_operatorBIndptr[static_cast<size_t>(a)];
             k < m_operatorBIndptr[static_cast<size_t>(a) + 1]; ++k)
            charge[static_cast<size_t>(a)] += m_operatorBData[static_cast<size_t>(k)]
                * magnetization[static_cast<size_t>(m_operatorBIndices[static_cast<size_t>(k)])];

    std::vector<double> positions;
    std::vector<double> strengths;
    auto append_cell_point = [&](int host, const double xi[2], double weight) {
        double density = 0.0;
        for (int a : m_cellCharges[static_cast<size_t>(host)])
            density += charge[static_cast<size_t>(a)]
                     * D2MonoCell(&m_expo[static_cast<size_t>(3*a)], xi);
        double X[2];
        const double* map = &m_d2CellMap[static_cast<size_t>(host)*m_d2CellMapStride];
        D2CellMap(m_d2CellType[static_cast<size_t>(host)], map, xi, X);
        positions.push_back(X[0]);
        positions.push_back(X[1]);
        strengths.push_back(weight*density);
    };

    for (int host = 0; host < m_n_el; ++host) {
        if (m_d2CellType[static_cast<size_t>(host)] == 1) {
            // Tensor Gauss is invariant under the signed reference-axis permutations used by IMA.
            for (size_t j = 0; j < m_glIn.size(); ++j)
                for (size_t i = 0; i < m_glIn.size(); ++i) {
                    const double xi[2] = {m_glIn[i], m_glIn[j]};
                    append_cell_point(host, xi, m_gwIn[i]*m_gwIn[j]);
                }
            continue;
        }

        // Four congruent sub-triangles times a symmetric Dunavant rule retain every vertex permutation.
        const double V[3][2] = {
            {D2_TRIREF_V[0][0], D2_TRIREF_V[0][1]},
            {D2_TRIREF_V[1][0], D2_TRIREF_V[1][1]},
            {D2_TRIREF_V[2][0], D2_TRIREF_V[2][1]}
        };
        const double M01[2] = {0.5*(V[0][0]+V[1][0]), 0.5*(V[0][1]+V[1][1])};
        const double M12[2] = {0.5*(V[1][0]+V[2][0]), 0.5*(V[1][1]+V[2][1])};
        const double M20[2] = {0.5*(V[2][0]+V[0][0]), 0.5*(V[2][1]+V[0][1])};
        const double sub[4][3][2] = {
            {{V[0][0],V[0][1]}, {M01[0],M01[1]}, {M20[0],M20[1]}},
            {{M01[0],M01[1]}, {V[1][0],V[1][1]}, {M12[0],M12[1]}},
            {{M20[0],M20[1]}, {M12[0],M12[1]}, {V[2][0],V[2][1]}},
            {{M01[0],M01[1]}, {M12[0],M12[1]}, {M20[0],M20[1]}}
        };
        for (const auto& T : sub) {
            const double e1[2] = {T[1][0]-T[0][0], T[1][1]-T[0][1]};
            const double e2[2] = {T[2][0]-T[0][0], T[2][1]-T[0][1]};
            const double scale = std::fabs(e1[0]*e2[1]-e1[1]*e2[0]);
            for (size_t q = 0; q < m_d2FarTriW.size(); ++q) {
                const double l1 = m_d2FarTriP[2*q], l2 = m_d2FarTriP[2*q+1];
                const double xi[2] = {T[0][0]+l1*e1[0]+l2*e2[0],
                                      T[0][1]+l1*e1[1]+l2*e2[1]};
                append_cell_point(host, xi, scale*m_d2FarTriW[q]);
            }
        }
    }

    for (int host = 0; host < m_d2_n_be; ++host) {
        const double* map = &m_d2EdgeMap[static_cast<size_t>(host)*m_d2EdgeMapStride];
        for (size_t q = 0; q < m_d2GlE.size(); ++q) {
            const double t = m_d2GlE[q];
            double density = 0.0;
            for (int a : m_faceCharges[static_cast<size_t>(host)]) {
                const int exponent = m_expo[static_cast<size_t>(3*a)];
                density += charge[static_cast<size_t>(a)]*D2Pow(t, exponent);
            }
            double X[2];
            D2EdgeMap(map, t, X);
            positions.push_back(X[0]);
            positions.push_back(X[1]);
            strengths.push_back(m_d2GwE[q]*density);
        }
    }
    if (strengths.empty()) {
        positions = {0.0, 0.0};
        strengths = {0.0};
    }
    if (!m_image_rot_angle.empty())
        throw std::invalid_argument(
            "2D planar field evaluation does not implement cyclic image rotations yet "
            "(the 2D charge Gram does); drop the rotations or evaluate the field in 3D");
    return std::make_shared<rad_planar_charges::PlanarFieldEvaluator>(
        std::move(positions), std::move(strengths), m_image_masks, m_image_signs);
}

std::vector<std::pair<std::string, double>> RadHACApKChargeGram::LastSolveTimings() const
{
    const SolveTiming& t = m_lastSolveTiming;
    return {
        {"solve_total_s", t.total_s},
        {"solve_factor_s", t.factor_s},
        {"solve_prec_s", t.prec_s},
        {"solve_bx_s", t.bx_s},
        {"solve_gmatvec_s", t.gmatvec_s},
        {"solve_btx_s", t.btx_s},
        {"solve_mass_s", t.mass_s},
        {"solve_dot_s", t.dot_s},
        {"solve_ax_total_s", t.ax_total_s},
        {"solve_ax_other_s", t.ax_other_s},
        {"solve_pcg_update_s", t.pcg_update_s},
        {"solve_apply_count", (double)t.apply_count},
        {"solve_prec_count", (double)t.prec_count},
        {"solve_dot_count", (double)t.dot_count},
        {"hmatvec_total_s", t.hmatvec_total_s},
        {"hmatvec_zero_s", t.hmatvec_zero_s},
        {"hmatvec_permute_s", t.hmatvec_permute_s},
        {"hmatvec_leaf_s", t.hmatvec_leaf_s},
        {"hmatvec_reduce_s", t.hmatvec_reduce_s},
        {"hmatvec_meta_s", t.hmatvec_meta_s},
        {"hmatvec_lowrank_flop_est", t.hmatvec_lowrank_flop_est},
        {"hmatvec_dense_flop_est", t.hmatvec_dense_flop_est},
        {"hmatvec_calls", t.hmatvec_calls},
        {"hmatvec_lowrank_leaves", t.hmatvec_lowrank_leaves},
        {"hmatvec_dense_leaves", t.hmatvec_dense_leaves},
        {"hmatvec_mirrored_upper_leaves", t.hmatvec_mirrored_upper_leaves},
        {"hmatvec_diagonal_leaves", t.hmatvec_diagonal_leaves},
        {"hmatvec_skipped_lower_leaves", t.hmatvec_skipped_lower_leaves},
        {"hmatvec_last_nd", t.hmatvec_last_nd},
        {"hmatvec_last_nthr", t.hmatvec_last_nthr},
    };
}

RadHACApKChargeGram::PicardResult RadHACApKChargeGram::SolveNonlinearPicard(
    const std::vector<int>& B_indptr, const std::vector<int>& B_indices,
    const std::vector<double>& B_data, int n_face,
    const std::vector<int>& mI, const std::vector<int>& mJ, const std::vector<double>& mV,
    const std::vector<double>& Mmass_diag, const std::vector<double>& N_diag,
    const std::vector<double>& mu, double denom,
    double chi0, double Msat, double H0,
    int picard_iters, double cg_tol, int cg_maxit)
{
    const int n_charge = (int)B_indptr.size() - 1;
    // TaskManager self-wrap (AGENTS.md "Parallelization: NGSolve TaskManager"): one region around the
    // whole Picard loop (inner CG matvecs) -> parallel without a caller `with TaskManager()`.
    ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
    auto mmass_apply = [&](const std::vector<double>& x, std::vector<double>& y) {  // y = M_mass x
        y.assign((size_t)n_face, 0.0);
        ngcore::ParallelFor(ngcore::IntRange((int)mV.size()), [&](size_t k) {
            ngcore::AtomicAdd(y[mI[k]], mV[k] * x[mJ[k]]);
        });
    };
    auto N_apply = [&](const std::vector<double>& x, std::vector<double>& y) {        // y = B^T G (B x)
        std::vector<double> q((size_t)n_charge, 0.0), Gq((size_t)n_charge, 0.0);
        ngcore::ParallelFor(ngcore::IntRange(n_charge), [&](size_t a) {
            double s = 0.0;
            for (int k = B_indptr[a]; k < B_indptr[a + 1]; ++k) s += B_data[k] * x[B_indices[k]];
            q[a] = s;
        });
        MatVec(q, Gq);
        y.assign((size_t)n_face, 0.0);
        ngcore::ParallelFor(ngcore::IntRange(n_charge), [&](size_t a) {
            double ga = Gq[a];
            for (int k = B_indptr[a]; k < B_indptr[a + 1]; ++k) ngcore::AtomicAdd(y[B_indices[k]], B_data[k] * ga);
        });
    };
    auto dot = [&](const std::vector<double>& a, const std::vector<double>& b) {
        double s = 0.0;
        ngcore::ParallelForRange(ngcore::IntRange(n_face), [&](ngcore::IntRange r) {
            double local = 0.0;
            for (auto f : r) local += a[f] * b[f];
            ngcore::AtomicAdd(s, local);
        });
        return s;
    };
    // b0 = M_mass mu ; Dscal = mu.(N mu)/denom (the uniform-mode demag factor, Rayleigh quotient).
    std::vector<double> b0, Nmu, mass_m, rhs((size_t)n_face), prec((size_t)n_face);
    mmass_apply(mu, b0);
    N_apply(mu, Nmu);
    double Dscal = dot(mu, Nmu);
    Dscal /= denom;

    std::vector<double> m((size_t)n_face, 0.0);
    double chi = chi0, Mavg = 0.0, Mprev = 0.0;
    int it = 0, done = 0;
    for (; it < picard_iters; ++it) {
        const double inv_chi = 1.0 / chi;
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) {
            prec[f] = inv_chi * Mmass_diag[f] + N_diag[f];
            rhs[f]  = H0 * b0[f];
        });
        int cg_iters = 0;
        m = SolveLinearMaterial(B_indptr, B_indices, B_data, n_face, mI, mJ, mV,
                                inv_chi, prec, rhs, cg_tol, cg_maxit, cg_iters,
                                /*mass_riesz=*/true, /*symmetric=*/true);
        mmass_apply(m, mass_m);
        Mavg = dot(mu, mass_m);
        Mavg /= denom;
        const double Hi = H0 - Dscal * Mavg;
        const double chi_sec = chi0 / (1.0 + chi0 * std::fabs(Hi) / Msat);   // M(H)=chi0 H/(1+chi0|H|/Msat)
        chi = 0.5 * chi + 0.5 * chi_sec;
        done = it + 1;
        if (it > 0 && std::fabs(Mavg - Mprev) < 1e-10 * (std::fabs(Mavg) + 1e-30)) break;
        Mprev = Mavg;
    }
    PicardResult r;
    r.m = m; r.Mavg = Mavg; r.chi = chi; r.Dscal = Dscal; r.iters = done;
    return r;
}
