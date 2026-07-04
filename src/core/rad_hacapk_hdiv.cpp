/* rad_hacapk_hdiv.cpp -- HACApK H-matrix for the symmetric HDiv-type VIM demag operator.
 * See rad_hacapk_hdiv.h.  The build/matvec/stats lifecycle is inherited from RadHACApKBase;
 * this file supplies only the HDiv kernel hooks (coordinates = face centroids, entry = the
 * charge-cluster Coulomb sum N[i][j] = sum_a sum_b B[a][i] G[a][b] B[b][j]). */
#include "rad_hacapk_hdiv.h"
#include "rad_parallel.h"
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

#ifdef HAVE_LAPACK
#include "mkl_pardiso.h"          // PARDISO sparse-direct factor of the RT0 mass for the MASS RIESZ precond
namespace {
// RAII PARDISO SPD (mtype=2 real symmetric positive definite) factor of the RT0 H(div) mass M_mass,
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
        MKL_INT phase = 33, nrhs = 1, idum = 0, error = 0;
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a.data(), ia.data(), ja.data(),
                &idum, &nrhs, iparm, &msglvl, const_cast<double*>(rhs), x, &error);
        if (error != 0)
            throw std::runtime_error("MassRieszPardiso: PARDISO solve phase failed");
    }
};
} // namespace
#endif // HAVE_LAPACK
// Per-thread memo for the high-order QuadDot lives inside QuadDot (a function-local static thread_local map):
// PhiAtHO(src, m_qp[tgt][k]) (the expensive analytic base + inner subtraction loop) depends ONLY on
// (kind,host of tgt, src) -- IDENTICAL across the co-located monomials that share a host's outer points -- so
// the H-matrix fill otherwise recomputes it n_mono(host) times per source.  See QuadDot for the rationale.

RadHACApKHDivManager::RadHACApKHDivManager(int nx, int ny, int nz,
                                           double h, double distort, int nsub)
    : m_nx(nx), m_ny(ny), m_nz(nz), m_nsub(nsub), m_h(h), m_distort(distort)
{
}

void RadHACApKHDivManager::ExtractCoordinates()
{
    m_mesh = rad_hdiv::BuildStructuredRT0(m_nx, m_ny, m_nz, m_h, m_distort);
    const int nf = m_mesh.n_face();
    m_n_elem = nf;   // one "element" per face (the clustering granularity)
    m_ndof   = nf;   // one normal-flux DOF per face
    m_coordinates.assign((size_t)nf * 3, 0.0);
    for (int f = 0; f < nf; ++f) {
        const rad_hdiv::Vec3& c = m_mesh.faces[f].c;
        m_coordinates[(size_t)f * 3 + 0] = c[0];
        m_coordinates[(size_t)f * 3 + 1] = c[1];
        m_coordinates[(size_t)f * 3 + 2] = c[2];
    }
}

void RadHACApKHDivManager::OnBeforeBuild()
{
    rad_hdiv::BuildChargeMapCSC(m_mesh, m_csc);
    rad_hdiv::BuildChargeQuad(m_mesh, m_nsub, m_quad);
    rad_hdiv::BuildMassCOO(m_mesh, m_mI, m_mJ, m_mV, m_mass_diag);
    // System-A H-LU mode: build the O(1) (i,j)->M_mass[i][j] lookup from the COO (RT0 mass couples
    // a face with the other faces of its 1-2 incident cells -> sparse, ~12 nnz/row).
    if (m_system_mode) {
        m_mass_map.clear();
        m_mass_map.reserve(m_mV.size() * 2);
        for (size_t k = 0; k < m_mV.size(); ++k) {
            long long key = (long long)m_mI[k] * (long long)m_ndof + (long long)m_mJ[k];
            m_mass_map[key] += m_mV[k];
        }
    }
}

double RadHACApKHDivManager::ComputeSystemEntry(int dof_i, int dof_j) const
{
    // Default (+N) for the matvec path; system-A (M_mass + chi*N) when SetSystemMode(chi>0) was called.
    double N = GetInteractionMatrixElement(dof_i, dof_j);
    if (!m_system_mode) return N;
    double mass = 0.0;
    auto it = m_mass_map.find((long long)dof_i * (long long)m_ndof + (long long)dof_j);
    if (it != m_mass_map.end()) mass = it->second;
    return mass + m_system_chi * N;
}

void RadHACApKHDivManager::InitializeInvChi()
{
    // The +N H-matrix does not fold in 1/chi (ComputeSystemEntry = default = +N); the material
    // system A = (1/chi) M_mass - N is applied by the caller.  Still must size m_inv_chi.
    m_inv_chi.assign(m_ndof, 0.0);
}

void RadHACApKHDivManager::ApplySystem(const std::vector<double>& x, double inv_chi,
                                       std::vector<double>& y)
{
    // y = N x  (O(N log N) H-matvec via the base)
    y.assign(m_ndof, 0.0);
    MatVec(x, y);
    // y = inv_chi * (M_mass x) - N x
    for (int f = 0; f < m_ndof; ++f) y[f] = -y[f];
    for (size_t k = 0; k < m_mV.size(); ++k)
        y[m_mI[k]] += inv_chi * m_mV[k] * x[m_mJ[k]];
}

std::vector<double> RadHACApKHDivManager::DiagSystem(double inv_chi) const
{
    std::vector<double> d((size_t)m_ndof, 0.0);
    for (int f = 0; f < m_ndof; ++f)
        d[f] = inv_chi * m_mass_diag[f] - GetInteractionMatrixElement(f, f);
    return d;
}

double RadHACApKHDivManager::GetInteractionMatrixElement(int dof_i, int dof_j) const
{
    // N[i][j] = sum_{a in supp(i)} sum_{b in supp(j)} B[a][i] G[a][b] B[b][j].  Each face's
    // support has <= 2 charges -> <= 4 Gram evaluations per matrix entry.
    double acc = 0.0;
    const std::array<int, 2>&    ri = m_csc.rows[dof_i];
    const std::array<double, 2>& ci = m_csc.coef[dof_i];
    const std::array<int, 2>&    rj = m_csc.rows[dof_j];
    const std::array<double, 2>& cj = m_csc.coef[dof_j];
    for (int p = 0; p < 2; ++p) {
        int a = ri[p];
        if (a < 0) continue;
        double ca = ci[p];
        for (int q = 0; q < 2; ++q) {
            int b = rj[q];
            if (b < 0) continue;
            acc += ca * cj[q] * rad_hdiv::CoulombGramEntry(m_quad, a, b);
        }
    }
    return acc;
}

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
// charge-Gram matches the dense Python reference).  Shared by the tet analytic ctor (outer quad on the
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

RadHACApKChargeGram::RadHACApKChargeGram(std::vector<double> centroids,
                                         std::vector<double> measures,
                                         std::vector<double> self_energy)
    : m_cent(std::move(centroids)), m_meas(std::move(measures)), m_self(std::move(self_energy))
{
    m_n = (int)m_meas.size();
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
// sum-of-sub-triangle Wilton potential (face), evaluated from m_srcTris.  Matches the dense Python
// radia.vim._core.analytic_charge_gram polytope path entry-by-entry (same tris, same quad rules).
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
// CurvedTriPotential).  The cell volume charge is DOMINANT (curved RT0 cannot represent uniform M exactly,
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

// ---- CURVED HIGH-ORDER (isoparametric P2) constructor: monomial-charge Gram on a mesh.Curve(2) geometry. ----
// Mirrors the flat HO build but uses the curved P2 map + curved measure for the OUTER quad (xi^expo folded at
// the REFERENCE point, no affine inverse) and the curved Duffy for the INNER potential (PhiInner -> PhiAtHO_
// Curved).  No analytic moments / inner-subtraction table / near-far split (m_ho_far_factor stays 1e30).
RadHACApKChargeGram::RadHACApKChargeGram(
    std::vector<double> cell_nodes, std::vector<double> face_nodes, int n_el, int curve_order,
    std::vector<int> charge_host, std::vector<int> charge_kind, std::vector<int> charge_expo,
    std::vector<double> ref_tet_pts, std::vector<double> ref_tet_w,
    std::vector<double> ref_tri_pts, std::vector<double> ref_tri_w,
    std::vector<double> curve_gl, std::vector<double> curve_gw)
    : m_n_el(n_el), m_curved(true), m_curve_order(curve_order),
      m_cellNodes(std::move(cell_nodes)), m_faceNodes(std::move(face_nodes)),
      m_gl(std::move(curve_gl)), m_gw(std::move(curve_gw)),
      m_highorder(true),
      m_host(std::move(charge_host)), m_kind(std::move(charge_kind)), m_expo(std::move(charge_expo))
{
    const int n_cell = n_el;
    const int n_bf   = (int)(m_faceNodes.size() / 18);
    m_n = (int)m_host.size();
    m_build_id = NextChargeGramBuildId();           // GLOBAL unique id (shared with the high-order ctor)
    m_nmono.assign(m_n, 1);
    {
        std::unordered_map<long long, int> cnt;
        for (int a = 0; a < m_n; ++a) cnt[(long long)m_host[a]*2 + m_kind[a]]++;
        for (int a = 0; a < m_n; ++a) m_nmono[a] = cnt[(long long)m_host[a]*2 + m_kind[a]];
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
            cellM[c][q]  = ref_tet_w[q] * dV;
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
            faceM[f][q]  = ref_tri_w[q] * dA;
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
    const int* e = &m_expo[(size_t)3*src];
    const int deg = e[0] + e[1] + e[2];
    // (1) the host barycentric gradients beta_i (l_i = beta_i . (y - V0)) and V0
    double beta[3][3] = {{0,0,0},{0,0,0},{0,0,0}}, V0[3];
    int ncoord;
    if (m_kind[src] == 0) {                                  // tet cell: l_i = Inv_i . (y - V0)
        const double* V = &m_cellV[(size_t)host*12];
        const double* Inv = &m_cellInv[(size_t)host*9];
        for (int i = 0; i < 3; ++i) { beta[i][0]=Inv[3*i]; beta[i][1]=Inv[3*i+1]; beta[i][2]=Inv[3*i+2]; }
        V0[0]=V[0]; V0[1]=V[1]; V0[2]=V[2]; ncoord = 3;
    } else {                                                 // tri face: l_i = Gi-combination of a_k . (y - V0)
        const double* V = &m_faceV[(size_t)host*9];
        const double* Gi = &m_faceGinv[(size_t)host*4];
        double a1[3], a2[3];
        for (int k=0;k<3;++k){ a1[k]=V[3+k]-V[k]; a2[k]=V[6+k]-V[k]; }
        for (int k=0;k<3;++k){ beta[0][k]=Gi[0]*a1[k]+Gi[1]*a2[k]; beta[1][k]=Gi[2]*a1[k]+Gi[3]*a2[k]; }
        V0[0]=V[0]; V0[1]=V[1]; V0[2]=V[2]; ncoord = 2;
    }
    // (2) collect the (at most 2 for deg<=2) affine factors l_i = alpha_i + beta_i . y
    double facA[2], facB[2][3]; int nf = 0;
    for (int i = 0; i < ncoord; ++i) {
        for (int c = 0; c < e[i]; ++c) {
            if (nf < 2) {
                facB[nf][0]=beta[i][0]; facB[nf][1]=beta[i][1]; facB[nf][2]=beta[i][2];
                facA[nf] = -(beta[i][0]*V0[0]+beta[i][1]*V0[1]+beta[i][2]*V0[2]);
            }
            ++nf;
        }
    }
    // (3) multiply the affine factors -> physical polynomial A + B.y + y^T C y  (nf <= 2)
    double A = 1.0, B[3] = {0,0,0}, C[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int f = 0; f < nf; ++f) {
        const double al = facA[f]; const double* be = facB[f];
        const double nA = A*al; double nB[3], nC[3][3];
        for (int k=0;k<3;++k) nB[k] = A*be[k] + al*B[k];
        for (int k=0;k<3;++k) for (int l=0;l<3;++l) nC[k][l] = al*C[k][l] + 0.5*(B[k]*be[l] + be[k]*B[l]);
        A = nA;
        for (int k=0;k<3;++k) { B[k]=nB[k]; for (int l=0;l<3;++l) C[k][l]=nC[k][l]; }
    }
    // (4) contract with the exact moment potentials
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
        for (int k=0;k<3;++k) for (int l=0;l<3;++l) res += C[k][l]*M2[k][l];
    }
    return res;
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
        return rad_hdiv::CurvedTetPotential(nd, e[0], e[1], e[2], p, m_gl.data(), m_gw.data(), nq);
    }
    const double (*nd)[3] = (const double(*)[3])&m_faceNodes[(size_t)host*18];
    return rad_hdiv::CurvedTriPotential(nd, e[0], e[1], p, m_gl.data(), m_gw.data(), nq);
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
            return tot;          // constant RT0 charge -> monomial exponent 0
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

double RadHACApKChargeGram::QuadDotRefl(int tgt, int src, int mask) const
{
    // IMA image entry G_img(tgt,src) = (1/4pi) INT_tgt Phi_{R(src)} = (1/4pi) INT_tgt Phi_src(R(x))
    // (reflection isometry |x - R(y)| = |R(x) - y|), so we mirror tgt's outer points on the mask axes and
    // evaluate the UNreflected source potential there.  Always the full analytic PhiAt (the image of a
    // charge straddling the mirror is singular at the plane -> needs the exact through-singularity potential).
    const std::vector<rad_hdiv::Vec3>& P = m_qp[tgt];
    const std::vector<double>&         W = m_qw[tgt];
    double s = 0.0;
    for (size_t k = 0; k < P.size(); ++k) {
        double p[3] = {P[k][0], P[k][1], P[k][2]};
        if (mask & 1) p[0] = -p[0];
        if (mask & 2) p[1] = -p[1];
        if (mask & 4) p[2] = -p[2];
        // m_highorder: the monomial-charge inner potential (host-agnostic potential-at-p, PhiAtHO_*);
        // m_analytic/polytope: the constant-charge PhiAt.  Mirror-image charge -> reflected eval point.
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

// ===================================================================== HEX RT1 mode (2026-07-02)
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

// ============================================ WEDGE (PRISM) RT1 geometry (2026-07-04) ===================
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
    const int* e = &m_expo[(size_t)3*charge];       // Q1: e in {0,1} -> a plain conditional product
    double v = 1.0;
    if (e[0]) v *= xi[0];
    if (e[1]) v *= xi[1];
    if (e[2]) v *= xi[2];                            // face charges carry e[2] = 0
    return v;
}

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
    double near_grade, double far_inner_factor)
    : m_n_el(n_el), m_hexmode(true), m_hex_n_bf(n_bf),
      m_hexNodes(std::move(hex_cell_nodes)), m_quadNodes(std::move(quad_face_nodes)),
      m_symTetP(std::move(sym_tet_pts)), m_symTetW(std::move(sym_tet_w)),
      m_symTriP(std::move(sym_tri_pts)), m_symTriW(std::move(sym_tri_w)),
      m_glOut(std::move(gl_out)), m_gwOut(std::move(gw_out)),
      m_glIn(std::move(gl_in)), m_gwIn(std::move(gw_in)),
      m_farTetP(std::move(far_tet_pts)), m_farTetW(std::move(far_tet_w)),
      m_farTriP(std::move(far_tri_pts)), m_farTriW(std::move(far_tri_w)),
      m_near_grade(near_grade), m_far_inner_factor(far_inner_factor),
      m_host(std::move(charge_host)), m_kind(std::move(charge_kind)), m_expo(std::move(charge_expo))
{
    m_n = (int)m_host.size();
    m_build_id = NextChargeGramBuildId();
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
    BuildHexSiteTables();   // static-site radial tables (non-self near inner) + mapped site positions
}

// WEDGE (PRISM) RT1 ctor -- mirror of the hex ctor with 3-sub-tet prism cells + mixed tri/quad faces (see
// the header doc).  Reuses the hex-mode quadrature-table + block-serving members; fills only the wedge
// geometry.  Initializer list is in member DECLARATION order (m_n_el, the shared quad tables, the wedge
// nodes, then m_host/m_kind/m_expo) to avoid -Wreorder.
RadHACApKChargeGram::RadHACApKChargeGram(
    std::vector<double> wedge_cell_nodes, std::vector<double> face_nodes, std::vector<int> face_type,
    int n_el, int n_bf,
    std::vector<int> charge_host, std::vector<int> charge_kind, std::vector<int> charge_expo,
    std::vector<double> sym_tet_pts, std::vector<double> sym_tet_w,
    std::vector<double> sym_tri_pts, std::vector<double> sym_tri_w,
    std::vector<double> gl_out, std::vector<double> gw_out,
    std::vector<double> gl_in, std::vector<double> gw_in,
    std::vector<double> far_tet_pts, std::vector<double> far_tet_w,
    std::vector<double> far_tri_pts, std::vector<double> far_tri_w,
    double near_grade, double far_inner_factor)
    : m_n_el(n_el),
      m_symTetP(std::move(sym_tet_pts)), m_symTetW(std::move(sym_tet_w)),
      m_symTriP(std::move(sym_tri_pts)), m_symTriW(std::move(sym_tri_w)),
      m_glOut(std::move(gl_out)), m_gwOut(std::move(gw_out)),
      m_glIn(std::move(gl_in)), m_gwIn(std::move(gw_in)),
      m_farTetP(std::move(far_tet_pts)), m_farTetW(std::move(far_tet_w)),
      m_farTriP(std::move(far_tri_pts)), m_farTriW(std::move(far_tri_w)),
      m_near_grade(near_grade), m_far_inner_factor(far_inner_factor),
      m_wedgemode(true), m_wedge_n_bf(n_bf),
      m_wCellNodes(std::move(wedge_cell_nodes)), m_wFaceNodes(std::move(face_nodes)),
      m_wFaceType(std::move(face_type)),
      m_host(std::move(charge_host)), m_kind(std::move(charge_kind)), m_expo(std::move(charge_expo))
{
    m_n = (int)m_host.size();
    m_build_id = NextChargeGramBuildId();
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
                            const double m1 = y[0], m2 = y[1], m4 = y[2];   // Q1 monomials, idx = e0+2e1+4e2
                            R.M.push_back(1.0);   R.M.push_back(m1);    R.M.push_back(m2);    R.M.push_back(m1*m2);
                            R.M.push_back(m4);    R.M.push_back(m1*m4); R.M.push_back(m2*m4); R.M.push_back(m1*m2*m4);
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
                        R.M.push_back(1.0); R.M.push_back(yu); R.M.push_back(yv); R.M.push_back(yu*yv);
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
    add("d2CellNodes", m_d2CellNodes); add("d2EdgeNodes", m_d2EdgeNodes);
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
    const int nm = cell ? 8 : 4;
    const int nS = (int)srcG.size();
    int col[8];
    for (int ls = 0; ls < nS; ++ls) {
        const int* e = &m_expo[(size_t)3*srcG[ls]];
        col[ls] = e[0] + 2*e[1] + (cell ? 4*e[2] : 0);
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

void RadHACApKChargeGram::PhiInnerHexSubVec(int kindS, int hS, int subB, const double p[3],
                                            const std::vector<int>& srcG, double* inn) const
{
    const bool cell = (kindS == 0);
    const size_t sid = cell ? ((size_t)hS*6 + subB) : ((size_t)hS*2 + subB);
    const double* cs = cell ? &m_cellSubC[sid*3] : &m_faceSubC[sid*3];
    const double  sz = cell ? m_cellSubS[sid] : m_faceSubS[sid];
    const double dxc = p[0]-cs[0], dyc = p[1]-cs[1], dzc = p[2]-cs[2];
    const bool far_pt = std::sqrt(dxc*dxc + dyc*dyc + dzc*dzc) > m_far_inner_factor*sz;
    if (!far_pt) {
        PhiInnerHexSiteVec(kindS, hS, subB, p, srcG, inn);
        return;
    }
    const double* nd = cell ? &m_hexNodes[(size_t)hS*81] : &m_quadNodes[(size_t)hS*27];
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
    const double* nd = cell ? &m_hexNodes[(size_t)hS*81] : &m_quadNodes[(size_t)hS*27];
    const int nR = (int)m_glIn.size();
    const double* GL = m_glIn.data();
    const double* GW = m_gwIn.data();
    const int nS = (int)srcG.size();
    double acc[8];
    for (int ls = 0; ls < nS; ++ls) acc[ls] = 0.0;

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

// The whole DIRECTED host-pair block (target host (kindT,hT) outer x source host (kindS,hS) inner) for every
// local charge pair, computed in ONE pass.  All near/far/grading decisions are host+sub geometric (identical
// across the block), so the per-entry value is served from ONE block computation -- the expensive kernel
// work on each (outer pt, inner pt) is shared across all nT*nS monomial combos.  Returns [nT*nS]
// row-major, INV4PI folded.
std::vector<double> RadHACApKChargeGram::QuadBlockHex(int kindT, int hT, int kindS, int hS) const
{
    const std::vector<int>& tgtG = (kindT == 0) ? m_cellCharges[hT] : m_faceCharges[hT];
    const std::vector<int>& srcG = (kindS == 0) ? m_cellCharges[hS] : m_faceCharges[hS];
    const int nT = (int)tgtG.size(), nS = (int)srcG.size();
    std::vector<double> blk((size_t)nT*nS, 0.0);
    if (nT == 0 || nS == 0) return blk;
    const bool cellT = (kindT == 0), cellS = (kindS == 0);
    const int nsubT = cellT ? 6 : 2, nsubS = cellS ? 6 : 2;
    const int rt = tgtG[0], rs = srcG[0];      // representative charges (host-level cent/size)
    const double dxh = m_cent[3*rt]-m_cent[3*rs], dyh = m_cent[3*rt+1]-m_cent[3*rs+1],
                 dzh = m_cent[3*rt+2]-m_cent[3*rs+2];
    const double r_h = std::sqrt(dxh*dxh + dyh*dyh + dzh*dzh);
    const bool near_hosts = (kindT == kindS && hT == hS)
                            || r_h <= m_near_grade*(m_size[rt] + m_size[rs]);
    // SELF host pair: the inner takes the RADIAL decomposition with the EXACT anchor xiT (the outer
    // point's own ref coords -- no Newton).  The OUTER grading below is UNCHANGED -- it is required by the
    // Q1 charge degree regardless of how the inner is computed (exact inner + regular outer -> eig 1.088).
    const bool self_pair = (kindT == kindS && hT == hS);
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
            const double* cB = cellS ? &m_cellSubC[sidB*3] : &m_faceSubC[sidB*3];
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
                if (self_pair) PhiInnerHexRadialVec(kindS, hS, sB, pq, xiT, srcG, inn.data());  // radial, exact anchor
                else           PhiInnerHexSubVec(kindS, hS, sB, pq, srcG, inn.data());          // far cloud / radial

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

// thread_local block cache (build_id-guarded, same discipline as the cloud cache).  Keyed by the directed
// (kindT,hT,kindS,hS); a HACApK dense leaf touches all nT*nS entries of a host pair -> computed once, reused.
struct HexBlockKey {
    int kindT;
    int hT;
    int kindS;
    int hS;
    bool operator==(const HexBlockKey& o) const
    {
        return kindT == o.kindT && hT == o.hT && kindS == o.kindS && hS == o.hS;
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
        mix(k.kindT); mix(k.hT); mix(k.kindS); mix(k.hS);
        return h;
    }
};

static thread_local long long s_hex_block_owner = -1;
static thread_local std::unordered_map<HexBlockKey, std::vector<double>, HexBlockKeyHash> s_hex_block_cache;

const std::vector<double>& RadHACApKChargeGram::GetHexBlock(int kindT, int hT, int kindS, int hS) const
{
    if (s_hex_block_owner != m_build_id) { s_hex_block_cache.clear(); s_hex_block_owner = m_build_id; }
    const HexBlockKey key{kindT, hT, kindS, hS};
    auto it = s_hex_block_cache.find(key);
    if (it == s_hex_block_cache.end()) {
        if (s_hex_block_cache.size() > 200000u) s_hex_block_cache.clear();
        it = s_hex_block_cache.emplace(key,
                   m_d2        ? QuadBlock2D(kindT, hT, kindS, hS)
                 : m_wedgemode ? QuadBlockWedge(kindT, hT, kindS, hS)
                 :               QuadBlockHex(kindT, hT, kindS, hS)).first;
    }
    return it->second;
}

// ===================================================================== 2D PLANAR mode (2026-07-03)
// Motor cross-section layer (memory hdiv-vim-tri-quad-motor): kernel -ln(r)/(2pi), charges = -div M on
// tri/quad cells + M.n on boundary edges, Piola-exact REF measures, regular symmetric outer (the log
// kernel's single-layer potentials are continuous -- numpy-validated that NO graded outer is needed),
// radial-cone inner for near/self, cheap far cloud otherwise.  See the header ctor doc.
static const double D2_TRIREF_V[3][2] = {{1, 0}, {0, 1}, {0, 0}};   // NGSolve trig reference

static inline double D2MonoCell(const int* e, const double xi[2])
{
    double v = 1.0;
    if (e[0]) v *= xi[0];
    if (e[1]) v *= xi[1];
    return v;
}

void RadHACApKChargeGram::Tri6Map(const double* nd12, const double xi[2], double X[2])
{
    const double l0 = xi[0], l1 = xi[1], l2 = 1.0 - xi[0] - xi[1];
    const double s[6] = {l0*(2*l0 - 1), l1*(2*l1 - 1), l2*(2*l2 - 1), 4*l0*l1, 4*l1*l2, 4*l2*l0};
    X[0] = X[1] = 0.0;
    for (int k = 0; k < 6; ++k) { X[0] += s[k]*nd12[2*k]; X[1] += s[k]*nd12[2*k + 1]; }
}

void RadHACApKChargeGram::Quad9Map(const double* nd18, const double xi[2], double X[2])
{
    double vx[3], dx[3], vy[3], dy[3];
    HexLag3(xi[0], vx, dx); HexLag3(xi[1], vy, dy);
    X[0] = X[1] = 0.0;
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i) {
            const double s = vx[i]*vy[j];
            const double* nd = &nd18[2*(i + 3*j)];
            X[0] += s*nd[0]; X[1] += s*nd[1];
        }
}

void RadHACApKChargeGram::Edge3Map(const double* nd6, double t, double X[2])
{
    double v[3], d[3];
    HexLag3(t, v, d);
    X[0] = v[0]*nd6[0] + v[1]*nd6[2] + v[2]*nd6[4];
    X[1] = v[0]*nd6[1] + v[1]*nd6[3] + v[2]*nd6[5];
}

// sub-tri ref vertices of cell host h, sub s (tri: itself; quad: 2 sub-tris of [0,1]^2)
static void D2SubTri(int cell_type, int s, double V[3][2])
{
    if (cell_type == 0) {
        for (int i = 0; i < 3; ++i) { V[i][0] = D2_TRIREF_V[i][0]; V[i][1] = D2_TRIREF_V[i][1]; }
    } else {
        const int* tv = QUADREF_TRIS[s];
        for (int i = 0; i < 3; ++i) { V[i][0] = QUADREF_V[tv[i]][0]; V[i][1] = QUADREF_V[tv[i]][1]; }
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

RadHACApKChargeGram::RadHACApKChargeGram(int /*dim2_tag*/,
    std::vector<double> cell_nodes9, std::vector<int> cell_type, std::vector<double> edge_nodes3,
    int n_el, int n_be,
    std::vector<int> charge_host, std::vector<int> charge_kind, std::vector<int> charge_expo,
    std::vector<double> sym_tri_pts, std::vector<double> sym_tri_w,
    std::vector<double> gl_edge, std::vector<double> gw_edge,
    std::vector<double> gl_in, std::vector<double> gw_in,
    std::vector<double> far_tri_pts, std::vector<double> far_tri_w,
    double near_grade, double far_inner_factor)
    : m_n_el(n_el),
      m_glIn(std::move(gl_in)), m_gwIn(std::move(gw_in)),
      m_near_grade(near_grade), m_far_inner_factor(far_inner_factor),
      m_host(std::move(charge_host)), m_kind(std::move(charge_kind)), m_expo(std::move(charge_expo))
{
    m_d2 = true;
    m_d2_n_be = n_be;
    m_d2CellNodes = std::move(cell_nodes9);
    m_d2CellType  = std::move(cell_type);
    m_d2EdgeNodes = std::move(edge_nodes3);
    m_d2SymTriP = std::move(sym_tri_pts); m_d2SymTriW = std::move(sym_tri_w);
    m_d2GlE = std::move(gl_edge); m_d2GwE = std::move(gw_edge);
    m_d2FarTriP = std::move(far_tri_pts); m_d2FarTriW = std::move(far_tri_w);
    m_n = (int)m_host.size();
    m_build_id = NextChargeGramBuildId();
    // ---- per-sub geometry: centroid/size (near test) + mapped anchor sites ----
    m_d2CellSubC.assign((size_t)n_el*2*2, 0.0); m_d2CellSubS.assign((size_t)n_el*2, 0.0);
    m_d2CellSiteX.assign((size_t)n_el*2*7*2, 0.0);
    for (int c = 0; c < n_el; ++c) {
        const double* nd = &m_d2CellNodes[(size_t)c*18];
        const int ct = m_d2CellType[c];
        const int nsub = (ct == 1) ? 2 : 1;
        for (int s = 0; s < nsub; ++s) {
            double V[3][2];
            D2SubTri(ct, s, V);
            double cen[2] = {0, 0}, P[3][2];
            for (int i = 0; i < 3; ++i) {
                if (ct == 0) Tri6Map(nd, V[i], P[i]); else Quad9Map(nd, V[i], P[i]);
                cen[0] += P[i][0]/3.0; cen[1] += P[i][1]/3.0;
            }
            double* pc = &m_d2CellSubC[((size_t)c*2 + s)*2];
            pc[0] = cen[0]; pc[1] = cen[1];
            double sz = 0.0;
            for (int i = 0; i < 3; ++i) {
                const double dx = P[i][0]-cen[0], dy = P[i][1]-cen[1];
                sz = std::max(sz, std::sqrt(dx*dx + dy*dy));
            }
            m_d2CellSubS[(size_t)c*2 + s] = sz;
            for (int k = 0; k < 7; ++k) {
                double x0[2], X[2];
                D2SiteRef(V, k, x0);
                if (ct == 0) Tri6Map(nd, x0, X); else Quad9Map(nd, x0, X);
                double* out = &m_d2CellSiteX[(((size_t)c*2 + s)*7 + k)*2];
                out[0] = X[0]; out[1] = X[1];
            }
        }
    }
    m_d2EdgeC.assign((size_t)n_be*2, 0.0); m_d2EdgeS.assign((size_t)n_be, 0.0);
    m_d2EdgeSiteX.assign((size_t)n_be*3*2, 0.0);
    for (int f = 0; f < n_be; ++f) {
        const double* nd = &m_d2EdgeNodes[(size_t)f*6];
        double P0[2], P1[2], Pm[2];
        Edge3Map(nd, 0.0, P0); Edge3Map(nd, 1.0, P1); Edge3Map(nd, 0.5, Pm);
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
            const int nsub = (ct == 1) ? 2 : 1;
            double cen[2] = {0, 0}, sz = 0.0;
            for (int s = 0; s < nsub; ++s) {
                cen[0] += m_d2CellSubC[((size_t)h*2 + s)*2] / nsub;
                cen[1] += m_d2CellSubC[((size_t)h*2 + s)*2 + 1] / nsub;
            }
            for (int s = 0; s < nsub; ++s) {
                const double dx = m_d2CellSubC[((size_t)h*2 + s)*2] - cen[0];
                const double dy = m_d2CellSubC[((size_t)h*2 + s)*2 + 1] - cen[1];
                sz = std::max(sz, m_d2CellSubS[(size_t)h*2 + s] + std::sqrt(dx*dx + dy*dy));
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
        const double* nd = &m_d2CellNodes[(size_t)hS*18];
        double V[3][2];
        D2SubTri(ct, subB, V);
        const double* cs = &m_d2CellSubC[((size_t)hS*2 + subB)*2];
        const double sz = m_d2CellSubS[(size_t)hS*2 + subB];
        const double dxc = p[0]-cs[0], dyc = p[1]-cs[1];
        if (std::sqrt(dxc*dxc + dyc*dyc) > m_far_inner_factor*sz) {
            // FAR: smooth -ln(r), the fixed bary tri rule mapped on the fly (2D is cheap)
            const int nq = (int)m_d2FarTriW.size();
            for (int q = 0; q < nq; ++q) {
                const double l1 = m_d2FarTriP[2*q], l2 = m_d2FarTriP[2*q + 1];
                const double xi[2] = {V[0][0] + l1*(V[1][0]-V[0][0]) + l2*(V[2][0]-V[0][0]),
                                      V[0][1] + l1*(V[1][1]-V[0][1]) + l2*(V[2][1]-V[0][1])};
                double X[2];
                if (ct == 0) Tri6Map(nd, xi, X); else Quad9Map(nd, xi, X);
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
            const double* sx = &m_d2CellSiteX[(((size_t)hS*2 + subB)*7)*2];
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
                    if (ct == 0) Tri6Map(nd, xi, X); else Quad9Map(nd, xi, X);
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
    const double* nd = &m_d2EdgeNodes[(size_t)hS*6];
    const double* ec = &m_d2EdgeC[(size_t)hS*2];
    const double es = m_d2EdgeS[hS];
    const double dxc = p[0]-ec[0], dyc = p[1]-ec[1];
    const bool far_pt = !xiT && std::sqrt(dxc*dxc + dyc*dyc) > m_far_inner_factor*es;
    const int nq = (int)m_d2GwE.size();
    auto accum = [&](double t, double w) {
        double X[2];
        Edge3Map(nd, t, X);
        const double dx = p[0]-X[0], dy = p[1]-X[1];
        const double r = std::sqrt(dx*dx + dy*dy);
        if (r < 1e-300) return;
        const double g2 = w*(-std::log(r));
        for (int ls = 0; ls < nS; ++ls) {
            const int* e = &m_expo[(size_t)3*srcG[ls]];
            inn[ls] += g2*(e[0] ? t : 1.0);
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
        Edge3Map(nd, 0.0, P0); Edge3Map(nd, 1.0, P1);
        const double du = P1[0]-P0[0], dv = P1[1]-P0[1];
        const double L2 = du*du + dv*dv;
        ts = (L2 > 1e-300) ? ((p[0]-P0[0])*du + (p[1]-P0[1])*dv)/L2 : 0.5;
        ts = std::min(1.0, std::max(0.0, ts));
        for (int it = 0; it < 3; ++it) {                   // Newton polish on the quadratic map
            double v[3], d[3];
            HexLag3(ts, v, d);
            const double X0 = v[0]*nd[0] + v[1]*nd[2] + v[2]*nd[4];
            const double X1 = v[0]*nd[1] + v[1]*nd[3] + v[2]*nd[5];
            const double T0 = d[0]*nd[0] + d[1]*nd[2] + d[2]*nd[4];
            const double T1 = d[0]*nd[1] + d[1]*nd[3] + d[2]*nd[5];
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
std::vector<double> RadHACApKChargeGram::QuadBlock2D(int kindT, int hT, int kindS, int hS) const
{
    const std::vector<int>& tgtG = (kindT == 0) ? m_cellCharges[hT] : m_faceCharges[hT];
    const std::vector<int>& srcG = (kindS == 0) ? m_cellCharges[hS] : m_faceCharges[hS];
    const int nT = (int)tgtG.size(), nS = (int)srcG.size();
    std::vector<double> blk((size_t)nT*nS, 0.0);
    if (nT == 0 || nS == 0) return blk;
    const bool self_pair = (kindT == kindS && hT == hS);
    const int nsubS = (kindS == 0) ? ((m_d2CellType[hS] == 1) ? 2 : 1) : 1;
    std::vector<double> inn(nS);
    auto accumulate = [&](const double xiA[2], double wg, const double Xp[2]) {
        for (int sB = 0; sB < nsubS; ++sB) {
            for (int ls = 0; ls < nS; ++ls) inn[ls] = 0.0;
            PhiInner2DVec(kindS, hS, sB, Xp, self_pair ? xiA : nullptr, srcG, inn.data());
            for (int lt = 0; lt < nT; ++lt) {
                const int* e = &m_expo[(size_t)3*tgtG[lt]];
                const double ma = (kindT == 0) ? D2MonoCell(e, xiA) : (e[0] ? xiA[0] : 1.0);
                double* row = &blk[(size_t)lt*nS];
                for (int ls = 0; ls < nS; ++ls) row[ls] += wg*ma*inn[ls];
            }
        }
    };
    if (kindT == 0) {
        const int ct = m_d2CellType[hT];
        const double* nd = &m_d2CellNodes[(size_t)hT*18];
        const int nsubT = (ct == 1) ? 2 : 1;
        const int nq = (int)m_d2SymTriW.size();
        for (int sA = 0; sA < nsubT; ++sA) {
            double V[3][2];
            D2SubTri(ct, sA, V);
            const double e1u = V[1][0]-V[0][0], e1v = V[1][1]-V[0][1];
            const double e2u = V[2][0]-V[0][0], e2v = V[2][1]-V[0][1];
            const double sc = std::fabs(e1u*e2v - e1v*e2u);            // 2*A_sub(ref)
            for (int q = 0; q < nq; ++q) {
                const double l1 = m_d2SymTriP[2*q], l2 = m_d2SymTriP[2*q + 1];
                const double xiA[2] = {V[0][0] + l1*e1u + l2*e2u, V[0][1] + l1*e1v + l2*e2v};
                double Xp[2];
                if (ct == 0) Tri6Map(nd, xiA, Xp); else Quad9Map(nd, xiA, Xp);
                accumulate(xiA, m_d2SymTriW[q]*sc, Xp);                // W sums 1/2 -> x sc = ref area
            }
        }
    } else {
        const double* nd = &m_d2EdgeNodes[(size_t)hT*6];
        const int nq = (int)m_d2GwE.size();
        for (int q = 0; q < nq; ++q) {
            const double t = m_d2GlE[q];
            const double xiA[2] = {t, 0.0};
            double Xp[2];
            Edge3Map(nd, t, Xp);
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
    cHACApK_set_sym_fill(1);
    const bool ok = RadHACApKBase::BuildHMatrix(params);
    cHACApK_set_sym_fill(0);
    return ok;
}

// ============================================ WEDGE (PRISM) RT1 compute (2026-07-04) ===================
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
                            const double m1 = y[0], m2 = y[1], m4 = y[2];   // Q1 monomials, idx = e0+2e1+4e2
                            R.M.push_back(1.0);   R.M.push_back(m1);    R.M.push_back(m2);    R.M.push_back(m1*m2);
                            R.M.push_back(m4);    R.M.push_back(m1*m4); R.M.push_back(m2*m4); R.M.push_back(m1*m2*m4);
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
                            R.M.push_back(1.0); R.M.push_back(yu); R.M.push_back(yv); R.M.push_back(yu*yv);
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
    for (int k = 0; k < nsite; ++k) {
        const double dx = p[0]-sx[3*k], dy = p[1]-sx[3*k+1], dz = p[2]-sx[3*k+2];
        const double d = dx*dx + dy*dy + dz*dz;
        if (d < bd) { bd = d; best = k; }
    }
    const HexSiteRad& R = cell ? m_wCellSiteRad[(size_t)subB*15 + best]
                               : (ft == 0 ? m_wFaceSiteRadTri[(size_t)best]
                                          : m_wFaceSiteRadQuad[(size_t)subB*7 + best]);
    const double* nd = cell ? &m_wCellNodes[(size_t)hS*54] : &m_wFaceNodes[(size_t)hS*27];
    const int nn = cell ? 18 : (ft == 0 ? 6 : 9);
    const int nm = cell ? 8 : 4;
    const int nS = (int)srcG.size();
    int col[8];
    for (int ls = 0; ls < nS; ++ls) { const int* e = &m_expo[(size_t)3*srcG[ls]]; col[ls] = e[0] + 2*e[1] + (cell ? 4*e[2] : 0); }
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

// Far field point -> cheap cached far cloud; else -> the static-site radial.  Mirror of PhiInnerHexSubVec.
void RadHACApKChargeGram::PhiInnerWedgeSubVec(int kindS, int hS, int subB, const double p[3],
                                              const std::vector<int>& srcG, double* inn) const
{
    const bool cell = (kindS == 0);
    const int ft = cell ? -1 : m_wFaceType[hS];
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

// SELF inner: the exact-anchor (xiT) REF-frame radial decomposition.  Mirror of PhiInnerHexRadialVec.
void RadHACApKChargeGram::PhiInnerWedgeRadialVec(int kindS, int hS, int subB, const double p[3],
                                                 const double* xiT, const std::vector<int>& srcG,
                                                 double* inn) const
{
    if (!xiT) throw std::logic_error("PhiInnerWedgeRadialVec: xiT required (SELF-only)");
    const bool cell = (kindS == 0);
    const int ft = cell ? -1 : m_wFaceType[hS];
    const double* nd = cell ? &m_wCellNodes[(size_t)hS*54] : &m_wFaceNodes[(size_t)hS*27];
    const int nR = (int)m_glIn.size();
    const double* GL = m_glIn.data();
    const double* GW = m_gwIn.data();
    const int nS = (int)srcG.size();
    double acc[8];
    for (int ls = 0; ls < nS; ++ls) acc[ls] = 0.0;
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

// Directed host-pair block (mirror of QuadBlockHex) with mixed-face sub counts / node strides.
std::vector<double> RadHACApKChargeGram::QuadBlockWedge(int kindT, int hT, int kindS, int hS) const
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
    const int rt = tgtG[0], rs = srcG[0];
    const double dxh = m_cent[3*rt]-m_cent[3*rs], dyh = m_cent[3*rt+1]-m_cent[3*rs+1], dzh = m_cent[3*rt+2]-m_cent[3*rs+2];
    const double r_h = std::sqrt(dxh*dxh + dyh*dyh + dzh*dzh);
    const bool near_hosts = (kindT == kindS && hT == hS) || r_h <= m_near_grade*(m_size[rt] + m_size[rs]);
    const bool self_pair = (kindT == kindS && hT == hS);
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
            const double* cB = cellS ? &m_wCellSubC[sidB*3] : &m_wFaceSubC[sidB*3];
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
                for (int i = 0; i < nvT; ++i) {
                    const double ddx = subVA[3*i]-cB[0], ddy = subVA[3*i+1]-cB[1], ddz = subVA[3*i+2]-cB[2];
                    const double d = ddx*ddx + ddy*ddy + ddz*ddz;
                    if (d < best) { best = d; corner = i; }
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
                if (self_pair) PhiInnerWedgeRadialVec(kindS, hS, sB, pq, xiT, srcG, inn.data());
                else           PhiInnerWedgeSubVec(kindS, hS, sB, pq, srcG, inn.data());
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

double RadHACApKChargeGram::GetInteractionMatrixElement(int a, int b) const
{
    // Fail-loud bounds guard (2026-07-03 flake hunt): a HACApK-side index bug (1-based lod handling /
    // buffer overrun) would otherwise read garbage hosts and produce plausible-but-wrong blocks.
    if (a < 0 || a >= m_n || b < 0 || b >= m_n)
        throw std::out_of_range("ChargeGram entry index out of range: a=" + std::to_string(a)
                                + " b=" + std::to_string(b) + " n=" + std::to_string(m_n));
    if (m_d2) {
        // 2D planar mode: served block-wise like the hex mode, symmetrized 0.5*(AB + BA).  Each scalar
        // is read BEFORE the next GetHexBlock fetch -- the memo's capacity clear would otherwise leave a
        // dangling reference (the same use-after-free family as the cloud-cache n=10 crash).
        const int kA = m_kind[a], hA = m_host[a], kB = m_kind[b], hB = m_host[b];
        const int la = m_hexLocalOf[a], lb = m_hexLocalOf[b];
        const int nB = (kB == 0) ? (int)m_cellCharges[hB].size() : (int)m_faceCharges[hB].size();
        const int nA = (kA == 0) ? (int)m_cellCharges[hA].size() : (int)m_faceCharges[hA].size();
        const double vAB = GetHexBlock(kA, hA, kB, hB)[(size_t)la*nB + lb];
        const double vBA = GetHexBlock(kB, hB, kA, hA)[(size_t)lb*nA + la];
        return 0.5*(vAB + vBA);
    }
    if (m_hexmode || m_wedgemode) {
        // HEX / WEDGE RT1: the pair-graded scheme (near subs -> both-domains-graded Duffy outer; far -> the
        // regular symmetric outer; inner always graded/far-dispatched), symmetrized like the other modes.
        // The wedge mode shares this block-serving path verbatim (GetHexBlock -> QuadBlockWedge dispatch).
        // Served from the whole-host-pair block memo (the 64x co-location win) -- bit-identical to
        // the symmetrized 0.5*(block_AB + block_BA) per-entry value, kernel work shared per block.
        // Each scalar is read BEFORE the next GetHexBlock fetch: the memo's capacity clear (fires on
        // ~20k-charge meshes) would otherwise leave a dangling reference.
        const int kA = m_kind[a], hA = m_host[a], kB = m_kind[b], hB = m_host[b];
        const int la = m_hexLocalOf[a], lb = m_hexLocalOf[b];
        const int nB = (kB == 0) ? (int)m_cellCharges[hB].size() : (int)m_faceCharges[hB].size();
        const int nA = (kA == 0) ? (int)m_cellCharges[hA].size() : (int)m_faceCharges[hA].size();
        const double vAB = GetHexBlock(kA, hA, kB, hB)[(size_t)la*nB + lb];   // target A, source B
        const double vBA = GetHexBlock(kB, hB, kA, hA)[(size_t)lb*nA + la];   // target B, source A
        return 0.5*(vAB + vBA);
    }
    if (m_highorder) {
        // polynomial charges, symmetrized; the HACApK ACA compresses the well-separated low-rank blocks.
        // NEAR/FAR adaptive quadrature: a well-separated pair uses the cheap LOW-quad plain double-Gauss
        // (QuadDotFar) -- the kernel is smooth there so the expensive HIGH-quad singularity-subtraction is
        // unnecessary; NEAR/self pairs keep the full QuadDot.  This is NOT a monopole far (zero-mean modes
        // have zero monopole) -- it is just a lower quadrature order where the integrand is smooth.
        // m_ho_far_factor = 1e30 (no LOW rule supplied) => every pair NEAR => original all-high-quad path.
        double base;
        if (a == b) {
            base = QuadDot(a, a);                                // self: always the full high-quad subtraction
        } else if (m_ho_far_factor < 1e29 &&
                   [&]{ const double dx = m_cent[3*a]-m_cent[3*b], dy = m_cent[3*a+1]-m_cent[3*b+1],
                                     dz = m_cent[3*a+2]-m_cent[3*b+2];
                        return std::sqrt(dx*dx + dy*dy + dz*dz) > m_ho_far_factor * (m_size[a] + m_size[b]); }()) {
            base = 0.5 * (QuadDotFar(a, b) + QuadDotFar(b, a));   // FAR: cheap low-quad plain double-Gauss
        } else {
            base = 0.5 * (QuadDot(a, b) + QuadDot(b, a));         // NEAR: full high-quad subtraction
        }
        // IMA: fold in the mirror-image charge interactions (QuadDotRefl uses PhiInner in this mode) so a
        // reduced (1/2, 1/4, 1/8) symmetry model reproduces the full model -- G_IMA = G + sum_i sign_i *
        // 0.5*(refl(a,b)+refl(b,a)).  Empty image => plain highorder.
        for (size_t i = 0; i < m_image_masks.size(); ++i)
            base += m_image_signs[i] * 0.5 *
                    (QuadDotRefl(a, b, m_image_masks[i]) + QuadDotRefl(b, a, m_image_masks[i]));
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
            // build_demag(analytic_gram=True) golden); near_factor ~ 2 gives the fast split.
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
                    (QuadDotRefl(a, b, m_image_masks[i]) + QuadDotRefl(b, a, m_image_masks[i]));
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
    double tol, int maxit, int& iters_out, bool mass_riesz, bool symmetric)
{
    const int n_charge = (int)B_indptr.size() - 1;     // B is n_charge x n_face (CSR over charges)
    // TaskManager self-wrap (AGENTS.md "Parallelization: NGSolve TaskManager"): keep the pool up across
    // the whole CG loop so the Gram H-matvec is parallel without a caller `with TaskManager()`.
    ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
    // MASS RIESZ preconditioner (the default 'auto' path): z = M_mass^{-1} r via a single PARDISO SPD
    // factor of the RT0 mass (built once, applied per iteration).  ~3-5x fewer iters than the diagonal
    // Jacobi (the diag under-resolves the RT0 mass off-diagonal coupling) and nearly mu_r-flat.  When
    // mass_riesz is false the legacy diagonal Jacobi z = r/prec is used (linear_solver="cpp-cg").
#ifdef HAVE_LAPACK
    std::unique_ptr<MassRieszPardiso> mr;
    if (mass_riesz) {
        mr = std::make_unique<MassRieszPardiso>();
        if (!mr->Factor(mI, mJ, mV, n_face))
            throw std::runtime_error("SolveLinearMaterial: PARDISO SPD factor of the RT0 mass "
                                     "(mass Riesz preconditioner) failed");
    }
#else
    if (mass_riesz)
        throw std::runtime_error("SolveLinearMaterial: mass Riesz preconditioner requires MKL PARDISO "
                                 "(HAVE_LAPACK)");
#endif
    auto applyPrec = [&](const std::vector<double>& rr, std::vector<double>& zz) {
#ifdef HAVE_LAPACK
        if (mass_riesz) { mr->Solve(rr.data(), zz.data()); return; }
#endif
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { zz[f] = rr[f] / prec[f]; });
    };
    // A x = inv_chi*(M_mass x) + B^T (G (B x)), with G applied as the charge-Gram H-matvec.
    std::vector<double> q((size_t)n_charge), Gq((size_t)n_charge);
    auto applyA = [&](const std::vector<double>& x, std::vector<double>& y) {
        std::fill(q.begin(), q.end(), 0.0);
        ngcore::ParallelFor(ngcore::IntRange(n_charge), [&](size_t a) {
            double s = 0.0;
            for (int k = B_indptr[a]; k < B_indptr[a + 1]; ++k) s += B_data[k] * x[B_indices[k]];
            q[a] = s;
        });
        std::fill(Gq.begin(), Gq.end(), 0.0);
        if (symmetric) MatVecSym(q, Gq);               // EXACTLY symmetric -> CG-valid Gram apply
        else           MatVec(q, Gq);                  // shadowed: also MatVecSym (sym-fill leaves lower empty)
        y.assign((size_t)n_face, 0.0);
        ngcore::ParallelFor(ngcore::IntRange(n_charge), [&](size_t a) {
            double ga = Gq[a];
            for (int k = B_indptr[a]; k < B_indptr[a + 1]; ++k) ngcore::AtomicAdd(y[B_indices[k]], B_data[k] * ga);
        });
        ngcore::ParallelFor(ngcore::IntRange((int)mV.size()), [&](size_t k) {
            ngcore::AtomicAdd(y[mI[k]], inv_chi * mV[k] * x[mJ[k]]);
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
    // Preconditioned conjugate gradients (SPD system; M^{-1} = mass Riesz or 1/prec diagonal Jacobi).
    std::vector<double> x((size_t)n_face, 0.0), r = rhs, z((size_t)n_face), p((size_t)n_face), Ap;
    applyPrec(r, z);
    p = z;
    double rz = dot(r, z);
    double bnorm = dot(rhs, rhs);
    bnorm = std::sqrt(bnorm); if (bnorm == 0.0) bnorm = 1.0;
    int it = 0;
    for (; it < maxit; ++it) {
        double rnorm = dot(r, r);
        if (std::sqrt(rnorm) <= tol * bnorm) break;
        applyA(p, Ap);
        double pAp = dot(p, Ap);
        double alpha = rz / pAp;
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { x[f] += alpha * p[f]; r[f] -= alpha * Ap[f]; });
        applyPrec(r, z);
        double rz_new = dot(r, z);
        double beta = rz_new / rz;
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { p[f] = z[f] + beta * p[f]; });
        rz = rz_new;
    }
    iters_out = it;
    return x;
}

std::vector<double> RadHACApKChargeGram::SolveMaterialMINRES(
    const std::vector<int>& B_indptr, const std::vector<int>& B_indices,
    const std::vector<double>& B_data, int n_face,
    const std::vector<int>& mI, const std::vector<int>& mJ, const std::vector<double>& mV,
    double inv_chi, const std::vector<double>& prec, const std::vector<double>& rhs,
    double tol, int maxit, int& iters_out, bool mass_riesz, bool symmetric)
{
    const int n_charge = (int)B_indptr.size() - 1;
    // Stand up (or reuse the caller's) TaskManager pool so the HACApK H-matvec runs multi-threaded.
    ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
    // MASS RIESZ preconditioner: y = M_mass^{-1} r via a single PARDISO SPD factor of the RT0 mass
    // (built once, applied per iteration); the bounded -N spectrum (eigenvalues vs M_mass = inv_chi - d,
    // d in [0,1]) makes the mass Riesz especially effective.  mass_riesz=false -> diagonal Jacobi y=r/prec.
#ifdef HAVE_LAPACK
    std::unique_ptr<MassRieszPardiso> mr;
    if (mass_riesz) {
        mr = std::make_unique<MassRieszPardiso>();
        if (!mr->Factor(mI, mJ, mV, n_face))
            throw std::runtime_error("SolveMaterialMINRES: PARDISO SPD factor of the RT0 mass "
                                     "(mass Riesz preconditioner) failed");
    }
#else
    if (mass_riesz)
        throw std::runtime_error("SolveMaterialMINRES: mass Riesz preconditioner requires MKL PARDISO "
                                 "(HAVE_LAPACK)");
#endif
    auto applyPrec = [&](const std::vector<double>& rr, std::vector<double>& zz) {
#ifdef HAVE_LAPACK
        if (mass_riesz) { mr->Solve(rr.data(), zz.data()); return; }
#endif
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { zz[f] = rr[f] / prec[f]; });
    };

    // A x = inv_chi*(M_mass x) - B^T (G (B x))  -- symmetric INDEFINITE -> MINRES.
    std::vector<double> q((size_t)n_charge), Gq((size_t)n_charge);
    auto applyA = [&](const std::vector<double>& x, std::vector<double>& y) {
        std::fill(q.begin(), q.end(), 0.0);
        ngcore::ParallelFor(ngcore::IntRange(n_charge), [&](size_t a) {
            double s = 0.0;
            for (int k = B_indptr[a]; k < B_indptr[a + 1]; ++k) s += B_data[k] * x[B_indices[k]];
            q[a] = s;
        });
        std::fill(Gq.begin(), Gq.end(), 0.0);
        if (symmetric) MatVecSym(q, Gq);                           // EXACTLY symmetric -> MINRES-valid Gram apply
        else           MatVec(q, Gq);                              // shadowed: also MatVecSym (sym-fill leaves lower empty)
        y.assign((size_t)n_face, 0.0);
        ngcore::ParallelFor(ngcore::IntRange(n_charge), [&](size_t a) { // y = -B^T (G B x)
            double ga = Gq[a];
            for (int k = B_indptr[a]; k < B_indptr[a + 1]; ++k) ngcore::AtomicAdd(y[B_indices[k]], -B_data[k] * ga);
        });
        ngcore::ParallelFor(ngcore::IntRange((int)mV.size()), [&](size_t k) {
            ngcore::AtomicAdd(y[mI[k]], inv_chi * mV[k] * x[mJ[k]]);  // + inv_chi M_mass x
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

    // ---- Jacobi-preconditioned MINRES (Paige-Saunders 1975; scipy.sparse.linalg.minres recurrence) ----
    std::vector<double> x((size_t)n_face, 0.0), r1 = rhs, r2 = rhs, y((size_t)n_face);
    applyPrec(r1, y);                                               // y = M^{-1} b
    double beta1 = dot(r1, y);                                      // b . M^{-1} b
    iters_out = 0;
    if (beta1 <= 0.0) return x;                                     // b = 0 (or M not SPD) -> x = 0
    beta1 = std::sqrt(beta1);
    double oldb = 0.0, beta = beta1, dbar = 0.0, epsln = 0.0, phibar = beta1, cs = -1.0, sn = 0.0;
    std::vector<double> v((size_t)n_face), Av((size_t)n_face),
                        w((size_t)n_face, 0.0), w1((size_t)n_face, 0.0), w2((size_t)n_face, 0.0);
    int it = 0;
    for (; it < maxit; ++it) {
        const double s = 1.0 / beta;
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { v[f] = s * y[f]; }); // Lanczos vector
        applyA(v, Av);
        if (it >= 1) {
            const double c = beta / oldb;
            ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { Av[f] -= c * r1[f]; });
        }
        const double alfa = dot(v, Av);
        {
            const double c = alfa / beta;
            ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { Av[f] -= c * r2[f]; });
        }
        r1 = r2; r2 = Av;
        applyPrec(r2, y);                                          // y = M^{-1} r2
        oldb = beta; beta = dot(r2, y);
        if (beta < 0.0) break;                                      // preconditioner not SPD
        beta = std::sqrt(beta);
        // previous + next Givens rotation
        const double oldeps = epsln;
        const double delta  = cs * dbar + sn * alfa;
        const double gbar   = sn * dbar - cs * alfa;
        epsln = sn * beta;
        dbar  = -cs * beta;
        double gamma = std::sqrt(gbar * gbar + beta * beta);
        if (gamma < 1e-300) gamma = 1e-300;
        cs = gbar / gamma; sn = beta / gamma;
        const double phi = cs * phibar; phibar = sn * phibar;
        const double denom = 1.0 / gamma;
        w1 = w2; w2 = w;
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { w[f] = (v[f] - oldeps * w1[f] - delta * w2[f]) * denom; });
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { x[f] += phi * w[f]; });
        iters_out = it + 1;
        if (phibar <= tol * beta1) break;                          // relative preconditioned residual
    }
    return x;
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
    std::vector<double> b0, Nmu, Mmm, rhs((size_t)n_face), prec((size_t)n_face);
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
                                inv_chi, prec, rhs, cg_tol, cg_maxit, cg_iters);
        mmass_apply(m, Mmm);
        Mavg = dot(mu, Mmm);
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

//=========================================================================
// RadHACApKPointKernel / RadHACApKChargeGaussOperator
//=========================================================================

RadHACApKPointKernel::RadHACApKPointKernel(std::vector<double> points)
    : m_points(std::move(points))
{
}

void RadHACApKPointKernel::ExtractCoordinates()
{
    m_n_elem = (int)(m_points.size() / 3);
    m_ndof = m_n_elem;
    m_coordinates = m_points;
}

double RadHACApKPointKernel::GetInteractionMatrixElement(int i, int j) const
{
    if (i == j) return 0.0;  // self singularity is carried by the charge-level near correction.
    const double dx = m_points[(size_t)3*i]     - m_points[(size_t)3*j];
    const double dy = m_points[(size_t)3*i + 1] - m_points[(size_t)3*j + 1];
    const double dz = m_points[(size_t)3*i + 2] - m_points[(size_t)3*j + 2];
    const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
    return (r > 1e-300) ? (RAD_INV_FOUR_PI / r) : 0.0;
}

RadHACApKChargeGaussOperator::RadHACApKChargeGaussOperator(
    std::vector<double> point_coords,
    std::vector<int> P_pt,
    std::vector<int> P_chg,
    std::vector<double> P_coef,
    int n_charge,
    std::vector<int> corr_i,
    std::vector<int> corr_j,
    std::vector<double> corr_v)
    : m_ncharge(n_charge),
      m_point_coords(std::move(point_coords)),
      m_corr_i(std::move(corr_i)),
      m_corr_j(std::move(corr_j)),
      m_corr_v(std::move(corr_v))
{
    m_npoint = (int)(m_point_coords.size() / 3);
    if ((int)m_point_coords.size() != 3 * m_npoint)
        throw std::runtime_error("RadHACApKChargeGaussOperator: point_coords size not a multiple of 3");
    const size_t nnz = P_coef.size();
    if (P_pt.size() != nnz || P_chg.size() != nnz)
        throw std::runtime_error("RadHACApKChargeGaussOperator: inconsistent P-scatter array sizes");
    if ((int)m_corr_i.size() != (int)m_corr_j.size() || (int)m_corr_i.size() != (int)m_corr_v.size())
        throw std::runtime_error("RadHACApKChargeGaussOperator: inconsistent correction array sizes");
    // Build BOTH CSR orientations of the scatter P from the COO triple (counting-sort by point, by charge).
    m_pt_indptr.assign((size_t)m_npoint + 1, 0);
    m_chg_indptr.assign((size_t)m_ncharge + 1, 0);
    for (size_t k = 0; k < nnz; ++k) {
        const int p = P_pt[k], a = P_chg[k];
        if (p < 0 || p >= m_npoint) throw std::runtime_error("RadHACApKChargeGaussOperator: P_pt out of range");
        if (a < 0 || a >= m_ncharge) throw std::runtime_error("RadHACApKChargeGaussOperator: P_chg out of range");
        ++m_pt_indptr[(size_t)p + 1];
        ++m_chg_indptr[(size_t)a + 1];
    }
    for (int p = 0; p < m_npoint; ++p)  m_pt_indptr[(size_t)p + 1]  += m_pt_indptr[(size_t)p];
    for (int a = 0; a < m_ncharge; ++a) m_chg_indptr[(size_t)a + 1] += m_chg_indptr[(size_t)a];
    m_pt_charge.resize(nnz);  m_pt_coef.resize(nnz);
    m_chg_point.resize(nnz);  m_chg_coef.resize(nnz);
    std::vector<int> pcur(m_pt_indptr.begin(), m_pt_indptr.end() - 1);
    std::vector<int> ccur(m_chg_indptr.begin(), m_chg_indptr.end() - 1);
    for (size_t k = 0; k < nnz; ++k) {
        const int p = P_pt[k], a = P_chg[k];
        const double c = P_coef[k];
        const int ip = pcur[(size_t)p]++;  m_pt_charge[(size_t)ip] = a; m_pt_coef[(size_t)ip] = c;
        const int ic = ccur[(size_t)a]++;  m_chg_point[(size_t)ic] = p; m_chg_coef[(size_t)ic] = c;
    }
    m_corr_map.reserve(m_corr_v.size() * 2 + 1);
    for (size_t k = 0; k < m_corr_v.size(); ++k) {
        const int i = m_corr_i[k], j = m_corr_j[k];
        if (i < 0 || i >= m_ncharge || j < 0 || j >= m_ncharge)
            throw std::runtime_error("RadHACApKChargeGaussOperator: correction index out of range");
        m_corr_map[(long long)i * (long long)m_ncharge + (long long)j] += m_corr_v[k];
    }
}

bool RadHACApKChargeGaussOperator::BuildHMatrix(const RadHACApKParams& params)
{
    m_kernel.reset(new RadHACApKPointKernel(m_point_coords));
    return m_kernel->BuildHMatrix(params);
}

double RadHACApKChargeGaussOperator::PointDirectEntry(int a, int b) const
{
    // (1/4pi) sum_{p in supp(a)} sum_{q in supp(b), q!=p} coef_a(p) coef_b(q) / |x_p - x_q|.
    double s = 0.0;
    for (int ka = m_chg_indptr[(size_t)a]; ka < m_chg_indptr[(size_t)a + 1]; ++ka) {
        const int p = m_chg_point[(size_t)ka];
        const double ca = m_chg_coef[(size_t)ka];
        const double x0 = m_point_coords[(size_t)3*p];
        const double x1 = m_point_coords[(size_t)3*p + 1];
        const double x2 = m_point_coords[(size_t)3*p + 2];
        for (int kb = m_chg_indptr[(size_t)b]; kb < m_chg_indptr[(size_t)b + 1]; ++kb) {
            const int q = m_chg_point[(size_t)kb];
            if (p == q) continue;
            const double dx = x0 - m_point_coords[(size_t)3*q];
            const double dy = x1 - m_point_coords[(size_t)3*q + 1];
            const double dz = x2 - m_point_coords[(size_t)3*q + 2];
            const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (r > 1e-300) s += ca * m_chg_coef[(size_t)kb] * RAD_INV_FOUR_PI / r;
        }
    }
    return s;
}

double RadHACApKChargeGaussOperator::GetChargeEntry(int a, int b) const
{
    double v = PointDirectEntry(a, b);
    auto it = m_corr_map.find((long long)a * (long long)m_ncharge + (long long)b);
    if (it != m_corr_map.end()) v += it->second;
    return v;
}

void RadHACApKChargeGaussOperator::MatVec(const std::vector<double>& q, std::vector<double>& y)
{
    if (!m_kernel || !m_kernel->IsValid())
        throw std::runtime_error("RadHACApKChargeGaussOperator.MatVec: H-matrix is not built");
    if ((int)q.size() != m_ncharge)
        throw std::runtime_error("RadHACApKChargeGaussOperator.MatVec: q size mismatch");
    std::vector<double> point_rhs((size_t)m_npoint, 0.0), point_phi((size_t)m_npoint, 0.0);
    // scatter: point_rhs[p] = sum_{(a,coef) at p} coef * q[a]   (per-point CSR -> lock-free)
    ngcore::ParallelFor(ngcore::IntRange(m_npoint), [&](size_t p) {
        double s = 0.0;
        for (int k = m_pt_indptr[p]; k < m_pt_indptr[p + 1]; ++k) s += m_pt_coef[(size_t)k] * q[(size_t)m_pt_charge[(size_t)k]];
        point_rhs[p] = s;
    });
    m_kernel->MatVec(point_rhs, point_phi);
    y.assign((size_t)m_ncharge, 0.0);
    // gather: y[a] = sum_{(p,coef) of a} coef * point_phi[p]   (per-charge CSR -> lock-free)
    ngcore::ParallelFor(ngcore::IntRange(m_ncharge), [&](size_t a) {
        double s = 0.0;
        for (int k = m_chg_indptr[a]; k < m_chg_indptr[a + 1]; ++k) s += m_chg_coef[(size_t)k] * point_phi[(size_t)m_chg_point[(size_t)k]];
        y[a] = s;
    });
    ngcore::ParallelFor(ngcore::IntRange((int)m_corr_v.size()), [&](size_t k) {
        ngcore::AtomicAdd(y[(size_t)m_corr_i[k]], m_corr_v[k] * q[(size_t)m_corr_j[k]]);
    });
}

std::vector<double> RadHACApKChargeGaussOperator::SolveLinearMaterial(
    const std::vector<int>& B_indptr, const std::vector<int>& B_indices,
    const std::vector<double>& B_data, int n_face,
    const std::vector<int>& mI, const std::vector<int>& mJ, const std::vector<double>& mV,
    double inv_chi, const std::vector<double>& prec, const std::vector<double>& rhs,
    double tol, int maxit, int& iters_out, bool mass_riesz)
{
    const int n_charge = (int)B_indptr.size() - 1;
    if (n_charge != m_ncharge) throw std::runtime_error("ChargeGauss SolveLinearMaterial: B row count mismatch");
    ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
    // MASS RIESZ preconditioner (default 'auto'): z = M_mass^{-1} r via one PARDISO SPD factor of the RT0
    // mass; mass_riesz=false keeps the diagonal Jacobi z = r/prec.  Same path as RadHACApKChargeGram.
#ifdef HAVE_LAPACK
    std::unique_ptr<MassRieszPardiso> mr;
    if (mass_riesz) {
        mr = std::make_unique<MassRieszPardiso>();
        if (!mr->Factor(mI, mJ, mV, n_face))
            throw std::runtime_error("ChargeGauss SolveLinearMaterial: PARDISO SPD factor of the RT0 mass "
                                     "(mass Riesz preconditioner) failed");
    }
#else
    if (mass_riesz)
        throw std::runtime_error("ChargeGauss SolveLinearMaterial: mass Riesz preconditioner requires MKL "
                                 "PARDISO (HAVE_LAPACK)");
#endif
    auto applyPrec = [&](const std::vector<double>& rr, std::vector<double>& zz) {
#ifdef HAVE_LAPACK
        if (mass_riesz) { mr->Solve(rr.data(), zz.data()); return; }
#endif
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { zz[f] = rr[f] / prec[f]; });
    };
    std::vector<double> q((size_t)n_charge), Gq((size_t)n_charge);
    auto applyA = [&](const std::vector<double>& x, std::vector<double>& y) {
        std::fill(q.begin(), q.end(), 0.0);
        ngcore::ParallelFor(ngcore::IntRange(n_charge), [&](size_t a) {
            double s = 0.0;
            for (int k = B_indptr[a]; k < B_indptr[a + 1]; ++k) s += B_data[k] * x[B_indices[k]];
            q[a] = s;
        });
        MatVec(q, Gq);
        y.assign((size_t)n_face, 0.0);
        ngcore::ParallelFor(ngcore::IntRange(n_charge), [&](size_t a) {
            const double ga = Gq[a];
            for (int k = B_indptr[a]; k < B_indptr[a + 1]; ++k)
                ngcore::AtomicAdd(y[B_indices[k]], B_data[k] * ga);
        });
        ngcore::ParallelFor(ngcore::IntRange((int)mV.size()), [&](size_t k) {
            ngcore::AtomicAdd(y[mI[k]], inv_chi * mV[k] * x[mJ[k]]);
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
    std::vector<double> x((size_t)n_face, 0.0), r = rhs, z((size_t)n_face), p((size_t)n_face), Ap;
    applyPrec(r, z);
    p = z;
    double rz = dot(r, z);
    double bnorm = std::sqrt(dot(rhs, rhs)); if (bnorm == 0.0) bnorm = 1.0;
    int it = 0;
    for (; it < maxit; ++it) {
        if (std::sqrt(dot(r, r)) <= tol * bnorm) break;
        applyA(p, Ap);
        const double pAp = dot(p, Ap);
        const double alpha = rz / pAp;
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { x[f] += alpha * p[f]; r[f] -= alpha * Ap[f]; });
        applyPrec(r, z);
        const double rz_new = dot(r, z);
        const double beta = rz_new / rz;
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { p[f] = z[f] + beta * p[f]; });
        rz = rz_new;
    }
    iters_out = it;
    return x;
}

//=========================================================================
// RadHACApKHDivSystemTet -- unstructured face-DOF system A = M_mass + chi*N (Phase 2)
//=========================================================================

RadHACApKHDivSystemTet::RadHACApKHDivSystemTet(
    std::vector<double> face_centroids, double chi,
    std::vector<int> face_charge, std::vector<double> face_coef,
    std::vector<int> mI, std::vector<int> mJ, std::vector<double> mV,
    std::vector<double> cell_verts, std::vector<double> face_verts,
    int n_el, double gram_near_factor)
    : m_chi(chi),
      m_face_cent(std::move(face_centroids)),
      m_face_charge(std::move(face_charge)),
      m_face_coef(std::move(face_coef)),
      // embedded analytic charge Gram (constructed -> geometry ready -> G(a,b) via entry, NOT built)
      m_G(std::move(cell_verts), std::move(face_verts), n_el, gram_near_factor)
{
    m_nface = (int)(m_face_cent.size() / 3);
    // O(1) (i,j)->M_mass[i][j] lookup from the COO (RT0 mass is sparse: ~couples within incident cells)
    m_mass_map.reserve(mV.size() * 2);
    for (size_t k = 0; k < mV.size(); ++k)
        m_mass_map[(long long)mI[k] * (long long)m_nface + (long long)mJ[k]] += mV[k];
}

double RadHACApKHDivSystemTet::GetInteractionMatrixElement(int dof_i, int dof_j) const
{
    // N[i][j] = sum_{a in supp(i)} sum_{b in supp(j)} B[a][i] G[a][b] B[b][j]  (<=2 charges/face)
    double N = 0.0;
    for (int p = 0; p < 2; ++p) {
        int a = m_face_charge[(size_t)dof_i * 2 + p];
        if (a < 0) continue;
        double ca = m_face_coef[(size_t)dof_i * 2 + p];
        for (int q = 0; q < 2; ++q) {
            int b = m_face_charge[(size_t)dof_j * 2 + q];
            if (b < 0) continue;
            N += ca * m_face_coef[(size_t)dof_j * 2 + q] * m_G.GetInteractionMatrixElement(a, b);
        }
    }
    double mass = 0.0;
    auto it = m_mass_map.find((long long)dof_i * (long long)m_nface + (long long)dof_j);
    if (it != m_mass_map.end()) mass = it->second;
    return mass + m_chi * N;   // A = M_mass + chi*N
}
