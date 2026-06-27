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
#include <atomic>
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
                                         std::vector<int> image_masks, std::vector<double> image_signs)
    : m_n_el(n_el), m_analytic(true), m_near_factor(near_factor),
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
    std::vector<int> image_masks, std::vector<double> image_signs)
    : m_n_el(n_el), m_analytic(true), m_near_factor(near_factor), m_polytope(true),
      m_image_masks(std::move(image_masks)), m_image_signs(std::move(image_signs))
{
    const int n_cell = n_el;
    const int n_bf   = (int)face_meas.size();
    m_n = n_cell + n_bf;
    m_cent.assign((size_t)m_n * 3, 0.0);
    m_meas.assign((size_t)m_n, 0.0);
    m_size.assign((size_t)m_n, 0.0);
    m_qp.resize(m_n); m_qw.resize(m_n); m_srcTris.resize(m_n);

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
        }
    }
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
    std::vector<double> ref_tri_pts_in, std::vector<double> ref_tri_w_in)
    : m_n_el(n_el), m_highorder(true), m_ho_far_factor(ho_far_factor),
      m_cellV(std::move(cell_verts)), m_faceV(std::move(face_verts)),
      m_host(std::move(charge_host)), m_kind(std::move(charge_kind)), m_expo(std::move(charge_expo))
{
    const int n_cell = n_el;
    const int n_bf   = (int)(m_faceV.size() / 9);
    m_n = (int)m_host.size();                       // number of polynomial CHARGES (the H-matrix dofs)
    { static std::atomic<long long> s_id{0}; m_build_id = s_id.fetch_add(1) + 1; }   // unique id for the QuadDot memo
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
        double tot = 0.0;                          // FACE: sum of sub-triangle Wilton potentials
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
                v[k] = PhiAtHO(src, p);
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
        s += W[k] * (m_highorder ? PhiAtHO(src, p) : PhiAt(src, p));
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
        s += W[k] * PhiAt(src, p);
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

void RadHACApKChargeGram::ExtractCoordinates()
{
    m_n_elem = m_n;
    m_ndof   = m_n;
    m_coordinates = m_cent;   // [n*3] charge centroids (the cluster-tree points)
}

double RadHACApKChargeGram::GetInteractionMatrixElement(int a, int b) const
{
    if (m_highorder) {
        // polynomial charges, symmetrized; the HACApK ACA compresses the well-separated low-rank blocks.
        // NEAR/FAR adaptive quadrature: a well-separated pair uses the cheap LOW-quad plain double-Gauss
        // (QuadDotFar) -- the kernel is smooth there so the expensive HIGH-quad singularity-subtraction is
        // unnecessary; NEAR/self pairs keep the full QuadDot.  This is NOT a monopole far (zero-mean modes
        // have zero monopole) -- it is just a lower quadrature order where the integrand is smooth.
        // m_ho_far_factor = 1e30 (no LOW rule supplied) => every pair NEAR => original all-high-quad path.
        if (a == b) return QuadDot(a, a);                        // self: always the full high-quad subtraction
        if (m_ho_far_factor < 1e29) {
            const double dx = m_cent[3*a]     - m_cent[3*b];
            const double dy = m_cent[3*a + 1] - m_cent[3*b + 1];
            const double dz = m_cent[3*a + 2] - m_cent[3*b + 2];
            const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (r > m_ho_far_factor * (m_size[a] + m_size[b]))
                return 0.5 * (QuadDotFar(a, b) + QuadDotFar(b, a));   // FAR: cheap low-quad plain double-Gauss
        }
        return 0.5 * (QuadDot(a, b) + QuadDot(b, a));            // NEAR: full high-quad subtraction
    }
    if (m_analytic) {
        // Diagonal = the analytic self (the Wilton/phi_tet potential is exact through the 1/r singularity).
        double base;
        if (a == b) {
            base = QuadDot(a, a);
        } else {
            // NEAR/FAR split (build speedup): the analytic entry 0.5*(outer-quad_a . Phi_b + outer-quad_b .
            // Phi_a) is EXPENSIVE (PhiTet/TriPotential per outer point) and only matters for NEAR pairs
            // (the non-uniform-M / div M != 0 interaction); FAR pairs use the cheap centroid-monopole.
            // near_factor = 1e30 (default) => all pairs NEAR => all-analytic (matches the dense
            // build_demag(analytic_gram=True) golden); near_factor ~ 2 gives the fast split.
            const double dx = m_cent[3*a]     - m_cent[3*b];
            const double dy = m_cent[3*a + 1] - m_cent[3*b + 1];
            const double dz = m_cent[3*a + 2] - m_cent[3*b + 2];
            const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (r <= m_near_factor * (m_size[a] + m_size[b]))
                base = 0.5 * (QuadDot(a, b) + QuadDot(b, a));    // NEAR: exact analytic
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
    double tol, int maxit, int& iters_out)
{
    const int n_charge = (int)B_indptr.size() - 1;     // B is n_charge x n_face (CSR over charges)
    // TaskManager self-wrap (AGENTS.md "Parallelization: NGSolve TaskManager"): keep the pool up across
    // the whole CG loop so the Gram H-matvec is parallel without a caller `with TaskManager()`.
    ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
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
        MatVec(q, Gq);                                 // O(N log N) Gram H-matvec
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
    // Jacobi-preconditioned conjugate gradients (the SPD system M^{-1} = 1/prec).
    std::vector<double> x((size_t)n_face, 0.0), r = rhs, z((size_t)n_face), p((size_t)n_face), Ap;
    ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { z[f] = r[f] / prec[f]; });
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
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { z[f] = r[f] / prec[f]; });
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
    double tol, int maxit, int& iters_out)
{
    const int n_charge = (int)B_indptr.size() - 1;
    // Stand up (or reuse the caller's) TaskManager pool so the HACApK H-matvec runs multi-threaded.
    ngcore::RegionTaskManager rtm(radia::GetMaxThreads());

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
        MatVec(q, Gq);                                              // O(N log N) HACApK H-matvec (parallel)
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
    ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { y[f] = r1[f] / prec[f]; }); // y = M^{-1} b
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
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { y[f] = r2[f] / prec[f]; }); // y = M^{-1} r2
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
    double tol, int maxit, int& iters_out)
{
    const int n_charge = (int)B_indptr.size() - 1;
    if (n_charge != m_ncharge) throw std::runtime_error("ChargeGauss SolveLinearMaterial: B row count mismatch");
    ngcore::RegionTaskManager rtm(radia::GetMaxThreads());
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
    ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { z[f] = r[f] / prec[f]; });
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
        ngcore::ParallelFor(ngcore::IntRange(n_face), [&](size_t f) { z[f] = r[f] / prec[f]; });
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
