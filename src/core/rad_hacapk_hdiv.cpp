/* rad_hacapk_hdiv.cpp -- HACApK H-matrix for the symmetric HDiv-type VIM demag operator.
 * See rad_hacapk_hdiv.h.  The build/matvec/stats lifecycle is inherited from RadHACApKBase;
 * this file supplies only the HDiv kernel hooks (coordinates = face centroids, entry = the
 * charge-cluster Coulomb sum N[i][j] = sum_a sum_b B[a][i] G[a][b] B[b][j]). */
#include "rad_hacapk_hdiv.h"
#include <cmath>
#include <utility>
#include <algorithm>

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

RadHACApKChargeGram::RadHACApKChargeGram(std::vector<double> centroids,
                                         std::vector<double> measures,
                                         std::vector<double> self_energy)
    : m_cent(std::move(centroids)), m_meas(std::move(measures)), m_self(std::move(self_energy))
{
    m_n = (int)m_meas.size();
}

RadHACApKChargeGram::RadHACApKChargeGram(std::vector<double> cell_verts,
                                         std::vector<double> face_verts,
                                         int n_el)
    : m_n_el(n_el), m_analytic(true),
      m_cellV(std::move(cell_verts)), m_faceV(std::move(face_verts))
{
    const int n_bf = (int)(m_faceV.size() / 9);
    m_n = n_el + n_bf;
    m_cent.assign((size_t)m_n * 3, 0.0);
    m_qp.resize(m_n);
    m_qw.resize(m_n);

    // Outer-quad rule on a CELL: tet barycentric sub-points _bary_tet(3) (10 nodes, equal weights,
    // sum = cell volume).  Matches radia.hdiv_vim._core.analytic_charge_gram nsub_tet=3.
    double lamT[10][4];
    int nT = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3 - i; ++j)
            for (int k = 0; k < 3 - i - j; ++k) {
                int l = 2 - i - j - k;
                lamT[nT][0] = (i + 0.25) / 3.0; lamT[nT][1] = (j + 0.25) / 3.0;
                lamT[nT][2] = (k + 0.25) / 3.0; lamT[nT][3] = (l + 0.25) / 3.0;
                ++nT;
            }
    // Outer-quad rule on a FACE: Dunavant degree-5 symmetric triangle rule (7 nodes; weights sum to 1).
    static const double DUN[7][4] = {
        {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.225},
        {0.0597158717, 0.4701420641, 0.4701420641, 0.1323941527},
        {0.4701420641, 0.0597158717, 0.4701420641, 0.1323941527},
        {0.4701420641, 0.4701420641, 0.0597158717, 0.1323941527},
        {0.7974269853, 0.1012865073, 0.1012865073, 0.1259391805},
        {0.1012865073, 0.7974269853, 0.1012865073, 0.1259391805},
        {0.1012865073, 0.1012865073, 0.7974269853, 0.1259391805},
    };

    for (int a = 0; a < n_el; ++a) {
        const double* V = &m_cellV[(size_t)a * 12];   // 4 x 3
        double cx = 0, cy = 0, cz = 0;
        for (int i = 0; i < 4; ++i) { cx += V[3*i]; cy += V[3*i+1]; cz += V[3*i+2]; }
        m_cent[3*a] = cx / 4; m_cent[3*a+1] = cy / 4; m_cent[3*a+2] = cz / 4;
        double e1[3], e2[3], e3[3];
        for (int k = 0; k < 3; ++k) { e1[k] = V[3+k]-V[k]; e2[k] = V[6+k]-V[k]; e3[k] = V[9+k]-V[k]; }
        double cr[3] = {e2[1]*e3[2]-e2[2]*e3[1], e2[2]*e3[0]-e2[0]*e3[2], e2[0]*e3[1]-e2[1]*e3[0]};
        double vol = std::fabs(e1[0]*cr[0] + e1[1]*cr[1] + e1[2]*cr[2]) / 6.0;
        m_qp[a].resize(nT);
        m_qw[a].assign(nT, vol / nT);
        for (int q = 0; q < nT; ++q) {
            rad_hdiv::Vec3 P = {0, 0, 0};
            for (int i = 0; i < 4; ++i) for (int k = 0; k < 3; ++k) P[k] += lamT[q][i] * V[3*i+k];
            m_qp[a][q] = P;
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

double RadHACApKChargeGram::PhiAt(int src, const double p[3]) const
{
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
    double s = 0.0;
    for (size_t k = 0; k < P.size(); ++k) {
        const double p[3] = {P[k][0], P[k][1], P[k][2]};
        s += W[k] * PhiAt(src, p);
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
    if (m_analytic) {
        // Exact analytic charge Gram: G[a][b] = (1/4pi) INT_a Phi_b (Phi_b = PhiTet/TriPotential of
        // source b).  Symmetrize the quadrature: 0.5*(outer-quad on a . Phi_b + outer-quad on b . Phi_a).
        // Diagonal = the analytic self (the Wilton/phi_tet potential is exact through the 1/r singularity).
        if (a == b) return QuadDot(a, a);
        return 0.5 * (QuadDot(a, b) + QuadDot(b, a));
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
    // A x = inv_chi*(M_mass x) + B^T (G (B x)), with G applied as the charge-Gram H-matvec.
    std::vector<double> q((size_t)n_charge), Gq((size_t)n_charge);
    auto applyA = [&](const std::vector<double>& x, std::vector<double>& y) {
        std::fill(q.begin(), q.end(), 0.0);
        for (int a = 0; a < n_charge; ++a) {
            double s = 0.0;
            for (int k = B_indptr[a]; k < B_indptr[a + 1]; ++k) s += B_data[k] * x[B_indices[k]];
            q[a] = s;
        }
        std::fill(Gq.begin(), Gq.end(), 0.0);
        MatVec(q, Gq);                                 // O(N log N) Gram H-matvec
        y.assign((size_t)n_face, 0.0);
        for (int a = 0; a < n_charge; ++a) {
            double ga = Gq[a];
            for (int k = B_indptr[a]; k < B_indptr[a + 1]; ++k) y[B_indices[k]] += B_data[k] * ga;
        }
        for (size_t k = 0; k < mV.size(); ++k) y[mI[k]] += inv_chi * mV[k] * x[mJ[k]];
    };
    // Jacobi-preconditioned conjugate gradients (the SPD system M^{-1} = 1/prec).
    std::vector<double> x((size_t)n_face, 0.0), r = rhs, z((size_t)n_face), p((size_t)n_face), Ap;
    for (int f = 0; f < n_face; ++f) z[f] = r[f] / prec[f];
    p = z;
    double rz = 0.0; for (int f = 0; f < n_face; ++f) rz += r[f] * z[f];
    double bnorm = 0.0; for (int f = 0; f < n_face; ++f) bnorm += rhs[f] * rhs[f];
    bnorm = std::sqrt(bnorm); if (bnorm == 0.0) bnorm = 1.0;
    int it = 0;
    for (; it < maxit; ++it) {
        double rnorm = 0.0; for (int f = 0; f < n_face; ++f) rnorm += r[f] * r[f];
        if (std::sqrt(rnorm) <= tol * bnorm) break;
        applyA(p, Ap);
        double pAp = 0.0; for (int f = 0; f < n_face; ++f) pAp += p[f] * Ap[f];
        double alpha = rz / pAp;
        for (int f = 0; f < n_face; ++f) { x[f] += alpha * p[f]; r[f] -= alpha * Ap[f]; }
        for (int f = 0; f < n_face; ++f) z[f] = r[f] / prec[f];
        double rz_new = 0.0; for (int f = 0; f < n_face; ++f) rz_new += r[f] * z[f];
        double beta = rz_new / rz;
        for (int f = 0; f < n_face; ++f) p[f] = z[f] + beta * p[f];
        rz = rz_new;
    }
    iters_out = it;
    return x;
}
