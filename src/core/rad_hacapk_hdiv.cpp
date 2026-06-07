/* rad_hacapk_hdiv.cpp -- HACApK H-matrix for the symmetric HDiv-type VIM demag operator.
 * See rad_hacapk_hdiv.h.  The build/matvec/stats lifecycle is inherited from RadHACApKBase;
 * this file supplies only the HDiv kernel hooks (coordinates = face centroids, entry = the
 * charge-cluster Coulomb sum N[i][j] = sum_a sum_b B[a][i] G[a][b] B[b][j]). */
#include "rad_hacapk_hdiv.h"

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
