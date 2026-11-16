// HACApK H-matrix of the tetrahedral volume PEEC partial inductance (see rad_hacapk_peec_tet.h).
#include "rad_hacapk_peec_tet.h"

#include <stdexcept>
#include <utility>

RadHACApKTetPEEC::RadHACApKTetPEEC(std::unique_ptr<radia::peec_tet::InductanceKernel> kernel)
    : RadHACApKBase(), m_kernel(std::move(kernel)) {
    if (!m_kernel || m_kernel->n_dof() == 0)
        throw std::invalid_argument("RadHACApKTetPEEC: empty inductance kernel");
}

double RadHACApKTetPEEC::GetInteractionMatrixElement(int dof_i, int dof_j) const {
    return m_kernel->Entry(dof_i, dof_j);
}

void RadHACApKTetPEEC::ExtractCoordinates() {
    const int n = static_cast<int>(m_kernel->n_dof());
    m_n_elem = n;
    m_ndof = n;
    m_coordinates.resize(static_cast<std::size_t>(n) * 3);
    m_pointRadius.resize(static_cast<std::size_t>(n));
    m_dof_offset.resize(static_cast<std::size_t>(n) + 1);
    for (int i = 0; i < n; ++i) {
        m_kernel->Support(i, &m_coordinates[3 * static_cast<std::size_t>(i)], m_pointRadius[i]);
        m_dof_offset[i] = i;
    }
    m_dof_offset[n] = n;
}

void RadHACApKTetPEEC::InitializeInvChi() {
    // The stored operator is L itself; resistance and frequency are applied by the circuit solver.
    m_inv_chi.assign(static_cast<std::size_t>(m_ndof), 0.0);
}

void RadHACApKTetPEEC::Apply(const std::vector<double>& x, int nrhs, std::vector<double>& y) {
    if (!IsValid())
        throw std::runtime_error("RadHACApKTetPEEC::Apply: the H-matrix is not built");
    MatVecSymMany(x, nrhs, y);
}
