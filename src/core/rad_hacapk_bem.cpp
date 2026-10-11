/*-------------------------------------------------------------------------
*
* File name:      rad_hacapk_bem.cpp
*
* Project:        RADIA
*
* Description:    Implementation of RadHACApKBEMManager.
*
* First release:  2026-04-30
*
-------------------------------------------------------------------------*/

#include "rad_hacapk_bem.h"

RadHACApKBEMManager::RadHACApKBEMManager(const double* coordinates,
                                          const double* dense_entries,
                                          int n_v)
    : RadHACApKBase()
    , m_coords_ext(coordinates)
    , m_entries_ext(dense_entries)
    , m_n_v(n_v)
{
}

RadHACApKBEMManager::RadHACApKBEMManager(std::shared_ptr<radia::bem::P1EntryProvider> provider, bool double_layer)
    : m_coords_ext(provider->Coordinates()), m_entries_ext(nullptr), m_n_v(provider->Size()),
      m_provider(std::move(provider)), m_double_layer(double_layer) {}

void RadHACApKBEMManager::ExtractCoordinates() {
    m_n_elem = m_n_v;
    m_ndof = m_n_v;
    m_coordinates.resize(static_cast<size_t>(m_n_v) * 3);
    m_dof_offset.resize(static_cast<size_t>(m_n_v) + 1);

    for (int i = 0; i < m_n_v; ++i) {
        m_coordinates[i * 3 + 0] = m_coords_ext[i * 3 + 0];
        m_coordinates[i * 3 + 1] = m_coords_ext[i * 3 + 1];
        m_coordinates[i * 3 + 2] = m_coords_ext[i * 3 + 2];
        m_dof_offset[i] = i;
    }
    m_dof_offset[m_n_v] = m_ndof;
}

void RadHACApKBEMManager::OnBeforeBuild() {
    // No kernel-specific precomputation: entries are already in the
    // dense table.
}

void RadHACApKBEMManager::InitializeInvChi() {
    // Scalar BEM has no "1/chi" material diagonal; the system matrix
    // stored by HACApK is M directly.  Frequency-dependent terms
    // (e.g. Z_s scaling for SIBC) are applied outside HACApK by the
    // enclosing solver loop.
    m_inv_chi.assign(static_cast<size_t>(m_ndof), 0.0);
}

double RadHACApKBEMManager::GetInteractionMatrixElement(int dof_i, int dof_j) const {
    if (m_provider) {
        auto entry=m_provider->Entry(dof_i,dof_j);
        return m_double_layer ? entry.dl : entry.sl;
    }
    // Legacy dense-entry reference route.
    return m_entries_ext[static_cast<size_t>(dof_i) * m_n_v + dof_j];
}
