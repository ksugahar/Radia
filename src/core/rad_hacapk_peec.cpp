/*-------------------------------------------------------------------------
*
* File name:      rad_hacapk_peec.cpp
*
* Project:        RADIA
*
* Description:    Implementation of RadHACApKPEECManager.
*
* First release:  2026-04-16
*
-------------------------------------------------------------------------*/

#include "rad_hacapk_peec.h"

RadHACApKPEECManager::RadHACApKPEECManager(const radia::PEECMatrixBuilder& builder)
    : RadHACApKBase()
    , m_builder(builder)
{
}

void RadHACApKPEECManager::ExtractCoordinates() {
    // One DOF per filament; coordinate = filament center point.
    const auto& segments = m_builder.GetSegments();
    m_n_elem = static_cast<int>(segments.size());
    m_ndof = m_n_elem;
    m_coordinates.resize(static_cast<size_t>(m_n_elem) * 3);
    m_dof_offset.resize(static_cast<size_t>(m_n_elem) + 1);

    for (int i = 0; i < m_n_elem; ++i) {
        m_coordinates[i * 3 + 0] = segments[i].center.x;
        m_coordinates[i * 3 + 1] = segments[i].center.y;
        m_coordinates[i * 3 + 2] = segments[i].center.z;
        m_dof_offset[i] = i;
    }
    m_dof_offset[m_n_elem] = m_ndof;
}

void RadHACApKPEECManager::OnBeforeBuild() {
    // No kernel-specific precomputation needed: Ruehli entries are
    // computed on demand from the filament geometry held by the builder.
    // The HACApK ACA+ loop will call GetInteractionMatrixElement /
    // ComputeSystemEntry directly.
}

void RadHACApKPEECManager::InitializeInvChi() {
    // PEEC has no "1/chi" material diagonal. The system matrix stored
    // by HACApK is L directly; frequency-dependent diagonal (R, jωL,
    // surface impedance Zs) is applied outside HACApK by the enclosing
    // BiCGSTAB solver. Leave m_inv_chi zero to keep UpdateDiagonal
    // a no-op for PEEC (subsequent R/Zs updates use a different path).
    m_inv_chi.assign(static_cast<size_t>(m_ndof), 0.0);
}

double RadHACApKPEECManager::GetInteractionMatrixElement(int dof_i, int dof_j) const {
    return m_builder.InductanceEntry(dof_i, dof_j);
}

//=========================================================================
// Sanity check: N parallel filaments, compare HACApK MatVec against
// dense L matrix multiplication. Returns max relative error.
//=========================================================================

#include <cmath>
#include <cstdio>

double RadHACApKPEECSanityCheck(int n_filaments) {
    if (n_filaments < 2) return -1.0;

    // Build a helical solenoid discretization: representative of a
    // realistic IH work coil. This gives ACA a 3D-spread geometry
    // (good compressibility) rather than a degenerate 1D arrangement.
    // Parameters: radius 50 mm, pitch 6 mm per turn, 25 segments per turn.
    radia::PEECMatrixBuilder builder;
    const double radius = 0.05;
    const double pitch = 0.006;
    const int seg_per_turn = 25;
    const double wire_w = 3.0e-3;
    const double wire_h = 3.0e-3;
    const double two_pi = 6.283185307179586;

    for (int k = 0; k < n_filaments; ++k) {
        double theta0 = two_pi * (double)k / seg_per_turn;
        double theta1 = two_pi * (double)(k + 1) / seg_per_turn;
        double z0 = pitch * (double)k / seg_per_turn;
        double z1 = pitch * (double)(k + 1) / seg_per_turn;
        TVector3d p0(radius * std::cos(theta0), radius * std::sin(theta0), z0);
        TVector3d p1(radius * std::cos(theta1), radius * std::sin(theta1), z1);

        radia::PEECSegment seg;
        seg.center = TVector3d(0.5 * (p0.x + p1.x),
                               0.5 * (p0.y + p1.y),
                               0.5 * (p0.z + p1.z));
        TVector3d d(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
        double dlen = std::sqrt(d.x*d.x + d.y*d.y + d.z*d.z);
        seg.length = dlen;
        seg.direction = TVector3d(d.x / dlen, d.y / dlen, d.z / dlen);
        seg.width = wire_w;
        seg.height = wire_h;
        seg.sigma = 5.8e7;
        builder.AddSegment(seg);
    }

    // Dense reference L via ComputeL.
    radia::PEECMatrices dense;
    dense.n_loop = n_filaments;
    dense.L.assign(static_cast<size_t>(n_filaments) * n_filaments, 0.0);
    builder.ComputeL(dense);

    // HACApK L via the PEEC adapter.
    RadHACApKPEECManager mgr(builder);
    // Use the PEEC-specific default parameters.
    RadHACApKParams params = RadHACApKPEECDefaultParams();
    params.print_level = 0;
    bool ok = mgr.BuildHMatrix(params);
    if (!ok) {
        std::fprintf(stderr, "[PEEC sanity] BuildHMatrix failed\n");
        return -2.0;
    }
    if (mgr.GetNDOF() != n_filaments) {
        std::fprintf(stderr, "[PEEC sanity] ndof mismatch: %d vs %d\n",
                     mgr.GetNDOF(), n_filaments);
        return -3.0;
    }

    // Deterministic test vector (keeps reproducibility across runs).
    std::vector<double> x(n_filaments), y_hmat(n_filaments), y_dense(n_filaments);
    for (int i = 0; i < n_filaments; ++i) {
        x[i] = std::sin(0.1 * i + 0.7) + 0.3 * std::cos(0.27 * i);
    }

    // y_dense = L * x
    for (int i = 0; i < n_filaments; ++i) {
        double s = 0.0;
        for (int j = 0; j < n_filaments; ++j) {
            s += dense.L[(size_t)i * n_filaments + j] * x[j];
        }
        y_dense[i] = s;
    }

    // y_hmat = HACApK_L * x
    mgr.MatVec(x, y_hmat);

    // Max relative error, protected against near-zero denominators.
    double max_abs_dense = 0.0;
    for (double v : y_dense) max_abs_dense = std::max(max_abs_dense, std::fabs(v));
    if (max_abs_dense == 0.0) return 0.0;

    double max_err = 0.0;
    for (int i = 0; i < n_filaments; ++i) {
        double err = std::fabs(y_hmat[i] - y_dense[i]) / max_abs_dense;
        if (err > max_err) max_err = err;
    }
    return max_err;
}
