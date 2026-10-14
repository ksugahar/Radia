#pragma once

#include <algorithm>
#include <cmath>
#include <vector>
#include <mkl_cblas.h>
#include <mkl_service.h>

namespace radia_detail {

// Contract the unchanged product quadrature in bounded target-point tiles.
// The caller supplies image-transformed target positions and a zeroed block.
inline void HexFarProductBlas(
    int nt, int ns, int qt, int qs,
    const double* tx, const double* tw, const double* tv,
    const double* sx, const double* sw, const double* sv, double* block)
{
    constexpr int tile = 32;
    const int capacity = std::min(tile, qt);
    std::vector<double> kernel(static_cast<size_t>(capacity)*qs);
    std::vector<double> inner(static_cast<size_t>(capacity)*ns);
    std::vector<double> outer(static_cast<size_t>(capacity)*nt);
    // NGSolve owns worker parallelism; restore this thread's MKL setting.
    struct LocalThreads {
        int previous = mkl_set_num_threads_local(1);
        ~LocalThreads() { mkl_set_num_threads_local(previous); }
    } threads;
    for (int begin = 0; begin < qt; begin += tile) {
        const int count = std::min(tile, qt-begin);
        for (int i = 0; i < count; ++i) {
            const int target = begin+i;
            for (int j = 0; j < qs; ++j) {
                const double dx = tx[3*target]-sx[3*j];
                const double dy = tx[3*target+1]-sx[3*j+1];
                const double dz = tx[3*target+2]-sx[3*j+2];
                const double radius = std::sqrt(dx*dx+dy*dy+dz*dz);
                kernel[static_cast<size_t>(i)*qs+j] = radius < 1e-300 ? 0.0 : sw[j]/radius;
            }
            for (int j = 0; j < nt; ++j)
                outer[static_cast<size_t>(i)*nt+j] = tw[target]*tv[static_cast<size_t>(target)*nt+j];
        }
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            count, ns, qs, 1.0, kernel.data(), qs, sv, ns,
            0.0, inner.data(), ns);
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
            nt, ns, count, 1.0, outer.data(), nt, inner.data(), ns,
            1.0, block, ns);
    }
}

} // namespace radia_detail
