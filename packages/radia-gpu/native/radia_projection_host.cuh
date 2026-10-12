#pragma once
// mexcuda host operation recovered from the retired MEX interpolation path.
// NGSolve owns FE interpolation. No TaskManager region is created here.
#include <gridfunction.hpp>
#include <postproc.hpp>
#include <chrono>
#include <memory>
#include <stdexcept>

namespace radia::cuda_projection {
struct Profile {
    double interpolation_seconds;
};
inline Profile Interpolate(
    const std::shared_ptr<ngfem::CoefficientFunction>& coefficient,
    ngcomp::GridFunction& gridfunction) {
    if (!coefficient) throw std::invalid_argument("Projection coefficient is null");
    ngstd::LocalHeap local_heap(1 << 20, "radia_cuda_gridfunction_interpolate");
    const auto begin = std::chrono::steady_clock::now();
    ngcomp::SetValues(coefficient, gridfunction, ngfem::VOL, nullptr, local_heap);
    return {std::chrono::duration<double>(
        std::chrono::steady_clock::now() - begin).count()};
}
} // namespace radia::cuda_projection
