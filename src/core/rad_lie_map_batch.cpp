#include "rad_lie_map_batch.h"

#include <core/taskmanager.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <stdexcept>

namespace rad_lie {
namespace {

constexpr int kDim = 6;

struct Generator {
    const double* tensor;
    int degree;
};

inline double DotRow(const double* row, const double* z) {
    return row[0]*z[0] + row[1]*z[1] + row[2]*z[2]
         + row[3]*z[3] + row[4]*z[4] + row[5]*z[5];
}

// Contract the trailing tensor indices with z down to the (6,6) partial
// contraction C2, from which both the gradient and the Hessian follow with
// the degree-specific combinatorial factors of the Python reference.
void GradientHessian(const Generator& generator, const double* z,
                     double* gradient, double* hessian) {
    double c4[1296];
    double c3[216];
    double c2[36];
    const double* source = generator.tensor;
    if (generator.degree == 5) {
        for (int a = 0; a < 1296; ++a) c4[a] = DotRow(source + a*kDim, z);
        for (int a = 0; a < 216; ++a) c3[a] = DotRow(c4 + a*kDim, z);
        for (int a = 0; a < 36; ++a) c2[a] = DotRow(c3 + a*kDim, z);
        for (int a = 0; a < 36; ++a) hessian[a] = c2[a] / 6.0;
        for (int i = 0; i < kDim; ++i)
            gradient[i] = DotRow(c2 + i*kDim, z) / 24.0;
        return;
    }
    if (generator.degree == 4) {
        for (int a = 0; a < 216; ++a) c3[a] = DotRow(source + a*kDim, z);
        for (int a = 0; a < 36; ++a) c2[a] = DotRow(c3 + a*kDim, z);
        for (int a = 0; a < 36; ++a) hessian[a] = 0.5 * c2[a];
        for (int i = 0; i < kDim; ++i)
            gradient[i] = DotRow(c2 + i*kDim, z) / 6.0;
        return;
    }
    for (int a = 0; a < 36; ++a) c2[a] = DotRow(source + a*kDim, z);
    for (int a = 0; a < 36; ++a) hessian[a] = c2[a];
    for (int i = 0; i < kDim; ++i)
        gradient[i] = 0.5 * DotRow(c2 + i*kDim, z);
}

// Solve the 6x6 system a x = b by Gaussian elimination with partial
// pivoting; a and b are destroyed.  Returns false on a singular pivot.
bool Solve6(double* a, double* b, double* x) {
    int order[kDim];
    for (int i = 0; i < kDim; ++i) order[i] = i;
    for (int column = 0; column < kDim; ++column) {
        int pivot = column;
        double best = std::fabs(a[order[column]*kDim + column]);
        for (int row = column+1; row < kDim; ++row) {
            const double candidate = std::fabs(a[order[row]*kDim + column]);
            if (candidate > best) { best = candidate; pivot = row; }
        }
        if (!(best > 0.0)) return false;
        std::swap(order[column], order[pivot]);
        const int lead = order[column];
        const double inverse_pivot = 1.0 / a[lead*kDim + column];
        for (int row = column+1; row < kDim; ++row) {
            const int target = order[row];
            const double factor = a[target*kDim + column] * inverse_pivot;
            if (factor == 0.0) continue;
            a[target*kDim + column] = 0.0;
            for (int k = column+1; k < kDim; ++k)
                a[target*kDim + k] -= factor * a[lead*kDim + k];
            b[target] -= factor * b[lead];
        }
    }
    for (int row = kDim-1; row >= 0; --row) {
        const int lead = order[row];
        double accumulated = b[lead];
        for (int k = row+1; k < kDim; ++k)
            accumulated -= a[lead*kDim + k] * x[k];
        x[row] = accumulated / a[lead*kDim + row];
    }
    return true;
}

// One generator flow on a single state; mirrors
// apply_homogeneous_lie_generator (predictor, Newton residual, implicit
// midpoint, infinity-norm convergence test) step for step.
bool ApplyGeneratorFlow(const Generator& generator, const double* poisson,
                        double* value, int substeps, double tolerance,
                        int iteration_cap) {
    const double step = 1.0 / static_cast<double>(substeps);
    double gradient[kDim];
    double hessian[kDim*kDim];
    double initial[kDim];
    double trial[kDim];
    double midpoint[kDim];
    double system[kDim*kDim];
    double residual[kDim];
    double update[kDim];
    for (int sub = 0; sub < substeps; ++sub) {
        for (int i = 0; i < kDim; ++i) initial[i] = value[i];
        GradientHessian(generator, initial, gradient, hessian);
        for (int i = 0; i < kDim; ++i)
            trial[i] = initial[i] + step * DotRow(poisson + i*kDim, gradient);
        bool converged = false;
        for (int iteration = 0; iteration < iteration_cap; ++iteration) {
            for (int i = 0; i < kDim; ++i)
                midpoint[i] = 0.5 * (initial[i] + trial[i]);
            GradientHessian(generator, midpoint, gradient, hessian);
            for (int i = 0; i < kDim; ++i)
                residual[i] = trial[i] - initial[i]
                    - step * DotRow(poisson + i*kDim, gradient);
            for (int i = 0; i < kDim; ++i) {
                const double* p_row = poisson + i*kDim;
                for (int j = 0; j < kDim; ++j) {
                    double product = 0.0;
                    for (int k = 0; k < kDim; ++k)
                        product += p_row[k] * hessian[k*kDim + j];
                    system[i*kDim + j] = (i == j ? 1.0 : 0.0)
                        - 0.5 * step * product;
                }
            }
            if (!Solve6(system, residual, update)) return false;
            double worst = 0.0;
            double reference = 1.0;
            for (int i = 0; i < kDim; ++i) {
                trial[i] -= update[i];
                worst = std::max(worst, std::fabs(update[i]));
                reference = std::max(reference, std::fabs(trial[i]));
            }
            if (worst <= tolerance * reference) { converged = true; break; }
        }
        if (!converged) return false;
        for (int i = 0; i < kDim; ++i) value[i] = trial[i];
    }
    return true;
}

}  // namespace

void ApplyDragtFinnMapBatch(
    const double* R,
    const double* f3,
    const double* f4,
    const double* f5,
    const double* poisson,
    const double* states_in,
    double* states_out,
    std::size_t n_states,
    int substeps,
    double newton_tolerance,
    int maximum_newton_iterations) {
    if (!R || !f3 || !f4 || !poisson || (!states_in && n_states) ||
        (!states_out && n_states))
        throw std::invalid_argument(
            "lie map batch: required array pointer is null");
    if (substeps < 1 || maximum_newton_iterations < 1 ||
        !(newton_tolerance > 0.0) || !std::isfinite(newton_tolerance))
        throw std::invalid_argument(
            "Lie-flow integration controls are invalid");
    if (n_states == 0) return;

    Generator generators[3];
    int generator_count = 0;
    if (f5) generators[generator_count++] = {f5, 5};
    generators[generator_count++] = {f4, 4};
    generators[generator_count++] = {f3, 3};

    std::atomic<bool> failed{false};
    auto advance = [&](std::size_t index) {
        if (failed.load(std::memory_order_relaxed)) return;
        double value[kDim];
        const double* source = states_in + kDim*index;
        for (int i = 0; i < kDim; ++i) value[i] = source[i];
        for (int g = 0; g < generator_count; ++g) {
            if (!ApplyGeneratorFlow(generators[g], poisson, value, substeps,
                                    newton_tolerance,
                                    maximum_newton_iterations)) {
                failed.store(true, std::memory_order_relaxed);
                return;
            }
        }
        double* target = states_out + kDim*index;
        double mapped[kDim];
        for (int i = 0; i < kDim; ++i)
            mapped[i] = DotRow(R + i*kDim, value);
        for (int i = 0; i < kDim; ++i) target[i] = mapped[i];
    };

    constexpr std::size_t kSerialStates = 16;
    if (n_states <= kSerialStates) {
        for (std::size_t index = 0; index < n_states; ++index) advance(index);
    } else {
        // Callers arrive from Python with no active TaskManager job; stand
        // up the region here (no-op when one is already active).
        ngcore::RegionTaskManager task_manager;
        ngcore::ParallelFor(ngcore::IntRange(n_states), advance);
    }
    if (failed.load())
        throw std::runtime_error(
            "implicit-midpoint Lie flow did not converge");
}

}  // namespace rad_lie
