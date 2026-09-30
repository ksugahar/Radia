#include "radia_ih_thermal.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace radia { namespace ih {
namespace {

void check_matrix(const CSRMatrix& a, int n, const char* name) {
    if (a.n != n || a.row_ptr.size() != static_cast<std::size_t>(n + 1) ||
        a.col.size() != a.value.size() ||
        a.row_ptr.front() != 0 || a.row_ptr.back() != static_cast<int>(a.col.size()))
        throw std::invalid_argument(std::string("invalid ") + name + " CSR matrix");
    for (int i = 0; i < n; ++i) {
        if (a.row_ptr[i] > a.row_ptr[i + 1])
            throw std::invalid_argument(std::string("non-monotone ") + name + " CSR rows");
        for (int k = a.row_ptr[i]; k < a.row_ptr[i + 1]; ++k)
            if (a.col[k] < 0 || a.col[k] >= n || !std::isfinite(a.value[k]))
                throw std::invalid_argument(std::string("invalid ") + name + " CSR entry");
    }
}

void matvec(const CSRMatrix& a, const std::vector<double>& x,
            std::vector<double>& y) {
    y.assign(static_cast<std::size_t>(a.n), 0.0);
    for (int i = 0; i < a.n; ++i)
        for (int k = a.row_ptr[i]; k < a.row_ptr[i + 1]; ++k)
            y[static_cast<std::size_t>(i)] +=
                a.value[static_cast<std::size_t>(k)] * x[static_cast<std::size_t>(a.col[static_cast<std::size_t>(k)])];
}

void add_scaled(const CSRMatrix& a, double scale, CSRMatrix& out) {
    out = a;
    for (double& value : out.value) value *= scale;
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) result += a[i] * b[i];
    return result;
}

double norm(const std::vector<double>& values) {
    return std::sqrt(std::max(0.0, dot(values, values)));
}

void true_residual(const CSRMatrix& a, const std::vector<double>& b,
                   const std::vector<double>& x, std::vector<double>& residual) {
    matvec(a, x, residual);
    for (int i = 0; i < a.n; ++i)
        residual[static_cast<std::size_t>(i)] =
            b[static_cast<std::size_t>(i)] - residual[static_cast<std::size_t>(i)];
}

void cg(const CSRMatrix& a, const std::vector<double>& b,
        double tolerance, int max_iterations, std::vector<double>& x) {
    const int n = a.n;
    if (!(tolerance > 0.0) || !std::isfinite(tolerance) || max_iterations < 0)
        throw std::invalid_argument("invalid IH thermal CG options");
    const double relative_limit = std::min(tolerance, 1.0e-6);

    std::vector<double> diagonal(static_cast<std::size_t>(n), 0.0);
    for (int i = 0; i < n; ++i)
        for (int k = a.row_ptr[i]; k < a.row_ptr[i + 1]; ++k)
            if (a.col[static_cast<std::size_t>(k)] == i)
                diagonal[static_cast<std::size_t>(i)] +=
                    a.value[static_cast<std::size_t>(k)];
    for (double value : diagonal)
        if (!(value > 0.0) || !std::isfinite(value))
            throw std::runtime_error(
                "IH thermal Jacobi preconditioner requires a positive finite diagonal");

    std::vector<double> r, z(static_cast<std::size_t>(n)), p, ap;
    true_residual(a, b, x, r); // Warm start from the previous accepted temperature.
    // The fixed effective load is b-A*x_initial.  Scaling by ||b|| would hide
    // an unresolved small heat increment behind the much larger stored-energy
    // term M*T_previous.
    const double effective_load_norm = norm(r);
    if (effective_load_norm == 0.0) return;
    const double target = relative_limit *
        std::max(effective_load_norm, 1.0e-300);
    if (norm(r) <= target) return;
    for (int i = 0; i < n; ++i)
        z[static_cast<std::size_t>(i)] =
            r[static_cast<std::size_t>(i)] / diagonal[static_cast<std::size_t>(i)];
    p = z;
    double rz = dot(r, z);
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        matvec(a, p, ap);
        const double pap = dot(p, ap);
        if (!(pap > 0.0) || !std::isfinite(pap))
            throw std::runtime_error("IH thermal matrix is not positive definite");
        const double alpha = rz / pap;
        for (int i = 0; i < n; ++i) {
            x[static_cast<std::size_t>(i)] += alpha * p[static_cast<std::size_t>(i)];
            r[static_cast<std::size_t>(i)] -= alpha * ap[static_cast<std::size_t>(i)];
        }
        if (norm(r) <= target) {
            true_residual(a, b, x, r);
            if (norm(r) <= target) return;
        }
        for (int i = 0; i < n; ++i)
            z[static_cast<std::size_t>(i)] =
                r[static_cast<std::size_t>(i)] / diagonal[static_cast<std::size_t>(i)];
        const double next_rz = dot(r, z);
        if (!(next_rz >= 0.0) || !std::isfinite(next_rz))
            throw std::runtime_error("IH thermal CG residual is not finite");
        const double beta = next_rz / rz;
        for (int i = 0; i < n; ++i)
            p[static_cast<std::size_t>(i)] = z[static_cast<std::size_t>(i)] +
                beta * p[static_cast<std::size_t>(i)];
        rz = next_rz;
    }
    true_residual(a, b, x, r);
    if (norm(r) > target)
        throw std::runtime_error("IH thermal CG true relative residual exceeds 1e-6");
}

}  // namespace

void advance_thermal(const CSRMatrix& mass, const CSRMatrix& stiffness,
                     const CSRMatrix* convection,
                     const std::vector<double>& source_W,
                     const std::vector<double>& cell_weights,
                     double angle_now_rad, const ThermalStepOptions& options,
                     ThermalState& state) {
    const int n = mass.n;
    if (n <= 0 || source_W.size() != static_cast<std::size_t>(n) ||
        cell_weights.size() != static_cast<std::size_t>(n) ||
        state.temperature_K.size() != static_cast<std::size_t>(n) ||
        !(options.dt_s > 0.0) || !std::isfinite(angle_now_rad))
        throw std::invalid_argument("invalid IH thermal step dimensions or time data");
    check_matrix(mass, n, "mass");
    check_matrix(stiffness, n, "stiffness");
    if (convection) check_matrix(*convection, n, "convection");
    const bool coefficients = !options.constant_coefficients.empty();
    if (coefficients && options.constant_coefficients.size() != static_cast<std::size_t>(n))
        throw std::invalid_argument("IH constant coefficient size mismatch");
    for (double w : cell_weights)
        if ((!coefficients && !(w > 0.0)) || !std::isfinite(w))
            throw std::invalid_argument("IH thermal cell weights must be finite and positive");

    // M + dt*K is assembled from the same workpiece mesh as the reference.
    // State and source are in material coordinates; no thermal-state transport.
    CSRMatrix system;
    add_scaled(stiffness, options.dt_s * options.conductivity_scale, system);
    if (system.row_ptr != mass.row_ptr || system.col != mass.col ||
        (convection && (system.row_ptr != convection->row_ptr || system.col != convection->col)))
        throw std::invalid_argument("mass and stiffness CSR sparsity must match");
    for (std::size_t k = 0; k < system.value.size(); ++k)
        system.value[k] += mass.value[k];

    std::vector<double> rhs;
    matvec(mass, state.temperature_K, rhs);
    for (int i = 0; i < n; ++i) {
        rhs[static_cast<std::size_t>(i)] += options.dt_s * source_W[static_cast<std::size_t>(i)];
        if (convection) {
            double row_sum = 0.0;
            for (int k = convection->row_ptr[i]; k < convection->row_ptr[i + 1]; ++k)
                row_sum += convection->value[static_cast<std::size_t>(k)] *
                    (coefficients ? options.constant_coefficients[
                        static_cast<std::size_t>(convection->col[static_cast<std::size_t>(k)])] : 1.0);
            rhs[static_cast<std::size_t>(i)] += options.dt_s * options.convection_W_per_m2K *
                options.ambient_temperature_K * row_sum;
        }
    }
    if (convection)
        for (std::size_t k = 0; k < system.value.size(); ++k)
            system.value[k] += options.dt_s * options.convection_W_per_m2K * convection->value[k];
    std::vector<double> accepted_temperature = state.temperature_K;
    cg(system, rhs, options.tolerance, options.max_iterations,
       accepted_temperature);
    state.temperature_K.swap(accepted_temperature);
    state.time_s += options.dt_s;
    state.previous_angle_rad = angle_now_rad;
}

}}  // namespace radia::ih
