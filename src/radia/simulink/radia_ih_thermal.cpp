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
        a.row_ptr.back() != static_cast<int>(a.col.size()))
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

void cg(const CSRMatrix& a, const std::vector<double>& b,
        double tolerance, int max_iterations, std::vector<double>& x) {
    const int n = a.n;
    std::vector<double> r = b, p = r, ap;
    x.assign(static_cast<std::size_t>(n), 0.0);
    double rr = dot(r, r);
    const double target = tolerance * tolerance * std::max(1.0, rr);
    for (int iteration = 0; iteration < max_iterations && rr > target; ++iteration) {
        matvec(a, p, ap);
        const double pap = dot(p, ap);
        if (!(pap > 0.0) || !std::isfinite(pap))
            throw std::runtime_error("IH thermal matrix is not positive definite");
        const double alpha = rr / pap;
        for (int i = 0; i < n; ++i) {
            x[static_cast<std::size_t>(i)] += alpha * p[static_cast<std::size_t>(i)];
            r[static_cast<std::size_t>(i)] -= alpha * ap[static_cast<std::size_t>(i)];
        }
        const double next_rr = dot(r, r);
        if (next_rr <= target) return;
        const double beta = next_rr / rr;
        for (int i = 0; i < n; ++i)
            p[static_cast<std::size_t>(i)] = r[static_cast<std::size_t>(i)] +
                beta * p[static_cast<std::size_t>(i)];
        rr = next_rr;
    }
    if (rr > target) throw std::runtime_error("IH thermal CG did not converge");
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
    for (double w : cell_weights)
        if (!(w > 0.0) || !std::isfinite(w))
            throw std::invalid_argument("IH thermal cell weights must be finite and positive");

    // M + dt*K is assembled from the same workpiece mesh as the Python
    // reference.  The right-hand side uses the previous workpiece state;
    // angle transport is performed by the Thermal S-Function before calling
    // this kernel, so dt ordering remains explicit and testable.
    CSRMatrix system;
    add_scaled(stiffness, options.dt_s * options.conductivity_scale, system);
    if (system.row_ptr.size() != mass.row_ptr.size() || system.col != mass.col ||
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
                row_sum += convection->value[static_cast<std::size_t>(k)];
            rhs[static_cast<std::size_t>(i)] += options.dt_s * options.convection_W_per_m2K *
                options.ambient_temperature_K * row_sum;
        }
    }
    if (convection)
        for (std::size_t k = 0; k < system.value.size(); ++k)
            system.value[k] += options.dt_s * options.convection_W_per_m2K * convection->value[k];
    cg(system, rhs, options.tolerance, options.max_iterations, state.temperature_K);
    state.time_s += options.dt_s;
    state.previous_angle_rad = angle_now_rad;
}

}}  // namespace radia::ih
