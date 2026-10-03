// Validation-only adapter. Not installed or exported by Radia.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <sparsesolv/sparsesolv.hpp>
#include <chrono>

namespace py = pybind11;
using Clock = std::chrono::steady_clock;
using namespace sparsesolv;
static double elapsed(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

template<class T> struct TimedIC : ICPreconditioner<T> {
    mutable double seconds = 0;
    mutable int calls = 0;
    bool profile;
    explicit TimedIC(bool enabled) : profile(enabled) {}
    T apply_fused_dot(const T* r, const T* x, T* y, index_t n, bool c) const override {
        if (!profile) return ICPreconditioner<T>::apply_fused_dot(r, x, y, n, c);
        auto t = Clock::now();
        T rho = ICPreconditioner<T>::apply_fused_dot(r, x, y, n, c);
        seconds += elapsed(t); ++calls;
        return rho;
    }
};

template<class T> using Array = py::array_t<T, py::array::c_style | py::array::forcecast>;
template<class T> py::dict run(Array<index_t> rp, Array<index_t> ci, Array<T> av,
                             Array<T> rhs, bool conjugate, bool abmc,
                             bool scaling, bool profile, double tolerance,
                             int maxiter, bool auto_shift, double shift) {
    const index_t n = static_cast<index_t>(rhs.size());
    if (rp.size() != n + 1 || ci.size() != av.size() || rp.data()[0] != 0
        || rp.data()[n] != av.size()) throw py::value_error("invalid CSR dimensions");
    for (index_t i = 0; i < n; ++i) {
        if (rp.data()[i] > rp.data()[i+1] || rp.data()[i] < 0 || rp.data()[i+1] > av.size())
            throw py::value_error("invalid CSR row pointer");
        bool diagonal = false;
        index_t last = -1;
        for (index_t k = rp.data()[i]; k < rp.data()[i+1]; ++k) {
            if (ci.data()[k] <= last || ci.data()[k] >= n) throw py::value_error("CSR must be sorted/unique");
            diagonal |= ci.data()[k] == i;
            last = ci.data()[k];
        }
        if (!diagonal) throw py::value_error("missing diagonal");
    }
    SparseMatrixView<T> a(n, n, rp.data(), ci.data(), av.data());
    SolverConfig config;
    config.conjugate = conjugate;
    config.use_abmc = abmc;
    config.diagonal_scaling = scaling;
    config.tolerance = tolerance;
    config.max_iterations = maxiter;
    config.auto_shift = auto_shift;
    config.shift_parameter = shift;
    config.save_residual_history = true;
    std::vector<T> x(n, T(0));
    double setup = 0, solve = 0, apply = 0;
    int calls = 0;
    auto start = Clock::now();
    // Intentionally bypass only the public conjugate rejection in this isolated
    // experiment. Scaling, stopping, best iterate and true residual stay native.
    auto result = detail::solve_scaled(a, rhs.data(), x.data(), n, config,
        [&](const SparseMatrixView<T>& m, const T* b, T* y, const SolverConfig& cfg) {
            TimedIC<T> pre(profile);
            pre.set_config(cfg);
            auto t = Clock::now(); pre.setup(m); setup += elapsed(t);
            CGSolver<T> solver; solver.set_config(cfg);
            t = Clock::now(); auto r = solver.solve(m, b, y, m.rows(), &pre); solve += elapsed(t);
            apply += pre.seconds; calls += pre.calls;
            r.actual_shift = pre.actual_shift();
            return r;
        });
    const double total = elapsed(start);
    std::vector<T> temp(n);
    start = Clock::now();
    for (int i = 0; i < 24; ++i) a.multiply(x.data(), temp.data());
    const double spmv = elapsed(start) / 24;
    py::array_t<T> solution(n);
    std::copy(x.begin(), x.end(), solution.mutable_data());
    py::dict d;
    d["x"] = solution; d["iterations"] = result.iterations;
    d["converged"] = result.converged; d["true_residual"] = result.true_residual;
    d["actual_shift"] = result.actual_shift; d["history"] = result.residual_history;
    d["setup_s"] = setup; d["krylov_s"] = solve; d["total_s"] = total;
    d["apply_and_dot_s"] = apply; d["apply_calls"] = calls; d["spmv_s"] = spmv;
    return d;
}

PYBIND11_MODULE(STUDY_MODULE, m) {
    for (auto name : {"real", "complex"}) {
        if (std::string(name) == "real")
            m.def(name, &run<double>);
        else m.def(name, &run<std::complex<double>>);
    }
}
