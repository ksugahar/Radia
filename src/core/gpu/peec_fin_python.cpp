// pybind11 entry for the device-resident fin PEEC solve (peec_fin_host.hpp).
#include "peec_fin_host.hpp"
#include "../peec_oriented_box.hpp"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <cmath>
#include <stdexcept>
namespace py = pybind11;
template <class T> using Array = py::array_t<T, py::array::c_style>;

template <class T> static std::vector<T> own(const Array<T>& a, py::ssize_t n, const char* name) {
    if (a.size() != n) throw py::value_error(std::string(name) + ": wrong length");
    std::vector<T> v(a.data(), a.data() + a.size());
    if constexpr (std::is_floating_point_v<T>)
        for (T x : v) if (!std::isfinite(x)) throw py::value_error(std::string(name) + ": non-finite entry");
    return v;
}

static py::tuple solve(Array<double> center, Array<double> direction, Array<double> width_dir,
                       Array<std::uint8_t> has_width_dir, Array<double> length, Array<double> width,
                       Array<double> height, Array<double> resistance, int order, double frequency, Array<int> branch_nodes,
                       Array<int> node_map, int n_reduced, Array<std::complex<double>> injection) {
    const py::ssize_t m = length.size();
    if (m <= 0) throw py::value_error("no bars");
    if (order < 1 || order > 32) throw py::value_error("order must be 1..32");
    const py::ssize_t nr = py::ssize_t(order) * (order + 1) / 2;
    auto c = own(center, 3 * m, "center"), d = own(direction, 3 * m, "direction"),
         w = own(width_dir, 3 * m, "width_dir"), l = own(length, m, "length"), wi = own(width, m, "width"),
         h = own(height, m, "height"), r = own(resistance, m, "resistance");
    std::vector<double> qx(nr), qw(nr);
    radia_peec_box::GaussRules(order, qx.data(), qw.data());
    auto hw = own(has_width_dir, m, "has_width_dir");
    auto bn = own(branch_nodes, 2 * m, "branch_nodes");
    const py::ssize_t nn = node_map.size();
    auto nm = own(node_map, nn, "node_map");
    for (int k : bn) if (k < 0 || k >= nn) throw py::value_error("branch node out of range");
    for (int k : nm) if (k < -1 || k >= n_reduced) throw py::value_error("node_map entry out of range");
    for (py::ssize_t i = 0; i < m; ++i)
        if (!(l[i] > 0 && wi[i] > 0 && h[i] > 0 && r[i] > 0)) throw py::value_error("bars need positive size and R");
    auto inj = own(injection, n_reduced, "injection");
    radia_peec_fin::Input in;
    in.m = int(m); in.center = c.data(); in.direction = d.data(); in.width_dir = w.data();
    in.has_width_dir = hw.data(); in.length = l.data(); in.width = wi.data(); in.height = h.data();
    in.resistance = r.data(); in.order = order; in.rule_x = qx.data(); in.rule_w = qw.data();
    in.frequency = frequency; in.branch_nodes = bn.data(); in.n_nodes = int(nn); in.node_map = nm.data();
    in.n_reduced = n_reduced; in.injection = inj.data();
    radia_peec_fin::Output o;
    {
        py::gil_scoped_release release;
        o = radia_peec_fin::solve_host(in);
    }
    Array<std::complex<double>> I(m), V(n_reduced);
    std::copy(o.current.begin(), o.current.end(), I.mutable_data());
    std::copy(o.voltage.begin(), o.voltage.end(), V.mutable_data());
    py::dict meta;
    meta["backend"] = "gpu"; meta["precision"] = "float64";
    meta["matrix_residency"] = "L and Z assembled, gated, factorised and solved on the device";
    meta["cholesky"] = o.cholesky; meta["min_eig_over_max"] = o.min_eig_over_max; meta["n_null"] = o.n_null;
    meta["max_abs_L"] = o.max_abs_L; meta["true_relative_residual"] = o.z_true_residual;
    meta["IRI"] = o.current_R_current; meta["ILI"] = o.current_L_current;
    meta["t_assemble_s"] = o.t_assemble_s; meta["t_energy_s"] = o.t_energy_s;
    meta["t_factor_s"] = o.t_factor_s; meta["t_solve_s"] = o.t_solve_s;
    meta["bytes_to_device"] = o.bytes_to_device; meta["bytes_to_host"] = o.bytes_to_host;
    return py::make_tuple(I, V, meta);
}

PYBIND11_MODULE(_peec_fin_cuda, m) {
    m.doc() = "Device-resident fin PEEC port solve; no import-time CUDA calls";
    m.def("solve", &solve, py::arg("center").noconvert(), py::arg("direction").noconvert(),
          py::arg("width_dir").noconvert(), py::arg("has_width_dir").noconvert(), py::arg("length").noconvert(),
          py::arg("width").noconvert(), py::arg("height").noconvert(), py::arg("resistance").noconvert(),
          py::arg("order"), py::arg("frequency"),
          py::arg("branch_nodes").noconvert(), py::arg("node_map").noconvert(), py::arg("n_reduced"),
          py::arg("injection").noconvert());
}
