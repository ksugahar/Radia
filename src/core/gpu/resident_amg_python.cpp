#ifndef RADIA_AMG_PYTHON_MODULE
#define RADIA_AMG_PYTHON_MODULE _gpu_resident_amg
#endif
// SPDX-License-Identifier: Apache-2.0
#include "resident_amg_host.hpp"
#include <algorithm>
#include <climits>
#include <cmath>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>
namespace py = pybind11;
template <class T> std::vector<T> array(py::handle o) {
  auto a = py::array_t<T, py::array::c_style>::ensure(o);
  if (!a || a.ndim() != 1)
    throw std::invalid_argument("contiguous typed vector required");
  return {a.data(), a.data() + a.size()};
}
AmgCSR csr(py::dict d) {
  AmgCSR a;
  a.n = d["n"].cast<int>();
  a.m = d["m"].cast<int>();
  a.row = array<int>(d["row"]);
  a.col = array<int>(d["col"]);
  a.val = array<double>(d["val"]);
  if (a.n < 1 || a.m < 1 || a.n > INT_MAX / 8 - 256 ||
      a.m > INT_MAX / 8 - 256 || a.row.size() != size_t(a.n) + 1 ||
      a.col.size() != a.val.size() || a.val.size() > 2147483647 ||
      a.row.front() != 0 || a.row.back() != a.val.size())
    throw std::invalid_argument("CSR dimensions");
  for (int i = 0; i < a.n; ++i) {
    if (a.row[i] < 0 || a.row[i + 1] < a.row[i] ||
        size_t(a.row[i + 1]) > a.val.size())
      throw std::invalid_argument("CSR rows");
    int last = -1;
    for (int j = a.row[i]; j < a.row[i + 1]; ++j) {
      if (a.col[j] <= last || a.col[j] >= a.m || !std::isfinite(a.val[j]))
        throw std::invalid_argument("CSR entries");
      last = a.col[j];
    }
  }
  return a;
}
struct Resident {
  void *p = nullptr;
  int n;
  Resident(py::list input, py::array_t<double, py::array::c_style> inverse) {
    std::vector<AmgLevel> levels;
    for (auto item : input) {
      auto d = py::cast<py::dict>(item);
      AmgLevel l;
      l.a = csr(d["a"].cast<py::dict>());
      if (l.a.n != l.a.m)
        throw std::invalid_argument("square level required");
      l.smoother = array<double>(d["smoother"]);
      if (l.smoother.size() != size_t(l.a.n))
        throw std::invalid_argument("smoother dimension");
      for (auto v : l.smoother)
        if (!std::isfinite(v) || v <= 0)
          throw std::invalid_argument("positive finite smoother");
      if (d.contains("p")) {
        l.p = csr(d["p"].cast<py::dict>());
        l.r = csr(d["r"].cast<py::dict>());
      }
      levels.push_back(std::move(l));
    }
    if (levels.empty() || levels.size() > 30)
      throw std::invalid_argument("hierarchy depth");
    n = levels.front().a.n;
    for (size_t i = 0; i + 1 < levels.size(); ++i)
      if (levels[i].p.n != levels[i].a.n ||
          levels[i].p.m != levels[i + 1].a.n ||
          levels[i].r.n != levels[i + 1].a.n || levels[i].r.m != levels[i].a.n)
        throw std::invalid_argument("transfer dimensions");
    int c = levels.back().a.n;
    if (c > 128 || inverse.size() != size_t(c) * c)
      throw std::invalid_argument("coarse inverse dimension");
    auto v = array<double>(inverse);
    for (auto q : v)
      if (!std::isfinite(q))
        throw std::invalid_argument("coarse nonfinite");
    p = amg_create(levels, v);
  }
  Resident(const Resident &) = delete;
  Resident &operator=(const Resident &) = delete;
  ~Resident() {
    if (p)
      amg_release(p);
  }
  void close() {
    if (p) {
      amg_release(p);
      p = nullptr;
    }
  }
  py::tuple solve(py::array_t<double, py::array::c_style> b, double tol,
                  int limit) {
    if (!p)
      throw std::runtime_error("released AMG context");
    if (b.ndim() != 1 || b.size() != n || !std::isfinite(tol) || tol <= 0 ||
        tol > 1e-6 || limit < 1 || limit > 100000)
      throw std::invalid_argument("solve inputs");
    for (auto v : array<double>(b))
      if (!std::isfinite(v))
        throw std::invalid_argument("RHS nonfinite");
    py::array_t<double> x(n);
    AmgStats s;
    s = amg_solve(p, b.data(), x.mutable_data(), tol, limit);
    py::dict info;
    info["backend"] = "gpu";
    info["preconditioner"] = "smoothed-aggregation AMG V-cycle";
    info["iterations"] = s.iterations;
    info["true_relative_residual"] = s.residual;
    info["certificate_gate"] = s.certificate_gate;
    info["owned_device_bytes"] = s.bytes;
    info["solve_count"] = s.solves;
    info["residual_refinements"] = s.residual_refinements;
    info["initial_original_residual"] = s.initial_original_residual;
    info["transfer_included"] = true;
    info["graph_batch_iterations"] = 8;
    info["iteration_host_syncs"] = s.graph_batches;
    info["matrix_and_work_vectors_resident"] = true;
    return py::make_tuple(x, info);
  }
};
PYBIND11_MODULE(RADIA_AMG_PYTHON_MODULE, m) {
  py::class_<Resident>(m, "Resident")
      .def(py::init<py::list, py::array_t<double, py::array::c_style>>())
      .def("solve", &Resident::solve, py::arg("rhs").noconvert(),
           py::arg("tolerance") = 1e-12, py::arg("max_iterations") = 10000)
      .def("close", &Resident::close);
}
