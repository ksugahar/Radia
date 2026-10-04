/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/// @file sparsesolv_python_export.hpp
/// @brief Pybind11 bindings: type registration + factory functions with auto-dispatch

#ifndef NGSOLVE_SPARSESOLV_PYTHON_EXPORT_HPP
#define NGSOLVE_SPARSESOLV_PYTHON_EXPORT_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "sparsesolv_precond.hpp"
#include "sparsesolv_solvers.hpp"
#include <comp.hpp>
#include <array>
#include <chrono>
#include <atomic>
#include <type_traits>

// Compact AMG/AMS (TaskManager-based, no external dependency)
#include "sparsesolv/preconditioners/compact_amg.hpp"
#include "sparsesolv/preconditioners/compact_ams.hpp"
#include "sparsesolv/preconditioners/complex_compact_ams.hpp"
#include "sparsesolv/ngsolve/ams_wirebasket.hpp"


namespace py = pybind11;

namespace ngla {

// ============================================================================
// Internal: SparseSolvResult (non-templated, called once)
// ============================================================================

inline void ExportSparseSolvResult_impl(py::module& m) {
  py::class_<SparseSolvResult>(m, "SparseSolvResult",
    "Result of a SparseSolv iterative solve.")
    .def_readonly("converged", &SparseSolvResult::converged,
        "Whether the recursive relative residual fell below tol")
    .def_readonly("iterations", &SparseSolvResult::iterations,
        "Number of iterations performed")
    .def_readonly("best_iteration", &SparseSolvResult::best_iteration,
        "Iteration of the returned iterate (0 = initial guess)")
    .def_readonly("final_residual", &SparseSolvResult::final_residual,
        "Recursive relative residual of the returned iterate (scaled system "
        "when diagonal_scaling is on)")
    .def_readonly("true_residual", &SparseSolvResult::true_residual,
        "||b - A x|| / ||b|| of the returned x on the original (free-DOF) system")
    .def_readonly("actual_shift", &SparseSolvResult::actual_shift,
        "IC shift used (0 when no IC factor was applied)")
    .def_readonly("residual_history", &SparseSolvResult::residual_history,
        "[initial, iteration 1, ...] recursive relative residuals "
        "(if save_residual_history enabled)")
    .def("__repr__", [](const SparseSolvResult& r) {
      return string("SparseSolvResult(converged=") +
             (r.converged ? "True" : "False") +
             ", iterations=" + std::to_string(r.iterations) +
             ", best_iteration=" + std::to_string(r.best_iteration) +
             ", residual=" + std::to_string(r.final_residual) +
             ", true_residual=" + std::to_string(r.true_residual) + ")";
    });
}

// ============================================================================
// Internal: Type registration (D/C suffix, no constructors — use factories)
// ============================================================================

template<typename SCAL>
void ExportSparseSolvTyped(py::module& m, const std::string& suffix) {

  // IC Preconditioner type registration (factory creates instances, methods accessible via downcast)
  {
    std::string cls_name = "ICPreconditioner" + suffix;
    py::class_<SparseSolvICPreconditioner<SCAL>,
               shared_ptr<SparseSolvICPreconditioner<SCAL>>,
               BaseMatrix>(m, cls_name.c_str())
      .def("Update", &SparseSolvICPreconditioner<SCAL>::Update)
      .def_property("shift",
          &SparseSolvICPreconditioner<SCAL>::GetShift,
          &SparseSolvICPreconditioner<SCAL>::SetShift)
      .def_property("use_abmc",
          &SparseSolvICPreconditioner<SCAL>::GetUseABMC,
          &SparseSolvICPreconditioner<SCAL>::SetUseABMC)
      .def_property("abmc_block_size",
          &SparseSolvICPreconditioner<SCAL>::GetABMCBlockSize,
          &SparseSolvICPreconditioner<SCAL>::SetABMCBlockSize)
      .def_property("abmc_num_colors",
          &SparseSolvICPreconditioner<SCAL>::GetABMCNumColors,
          &SparseSolvICPreconditioner<SCAL>::SetABMCNumColors)
      .def_property("diagonal_scaling",
          &SparseSolvICPreconditioner<SCAL>::GetDiagonalScaling,
          &SparseSolvICPreconditioner<SCAL>::SetDiagonalScaling);
  }

  // SparseSolv Solver type registration
  {
    std::string cls_name = "SparseSolvSolver" + suffix;
    py::class_<SparseSolvSolver<SCAL>,
               shared_ptr<SparseSolvSolver<SCAL>>,
               BaseMatrix>(m, cls_name.c_str())
      .def("Solve", [](SparseSolvSolver<SCAL>& self,
                       const BaseVector& rhs, BaseVector& sol) {
        return self.Solve(rhs, sol);
      }, py::arg("rhs"), py::arg("sol"))
      .def_property("method",
          &SparseSolvSolver<SCAL>::GetMethod,
          &SparseSolvSolver<SCAL>::SetMethod)
      .def_property("tol",
          &SparseSolvSolver<SCAL>::GetTolerance,
          &SparseSolvSolver<SCAL>::SetTolerance)
      .def_property("maxiter",
          &SparseSolvSolver<SCAL>::GetMaxIterations,
          &SparseSolvSolver<SCAL>::SetMaxIterations)
      .def_property("shift",
          &SparseSolvSolver<SCAL>::GetShift,
          &SparseSolvSolver<SCAL>::SetShift)
      .def_property("save_best_result",
          &SparseSolvSolver<SCAL>::GetSaveBestResult,
          &SparseSolvSolver<SCAL>::SetSaveBestResult)
      .def_property("save_residual_history",
          &SparseSolvSolver<SCAL>::GetSaveResidualHistory,
          &SparseSolvSolver<SCAL>::SetSaveResidualHistory)
      .def_property("printrates",
          &SparseSolvSolver<SCAL>::GetPrintRates,
          &SparseSolvSolver<SCAL>::SetPrintRates)
      .def_property("auto_shift",
          &SparseSolvSolver<SCAL>::GetAutoShift,
          &SparseSolvSolver<SCAL>::SetAutoShift)
      .def_property("diagonal_scaling",
          &SparseSolvSolver<SCAL>::GetDiagonalScaling,
          &SparseSolvSolver<SCAL>::SetDiagonalScaling)
      .def_property("divergence_check",
          &SparseSolvSolver<SCAL>::GetDivergenceCheck,
          &SparseSolvSolver<SCAL>::SetDivergenceCheck)
      .def_property("divergence_threshold",
          &SparseSolvSolver<SCAL>::GetDivergenceThreshold,
          &SparseSolvSolver<SCAL>::SetDivergenceThreshold)
      .def_property("divergence_count",
          &SparseSolvSolver<SCAL>::GetDivergenceCount,
          &SparseSolvSolver<SCAL>::SetDivergenceCount)
      .def_property("conjugate",
          &SparseSolvSolver<SCAL>::GetConjugate,
          &SparseSolvSolver<SCAL>::SetConjugate)
      .def_property("use_abmc",
          &SparseSolvSolver<SCAL>::GetUseABMC,
          &SparseSolvSolver<SCAL>::SetUseABMC)
      .def_property("abmc_block_size",
          &SparseSolvSolver<SCAL>::GetABMCBlockSize,
          &SparseSolvSolver<SCAL>::SetABMCBlockSize)
      .def_property("abmc_num_colors",
          &SparseSolvSolver<SCAL>::GetABMCNumColors,
          &SparseSolvSolver<SCAL>::SetABMCNumColors)
      .def_property("abmc_reorder_spmv",
          &SparseSolvSolver<SCAL>::GetABMCReorderSpMV,
          &SparseSolvSolver<SCAL>::SetABMCReorderSpMV)
      .def_property("abmc_use_rcm",
          &SparseSolvSolver<SCAL>::GetABMCUseRCM,
          &SparseSolvSolver<SCAL>::SetABMCUseRCM)
      .def_property_readonly("last_result",
          &SparseSolvSolver<SCAL>::GetLastResult);
  }
}

// ============================================================================
// Internal helpers
// ============================================================================

inline shared_ptr<BitArray> ExtractFreeDofs(py::object freedofs) {
  if (freedofs.is_none()) return nullptr;
  return py::cast<shared_ptr<BitArray>>(freedofs);
}

/// Preconditioned CG for an SPD system on free dofs, stopped on the true
/// relative residual |b - A x|_free / |b|_free. Iterates on the recurrence
/// residual and confirms with the true residual before stopping (continuing
/// from the true residual when they disagree), so it costs one matrix product
/// and one preconditioner application per iteration. Inner products are summed
/// over a fixed number of chunks in a fixed order: deterministic. The
/// preconditioner must return zero on constrained dofs (the AMS does).
class NativePCG {
public:
  NativePCG(shared_ptr<BaseMatrix> mat, shared_ptr<BaseMatrix> pre, shared_ptr<BitArray> freedofs)
      : mat_(mat), pre_(pre), free_(freedofs), r_(mat->CreateColVector()), z_(mat->CreateColVector()),
        p_(mat->CreateColVector()), q_(mat->CreateColVector())
  {
    n_ = mat_->Height();
    if (mat_->Width() != n_ || pre_->Height() != n_ || pre_->Width() != n_)
      throw py::value_error("NativePCG: square matrix and preconditioner of one size required");
    if (free_ && free_->Size() != n_) throw py::value_error("NativePCG: freedofs size mismatch");
  }

  py::tuple Solve(const BaseVector& b, BaseVector& x, double tolerance, int maxiter) {
    if (!(tolerance > 0.0 && tolerance < 1.0) || maxiter < 1)
      throw py::value_error("NativePCG.Solve: tolerance in (0, 1) and maxiter >= 1 required");
    auto fb = b.FVDouble(), fx = x.FVDouble();
    auto r = r_.FV<double>(), z = z_.FV<double>(), p = p_.FV<double>(), q = q_.FV<double>();
    const bool masked = bool(free_);
    auto is_free = [&](size_t i) { return !masked || free_->Test(i); };
    ParallelFor(n_, [&](size_t i) { fx[i] = 0.0; r[i] = is_free(i) ? fb[i] : 0.0; });
    const double bnorm = std::sqrt(Dot(r, r));
    if (bnorm == 0.0) return py::make_tuple(0, 0.0, true);
    pre_->Mult(r_, z_);
    ParallelFor(n_, [&](size_t i) { p[i] = z[i]; });
    double rz = Dot(r, z);
    if (!(rz > 0.0) || !std::isfinite(rz))
      throw std::runtime_error("NativePCG: r.z <= 0 (preconditioner not SPD on the free dofs)");
    const double target = tolerance * bnorm;
    double true_relative = 1.0;
    for (int it = 1; it <= maxiter; it++) {
      mat_->Mult(p_, q_);
      if (masked) ParallelFor(n_, [&](size_t i) { if (!free_->Test(i)) q[i] = 0.0; });
      const double pq = Dot(p, q);
      if (!(pq > 0.0) || !std::isfinite(pq))
        throw std::runtime_error("NativePCG: p.Ap <= 0 (matrix or preconditioner not SPD on the free dofs)");
      const double alpha = rz / pq;
      ParallelFor(n_, [&](size_t i) { fx[i] += alpha * p[i]; r[i] -= alpha * q[i]; });
      if (std::sqrt(Dot(r, r)) <= target) {
        // Confirm with the true residual; continue from it when it disagrees.
        mat_->Mult(x, q_);
        ParallelFor(n_, [&](size_t i) { r[i] = is_free(i) ? fb[i] - q[i] : 0.0; });
        true_relative = std::sqrt(Dot(r, r)) / bnorm;
        if (true_relative <= tolerance) return py::make_tuple(it, true_relative, true);
      }
      pre_->Mult(r_, z_);
      const double rz_new = Dot(r, z);
      if (!(rz_new > 0.0) || !std::isfinite(rz_new))
        throw std::runtime_error("NativePCG: r.z <= 0 (preconditioner not SPD on the free dofs)");
      const double beta = rz_new / rz;
      rz = rz_new;
      ParallelFor(n_, [&](size_t i) { p[i] = z[i] + beta * p[i]; });
    }
    mat_->Mult(x, q_);
    ParallelFor(n_, [&](size_t i) { r[i] = is_free(i) ? fb[i] - q[i] : 0.0; });
    return py::make_tuple(maxiter, std::sqrt(Dot(r, r)) / bnorm, false);
  }

private:
  static constexpr size_t chunks_ = 256;
  double Dot(FlatVector<double> a, FlatVector<double> b) const {
    std::array<double, chunks_> partial{};
    ParallelFor(chunks_, [&](size_t c) {
      const size_t begin = n_ * c / chunks_, end = n_ * (c + 1) / chunks_;
      double s = 0.0;
      for (size_t i = begin; i < end; i++) s += a[i] * b[i];
      partial[c] = s;
    });
    double s = 0.0;
    for (double v : partial) s += v;
    return s;
  }
  shared_ptr<BaseMatrix> mat_, pre_;
  shared_ptr<BitArray> free_;
  size_t n_ = 0;
  AutoVector r_, z_, p_, q_;
};

/// A-phi DC current of a closed conductor with one thick cut (see
/// radia.meshed_current.solve_closed_coil_current_phi for the method). Every
/// step works on the conductor only: its tets, faces, connectivity, the cut
/// and its completeness, the jump function, the P1 Laplacian on the conductor
/// vertices (assembled per row in ascending element order and solved by the
/// given NGSolve inverse type with one vertex gauged), the current density,
/// the weak-divergence gate and the cut-face flux.
inline py::dict ClosedCoilCurrentPhiImpl(shared_ptr<ngcomp::MeshAccess> ma,
                                         const std::vector<int>& materials, double current_A,
                                         std::array<double, 3> origin, std::array<double, 3> normal,
                                         double radius, const std::string& inverse)
{
  using Vec3 = std::array<double, 3>;
  py::dict timing;
  auto tic = std::chrono::steady_clock::now();
  auto lap = [&](const char* name) { auto now = std::chrono::steady_clock::now(); timing[name] = std::chrono::duration<double>(now - tic).count(); tic = now; };
  const size_t ne = ma->GetNE(ngfem::VOL);
  std::vector<char> wanted(*std::max_element(materials.begin(), materials.end()) + 1, 0);
  for (int m : materials) wanted[m] = 1;
  std::vector<int64_t> cells;
  for (size_t i = 0; i < ne; i++) {
    const int index = ma->GetElIndex(ngfem::ElementId(ngfem::VOL, i));
    if (index >= 0 && size_t(index) < wanted.size() && wanted[index]) cells.push_back(int64_t(i));
  }
  const size_t n = cells.size();
  if (n == 0) throw py::value_error("nonempty tetrahedral conductor required");
  std::vector<std::array<int, 4>> tets(n);
  std::atomic<size_t> bad{0};
  ParallelFor(n, [&](size_t c) {
    ngfem::ElementId ei(ngfem::VOL, cells[c]);
    auto vertices = ma->GetElVertices(ei);
    if (ma->GetElType(ei) != ngfem::ET_TET || vertices.Size() != 4) { bad++; return; }
    for (int k = 0; k < 4; k++) tets[c][k] = vertices[k];
  });
  if (bad) throw py::value_error("closed-coil current requires straight tet4 conductor cells");
  auto point = [&](int v) { auto p = ma->GetPoint<3>(v); return Vec3{p(0), p(1), p(2)}; };

  lap("conductor_s");
  // Faces: sorted vertex triples with their owners; unshared faces are walls.
  struct Face { std::array<int, 3> v; int64_t owner; };
  std::vector<Face> all(4 * n);
  static const int local[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
  ParallelFor(n, [&](size_t c) {
    for (int f = 0; f < 4; f++) {
      std::array<int, 3> t{tets[c][local[f][0]], tets[c][local[f][1]], tets[c][local[f][2]]};
      std::sort(t.begin(), t.end());
      all[4 * c + f] = Face{t, int64_t(c)};
    }
  });
  // Group equal triples within buckets of their smallest vertex (in parallel),
  // then list faces in bucket (= vertex) order: deterministic.
  const size_t nvert = ma->GetNV();
  TableCreator<int> face_creator(nvert);
  for (; !face_creator.Done(); face_creator++)
    for (size_t i = 0; i < all.size(); i++) face_creator.Add(all[i].v[0], int(i));
  Table<int> buckets = face_creator.MoveTable();
  std::vector<int8_t> kind(all.size(), 0);  // 1 shared (first owner), 2 wall, 0 skip
  std::vector<int64_t> partner(all.size(), -1);
  std::atomic<size_t> nonmanifold{0};
  ParallelFor(nvert, [&](size_t vertex) {
    auto bucket = buckets[vertex];
    std::sort(bucket.Data(), bucket.Data() + bucket.Size(), [&](int a, int b) {
      return all[a].v != all[b].v ? all[a].v < all[b].v : all[a].owner < all[b].owner; });
    for (size_t i = 0; i < bucket.Size();) {
      size_t j = i + 1;
      while (j < bucket.Size() && all[bucket[j]].v == all[bucket[i]].v) j++;
      if (j - i > 2) nonmanifold++;
      else if (j - i == 2) { kind[bucket[i]] = 1; partner[bucket[i]] = all[bucket[i + 1]].owner; }
      else kind[bucket[i]] = 2;
      i = j;
    }
  });
  if (nonmanifold) throw py::value_error("nonmanifold conductor face");
  std::vector<std::array<int, 3>> shared_faces, wall_faces;
  std::vector<int64_t> first, second;
  for (size_t vertex = 0; vertex < nvert; vertex++)
    for (int i : buckets[vertex]) {
      if (kind[i] == 1) { shared_faces.push_back(all[i].v); first.push_back(all[i].owner); second.push_back(partner[i]); }
      else if (kind[i] == 2) wall_faces.push_back(all[i].v);
    }
  {
    std::vector<int64_t> parent(n);
    for (size_t i = 0; i < n; i++) parent[i] = int64_t(i);
    auto find = [&](int64_t i) { while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; } return i; };
    for (size_t k = 0; k < first.size(); k++) {
      const int64_t ra = find(first[k]), rb = find(second[k]);
      if (ra != rb) parent[ra] = rb;
    }
    const int64_t root = find(0);
    for (size_t i = 1; i < n; i++)
      if (find(int64_t(i)) != root)
        throw py::value_error("conductor must be face-connected; split independent coils");
  }

  lap("faces_connectivity_s");
  // Cells near the cut: signed distance and in-plane radius of the centroid.
  std::vector<double> signed_distance(n);
  std::vector<char> near(n);
  ParallelFor(n, [&](size_t c) {
    Vec3 centroid{0, 0, 0};
    for (int k = 0; k < 4; k++) { auto p = point(tets[c][k]); for (int d = 0; d < 3; d++) centroid[d] += p[d]; }
    Vec3 offset;
    for (int d = 0; d < 3; d++) offset[d] = centroid[d] / 4.0 - origin[d];
    const double s = offset[0] * normal[0] + offset[1] * normal[1] + offset[2] * normal[2];
    double r2 = 0.0;
    for (int d = 0; d < 3; d++) { const double t = offset[d] - s * normal[d]; r2 += t * t; }
    signed_distance[c] = s;
    near[c] = std::sqrt(r2) < radius;
  });
  std::vector<size_t> crossing;
  for (size_t k = 0; k < shared_faces.size(); k++) {
    const auto a = first[k], b = second[k];
    if (near[a] && near[b] && ((signed_distance[a] < 0) != (signed_distance[b] < 0))) crossing.push_back(k);
  }
  if (crossing.empty())
    throw py::value_error("cut plane does not cross the conductor within cut_radius_m");
  // Completeness: every border edge of the cut sheet must lie on the wall.
  {
    std::vector<std::array<int, 2>> cut_edges, wall_edges;
    for (auto k : crossing) {
      const auto& f = shared_faces[k];
      cut_edges.push_back({f[0], f[1]}); cut_edges.push_back({f[0], f[2]}); cut_edges.push_back({f[1], f[2]});
    }
    for (const auto& f : wall_faces) {
      wall_edges.push_back({f[0], f[1]}); wall_edges.push_back({f[0], f[2]}); wall_edges.push_back({f[1], f[2]});
    }
    std::sort(cut_edges.begin(), cut_edges.end());
    std::sort(wall_edges.begin(), wall_edges.end());
    size_t loose = 0;
    for (size_t i = 0; i < cut_edges.size();) {
      size_t j = i + 1;
      while (j < cut_edges.size() && cut_edges[j] == cut_edges[i]) j++;
      if (j - i == 1 && !std::binary_search(wall_edges.begin(), wall_edges.end(), cut_edges[i])) loose++;
      i = j;
    }
    if (loose)
      throw py::value_error("cut does not span the conductor section (" + std::to_string(loose)
                            + " interior border edges)");
  }

  // Geometry, the jump function tau (one at cut vertices on the positive side)
  // and h = grad tau on its support cells.
  const size_t nv = ma->GetNV();
  std::vector<char> on_cut_vertex(nv, 0);
  for (auto k : crossing) for (int v : shared_faces[k]) on_cut_vertex[v] = 1;
  std::vector<double> volume(n);
  std::vector<std::array<Vec3, 4>> grads(n);
  std::vector<Vec3> jump(n);
  std::atomic<size_t> degenerate{0}, support_cells{0};
  ParallelFor(n, [&](size_t c) {
    const Vec3 x0 = point(tets[c][0]);
    double E[3][3];
    for (int i = 0; i < 3; i++) { const Vec3 xi = point(tets[c][i + 1]); for (int d = 0; d < 3; d++) E[i][d] = xi[d] - x0[d]; }
    const double det = E[0][0] * (E[1][1] * E[2][2] - E[1][2] * E[2][1])
                     - E[0][1] * (E[1][0] * E[2][2] - E[1][2] * E[2][0])
                     + E[0][2] * (E[1][0] * E[2][1] - E[1][1] * E[2][0]);
    volume[c] = std::abs(det) / 6.0;
    if (!(volume[c] > 0.0)) { degenerate++; return; }
    // grad lambda_i (i = 1..3) is column i of E^{-1}: the rows of E^{-T}.
    double inv[3][3];
    inv[0][0] = (E[1][1] * E[2][2] - E[1][2] * E[2][1]) / det;
    inv[0][1] = (E[0][2] * E[2][1] - E[0][1] * E[2][2]) / det;
    inv[0][2] = (E[0][1] * E[1][2] - E[0][2] * E[1][1]) / det;
    inv[1][0] = (E[1][2] * E[2][0] - E[1][0] * E[2][2]) / det;
    inv[1][1] = (E[0][0] * E[2][2] - E[0][2] * E[2][0]) / det;
    inv[1][2] = (E[0][2] * E[1][0] - E[0][0] * E[1][2]) / det;
    inv[2][0] = (E[1][0] * E[2][1] - E[1][1] * E[2][0]) / det;
    inv[2][1] = (E[0][1] * E[2][0] - E[0][0] * E[2][1]) / det;
    inv[2][2] = (E[0][0] * E[1][1] - E[0][1] * E[1][0]) / det;
    for (int i = 0; i < 3; i++)
      for (int d = 0; d < 3; d++) grads[c][i + 1][d] = inv[d][i];
    for (int d = 0; d < 3; d++) grads[c][0][d] = -(grads[c][1][d] + grads[c][2][d] + grads[c][3][d]);
    bool touches = false;
    for (int k = 0; k < 4; k++) touches = touches || on_cut_vertex[tets[c][k]];
    const bool support = near[c] && signed_distance[c] >= 0 && touches;
    Vec3 h{0, 0, 0};
    if (support) {
      support_cells++;
      for (int k = 0; k < 4; k++)
        if (on_cut_vertex[tets[c][k]]) for (int d = 0; d < 3; d++) h[d] += grads[c][k][d];
    }
    jump[c] = h;
  });
  if (degenerate) throw py::value_error("degenerate conductor tetrahedron");

  lap("cut_geometry_s");
  // P1 Laplacian on the conductor vertices (local numbering, ascending global).
  std::vector<int> local_of(nv, -1), global_of;
  for (size_t c = 0; c < n; c++) for (int k = 0; k < 4; k++) local_of[tets[c][k]] = 0;
  for (size_t v = 0; v < nv; v++) if (local_of[v] == 0) { local_of[v] = int(global_of.size()); global_of.push_back(int(v)); }
  const size_t m = global_of.size();
  Array<int> sizes(n);
  sizes = 4;
  Table<int> elements(sizes);
  for (size_t c = 0; c < n; c++) for (int k = 0; k < 4; k++) elements[c][k] = local_of[tets[c][k]];
  // The graph constructor sorts the element rows in place: after it, index the
  // element vertices through tets / local_of, never through this table.
  auto laplace = make_shared<SparseMatrix<double>>(m, m, elements, elements, false);
  laplace->AsVector() = 0.0;
  TableCreator<int> creator(m);
  for (; !creator.Done(); creator++)
    for (size_t c = 0; c < n; c++) for (int k = 0; k < 4; k++) creator.Add(local_of[tets[c][k]], int(4 * c + k));
  Table<int> vertex_slots = creator.MoveTable();
  VVector<double> rhs(m), phi(m);
  auto values = laplace->AsVector().FVDouble();
  auto dot = [](const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };
  ParallelFor(m, [&](size_t row) {
    double f = 0.0;
    for (int slot : vertex_slots[row]) {
      const size_t c = size_t(slot) / 4;
      const int a = slot % 4;
      for (int b = 0; b < 4; b++)
        values[laplace->GetPosition(row, local_of[tets[c][b]])] += volume[c] * dot(grads[c][a], grads[c][b]);
      f -= volume[c] * dot(jump[c], grads[c][a]);
    }
    rhs.FV()[row] = f;
  });
  auto free = make_shared<BitArray>(m);
  free->Set();
  free->Clear(local_of[tets[0][0]]);  // one potential gauge for the connected conductor
  lap("assembly_s");
  if (inverse == "iccg") {
    // IC(0)-preconditioned CG on the diagonally scaled, gauged SPD system to
    // a recursive relative residual of 1e-12; non-convergence fails loudly.
    SparseSolvSolver<double> solver(laplace, "ICCG", free, 1e-12, 20000, 1.0);
    phi.FV() = 0.0;
    const auto result = solver.Solve(rhs, phi);
    if (!result.converged)
      throw std::runtime_error("A-phi potential: ICCG did not converge (residual "
                               + std::to_string(result.final_residual) + ", true residual "
                               + std::to_string(result.true_residual) + ")");
  } else {
    laplace->SetInverseType(inverse);
    auto solver = laplace->InverseMatrix(free);
    solver->Mult(rhs, phi);
  }

  lap("solve_s");
  // Field E = grad phi + h, energy, current density J = scale E.
  py::array_t<double> density({py::ssize_t(n), py::ssize_t(3)});
  double* J = density.mutable_data();
  std::vector<Vec3> field(n);
  ParallelFor(n, [&](size_t c) {
    Vec3 e = jump[c];
    for (int k = 0; k < 4; k++) {
      const double value = phi.FV()[local_of[tets[c][k]]];
      for (int d = 0; d < 3; d++) e[d] += value * grads[c][k][d];
    }
    field[c] = e;
  });
  double energy = 0.0;
  for (size_t c = 0; c < n; c++) energy += volume[c] * dot(field[c], field[c]);
  if (!std::isfinite(energy) || !(energy > 0.0))
    throw py::value_error("cut does not drive a resolved closed current path");
  const double scale = -current_A / energy;  // flux of E through the cut along +n is -energy
  for (size_t c = 0; c < n; c++) for (int d = 0; d < 3; d++) J[3 * c + d] = scale * field[c][d];

  // Weak divergence against every conductor P1 function, relative to the drive.
  std::vector<double> divergence(m), drive(m);
  ParallelFor(m, [&](size_t row) {
    double dv = 0.0, dr = 0.0;
    for (int slot : vertex_slots[row]) {
      const size_t c = size_t(slot) / 4;
      const int a = slot % 4;
      const Vec3 jc{J[3 * c], J[3 * c + 1], J[3 * c + 2]};
      dv += volume[c] * dot(jc, grads[c][a]);
      dr += volume[c] * scale * dot(jump[c], grads[c][a]);
    }
    divergence[row] = dv;
    drive[row] = dr;
  });
  double dn = 0.0, rn = 0.0;
  for (size_t i = 0; i < m; i++) { dn += divergence[i] * divergence[i]; rn += drive[i] * drive[i]; }
  const double relative_divergence = rn > 0.0 ? std::sqrt(dn) / std::sqrt(rn) : 0.0;
  if (current_A != 0.0 && !(relative_divergence <= 1e-8))
    throw std::runtime_error("weak current divergence gate failed (" + std::to_string(relative_divergence) + ")");
  // Independent check: face-averaged flux through the cut sheet along +n.
  double face_flux = 0.0;
  for (auto k : crossing) {
    const auto& f = shared_faces[k];
    const Vec3 p0 = point(f[0]), p1 = point(f[1]), p2 = point(f[2]);
    const Vec3 u{p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]}, w{p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
    Vec3 area{0.5 * (u[1] * w[2] - u[2] * w[1]), 0.5 * (u[2] * w[0] - u[0] * w[2]), 0.5 * (u[0] * w[1] - u[1] * w[0])};
    const double sign = dot(area, normal) > 0 ? 1.0 : (dot(area, normal) < 0 ? -1.0 : 0.0);
    for (int d = 0; d < 3; d++) area[d] *= sign;
    const auto a = first[k], b = second[k];
    for (int d = 0; d < 3; d++) face_flux += 0.5 * (J[3 * a + d] + J[3 * b + d]) * area[d];
  }
  lap("current_checks_s");
  py::array_t<int64_t> element_numbers{std::vector<py::ssize_t>{py::ssize_t(n)}};
  std::copy(cells.begin(), cells.end(), element_numbers.mutable_data());
  py::array_t<int32_t> vertices{std::vector<py::ssize_t>{py::ssize_t(m)}};
  std::copy(global_of.begin(), global_of.end(), vertices.mutable_data());
  py::array_t<double> potential{std::vector<py::ssize_t>{py::ssize_t(m)}};
  for (size_t i = 0; i < m; i++) potential.mutable_data()[i] = phi.FV()[i];
  py::dict out;
  out["elements"] = element_numbers;
  out["density"] = density;
  out["vertices"] = vertices;
  out["phi"] = potential;
  out["cut_faces"] = crossing.size();
  out["support_cells"] = support_cells.load();
  out["phi_ndof"] = m - 1;
  out["relative_weak_divergence"] = relative_divergence;
  out["cut_face_flux_A"] = face_flux;
  out["timing"] = timing;
  return out;
}

/// Nonlinear magnetostatic state and residual of a lowest-order HCurl
/// curl-curl problem with element-constant flux, in one native pass.
/// Built from LowestOrderCurlSystem data of all elements and the nonlinear
/// ("iron") element numbers. Evaluate(x, ...) computes the element curl
/// c_e = C_e x[dofs_e]; on the iron B_e = c_e (+ the source mean), |B_e|,
/// nu_e and dH/dB by linear interpolation of the tabulated law (np.interp
/// semantics), q_e = (dH/dB - nu)/|B|^2 (0 below 1e-9 T); and the residual
/// r = sum_e vol_e C_e^T (nu_e c_e + [iron] (nu_e - nu0) Bs_e) - f,
/// gathered per dof over its elements in ascending order (deterministic).
class LowestOrderCurlResidual {
public:
  LowestOrderCurlResidual(py::array_t<int32_t, py::array::c_style | py::array::forcecast> dofs,
                          py::array_t<double, py::array::c_style | py::array::forcecast> curl,
                          py::array_t<double, py::array::c_style | py::array::forcecast> volume,
                          size_t ndof,
                          py::array_t<int64_t, py::array::c_style | py::array::forcecast> iron,
                          double nu0)
      : ndof_(ndof), nu0_(nu0)
  {
    ne_ = volume.ndim() == 1 ? size_t(volume.shape(0)) : 0;
    if (dofs.ndim() != 2 || size_t(dofs.shape(0)) != ne_ || dofs.shape(1) != 6
        || curl.ndim() != 3 || size_t(curl.shape(0)) != ne_ || curl.shape(1) != 3 || curl.shape(2) != 6
        || iron.ndim() != 1)
      throw py::value_error("LowestOrderCurlResidual: dofs (ne, 6), curl (ne, 3, 6), volume (ne,), iron (n,)");
    dofs_.assign(dofs.data(), dofs.data() + 6 * ne_);
    curl_.assign(curl.data(), curl.data() + 18 * ne_);
    volume_.assign(volume.data(), volume.data() + ne_);
    for (auto d : dofs_)
      if (d < 0 || size_t(d) >= ndof_) throw py::value_error("LowestOrderCurlResidual: dof outside ndof");
    iron_.assign(iron.data(), iron.data() + iron.shape(0));
    iron_slot_.assign(ne_, -1);
    for (size_t k = 0; k < iron_.size(); k++) {
      if (iron_[k] < 0 || size_t(iron_[k]) >= ne_ || iron_slot_[iron_[k]] >= 0)
        throw py::value_error("LowestOrderCurlResidual: iron elements must be distinct element numbers");
      iron_slot_[iron_[k]] = int(k);
    }
    TableCreator<int> creator(ndof_);
    for (; !creator.Done(); creator++)
      for (size_t e = 0; e < ne_; e++)
        for (int a = 0; a < 6; a++)
          creator.Add(dofs_[6 * e + a], int(6 * e + a));
    dof_slots_ = creator.MoveTable();
    element_curl_.resize(3 * ne_);
    loads_.resize(3 * ne_);
  }

  py::tuple Evaluate(py::array_t<double, py::array::c_style | py::array::forcecast> x,
                     py::object source_mean,
                     py::array_t<double, py::array::c_style | py::array::forcecast> grid,
                     py::array_t<double, py::array::c_style | py::array::forcecast> nu_table,
                     py::array_t<double, py::array::c_style | py::array::forcecast> dhdb_table,
                     py::array_t<double, py::array::c_style | py::array::forcecast> load,
                     py::array_t<double, py::array::c_style> residual)
  {
    const size_t n = iron_.size();
    if (x.ndim() != 1 || size_t(x.shape(0)) != ndof_ || load.ndim() != 1 || size_t(load.shape(0)) != ndof_
        || residual.ndim() != 1 || size_t(residual.shape(0)) != ndof_)
      throw py::value_error("LowestOrderCurlResidual.Evaluate: x, load and residual need ndof entries");
    const size_t m = grid.ndim() == 1 ? size_t(grid.shape(0)) : 0;
    if (m < 2 || nu_table.ndim() != 1 || size_t(nu_table.shape(0)) != m
        || dhdb_table.ndim() != 1 || size_t(dhdb_table.shape(0)) != m)
      throw py::value_error("LowestOrderCurlResidual.Evaluate: grid, nu and dH/dB tables of one length >= 2");
    const double* sm = nullptr;
    py::array_t<double, py::array::c_style | py::array::forcecast> sm_array;
    if (!source_mean.is_none()) {
      sm_array = source_mean.cast<py::array_t<double, py::array::c_style | py::array::forcecast>>();
      if (sm_array.ndim() != 2 || size_t(sm_array.shape(0)) != n || sm_array.shape(1) != 3)
        throw py::value_error("LowestOrderCurlResidual.Evaluate: source_mean (n_iron, 3)");
      sm = sm_array.data();
    }
    const double* xv = x.data();
    const double* g = grid.data();
    const double* nt = nu_table.data();
    const double* dt = dhdb_table.data();
    const double* f = load.data();
    double* r = residual.mutable_data();
    py::array_t<double> b({py::ssize_t(n), py::ssize_t(3)});
    py::array_t<double> magnitude{std::vector<py::ssize_t>{py::ssize_t(n)}};
    py::array_t<double> nu{std::vector<py::ssize_t>{py::ssize_t(n)}};
    py::array_t<double> q{std::vector<py::ssize_t>{py::ssize_t(n)}};
    double* bv = b.mutable_data();
    double* mv = magnitude.mutable_data();
    double* nv = nu.mutable_data();
    double* qv = q.mutable_data();
    auto interp = [&](double value, const double* table) {
      // numpy.interp for an ascending grid.
      if (!(value > g[0])) return table[0];
      if (!(value < g[m - 1])) return table[m - 1];
      const size_t j = size_t(std::upper_bound(g, g + m, value) - g) - 1;
      const double slope = (table[j + 1] - table[j]) / (g[j + 1] - g[j]);
      return slope * (value - g[j]) + table[j];
    };
    ParallelFor(ne_, [&](size_t e) {
      const double* C = curl_.data() + 18 * e;
      const int32_t* d = dofs_.data() + 6 * e;
      for (int k = 0; k < 3; k++) {
        double s = 0.0;
        for (int a = 0; a < 6; a++) s += C[6 * k + a] * xv[d[a]];
        element_curl_[3 * e + k] = s;
      }
      const int slot = iron_slot_[e];
      if (slot < 0) {
        for (int k = 0; k < 3; k++) loads_[3 * e + k] = volume_[e] * (nu0_ * element_curl_[3 * e + k]);
        return;
      }
      double be[3];
      for (int k = 0; k < 3; k++) be[k] = element_curl_[3 * e + k] + (sm ? sm[3 * slot + k] : 0.0);
      const double mag = std::sqrt(be[0] * be[0] + be[1] * be[1] + be[2] * be[2]);
      const double nue = interp(mag, nt);
      const double dhdb = interp(mag, dt);
      for (int k = 0; k < 3; k++) bv[3 * slot + k] = be[k];
      mv[slot] = mag;
      nv[slot] = nue;
      qv[slot] = mag > 1.0e-9 ? (dhdb - nue) / (mag * mag) : 0.0;
      for (int k = 0; k < 3; k++)
        loads_[3 * e + k] = volume_[e] * (nue * element_curl_[3 * e + k]
                                          + (sm ? (nue - nu0_) * sm[3 * slot + k] : 0.0));
    });
    ParallelFor(ndof_, [&](size_t row) {
      double s = 0.0;
      for (int slot : dof_slots_[row]) {
        const size_t e = size_t(slot) / 6;
        const int a = slot % 6;
        const double* C = curl_.data() + 18 * e;
        s += C[a] * loads_[3 * e] + C[6 + a] * loads_[3 * e + 1] + C[12 + a] * loads_[3 * e + 2];
      }
      r[row] = s - f[row];
    });
    return py::make_tuple(b, magnitude, nu, q);
  }

private:
  size_t ndof_, ne_ = 0;
  double nu0_;
  std::vector<int32_t> dofs_;
  std::vector<double> curl_, volume_, element_curl_, loads_;
  std::vector<int64_t> iron_;
  std::vector<int> iron_slot_;
  Table<int> dof_slots_;
};

/// Newton Jacobian refresh of a lowest-order HCurl curl-curl matrix on a set
/// of elements with the tangent reluctivity nu I + q b b^T (element-constant).
/// Built from LowestOrderCurlSystem data of those elements; the values of every
/// row they touch are saved at construction (the constant part) and each
/// Refresh rewrites exactly those rows as constant part + sum over the
/// row's elements, in ascending element order, of
/// vol (nu C_a.C_b + q (b.C_a)(b.C_b)). Entries (i, j) and (j, i) sum
/// identical terms in the same order, so the matrix stays exactly symmetric.
class LowestOrderCurlJacobian {
public:
  LowestOrderCurlJacobian(shared_ptr<BaseMatrix> matrix,
                          py::array_t<int32_t, py::array::c_style | py::array::forcecast> dofs,
                          py::array_t<double, py::array::c_style | py::array::forcecast> curl,
                          py::array_t<double, py::array::c_style | py::array::forcecast> volume,
                          py::array_t<int32_t, py::array::c_style | py::array::forcecast> positions)
  {
    matrix_ = dynamic_pointer_cast<SparseMatrix<double>>(matrix);
    if (!matrix_) throw py::type_error("LowestOrderCurlJacobian: needs a real SparseMatrix");
    n_ = volume.ndim() == 1 ? size_t(volume.shape(0)) : 0;
    if (dofs.ndim() != 2 || size_t(dofs.shape(0)) != n_ || dofs.shape(1) != 6
        || curl.ndim() != 3 || size_t(curl.shape(0)) != n_ || curl.shape(1) != 3 || curl.shape(2) != 6
        || positions.ndim() != 2 || size_t(positions.shape(0)) != n_ || positions.shape(1) != 36)
      throw py::value_error("LowestOrderCurlJacobian: dofs (n, 6), curl (n, 3, 6), volume (n,), "
                            "positions (n, 36) of the same n");
    dofs_.assign(dofs.data(), dofs.data() + 6 * n_);
    curl_.assign(curl.data(), curl.data() + 18 * n_);
    volume_.assign(volume.data(), volume.data() + n_);
    positions_.assign(positions.data(), positions.data() + 36 * n_);
    const size_t height = matrix_->Height();
    const size_t nze = matrix_->NZE();
    for (size_t i = 0; i < 6 * n_; i++)
      if (dofs_[i] < 0 || size_t(dofs_[i]) >= height)
        throw py::value_error("LowestOrderCurlJacobian: dof outside the matrix");
    for (size_t i = 0; i < 36 * n_; i++)
      if (positions_[i] < 0 || size_t(positions_[i]) >= nze)
        throw py::value_error("LowestOrderCurlJacobian: position outside the matrix");
    // Positions must be the (dof_a, dof_b) entries: check the columns and rows.
    std::atomic<size_t> wrong{0};
    ParallelFor(n_, [&](size_t e) {
      for (int a = 0; a < 6; a++) {
        auto first = matrix_->First(dofs_[6 * e + a]);
        auto next = matrix_->First(dofs_[6 * e + a] + 1);
        auto cols = matrix_->GetRowIndices(dofs_[6 * e + a]);
        for (int b = 0; b < 6; b++) {
          const size_t p = size_t(positions_[36 * e + 6 * a + b]);
          if (p < first || p >= next || cols[p - first] != dofs_[6 * e + b]) wrong++;
        }
      }
    });
    if (wrong) throw py::value_error("LowestOrderCurlJacobian: positions do not match the dofs");
    TableCreator<int> creator(height);
    for (; !creator.Done(); creator++)
      for (size_t e = 0; e < n_; e++)
        for (int a = 0; a < 6; a++)
          creator.Add(dofs_[6 * e + a], int(6 * e + a));
    row_slots_ = creator.MoveTable();
    for (size_t row = 0; row < height; row++)
      if (row_slots_[row].Size()) rows_.push_back(int(row));
    // Store each lower/upper pair once. Fast-math vectorization may round
    // separately gathered transposed entries differently; PCG's matrix
    // contract requires a symmetric operator on every supported build.
    for (int row : rows_) {
      const size_t first = matrix_->First(row);
      auto columns = matrix_->GetRowIndices(row);
      for (int k = 0; k < columns.Size(); k++) {
        const int col = columns[k];
        if (col >= row) continue;
        const size_t transpose = matrix_->GetPositionTest(col, row);
        if (transpose == std::numeric_limits<size_t>::max())
          throw py::value_error("LowestOrderCurlJacobian: nonsymmetric sparsity pattern");
        symmetric_entries_.emplace_back(first + k, transpose);
      }
    }
    auto values = matrix_->AsVector().FVDouble();
    base_offset_.resize(rows_.size() + 1);
    base_offset_[0] = 0;
    for (size_t r = 0; r < rows_.size(); r++)
      base_offset_[r + 1] = base_offset_[r] + (matrix_->First(rows_[r] + 1) - matrix_->First(rows_[r]));
    base_.resize(base_offset_.back());
    ParallelFor(rows_.size(), [&](size_t r) {
      const size_t first = matrix_->First(rows_[r]);
      for (size_t k = 0; k < base_offset_[r + 1] - base_offset_[r]; k++)
        base_[base_offset_[r] + k] = values[first + k];
    });
  }

  void Refresh(py::array_t<double, py::array::c_style | py::array::forcecast> nu,
               py::array_t<double, py::array::c_style | py::array::forcecast> q,
               py::array_t<double, py::array::c_style | py::array::forcecast> b)
  {
    if (nu.ndim() != 1 || size_t(nu.shape(0)) != n_ || q.ndim() != 1 || size_t(q.shape(0)) != n_
        || b.ndim() != 2 || size_t(b.shape(0)) != n_ || b.shape(1) != 3)
      throw py::value_error("LowestOrderCurlJacobian.Refresh: nu (n,), q (n,), b (n, 3)");
    const double* nuv = nu.data();
    const double* qv = q.data();
    const double* bv = b.data();
    auto values = matrix_->AsVector().FVDouble();
    ParallelFor(rows_.size(), [&](size_t r) {
      const int row = rows_[r];
      const size_t first = matrix_->First(row);
      for (size_t k = 0; k < base_offset_[r + 1] - base_offset_[r]; k++)
        values[first + k] = base_[base_offset_[r] + k];
      for (int slot : row_slots_[row]) {
        const size_t e = size_t(slot) / 6;
        const int a = slot % 6;
        const double* C = curl_.data() + 18 * e;
        const double* be = bv + 3 * e;
        const double vol = volume_[e], n = nuv[e], qq = qv[e];
        const double ta = be[0] * C[a] + be[1] * C[6 + a] + be[2] * C[12 + a];
        for (int c = 0; c < 6; c++) {
          const double s = C[a] * C[c] + C[6 + a] * C[6 + c] + C[12 + a] * C[12 + c];
          const double tc = be[0] * C[c] + be[1] * C[6 + c] + be[2] * C[12 + c];
          values[positions_[36 * e + 6 * a + c]] += vol * (n * s + qq * (ta * tc));
        }
      }
    });
    // The row gather is complete; upper entries are read-only in this pass
    // and every lower destination is unique, so there are no write races.
    ParallelFor(symmetric_entries_.size(), [&](size_t k) {
      const auto& pair = symmetric_entries_[k];
      values[pair.first] = values[pair.second];
    });
  }

  size_t NumElements() const { return n_; }
  size_t NumRows() const { return rows_.size(); }

private:
  shared_ptr<SparseMatrix<double>> matrix_;
  size_t n_ = 0;
  std::vector<int32_t> dofs_, positions_;
  std::vector<double> curl_, volume_, base_;
  std::vector<size_t> base_offset_;
  std::vector<int> rows_;
  std::vector<std::pair<size_t, size_t>> symmetric_entries_;
  Table<int> row_slots_;
};

// ============================================================================
// Internal: Factory functions with auto-dispatch via mat->IsComplex()
// ============================================================================

inline void ExportSparseSolvFactories(py::module& m) {

  // ---- ICPreconditioner factory ----
  m.def("ICPreconditioner", [](shared_ptr<BaseMatrix> mat,
                                py::object freedofs, double shift) {
    auto sp_freedofs = ExtractFreeDofs(freedofs);
    shared_ptr<BaseMatrix> result;
    if (mat->IsComplex()) {
      auto sp = dynamic_pointer_cast<SparseMatrix<Complex>>(mat);
      if (!sp) throw py::type_error("ICPreconditioner: expected SparseMatrix");
      auto p = make_shared<SparseSolvICPreconditioner<Complex>>(sp, sp_freedofs, shift);
      p->Update();
      result = p;
    } else {
      auto sp = dynamic_pointer_cast<SparseMatrix<double>>(mat);
      if (!sp) throw py::type_error("ICPreconditioner: expected SparseMatrix");
      auto p = make_shared<SparseSolvICPreconditioner<double>>(sp, sp_freedofs, shift);
      p->Update();
      result = p;
    }
    return result;
  },
  py::arg("mat"), py::arg("freedofs") = py::none(), py::arg("shift") = 1.05,
  R"raw_string(
Incomplete Cholesky (IC) Preconditioner.

Parameters:

mat : SparseMatrix
  Real SPD or complex-symmetric matrix (transpose IC, auto-detected type).
  For Hermitian positive-definite systems use SparseSolvSolver with
  method="ICCG", conjugate=True instead.
freedofs : BitArray, optional
  Free DOFs. Constrained DOFs treated as identity.
shift : float
  Shift parameter (default: 1.05).
)raw_string");

  // ---- SparseSolvSolver factory ----
  m.def("SparseSolvSolver", [](shared_ptr<BaseMatrix> mat,
                                 const string& method, py::object freedofs,
                                 double tol, int maxiter, double shift,
                                 bool save_best_result, bool save_residual_history,
                                 bool printrates, bool conjugate,
                                 bool use_abmc, int abmc_block_size, int abmc_num_colors,
                                 bool abmc_reorder_spmv, bool abmc_use_rcm,
                                 bool auto_shift, bool diagonal_scaling,
                                 bool divergence_check, double divergence_threshold,
                                 int divergence_count) {
    auto sp_freedofs = ExtractFreeDofs(freedofs);

    auto configure = [&](auto& solver) {
      solver->SetAutoShift(auto_shift);
      solver->SetDiagonalScaling(diagonal_scaling);
      solver->SetDivergenceCheck(divergence_check);
      solver->SetDivergenceThreshold(divergence_threshold);
      solver->SetDivergenceCount(divergence_count);
      solver->SetConjugate(conjugate);
      solver->SetUseABMC(use_abmc);
      solver->SetABMCBlockSize(abmc_block_size);
      solver->SetABMCNumColors(abmc_num_colors);
      solver->SetABMCReorderSpMV(abmc_reorder_spmv);
      solver->SetABMCUseRCM(abmc_use_rcm);
    };

    shared_ptr<BaseMatrix> result;
    if (mat->IsComplex()) {
      auto sp = dynamic_pointer_cast<SparseMatrix<Complex>>(mat);
      if (!sp) throw py::type_error("SparseSolvSolver: expected SparseMatrix");
      auto solver = make_shared<SparseSolvSolver<Complex>>(
          sp, method, sp_freedofs, tol, maxiter, shift,
          save_best_result, save_residual_history, printrates);
      configure(solver);
      result = solver;
    } else {
      auto sp = dynamic_pointer_cast<SparseMatrix<double>>(mat);
      if (!sp) throw py::type_error("SparseSolvSolver: expected SparseMatrix");
      auto solver = make_shared<SparseSolvSolver<double>>(
          sp, method, sp_freedofs, tol, maxiter, shift,
          save_best_result, save_residual_history, printrates);
      configure(solver);
      result = solver;
    }
    return result;
  },
  py::arg("mat"),
  py::arg("method") = "ICCG",
  py::arg("freedofs") = py::none(),
  py::arg("tol") = 1e-8,
  py::arg("maxiter") = 0,
  py::arg("shift") = 1.0,
  py::arg("save_best_result") = true,
  py::arg("save_residual_history") = false,
  py::arg("printrates") = false,
  py::arg("conjugate") = false,
  py::arg("use_abmc") = false,
  py::arg("abmc_block_size") = 4,
  py::arg("abmc_num_colors") = 4,
  py::arg("abmc_reorder_spmv") = false,
  py::arg("abmc_use_rcm") = false,
  py::arg("auto_shift") = true,
  py::arg("diagonal_scaling") = true,
  py::arg("divergence_check") = true,
  py::arg("divergence_threshold") = 10.0,
  py::arg("divergence_count") = 10,
  R"raw_string(
Iterative solver (ICCG / CG / COCR). Auto-detects real/complex.

Can be used as BaseMatrix (inverse operator, zero initial guess) or via
Solve(rhs, sol) with sol as the initial guess for detailed results.

Parameters:

mat : SparseMatrix
  System matrix (real or complex), full (non-symmetric) storage.
method : str
  "ICCG", "CG", or "COCR".
  COCR: Conjugate Orthogonal Conjugate Residual for complex-symmetric A^T=A.
  For non-symmetric systems, use GMRESSolver instead.
freedofs : BitArray, optional
  Free DOFs.
tol : float
  Stop when the recursive relative residual ||r||/||b|| < tol (default 1e-8).
  With diagonal_scaling the test uses the scaled system.
maxiter : int
  Iteration limit; 0 means 2*n (default 0).
shift : float
  IC diagonal shift alpha >= 1 (alpha*a_ii on rows with Re(a_ii) > 0);
  the start value when auto_shift is on (default 1.0).
save_best_result : bool
  Return the iterate with the smallest residual, the initial guess
  included (default True); False returns the last iterate.
save_residual_history : bool
  Record residual history (default: False).
printrates : bool
  Print convergence info (default: False).
conjugate : bool
  Hermitian products for CG/ICCG (default False); ICCG uses adjoint IC.
  ICCG checks Hermitian structure; the caller must ensure positive definiteness.
  Encountered non-positive curvature raises. COCR rejects True.
auto_shift : bool
  Restart the IC factorization with shift + 0.01 while a pivot has
  Re(d) < 1e-6 |a_ii| and shift < 5; a pivot still below it at the limit
  is an error (default True).
diagonal_scaling : bool
  Solve (S A S) y = S b, x = S y with S = diag(1/sqrt|a_ii|) (default True).
divergence_check, divergence_threshold, divergence_count : bool, float, int
  Stop when more than divergence_count consecutive iterations neither set a
  new best residual nor stay below best*divergence_threshold
  (default True, 10, 10).
use_abmc, abmc_block_size, abmc_num_colors, abmc_reorder_spmv, abmc_use_rcm :
  Parallel IC ordering (changes the IC factor).
)raw_string");
}

// ============================================================================
// Compact AMG/AMS factory (TaskManager-based)
// ============================================================================

inline void ExportHypreBasedAMS(py::module& m) {

  // Type registration for CompactAMG (needed for correct virtual dispatch in Python)
  py::class_<CompactAMG, shared_ptr<CompactAMG>, BaseMatrix>
      (m, "CompactAMGPreconditionerImpl");

  // Type registrations for HypreBasedAMS / ComplexHypreBasedAMS (enables Update() from Python)
  py::class_<HypreBasedAMS, shared_ptr<HypreBasedAMS>, BaseMatrix>
      (m, "HypreBasedAMSPreconditionerImpl")
      .def_property_readonly("beta_zero", &HypreBasedAMS::GetBetaZero)
      .def_property_readonly("reuse_hierarchy", &HypreBasedAMS::GetReuseHierarchy)
      .def_property_readonly("mixed_precision", &HypreBasedAMS::GetMixedPrecision)
      .def_property_readonly("hierarchy_refreshes", &HypreBasedAMS::GetHierarchyRefreshes,
          "Updates that refreshed the frozen AMG hierarchies instead of rebuilding them.")
      .def_property_readonly("in_place_updates", &HypreBasedAMS::GetInPlaceUpdates,
          "Updates (reuse_hierarchy, unchanged sparsity) whose Galerkin products ran "
          "numerically on the previous patterns, without symbolic products or allocation.")
      .def_property_readonly("setup_workers", &HypreBasedAMS::GetSetupWorkers,
          "Observed workers in AMG strength-row construction (zero when no coarsening occurs).")
      .def("Update", py::overload_cast<>(&HypreBasedAMS::Update),
           "Rebuild preconditioner with current matrix values (geometry preserved).")
      .def("Update", py::overload_cast<shared_ptr<SparseMatrix<double>>>(&HypreBasedAMS::Update),
           py::arg("new_mat"),
           "Update with a new system matrix, then rebuild.");

  py::class_<ComplexHypreBasedAMS, shared_ptr<ComplexHypreBasedAMS>, BaseMatrix>
      (m, "ComplexHypreBasedAMSPreconditionerImpl")
      .def("Update", py::overload_cast<>(&ComplexHypreBasedAMS::Update),
           "Rebuild preconditioner with current matrix values (geometry preserved).")
      .def("Update", py::overload_cast<shared_ptr<SparseMatrix<double>>>(&ComplexHypreBasedAMS::Update),
           py::arg("new_a_real"),
           "Update with a new real auxiliary matrix, then rebuild.");

  m.def("CompactAMGPreconditioner",
    [](shared_ptr<BaseMatrix> mat,
       py::object freedofs_obj,
       double theta,
       int max_levels,
       int min_coarse,
       int num_smooth,
       int print_level) -> shared_ptr<CompactAMG>
    {
      auto sp_mat = dynamic_pointer_cast<SparseMatrix<double>>(mat);
      if (!sp_mat)
        throw py::type_error("CompactAMGPreconditioner: expected real SparseMatrix<double>");
      auto sp_freedofs = ExtractFreeDofs(freedofs_obj);

      auto amg = make_shared<CompactAMG>(sp_mat, sp_freedofs, theta,
                                          max_levels, min_coarse, num_smooth,
                                          print_level);
      amg->Setup();
      return amg;
    },
    py::arg("mat"),
    py::arg("freedofs") = py::none(),
    py::arg("theta") = 0.25,
    py::arg("max_levels") = 25,
    py::arg("min_coarse") = 50,
    py::arg("num_smooth") = 1,
    py::arg("print_level") = 0,
    R"raw_string(
Compact Algebraic Multigrid (AMG) Preconditioner.

TaskManager-parallel AMG for scalar H1 problems. No external dependency.
Uses PMIS coarsening + classical interpolation + l1-Jacobi smoother.

Parameters:

mat : SparseMatrix (real)
  H1 system matrix.
freedofs : BitArray, optional
  Free DOFs mask.
theta : float
  Strength threshold (default=0.25 for 3D).
max_levels : int
  Maximum AMG levels (default=25).
min_coarse : int
  Minimum DOFs for direct solve (default=50).
num_smooth : int
  Smoother sweeps per level (default=1).
print_level : int
  Verbosity (0=silent, default=0).
)raw_string");

  m.def("TaskManagerActive", []() { return ngcore::GetTaskManager() != nullptr; },
    "True inside an ngsolve.TaskManager context: AMS construction and Update refuse to run there.");

  m.def("LowestOrderGradient",
    [](shared_ptr<ngcomp::FESpace> fes) -> shared_ptr<BaseMatrix>
    {
      // Discrete gradient H1(order 1) -> lowest-order HCurl straight from the
      // mesh edge table: row e holds -1 / +1 at the edge's first / second
      // vertex, exactly as FESpace.CreateGradient() for this space.
      auto ma = fes->GetMeshAccess();
      const size_t nedge = ma->GetNEdges();
      const size_t nv = ma->GetNV();
      if (fes->GetNDof() != nedge)
        throw py::value_error("LowestOrderGradient: needs one HCurl dof per edge "
                              "(HCurl(order=1, nograds=True) or order=0)");
      std::atomic<size_t> mismatched{0};
      ParallelFor(nedge, [&](size_t e) {
        Array<ngcomp::DofId> dnums;
        fes->GetDofNrs(ngfem::NodeId(ngfem::NT_EDGE, e), dnums);
        if (dnums.Size() != 1 || size_t(dnums[0]) != e) mismatched++;
      });
      if (mismatched)
        throw py::value_error("LowestOrderGradient: HCurl dofs are not numbered by edge ("
                              + std::to_string(mismatched.load()) + " edges differ)");
      (void)nv;
      return shared_ptr<BaseMatrix>(BuildLowestOrderGradient(*ma));
    },
    py::arg("fes"),
    R"raw_string(
Discrete gradient of the vertex space into a lowest-order HCurl space.

Built directly from the mesh edge table in parallel; the matrix equals the
first return value of fes.CreateGradient() for HCurl(order=1, nograds=True)
(or order=0). Fails for any space whose dofs are not one per edge in edge
order.
)raw_string");

  m.def("LowestOrderCurlSystem",
    [](shared_ptr<ngcomp::FESpace> fes, py::object coefficient) -> py::dict
    {
      // One pass over the volume elements of a lowest-order HCurl space on
      // straight tetrahedra: element dofs, the element-constant curl of the
      // space's own basis, the volume, the element-graph matrix and every
      // element's 36 value positions in it.
      auto ma = fes->GetMeshAccess();
      if (ma->GetDimension() != 3)
        throw py::value_error("LowestOrderCurlSystem: needs a three-dimensional mesh");
      const size_t ne = ma->GetNE(ngfem::VOL);
      const size_t ndof = fes->GetNDof();
      py::array_t<int32_t> dofs({py::ssize_t(ne), py::ssize_t(6)});
      py::array_t<double> curl({py::ssize_t(ne), py::ssize_t(3), py::ssize_t(6)});
      py::array_t<double> volume{std::vector<py::ssize_t>{py::ssize_t(ne)}};
      int32_t* dof_out = dofs.mutable_data();
      double* curl_out = curl.mutable_data();
      double* volume_out = volume.mutable_data();
      std::atomic<size_t> rejected{0};
      ParallelForRange(ne, [&](IntRange range) {
        LocalHeap lh(100000, "LowestOrderCurlSystem", true);
        Array<ngcomp::DofId> dnums;
        for (auto i : range) {
          HeapReset hr(lh);
          ngfem::ElementId ei(ngfem::VOL, i);
          if (ma->GetElType(ei) != ngfem::ET_TET) { rejected++; continue; }
          fes->GetDofNrs(ei, dnums);
          const auto& fel = fes->GetFE(ei, lh);
          auto hcurl = dynamic_cast<const ngfem::HCurlFiniteElement<3>*>(&fel);
          if (!hcurl || fel.GetNDof() != 6 || dnums.Size() != 6) { rejected++; continue; }
          const auto& trafo = ma->GetTrafo(ei, lh);
          if (trafo.IsCurvedElement()) { rejected++; continue; }
          ngfem::IntegrationPoint ip(0.25, 0.25, 0.25);
          ngfem::MappedIntegrationPoint<3, 3> mip(ip, trafo);
          FlatMatrix<double> shape(6, 3, lh);
          hcurl->CalcMappedCurlShape(mip, shape);
          // Sort the six dofs ascending and carry the curl columns along.
          int order[6] = {0, 1, 2, 3, 4, 5};
          std::sort(order, order + 6, [&](int a, int b) { return dnums[a] < dnums[b]; });
          for (int a = 0; a < 6; a++) {
            if (!ngcomp::IsRegularDof(dnums[order[a]])) { rejected++; break; }
            dof_out[6 * i + a] = int32_t(dnums[order[a]]);
            for (int k = 0; k < 3; k++)
              curl_out[18 * i + 6 * k + a] = shape(order[a], k);
          }
          volume_out[i] = std::abs(mip.GetJacobiDet()) / 6.0;
        }
      });
      if (rejected)
        throw py::value_error("LowestOrderCurlSystem: needs straight tetrahedra with six regular "
                              "HCurl dofs each (HCurl(order=1, nograds=True)); "
                              + std::to_string(rejected.load()) + " elements differ");
      Array<int> sizes(ne);
      sizes = 6;
      Table<int> elements(sizes);
      ParallelFor(ne, [&](size_t i) {
        for (int a = 0; a < 6; a++) elements[i][a] = dof_out[6 * i + a];
      });
      auto matrix = make_shared<SparseMatrix<double>>(ndof, ndof, elements, elements, false);
      matrix->AsVector() = 0.0;
      py::array_t<int32_t> positions({py::ssize_t(ne), py::ssize_t(36)});
      int32_t* position_out = positions.mutable_data();
      if (matrix->NZE() > size_t(std::numeric_limits<int32_t>::max()))
        throw py::value_error("LowestOrderCurlSystem: matrix too large for int32 positions");
      std::atomic<size_t> missing{0};
      ParallelFor(ne, [&](size_t i) {
        const int32_t* d = dof_out + 6 * i;
        for (int a = 0; a < 6; a++)
          for (int b = 0; b < 6; b++) {
            const size_t position = matrix->GetPositionTest(d[a], d[b]);
            if (position == numeric_limits<size_t>::max()) { missing++; continue; }
            position_out[36 * i + 6 * a + b] = int32_t(position);
          }
      });
      if (missing)
        throw std::runtime_error("LowestOrderCurlSystem: element entries outside the graph");
      if (!coefficient.is_none()) {
        // matrix = sum_e c_e vol_e C_e^T C_e, gathered row by row over the
        // elements of each dof in ascending element order: deterministic, and
        // entries (i, j) and (j, i) sum identical terms in the same order, so
        // the matrix is exactly symmetric.
        auto c = coefficient.cast<py::array_t<double, py::array::c_style | py::array::forcecast>>();
        if (c.ndim() != 1 || size_t(c.shape(0)) != ne)
          throw py::value_error("LowestOrderCurlSystem: coefficient needs one value per element");
        const double* cf = c.data();
        TableCreator<int> creator(ndof);
        for (; !creator.Done(); creator++)
          for (size_t i = 0; i < ne; i++)
            for (int a = 0; a < 6; a++)
              creator.Add(dof_out[6 * i + a], int(6 * i + a));
        Table<int> dof_elements = creator.MoveTable();
        auto values = matrix->AsVector().FVDouble();
        ParallelFor(ndof, [&](size_t row) {
          for (int slot : dof_elements[row]) {
            const size_t i = size_t(slot) / 6;
            const int a = slot % 6;
            const double scale = cf[i] * volume_out[i];
            if (scale == 0.0) continue;
            const double* C = curl_out + 18 * i;
            for (int b = 0; b < 6; b++)
              values[position_out[36 * i + 6 * a + b]] +=
                  scale * (C[a] * C[b] + C[6 + a] * C[6 + b] + C[12 + a] * C[12 + b]);
          }
        });
      }
      py::dict result;
      result["dofs"] = dofs;
      result["curl"] = curl;
      result["volume"] = volume;
      result["matrix"] = shared_ptr<BaseMatrix>(matrix);
      result["positions"] = positions;
      return result;
    },
    py::arg("fes"), py::arg("coefficient") = py::none(),
    R"raw_string(
Element data of a lowest-order HCurl space on straight tetrahedra, in one pass.

Returns a dict: ``dofs`` (ne, 6) int32 ascending per element; ``curl``
(ne, 3, 6), the element-constant curl of the space's basis functions in that
dof order; ``volume`` (ne,); ``matrix``, a zero SparseMatrix with the element
graph (full storage, sorted columns), i.e. the pattern of any bilinear form on
the space; ``positions`` (ne, 36) int32, the index into ``matrix`` values of
entry (dofs[e, a], dofs[e, b]) at 6 a + b. With `coefficient` (ne,), the
matrix holds sum_e coefficient[e] vol_e curl_e^T curl_e (the curl-curl form
with an element-constant coefficient), summed deterministically and exactly
symmetric; otherwise its values are zero. Parallel under an ngsolve
TaskManager. Fails for curved, non-tetrahedral or higher-order elements.
)raw_string");

  py::class_<LowestOrderCurlJacobian>(m, "LowestOrderCurlJacobian", R"raw_string(
Newton Jacobian refresh for a lowest-order HCurl curl-curl matrix.

Constructed from the ``matrix`` of LowestOrderCurlSystem (holding the constant
part) and the ``dofs``, ``curl``, ``volume`` and ``positions`` rows of the
nonlinear elements. ``Refresh(nu, q, b)`` rewrites every row those elements
touch as its constant part plus sum_e vol_e C_e^T (nu_e I + q_e b_e b_e^T) C_e,
gathered per row in ascending element order (deterministic, exactly symmetric).
Parallel under an ngsolve TaskManager.
)raw_string")
    .def(py::init<shared_ptr<BaseMatrix>,
                  py::array_t<int32_t, py::array::c_style | py::array::forcecast>,
                  py::array_t<double, py::array::c_style | py::array::forcecast>,
                  py::array_t<double, py::array::c_style | py::array::forcecast>,
                  py::array_t<int32_t, py::array::c_style | py::array::forcecast>>(),
         py::arg("matrix"), py::arg("dofs"), py::arg("curl"), py::arg("volume"), py::arg("positions"))
    .def("Refresh", &LowestOrderCurlJacobian::Refresh, py::arg("nu"), py::arg("q"), py::arg("b"),
         "Rewrite the touched rows for element reluctivity nu, rank-one q and flux b (n, 3).")
    .def_property_readonly("elements", &LowestOrderCurlJacobian::NumElements)
    .def_property_readonly("rows", &LowestOrderCurlJacobian::NumRows);

  py::class_<NativePCG>(m, "NativePCG", R"raw_string(
Preconditioned CG stopped on the true relative residual over free dofs.

``NativePCG(mat, pre, freedofs)``; ``Solve(b, x, tolerance, maxiter)`` starts
from x = 0 and returns ``(iterations, true_relative_residual, converged)``.
Each iteration is one matrix product and one preconditioner application; the
true residual is computed only to confirm convergence. Deterministic inner
products. Parallel under an ngsolve TaskManager.
)raw_string")
    .def(py::init<shared_ptr<BaseMatrix>, shared_ptr<BaseMatrix>, shared_ptr<BitArray>>(),
         py::arg("mat"), py::arg("pre"), py::arg("freedofs"))
    .def("Solve", &NativePCG::Solve, py::arg("b"), py::arg("x"), py::arg("tolerance"),
         py::arg("maxiter"));

  m.def("ClosedCoilCurrentPhi",
    [](shared_ptr<ngcomp::MeshAccess> mesh, std::vector<int> materials, double current_A,
       std::array<double, 3> origin, std::array<double, 3> normal, double radius,
       std::string inverse) {
      if (materials.empty()) throw py::value_error("conductor material labels are missing");
      return ClosedCoilCurrentPhiImpl(mesh, materials, current_A, origin, normal, radius, inverse);
    },
    py::arg("mesh"), py::arg("materials"), py::arg("current_A"), py::arg("origin"),
    py::arg("normal"), py::arg("radius"), py::arg("inverse") = "sparsecholesky",
    R"raw_string(
A-phi DC current of a closed conductor with one thick cut, natively.

``materials`` are 0-based material indices of the conductor; ``normal`` must
be a unit vector. Returns ``elements`` (conductor element numbers),
``density`` (n, 3) in A/m^2, ``vertices`` and ``phi`` (conductor vertices and
the potential), ``cut_faces``, ``support_cells``, ``phi_ndof``,
``relative_weak_divergence`` and ``cut_face_flux_A``. Raises the same errors as
radia.meshed_current.solve_closed_coil_current_phi. Parallel under a
TaskManager.
)raw_string");

  py::class_<LowestOrderCurlResidual>(m, "LowestOrderCurlResidual", R"raw_string(
Nonlinear state and residual of a lowest-order HCurl curl-curl problem.

Constructed from the ``dofs``, ``curl`` and ``volume`` of LowestOrderCurlSystem
(all elements), ``ndof``, the nonlinear element numbers ``iron`` and ``nu0``.
``Evaluate(x, source_mean, grid, nu, dhdb, load, residual)`` returns
``(b, magnitude, nu, q)`` on the iron (b = element curl + source_mean, which
may be None; nu and dH/dB interpolated like numpy.interp on the ascending
grid; q = (dH/dB - nu)/|b|^2, 0 below 1e-9) and writes ``residual`` =
sum_e vol_e C_e^T (nu_e c_e + [iron] (nu_e - nu0) Bs_e) - load, with nu0 off
the iron, gathered per dof (deterministic). Parallel under a TaskManager.
)raw_string")
    .def(py::init<py::array_t<int32_t, py::array::c_style | py::array::forcecast>,
                  py::array_t<double, py::array::c_style | py::array::forcecast>,
                  py::array_t<double, py::array::c_style | py::array::forcecast>, size_t,
                  py::array_t<int64_t, py::array::c_style | py::array::forcecast>, double>(),
         py::arg("dofs"), py::arg("curl"), py::arg("volume"), py::arg("ndof"), py::arg("iron"),
         py::arg("nu0"))
    .def("Evaluate", &LowestOrderCurlResidual::Evaluate, py::arg("x"), py::arg("source_mean"),
         py::arg("grid"), py::arg("nu"), py::arg("dhdb"), py::arg("load"), py::arg("residual"));

  m.def("HypreBasedAMSPreconditioner",
    [](shared_ptr<BaseMatrix> mat,
       shared_ptr<BaseMatrix> grad_mat,
       py::object freedofs_obj,
       py::list coord_x_list,
       py::list coord_y_list,
       py::list coord_z_list,
       int cycle_type,
       int print_level,
       int subspace_solver,
       int num_smooth,
       bool beta_zero,
       bool reuse_hierarchy,
       bool mixed_precision) -> shared_ptr<HypreBasedAMS>
    {
      auto sp_mat = dynamic_pointer_cast<SparseMatrix<double>>(mat);
      if (!sp_mat)
        throw py::type_error("HypreBasedAMSPreconditioner: expected real SparseMatrix<double>");

      auto sp_grad = dynamic_pointer_cast<SparseMatrix<double>>(grad_mat);
      if (!sp_grad)
        throw py::type_error("HypreBasedAMSPreconditioner: grad_mat must be real SparseMatrix<double>");

      auto sp_freedofs = ExtractFreeDofs(freedofs_obj);

      auto to_vec = [](py::list lst) {
        std::vector<double> v(lst.size());
        for (size_t i = 0; i < lst.size(); i++)
          v[i] = lst[i].cast<double>();
        return v;
      };

      return make_shared<HypreBasedAMS>(
          sp_mat, sp_grad, sp_freedofs,
          to_vec(coord_x_list), to_vec(coord_y_list), to_vec(coord_z_list),
          cycle_type, num_smooth, 0.25, print_level, 1.0, subspace_solver, beta_zero,
          reuse_hierarchy, mixed_precision);
    },
    py::arg("mat"),
    py::arg("grad_mat"),
    py::arg("freedofs") = py::none(),
    py::arg("coord_x") = py::list(),
    py::arg("coord_y") = py::list(),
    py::arg("coord_z") = py::list(),
    py::arg("cycle_type") = 1,
    py::arg("print_level") = 0,
    py::arg("subspace_solver") = 0,
    py::arg("num_smooth") = 1,
    py::arg("beta_zero") = false,
    py::arg("reuse_hierarchy") = false,
    py::arg("mixed_precision") = false,
    R"raw_string(
Compact AMS (Auxiliary-space Maxwell Solver) Preconditioner.

TaskManager-parallel AMS for HCurl curl-curl + mass systems. No external dependency.
Uses CompactAMG as sub-solver for gradient and nodal auxiliary spaces.
With beta_zero=True, omit the gradient correction and hierarchy, as in the
two-level beta=0 AMS mode. Use for compatible pure curl-curl systems; this
does not add a mass penalty, shift the matrix, or project the right-hand side.

Construct and call Update outside ngsolve.TaskManager. An active context raises
RuntimeError; applying an already-built preconditioner may run inside TaskManager.

Requires a lowest-order HCurl space (order=1, nograds=True); any other space
raises RuntimeError. For order >= 2 use NGSolve's bddc preconditioner.

Parameters:

mat : SparseMatrix (real)
  HCurl system matrix (must be nonsymmetric storage).
grad_mat : SparseMatrix (real)
  Discrete gradient matrix (HCurl -> H1).
freedofs : BitArray, optional
  Free DOFs mask.
coord_x, coord_y, coord_z : list of float
  Vertex coordinates (length = number of H1 DOFs).
cycle_type : int
  AMS cycle type (1=01210, 7=0201020, default=1).
print_level : int
  Verbosity (0=silent, default=0).
beta_zero : bool
  Pure curl-curl mode without the gradient subspace solver (default=False).
reuse_hierarchy : bool
  Update() keeps the AMG coarsening and interpolation of the first build and
  refreshes only the Galerkin coarse matrices (frozen interpolation), for a
  sequence of matrices on one sparsity pattern, e.g. Newton (default=False).
mixed_precision : bool
  The residual products inside the cycle (fine level and every AMG level) read
  float32 copies of the matrix values; vectors and sums stay double. The
  preconditioner stays a fixed symmetric operator, so CG remains valid; use it
  with a Krylov method that checks the double-precision residual. Requires
  beta_zero: with a gauge the float32 values lose the gauge-scale gradient
  components the gradient correction amplifies (default=False).
)raw_string");

  m.def("ComplexHypreBasedAMSPreconditioner",
    [](shared_ptr<BaseMatrix> a_real_mat,
       shared_ptr<BaseMatrix> grad_mat,
       py::object freedofs_obj,
       py::list coord_x_list,
       py::list coord_y_list,
       py::list coord_z_list,
       int ndof_complex,
       int cycle_type,
       int print_level,
       double correction_weight,
       int subspace_solver,
       int num_smooth) -> shared_ptr<ComplexHypreBasedAMS>
    {
      auto sp_mat = dynamic_pointer_cast<SparseMatrix<double>>(a_real_mat);
      if (!sp_mat)
        throw py::type_error("ComplexHypreBasedAMSPreconditioner: a_real_mat must be real SparseMatrix<double>");

      auto sp_grad = dynamic_pointer_cast<SparseMatrix<double>>(grad_mat);
      if (!sp_grad)
        throw py::type_error("ComplexHypreBasedAMSPreconditioner: grad_mat must be real SparseMatrix<double>");

      auto sp_freedofs = ExtractFreeDofs(freedofs_obj);

      auto to_vec = [](py::list lst) {
        std::vector<double> v(lst.size());
        for (size_t i = 0; i < lst.size(); i++)
          v[i] = lst[i].cast<double>();
        return v;
      };

      // Auto-derive ndof_complex from matrix if not specified (0 = auto)
      int ndof = ndof_complex;
      if (ndof <= 0) {
        ndof = static_cast<int>(sp_mat->VHeight());
      }

      return make_shared<ComplexHypreBasedAMS>(
          sp_mat, sp_grad, sp_freedofs,
          to_vec(coord_x_list), to_vec(coord_y_list), to_vec(coord_z_list),
          ndof, cycle_type, print_level, correction_weight,
          subspace_solver, num_smooth);
    },
    py::arg("a_real_mat"),
    py::arg("grad_mat"),
    py::arg("freedofs") = py::none(),
    py::arg("coord_x") = py::list(),
    py::arg("coord_y") = py::list(),
    py::arg("coord_z") = py::list(),
    py::arg("ndof_complex") = 0,
    py::arg("cycle_type") = 1,
    py::arg("print_level") = 0,
    py::arg("correction_weight") = 1.0,
    py::arg("subspace_solver") = 0,
    py::arg("num_smooth") = 1,
    R"raw_string(
Complex Compact AMS preconditioner with TaskManager Re/Im parallelism.

For complex eddy current problems (A = K + jw*sigma*M). Creates TWO
independent HypreBasedAMS solver instances and applies them to the real
and imaginary parts in parallel via NGSolve TaskManager.

No external dependency (pure C++ header-only).
Use with COCRSolver (complex symmetric) or GMRESSolver (general).

Construct and call Update outside ngsolve.TaskManager. An active context raises
RuntimeError; applying an already-built preconditioner may run inside TaskManager.

Requires a lowest-order HCurl space (order=1, nograds=True); any other space
raises RuntimeError. For order >= 2 use NGSolve's bddc preconditioner.

Parameters:

a_real_mat : SparseMatrix (real)
  Real auxiliary matrix (K + eps*M + |omega|*sigma*M_cond).
grad_mat : SparseMatrix (real)
  Discrete gradient matrix (HCurl -> H1).
freedofs : BitArray, optional
  Free DOFs mask for HCurl space.
coord_x, coord_y, coord_z : list of float
  Vertex coordinates (length = number of H1 DOFs).
ndof_complex : int
  Number of complex DOFs (= fes.ndof for complex HCurl space).
cycle_type : int
  AMS cycle type (1=01210, 7=0201020, default=1).
print_level : int
  Verbosity (0=silent, default=0).
)raw_string");

  RegisterAMSWirebasket();
  m.def("AMSCoarseStats", []() {
      const auto& s = GetAMSWirebasketStats();
      py::dict d;
      d["builds"] = s.builds;
      d["n_edges"] = s.n_edges;
      d["n_extra"] = s.n_extra;
      d["cycles"] = s.cycles;
      d["n_free"] = s.n_free;
      d["n_vertices"] = s.n_vertices;
      d["complex"] = s.complex;
      d["extract_s"] = s.extract_s;
      d["setup_s"] = s.setup_s;
      d["applies"] = s.applies;
      d["apply_s"] = s.apply_s;
      if (auto ams = s.complex_ams.lock()) {
        const char* names[10] = {"split", "grad_residual", "grad_restrict", "grad_amg",
                                 "grad_prolong", "nodal_residual", "nodal_restrict",
                                 "nodal_amg", "nodal_prolong", "final_sweep"};
        auto sec = ams->StageSeconds();
        py::dict stages;
        for (int k = 0; k < 10; k++) stages[names[k]] = sec[k];
        d["cycle_stage_s"] = stages;
        d["cycles_run"] = ams->CyclesRun();
        py::dict amg_levels;
        const auto& real = ams->Real();
        const std::pair<const char*, CompactAMG*> amgs[4] = {
            {"G", real.GetBGAsAMG()}, {"Px", real.GetBPixAsAMG()},
            {"Py", real.GetBPiyAsAMG()}, {"Pz", real.GetBPizAsAMG()}};
        for (const auto& [name, amg] : amgs) {
          if (!amg) continue;
          py::list rows;
          for (const auto& row : amg->DualLevelProfile())
            rows.append(py::make_tuple(int(row[0]), int(row[1]), row[2], int(row[3])));
          amg_levels[name] = rows;
        }
        d["amg_levels"] = amg_levels;
      }
      return d;
    },
    R"raw_string(
Counters of the most recent BDDC wirebasket AMS (coarsetype="sparsesolv_ams").

Importing this module registers "sparsesolv_ams" with NGSolve's preconditioner
classes. For an HCurl BilinearForm, Preconditioner(a, "bddc",
coarsetype="sparsesolv_ams", coarseflags={...}) replaces the direct wirebasket
inverse by Compact AMS on the lowest-order edge block. coarseflags: cycles
(k AMS cycles as k stationary steps on the wirebasket system, default 1),
lean_coarse (default 1: auxiliary AMGs stop where coarsening stalls and solve
a coarsest level of at most 1024 rows densely when that reproduces a test
vector; 0 keeps sparse Cholesky),
cycle_type, num_smooth, print_level, eps (relative diagonal shift of the AMS
surrogate), beta_zero and mixed_precision (real systems only). A complex
wirebasket matrix S uses the real surrogate Re S + Im S.
)raw_string");

  m.def("has_compact_ams", []() { return true; },
    "Returns True if Compact AMG/AMS support is available.");

  // Back-compat aliases (renamed 2026-06-27: CompactAMS -> HypreBasedAMS, since the algorithm IS HYPRE's
  // AMS / Kolev-Vassilevski, just HYPRE-free; "Compact" misleadingly implied a distinct reduced variant).
  // Old names keep working so panels / MCP recipes / external code do not break.
  m.attr("CompactAMSPreconditioner") = m.attr("HypreBasedAMSPreconditioner");
  m.attr("ComplexCompactAMSPreconditioner") = m.attr("ComplexHypreBasedAMSPreconditioner");
}

// ============================================================================
// COCR solver (NGSolve BaseMatrix interface, accepts external preconditioner)
// ============================================================================

inline void ExportCOCRSolver(py::module& m) {
  // Register C++ types so .iterations property works
  py::class_<COCRSolverNGS<double>, shared_ptr<COCRSolverNGS<double>>, BaseMatrix>
      (m, "COCRSolverD")
      .def_property_readonly("iterations", &COCRSolverNGS<double>::GetIterations);

  py::class_<COCRSolverNGS<Complex>, shared_ptr<COCRSolverNGS<Complex>>, BaseMatrix>
      (m, "COCRSolverC")
      .def_property_readonly("iterations", &COCRSolverNGS<Complex>::GetIterations);

  m.def("COCRSolver",
    [](shared_ptr<BaseMatrix> mat,
       shared_ptr<BaseMatrix> pre,
       py::object freedofs_obj,
       int maxiter,
       double tol,
       bool printrates) -> shared_ptr<BaseMatrix>
    {
      auto sp_freedofs = ExtractFreeDofs(freedofs_obj);
      if (mat->IsComplex()) {
        return make_shared<COCRSolverNGS<Complex>>(
            mat, pre, sp_freedofs, maxiter, tol, printrates);
      } else {
        return make_shared<COCRSolverNGS<double>>(
            mat, pre, sp_freedofs, maxiter, tol, printrates);
      }
    },
    py::arg("mat"),
    py::arg("pre"),
    py::arg("freedofs") = py::none(),
    py::arg("maxiter") = 500,
    py::arg("tol") = 1e-8,
    py::arg("printrates") = false,
    R"raw_string(
COCR (Conjugate Orthogonal Conjugate Residual) solver for complex-symmetric systems.

For A^T = A (NOT Hermitian). Uses unconjugated inner products (x^T y).
Minimizes ||A r~||_2 for smoother convergence than COCG/CG.

When to use COCRSolver vs SparseSolvSolver(method="COCR"):
  - COCRSolver(mat, pre): accepts any external BaseMatrix preconditioner
    (e.g., IC, Compact AMS). Same interface as NGSolve CGSolver.
  - SparseSolvSolver(method="COCR"): uses internal IC preconditioner with
    auto-shift, ABMC ordering, divergence detection. Unified solver interface.

Usage (same as NGSolve CGSolver):
  inv = COCRSolver(mat, pre, maxiter=500, tol=1e-8)
  gfu.vec.data = inv * rhs.vec

For COCG, use CGSolver(mat, pre, conjugate=False) instead.

Parameters:

mat : BaseMatrix
  System matrix (real or complex).
pre : BaseMatrix
  Preconditioner (must be symmetric for COCR).
maxiter : int
  Maximum iterations (default: 500).
tol : float
  Relative convergence tolerance (default: 1e-8).
printrates : bool
  Print convergence info (default: False).

Reference: Sogabe & Zhang (2007), J. Comput. Appl. Math., 199(2), 297-303.
)raw_string");
}

// ============================================================================
// Public API: Single entry point for NGSolve integration
// ============================================================================

// ============================================================================
// GMRES solver (NGSolve BaseMatrix interface, accepts external preconditioner)
// ============================================================================

inline void ExportGMRESSolver(py::module& m) {
  py::class_<GMRESSolverNGS<double>, shared_ptr<GMRESSolverNGS<double>>, BaseMatrix>
      (m, "GMRESSolverD")
      .def_property_readonly("iterations", &GMRESSolverNGS<double>::GetIterations);

  py::class_<GMRESSolverNGS<Complex>, shared_ptr<GMRESSolverNGS<Complex>>, BaseMatrix>
      (m, "GMRESSolverC")
      .def_property_readonly("iterations", &GMRESSolverNGS<Complex>::GetIterations);

  m.def("GMRESSolver",
    [](shared_ptr<BaseMatrix> mat,
       shared_ptr<BaseMatrix> pre,
       py::object freedofs_obj,
       int maxiter,
       double tol,
       int restart,
       bool printrates) -> shared_ptr<BaseMatrix>
    {
      auto freedofs = ExtractFreeDofs(freedofs_obj);
      if (mat->IsComplex()) {
        return make_shared<GMRESSolverNGS<Complex>>(
            mat, pre, freedofs, maxiter, tol, restart, printrates);
      } else {
        return make_shared<GMRESSolverNGS<double>>(
            mat, pre, freedofs, maxiter, tol, restart, printrates);
      }
    },
    py::arg("mat"),
    py::arg("pre"),
    py::arg("freedofs") = py::none(),
    py::arg("maxiter") = 500,
    py::arg("tol") = 1e-8,
    py::arg("restart") = 0,
    py::arg("printrates") = false,
    R"raw_string(
Left-preconditioned GMRES solver for non-symmetric linear systems.

1 SpMV + 1 preconditioner application per iteration.
Optimal for AMS preconditioned eddy current problems.

Parameters:

mat : BaseMatrix
  System matrix (real or complex, auto-detected).
pre : BaseMatrix
  Preconditioner (any BaseMatrix, need not be symmetric).
freedofs : BitArray, optional
  Free DOFs. Constrained DOFs are zeroed out during iteration.
maxiter : int
  Maximum iterations (default: 500).
tol : float
  Relative convergence tolerance (default: 1e-8).
restart : int
  Restart after this many iterations (0 = no restart, default: 0).
printrates : bool
  Print convergence info (default: False).
)raw_string");
}

/// Register all SparseSolv Python bindings (type registration + factory functions)
inline void ExportSparseSolvBindings(py::module& m) {
  ExportSparseSolvResult_impl(m);
  ExportSparseSolvTyped<double>(m, "D");
  ExportSparseSolvTyped<Complex>(m, "C");
  ExportSparseSolvFactories(m);
  ExportHypreBasedAMS(m);
  ExportCOCRSolver(m);
  ExportGMRESSolver(m);
}

} // namespace ngla

#endif // NGSOLVE_SPARSESOLV_PYTHON_EXPORT_HPP
