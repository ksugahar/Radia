/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/// @file sparsesolv_python_export.hpp
/// @brief Pybind11 bindings: type registration + factory functions with auto-dispatch

#ifndef NGSOLVE_SPARSESOLV_PYTHON_EXPORT_HPP
#define NGSOLVE_SPARSESOLV_PYTHON_EXPORT_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "sparsesolv_precond.hpp"
#include "sparsesolv_solvers.hpp"
#include <comp.hpp>
#include <type_traits>

// Compact AMG/AMS (TaskManager-based, no external dependency)
#include "sparsesolv/preconditioners/compact_amg.hpp"
#include "sparsesolv/preconditioners/compact_ams.hpp"
#include "sparsesolv/preconditioners/complex_compact_ams.hpp"


namespace py = pybind11;

namespace ngla {

// ============================================================================
// Internal: SparseSolvResult (non-templated, called once)
// ============================================================================

inline void ExportSparseSolvResult_impl(py::module& m) {
  py::class_<SparseSolvResult>(m, "SparseSolvResult",
    "Result of a SparseSolv iterative solve (converged, iterations, final_residual, residual_history).")
    .def_readonly("converged", &SparseSolvResult::converged,
        "Whether the solver converged within tolerance")
    .def_readonly("iterations", &SparseSolvResult::iterations,
        "Number of iterations performed")
    .def_readonly("final_residual", &SparseSolvResult::final_residual,
        "Final relative residual (or best residual if save_best_result enabled)")
    .def_readonly("residual_history", &SparseSolvResult::residual_history,
        "Residual at each iteration (if save_residual_history enabled)")
    .def("__repr__", [](const SparseSolvResult& r) {
      return string("SparseSolvResult(converged=") +
             (r.converged ? "True" : "False") +
             ", iterations=" + std::to_string(r.iterations) +
             ", residual=" + std::to_string(r.final_residual) + ")";
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

  // SGS Preconditioner type registration
  {
    std::string cls_name = "SGSPreconditioner" + suffix;
    py::class_<SparseSolvSGSPreconditioner<SCAL>,
               shared_ptr<SparseSolvSGSPreconditioner<SCAL>>,
               BaseMatrix>(m, cls_name.c_str())
      .def("Update", &SparseSolvSGSPreconditioner<SCAL>::Update);
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
  SPD matrix (real or complex, auto-detected).
freedofs : BitArray, optional
  Free DOFs. Constrained DOFs treated as identity.
shift : float
  Shift parameter (default: 1.05).
)raw_string");

  // ---- SGSPreconditioner factory ----
  m.def("SGSPreconditioner", [](shared_ptr<BaseMatrix> mat,
                                  py::object freedofs) {
    auto sp_freedofs = ExtractFreeDofs(freedofs);
    shared_ptr<BaseMatrix> result;
    if (mat->IsComplex()) {
      auto sp = dynamic_pointer_cast<SparseMatrix<Complex>>(mat);
      if (!sp) throw py::type_error("SGSPreconditioner: expected SparseMatrix");
      auto p = make_shared<SparseSolvSGSPreconditioner<Complex>>(sp, sp_freedofs);
      p->Update();
      result = p;
    } else {
      auto sp = dynamic_pointer_cast<SparseMatrix<double>>(mat);
      if (!sp) throw py::type_error("SGSPreconditioner: expected SparseMatrix");
      auto p = make_shared<SparseSolvSGSPreconditioner<double>>(sp, sp_freedofs);
      p->Update();
      result = p;
    }
    return result;
  },
  py::arg("mat"), py::arg("freedofs") = py::none(),
  R"raw_string(
Symmetric Gauss-Seidel (SGS) Preconditioner.

Parameters:

mat : SparseMatrix
  SPD matrix (real or complex, auto-detected).
freedofs : BitArray, optional
  Free DOFs. Constrained DOFs treated as identity.
)raw_string");

  // ---- SparseSolvSolver factory ----
  m.def("SparseSolvSolver", [](shared_ptr<BaseMatrix> mat,
                                 const string& method, py::object freedofs,
                                 double tol, int maxiter, double shift,
                                 bool save_best_result, bool save_residual_history,
                                 bool printrates, bool conjugate,
                                 bool use_abmc, int abmc_block_size, int abmc_num_colors,
                                 bool abmc_reorder_spmv, bool abmc_use_rcm) {
    auto sp_freedofs = ExtractFreeDofs(freedofs);

    auto configure = [&](auto& solver) {
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
  py::arg("tol") = 1e-10,
  py::arg("maxiter") = 1000,
  py::arg("shift") = 1.05,
  py::arg("save_best_result") = true,
  py::arg("save_residual_history") = false,
  py::arg("printrates") = false,
  py::arg("conjugate") = false,
  py::arg("use_abmc") = false,
  py::arg("abmc_block_size") = 4,
  py::arg("abmc_num_colors") = 4,
  py::arg("abmc_reorder_spmv") = false,
  py::arg("abmc_use_rcm") = false,
  R"raw_string(
Iterative solver (ICCG / SGSMRTR / CG / COCR). Auto-detects real/complex.

Can be used as BaseMatrix (inverse operator) or via Solve() for detailed results.

Parameters:

mat : SparseMatrix
  System matrix (real or complex).
method : str
  "ICCG", "SGSMRTR", "CG", or "COCR".
  COCR: Conjugate Orthogonal Conjugate Residual for complex-symmetric A^T=A.
  For non-symmetric systems, use GMRESSolver instead.
freedofs : BitArray, optional
  Free DOFs.
tol : float
  Convergence tolerance (default: 1e-10).
maxiter : int
  Max iterations (default: 1000).
shift : float
  IC shift parameter (default: 1.05).
save_best_result : bool
  Track best solution (default: True).
save_residual_history : bool
  Record residual history (default: False).
printrates : bool
  Print convergence info (default: False).
conjugate : bool
  Conjugated inner product for Hermitian systems (default: False).

Properties (set after construction):
  auto_shift, diagonal_scaling, divergence_check, divergence_threshold,
  divergence_count, use_abmc, abmc_block_size, abmc_num_colors,
  abmc_reorder_spmv, abmc_use_rcm.
)raw_string");
}

// ============================================================================
// Compact AMG/AMS factory (TaskManager-based)
// ============================================================================

inline void ExportCompactAMS(py::module& m) {

  // Type registration for CompactAMG (needed for correct virtual dispatch in Python)
  py::class_<CompactAMG, shared_ptr<CompactAMG>, BaseMatrix>
      (m, "CompactAMGPreconditionerImpl");

  // Type registrations for CompactAMS / ComplexCompactAMS (enables Update() from Python)
  py::class_<CompactAMS, shared_ptr<CompactAMS>, BaseMatrix>
      (m, "CompactAMSPreconditionerImpl")
      .def("Update", py::overload_cast<>(&CompactAMS::Update),
           "Rebuild preconditioner with current matrix values (geometry preserved).")
      .def("Update", py::overload_cast<shared_ptr<SparseMatrix<double>>>(&CompactAMS::Update),
           py::arg("new_mat"),
           "Update with a new system matrix, then rebuild.");

  py::class_<ComplexCompactAMS, shared_ptr<ComplexCompactAMS>, BaseMatrix>
      (m, "ComplexCompactAMSPreconditionerImpl")
      .def("Update", py::overload_cast<>(&ComplexCompactAMS::Update),
           "Rebuild preconditioner with current matrix values (geometry preserved).")
      .def("Update", py::overload_cast<shared_ptr<SparseMatrix<double>>>(&ComplexCompactAMS::Update),
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

  m.def("CompactAMSPreconditioner",
    [](shared_ptr<BaseMatrix> mat,
       shared_ptr<BaseMatrix> grad_mat,
       py::object freedofs_obj,
       py::list coord_x_list,
       py::list coord_y_list,
       py::list coord_z_list,
       int cycle_type,
       int print_level,
       int subspace_solver,
       int num_smooth) -> shared_ptr<CompactAMS>
    {
      auto sp_mat = dynamic_pointer_cast<SparseMatrix<double>>(mat);
      if (!sp_mat)
        throw py::type_error("CompactAMSPreconditioner: expected real SparseMatrix<double>");

      auto sp_grad = dynamic_pointer_cast<SparseMatrix<double>>(grad_mat);
      if (!sp_grad)
        throw py::type_error("CompactAMSPreconditioner: grad_mat must be real SparseMatrix<double>");

      auto sp_freedofs = ExtractFreeDofs(freedofs_obj);

      auto to_vec = [](py::list lst) {
        std::vector<double> v(lst.size());
        for (size_t i = 0; i < lst.size(); i++)
          v[i] = lst[i].cast<double>();
        return v;
      };

      return make_shared<CompactAMS>(
          sp_mat, sp_grad, sp_freedofs,
          to_vec(coord_x_list), to_vec(coord_y_list), to_vec(coord_z_list),
          cycle_type, num_smooth, 0.25, print_level, 1.0, subspace_solver);
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
    R"raw_string(
Compact AMS (Auxiliary-space Maxwell Solver) Preconditioner.

TaskManager-parallel AMS for HCurl curl-curl + mass systems. No external dependency.
Uses CompactAMG as sub-solver for gradient and nodal auxiliary spaces.

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
)raw_string");

  m.def("ComplexCompactAMSPreconditioner",
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
       int num_smooth) -> shared_ptr<ComplexCompactAMS>
    {
      auto sp_mat = dynamic_pointer_cast<SparseMatrix<double>>(a_real_mat);
      if (!sp_mat)
        throw py::type_error("ComplexCompactAMSPreconditioner: a_real_mat must be real SparseMatrix<double>");

      auto sp_grad = dynamic_pointer_cast<SparseMatrix<double>>(grad_mat);
      if (!sp_grad)
        throw py::type_error("ComplexCompactAMSPreconditioner: grad_mat must be real SparseMatrix<double>");

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

      return make_shared<ComplexCompactAMS>(
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
independent CompactAMS solver instances and applies them to the real
and imaginary parts in parallel via NGSolve TaskManager.

No external dependency (pure C++ header-only).
Use with COCRSolver (complex symmetric) or GMRESSolver (general).

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

  m.def("has_compact_ams", []() { return true; },
    "Returns True if Compact AMG/AMS support is available.");
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
  ExportCompactAMS(m);
  ExportCOCRSolver(m);
  ExportGMRESSolver(m);
}

} // namespace ngla

#endif // NGSOLVE_SPARSESOLV_PYTHON_EXPORT_HPP
