/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file python_module.cpp
 * @brief pybind11 module for SparseSolv NGSolve integration
 *
 * Built into the radia wheel as `radia/sparsesolv_ngsolve.pyd` via the
 * top-level Radia CMake target `add_ngsolve_python_module(sparsesolv_ngsolve ...)`.
 * Usage: `import radia.sparsesolv_ngsolve as ssn`
 *
 * (`import radia` first to register the NGSolve / MKL DLL search paths
 * on Windows; otherwise the .pyd will fail to load.)
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "sparsesolv/ngsolve/sparsesolv_python_export.hpp"

namespace py = pybind11;

PYBIND11_MODULE(sparsesolv_ngsolve, m) {
    // Import NGSolve modules so pybind11 knows about BaseMatrix, BilinearForm, etc.
    py::module_::import("ngsolve.la");
    py::module_::import("ngsolve.comp");

    m.doc() = "SparseSolv iterative solvers and preconditioners for NGSolve\n\n"
              "Provides IC preconditioner, HypreBasedAMS/ComplexHypreBasedAMS (HCurl AMS, HYPRE-free),\n"
              "CompactAMG, and ICCG/CG/COCR/GMRES iterative solvers\n"
              "for use with NGSolve's sparse linear algebra.\n\n"
              "Based on JP-MARs/SparseSolv (https://github.com/JP-MARs/SparseSolv)";

    ngla::ExportSparseSolvBindings(m);
}
