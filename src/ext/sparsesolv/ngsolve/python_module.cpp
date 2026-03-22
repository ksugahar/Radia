/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file python_module.cpp
 * @brief Standalone pybind11 module for SparseSolv NGSolve integration
 *
 * Builds as an independent .pyd/.so that links against NGSolve.
 * Usage: import sparsesolv_ngsolve
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
              "Provides IC/SGS preconditioners, CompactAMS/ComplexCompactAMS (HCurl),\n"
              "CompactAMG, and ICCG/SGSMRTR/CG/COCR/GMRES iterative solvers\n"
              "for use with NGSolve's sparse linear algebra.\n\n"
              "Based on JP-MARs/SparseSolv (https://github.com/JP-MARs/SparseSolv)";

    ngla::ExportSparseSolvBindings(m);
}
