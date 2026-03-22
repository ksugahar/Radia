/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/// @file sparsesolv.hpp
/// @brief Main header — includes all SparseSolv components

#ifndef SPARSESOLV_HPP
#define SPARSESOLV_HPP

// Core components
#include "core/types.hpp"
#include "core/constants.hpp"
#include "core/solver_config.hpp"
#include "core/sparse_matrix_view.hpp"
#include "core/preconditioner.hpp"
#include "core/abmc_ordering.hpp"

// Preconditioners
#include "preconditioners/ic_preconditioner.hpp"
#include "preconditioners/sgs_preconditioner.hpp"
// Solvers
#include "solvers/iterative_solver.hpp"
#include "solvers/cg_solver.hpp"
#include "solvers/cocr_solver.hpp"
#include "solvers/sgs_mrtr_solver.hpp"

namespace sparsesolv {

struct Version {
    static constexpr int major = 2;
    static constexpr int minor = 7;
    static constexpr int patch = 0;

    static std::string string() {
        return std::to_string(major) + "." +
               std::to_string(minor) + "." +
               std::to_string(patch);
    }
};

/// Convenience function: Solve Ax=b using ICCG
template<typename Scalar = double>
[[nodiscard]] inline SolverResult solve_iccg(
    const SparseMatrixView<Scalar>& A,
    const Scalar* b,
    Scalar* x,
    index_t size,
    const SolverConfig& config = SolverConfig()
) {
    // Create IC preconditioner with full config (for auto-shift, diagonal scaling, etc.)
    ICPreconditioner<Scalar> precond(config.shift_parameter);
    precond.set_config(config);
    precond.setup(A);

    // Create CG solver
    CGSolver<Scalar> solver;
    solver.set_config(config);

    if (precond.has_reordered_matrix()) {
        // Path 1: ABMC with abmc_reorder_spmv=true
        // CG runs entirely in ABMC-reordered space (legacy behavior)
        const auto& ord = precond.abmc_ordering();

        std::vector<Scalar> b_reord(size);
        for (index_t i = 0; i < size; ++i) b_reord[ord[i]] = b[i];

        std::vector<Scalar> x_reord(size);
        for (index_t i = 0; i < size; ++i) x_reord[ord[i]] = x[i];

        auto A_reord = precond.reordered_matrix_view();
        ICPrecondReorderedAdapter<Scalar> adapter(precond);
        auto result = solver.solve(A_reord, b_reord.data(), x_reord.data(),
                                   size, &adapter);

        for (index_t i = 0; i < size; ++i) x[i] = x_reord[ord[i]];
        return result;
    }

    if (precond.has_rcm_matrix()) {
        // Path 2: RCM+ABMC split mode
        // SpMV uses RCM-reordered matrix (better cache locality than original)
        // Preconditioner handles RCM->ABMC permutation internally
        const auto& rcm_ord = precond.rcm_perm();

        std::vector<Scalar> b_rcm(size);
        for (index_t i = 0; i < size; ++i) b_rcm[rcm_ord[i]] = b[i];

        std::vector<Scalar> x_rcm(size);
        for (index_t i = 0; i < size; ++i) x_rcm[rcm_ord[i]] = x[i];

        auto A_rcm = precond.rcm_matrix_view();
        ICPrecondRCMABMCAdapter<Scalar> rcm_adapter(precond);
        auto result = solver.solve(A_rcm, b_rcm.data(), x_rcm.data(),
                                   size, &rcm_adapter);

        for (index_t i = 0; i < size; ++i) x[i] = x_rcm[rcm_ord[i]];
        return result;
    }

    // Path 3: Standard (default)
    // SpMV uses original matrix A (preserves FEM mesh cache locality)
    // ABMC permutation handled internally by preconditioner::apply_abmc()
    return solver.solve(A, b, x, size, &precond);
}

/// Convenience function: Solve Ax=b using ICCG with std::vector
template<typename Scalar = double>
[[nodiscard]] inline SolverResult solve_iccg(
    const SparseMatrixView<Scalar>& A,
    const std::vector<Scalar>& b,
    std::vector<Scalar>& x,
    const SolverConfig& config = SolverConfig()
) {
    if (x.size() != b.size()) {
        x.resize(b.size());
    }
    return solve_iccg(A, b.data(), x.data(), static_cast<index_t>(b.size()),
                      config);
}

/// Convenience function: Solve Ax=b using SGS-MRTR
template<typename Scalar = double>
[[nodiscard]] inline SolverResult solve_sgsmrtr(
    const SparseMatrixView<Scalar>& A,
    const Scalar* b,
    Scalar* x,
    index_t size,
    const SolverConfig& config = SolverConfig()
) {
    // Use specialized SGS-MRTR solver with split formula
    SGSMRTRSolver<Scalar> solver;
    solver.set_config(config);

    // Solve
    return solver.solve(A, b, x, size);
}

/// Convenience function: Solve Ax=b using SGS-MRTR with std::vector
template<typename Scalar = double>
[[nodiscard]] inline SolverResult solve_sgsmrtr(
    const SparseMatrixView<Scalar>& A,
    const std::vector<Scalar>& b,
    std::vector<Scalar>& x,
    const SolverConfig& config = SolverConfig()
) {
    if (x.size() != b.size()) {
        x.resize(b.size());
    }
    return solve_sgsmrtr(A, b.data(), x.data(), static_cast<index_t>(b.size()), config);
}

/// Convenience function: Solve Ax=b using COCR (IC preconditioned)
template<typename Scalar = double>
[[nodiscard]] inline SolverResult solve_cocr(
    const SparseMatrixView<Scalar>& A,
    const Scalar* b,
    Scalar* x,
    index_t size,
    const SolverConfig& config = SolverConfig()
) {
    ICPreconditioner<Scalar> precond(config.shift_parameter);
    precond.set_config(config);
    precond.setup(A);

    COCRSolver<Scalar> solver;
    solver.set_config(config);
    return solver.solve(A, b, x, size, &precond);
}

/// Convenience function: Solve Ax=b using COCR with std::vector
template<typename Scalar = double>
[[nodiscard]] inline SolverResult solve_cocr(
    const SparseMatrixView<Scalar>& A,
    const std::vector<Scalar>& b,
    std::vector<Scalar>& x,
    const SolverConfig& config = SolverConfig()
) {
    if (x.size() != b.size()) {
        x.resize(b.size());
    }
    return solve_cocr(A, b.data(), x.data(), static_cast<index_t>(b.size()), config);
}

} // namespace sparsesolv

#endif // SPARSESOLV_HPP
