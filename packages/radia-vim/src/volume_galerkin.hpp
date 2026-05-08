// volume_galerkin.hpp — full K-matrix assembler for the HDiv div-free
// hex BEM (Phase F-4).
//
// Loops over all unique (i, j) basis pairs and calls
// `compute_K_entry_spherical_duffy` for each. Returns the dense, symmetric
// K matrix as a flat std::vector<T> in row-major order with shape
// (n_dofs, n_dofs). The pybind11 binding wraps it in a NumPy array.
//
// Symmetry: the kernel 1/|r-r'| is symmetric under (r ↔ r'), so K is
// symmetric. We compute K[i,j] for i ≤ j (n_dofs * (n_dofs + 1) / 2 unique
// entries) and mirror the result.
//
// Parallelism: OpenMP over a flattened pair-index `pair_id` ∈ [0, n_pairs).
// Each entry is independent — no shared state writes inside the loop.
//
// The bare quadrature does NOT include the (mu0 / 4π) prefactor; downstream
// Python code multiplies after the fact (matches the convention of
// `compute_K_entry_spherical_duffy` and `K_reference_values.json`).

#pragma once

#include <vector>
#include <utility>
#include <cstddef>

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "hdiv_divfree_hex_basis.hpp"
#include "sauter_schwab_hex.hpp"

namespace radia_vim {

// Pair-index helpers: flatten the upper triangle (i <= j) of an n × n
// matrix. There are n*(n+1)/2 unique pairs. Mapping:
//   pair_id = i * n - i*(i-1)/2 + (j - i)         for 0 <= i <= j < n
inline std::pair<int, int> pair_index_to_ij(std::size_t pair_id, int n) {
    // Solve quadratic: pair_id = i*(2n - i + 1)/2 + (j - i) where
    // (j - i) ∈ [0, n - i). Equivalent: find largest i s.t. T(i) <= pair_id
    // where T(i) = i*(2n - i + 1)/2 = sum_{k=0}^{i-1} (n - k).
    // Use float estimate then linear adjust.
    long long pid = static_cast<long long>(pair_id);
    long long N   = static_cast<long long>(n);
    // Rough estimate: i ≈ n - sqrt(n^2 - 2*pid)
    double est = static_cast<double>(N) - 0.5
               - 0.5 * std::sqrt(4.0 * N * (N + 1.0) - 8.0 * static_cast<double>(pid) - 7.0);
    long long i = static_cast<long long>(est);
    if (i < 0)         i = 0;
    if (i >= N)        i = N - 1;
    auto T = [&](long long ii) { return ii * (2 * N - ii + 1) / 2; };
    while (i > 0 && T(i) > pid) --i;
    while (i + 1 < N && T(i + 1) <= pid) ++i;
    long long j = i + (pid - T(i));
    if (j >= N) j = N - 1;
    return {static_cast<int>(i), static_cast<int>(j)};
}

inline std::size_t n_unique_pairs(int n) {
    return static_cast<std::size_t>(n) * static_cast<std::size_t>(n + 1) / 2;
}

// -----------------------------------------------------------------------
// Assemble the bare (i.e. no mu0/4pi prefactor) K matrix for the cuboid
// (a × b × c). Returns row-major flat vector of length n_dofs * n_dofs.
// -----------------------------------------------------------------------
template <typename T>
std::vector<T> assemble_K_bare(
    const HDivDivFreeHexBasis<T>& basis,
    T a, T b, T c,
    int n_omega_theta = 16,
    int n_omega_phi   = 32,
    int n_rho         = 16,
    int n_c           = 8,
    int verbose       = 0) {

    const int n = basis.n_dofs();
    std::vector<T> K(static_cast<std::size_t>(n) * static_cast<std::size_t>(n),
                     T(0));
    const std::size_t n_pairs = n_unique_pairs(n);

    // Counter for progress reporting under OpenMP.
    std::size_t done = 0;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (long long pid = 0; pid < static_cast<long long>(n_pairs); ++pid) {
        auto [i, j] = pair_index_to_ij(static_cast<std::size_t>(pid), n);
        T val = compute_K_entry_spherical_duffy<T>(
            basis, i, j, a, b, c,
            n_omega_theta, n_omega_phi, n_rho, n_c);
        K[static_cast<std::size_t>(i) * n + j] = val;
        if (i != j)
            K[static_cast<std::size_t>(j) * n + i] = val;

        if (verbose > 0) {
#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                ++done;
                if (verbose >= 2 ||
                    (done % static_cast<std::size_t>(verbose) == 0) ||
                    done == n_pairs) {
                    // We use printf rather than std::cout to avoid std::sync overhead.
                    std::fprintf(stderr,
                        "  K assembler: %zu / %zu  (%.1f%%)\n",
                        done, n_pairs,
                        100.0 * static_cast<double>(done) /
                                static_cast<double>(n_pairs));
                    std::fflush(stderr);
                }
            }
        }
    }

    return K;
}

}  // namespace radia_vim
