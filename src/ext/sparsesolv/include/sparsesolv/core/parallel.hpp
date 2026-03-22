/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file parallel.hpp
 * @brief Portable parallel primitives for SparseSolv
 *
 * Provides compile-time dispatch to:
 * - NGSolve TaskManager (when SPARSESOLV_USE_NGSOLVE_TASKMANAGER is defined)
 * - OpenMP (when _OPENMP is defined)
 * - Serial fallback (otherwise)
 *
 * This allows SparseSolv to work both as a standalone library and
 * integrated into NGSolve without OpenMP/TaskManager thread pool conflicts.
 */

#ifndef SPARSESOLV_CORE_PARALLEL_HPP
#define SPARSESOLV_CORE_PARALLEL_HPP

#include "types.hpp"
#include <type_traits>
#include <complex>
#include <atomic>
#include <thread>

#ifdef SPARSESOLV_USE_NGSOLVE_TASKMANAGER
  #include <core/taskmanager.hpp>
#elif defined(_OPENMP)
  #include <omp.h>
#endif

namespace sparsesolv {

/**
 * @brief Parallel for loop over [0, n)
 *
 * Dispatches to ngcore::ParallelFor, OpenMP, or serial loop.
 */
template<typename FUNC>
inline void parallel_for(index_t n, FUNC f) {
#ifdef SPARSESOLV_USE_NGSOLVE_TASKMANAGER
    if (n > 0)
        ngcore::ParallelFor(static_cast<size_t>(n),
            [&](size_t i) { f(static_cast<index_t>(i)); });
#elif defined(_OPENMP)
    #pragma omp parallel for schedule(static)
    for (index_t i = 0; i < n; ++i)
        f(i);
#else
    for (index_t i = 0; i < n; ++i)
        f(i);
#endif
}

/**
 * @brief Parallel reduction (summation) over [0, n)
 *
 * Computes sum of f(i) for i in [0, n).
 * Uses ngcore::ParallelReduce (supports all types including complex),
 * OpenMP reduction (double only on MSVC), or serial loop.
 *
 * @tparam T Result type (double, complex<double>, etc.)
 * @tparam FUNC Function type: index_t -> T
 * @param n Loop bound
 * @param f Function to evaluate at each index
 * @param init Initial value for reduction
 * @return Sum of f(i) for i in [0, n)
 */
template<typename T, typename FUNC>
inline T parallel_reduce_sum(index_t n, FUNC f, T init = T(0)) {
#ifdef SPARSESOLV_USE_NGSOLVE_TASKMANAGER
    if (n <= 0) return init;
    return ngcore::ParallelReduce(static_cast<size_t>(n),
        [&](size_t i) { return f(static_cast<index_t>(i)); },
        [](T a, T b) { return a + b; }, init);
#elif defined(_OPENMP)
    T sum = init;
    // MSVC OpenMP (v2.0) only supports reduction on scalar arithmetic types
    if constexpr (std::is_same_v<T, double>) {
        #pragma omp parallel for reduction(+:sum)
        for (index_t i = 0; i < n; ++i)
            sum += f(i);
    } else {
        // Complex or other types: serial fallback on MSVC
        for (index_t i = 0; i < n; ++i)
            sum += f(i);
    }
    return sum;
#else
    T sum = init;
    for (index_t i = 0; i < n; ++i)
        sum += f(i);
    return sum;
#endif
}

/**
 * @brief Get the current thread index within a parallel region
 *
 * Returns a value in [0, get_num_threads()) that identifies the calling thread.
 */
inline int get_thread_id() {
#ifdef SPARSESOLV_USE_NGSOLVE_TASKMANAGER
    return ngcore::TaskManager::GetThreadId();
#elif defined(_OPENMP)
    return omp_get_thread_num();
#else
    return 0;
#endif
}

/**
 * @brief Get the number of available parallel threads
 */
inline int get_num_threads() {
#ifdef SPARSESOLV_USE_NGSOLVE_TASKMANAGER
    return ngcore::TaskManager::GetNumThreads();
#elif defined(_OPENMP)
    return omp_get_max_threads();
#else
    return 1;
#endif
}

/**
 * @brief Spin-barrier for synchronizing threads within a persistent parallel region
 *
 * Uses sense-reversing algorithm with std::atomic for C++17 compatibility.
 * Designed for use inside a single parallel_for(nthreads, ...) where each
 * index maps to a separate thread, enabling level-by-level synchronization
 * without repeated parallel_for dispatch overhead.
 *
 * Requirements:
 * - Must be called with exactly num_threads concurrent threads
 * - All threads must call wait() the same number of times
 */
class SpinBarrier {
public:
    explicit SpinBarrier(int num_threads)
        : num_threads_(num_threads) {
        if (num_threads < 1)
            throw std::invalid_argument("SpinBarrier: num_threads must be >= 1, got "
                + std::to_string(num_threads));
        count_.store(0, std::memory_order_relaxed);
        sense_.store(0, std::memory_order_relaxed);
    }

    void wait() {
        int my_sense = sense_.load(std::memory_order_acquire);
        if (count_.fetch_add(1, std::memory_order_acq_rel) == num_threads_ - 1) {
            // Last thread to arrive: reset counter and flip sense
            count_.store(0, std::memory_order_relaxed);
            sense_.fetch_add(1, std::memory_order_release);
        } else {
            // Spin until the last thread flips the sense
            int spins = 0;
            while (sense_.load(std::memory_order_acquire) == my_sense) {
                if (++spins > 4096) {
                    std::this_thread::yield();
                    spins = 0;
                }
            }
        }
    }

private:
    // Separate cache lines to avoid false sharing
    alignas(64) std::atomic<int> count_;
    alignas(64) std::atomic<int> sense_;
    int num_threads_;
};

} // namespace sparsesolv

#endif // SPARSESOLV_CORE_PARALLEL_HPP
