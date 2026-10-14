/* SPDX-License-Identifier: MPL-2.0 */
#pragma once
#include <core/taskmanager.hpp>
#include <algorithm>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <utility>
#ifdef SPARSESOLV_OPENMP_AMS_SETUP
#ifndef SPARSESOLV_SERIAL_AMS_SETUP
#error OpenMP snapshot setup must disable internally owned TaskManager regions
#endif
#ifndef _OPENMP
#error SPARSESOLV_OPENMP_AMS_SETUP requires OpenMP compiler support
#endif
#include <omp.h>
#endif

namespace ngla {
// Only the Python snapshot adapter defines this macro. CUDA MEX sources
// include neither this header nor NGSolve. Normal CPU dispatch is unchanged.
inline int& CompactSetupThreads() { static thread_local int value = 1; return value; }
struct CompactSetupThreadScope {
    int previous;
    explicit CompactSetupThreadScope(int threads) : previous(CompactSetupThreads()) {
        if (threads < 1 || threads > 64) throw std::invalid_argument("AMS setup threads must be in [1,64]");
#ifndef SPARSESOLV_OPENMP_AMS_SETUP
        if (threads != 1) throw std::invalid_argument("AMS snapshot was built without OpenMP setup");
#endif
        CompactSetupThreads() = threads;
    }
    ~CompactSetupThreadScope() { CompactSetupThreads() = previous; }
};
inline int CompactWorkerCount() {
#ifdef SPARSESOLV_OPENMP_AMS_SETUP
    return CompactSetupThreads();
#else
    return ngcore::TaskManager::GetNumThreads();
#endif
}
inline int CompactWorkerId() {
#ifdef SPARSESOLV_OPENMP_AMS_SETUP
    return omp_in_parallel() ? omp_get_thread_num() : 0;
#else
    return ngcore::TaskManager::GetThreadId();
#endif
}
template<class F> inline void CompactParallelFor(size_t n, F&& f, size_t work = 0) {
#ifdef SPARSESOLV_OPENMP_AMS_SETUP
    const int threads = CompactSetupThreads();
    if (threads == 1 || (work ? work : n) < 1024 || omp_in_parallel()) {
        for (size_t i = 0; i < n; ++i) f(i);
        return;
    }
    std::exception_ptr error;
    std::mutex error_mutex;
    // Row reductions stay in their original order. Disjoint rows/chunks
    // alone run concurrently; exceptions are rethrown on the calling thread.
#pragma omp parallel for schedule(static) num_threads(threads)
    for (long long i = 0; i < static_cast<long long>(n); ++i) {
        try { f(static_cast<size_t>(i)); }
        catch (...) {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!error) error = std::current_exception();
        }
    }
    if (error) std::rethrow_exception(error);
#else
    ngcore::ParallelFor(n, std::forward<F>(f));
#endif
}
template<class F> inline void CompactParallelForRange(size_t n, F&& f) {
#ifdef SPARSESOLV_OPENMP_AMS_SETUP
    const size_t chunks = std::max<size_t>(1, std::min<size_t>(n, 256));
    CompactParallelFor(chunks, [&](size_t c) {
        f(ngcore::IntRange(n * c / chunks, n * (c + 1) / chunks));
    }, n);
#else
    ngcore::ParallelForRange(n, std::forward<F>(f));
#endif
}
}
