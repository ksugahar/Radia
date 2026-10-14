#pragma once

/*-------------------------------------------------------------------------
 * Radia parallel abstraction layer
 *
 * Uses NGSolve TaskManager for all parallelism.
 * NGSolve is a hard prerequisite for this fork — no OpenMP fallback.
 *
 * Replaces all #pragma omp directives with TaskManager equivalents:
 *   #pragma omp parallel for       -> ParallelFor(IntRange(n), lambda)
 *   #pragma omp parallel for       -> ParallelForRange(IntRange(n), lambda)
 *     schedule(dynamic)               (TaskManager uses work-stealing by default)
 *   #pragma omp critical           -> std::mutex / AtomicMax / AtomicMin
 *   #pragma omp atomic             -> AtomicAdd
 *   #pragma omp parallel for       -> ParallelForRange + local accum + AtomicAdd
 *     reduction(+:sum)
 *   omp_get_max_threads()          -> TaskManager::GetNumThreads()
 *   omp_get_thread_num()           -> TaskManager::GetThreadId()
 *-------------------------------------------------------------------------*/

#include <core/taskmanager.hpp>    // ParallelFor, ParallelForRange, ParallelJob
#include <core/utils.hpp>          // AtomicAdd, AtomicMax, AtomicMin, SuspendTaskManager
#include <algorithm>
#include <cstdlib>

#ifdef HAVE_LAPACK
#include "mkl_service.h"           // mkl_set_num_threads, mkl_get_max_threads
#endif

namespace radia {

// Re-export commonly used TaskManager names into radia namespace
using ngcore::ParallelFor;
using ngcore::ParallelForRange;
using ngcore::ParallelJob;
using ngcore::AtomicAdd;
using ngcore::AtomicMax;
using ngcore::AtomicMin;
using ngcore::TaskManager;
using ngcore::SuspendTaskManager;
using ngcore::IntRange;
using ngcore::TaskInfo;

#ifdef HAVE_LAPACK
// RAII guard: temporarily sets MKL threads for LAPACK calls (dgesv_ etc.),
// restoring the original value on scope exit.
// NGSolve sets mkl_set_num_threads(1) globally; this guard re-enables
// multi-threaded MKL for computationally intensive LAPACK operations.
class MKLThreadGuard {
    int saved_;
public:
    explicit MKLThreadGuard(int n) : saved_(mkl_get_max_threads()) { mkl_set_num_threads(n); }
    ~MKLThreadGuard() { mkl_set_num_threads(saved_); }
    MKLThreadGuard(const MKLThreadGuard&) = delete;
    MKLThreadGuard& operator=(const MKLThreadGuard&) = delete;
};
#endif

// Convenience wrappers.  RADIA_NUM_THREADS is intentionally read on demand so
// benchmark drivers can set it before calling rad.Solve without rebuilding or
// relying on the host Python's NGSolve default.
inline int GetMaxThreads()
{
    if(const char* env = std::getenv("RADIA_NUM_THREADS"))
    {
        int n = std::atoi(env);
        if(n > 0) return n;
    }
    return std::max(1, TaskManager::GetMaxThreads());
}

inline int GetNumThreads()
{
    int n = TaskManager::GetNumThreads();
    return (n > 0) ? n : GetMaxThreads();
}
inline int GetThreadId()   { return TaskManager::GetThreadId(); }

} // namespace radia
