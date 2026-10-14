/*-------------------------------------------------------------------------
 * C++ implementation of TaskManager wrappers for HACApK C code
 *
 * These functions are called from cHACApK_cpp_impl.c and cHACApK_base.c
 * via the C-compatible interface in rad_hacapk_parallel.h.
 *-------------------------------------------------------------------------*/

#include "rad_hacapk_parallel.h"
#include "rad_parallel.h"
#include <algorithm>

// Tasks per worker for dynamic load balancing of the leaf fill (see hacapk_parallel_for).
static constexpr int HACAPK_TASKS_PER_THREAD = 32;

extern "C" {

int hacapk_get_num_threads(void)
{
    return ngcore::TaskManager::GetNumThreads();
}

int hacapk_get_thread_id(void)
{
    return ngcore::TaskManager::GetThreadId();
}

void hacapk_parallel_for(int n,
    void (*func)(int idx, void* data), void* data)
{
    if (n <= 0) return;
    // Many tasks per thread, not one: ngcore::ParallelFor's default splits the range into ONE
    // contiguous chunk per thread, which is a static schedule.  The H-matrix leaf fill is the main
    // caller and its leaves are wildly uneven (a dense near leaf of Duffy pair blocks costs seconds,
    // a far low-rank leaf microseconds) and sorted by row block, so the static split left most
    // workers idle while a few finished their chunk: on the CEFC 2020 quadrupole build the
    // quadrature branches summed to a quarter of the thread capacity.  Tasks are pulled from a
    // shared counter, so a fine split balances the load; the per-task overhead is microseconds.
    const int nthr = std::max(1, ngcore::TaskManager::GetNumThreads());
    const int ntasks = std::max(1, std::min(n, HACAPK_TASKS_PER_THREAD * nthr));
    ngcore::ParallelFor(ngcore::IntRange(n), [&](size_t i) {
        func((int)i, data);
    }, ntasks);
}

void hacapk_parallel_job(
    void (*func)(int tid, int nthr, void* data), void* data)
{
    int nthr = ngcore::TaskManager::GetNumThreads();
    if (nthr <= 1 || !ngcore::GetTaskManager()) {
        // No task manager or single thread → run sequentially
        func(0, 1, data);
        return;
    }
    // Use ParallelFor(nthr) instead of ParallelJob to ensure each "thread"
    // gets a unique tid and nthr matches the buffer allocation.
    // ParallelJob may report different ti.nthreads than GetNumThreads().
    ngcore::ParallelFor(ngcore::IntRange(nthr), [&](size_t i) {
        func((int)i, nthr, data);
    });
}

} // extern "C"
