/*-------------------------------------------------------------------------
 * C++ implementation of TaskManager wrappers for HACApK C code
 *
 * These functions are called from cHACApK_cpp_impl.c and cHACApK_base.c
 * via the C-compatible interface in rad_hacapk_parallel.h.
 *-------------------------------------------------------------------------*/

#include "rad_hacapk_parallel.h"
#include "rad_parallel.h"

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
    ngcore::ParallelFor(ngcore::IntRange(n), [&](size_t i) {
        func((int)i, data);
    });
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
