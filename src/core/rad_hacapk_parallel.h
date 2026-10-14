#ifndef RAD_HACAPK_PARALLEL_H
#define RAD_HACAPK_PARALLEL_H

/*-------------------------------------------------------------------------
 * C-compatible wrappers for NGSolve TaskManager parallelization
 *
 * HACApK is written in C and cannot directly call C++ TaskManager API.
 * These wrapper functions bridge the C/C++ boundary.
 *
 * Implementation: rad_hacapk_parallel.cpp (compiled as C++)
 *-------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/* Get number of available threads from TaskManager */
int hacapk_get_num_threads(void);

/* Get current thread ID from TaskManager */
int hacapk_get_thread_id(void);

/* Parallel for loop: calls func(idx, data) for idx in [0, n)
 * Replacement for: #pragma omp parallel for */
void hacapk_parallel_for(int n,
    void (*func)(int idx, void* data), void* data);

/* Parallel job: calls func(tid, nthr, data) on each thread
 * Each thread gets a unique tid in [0, nthr).
 * Replacement for: #pragma omp parallel { ... omp_get_thread_num() ... }
 * Thread function should internally split work using tid/nthr. */
void hacapk_parallel_job(
    void (*func)(int tid, int nthr, void* data), void* data);

#ifdef __cplusplus
}
#endif

#endif /* RAD_HACAPK_PARALLEL_H */
