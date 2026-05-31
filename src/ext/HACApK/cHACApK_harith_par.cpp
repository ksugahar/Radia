/* cHACApK_harith_par.cpp
 *
 * C-callable bridge from the (C) H-LU arithmetic in cHACApK_harith.c to
 * NGSolve's ngcore TaskManager, for BLOCK-LEVEL parallelism of the
 * block-recursive H-matrix multiply (h_addmul).
 *
 * Rationale (2026-05-31): the H-LU leaf BLAS are too small (~60 DOF blocks)
 * to trigger MKL's internal threading, so the factorization runs effectively
 * single-threaded while MKL dense LU uses all cores. The only way to use the
 * cores is task parallelism over INDEPENDENT block-tree operations. HACApK
 * clusters are binary (s=2), so the parallelism lives in the recursive
 * h_addmul descent: the output sub-blocks C(i,k) are independent (the inner
 * j-sum that accumulates into one C(i,k) stays serial), so we parallelize
 * over the (i,k) grid with a size cutoff and let ngcore's work-stealing
 * handle the nested recursion.
 *
 * Policy note: this honours the TaskManager-only policy (no raw OpenMP); the
 * ngcore ParallelFor uses the same NGSolve threadpool. Without an enclosing
 * RegionTaskManager, ParallelFor falls back to a single-threaded loop, so the
 * decomp wraps hlu_rec in chacapk_par_region to bring the pool up once. */

#include <core/taskmanager.hpp>

extern "C" {

int chacapk_max_threads(void)
{
    int n = ngcore::TaskManager::GetMaxThreads();
    return n > 0 ? n : 1;
}

/* Run body(ctx) inside an active TaskManager region so nested chacapk_par_for
 * calls actually parallelize (RAII brings the pool up for the region). */
void chacapk_par_region(void (*body)(void*), void *ctx)
{
    int nthr = chacapk_max_threads();
    ngcore::RegionTaskManager rtm(nthr);
    body(ctx);
}

/* Parallel for over [0,n): body(i, ctx) for each i. Serial fallback if no
 * TaskManager region is active. */
void chacapk_par_for(int n, void (*body)(int, void*), void *ctx)
{
    if (n <= 0) return;
    ngcore::ParallelFor(ngcore::IntRange((size_t)n), [&](size_t i) {
        body((int)i, ctx);
    });
}

} /* extern "C" */
