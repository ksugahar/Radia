/* Deliberately serial MATLAB host scheduler for native HACApK assembly.
 * Link instead of rad_hacapk_parallel.cpp, never alongside it. Keeping the
 * entry callback on its owning thread preserves scoped CUDA binding and
 * avoids oversubscription while CUDA performs the expensive integration.
 */
#include "rad_hacapk_parallel.h"
int hacapk_get_num_threads(void) { return 1; }
int hacapk_get_thread_id(void) { return 0; }
void hacapk_parallel_for(int n, void (*func)(int, void*), void* data) {
    int i;
    for(i=0;i<n;++i) func(i,data);
}
void hacapk_parallel_job(void (*func)(int,int,void*), void* data) {
    func(0,1,data);
}
