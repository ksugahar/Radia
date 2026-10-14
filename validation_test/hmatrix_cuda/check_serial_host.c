#include "rad_hacapk_parallel.h"
#include <stdio.h>
struct State { int next, jobs, failed; };
static void step(int i,void* data) {
    struct State* s=(struct State*)data;
    if(i!=s->next++ || hacapk_get_thread_id()!=0 || hacapk_get_num_threads()!=1) s->failed=1;
}
static void job(int tid,int threads,void* data) {
    struct State* s=(struct State*)data;
    if(tid!=0 || threads!=1) s->failed=1;
    ++s->jobs;
}
int main(void) {
    struct State s={0,0,0};
    hacapk_parallel_for(0,step,&s);
    if(s.next) return 1;
    hacapk_parallel_for(257,step,&s);
    hacapk_parallel_job(job,&s);
    if(s.failed || s.next!=257 || s.jobs!=1) return 2;
    puts("{\"iterations\":257,\"jobs\":1,\"threads\":1,\"taskmanager_used\":false}");
    return 0;
}
