/* Own only the serial C build invocation. No MATLAB API or C++ frames here. */
#include "hmatrix_c_memory.h"
#include <setjmp.h>
#include <string.h>
#include <stdint.h>
#include "cHACApK_base.h"
#ifdef _MSC_VER
#define RAD_TLS __declspec(thread)
#else
#define RAD_TLS _Thread_local
#endif
typedef struct allocation { void *pointer; size_t bytes; struct allocation *next; } allocation;
typedef struct { jmp_buf jump; allocation *head; rad_hmatrix_c_memory_report report; } transaction;
static RAD_TLS transaction *active;
static RAD_TLS size_t requested_failure;
static RAD_TLS rad_hmatrix_c_memory_report last_report;
extern int HACApK_build_hmatrix_wrapper(void*,void*,double*,int,int,int,double,int,double,int,int);

void rad_hmatrix_c_memory_configure(size_t fail_at) { requested_failure=fail_at; }
rad_hmatrix_c_memory_report rad_hmatrix_c_memory_last(void) { return last_report; }
static void fail(int code,const char *file,int line) {
    active->report.failure_code=code;
    active->report.last_file=file; active->report.last_line=line;
    longjmp(active->jump,1);
}
static int attempt(const char *file,int line,const char *kind){
    active->report.last_file=file;active->report.last_line=line;active->report.last_kind=kind;
    ++active->report.allocations;
    return !(active->report.fail_at && active->report.allocations==active->report.fail_at);
}
static void account(size_t bytes){
    active->report.live_bytes+=bytes;
    if(active->report.live_bytes>active->report.peak_bytes)
        active->report.peak_bytes=active->report.live_bytes;
}
static void *allocate(size_t bytes,int zero,const char *file,int line) {
    void *p; allocation *node;
    if(!bytes)bytes=1;
    if(!active) return zero?calloc(1,bytes):malloc(bytes);
    if(active->report.live_bytes>SIZE_MAX-sizeof(*node)||bytes>SIZE_MAX-active->report.live_bytes-sizeof(*node))fail(2,file,line);
    p=attempt(file,line,"payload")?(zero?calloc(1,bytes?bytes:1):malloc(bytes?bytes:1)):NULL;
    if(!p) fail(1,file,line);
    account(bytes);
    node=attempt(file,line,"ledger")?(allocation*)malloc(sizeof(*node)):NULL;
    if(!node){free(p);active->report.live_bytes-=bytes;++active->report.frees;fail(1,file,line);}
    account(sizeof(*node));
    node->pointer=p;node->bytes=bytes;node->next=active->head;active->head=node;
    return p;
}
void *rad_hmatrix_c_malloc(size_t n,const char *f,int l){return allocate(n,0,f,l);}
void *rad_hmatrix_c_calloc(size_t n,size_t s,const char *f,int l){
    if(s && n>SIZE_MAX/s){if(active)fail(2,f,l);return NULL;}
    return allocate(n*s,1,f,l);
}
void rad_hmatrix_c_free(void *p){
    allocation **at;
    if(!p)return;
    if(active){
        at=&active->head;
        while(*at && (*at)->pointer!=p)at=&(*at)->next;
        if(*at){allocation *node=*at;*at=node->next;
            active->report.live_bytes-=node->bytes+sizeof(*node);active->report.frees+=2;free(node);}
    }
    free(p);
}
void *rad_hmatrix_c_realloc(void *p,size_t n,const char *f,int l){
    allocation *node;void *q;
    if(!active)return realloc(p,n);
    if(!p)return allocate(n,0,f,l);
    if(!n){rad_hmatrix_c_free(p);return NULL;}
    for(node=active->head;node && node->pointer!=p;node=node->next){}
    /* Never guess the extent of an allocation owned outside this transaction. */
    if(!node)fail(3,f,l);
    q=allocate(n,0,f,l);memcpy(q,p,n<node->bytes?n:node->bytes);
    rad_hmatrix_c_free(p);return q;
}
void rad_hmatrix_c_exit(int status,const char *f,int l){
    if(active)fail(status?4:5,f,l);
    /* No intercepted build function is allowed outside the C boundary. */
    abort();
}
int rad_hmatrix_guarded_build(void *m,void *c,double *xyz,int n,int nffc,int dim,
                             double eps,int leaf,double eta,int rank,int print){
    transaction *tx;int result;allocation *node;
    if(active)return -91;
    /* Heap storage avoids indeterminate automatic variables after longjmp. */
    memset(&last_report,0,sizeof(last_report));
    tx=requested_failure==1?NULL:(transaction*)calloc(1,sizeof(*tx));
    if(!tx){last_report.allocations=1;last_report.fail_at=requested_failure;
        last_report.failure_code=1;last_report.last_file=__FILE__;last_report.last_line=__LINE__;
        last_report.last_kind="transaction";return -92;}
    tx->report.fail_at=requested_failure;tx->report.allocations=1;active=tx;account(sizeof(*tx));
    if(setjmp(tx->jump)==0)
        result=HACApK_build_hmatrix_wrapper(m,c,xyz,n,nffc,dim,eps,leaf,eta,rank,print);
    else result=-93;
    if(result){
        while(tx->head)rad_hmatrix_c_free(tx->head->pointer);
        /* The caller owns these two outer objects; all their children were
         * rolled back, so its ordinary C++ destructor must see empty objects. */
        memset(m,0,sizeof(st_cHACApK_leafmtxp_t));
        memset(c,0,sizeof(st_cHACApK_lcontrol_t));
    }else{
        /* Successful allocations transfer to native ownership, without copying.
         * Only bookkeeping is discarded; normal HACApK free remains valid. */
        while(tx->head){node=tx->head;tx->head=node->next;
            tx->report.live_bytes-=sizeof(*node);++tx->report.frees;free(node);}
        tx->report.transferred_bytes=tx->report.live_bytes-sizeof(*tx);
        tx->report.live_bytes=sizeof(*tx);
    }
    tx->report.live_bytes-=sizeof(*tx);++tx->report.frees;
    last_report=tx->report;active=NULL;free(tx);return result;
}
