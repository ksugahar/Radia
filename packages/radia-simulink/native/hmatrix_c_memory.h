/* MEX-only C allocation transaction. Never include its macros in C++ code. */
#ifndef RADIA_HMATRIX_C_MEMORY_H
#define RADIA_HMATRIX_C_MEMORY_H
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {
    size_t allocations, frees, peak_bytes, live_bytes, fail_at, transferred_bytes;
    const char *last_file;
    const char *last_kind;
    int last_line, failure_code;
} rad_hmatrix_c_memory_report;
void rad_hmatrix_c_memory_configure(size_t fail_at);
rad_hmatrix_c_memory_report rad_hmatrix_c_memory_last(void);
int rad_hmatrix_guarded_build(void*,void*,double*,int,int,int,double,int,double,int,int);
void *rad_hmatrix_c_malloc(size_t,const char*,int);
void *rad_hmatrix_c_calloc(size_t,size_t,const char*,int);
void *rad_hmatrix_c_realloc(void*,size_t,const char*,int);
void rad_hmatrix_c_free(void*);
void rad_hmatrix_c_exit(int,const char*,int);
#ifdef __cplusplus
}
#elif defined(RAD_HACAPK_C_MEMORY)
/* This opt-in is passed only to the serial HACApK C translation units.
 * setjmp/longjmp never spans the C++ owner or a CUDA entry callback. */
#define malloc(n) rad_hmatrix_c_malloc((n),__FILE__,__LINE__)
#define calloc(n,s) rad_hmatrix_c_calloc((n),(s),__FILE__,__LINE__)
#define realloc(p,n) rad_hmatrix_c_realloc((p),(n),__FILE__,__LINE__)
#define free(p) rad_hmatrix_c_free(p)
#define exit(n) rad_hmatrix_c_exit((n),__FILE__,__LINE__)
#endif
#endif
