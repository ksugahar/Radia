#pragma once
/* Legacy transform interface over the selected NVPL FFTW3 float API.
 * Forward and backward transforms are both unnormalised, matching the
 * unchanged MKL descriptor path; callers in gmfft own physical scaling. */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <vector>
#include <stdexcept>
#include <limits>
namespace radia_nvpl_fft {
#include <nvpl_fftw/fftw3.h>
constexpr unsigned plan_flags = FFTW_ESTIMATE | FFTW_UNALIGNED;
}
#undef FFTW_FORWARD
#undef FFTW_BACKWARD
#undef FFTW_ESTIMATE
#undef FFTW_MEASURE
#define FFTW_ESTIMATE 0
#define FFTW_MEASURE 1
#define FFTW_OUT_OF_PLACE 0
#define FFTW_IN_PLACE 8
#define FFTW_USE_WISDOM 16
#define FFTW_THREADSAFE 128
#define FFTWND_FORCE_BUFFERED 256
#define FFTW_NO_VECTOR_RECURSE 512
#define FFTW_ENABLE_FLOAT 1
typedef float fftw_real;
typedef struct { fftw_real re, im; } fftw_complex;
#define c_re(c) ((c).re)
#define c_im(c) ((c).im)
typedef enum { FFTW_FORWARD=-1, FFTW_BACKWARD=1 } fftw_direction;
typedef enum { FFTW_SUCCESS=0, FFTW_FAILURE=-1 } fftw_status;
typedef fftw_complex FFTW_COMPLEX;
typedef fftw_real FFTW_REAL;
struct radia_nvpl_fft_plan {
    std::vector<int> dims;
    fftw_direction direction;
    bool inplace;
};
typedef radia_nvpl_fft_plan* fftw_plan;
typedef radia_nvpl_fft_plan* fftwnd_plan;

static inline fftwnd_plan fftwnd_create_plan(int rank, const int* dims, fftw_direction dir, int flags) {
    if(rank<=0 || !dims || (dir!=FFTW_FORWARD && dir!=FFTW_BACKWARD))
        throw std::invalid_argument("Invalid NVPL FFT plan");
    size_t size=1;
    for(int i=0;i<rank;++i) {
        if(dims[i]<=0 || size>std::numeric_limits<size_t>::max()/static_cast<size_t>(dims[i]))
            throw std::invalid_argument("Invalid NVPL FFT dimensions");
        size*=static_cast<size_t>(dims[i]);
    }
    return new radia_nvpl_fft_plan{std::vector<int>(dims,dims+rank),dir,(flags & FFTW_IN_PLACE)!=0};
}
static inline fftw_plan fftw_create_plan(int n, fftw_direction dir, int flags) {
    return fftwnd_create_plan(1,&n,dir,flags);
}
static inline fftwnd_plan fftw2d_create_plan(int nx,int ny,fftw_direction dir,int flags) {
    const int dims[]={nx,ny}; return fftwnd_create_plan(2,dims,dir,flags);
}
static inline fftwnd_plan fftw3d_create_plan(int nx,int ny,int nz,fftw_direction dir,int flags) {
    const int dims[]={nx,ny,nz}; return fftwnd_create_plan(3,dims,dir,flags);
}
static inline void fftwnd_destroy_plan(fftwnd_plan p) { delete p; }
static inline void fftw_destroy_plan(fftw_plan p) { delete p; }
static inline void fftwnd(fftwnd_plan p,int howmany,fftw_complex* in,int istride,int idist,
                          fftw_complex* out,int ostride,int odist) {
    if(!p || !in || howmany<1 || istride<1 || ostride<1 || idist<0 || odist<0)
        throw std::invalid_argument("Invalid NVPL FFT execution");
    if(p->inplace) out=in;
    if(!out) throw std::invalid_argument("NVPL FFT output is required");
    static_assert(sizeof(fftw_complex)==sizeof(radia_nvpl_fft::fftwf_complex),"Float interleaved complex ABI");
    // ESTIMATE planning leaves input untouched; no measured plan can overwrite
    // caller data. Actual strides/batch distances are owned by this invocation.
    auto plan=radia_nvpl_fft::fftwf_plan_many_dft(static_cast<int>(p->dims.size()),p->dims.data(),howmany,
        reinterpret_cast<radia_nvpl_fft::fftwf_complex*>(in),nullptr,istride,idist,
        reinterpret_cast<radia_nvpl_fft::fftwf_complex*>(out),nullptr,ostride,odist,p->direction,radia_nvpl_fft::plan_flags);
    if(!plan) throw std::runtime_error("NVPL FFT plan creation failed");
    radia_nvpl_fft::fftwf_execute(plan);
    radia_nvpl_fft::fftwf_destroy_plan(plan);
}
static inline void fftw(fftw_plan p,int n,fftw_complex* in,int is,int id,fftw_complex* out,int os,int od) {
    fftwnd(p,n,in,is,id,out,os,od);
}
static inline void fftw_one(fftw_plan p,fftw_complex* in,fftw_complex* out) { fftw(p,1,in,1,0,out,1,0); }
static inline void fftwnd_one(fftwnd_plan p,fftw_complex* in,fftw_complex* out) { fftwnd(p,1,in,1,0,out,1,0); }
static inline void fftwnd_fprint_plan(FILE* f,fftwnd_plan p) {
    if(!f || !p) throw std::invalid_argument("Invalid FFT plan print");
    fprintf(f,"NVPL FFT rank=%zu direction=%d inplace=%d dims=",p->dims.size(),p->direction,p->inplace);
    for(auto n:p->dims) fprintf(f,"%d ",n);
    fprintf(f,"\n");
}
static inline void fftw_fprint_plan(FILE* f,fftw_plan p) { fftwnd_fprint_plan(f,p); }
static inline void fftw_print_plan(fftw_plan p) { fftw_fprint_plan(stdout,p); }
static inline void fftwnd_print_plan(fftwnd_plan p) { fftwnd_fprint_plan(stdout,p); }
static inline void* fftw_malloc(size_t n) { return malloc(n); }
static inline void fftw_free(void* p) { free(p); }
static inline size_t fftw_sizeof_fftw_real(void) { return sizeof(fftw_real); }
static inline void fftw_die(const char* message) { fprintf(stderr,"NVPL FFT: %s\n",message); abort(); }
static inline void fftw_forget_wisdom() { radia_nvpl_fft::fftwf_forget_wisdom(); }
static inline void fftw_export_wisdom(void (*fn)(char,void*),void* data) { radia_nvpl_fft::fftwf_export_wisdom(fn,data); }
static inline fftw_status fftw_import_wisdom(int (*fn)(void*),void* data) {
    return radia_nvpl_fft::fftwf_import_wisdom(fn,data) ? FFTW_SUCCESS : FFTW_FAILURE;
}
static inline void fftw_export_wisdom_to_file(FILE* f) { radia_nvpl_fft::fftwf_export_wisdom_to_file(f); }
static inline fftw_status fftw_import_wisdom_from_file(FILE* f) {
    return radia_nvpl_fft::fftwf_import_wisdom_from_file(f) ? FFTW_SUCCESS : FFTW_FAILURE;
}
static inline char* fftw_export_wisdom_to_string() { return radia_nvpl_fft::fftwf_export_wisdom_to_string(); }
static inline fftw_status fftw_import_wisdom_from_string(const char* text) {
    return radia_nvpl_fft::fftwf_import_wisdom_from_string(text) ? FFTW_SUCCESS : FFTW_FAILURE;
}
#define FFTW_HAS_WISDOM
#define FFTW_HAS_PLAN_SPECIFIC
#define FFTW_HAS_FPRINT_PLAN
#define FFTWND_HAS_PRINT_PLAN
