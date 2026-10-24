#pragma once
/* Build-selected LP64 provider. There is no runtime provider selection.
 * This C-compatible boundary is also used by the shared HACApK C sources. */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#if defined(RADIA_BLAS_MKL) == defined(RADIA_BLAS_NVPL)
#error "Select exactly one Radia BLAS provider"
#endif
#if defined(MKL_ILP64) || defined(NVPL_ILP64)
#error "Radia requires LP64 BLAS/LAPACK"
#endif

#if defined(RADIA_BLAS_MKL)
#if !(defined(__x86_64__) || defined(_M_X64))
#error "MKL is only admitted on x86_64"
#endif
#include <mkl.h>
typedef MKL_Complex16 radia_blas_complex_double;
#define RADIA_BLAS_NAME "mkl"
#define RADIA_BLAS_THREAD_POLICY "mkl-caller-controlled"
static inline int radia_blas_get_max_threads(void) { return mkl_get_max_threads(); }
static inline void radia_blas_set_num_threads(int n) { mkl_set_num_threads(n); }
static inline int radia_blas_set_num_threads_local(int n) { return mkl_set_num_threads_local(n); }
#else
#if !defined(__aarch64__) || !defined(__linux__)
#error "NVPL is only admitted on Linux aarch64"
#endif
/* Use the documented two-real struct representation for the Fortran ABI,
 * including C++ translation units (never compiler-specific C99 _Complex). */
#define NVPL_LAPACK_COMPLEX_STRUCTURE
#include <nvpl_lapacke.h>
#include <nvpl_blas_cblas.h>
typedef nvpl_dcomplex_t radia_blas_complex_double;
#define RADIA_BLAS_NAME "nvpl"
#define RADIA_BLAS_THREAD_POLICY "lp64-sequential"
static inline int radia_blas_get_max_threads(void) { return 1; }
static inline void radia_blas_set_num_threads(int n) {
    if(n != 1) { fprintf(stderr, "Radia NVPL requires sequential BLAS (threads=1)\n"); abort(); }
}
static inline int radia_blas_set_num_threads_local(int n) {
    if(n != 0 && n != 1) { fprintf(stderr, "Radia NVPL requires sequential BLAS (threads=1)\n"); abort(); }
    return 1; /* restoring 1 or clearing a local override both retain seq */
}
#endif

typedef char radia_lp64_integer_check[(sizeof(lapack_int)==4) ? 1 : -1];
typedef char radia_complex_pair_check[
    (sizeof(radia_blas_complex_double)==2*sizeof(double) &&
     offsetof(radia_blas_complex_double,real)==0 &&
     offsetof(radia_blas_complex_double,imag)==sizeof(double)) ? 1 : -1];

#ifdef __cplusplus
#include <complex>
#include <stdexcept>
#include <type_traits>
static_assert(sizeof(lapack_int) == 4, "LP64 LAPACK integer ABI required");
static_assert(sizeof(radia_blas_complex_double) == 2*sizeof(double), "Complex pair ABI required");
static_assert(sizeof(std::complex<double>) == sizeof(radia_blas_complex_double), "Complex layout size");
static_assert(alignof(std::complex<double>) == alignof(radia_blas_complex_double), "Complex layout alignment");
static_assert(std::is_trivially_copyable<std::complex<double>>::value, "Complex byte representation");
namespace radia {
class DenseThreadGuard {
    int saved_;
public:
    explicit DenseThreadGuard(int n) : saved_(radia_blas_get_max_threads()) {
#ifdef RADIA_BLAS_NVPL
        // TaskManager may own several workers; each uses the selected seq BLAS.
        if(n < 1) throw std::invalid_argument("Positive caller thread budget required");
        radia_blas_set_num_threads(1);
#else
        radia_blas_set_num_threads(n);
#endif
    }
    ~DenseThreadGuard() { radia_blas_set_num_threads(saved_); }
    DenseThreadGuard(const DenseThreadGuard&) = delete;
    DenseThreadGuard& operator=(const DenseThreadGuard&) = delete;
};
}
#endif
