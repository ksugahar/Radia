#ifndef exafmm_t_h
#define exafmm_t_h

// Prevent Windows min/max macro conflicts
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <set>
#include <vector>
#include "align.h"
#include "vec.h"

// Intel MKL DFTI for FFT operations (replacing FFTW3)
#include <mkl_dfti.h>

// M_PI may not be defined on MSVC without _USE_MATH_DEFINES
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// FFTW3-compatible constants
#define FFTW_FORWARD (-1)
#define FFTW_BACKWARD (+1)
#define FFTW_ESTIMATE (1U << 6)
#define FFTW_MEASURE (0U)

// Parallelization: NGSolve TaskManager (OpenMP removed)

namespace exafmm_t {
  const int MEM_ALIGN = 64;
  const int CACHE_SIZE = 512;
  const int NCHILD = 8;

  // Always use double precision (single precision dropped per user request)
  typedef double real_t;                       //!< Real number type
  const real_t EPS = 1e-16;

  // MKL DFTI types for FFT operations
  typedef DFTI_DESCRIPTOR_HANDLE fft_plan;

  // Complex type compatible with MKL (interleaved real/imag)
  typedef struct { double real; double imag; } fft_complex;

  const real_t PI = M_PI;
  typedef std::complex<real_t> complex_t;       //!< Complex number type
  typedef vec<3,int> ivec3;                     //!< Vector of 3 int types
  typedef vec<3,real_t> vec3;                   //!< Vector of 3 real_t types
  typedef vec<3,complex_t> cvec3;               //!< Vector of 3 complex_t types

  //! SIMD vector types for AVX512, AVX, and SSE
  const int NSIMD = SIMD_BYTES / int(sizeof(real_t));  //!< SIMD vector length (SIMD_BYTES defined in vec.h)
  typedef vec<NSIMD, real_t> simdvec;           //!< SIMD vector of NSIMD real_t types

  typedef std::vector<real_t> RealVec;          //!< Vector of real_t types
  typedef std::vector<complex_t> ComplexVec;    //!< Vector of complex_t types
  typedef AlignedAllocator<real_t, MEM_ALIGN> AlignAllocator;   //!< Allocator for memory alignment
  typedef std::vector<real_t, AlignAllocator> AlignedVec;       //!< Aligned vector of real_t types

  //! Interaction types that need to be pre-computed.
  typedef enum {
    M2M_Type = 0,
    L2L_Type = 1,
    M2L_Helper_Type = 2,
    M2L_Type = 3,
    Type_Count = 4
  } Precompute_Type;

  /**
   * @brief Structure of bodies.
   *
   * @tparam T Value type of sources and targets (real or complex).
   */
  template <typename T>
  struct Body {
    int ibody;                             //!< Initial body numbering for sorting back
    vec3 X;                                //!< Coordinates
    T q;                                   //!< Charge
    T p;                                   //!< Potential
    vec<3,T> F;                            //!< Gradient
  };
  template <typename T> using Bodies = std::vector<Body<T>>;     //!< Vector of nodes

  /**
   * @brief Structure of nodes.
   *
   * @tparam Value type of sources and targets (real or complex).
   */
  template <typename T>
  struct Node {
    size_t idx;                                 //!< Index in the octree
    size_t idx_M2L;                             //!< Index in global M2L interaction list
    bool is_leaf;                               //!< Whether the node is leaf
    int ntrgs;                                  //!< Number of targets
    int nsrcs;                                  //!< Number of sources
    vec3 x;                                     //!< Coordinates of the center of the node
    real_t r;                                   //!< Radius of the node
    uint64_t key;                               //!< Morton key
    int level;                                  //!< Level in the octree
    int octant;                                 //!< Octant
    Node* parent;                               //!< Pointer to parent
    std::vector<Node*> children;                //!< Vector of pointers to child nodes
    std::vector<Node*> P2L_list;                //!< Vector of pointers to nodes in P2L interaction list
    std::vector<Node*> M2P_list;                //!< Vector of pointers to nodes in M2P interaction list
    std::vector<Node*> P2P_list;                //!< Vector of pointers to nodes in P2P interaction list
    std::vector<Node*> M2L_list;                //!< Vector of pointers to nodes in M2L interaction list
    std::vector<int> isrcs;                     //!< Vector of initial source numbering
    std::vector<int> itrgs;                     //!< Vector of initial target numbering
    RealVec src_coord;                          //!< Vector of coordinates of sources in the node
    RealVec trg_coord;                          //!< Vector of coordinates of targets in the node
    std::vector<T> src_value;                   //!< Vector of charges of sources in the node
    std::vector<T> trg_value;                   //!< Vector of potentials and gradients of targets in the node
    std::vector<T> up_equiv;                    //!< Upward check potentials / Upward equivalent densities
    std::vector<T> dn_equiv;                    //!< Downward check potentials / Downward equivalent densites
  };

  // alias template
  template <typename T> using Nodes = std::vector<Node<T>>;        //!< Vector of nodes
  template <typename T> using NodePtrs = std::vector<Node<T>*>;    //!< Vector of Node pointers
  using Keys = std::vector<std::set<uint64_t>>;                    //!< Vector of Morton keys of each level

  //! M2L setup data
  struct M2LData {
    std::vector<size_t> fft_offset;   // source's first child's upward_equiv's displacement
    std::vector<size_t> ifft_offset;  // target's first child's dnward_equiv's displacement
    RealVec ifft_scale;
    std::vector<size_t> interaction_offset_f;
    std::vector<size_t> interaction_count_offset;
  };

  // Relative coordinates and interaction lists (C++17 inline variables)
  inline std::vector<std::vector<ivec3>> REL_COORD;  //!< Vector of possible relative coordinates (inner) of each interaction type (outer)
  inline std::vector<std::vector<int>> HASH_LUT;     //!< Vector of hash Lookup tables (inner) of relative positions for each interaction type (outer)
  inline std::vector<std::vector<int>> M2L_INDEX_MAP;  //!< [M2L_relpos_idx][octant] -> M2L_Helper_relpos_idx

  // ========================================================================
  // MKL DFTI FFT wrapper functions (replacing FFTW3 calls)
  // Direct MKL integration - no FFTW3 wrapper dependency
  // ========================================================================

  // Complex-to-Complex FFT plan
  inline fft_plan fft_plan_dft(int rank, const int* n, fft_complex* in, fft_complex* out, int sign, unsigned flags) {
    (void)sign;   // MKL handles direction in execute call
    (void)flags;  // MKL doesn't use planning flags
    DFTI_DESCRIPTOR_HANDLE handle = nullptr;
    MKL_LONG dims[3];
    for (int i = 0; i < rank; i++) dims[i] = n[i];

    DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, rank, dims);
    if (in == out) {
      DftiSetValue(handle, DFTI_PLACEMENT, DFTI_INPLACE);
    } else {
      DftiSetValue(handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    }
    DftiCommitDescriptor(handle);
    return handle;
  }

  // Real-to-Complex FFT plan (R2C)
  inline fft_plan fft_plan_dft_r2c(int rank, const int* n, real_t* in, fft_complex* out, unsigned flags) {
    (void)in;
    (void)out;
    (void)flags;
    DFTI_DESCRIPTOR_HANDLE handle = nullptr;
    MKL_LONG dims[3];
    for (int i = 0; i < rank; i++) dims[i] = n[i];

    DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_REAL, rank, dims);
    DftiSetValue(handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    DftiSetValue(handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    DftiCommitDescriptor(handle);
    return handle;
  }

  // Batched Complex-to-Complex FFT plan
  inline fft_plan fft_plan_many_dft(int rank, const int* n, int howmany,
                                     fft_complex* in, const int* inembed, int istride, int idist,
                                     fft_complex* out, const int* onembed, int ostride, int odist,
                                     int sign, unsigned flags) {
    (void)in; (void)out;
    (void)inembed; (void)onembed;
    (void)istride; (void)ostride;
    (void)sign; (void)flags;

    DFTI_DESCRIPTOR_HANDLE handle = nullptr;
    MKL_LONG dims[3];
    for (int i = 0; i < rank; i++) dims[i] = n[i];

    DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, rank, dims);
    DftiSetValue(handle, DFTI_NUMBER_OF_TRANSFORMS, (MKL_LONG)howmany);
    DftiSetValue(handle, DFTI_INPUT_DISTANCE, (MKL_LONG)idist);
    DftiSetValue(handle, DFTI_OUTPUT_DISTANCE, (MKL_LONG)odist);
    DftiSetValue(handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    DftiCommitDescriptor(handle);
    return handle;
  }

  // Batched Real-to-Complex FFT plan (R2C)
  inline fft_plan fft_plan_many_dft_r2c(int rank, const int* n, int howmany,
                                         real_t* in, const int* inembed, int istride, int idist,
                                         fft_complex* out, const int* onembed, int ostride, int odist,
                                         unsigned flags) {
    (void)in; (void)out;
    (void)inembed; (void)onembed;
    (void)istride; (void)ostride;
    (void)flags;

    DFTI_DESCRIPTOR_HANDLE handle = nullptr;
    MKL_LONG dims[3];
    for (int i = 0; i < rank; i++) dims[i] = n[i];

    DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_REAL, rank, dims);
    DftiSetValue(handle, DFTI_NUMBER_OF_TRANSFORMS, (MKL_LONG)howmany);
    DftiSetValue(handle, DFTI_INPUT_DISTANCE, (MKL_LONG)idist);
    DftiSetValue(handle, DFTI_OUTPUT_DISTANCE, (MKL_LONG)odist);
    DftiSetValue(handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    DftiSetValue(handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    DftiCommitDescriptor(handle);
    return handle;
  }

  // Batched Complex-to-Real FFT plan (C2R)
  inline fft_plan fft_plan_many_dft_c2r(int rank, const int* n, int howmany,
                                         fft_complex* in, const int* inembed, int istride, int idist,
                                         real_t* out, const int* onembed, int ostride, int odist,
                                         unsigned flags) {
    (void)in; (void)out;
    (void)inembed; (void)onembed;
    (void)istride; (void)ostride;
    (void)flags;

    DFTI_DESCRIPTOR_HANDLE handle = nullptr;
    MKL_LONG dims[3];
    for (int i = 0; i < rank; i++) dims[i] = n[i];

    DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_REAL, rank, dims);
    DftiSetValue(handle, DFTI_NUMBER_OF_TRANSFORMS, (MKL_LONG)howmany);
    DftiSetValue(handle, DFTI_INPUT_DISTANCE, (MKL_LONG)idist);
    DftiSetValue(handle, DFTI_OUTPUT_DISTANCE, (MKL_LONG)odist);
    DftiSetValue(handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    DftiSetValue(handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    DftiCommitDescriptor(handle);
    return handle;
  }

  // Execute Complex-to-Complex FFT (forward)
  inline void fft_execute_dft(fft_plan plan, fft_complex* in, fft_complex* out) {
    if (in == out) {
      DftiComputeForward(plan, in);
    } else {
      DftiComputeForward(plan, in, out);
    }
  }

  // Execute Real-to-Complex FFT (forward)
  inline void fft_execute_dft_r2c(fft_plan plan, real_t* in, fft_complex* out) {
    DftiComputeForward(plan, in, out);
  }

  // Execute Complex-to-Real FFT (backward/inverse)
  inline void fft_execute_dft_c2r(fft_plan plan, fft_complex* in, real_t* out) {
    DftiComputeBackward(plan, in, out);
  }

  // Destroy FFT plan
  inline void fft_destroy_plan(fft_plan plan) {
    if (plan) DftiFreeDescriptor(&plan);
  }

  // Get FFT flop count (MKL doesn't provide this, return 0)
  inline void fft_flops(fft_plan plan, double* add, double* mul, double* fma) {
    (void)plan;
    *add = 0;
    *mul = 0;
    *fma = 0;
  }

}  // namespace exafmm_t

#endif  // exafmm_t_h
