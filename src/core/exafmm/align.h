#ifndef align_h
#define align_h
#include <cstdlib>
#include <memory>

// Windows compatibility for aligned memory allocation
#ifdef _WIN32
#include <malloc.h>  // For _aligned_malloc, _aligned_free
#endif

template <typename T, size_t NALIGN>
struct AlignedAllocator : public std::allocator<T> {
  template <typename U>
  struct rebind {
    typedef AlignedAllocator<U, NALIGN> other;
  };

  T * allocate(size_t n) {
    void *ptr = nullptr;
#ifdef _WIN32
    // Windows: use _aligned_malloc
    ptr = _aligned_malloc(n * sizeof(T), NALIGN);
    if (ptr == nullptr) throw std::bad_alloc();
#else
    // POSIX: use posix_memalign
    int rc = posix_memalign(&ptr, NALIGN, n * sizeof(T));
    if (rc != 0) return nullptr;
    if (ptr == nullptr) throw std::bad_alloc();
#endif
    return reinterpret_cast<T*>(ptr);
  }

  void deallocate(T * p, size_t) {
#ifdef _WIN32
    // Windows: use _aligned_free
    _aligned_free(p);
#else
    // POSIX: use free
    free(p);
#endif
  }
};
#endif
