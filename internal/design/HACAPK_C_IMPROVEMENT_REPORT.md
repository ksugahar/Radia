# HACApK C Version Improvement Report

**Date**: 2025-12-29
**For**: HACApK Developers (Akihiro Ida, Takeshi Iwashita)
**Reviewed by**: Radia Development Team

## Executive Summary

We reviewed the HACApK C version (v1.3.0) used in Radia and compared it with the original Fortran version (v1.0.1) used in ELF_MAGIC. This report identifies potential improvements to reduce H-matrix construction time.

**Current Performance Gap**:
- Radia C version: `t_hmatrix_build = 0.42s` (tetra maxh=0.20, 627 elements)
- ELF Fortran version: `t_hmatrix_build = 0.165s`
- **Gap: 2.5x slower**

## Key Findings

### 1. OpenMP Parallelization Strategy Difference

**Fortran version** (m_HACApK_base.f90:755-822):
```fortran
!$OMP parallel default(none) ...
 ith = omp_get_thread_num()
 nths=lthr(ith); nthe=lthr(ith1)-1
 do ip=nths,nthe  ! Static thread-to-block assignment
   ...
 enddo
!$omp end parallel
```

**C version** (cHACApK_base.c:964-1052 in ELF, 976-1052 in Radia):
```c
// ELF's C version: Same as Fortran (static assignment)
#pragma omp parallel default(none) ...
{
  ith = omp_get_thread_num();
  nths=lthr[ith]; nthe=lthr[ith1]-1;
  for (ip=nths; ip<=nthe; ip++) { ... }
}

// Radia's C version: Uses dynamic scheduling
#pragma omp parallel for schedule(dynamic, 8) ...
for (ip = 1; ip <= nlf; ip++) { ... }
```

**Observation**: Radia uses `schedule(dynamic, 8)` which should improve load balancing but adds scheduling overhead. The Fortran/ELF versions use static pre-computed thread assignment via `lthr[]` array.

**Recommendation**: Test both approaches. For well-balanced workloads, static scheduling may be faster.

### 2. Memory Allocation Pattern in ACA+ Loop

**Issue**: In `cHACApK_calc_vec()` (line 696-708), temporary memory is allocated/freed inside the critical path:

```c
void cHACApK_calc_vec(...) {
  ...
  if(k==0) return;
  zz = (double *) calloc(k,sizeof(double));  // Allocation in hot path
  ...
  free(zz);  // Free in hot path
}
```

**Recommendation**: Pre-allocate workspace at the start of each thread's work and reuse:
```c
// In cHACApK_acaplus, allocate once:
double *zz_workspace = (double *) calloc(kmax, sizeof(double));
// Pass to cHACApK_calc_vec, reuse across calls
```

### 3. ACA+ vs ACA Algorithm Selection

**Fortran version** uses ACA (param[60]=1):
```fortran
if(param(60)==1)then
  kt=HACApK_aca(...)
else
  print*,'Only ACA is avairable!'
```

**C version** uses ACA+ (param[60]=2):
```c
if(param[60]==2) {
  kt=cHACApK_acaplus(...);
}
```

**Observation**: ACA+ is more robust but computationally more expensive than basic ACA. The performance gap may partially be due to algorithm choice.

**Recommendation**: Implement basic ACA (`cHACApK_aca`) in C for comparison. It may be faster for problems where ACA+ robustness is not needed.

### 4. BLAS Optimization Coverage

**Already optimized** (good):
- `cHACApK_unrm_d()` uses `cblas_dnrm2()` ✓
- `cHACApK_adotsub_dsm()` uses `cblas_dgemv()` ✓
- `cHACApK_extract_col()` uses `cblas_dcopy()` ✓

**Not using BLAS** (potential improvement):
- `cHACApK_maxabsvalloc_d()`: Could use `cblas_idamax()` for finding max absolute value
- `cHACApK_minabsvalloc_d()`: No BLAS equivalent, loop is necessary

**Recommendation**: Replace max/min value search with BLAS:
```c
void cHACApK_maxabsvalloc_d(double *vec, double *maxval, int *loc, int n) {
#ifdef HAVE_LAPACK
  int idx = cblas_idamax(n, vec, 1);  // Returns 0-indexed location
  *loc = idx;
  *maxval = fabs(vec[idx]);
#else
  // fallback loop
#endif
}
```

### 5. Dense Block Fill Parallelization

**Current code** (cHACApK_base.c:1043-1050):
```c
for (int il=0; il<ndl; il++) {
  int ill=il+nstrtl;
  for (int it=0; it<ndt; it++) {
    int itt=it+nstrtt;
    double val = cHACApK_entry_ij(lodl[ill],lodt[itt],i_bemv);
    st_lf[ip]->a1[it+ndt*il] = val;
  }
}
```

**Observation**: The inner loop over `ndt` elements is sequential within each block.

**Recommendation**: For large dense blocks, consider nested parallelism or SIMD:
```c
#pragma omp simd
for (int it=0; it<ndt; it++) {
  ...
}
```

### 6. Cluster Tree Construction

The cluster tree construction in `cHACApK_generate_cbitree()` and `cHACApK_count_blrleaf()` uses recursive function calls. This is memory-efficient but may have overhead.

**Recommendation**: Consider iterative implementation using explicit stack for very deep trees.

## Performance Comparison Summary

| Component | Radia C | ELF Fortran | Ratio | Notes |
|-----------|---------|-------------|-------|-------|
| H-matrix build | 0.42s | 0.165s | 2.5x | Algorithm + parallelization |
| Linear solve | 0.23s | 0.23s | 1.0x | Similar performance |
| Matrix build | 0.16s | 0.16s | 1.0x | Similar performance |

## Recommended Priority

1. **High Priority**: Test static vs dynamic OpenMP scheduling
2. **High Priority**: Pre-allocate workspace in ACA+ to avoid malloc/free in hot path
3. **Medium Priority**: Implement basic ACA for comparison
4. **Medium Priority**: Use `cblas_idamax()` for max value search
5. **Low Priority**: SIMD optimization for dense block fill

## Files Modified in Radia

Radia's HACApK C version includes the following enhancements over ELF's C version:
- [cHACApK_base.c](../../src/ext/HACApK/cHACApK_base.c): Dynamic scheduling, BLAS optimization
- [cHACApK_lib.c](../../src/ext/HACApK/cHACApK_lib.c): MKL BLAS integration

## Conclusion

The 2.5x performance gap between Radia C and ELF Fortran is primarily due to:
1. ACA+ vs ACA algorithm difference
2. Dynamic vs static OpenMP scheduling overhead
3. Repeated memory allocation in hot paths

Addressing items 1-2 in the priority list should significantly close the performance gap.
