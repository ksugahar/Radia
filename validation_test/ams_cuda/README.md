# Experimental AMS CUDA MEX

This experiment applies the existing real FP64 Radia AMS preconditioner on a
CUDA GPU. It is an assembled-operator bridge, not NGSolve's upstream
matrix-free GPU backend. It is not a released MATLAB solver or a claim of
production AMS GPU support.

Python/NGSolve owns the HCurl space, numbering, discrete gradient, assembly
and boundary conditions. The shared `HypreBasedAMS`/`CompactAMG` C++ code
owns interpolation, PMIS coarsening, Galerkin products and coarse
factorization. A validation-only Python extension copies that existing
numeric hierarchy and invokes its actual CPU `Mult` for probe comparisons.
There is no second MATLAB CPU solver implementation.

CUDA MEX owns the FP64 AMS cycle, l1-Jacobi sweeps, residual products,
restriction, prolongation, auxiliary V-cycles and PCG vector operations.
Both real AMS cycle types 1 and 7 are supported. Mixed precision and
non-CompactAMG auxiliary solvers fail explicitly. Coarse sparse solves use
NGSolve's existing ordered sparse LDL factors on the GPU; this initial
small-coarse-level kernel is serial. Existing dense coarse inverses and the
explicit ten-sweep coarse treatment are preserved when supplied by setup.
Coarse factorization and hierarchy creation remain CPU setup work. The
host computes PCG scalar control from cuBLAS GPU dot-product results; no
CPU vector or coarse solve is used as a fallback.

The gateway takes a versioned struct containing GPU arrays, a GPU RHS,
tolerance, maximum iteration count and a logical apply-only flag. It has
no persistent handles. Every call validates arrays and allocates temporary
workspace; timings include that overhead. Reusing the input hierarchy on
the GPU therefore does not mean reusing all gateway workspace. This is a
correctness-first prototype and can be slower than CPU AMS.

Build on 100/T7920 in an isolated checkout with the configured MATLAB
Engine Python environment, MATLAB R2026a and its compatible VS2022 toolchain:

```powershell
./Build.ps1 -AmsCudaOnly -RequireNativeProvenance
python validation_test/ams_cuda/run.py --out C:/temp/ams-cuda-owned-run --maxh 0.5 0.2 0.12
```

The CUDA binary includes sm75, sm80 and compute80 PTX. A different OS needs
its own MEX build. This initial build script targets the confirmed Windows
100 environment. A100 performance is unverified until separately measured
on an idle mdx GPU; no A100 result is inferred from RTX measurements.

Acceptance compares random preconditioner applications with shared C++ AMS,
and the resulting solution with an independent NGSolve `sparsecholesky`
solve of the same manufactured HCurl curl-curl plus mass problem. The
FP64 true residual must be at most 1e-6; this small run requests 1e-9.
Results include iterations, first-call and warmed times, input upload,
transfer-inclusive time, GPU memory observations, versions and hashes.
CPU timings use caller-owned TaskManager. Invalid CSR, CPU input, NaN,
invalid smoother, nonconvergence, zero RHS and clear/reload checks run on
the actual gateway. Only sanitized JSON evidence belongs in Git; generated
MAT inputs, raw logs and binaries stay outside tracked source.

Source of the sparse factor layout: NGSolve's public `SparseCholeskyTM`
accessors in `sparsecholesky.hpp` for 6.2.2607. The GPU kernel evaluates the
standard ordered unit-triangular/diagonal solves from those factors; it
does not recompute or replace NGSolve's factorization.
