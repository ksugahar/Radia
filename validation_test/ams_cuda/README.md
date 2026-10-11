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
NGSolve's existing ordered sparse LDL factors on the GPU; row products
run within one warp while triangular dependencies stay ordered. Existing dense coarse inverses and the
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

The first gateway call is measured after MATLAB and `gpuDevice`
initialization; MATLAB startup is reported separately. Available-memory
deltas and live-workspace samples are observations, not allocator peaks.
CPU CG stops on its preconditioned norm and counts the initial residual
check; GPU PCG stops on the Euclidean norm and counts updates. Both undergo
the independent true-residual gate, but iteration counts are not directly
interchangeable.

## RTX 5000 result, 2026-10-11

[Saved numerical evidence](rtx5000_20261011.json) covers six small real
curl-curl plus positive-mass cases, with `beta_zero=False`, 96 to 5,326 DOFs,
both AMS cycles and one/two auxiliary levels. Pure curl-curl, complex
systems, A100 scaling and production model integration remain unverified.

At 5,326 DOFs, on the recorded four-thread CPU comparison:

| AMS cycle | CPU warm median | GPU warm median | GPU with transfers |
| --- | ---: | ---: | ---: |
| 1 | 15.80 ms | 41.14 ms | 46.55 ms |
| 7 | 17.22 ms | 34.23 ms | 49.25 ms |

All six true residual and preconditioner-parity checks pass. The prototype
does **not** beat CPU AMS at these sizes. Parallelizing coarse row products
within a warp reduced the earlier ordered single-thread GPU cost, but
per-call validation/allocation, kernel launches, host scalar synchronization
and small triangular dependencies still need measurement and optimization
before claiming fast GPU AMS. A100 must be tested directly on the correct
coordinated host.

Source of the sparse factor layout: NGSolve's public `SparseCholeskyTM`
accessors in `sparsecholesky.hpp` for 6.2.2607. The GPU kernel evaluates the
standard ordered unit-triangular/diagonal solves from those factors; it
does not recompute or replace NGSolve's factorization.
