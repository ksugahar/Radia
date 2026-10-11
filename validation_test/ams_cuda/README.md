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

Build on the designated Windows build host in an isolated checkout with the configured MATLAB
Engine Python environment, MATLAB R2026a and its compatible VS2022 toolchain:

```powershell
./Build.ps1 -AmsCudaOnly -RequireNativeProvenance
python validation_test/ams_cuda/run.py --out "$env:TEMP/ams-cuda-owned-run" --maxh 0.5 0.2 0.12
```

The CUDA binary includes sm75, sm80 and compute80 PTX. A different OS needs
its own MEX build. This initial build script targets the confirmed Windows
build environment. A100 performance is unverified until separately measured
on an idle compute GPU; no A100 result is inferred from RTX measurements.

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

The launch-reduction candidate fuses CSR Jacobi products and updates with
separate old/new buffers, fuses prolongation and addition, and combines PCG
solution/residual updates. GPU-to-GPU copies are queued on the same default
stream as kernels and cuBLAS. Device-to-device cudaMemcpy is already host
asynchronous according to NVIDIA's API synchronization contract; replacing it
with cudaMemcpyAsync alone is not evidence of fewer synchronization points.
PCG dot products now return device scalars; rho and pAp are read together once
per iteration, while next-rho and beta stay on device. This removes the separate
host wait for next-rho. Host scalar convergence and breakdown checks and the final FP64 true-residual
gate remain in place. No CPU solve fallback or CUDA Graph is introduced.

Before uploading, the benchmark estimates numeric inputs, AMS/PCG workspace,
duplicate uploads and diagnostic arrays, plus library headroom. It rejects
cases that would consume the available GPU memory minus the larger of 2 GiB
or 20% of total device memory. The estimate is conservative admission control,
not a measured allocator peak. Develop within GPU memory; increase case sizes
gradually. A 10 GB matrix alone is not a 10 GB solver working set.

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

## Residual and validation overhead experiment

The next implementation combines CSR products with residual output, combines
the smoother's residual subtraction and update, and starts a zero-initialized
Jacobi sweep directly with `b/l1`. Smoother updates still read a completed
matrix product: updating `x` inside the CSR kernel would introduce row races.
Independent input checks share a device error flag. Coarse row-pointer bounds
are checked before factor validation dereferences them, and all pending
checks complete before any solve. The hierarchy, factors, FP64 precision,
iteration rules and acceptance thresholds are unchanged.

`compare.py` compares two hash-verified MEX builds on the **same saved MAT
inputs**, with before/after/after/before order repeated three times. Each
block discards its first call and measures five warmed solves, giving 30
samples per build and case. It also records transfer-inclusive samples,
CPU-preconditioner probe errors, true residuals and gateway checks in JSON.
Run it only when competing CI and numerical work have finished:

```powershell
python validation_test/ams_cuda/compare.py --before $beforeBuild --after $afterBuild --inputs $inputCycle1 $inputCycle7 --out $comparisonOutput
```

Run `run.py` separately for independent direct-solution checks and CPU timings.
MEX load/clear between blocks and MATLAB's cached GPU allocator can affect
individual samples; retain all samples and compare medians rather than a
single fastest run. These are small acceptance measurements, not GPU scaling
evidence.

[Saved overhead experiment](rtx5000_fusion_20261011.json) records all samples
and six numerical-acceptance cases. The paired measurements use identical
5,326-DOF inputs with 30 warmed samples per build and cycle:

| AMS cycle | Before warm median | After warm median | Reduction | Before/after with transfers |
| --- | ---: | ---: | ---: | ---: |
| 1 | 37.56 ms | 31.17 ms | 17.0% | 43.84 / 37.70 ms |
| 7 | 33.45 ms | 27.49 ms | 17.8% | 39.73 / 34.11 ms |

This is a **shared-host relative comparison**: unrelated CPU CI environment
preparation/build was active, and the measurement process used BelowNormal
priority. The JSON records the observed CI steps. It is not an isolated-host
benchmark or evidence of beating CPU AMS; do not combine its samples with
the earlier CPU table to claim a CPU speedup. Preliminary switching-check
samples are excluded. A100 was neither started nor used.

All six numerical cases pass: maximum preconditioner probe relative error
is 3.45e-14, maximum FP64 true relative residual is 9.36e-10, and maximum
relative solution difference from independent CPU sparsecholesky is 4.76e-9.
All eleven gateway checks pass. The numerical-only run's timing samples
are excluded because it overlapped CPU CI as well.

## Idle-host repeat, 2026-10-11

[Saved idle-host evidence](rtx5000_idle_20261011.json) repeats the same six
small cases after repository CI finished. No CI was active or queued at the
start, midpoint or end. Two process/load probes five seconds apart preceded
measurement; total CPU utilization was 5% before and 9% afterward, with 2%
GPU utilization at both boundaries. Existing idle desktop/MATLAB sessions
were retained, so this is operational idle rather than exclusive reservation.
A100 was not started or used.

CPU uses Python/pybind11 + NGSolve with four threads; GPU uses the unchanged
FP64 mexcuda binary. These are the five-sample warmed solve medians from the
same acceptance run, excluding hierarchy setup for both routes. GPU gateway
validation/allocation is included; the last column also includes transfers.

| AMS cycle, 5,326 DOFs | CPU | GPU resident input | GPU with transfers |
| --- | ---: | ---: | ---: |
| 1 | 20.44 ms | 32.21 ms | 42.41 ms |
| 7 | 17.76 ms | 27.40 ms | 33.36 ms |

The GPU is still 1.58x / 1.54x slower than CPU for these cases. All six CPU
and GPU true residuals are below 1e-9; maximum CPU-preconditioner probe error
is 3.29e-14 and maximum solution difference from independent sparsecholesky
is 4.76e-9. All eleven gateway checks pass. Stopping-norm and iteration-count
semantics remain as described above.

A separate identical-input ABBA comparison retains 30 warmed samples per
build/cycle: cycle 1 improves from 39.01 to 32.52 ms (16.6%), and cycle 7
from 34.77 to 28.60 ms (17.7%). These medians belong to that separate
comparison, not the five-sample CPU/GPU table.

Zero-RHS calls skip PCG iterations but retain gateway validation/allocation,
cuBLAS setup and the RHS norm. Their paired medians drop from about 7.6 ms
to 3.2 ms. Apply-only calls, including one AMS cycle and their own setup,
take about 3.5--3.6 ms after the change. These wall-time diagnostics are not
exact kernel/allocation counters and use different workspace paths. Do not
subtract zero-RHS time from apply-only time. Gateway reuse alone therefore
cannot explain or eliminate the full CPU/GPU gap; iterative launch/sync and
coarse-solve costs need further profiling before selecting an optimization.

Source of the sparse factor layout: NGSolve's public `SparseCholeskyTM`
accessors in `sparsecholesky.hpp` for 6.2.2607. The GPU kernel evaluates the
standard ordered unit-triangular/diagonal solves from those factors; it
does not recompute or replace NGSolve's factorization.
