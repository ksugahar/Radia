# NGSolve Poisson / CUDA MEX acceptance

CPU solver execution belongs to Python/NGSolve; MATLAB solver execution uses
CUDA MEX. MATLAB is used for GPU execution independently of Simulink.
This first bridge solves **NGSolve-assembled CSR systems** with FP64 Jacobi
preconditioned CG. It does not implement finite-element bases, upstream
`ngscuda` / matrix-free operators, ICCG, AMS, or HDiv-MMM.

The manufactured problem is `-Delta u = 3*pi^2*u` on the unit cube,
`u = sin(pi*x)*sin(pi*y)*sin(pi*z)`, with homogeneous Dirichlet data and H1
order two. NGSolve generates the mesh, assembles the form, and eliminates
Dirichlet unknowns. The exact same reduced matrix and RHS go to CUDA MEX.
Current NGSolve preprocessing and CPU references run serially without
TaskManager. Historical JSON preserves its original threading configuration.
Python NGSolve sparsecholesky is the independent discrete reference; NGSolve
Jacobi CG is the iterative CPU timing comparison. No CPU MEX is built.

## Build on 100

Use an isolated Python 3.12 environment with the NGSolve/Netgen pins from
package metadata, NumPy, SciPy, psutil and the official MATLAB Engine for the installed
MATLAB. For R2026a, install Engine from `matlabroot/extern/engines/python`.
R2026a supplies CUDA 12.8. Use its supported Visual Studio 2022 C++ toolchain;
do not bypass CUDA's unsupported-compiler checks to use VS 2026. Both Visual
Studio installations can coexist. Parallel Computing Toolbox is required.

```powershell
# Put the isolated Python Scripts directory first on PATH for this invocation.
./Build.ps1 -CudaMexOnly
python validation_test/ngsolve_cuda_mex/run_poisson.py `
  --binary-dir ./build-cuda-mex --out "$env:TEMP/poisson-cuda-acceptance" --threads 1
```

The build stages `radia_cuda_csr_pcg.mexw64` and `build.json` in
`build-cuda-mex`; it does not replace an installed or loaded gateway. Targets
are native sm_75 (RTX 5000), sm_80 (A100), and compute_80 PTX. An A100 target
in a binary is not evidence of execution on A100.

The MEX accepts resident real `gpuArray` int32 zero-based CSR indices and
double values/RHS/positive Jacobi diagonal inverse. It returns a GPU solution
and diagnostic struct, and rejects CPU inputs, invalid CSR, nonfinite values,
CG breakdown or failed true residual. It is an SPD-only experimental entry
point, not a general sparse solver. All temporary GPU allocations are scoped
to a call. There is no hidden Python call or CPU solver fallback.

## Evidence and limits

Every run saves JSON: versions, source/binary identity, GPU, CPU thread count,
assembly/setup, initial GPU solve, five warmed solves, upload/download cost,
memory, CG iteration counts, analytic L2 error, true residual and error
against sparsecholesky. Gateway checks include invalid inputs, zero RHS,
nonconvergence, repeated calls and clear/reload.

GPU timing schema `radia.cuda-poisson-gpu.v2` measures five complete warmed
upload/solve/download/release calls. `warm_with_transfers_seconds` is their
median; `warm_full_transfer_seconds` retains the individual samples.
`estimated_first_transfer_plus_warm_kernel_seconds` retains the earlier
sum of first upload, resident-kernel median and download as an estimate.
Historical v1 records keep their original values and must not be interpreted
as measured end-to-end timings. Neither timing includes FE assembly or disk I/O.

CPU/GPU CG use different stopping norms; both must pass an independently
recomputed true residual <= 1e-6, and the GPU solution must agree with the
direct reference to relative vector L2 <= 1e-8. The direct factorization and
reuse times are reported separately. Small 100 acceptance cases establish
operation, not peak GPU speed or A100 scaling; unfavorable speed results must
remain visible. Large measurements belong on the coordinated gpu evaluation worker/A100 host.

Official references:
- [MathWorks CUDA MEX](https://www.mathworks.com/help/parallel-computing/run-mex-functions-containing-cuda-code.html)
- [NGSolve matrix-free GPU operators](https://docu.ngsolve.org/latest/i-tutorials/unit-5.6-gpu/matrixfree.html)
