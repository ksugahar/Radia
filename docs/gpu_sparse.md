# Optional Python sparse CUDA operations

`radia.gpu` provides explicit CPU/GPU selection for real SPD systems. Importing
the module does not import CUDA libraries. Optional GPU modules are built
separately with CUDA 12.8 and use the exact pip runtime versions in Radia's
`gpu` dependency extra. A CPU or GPU binary from a different source revision
is rejected. CUDA failures propagate; a GPU request never selects a CPU solve.

The current implementation is experimental. The recorded GB200 MIG result
covers the standalone two-case CSR Jacobi-PCG component only. It does not
qualify this public package entry, the new AMG/AMS bridge, full-GPU execution,
all Radia GPU features, performance, or release readiness. Those stages need
their own built-artifact and runtime evidence.

## Explicit SPD solve

NGSolve owns matrix assembly, orientation, numbering and the caller's normal
TaskManager policy. Given an assembled real matrix `a` and NGSolve vector `b`:

```python
from radia.gpu import solve_spd

x, info = solve_spd(a, b, backend="gpu", tolerance=1e-12, max_iterations=1000)
```

The GPU route uploads the CSR operator and uses the shared FP64 Jacobi-PCG
kernel. After the solve, Python recomputes the original residual through the
actual NGSolve matrix operator, rather than trusting exported CSR storage or
native residual metadata. A symmetric triangular storage export cannot certify
an incorrect solution to the full operator.

`backend="cpu"` explicitly selects the independent NGSolve `sparsecholesky`
route. It does not load a CUDA module. This is a separate requested operation,
not a fallback or an accuracy test that GPU users must run before solving.

## Existing CompactAMG and AMS hierarchies

Construct a real CPU hierarchy through the existing Radia/NGSolve setup workflow
first. The GPU adapter consumes an owned algebraic snapshot; it does not create
a finite-element space, replace a physical formulation or take over TaskManager.

```python
from radia.gpu import AMG, AMS

with AMG(amg_hierarchy, backend="gpu") as solver:
    x, info = solver.solve(rhs_array, tolerance=1e-12, max_iterations=1000)

with AMS(ams_hierarchy, backend="gpu") as solver:
    applied, metadata = solver.apply(rhs_array)
    x, info = solver.solve(rhs_array, tolerance=1e-12, max_iterations=1000)
```

`apply` performs one preconditioner cycle and reports
`operation="preconditioner_apply", residual_checked=False`. It is not a solved
system. `solve` runs GPU PCG, checks the FP64 true residual, and independently
recomputes the snapshot-system residual on the host. Numerical arrays require
finite float64 values and exact int32 indices; no precision conversion is
performed by the native snapshot parser.

The AMG adapter preserves the existing CompactAMG smoothing counts, transfers
and coarse action. The AMS adapter supports existing real FP64 cycles 1 and 7,
including their existing auxiliary AMG hierarchies. Mixed-precision and complex
AMS are unsupported. Constrained RHS entries must be zero for a solve. Matrix
and hierarchy snapshots are frozen: passing a different operator to
`solve_spd(..., preconditioner=solver)` is rejected. An update needs a new owner.

Both adapters own their uploaded hierarchy until `close` or context exit.
`close` is idempotent; operations after close fail. Native operations serialize
against close without holding the Python GIL. A context belongs to the CUDA
device selected at construction and refuses execution after a device change.

The separate `ResidentAMG(levels, coarse_inverse, backend="gpu")` interface
retains the existing caller-supplied algebraic hierarchy contract and its native
cycle. It is distinct from uploading an actual CompactAMG snapshot; equality
between those two hierarchy constructions is not assumed.

## Shared source and build boundary

`src/core/gpu` is the single owner of shared CUDA numerical bodies. Radia's
Python modules and the Windows Simulink MEX gateways consume that source;
normal Radia code does not import `radia_simulink` or MATLAB. MEX-specific array
parsing and handles remain in the Simulink package. HDiv-MMM is outside this
sparse GPU implementation.

`tools/build_python_gpu.py` builds `_gpu_csr`, `_gpu_resident_amg` and `_gpu_ams`
from a clean exact-source checkout into a job-owned directory outside source.
It requires an exact-source build admission and does not install packages,
stage binaries into a running runtime, execute GPU kernels or dispatch CI.
Windows uses the supported MSVC/MATLAB-bundled toolchain; Linux uses explicit
CUDA 12.8 and host compiler paths. Resource admission and runtime qualification
remain separate from source preparation.

Public companion `.build.json` receipts contain source/artifact identity and
a hash of the private operational receipt. Raw commands, admission records and
absolute workstation/toolchain paths stay in private recovery. Pip libraries
are verified against their SHA256 RECORD entries and actual loaded DLL/ELF
paths before and after import and on reuse. No toolkit runtime is bundled.
