# Sparse Python GPU surface audit — 2026-10-28

This targeted audit covers `radia.gpu` and its optional sparse CUDA boundary.
It is not a new inventory or certification of the full Radia package.

| Family | Bucket | Decision |
|---|---|---|
| Explicit `solve_spd(..., backend=...)` workflow | user-intent | Keep one entry. Delegate CPU factorization to NGSolve sparsecholesky; use the shared CUDA kernel only after explicit selection. |
| Existing CompactAMG/AMS CUDA cycles | method-keep | Preserve the existing Radia hierarchy/cycle mathematics. NGSolve owns FE construction, ordering and CPU hierarchy setup. |
| Caller-supplied `ResidentAMG` algebraic hierarchy | method-demote | Keep the existing expert/reproducibility route in `radia.gpu`; do not re-export it from the package root or treat it as an application UI. |
| `_accelerator_snapshot`, native `_gpu_*` modules and runtime admission | method-demote | Private implementation boundary. Do not promote generic array/device helpers into public API. |
| Mesh/CAD/FE assembly and generic CPU sparse algebra | plumbing-delete | Introduce no competing API here. Delegate to NGSolve/Netgen and existing CAD/mesh owners. |

Counts for these five families: user-intent 1, method-keep 1, method-demote 2,
plumbing-delete 1, deprecated-drop 0.

## Safe in this candidate

- Keep CUDA imports lazy and backend choice explicit; reject a failed GPU
  request instead of selecting a CPU substitute.
- Add `AMG`/`AMS` owners around the existing CPU hierarchy snapshots and shared
  CUDA cycles, with explicit apply-versus-solve and lifetime contracts.
- Preserve one canonical numerical source under `src/core/gpu` for Python and
  the Windows Simulink MEX consumers. Keep gateway parsing in its owner.
- Recompute solve residuals through the actual original NGSolve operator, so
  triangular CSR storage cannot falsely certify a full-system solution.

## Deferred actions and qualification

Do not add top-level package re-exports, generic GPU array arithmetic, CAD/mesh
wrappers, compatibility aliases or a second CPU solver. Broader application
intent interfaces require the corresponding physical workflow to be ready;
HDiv-MMM is explicitly outside this candidate.

Source contracts and mocked loader tests do not qualify new native artifacts.
The public package entry and AMG/AMS bridge need exact-source native builds,
independent CPU-reference comparisons and recorded CUDA runtime/device evidence.
The existing standalone GB200 MIG CSR component result does not establish
full-GPU, all-package GPU, performance or release acceptance.
