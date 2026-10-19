# MATLAB/Simulink distribution and explicit backend migration

Status: breaking source/metadata/import/CLI/CI-reference rename implemented in
an isolated checkout. Not published or installed-wheel accepted. MATLAB CPU
execution and direct Python GPU numerical bindings remain incomplete.
The canonical PyPI distribution name is `radia-simulink`, and its Python
package is `radia_simulink`. It owns standalone MATLAB, Simulink, LTspice
and direct Python GPU integration with explicit CPU/GPU execution. CUDA remains
an acceleration backend rather than the definition of the entire interface.

## Verified present state

- `packages/radia-simulink/pyproject.toml` names the distribution `radia-simulink`.
- Its README assigns CPU electromagnetic execution to Python `radia`, excludes
  CPU solver fallback, and explicitly excludes complete migration of older CPU
  MEX application models from the published feature claims.
- `matlab/python_api_parity_manifest.json` marks CPU MEX unsupported and does
  not promise MATLAB equivalents for all Python APIs.
- Native projection uses NGSolve `SetValues` on the mexcuda host. This proves
  that CPU host FE operations can exist in the MATLAB native bridge; it does
  not prove that all legacy Mesh/space/form commands or CPU solvers exist.
- RTX numerical acceptance and A100 measurements apply to their recorded
  binaries and cases, not arbitrary installed MATLAB paths or every API.
- A legacy `radia:mex:Retired` failure is not evidence that every CUDA entry
  point fails, nor evidence that native MATLAB FE assembly has been restored.

## Revised ownership and execution contract

The user's requested direction supersedes CPU-MEX prohibition for this
migration. Keep CPU Python/pybind11 in `radia`; let the MATLAB/Simulink
owner also expose explicit CPU native operations where supported. Both MATLAB CPU and GPU builds use
`mexcuda`; do not introduce a separate `mex` CPU build. CPU execution must
not require GPU initialization. MATLAB owns host scheduling (no NGSolve
TaskManager in these modules). Python builds, ABI and execution follow
NGSolve; caller-owned TaskManager remains the CPU parallel boundary. GPU setup
and execution remain native MATLAB/mexcuda, without mandatory Python GPU
adapters. Reuse NGSolve FE and shared Radia mathematics; do not recreate FE
semantics in MATLAB or call Python implicitly on each Simulink step.

Python GPU execution is explicitly authorized. Share CUDA domain kernels
through thin Python native and MEX adapters; MATLAB Engine is not a required
Python GPU user route. Small independent CPU validation remains mandatory.

Public operations select `Backend="cpu"` or `Backend="gpu"` explicitly where
both are implemented. Unsupported operations fail with a capability error.
No silent fallback is introduced. Native host FE assembly and a full CPU
numerical solver are separate capabilities. A no-GPU CPU route must be
validated on a machine/session without device initialization; merely building
host code with mexcuda does not establish that capability.

## Migration sequence

1. Inventory required MATLAB operations and old Mesh/space/form contracts.
   Record CPU host, CPU solve and GPU support separately. Do not restore a
   historical binary or redirect a missing command to different mathematics.
2. Implement the smallest needed shared-C++/NGSolve CPU native boundary,
   checked handles, errors and lifecycle. Validate assembly from the original
   mesh independently against Python, not by replaying Python-assembled
   matrices. Preserve `sparsecholesky` and original-equation residual gates.
3. Change root policies, distribution metadata, build paths, CI impact scopes,
   release tooling and manifests together in a clean integration checkout.
   Preserve concurrent policy edits and all existing solver research.
4. Publish `radia-simulink` as the sole canonical distribution after exact-wheel
   acceptance. Rename source directory, Python imports, CLI, metadata, build,
   CI, manifests and documentation together. Stop new `radia-simulink` releases.
   The user explicitly rejects backward compatibility: no compatibility wheel,
   aliases, forwarding imports, old-name CLI or silent migration/fallback.
5. Validate a clean installation of the new name, actual import/MATLAB paths
   and hashes, CPU-only operation, Python GPU and MATLAB RTX execution,
   Simulink save/reopen/repeat/error release. Historical immutable releases
   and numerical evidence remain historical records, not supported aliases.

The public PyPI JSON endpoint returned HTTP 404 for radia-simulink on
2026-10-19; that does not reserve the name or establish publication rights.
No new-name wheel publication, installed-wheel acceptance, new CPU solver or
direct Python GPU numerical implementation is claimed. Tracked SLX models
must be edited/checked/saved/reopened through official MATLAB tooling where
renamed entry points affect them; binary model content has not been patched.
