# Four-property producer and Simulink connection gap

Status: implementation plan, not numerical or model acceptance. CPU reference result `four_property_cpu_3c75fc807.json` tests the Python thermal law and energy arithmetic only. Tested implementation: `3c75fc8078c20596ee872f46c4316319b4130e33`; evidence-only revision: `f764add37f0d14fc1eeb165fc0d8e246f3295f0a`.

## Observed runtime boundary

- `makeIHNativeConfig.m` and `validateIHNativeConfig.m` explicitly reject non-linear `bh_mode`; the former still builds the legacy preassembled scalar thermal route.
- `ihEddySFunction.m` owns an `ih.eddy` object and emits heat from externally supplied accepted temperature. `ihThermalSFunction.m` owns a separate `ih.thermal` object and advances temperature in Update. Both call `radia_mex`.
- The new candidate instead lives in `native/radia_cuda_ih.cu`, with a coupled context exposing `create`, `update(current, ambient)`, `outputs`, and `electromagnetic_state`. Replacing only the thermal command would either duplicate an EM update or leave the Eddy output on the old material law.
- `src/runtime/ih/radia_ih_runtime.cpp` is a different runtime implementation. Neither CPU reference results nor compiling the new standalone CUDA entry point proves that the legacy S-Functions use it.

## Required next implementation

1. Keep FE basis, gradients, quadrature and element numbering owned by NGSolve. Extend the physical producer with explicit volume/boundary evaluation maps and weights, mass density, Kelvin cp/kappa/emissivity tables, per-element electrical resistivity/BH data and an independently assembled port-voltage projection. Reject missing fields and unsupported size/geometry; never infer maps from lumped weights or heat being tested.
2. Agree the backend/handle boundary with the shared-runtime owner before edits. A single accepted state must drive both block outputs and exactly one EM-to-thermal update per sample. A shared handle/registry change belongs to the common owner. The target is one latest coupled backend with a documented two-block ownership and schedule contract; changing only a command name is insufficient. Do not retain the old runtime as a selectable alternative or fallback.
3. Extend the IH configuration producer and validator together. Replace the legacy linear-preview contract with the required four-property contract in the same migration; do not add a mode switch that preserves the old fixed-property execution path. Validate Kelvin units, dimensions, FE constant preservation, BH tensor ordering and independent power observation before native allocation.
4. Connect the readable Eddy and Thermal wrappers using the agreed lifecycle: create once, output accepted state, advance once, atomically commit, reset together and destroy once. Preserve previous accepted state on material-domain, nonlinear, residual or energy failures.
5. Validate producer roundtrip and failure paths, then native numerical parity, lifecycle/reset/repeated run, timing and independent FEM. Perform formal SLX editing/check/save/close/reopen through supported MATLAB tooling after resource authorization.

## Latest-only migration and removal gate

User direction: keep no legacy execution route. Migrate every IH producer, validator, readable S-Function, monitor and generated-model caller to the same coupled-state contract. A single owner performs create/reset/destroy and one scheduled update; consumers read the same accepted state. Reject unsupported old configurations rather than converting them to a numerically different route. Remove obsolete IH runtime implementations, scalar-only dispatch branches, smoke configurations and tests only after identifying all repository consumers and replacing each required failure signal in the current route. Shared dispatch/registry/build code remains under its assigned owner; preserve other packages and active external consumers until the coordinated switch. Do not delete another owner's WIP or a running installation. No release may claim this migration complete while a production model can still select the old path.

## Acceptance not yet established

NGSolve-generated physical producer maps, native BH/rho temperature updates, both S-Functions executing the same coupled state, GPU numerical agreement, independent FEM, and formal SLX/package acceptance remain NOT VERIFIED. The existing scalar preview and the independent Python CPU reference cannot satisfy these gates.
