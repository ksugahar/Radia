---
name: pybind-mex-bridge
description: Expose an existing Radia C++/pybind11 capability through a standalone MATLAB MEX ABI and MATLAB wrapper, with numerical-parity, handle-lifecycle, error, and performance validation. Use when adding or reviewing MATLAB parity for a pybind11-backed Radia or NGSolve feature; not for MATLAB-only algorithms or Simulink model layout.
---

# pybind to MEX bridge

Start from the shared C++ numerical kernel, not from the Python wrapper. Keep
pybind11 and MEX as thin adapters over one implementation. If the pybind11
surface contains Python object identity or callbacks, choose a stable numeric,
sparse, file, struct, or checked-handle boundary instead of emulating Python
objects in MATLAB.

## Boundary design

- Prefer numeric arrays, MATLAB sparse matrices, versioned structs, `.vol`
  paths, and independently callable commands.
- Use an opaque `uint64` handle when copying or reconstructing the native object
  would be material. Validate registry membership and object type on every
  command; `create` owns one `mexLock`, and `destroy` releases it exactly once.
- Never expose a raw pointer. MATLAB handle classes own tokens and clear them in
  `delete`; Simulink Level-2 S-Functions release them in `Terminate`.
- Preserve NGSolve ownership of spaces, element transforms, GridFunctions,
  forms, matrices, vectors, and TaskManager work. Do not reconstruct FE basis
  semantics in MATLAB.
- Keep the standalone MEX command useful without Simulink. A readable Level-2
  MATLAB S-Function may compose it later; do not bury a reusable kernel inside
  a MEX S-Function.

## Parity workflow

1. Inventory the public pybind11 operation and identify the shared C++ call it
   already uses. Record any deliberately retained Python-only boundary.
2. Add the MEX command, command inventory entry, strict arity/type/shape checks,
   MATLAB wrapper, and Python-to-MATLAB parity-manifest entry together.
3. Build with `Build.ps1 -MatlabMexOnly`. Test the staged MEX before replacing a
   loaded production binary; never work around a MATLAB file lock by weakening
   the ABI check.
4. Generate expected numerical values by the public Python/pybind11 path on the
   same input. Compare MATLAB MEX results with explicit relative and absolute
   tolerances; do not use handwritten MATLAB goldens for shared behavior.
5. Test invalid arity, scalar type, dimensions, stale/wrong handles, double
   destroy, repeated create/use/destroy cycles, `clear mex`, and final live
   handle count. A failed command must not leak a lock or mutate accepted state.
6. Measure cold load separately from warmed execution and include transfer
   costs. Do not assume MEX or pybind11 is faster without end-to-end evidence.

Small deterministic regressions belong in `tests/`. Multi-mesh, scale,
manufactured-solution, and published numerical evidence belongs in
`validation_test/` with a versioned JSON result. Long MATLAB runs execute on
HIBINO or mdx; routine CI checks the artifact contract and focused fast lane.
