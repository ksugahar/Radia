# Existing RadFld MEX migration audit

The retired gateway at `b59efa879^:src/matlab/radia_mex.cpp` already implements:

- `ngsolve.radia_field.create`: shared RadiaFieldCoefficient constructor.
- `ngsolve.radia_field.prepare_cache`, `set_memoize`, `clear_cache`,
  `cache_stats`: the same shared CPU coefficient cache implementation.
- `ngsolve.grid_function.interpolate`: NGSolve SetValues with VOL and LocalHeap.
- Source object creation and field command dispatch in the same native module.

Current MATLAB RadiaField/GridFunction wrappers preserve these command names.
The old gateway source was removed in b59efa879; current CMake explicitly
rejects RADIA_BUILD_MATLAB_MEX. Current CUDA-only builders do not implement this
whole command surface. Missing radia_mex on the current MATLAB path was not
proof of missing historical functionality. Existing old binaries found in
other job directories are evidence only, not acceptable migration artifacts.

## Reuse and migration

`packages/radia-gpu/native/radia_projection_host.cuh` preserves the old
interpolation call as a mexcuda host helper with native wall-time profiling.
It deliberately removes the old explicit RegionTaskManager startup; NGSolve
still owns FE mapping and SetValues. The helper is not yet compiled or wired
into a checked native command registry and is not a validated gateway.

Reuse the shared RadiaFieldCoefficient and existing MATLAB commands instead
of inventing a new FE projection. Preserve source/CF/GridFunction checked-handle
ownership in one module; handles cannot be passed between unrelated MEX
registries. The existing numeric CUDA field evaluator must feed values through
a native coefficient compatible with SetValues, with explicit source revision
and cache invalidation. GPU execution is not established by compiling the old
host-only source with mexcuda. Record host setup and device work separately.

The old gateway starts TaskManager in many unrelated commands. Do not simply
copy it wholesale and call the migration accepted. Port and verify those host
regions, supported APIs, dependencies, numerical contracts and lifecycle in the
radia-gpu owner. Python DLL linkage is allowed by the user; a runtime Python
adapter is not. SparseCholesky and original FE semantics remain required.

Next build a narrow owned mexcuda gateway containing source construction,
RadiaField CF and GridFunction ownership plus existing projection/cache
commands. Verify direct and cached projection with the same physical inputs,
then integrate CUDA batching. Keep physical 1e-4 and independent CPU parity
pending until measured. A100 is not operated.
