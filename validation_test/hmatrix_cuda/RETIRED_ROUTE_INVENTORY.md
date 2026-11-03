# Retired HDiv host-upload route: inventory

The matrix residency policy forbids moving a matrix between CPU and GPU, so a
route where the caller or host assembles the mass, charge map or H-matrix leaves
and uploads them is not a GPU implementation. This file records what left the
repository, who used it, what it covered and what stands in its place. Two
removal steps are recorded: `e16e13f26` (scripts bound to the removed mesh-path
signatures) and this one (the gateway itself).

## Production entries removed

| Entry | Consumers at removal | Replacement |
|---|---|---|
| `radia.gpu.solveHDivRT0` | None but `setup.py` exclusion and the removed `compareHDivAssembly.m`, `validate_native_rt0.py` | `prepareHDivDeviceSetup(mesh, 0, ...)` + `solveHDivElements` |
| `radia.gpu.solveCurvedHDiv` | None but `setup.py` exclusion | `prepareHDivDeviceSetup(mesh, order, CurveOrder=2, ...)` + `solveHDivElements` |
| `radia.gpu.HMatrixOperator` | Tests and validation only, listed below; no block, wrapper or MCP tool | HDiv: `HDivDeviceOperator`, `HDivDeviceMass`. Generic caller-supplied leaves: none, see below |
| `radia.gpu.evaluateHDivBHMaterial`, `evaluateHDivBHChunks` | Removed scripts only | `HDivBHMaterial` over the device setup |
| `radia.gpu.affineHDivQuadrature`, `curvedHDivQuadrature` | Configuration of the removed gateway and `test_affine_hdiv_quadrature.m` | Rules staged by `radia_simulink.hdiv` from the CPU route |
| `radia_cuda_hmatrix` gateway, `tools/build_hmatrix_cuda.py`, `Build.ps1 -HMatrixCudaOnly`, its native `hmatrix_*` and `hdiv_rt0_host`, `hdiv_element_setup` headers | The entries above and the wheel payload | `radia_hdiv_device` |
| HDiv commands of `radia_cuda_projection` (`assemble_hdiv_*`, `inspect_hdiv_elements`, `hdiv_material_*`) | The entries above | `python -m radia_simulink.hdiv`, which stages NGSolve tables only |

### Effect on MATLAB users

- **MATLAB without Python.** The removed solves ran NGSolve in the MEX. Their
  replacement has one Python step at the initialization boundary, which may run
  beforehand on any machine with `radia[gpu]` and `h5py`:
  `python -m radia_simulink.hdiv request.json setup.h5`. MATLAB then needs only
  `radia.gpu.loadHDivDeviceSetup('setup.h5')`. Staging builds no matrix.
- **Caller-supplied leaves, mass or charge map.** No replacement. No production
  consumer existed, and this is the input the residency policy forbids.

## What stays, and is not part of the HDiv route

`radia_cuda_projection` keeps `create`, `targets`, `project`, `destroy`,
`live_count` and its lifecycle diagnostics. It maps a batched field onto NGSolve
targets and owns no matrix. Consumers: `validation_test/radfld_cuda/`,
`checkProjection*.m`, `makeGatewayMemoryWorkload.m`. It is still built in CI and
is no longer staged in the wheel, because the removed solves were its only
shipped caller.

Pyramid and mixed meshes: the CPU charge Gram refuses them
(`radia.vim` "Radia pyramid and MIXED-element charge-Gram are not implemented"),
and the device route refuses them through `_mesh_family`. No solve consumer
existed. `tests/test_hdiv_pyramid_material.py` covers the NGSolve material
evaluation only and is unchanged; the status is unchanged.

## Consumers of `HMatrixOperator` at removal

| Consumer | Kind | Disposition |
|---|---|---|
| `tests/matlab/test_hmatrix_normalized_leaves.m` | Test of the generic leaf boundary | Removed with its subject; manifest entry `normalized-leaf-resident-boundary` removed |
| `checkResidentHMatrix.m`, `profileResidentHMatrix.m`, `validate_resident.py`, `measure_resident.py`, `validate_same_artifact.py` | Generic resident apply | Removed with their subject |
| `checkChebyshevLinear.m`, `checkHDivSparseReuse.m`, `checkPolynomialLinearGraph.m` | HDiv solve options | Covered on the device by `test_hdiv_gpu_linear_material.py` (every preconditioner) |
| `checkNormalizedCoverageFault.m` | Partition fault | `tests/cpp/test_hmatrix_leaf_coverage.cpp` keeps the shared check |

## Validation scripts removed, and their coverage

No removed script was referenced by a record in `validation_records/`; no
adopted case was lost. Recorded JSON stays under `results/` with its original
provenance. It certifies the removed route only.

| Coverage | Removed scripts | Standing on the device route |
|---|---|---|
| Whole linear solve vs CPU | `checkHDivHacapk*.m`, `checkHDivCachedEntries.m`, `checkHDivIndexedEntries.m`, `checkHDivFusedMoments.m`, `checkHDivDeviceDense.m`, `validate_native_rt0.py`, `compareHDivAssembly.m`, `measure_hdiv_assembly.py` | `validate_curved_hdiv.py`, lanes `test_hdiv_gpu_linear_material.py`, `test_hdiv_gpu_curved_device.py`, `test_hdiv_gpu_hex_device.py` |
| Host staging of sparse mass and charge map | `checkMaterialCSC.m`, `checkMaterialContext.m`, `checkSparseStaging.m`, `validate_element_setup.py`, `checkHDivMassPreconditioner.m` | No host copy exists. Device mass and map are compared in the same lanes |
| Nonlinear whole solve vs independent CPU coefficients | `checkNonlinearSolve.m`, `checkPreparedNonlinear.m`, `prepare_nonlinear_solve.py` | `validate_device_nonlinear.py`; ideal sphere physics: `validate_nonlinear_sphere.py` |
| BH energy, gradient, tangent vs CPU and finite differences | `checkNonlinearMaterial.m`, `checkNonlinearChunks.m`, `prepare_nonlinear_material.py` | `validate_device_nonlinear.py`. Device chunking does not exist, so per-chunk accumulation has no subject |
| Near-root behaviour | `diagnoseNearRootNonlinear.m` | Recorded as a diagnostic in `validate_device_nonlinear.py` |
| 100-run process and device memory slope | `checkHostSetupMemory.m` | `checkHDivDeviceMemory.m` (`makeGatewayMemoryWorkload.m` case `hdiv-device`) under `gateway_memory_protocol.json` |
| Scale and timing without a dense oracle | `measurePreparedNonlinear.m`, `measureStrongNonlinear.m`, `measure_hdiv_cpu_gpu.py`, `validate_scaled_hdiv.py`, `prepare_hdiv_scale.py`, `prepareCurvedMxOutputFixture.m` | `device_scale.py`. Records only; no speed or scale threshold |
| Native allocation and MEX fault injection of the gateway | `checkAllocationFailure.m`, `checkCurvedAllocationFailure.m`, `checkMxCommandFailure.m`, `checkMxOutputFailure.m`, `checkResourceFailure.m`, `gates/p1-*.json`, `check_c_memory_transaction.c`, `check_host_fault_allocator.cpp`, `check_entry_batch.*`, `check_full_host.*`, `check_leaf_bounds.*`, `check_serial_host.c`, `check_aca_cpu.c`, `check_aca_regression.c` | Subject removed. No device-route counterpart was written: a 100-run memory slope covers retention, not allocation failure paths |

### Correction to the record of `e16e13f26`

Its message counts 23 scripts and 20 research scripts. The commit deleted 24
files: 21 MATLAB scripts and 3 Python runners.

Its reason, that the scripts "could no longer run" because the mesh-path and
positional-option signatures of `solveHDivElements` and `PreparedHDivNonlinear`
were removed, is false for three of them:

- `checkNonlinearMaterial.m` and `checkNonlinearChunks.m` called
  `radia_cuda_projection` and `evaluateHDivBH*`, which still existed at that commit.
- `checkHDivCPUActions.m` called `radia_cuda_hmatrix` and `radia_cuda_projection`,
  which still existed at that commit.

They still ran. Their coverage was not obsolete, it was unique: the first two
have no replacement until `validate_device_nonlinear.py` in this change. The
history is kept as it is; this paragraph corrects the reading of it.
