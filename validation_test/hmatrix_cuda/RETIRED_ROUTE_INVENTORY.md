# Retired HDiv host-upload route: inventory

The matrix residency policy forbids moving a matrix between CPU and GPU, so a
route where the caller or host assembles the mass, charge map or H-matrix leaves
and uploads them is not a GPU implementation. This file records what left the
repository since `93efc9a60`, who used it, what it covered, what stands in its
place and what was lost.

Removal steps:

| Commit | What it removed |
|---|---|
| `0bab6d074` | Host-supplied P1 diagonal route, its plan and criterion, `validate_python_hdiv_resident.py`, `validate_hdiv_nonlinear_manufactured.py` |
| `8fc81f915` | `packages/radia-simulink/tests/test_hdiv_host_publisher.py` |
| `d38eca2d9` | `tests/matlab/test_hdiv_shared_newton_api.m`, `test_prepared_hdiv_host_setup.m`, `test_prepared_hdiv_shared_workflow.m` |
| `e16e13f26` | 24 scripts bound to the mesh-path signatures, see the correction below |
| `492612d92` | The gateway, its MATLAB entries and native code, the remaining scripts |
| The commit adding this table | `preflight_hdiv_scale.py`, `curved_projection_validation.contract.json`, `hmatrix_leaf_coverage.hpp`, `test_hmatrix_leaf_coverage.cpp`, the Python fault-injection hooks and four retired-route plans |

## Production entries removed

| Entry | Consumers at removal | Replacement |
|---|---|---|
| `radia.gpu.HMatrixOperator` | Shipped: a wheel wrapper in `hdiv_payload.WRAPPERS`, so a released public class. Also `verify_wheel.py`, `test_hdiv_payload.py`, both parity manifests and the scripts listed below. No block or MCP tool | HDiv: `HDivDeviceOperator`, `HDivDeviceMass`. Generic caller-supplied leaves: none, see below |
| `radia.gpu.affineHDivQuadrature`, `curvedHDivQuadrature` | Shipped wheel wrappers. Configuration of the removed gateway, `test_affine_hdiv_quadrature.m`, `makeGatewayMemoryWorkload.m` | Rules staged by `radia_simulink.hdiv` from the CPU route |
| `radia.gpu.solveHDivRT0` | Excluded from the wheel. `compareHDivAssembly.m`, `validate_native_rt0.py`, `measure_hdiv_assembly.py`, `parity-manifest.json` | `prepareHDivDeviceSetup(mesh, 0, ...)` + `solveHDivElements` |
| `radia.gpu.solveCurvedHDiv` | Excluded from the wheel | `prepareHDivDeviceSetup(mesh, order, CurveOrder=2, ...)` + `solveHDivElements` |
| `radia.gpu.evaluateHDivBHMaterial`, `evaluateHDivBHChunks` | Removed scripts only | `HDivBHMaterial` over the device setup. Chunked evaluation has no replacement: the device route does not chunk |
| `radia_cuda_hmatrix` gateway, `tools/build_hmatrix_cuda.py`, `Build.ps1 -HMatrixCudaOnly`, its native `hmatrix_*` and `hdiv_rt0_host`, `hdiv_element_setup` headers | The entries above, the wheel payload, `makeGatewayMemoryWorkload.m` | `radia_hdiv_device` |
| HDiv commands of `radia_cuda_projection` (`assemble_hdiv_*`, `inspect_hdiv_elements`, `hdiv_material_*`) | The entries above | `python -m radia_simulink.hdiv`, which stages NGSolve tables only |

### Effect on MATLAB users

- **Physics.** RT0, BDM1 and BDM2 on affine and P2 tetrahedra and hexahedra, and
  the BH law at orders 1 and 2, all have a device route.
- **MATLAB without Python.** The removed solves ran NGSolve in the MEX. Their
  replacement has one Python step at the initialization boundary, which may run
  beforehand: `python -m radia_simulink.hdiv request.json setup.h5`. That Python
  needs `radia[gpu]`, `radia-simulink` and `h5py`. It loads the pip CUDA runtime
  because staging imports `_hdiv_cuda`; running it without a GPU has not been
  tested. MATLAB then needs only `radia.gpu.loadHDivDeviceSetup('setup.h5')`.
  Staging builds no matrix.
- **Checked input.** The old entries took a bare mesh path. The route now requires
  a versioned label contract, a `check-vol` receipt and its provenance.
- **Caller-supplied leaves, mass or charge map.** No replacement. This is the
  input the residency policy forbids. `HMatrixOperator` was a released public
  class, so this removal is a user decision and is flagged as such.

## What stays, and is not part of the HDiv route

`radia_cuda_projection` keeps `create`, `targets`, `project`, `destroy`,
`live_count` and its lifecycle diagnostics. It maps a batched field onto NGSolve
targets and owns no matrix. Consumers: `validation_test/radfld_cuda/`,
`checkProjection*.m`, `makeGatewayMemoryWorkload.m`. It is still built in CI and
is no longer staged in the wheel: no shipped wrapper calls it now. Its last
shipped caller was the mesh-path route of `solveHDivElements`.

Pyramid and mixed meshes: the CPU charge Gram refuses them
(`radia.vim` "Radia pyramid and MIXED-element charge-Gram are not implemented"),
and the device route refuses them through `_mesh_family`. No solve consumer
existed. `tests/test_hdiv_pyramid_material.py` covers the NGSolve material
evaluation only and is unchanged; the status is unchanged.

## Script consumers of `HMatrixOperator` at removal

| Consumer | Kind | Disposition |
|---|---|---|
| `tests/matlab/test_hmatrix_normalized_leaves.m` | Test of the generic leaf boundary | Removed with its subject; manifest entry `normalized-leaf-resident-boundary` removed |
| `checkResidentHMatrix.m`, `profileResidentHMatrix.m`, `validate_resident.py`, `measure_resident.py`, `validate_same_artifact.py` | Generic resident apply | Removed with their subject |
| `checkChebyshevLinear.m`, `checkPolynomialLinearGraph.m` | HDiv solve options | `test_hdiv_gpu_linear_material.py` runs every preconditioner on the pybind11 route |
| `checkHDivSparseReuse.m` | Operator reuse across RHS scalings | `test_device_solve_reuses_its_operator_across_rhs_scalings`. The invalid-mass rejection has no subject: no caller mass exists |
| `checkNormalizedCoverageFault.m` | Partition and allocation fault of the gateway | No counterpart. The partition check it exercised had no other caller and was removed |

## Validation scripts removed, and their coverage

No removed script was referenced by a record in `validation_records/`; no
adopted case was lost. Recorded JSON stays under `results/` with its original
provenance and certifies the removed route only. `gates/p1-*.json` was deleted
and remains in git history.

"Pybind11" below means the entry runs `_hdiv_cuda` through Python, not the MEX.
The two frontends share one kernel object; their bitwise equality is recorded for
`device-flat-tet` only (`PYBIND_MEX_BITWISE.md`).

| Coverage | Removed scripts | Standing on the device route |
|---|---|---|
| Whole linear solve vs CPU, through MATLAB | `checkHDivHacapk*.m`, `checkHDivCachedEntries.m`, `checkHDivIndexedEntries.m`, `checkHDivFusedMoments.m`, `checkHDivDeviceDense.m`, `validate_native_rt0.py`, `compareHDivAssembly.m`, `measure_hdiv_assembly.py`, `validate_hdiv_artifact.py` | `validate_curved_hdiv.py`: orders 0, 1, 2 with `--geometry affine`, `curved` or `hex` |
| Same, on pybind11 | | Lanes `test_hdiv_gpu_linear_material.py`, `test_hdiv_gpu_curved_device.py`, `test_hdiv_gpu_hex_device.py` |
| Host staging of sparse mass and charge map | `checkMaterialCSC.m`, `checkMaterialContext.m`, `checkSparseStaging.m`, `validate_element_setup.py`, `checkHDivMassPreconditioner.m`, `checkHDivSparseMass.m` | No host copy exists. Device mass and map are compared directly in the curved and HEX lanes, and through the solve in the affine lane |
| Gram, B and mass actions vs CPU | `checkHDivCPUActions.m`, `validate_python_hdiv_resident.py` | Same lanes, on pybind11 |
| Chebyshev and polynomial settings | `checkChebyshev.m`, `checkChebyshevSettings.m`, `screenChebyshevSettings.m`, `checkPolynomialGraph.m` | Fixed settings in the linear lane. The settings screen was a one-off and has no replacement |
| Nonlinear whole solve vs CPU coefficients | `checkNonlinearSolve.m`, `checkPreparedNonlinear.m`, `prepare_nonlinear_solve.py`, `checkFinalResidual.m`, `checkFrameAdmission.m` | `validate_device_nonlinear.py` against `radia.vim.Solve` on pybind11. It shares the near/far profile, so it is implementation parity, not independent accuracy. Ideal sphere physics: `validate_nonlinear_sphere.py`. Frame coefficient telemetry has no replacement |
| Known-root nonlinear solve | `validate_hdiv_nonlinear_manufactured.py` | **Lost.** Its analytic cubic material and manufactured root were independent of the Newton core. `radia.vim.Solve` is the nearest independent check |
| BH energy, gradient, tangent vs CPU and finite differences | `checkNonlinearMaterial.m`, `checkNonlinearChunks.m`, `prepare_nonlinear_material.py` | `validate_device_nonlinear.py` on pybind11. Through the MEX only `device-flat-tet` bitwise parity. Device chunking does not exist, so per-chunk accumulation has no subject |
| Near-root behaviour | `diagnoseNearRootNonlinear.m` | Recorded as a diagnostic in `validate_device_nonlinear.py` |
| MATLAB shared workflow | `test_hdiv_shared_newton_api.m`, `test_prepared_hdiv_host_setup.m`, `test_prepared_hdiv_shared_workflow.m` | `tests/matlab/test_prepared_hdiv_device_setup.m`. **No workflow invokes it yet**, it was run by hand on the candidate |
| 100-run process and device memory slope | `checkHostSetupMemory.m` | `checkHDivDeviceMemory.m` (`makeGatewayMemoryWorkload.m` case `hdiv-device`) under `gateway_memory_protocol.json` |
| Linear and nonlinear device scale, GPU only | `measurePreparedNonlinear.m`, `measureStrongNonlinear.m`, `validate_scaled_hdiv.py`, `prepare_hdiv_scale.py`, `preflight_hdiv_scale.py` | `device_scale.py`, with `--nonlinear`. Records only; no speed or scale threshold |
| CPU/GPU and ABBA timing | `measure_hdiv_cpu_gpu.py`, `measureChebyshev.m`, `measureChebyshevSettings.m`, `measureFinalResidual.m`, `measurePolynomialGraph.m`, `measureMaterialCSC.m` | None. One-off comparisons are not repeated |
| Native allocation and MEX fault injection of the gateway | `checkAllocationFailure.m`, `checkCurvedAllocationFailure.m`, `checkMxCommandFailure.m`, `checkMxOutputFailure.m`, `checkResourceFailure.m`, `gates/p1-*.json`, `check_c_memory_transaction.c`, `check_host_fault_allocator.cpp`, `check_entry_batch.*`, `check_full_host.*`, `check_leaf_bounds.*`, `check_serial_host.c`, `check_aca_cpu.c`, `check_aca_regression.c` | Subject removed. No device-route counterpart was written: a 100-run memory slope covers retention, not allocation failure paths |
| MATLAB output fixture | `prepareCurvedMxOutputFixture.m` | Subject removed |
| Scale estimate formulas | `tests/test_hdiv_scale_estimates.py` | Subject removed |

The retired-route plans `large-memory-plan.md`, `nonlinear-chunking-plan.md`,
`nonlinear-bh-cuda-plan.md` and `research-claims-20261018.md` described the host
route as current and were removed; they remain in git history.

### Correction to the record of `e16e13f26`

The commit deleted 24 files: 21 MATLAB scripts and 3 Python runners. Its message
counts them differently but consistently (20 research scripts, 3 runners and one
helper).

Its reason, that the scripts "could no longer run" because the mesh-path and
positional-option signatures of `solveHDivElements` and `PreparedHDivNonlinear`
were removed, is false for three of them:

- `checkNonlinearMaterial.m` and `checkNonlinearChunks.m` called
  `radia_cuda_projection` and `evaluateHDivBH*`, which still existed at that commit.
- `checkHDivCPUActions.m` called `radia_cuda_hmatrix` and `radia_cuda_projection`,
  which still existed at that commit.

Their calls still resolved; this was read from the source, not executed. Their
coverage was not obsolete, it was unique: the first two have no replacement until
`validate_device_nonlinear.py` in this change. The history is kept as it is; this
paragraph corrects the reading of it.
