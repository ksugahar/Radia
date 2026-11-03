## HDiv-MMM test tiers

The required mesh matrix is **tet only**, **hex only**, and **tet-hex mixed through
pyramids**. See the canonical [HDiv-MMM CPU/GPU test policy](../../tests/README.md#hdiv-mmm-mesh-families-and-cpugpu-test-tiers).
Small numerical and lifecycle regressions belong in the lightweight test tier;
large physical, nonlinear, mixed-mesh and GPU performance studies belong here
with JSON evidence. Existing tetrahedral results do not certify the other two
families. Pyramid-specific upstream defects do not suspend the other GPU work.

## Current device route

MATLAB HDiv solves take a checked `radia.hdiv-device-setup.v1`; every Gram, charge,
evaluation and mass operator is built on the device:

```matlab
setup = radia.gpu.prepareHDivDeviceSetup(checkedMeshPath, 2, LabelContract=contractPath, ...
    MeshCheckReceipt=receiptPath, MeshCheckReceiptSha256=receiptSha256, ...
    MeshCheckProvenance=provenance);
solution = radia.gpu.solveHDivElements(setup, [1 .2 .1], 2.5, Tolerance=1e-11, ...
    MaxIterations=1000, Preconditioner='mass-cholesky');
```

`validate_curved_hdiv.py` compares that whole solve with independent CPU pybind11;
`PYBIND_MEX_BITWISE.md` describes the same-kernel frontend parity.

Scripts of the retired host-normalized MATLAB route, such as
`validate_scaled_hdiv.py`, `validate_hdiv_artifact.py`, `measure_hdiv_cpu_gpu.py`
and the `check*`/`measure*` nonlinear scripts, were removed with it. Their recorded
JSON stays under `results/` and the sources stay in history. Those results do not
certify the device route.

## Retired host-upload route

The gateway `radia_cuda_hmatrix`, the HDiv commands of `radia_cuda_projection`,
`radia.gpu.HMatrixOperator`, `solveHDivRT0`, `solveCurvedHDiv` and their scripts
took caller-built matrices or host-assembled mass and charge maps, which the
matrix residency policy forbids. They are removed. Their recorded JSON stays
under `results/` as provenance for that route only and certifies nothing about
the device route. The dated narrative of those campaigns is in git history at
`93efc9a60`. [RETIRED_ROUTE_INVENTORY.md](RETIRED_ROUTE_INVENTORY.md) lists the
removed entries, their consumers, their coverage and what replaces it or what was lost.

## Device-route entries

| Entry | Scope | Output |
|---|---|---|
| `validate_curved_hdiv.py` | Whole linear solve through MATLAB against the CPU route on a checked mesh: orders 0, 1, 2; affine, curved or HEX | JSON with the matched-profile label |
| `validate_device_nonlinear.py` | Nonlinear solve and material energy, gradient and tangent on pybind11 against the CPU route, with finite differences | JSON |
| `device_scale.py` | Compressed device H-matrix setup and linear solve, with `--nonlinear` also the BH solve, at one or more meshes; no threshold | JSON |
| `check_pybind_mex_bitwise.py`, `hdiv_device_parity.py`, `produceHDivDeviceMexRecord.m` | Same-kernel frontend parity, see [PYBIND_MEX_BITWISE.md](PYBIND_MEX_BITWISE.md) | JSON |
| `check_analytic_gram_bitwise.py` | Analytic charge Gram action and inverse-BH samples, pybind11 against MEX, bitwise on one kernel object | JSON |
| `checkHDivDeviceMemory.m`, `makeGatewayMemoryWorkload.m`, `sampleGatewayMemory.m` | 100-run memory slope under `gateway_memory_protocol.json` (`hdiv-device` workload) | JSON |

`tests/test_hdiv_gpu_analytic_gram.py` and `tests/test_hdiv_gpu_device_guards.py` are the
surrogate tests of the analytic Gram and of the material and Newton guards on the device.

The `hdiv-device` memory workload does not meet the plateau limit of the protocol. Over 100
runs of the gateway, private bytes rise 4.26 MB, at 38.9 kB per run, against limits of 16 MiB
and 1 MiB per run. The last 25 runs rise 2.49 MB at 49.1 kB per run, against limits of 1 MiB
and 64 KiB per run. Device memory returns to its baseline and no handle stays live. A control
of the same sampling loop without the gateway rises 5.72 MB over its last 25 runs, also above
the limit, so this measurement does not separate the gateway from the sampling loop. The cause
is not identified and the limits are unchanged.
