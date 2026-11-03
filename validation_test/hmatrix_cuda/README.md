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

`validate_curved_hdiv.py` compares that whole solve with CPU pybind11, on curved geometry with the device mass rule ([record](CURVED_GATE_RECORD.md));
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
| `validate_curved_hdiv.py` | Whole linear solve through MATLAB on a checked mesh against the CPU route, on curved geometry against the CPU Gram with the device mass rule: orders 0, 1, 2; affine, curved or HEX | JSON with the matched-profile label |
| `validate_device_nonlinear.py` | Nonlinear solve and material energy, gradient and tangent on pybind11 against the CPU route, with finite differences | JSON |
| `device_scale.py` | Compressed device H-matrix setup and linear solve, with `--nonlinear` also the BH solve, at one or more meshes; no threshold | JSON |
| `check_pybind_mex_bitwise.py`, `hdiv_device_parity.py`, `produceHDivDeviceMexRecord.m` | Same-kernel frontend parity, see [PYBIND_MEX_BITWISE.md](PYBIND_MEX_BITWISE.md) | JSON |
| `check_analytic_gram_bitwise.py` | Analytic charge Gram action and inverse-BH samples, pybind11 against MEX, bitwise on one kernel object | JSON |
| `checkHDivDeviceMemory.m`, `makeGatewayMemoryWorkload.m`, `sampleGatewayMemory.m` | 100-run memory slope under `gateway_memory_protocol.json` (`hdiv-device` workload) | JSON |

`tests/test_hdiv_gpu_analytic_gram.py` and `tests/test_hdiv_gpu_device_guards.py` are the
surrogate tests of the analytic Gram and of the material and Newton guards on the device.

The `hdiv-device` memory workload meets the limits of the protocol on both device hosts, in two
fresh MATLAB runs each. Over 100 runs of the gateway, the limits are 16 MiB and 1 MiB per run; over
the last 25 runs they are 1 MiB and 64 KiB per run.

| Device | 100 runs | Per run | Last 25 runs | Per run |
|---|---|---|---|---|
| RTX 5000 | 1.20, 1.26 MB | 6.1, 5.0 kB | 0.00, 0.14 MB | -0.3, 7.9 kB |
| A100 | 2.58, 2.35 MB | 32.0, 32.0 kB | 0.00, 0.00 MB | 0.0, 0.0 kB |

Device memory returns to its baseline and no handle stays live on both. Passing the limits is not
proof of zero leaks.

The workload failed the last-25-run limits before `readHDivDeviceSetup` stopped checking an unchanged
payload again (RTX 5000 2.70 MB, A100 2.06 MB). On the A100, runs of one layer at a time showed late
one-off steps when the checks ran on every solve and none when the same solves skipped them, so a
changed payload is now the only one checked again. What inside MATLAB makes those steps is not
identified. Reading the process counter also adds 64 KiB per read for about the first 66 reads, which
ends before the last 25 runs. On the RTX 5000, private memory still alternates by about 0.52 MB per
run, cause not identified. An idle control on that shared desktop host moved by up to 5.77 MB in one
step, so a later result there can be disturbed by the host.
