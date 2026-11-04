## HDiv-MMM test tiers

The required mesh matrix is **tet only**, **hex only**, and **tet-hex mixed through
pyramids**. See the canonical [HDiv-MMM CPU/GPU test policy](../../tests/README.md#hdiv-mmm-mesh-families-and-cpugpu-test-tiers).
Small numerical and lifecycle regressions belong in the lightweight test tier;
large physical, nonlinear, mixed-mesh and GPU performance studies belong here
with JSON evidence. Existing tetrahedral results do not certify the other two
families. The pyramid map of Netgen 6.2.2608 scales every pyramid `z` by `1 - 1e-10`, so
pyramid results on the stock wheel differ from an exact map at that level. The installed-wheel
GPU lane runs on that source plus `tools/patches/netgen-pyramid-apex.patch`, built and placed in
its run-local environment by `tools/netgen_overlay.py`.

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
| `validate_mixed_sphere.py` | Ideal sphere on connected tet, pyramid and hex meshes: mesh and order convergence of the device route against the analytic magnetization, with the CPU route on the coarse levels; the limit is proposed | JSON |
| `validate_hdiv_polynomial_charges.py`, `polynomial_box_charges.py` | Charge Gram of polynomial magnetization on a cube against its closed box integrals, independent of the mesh: typed and flat routes, operator, solution and mass differences by profile; no threshold | JSON |
| `check_pybind_mex_bitwise.py`, `hdiv_device_parity.py`, `produceHDivDeviceMexRecord.m` | Same-kernel frontend parity, see [PYBIND_MEX_BITWISE.md](PYBIND_MEX_BITWISE.md) | JSON |
| `check_analytic_gram_bitwise.py` | Analytic charge Gram action and inverse-BH samples, pybind11 against MEX, bitwise on one kernel object | JSON |
| `run_device_boundary.py`, `runDeviceBoundary.m` | That frontend parity, the setup read controls, the MATLAB device workflow and the fault sweep on one exact build; runs in the GPU lane of `radia-simulink.yml` | JSON |
| `device_fault/` | Failure paths of the device commands, see below | JSON |
| `checkHDivDeviceMemory.m`, `makeGatewayMemoryWorkload.m`, `sampleGatewayMemory.m` | 100-run memory slope under `gateway_memory_protocol.json` (`hdiv-device` workload) | JSON |

`tests/test_hdiv_gpu_analytic_gram.py` and `tests/test_hdiv_gpu_device_guards.py` are the
surrogate tests of the analytic Gram and of the material and Newton guards on the device.

The `hdiv-device` memory workload meets the limits of the protocol on both device hosts, in two
fresh MATLAB runs each. Over 100 runs of the gateway, the limits are 16 MiB and 1 MiB per run; over
the last 25 runs they are 1 MiB and 64 KiB per run.

| Device | 100 runs | Per run | Last 25 runs | Per run |
|---|---|---|---|---|
| RTX 5000 | 1.09, 6.93 MB | 4.7, 55.3 kB | 0.00, 0.00 MB | 0.0, -0.2 kB |
| A100 | 2.39, 2.43 MB | 32.0, 32.0 kB | 0.00, 0.00 MB | 0.0, 0.0 kB |

Device memory returns to its baseline and no handle stays live on both. Passing the limits is not
proof of zero leaks.

The workload failed the last-25-run limits before `readHDivDeviceSetup` stopped checking an unchanged
payload again (RTX 5000 2.70 MB, A100 2.06 MB). On the A100, runs of one layer at a time showed late
one-off steps when the checks ran on every solve and none when the same solves skipped them, so a
changed payload is now the only one checked again. What inside MATLAB makes those steps is not
identified. The comparison uses builtins only: an earlier version written as a recursive MATLAB
function failed one of two RTX 5000 runs (16.78 MB, last 25 runs 5.70 MB) and added a one-off step of
6 to 9 MB on both hosts.

Reading the process counter also adds 64 KiB per read for about the first 66 reads, which ends before
the last 25 runs. On the RTX 5000, private memory can alternate by up to 1.57 MB per run, cause not
identified; it returns to the same value every other run, so the last-25-run difference does not see
it. An idle control on that shared desktop host moved by up to 5.77 MB in one step, so a later result
there can be disturbed by the host.

## Connected tet, pyramid and hex meshes

- Charge Gram: every pair takes the typed route on the trilinear image of the unit cube, on CPU and device, from one
  source (`src/core/rad_hdiv_field_*.inc`, `src/core/gpu/hdiv/hdiv_typed_host_math.inc`). Its 1D rules come from
  `_mixed_gram_rules` in `src/radia/vim/_vim.py` as kernel arguments; the near band and the near and pyramid counts are
  proposed values, not an adopted criterion.
- Rational pyramid charge: on an affine pyramid each level `z` of `xy/(1-z)` is a parallelogram carrying a polynomial,
  whose potential is exact, and the rule in `z` is split at the level of the point. The radial rule it replaced did not
  converge at order 2 on the small meshes of the record. A non-parallelogram pyramid keeps a radial rule on all charges.
- Measured limits of the counts, on small meshes against the finest profile of the same route, so they say how far the
  rule converged and nothing about the formulation. Order 2 on affine pyramids converges geometrically when the pyramid
  and near counts rise together; at the proposed (outer 6, pyramid 14, near 10) the field difference is about `5e-5`.
  Order 2 on non-parallelogram pyramids is slower and not monotone: the field difference stays near `3e-5` to `4e-5` from
  (6, 14, 16) to (6, 20, 16). Order 1 on the same mesh reaches about `2e-6` at (5, 14, 16), and its proposed (5, 6, 10)
  gives about `2e-4`. The flat tet route at its default is about `1e-3` in the field at order 2 on the cube of the record,
  and it needs intorder 60 and far quadrature 14 for about `3e-6`; the typed route is not its default.
- Reference: a polynomial magnetization has no jump across interior faces, so its charge on any mesh of a cube is that
  of the cube, whose box integrals have closed forms after a Duffy map at the corner (Duffy, SIAM J. Numer. Anal. 19, 1982;
  Bailey, Borwein and Crandall, J. Comput. Appl. Math. 206, 2007, for the cube and square values it reproduces). It covers
  polynomial directions only. The rational charge lies outside it. On one pyramid, `rational_pyramid_reference.py` reduces
  its energies to separated one-dimensional weights without the Gram; on a mesh the record only compares the typed route
  with its own finer profiles.
- Mass: the default rule of `u*v*dx` is too short on non-parallelogram pyramids and hexahedra. The device mass, the
  nonlinear material and the mixed CPU route use `IntegrationRuleSpace(order=max(6, 2p+2))`.
- Newton: each linearization builds the frozen tangent `W` on the mass pattern on the device, and `evaluate` keeps the
  matrix-free action; the two share one element table, so their agreement is an implementation check. The defaults are
  `reference-diagonal`, inner tolerance `1e-6`, inner cap 30000. A preconditioner from the tangent was measured on the
  ball at that inner tolerance and is not provided: its diagonal failed at order 1 and took longer at order 2.

## Device failure paths

`device_fault/` builds a test-only MEX, `radia_hdiv_device_fault`, from the unmodified adapter source and the same
kernel object as the shipped MEX. A shim takes the place of `cudaMalloc`, `cudaFree` and `cudaStreamSynchronize` at
link time, and a force-included header wraps the mx allocations of the adapter, so no shipped source carries a hook.
It also arms the host and resource counters the kernel already has.

`checkDeviceFaults.m` fails one site at a time on a low-rank affine input, the hexahedral pair and a connected tet,
pyramid and hex row. Each failure must raise `radia:hdiv:device`, publish no handle and return device and host memory to
their level before the call; the same call must then give the same bytes. A stream wait in a destructor cannot raise,
so it may pass silently if the answer is unchanged; those are counted. The solve, assembly and stats commands share one
kernel across mesh families, so the affine input sweeps them, each preconditioner path included, and the other two sweep
what differs by family: creation, apply, mass and material evaluation.

The first 40 sites of a call are all failed and later ones stepped to about 30 more, so this is not every site, and a
solve repeats its call only at every eighth site, as it is slow. CUDA
library workspaces of cublas and cusolver do not pass through the shim. Curved geometry is not covered: no small
curved fixture exists. The sweep works on small inputs only and says nothing about real out-of-memory behaviour at scale.
