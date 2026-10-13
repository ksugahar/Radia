# GPU IC(0) acceptance

`radia_cuda_ic` owns a checked real FP64 CSR operator, natural-order
fixed-shift IC(0) factors and reusable PCG workspace. CUDA performs numerical
weighted LDLT factorization, then converts it to Cholesky form. cuSPARSE SpSV
applies the forward and transpose triangular solves. Host code performs
symbolic indexing and iteration control, including scalar readbacks. There
is no CPU numerical factorization or CPU solve fallback.

`radia.gpu.ICContext` accepts GPU CSR inputs. `radia.gpu.HCurlICSystem` builds
the checked affine lowest-order HCurl operator directly within MATLAB/mexcuda
and retains its matrix, factor and workspace. Geometry, boundary conditions,
reluctivity and positive gauge mass are fixed; only current scale changes.
CPU/Python is used for independent validation, not as a GPU-user prerequisite.

An optional explicit GPU `free` mask identifies homogeneous constrained
DOFs. Only already-projected zero constrained couplings are accepted. Checked
zero constrained diagonals are extended by identity; constrained RHS entries
must remain zero. This preserves the free-DOF equations, with zero constrained
solution entries. Without a mask, every row must have a positive diagonal.
Missing/duplicate/unsorted CSR entries, asymmetric or nonfinite data,
nonpositive factor pivots, nonpositive PCG curvature, nonconvergence and stale
or foreign handles raise. The diagonal shift is fixed (default 1.05), never
silently increased.

`run.py --binary-dir <candidate> --matlab-dir <payload> --out <new-directory>`
runs the independent NumPy original-A solution, full-pattern shifted solve,
and separately written square-root sparse IC(0) application oracles. It also
checks the MATLAB context, explicit constraints, rejected operations, handle
generations and cleanup. Actual loaded module paths and SHA-256 hashes are
verified in the owned MATLAB Engine. Installed-wheel acceptance invokes the
same checks when the IC payload is present. NumPy/SciPy are test dependencies.

The six-route FEM study is
`validation_test/ams_cuda/measure_scale_installed.py --job <checked-job>
--wheel <exact-wheel> --ic-directory <IC-payload> --out <new-directory>
--cycles 2`. The checked unit cubes contain 10,460 / 32,233 / 106,272 total
DOFs (7,571 / 26,839 / 91,533 free DOFs). CPU IC and AMS run with 1 and 16
NGSolve threads; GPU IC and AMS use explicit separate routes. Fixed-shift
natural-order CPU IC has no ABMC or diagonal scaling. CPU CG uses NGSolve's
preconditioned residual stopping norm; GPU IC uses the unpreconditioned
relative residual. Both independently check the original free-DOF equations,
so small iteration-count differences are expected.

The study records fresh mesh/operator/RHS/factor construction, transfer-inclusive
first solve, and matrix/factor reuse at current scales 2, 0.5 and 1. Every
solve starts from zero and is checked against a separate Python NGSolve
`sparsecholesky` solution and original-A residual. Startup and first initialized
calls are separate; no OS/device cold-start claim is made. Process memory
sampling includes loaded runtimes and validation/serialization; CPU RSS also
includes the retained independent reference operator. Native owned buffers,
queried SpSV workspaces, symbolic transfer counts and device-free endpoints
are reported separately. These do not certify a transient whole-device peak.

100 is a numerical/build acceptance host. Its `--diagnostic` results do not
certify performance and retain explicit memory admission checks. Formal
same-host A100 timing keeps CPU/GPU idle admission thresholds. The A100 host
was powered off on 2026-10-13, so the new IC implementation has no A100
performance or same-binary runtime acceptance yet. The pending record must
remain pending until that host is available and the exact accepted candidate
artifact is loaded and verified. The published 0.1.4 wheel is unchanged.

Candidate numerical acceptance on 100 (2026-10-13): installed CI wheel
`0.1.5rc1`, run `38387288896`, passed all 360 independent FEM solution
checks. Maximum relative solution error was 1.122e-11;
maximum original-equation true residual was 5.621e-12.
The installed-module oracle and 24 rejection gates also passed.
The full phase, transfer, reuse, memory and artifact records are in
[100-ci-candidate-scaling.json](results/100-ci-candidate-scaling.json) and
[100-installed-acceptance.json](results/100-installed-acceptance.json).
All times below are diagnostic 100 values, not formal performance acceptance.

| Total DoFs | Route | Fresh median (s) | Reused median (s) |
| ---: | --- | ---: | ---: |
| 10,460 | cpu-ams-1 | 0.3140 | 0.0621 |
| 10,460 | cpu-ic-1 | 0.3851 | 0.2310 |
| 10,460 | gpu-ams-1 | 0.5074 | 0.0321 |
| 10,460 | cpu-ams-16 | 0.4150 | 0.0353 |
| 10,460 | cpu-ic-16 | 0.4886 | 0.3631 |
| 10,460 | gpu-ic-1 | 0.8602 | 0.3520 |
| 32,233 | cpu-ams-1 | 1.0019 | 0.1871 |
| 32,233 | cpu-ic-1 | 1.5506 | 1.0620 |
| 32,233 | gpu-ams-1 | 1.0225 | 0.0299 |
| 32,233 | cpu-ams-16 | 0.7345 | 0.0565 |
| 32,233 | cpu-ic-16 | 0.8487 | 0.5314 |
| 32,233 | gpu-ic-1 | 1.6980 | 0.7753 |
| 106,272 | cpu-ams-1 | 3.5753 | 0.8050 |
| 106,272 | cpu-ic-1 | 7.0548 | 5.3869 |
| 106,272 | gpu-ams-1 | 2.9077 | 0.0651 |
| 106,272 | cpu-ams-16 | 1.8470 | 0.1240 |
| 106,272 | cpu-ic-16 | 3.3028 | 2.3618 |
| 106,272 | gpu-ic-1 | 5.8940 | 3.4049 |

Natural-order GPU IC has many triangular dependency levels and scalar host
synchronizations. These diagnostic results do not establish an IC speedup;
A100 tuning and formal timings remain pending. AMS remains the preferred
candidate for the larger HCurl cases. The selected wheel and IC binary
hashes are pinned in [a100-pending.json](results/a100-pending.json).
