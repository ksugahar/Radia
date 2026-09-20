# HCurl Eddy Bubble Production Review, 2026-09-11

## Scope and Verdict

Reviewed source baseline: `b7e89ae53` (`codex/eddy-runtime-safety`).
This review is limited to epsilon treatment, pyramid projection tolerance,
actual reduced dimensions, and the evidence needed for HDiv-MMM/SIBC coupling.
Cubit, MCP deployment, and editable-install repair are separate work.

The focused regressions pass, but they do not close production acceptance of
the complete coupled solver. No numerical defaults were changed in this audit.

## Epsilon: Still Open for Coupled Blocks

`src/radia/vim/_hcurl_tet_interaction.py` constructs analytic-moment volume
interactions without kernel epsilon. The planar volume route also reports no
kernel epsilon. These are block-local properties, not a claim about every
block of a hybrid system.

In `NgsolveHybridVIMFromHCurl`, `_eddy_hybrid.py` constructs
`HACApKSampledLaplaceInteraction` or `HACApKSampledPlanarLogInteraction` around
the analytic volume block. The sampled operator receives either an explicit
epsilon or `_default_kernel_epsilon`. Only the volume diagonal is replaced;
bridge/surface diagonals and cross interactions still use sampled kernels.
HACApK compression does not remove this regularization error.

The nonlinear iron/ESIM validation explicitly supplies `kernel_epsilon=0.05`
and `coupling_kernel_epsilon=0.05`. It cannot certify epsilon-free coupling.

Required closure: enumerate actual active block backends; independently vary
epsilon, sampling/quadrature, and H-matrix tolerance; compare coupling, loss,
force and passivity against a singular-integrated reference. Do not merely set
epsilon to zero on coincident samples. Shared/adjacent volume-surface and
surface-surface interactions need appropriate singular integration.

## Pyramid: Projection Test Exists, Operator Convergence Is Not Closed

`NgsolveHCurlCellVolumeInteraction` currently defaults to projection tolerance
`1e-4` and geometry tolerance `1e-3` for every supported family. Its docstring
explicitly motivates the former by the p=6 pyramid at degree 18. This is not a
per-family error budget, nor a bound on loss or force error.

The current tests cover:

- p=6 HEX/WEDGE projection residual below `1e-8` without subdivision.
- p=6 PYRAMID projection residual below `1e-4` without subdivision.
- p=6 PYRAMID apex subdivision level 8: 114 subtets, projection residual below
  `1e-8`, and geometry residual below `1e-12`.

All these tests passed again in this review. The last test calls the projection
helper, not a complete refined interaction and solve. It therefore does not
establish converged inductance, loss, force, or affordable build cost.

Required closure: hold the parent response space fixed, sweep projection and
outer quadrature independently, and compare actual interaction/solution
observables. Then select explicit family/geometry budgets; do not relax tests
to make a blanket default pass. Include distorted/curved and mixed-cell cases.

## Actual DoF: Strong Narrow Evidence, Not General Accuracy Certification

The existing `maglev/team28_hcurl_vim_force_live_20260910.json` records a real
mdx2 solve, not an estimated rank:

| Quantity | Recorded value |
| --- | ---: |
| Frequency | 50 Hz |
| Parent HCurl order | 6 |
| Elements | 227 |
| Parent DoF | 22,814 |
| Actual response modes | 3 |
| Retained response fraction | 0.0131498% |
| Force magnitude relative error against stored reference | 0.512924% |
| Internal H-matrix charge count | 2,043 |
| Current quadrature samples | 49,032 |

Three response coordinates do not mean the complete computation stores or
processes only three unknowns. Report response dimension, charge dimension,
quadrature size, preprocessing time, and memory separately.

The runner trains three prescribed spatial loads with `steps=1`, then solves
the volume-only interaction. Surface modes are zero. The diagnostic graph's
130 cycle modes and estimated reduction ratio are not additional solved modes
in this test. There is no live HDiv-MMM/SIBC coupled solve in this evidence.
The reference is a stored force target, not a freshly solved unreduced system
on the same mesh. Thus the force difference is not isolated reduction error.

Required closure: compare reduced and unreduced solutions of the same operator
on tractable meshes, using held-out loads/frequencies. Separately compare to
independent FEM/physical references. Sweep retained rank, h, p, quadrature and
compression independently; include corner current and Joule-loss observables,
not only total force. Neither p=6 sufficiency nor minimality is established.

## Verification Performed Now

The isolated worktree's source and native extension were explicitly selected
and their imported paths verified; installed editable packages were untouched.

Command: `python -m pytest tests/test_eddy_native_required.py
tests/test_team28_live_diagnostics.py tests/test_hcurl_eddy_bubble_cells.py
tests/test_hcurl_tet_interaction.py -q` (one command, line-wrapped here).

Result: **65 passed in 35.19 seconds**, Windows, Python 3.12.10, pytest 9.0.2;
OMP and MKL thread limits were 2. The long TEAM28 solve was not rerun.

SHA-256 provenance:

- `_radia_pybind.pyd`:
  `060e75808a1b4863aaa80976e7ca043031471f35672f834535ef9dc1d481e1f3`
- `_eddy_hybrid.py`:
  `e3d0482d40341bb435853280a80e0b1312406807233e892594cc6c0306b589be`
- `_hcurl_tet_interaction.py`:
  `d8a1baa58cc5660cfa913a34c5894158ffec9b8efc67f61e635b3efa40b8fb9d`

This is source-plus-existing-binary verification, not a fresh C++ build,
MATLAB MEX parity run, or release-wheel acceptance.

## Next Production Gate

Prioritize the coupled-block epsilon sensitivity/reference test, followed by
the pyramid interaction convergence test and same-operator reduction-error
test. Keep the test setup and error budgets separate so a good total-force
number cannot hide cancellation between regularization and reduction errors.
After those pass, validate HDiv-MMM/SIBC coupling and the MATLAB entry point on
the same input artifacts, including native lifecycle and failure behavior.
