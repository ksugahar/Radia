# HDiv validation repair, 2026-09-28

The remaining FEEC failures were reproduced against `c566acfca` on mdx2
with Python 3.12.10 and NGSolve 6.2.2606. The selected baseline had 10
failures and 25 passes; tiny reduction-order differences explain the
different count from the earlier 11-failure triage.

The repaired selection, including the existing coupled multi-iron dispatch
test, passed all 37 cases. This is a focused result, not a claim that the
complete FEEC suite, ESRF campaign, or four-machine release has passed.

## Contract corrections

- WEDGE native fixtures now supply the required field triangle rule, using
  the same NGSolve order-5 rule as the production constructor.
- Direct entries and symmetrized matvecs allow a small floating-point
  reduction-order difference instead of requiring bitwise equality.
- The full-domain and image-reduced RT2 systems retain a roundoff-scale
  comparison: 64 machine epsilons for mean magnetization and 20 for each
  scaled field component. They have different system sizes and neither is
  an exact-arithmetic oracle. Geometry and evaluator contracts are unchanged.
- Order 0 remains supported by the topology operator; order 3 is the
  unsupported-order check. Production body-fitted solve restrictions remain.
- Order-6 outer integration expects the supported symmetric rules (81
  tetrahedron points and 25 triangle points), not the retired product rule.
- The obsolete multi-iron rejection test is removed. The existing coupled
  dispatch test compares the public route with `SolveCoupled`.
- The planar material knot check now asserts its values; its former
  `or True` expression could never fail.

## HEX compression

`results_hex_compression_support_radius_20260928.json` records a diagnostic
at 2 threads with production defaults and the HEX support-radius guard on.
The 4-cubed grid has no admissible far blocks and is correctly dense.
The 6-cubed and 8-cubed grids have 40 and 728 low-rank blocks, respectively;
their stored/dense memory ratios are 0.9542 and 0.6439. The computed cube
demagnetization factors remain within 7e-7 of 1/3.

The compression gate therefore uses 8 cubed and explicitly requires the
production radius guard and release-eligible settings. It does not disable
near-interaction protection to obtain a favorable compression ratio. Timings
in the JSON are diagnostics, not an accepted comparative speed benchmark.

Reproduce the repaired selection with `python -m pytest` on these files in
this directory: `test_hdiv_hacapk_gram_performance.py`,
`test_hdiv_vim_curved_ima_roundoff.py`, `test_hdiv_vim_hex_rt1_chargegram.py`,
`test_hdiv_vim_operator.py`, `test_hdiv_vim_symmetric_outer_quad.py`,
`test_hdiv_vim_wedge_rt1_chargegram_smoke.py`, `test_soft_iron_shapes.py`,
`test_planar_materials.py`, and `test_hdiv_radsolve_dispatch.py`.

## Mass Riesz sparse factor

The connected conforming mass now uses NGSolve `SparseCholesky`, with no
PARDISO fallback. Element-local dense Cholesky blocks remain unchanged.
Duplicate COO entries are combined in input order before parallel assembly;
the exact mass cache key and batched row-major right-hand-side contract remain.
Non-positive or non-finite factors fail explicitly.

`results_hdiv_sparsecholesky_20260928.json` records 7 focused contracts and
235 extended tests passing on mdx1. They cover actual mass residuals, cache
updates, multiple right-hand sides, non-SPD rejection, nonlinear solves,
multiple materials and topology optimization. These are working-tree build
checks, not exact release-artifact acceptance or production-scale timing.
Serial nonlinear repeat solves remain bitwise identical. Parallel reduction
order may differ; the repeated-solution vector difference must be below
1e-12 relative norm, compared with a 1e-9 nonlinear solve tolerance.

## Broad regression baseline

The committed repair baseline `e6b5f9cca` passed 397 of 398 selected HDiv
checks on mdx2. The remaining check is the explicit expected failure for
HDiv pyramids in the pinned NGSolve 6.2.2606 runtime; pyramid support is
not claimed. See `results_hdiv_regression_20260928.json` for the exact list.
This broad run predates the sparse-factor replacement above; its counts
must not be presented as verification of the replacement or a release wheel.
