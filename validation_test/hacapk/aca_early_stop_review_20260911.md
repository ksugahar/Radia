# ACA Compression-Path Investigation

## Findings

The Python/native wrapper selects ordinary ACA (`param[60]=1`), not ACA+,
in `src/ext/HACApK/cHACApK_cpp_impl.c`. The "use ACA+" comment in the leaf
filler does not describe the selected branch.

`cHACApK_aca` in `src/ext/HACApK/cHACApK_base.c` stops the whole leaf when
the largest residual entry of the selected row is at most `1e-20`. It does
not inspect a different unused row first. A zero residual row is not a
certificate that the entire residual block is zero.

The native probe `probe_aca_zero_row.c` calls that production C function on
the rank-two matrix with rows `(1,1,1)`, `(1,1,2)`, `(1,1,1)`.
Initial selection is the last row. The first rank-one update contains only
ones. Its next selected row is the first row, whose residual is zero, while
the middle row still has residual `(0,0,1)`.

Observed using freshly compiled repository C sources:

```text
rank=1 relative_frobenius_error=0.28867513459481287
exit code: 1 (intentional failed accuracy gate)
```

This proves an early-stop failure in the current source. It does not by itself
prove that every discrepancy in the spherical test takes this exact branch.
Neither the production algorithm nor the installed extension has been changed.

## Existing-Binary Spherical Cross-Block Evidence

The same volume-volume fixture as the epsilon sweep uses q=5, radius 0.5 m,
separation 3 m, epsilon 0.3 m and ACA tolerance 1e-10. Here errors are divided
by the direct *regularized* sampled sum, not the unregularized analytic value
used in the previous sweep; the slightly different reported numbers are
therefore expected.

| Leaf size | Native relative cross error |
| --- | ---: |
| 4 | 1.41173e-6 |
| 8 | 1.41173e-6 |
| 16 | -7.41533e-10 |
| 32 | -2.43457e-10 |
| 64 | -2.43457e-10 |
| 1024 (dense) | 1.53e-16 |

Scaling all coordinates, radii and epsilon together by 1e-3 or 1e3 leaves
these relative errors essentially unchanged. This does not support a simple
absolute physical-size cutoff as the main explanation for this fixture.

Rigidly rotating all sample coordinates about z, with no changes to weights
or pair distances, gives:

| Rotation (radians) | Relative cross error (leaf 8) |
| --- | ---: |
| 0 | 1.41173e-6 |
| 0.1 | 1.41173e-6 |
| 0.37 | 7.77110e-12 |

This localizes sensitivity to geometry-dependent clustering/pivot selection,
not a change to the physical kernel. The existing binary hash is
`060e75808a1b4863aaa80976e7ca043031471f35672f834535ef9dc1d481e1f3`.

## Standalone Probe Build

Run in an x64 MSVC developer environment from the repository root:

```text
cl /nologo /O2 /Gy /I src/ext/HACApK /I src/core
  validation_test/hacapk/probe_aca_zero_row.c
  validation_test/hacapk/aca_probe_unreachable.c
  src/ext/HACApK/cHACApK_base.c src/ext/HACApK/cHACApK_lib.c
  /FoC:/temp/ /FeC:/temp/probe_aca_zero_row.exe /link /OPT:REF
C:/temp/probe_aca_zero_row.exe
```

The compiler invocation above is line-wrapped for readability. This run used
Visual Studio Build Tools 18.5.1. The link-only guard file aborts if an unused
BLAS/tree/parallel entry point is unexpectedly called; it is not a numerical
substitute and must never enter a production target. Existing compiler warnings
in the shared source were not modified. No pybind11 or MEX build was performed.

## Next Fix Gate

Instrument the actual failing leaf to record rank and stop reason, then test
unused-row reselection/residual validation before accepting zero-pivot
termination. Test zero, exact-low-rank, hidden-residual, rescaled, permuted and
rectangular matrices with the actual C routine. Also examine the separate
small-update stopping condition; fixing only zero pivots does not provide a
global error guarantee.

The shared ACA kernel affects HDiv, HCurl, BEM and PEEC. A production change
must pass the original 54-case native sweep, including the strict expected
failure, plus representative existing compression regressions and build-cost
checks. Increasing leaf size or rotating the model is a diagnostic here, not
the proposed production repair.
