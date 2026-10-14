# Native Eddy Cross-Kernel Epsilon Sweep

## Scope

This is a native-kernel validation, not a coupled electromagnetic solve.
Uniform constant vector densities on disjoint balls and spherical shells test
volume-volume, volume-surface and surface-surface cross blocks. These fixtures
are not closed eddy-current modes or tangential SIBC modes.

For disjoint spheres the mean-value theorem gives the unregularized integral
exactly as `measure(left) * measure(right) / center_distance`. Both radii are
0.5 m and center separation is 3 m. Setting mu to 4*pi removes the kernel's
physical prefactor. The sphere measures are computed analytically, not from
the sampled weights used by the implementation.

The validator evaluates the actual C++ sampled HACApK operator. NumPy direct
summation is a validation-only oracle, not a production fallback. The sweep
holds geometry/densities fixed and varies quadrature order (3, 5), kernel
epsilon (0.3, 0.1, 0.03 m), ACA tolerance (1e-6, 1e-10), and dense-leaf versus
compressed native evaluation. The adjacent JSON records 54 evaluations,
loaded source/binary hashes and actual low-rank/dense leaf counts.

## Findings

At order 5, unregularized quadrature errors are respectively 1.49e-9, 3.99e-9,
and 6.48e-9 for the three pair types. Epsilon bias is therefore separately
observable:

| Pair | Relative epsilon shift at 0.3 m | At 0.03 m |
| --- | ---: | ---: |
| Volume-volume | -5.133e-3 | -5.175e-5 |
| Volume-surface | -5.193e-3 | -5.236e-5 |
| Surface-surface | -5.255e-3 | -5.300e-5 |

These dimensioned epsilon values are fixtures, not recommended solver
settings or a prediction of machine force/loss errors.

The native dense-leaf route agrees with direct regularized summation to
roundoff. However, the volume-volume compressed route has relative error
1.40851e-6 at ACA tolerance 1e-6 and 1.40449e-6 at tolerance 1e-10
(order 5, epsilon 0.3 m). Dense-leaf error is 1.52e-16. This localizes the
discrepancy to the compressed route but does not yet identify the ACA stopping
or pivoting mechanism responsible.

`tests/test_eddy_cross_epsilon.py` preserves this as a strict expected failure
against an explicit 1e-8 cross-block accuracy budget. That budget is a
validation target, not an assertion that ACA tolerance directly bounds every
matrix entry. This expected failure is an OPEN issue, not production success.
It must be removed once the native fix is independently verified.

## HDiv Pyramid Gate

The existing functional probe was rerun on installed NGSolve 6.2.2606.
The pyramid mesh passed H1/volume checks, but HDiv mass assembly raised
`HDivHighOrderFESpace: Pyramid elements not implemented yet!`.
HCurl pyramid projection support does not enable HDiv pyramid or its Radia
charge-Gram kernel. Keep this distinction when selecting coupled meshes.

## Remaining Acceptance

First diagnose the compressed-route discrepancy without changing kernel
epsilon to compensate. The current evidence uses the existing binary whose
SHA-256 is recorded, not a newly compiled C++ or MEX artifact.

Next cover touching/self-panel singular interactions and actual admissible
HCurl/BDM/SIBC modes; then compare coupled current, loss, force and passivity
with independent references. The present separated-support test must not be
used to approve those untested cases.
