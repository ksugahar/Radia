# NGSolve BEM near-surface evaluation

These analytic checks exercise the NGSolve 6.2.2607 layer-potential consumer
used by RadField. They do not establish acceptance of an induction-heating
model or a Radia solver release.

```powershell
python validation_test/ngsolve_bem_nearfield/run_square_sheet.py --output private-runtime-path
python validation_test/ngsolve_bem_nearfield/run_curved_sphere.py --output private-runtime-path
```

The square sheet checks the single-layer potential against the closed-form
rectangle integral, including separation 1e-6. Its saved result checks values.
Both 6.2.2606 and 6.2.2607 passed that flat-panel case during migration review;
it does not establish an accuracy gain from upgrading. NGSolve 2607 includes
near-element quadrature and tangent corrections in `bem/potentialcf.cpp`.
Both runners own their TaskManager and fix one thread. They require only an
isolated NGSolve environment, without a Radia native extension.

The curved test uses a unit sphere with unit surface density and kernel
1/(4*pi*distance). The exact potential is 1 inside and 1/r outside; its
gradient is zero inside and -x/r^3 outside. Target patches stay strictly on
one side, at distances 0.1, 0.01, 0.001 and 1e-5 near the north pole. This is
an off-surface test; it does not test a limiting trace or jump condition.

The sphere uses order-5 geometry and constant SurfaceL2 density. Both direct
potential evaluation and target-region local expansion are evaluated, along
with `grad(potential)`. The identical gradient results do not establish that
the gradient uses the same accelerated evaluation route; no speed claim is made.

Saved `curved-sphere2607.json` records all four mesh/quadrature combinations,
including the two coarse configurations that fail the declared limits:

| maxh | bonus order | maximum potential RMS | maximum gradient RMS | passes |
|---|---:|---:|---:|---|
| 0.35 | 12 | 2.13e-6 | 3.82e-4 | no |
| 0.35 | 24 | 1.93e-7 | 1.14e-4 | no |
| 0.175 | 12 | 2.60e-7 | 9.23e-5 | yes |
| 0.175 | 24 | 1.37e-7 | 7.19e-5 | yes |

Errors are absolute RMS normalized by the unit surface potential/field scale.
The zero interior field cannot define a relative field error. Limits are
1e-6 for potential and 1e-4 for gradient; the runner's exit status requires
the explicitly selected maxh=0.175, bonus=24 case to pass. Coarse-case failures
remain visible. Geometry radial RMS improves from 2.38e-7 to 3.92e-9 with
refinement. Increasing quadrature alone does not meet the gradient limit on
the coarse geometry. These results cover this sphere and density; arbitrary
curved meshes, nonconstant densities and IH coupling require their own gates.
