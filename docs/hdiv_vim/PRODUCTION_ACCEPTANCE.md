# HDiv-MMM production acceptance

User-approved acceptance policy, 2026-09-11.

Three-formulation comparison is a preferred independent validation method,
not an unconditional prerequisite for every HDiv-MMM production capability.
A resource limit or unresolved defect in a comparison FEM implementation must
not be reported as a demonstrated defect in HDiv-MMM. Conversely, inability
to run FEM is not evidence that HDiv-MMM is correct.

## Acceptance by supported scope

Each release assessment names the tested element family and order, geometry
order, material model, boundary/symmetry conditions, API operations, physical
observables, and operating range. Evidence for one scope does not establish
another: a linear gap-core comparison does not certify nonlinear response,
fringe fields, periodic FFAG models, or all high-order elements.

For scopes where the comparison FEM routes are validated and executable, use
HDiv-MMM, reduced-A, and mixed total/reduced Omega on matched physical inputs.
Retain each comparison's original convergence and discrepancy gates.

Where a comparison exceeds available resources, document the failed run and
resource measurements. Require HDiv-MMM mesh/order, quadrature, and H-matrix
tolerance studies, appropriate algebraic residuals, physical invariants and
symmetry checks, and independent evidence appropriate to the claimed scope.
Independent evidence may include analytic solutions, measurements, or a
validated alternative implementation. Internal convergence alone does not
establish absolute accuracy.

Where a FEM formulation or implementation remains unvalidated, do not use its
result as an acceptance reference. Investigate it separately and identify the
independent evidence replacing that reference. If evidence is insufficient,
retain the affected scope as unvalidated rather than granting a blanket pass.

## Separate outcomes

- **HDiv-MMM defect:** a reproducible failure in the claimed capability blocks
  that capability until corrected and regression-tested.
- **Comparator defect or resource limit:** preserve the diagnostic and failed
  comparison status; assess HDiv-MMM with the alternative evidence above.
- **Validation incomplete:** identify missing evidence and limit the published
  support claim. Do not silently relabel this as either a solver defect or pass.
- **Scoped production accepted:** record the supported scope, criteria,
  evidence, provenance, limitations, and release review decision.
- **Full comparison campaign complete:** separately report completion of the
  requested seven examples and periodic FFAG comparisons. This is not implied
  by scoped production acceptance.

## Evidence and reporting

Record meshes and hashes, CAD/material/source definitions, observation points,
solver and quadrature controls, implementation/native hashes, dependency
versions, convergence diagnostics, timings, and memory use. Define acceptance
criteria before interpreting results; do not loosen thresholds to turn a
failure into a pass. Near-zero block right-hand sides require a justified,
scale-aware residual assessment, not deletion of their existing diagnostics.

Keep raw field and symmetry diagnostics alongside any projected or gap-core
metric. Gap-only mesh refinement measures sensitivity along that refinement
sequence, not a rigorous whole-model error bound. Cross-formulation agreement
and mesh quality checks likewise do not by themselves prove absolute accuracy.

This policy changes the release decision framework, not numerical results or
existing runner pass/fail semantics. It does not itself certify a release.
