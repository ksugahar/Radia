# ESRF 6 algebraic and field diagnostics: HOLD

These are new solves on hibino, not a relabeling of the earlier energy audit.
The physical problem is linear iron, mu_r=1000, response p=1, source projection
order 3, eight threads, and the existing 134686-element curved TET mesh. This is
not the original nonlinear ESRF example or current-main release acceptance.

## Identity

The installed non-editable candidate is the same unpublished 4.95.91 wheel
recorded in `candidate_de7feea0/README.md`. Wheel and mesh hashes were checked
again on hibino before execution. No installation or native payload changed.

- Wheel SHA256: `4952ba92cc0502e81949c35745108f6da499d6c6eb5bf1055a46de1be3a9ad4a`.
- Native SHA256: `5d8fe9320b6776ed5712ca7158b9878994a827d99d5725278263734dbb8a5ac8`.
- Mesh SHA256: `dfc12b84f80db1fa17fb5012b6f87c072e8aa1fb3c085b61269f047a879a54ac`.
- Diagnostic runner source: `e565507e27b5d273742b458a1410e94cfb90c312`.
- Both JSON files contain the exact runner, diagnostics, factory, solver and
  installed Python-file hashes. Source and mesh identity were checked at saves.

`bonus4_algebraic_20260912.json` and `bonus8_algebraic_20260912.json` retain
the individual samples, original residuals and block-action norms. Both runs
completed the algebraic-only mode, whose missing energy/three-engine gates
remain false. No old failed gate or historical JSON was changed.

## Observed Results

| Quantity | Assembly bonus 4 | Assembly bonus 8 |
| --- | ---: | ---: |
| DOFs | 34924 | 34924 |
| Solve seconds | 212.785 | 708.973 |
| Solve plus short diagnostics seconds | 213.018 | 709.122 |
| Iron residual L2 | 6.3288e-17 | 7.7547e-17 |
| Iron RHS L2 | 1.8940e-15 | 1.8614e-6 |
| Iron block-action relative residual | 1.4155e-13 | 1.7269e-13 |

The source projection and startup are outside the displayed solve timings.
In bonus 4 the two iron-row column-block actions have norms about 2.2356e-4
each and cancel. Dividing the 6.33e-17 residual by the nearly zero RHS produces
0.0334. This establishes the cancellation scale, not a forward-error bound.
Hodge projection orthogonality is a possible explanation for the tiny source
load, but its actual quadrature and gauge contributions have not been isolated.

`field_delta_bonus4_8_20260912.json` compares identical source/native/mesh and
physical sample identities and recomputes all averages from the saved samples:

- 360 common samples, eight per each of 45 observation centres.
- Higher-rule average-field RMS: 0.09531886 T.
- Difference RMS: 4.57187e-6 T, relative 4.79640e-5 (0.0047964 percent).
- Maximum vector difference: 1.06508e-5 T.
- Raw-sample RMS and maximum differences are also retained; the average does
  not conceal a large raw-sample difference in this run.

Two assembly rules at fixed p and fixed source projection do not establish
quadrature convergence. Evaluation-energy convergence, source projection
accuracy, nonlinear behavior, and three-engine candidate-wheel comparison
remain separate. HOLD is unchanged.

## Frozen-System Residual Correction

`bonus4_correction_20260912.json` is a separate p=1/bonus=4 re-solve with the
same wheel and mesh, using the diagnostic instrumentation at `5bbdcbe9c`.
It records its own source hashes and is not merged into the two-rule comparison.
The original solution is restored after one PARDISO solve of `A delta = b-Ax`.

- Solve: 212.112 s; solve plus diagnostics/correction: 213.241 s.
- Correction coefficient L2: 2.21217e-8 (not a gauge-invariant field error).
- Relative RMS change of the 45 averaged B vectors: 9.51857e-16.
- Maximum B-vector change: 1.33292e-16 T.
- Iron residual L2 before/after: 8.65766e-17 / 8.76428e-17. The residual
  does not decrease beyond its roundoff-scale floor in this observation.

Thus the approximately 3 percent RHS-relative block residual is not a
3 percent field change under this correction in this linear p=1 case.
This is not a proof of forward-error bounds, source-load accuracy, or
nonlinear convergence. No residual threshold or release status was changed.
