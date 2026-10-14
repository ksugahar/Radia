# ESRF 5 three-engine numerical acceptance

The released 5.0.1 wheel passed the declared 1% relative vector RMS gate on
the 27-point core stencil. HDiv-MMM converged in 2 Newton iterations,
reduced-A in 30 damped Picard iterations, and mixed total/reduced Omega in
9 Newton iterations.

| Pair (left field supplies the normalization) | Core relative RMS |
| --- | ---: |
| HDiv-MMM / reduced-A | 0.224415% |
| HDiv-MMM / mixed Omega | 0.251845% |
| reduced-A / mixed Omega | 0.113532% |

`result.json` preserves all 45 field vectors per method, all-point and core
comparisons, convergence statistics, source settings, and mesh and
implementation hashes. The core is selected before solving by
`abs(x) <= 0.010 m`; all-point comparisons are reported separately.
Agreement between these discretizations does not establish mesh convergence
or absolute accuracy.

All three solves were fresh. Mixed Omega uses P2, Newton, material quadrature
bonus 4, assembly bonus 12, source projection order 2, and surface-flux loading.
Reduced-A uses relaxation 0.1 without Anderson acceleration.

Elapsed times are excluded because profiling and CI pauses affected this run.
The evidence supports no speed ranking or universal component-wise ordering
of HDiv-MMM between the other methods.

The 2026-09-29 driver defaults to `--mixed-source-load auto` and records both
the requested policy and the selected loads in its checkpoint diagnostics.
This historical result retains its explicit surface-flux setting and original
runtime identity; it has not been recomputed with the current default or 2607.
