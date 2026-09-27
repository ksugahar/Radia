# ESRF 6 three-engine numerical acceptance

The released 5.0.1 wheel passed the declared 1% relative vector RMS gate on
the 27-point core stencil. All three nonlinear solves converged: HDiv-MMM
in 10 Newton iterations, reduced-A in 11 Anderson-accelerated Picard
iterations, and mixed total/reduced Omega in 6 Newton iterations.

| Pair (left field supplies the normalization) | Core relative RMS |
| --- | ---: |
| HDiv-MMM / reduced-A | 0.333994% |
| HDiv-MMM / mixed Omega | 0.669053% |
| reduced-A / mixed Omega | 0.780721% |

`result.json` retains all 45 field vectors for each method, both all-point
and core comparisons, convergence statistics, source settings, and mesh and
implementation hashes. The all-point comparisons are reported separately;
the 1% gate applies only to the core stencil selected before solving.
This establishes agreement between formulations on these discretizations,
not mesh convergence or absolute field accuracy.

The HDiv result was computed in this fresh campaign and reused through the
runner's checkpoint contract when restarting reduced-A with Anderson depth
3. The reduced-A and mixed Omega results were freshly computed. Mixed Omega
uses P2, Newton, material quadrature bonus 4, assembly bonus 12, source
projection order 2, and surface-flux loading. The companion
[quadrature comparison](../mixed_quadrature_20260927/README.md) checks bonus 8.

No elapsed times are accepted as performance evidence: other jobs and
profiler sampling overlapped this campaign. These results do not establish
a speed ranking. Component-wise ordering is also not an acceptance gate;
HDiv-MMM does not lie between the other two methods at every sample.
