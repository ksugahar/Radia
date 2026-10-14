# Post-hoc one-percent core reassessment

This is a stricter replay of the recorded candidate-59b094d8 comparison,
not a new solve, not a predeclared one-percent experiment, and not acceptance
of a later rebuilt wheel. The original JSON and its three-percent threshold
remain unchanged.

Source: `three_engine_case6_bdm1_bonus12.json` in this directory.

| Pair | Core relative RMS |
| --- | --- |
| HDiv / reduced A | 0.3348911186840616 % |
| HDiv / mixed Omega | 0.6694965626678142 % |
| Reduced A / mixed Omega | 0.7818616219259033 % |

All three recorded core comparisons satisfy one percent. The current solver
contract replay independently accepts all three engines; historical HDiv
warm-start false-convergence records remain rejected regardless of their
small field discrepancy.

The original core stencil has 27 points. The full 45-point maximum is
1.407951932211685 %, which exceeds one percent and is reported, not gated.
Exclusion from the core does not establish the cause of that discrepancy or
whole-domain accuracy. This result does not validate iron-side nonlinear
observables, energy, or all excitations. Core definition was not changed.
