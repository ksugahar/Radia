# Radia 5.3.2 exact-wheel acceptance

Source: `0bdd66c9cf8c94410f14dc7fa5fe81e5b04c6619`, tag `v5.3.2`.
[Exact-tag native CI](https://github.com/ksugahar/Radia/actions/runs/38262424368) passed 512 tests without skips.

- Wheel SHA-256: `bb5337ea9aacf8f720634c5c1a3b8ea1722a53869bff6f26994e7eeb33eecd0a`.
- `_radia_pybind.pyd` SHA-256: `9020a4f700d381fe5388870ef3b77a43b67506f9f1610f9ef8b94906bdfc4996`.
- All 263 packaged Python sources match the release source. Fresh imports and installed RECORD hashes were checked in isolated wheel environments.
- The CPU wheel contains no MATLAB, MEX or Simulink artifacts.

| Execution role | Strong nonlinear cases | Focused regressions | Record |
| --- | --- | --- | --- |
| Build | 6 passed | 8 passed | [acceptance](build/acceptance.json), [full6](build/full6.json), [focused](build/focused.xml) |
| Independent validation | 6 passed | 8 passed | [acceptance](validation/acceptance.json), [full6](validation/full6.json), [focused](validation/focused.xml) |

The six cases cover TET, HEX and WEDGE at orders 1 and 2 with 300000 A/m excitation. Nonlinear true residual must be at most 2e-5; linear/BH parity must be below 1e-6. Every case must exercise tangent assembly. Focused regressions cover Jacobi diagonals, nonfinite inputs, nonpositive curvature and the frozen WEDGE true-residual/budget contract. No skips are accepted.

Each public export verifies current-machine identity and numerical gates before removing identifying runtime metadata. Records retain private-original hashes, a shared verification-batch identifier, distinct opaque machine tokens and canonical public-record digests. These are maintainer attestations, not remote hardware authentication. Originals and the batch key remain in private release storage.

Durations are diagnostic provenance, not performance acceptance. Numerical formulations are unchanged from 5.3.1; this acceptance was executed anew on the exact 5.3.2 wheel.
