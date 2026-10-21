# Radia 5.3.3 exact-wheel acceptance

Source: `8a98815230e8751c7ce7c21650f4fa68675718e2`, tag `v5.3.3`.
[Exact-tag native CI](https://github.com/ksugahar/Radia/actions/runs/39663524572) completed successfully.

- Wheel SHA-256: `0e2a43900f6b7331dd3bc4a2aa8ba917545e9dc56be50514de840b5e66d34f23`.
- `_radia_pybind.pyd` SHA-256: `5817b5a8ca57995c2c242bb59829e05a6e4753106680f3bb5f5dd4addc970a73`.
- All 265 packaged Python sources match the release source. Fresh imports and installed RECORD hashes were checked in isolated wheel environments.
- The CPU wheel contains no MATLAB, MEX or Simulink artifacts.

| Execution role | Strong nonlinear cases | Focused regressions | Record |
| --- | --- | --- | --- |
| Build | 6 passed | 8 passed | [acceptance](build/acceptance.json), [full6](build/full6.json), [focused](build/focused.xml) |
| Independent validation | 6 passed | 8 passed | [acceptance](validation/acceptance.json), [full6](validation/full6.json), [focused](validation/focused.xml) |

The six cases cover TET, HEX and WEDGE at orders 1 and 2 with 300000 A/m excitation. Nonlinear true residual must be at most 2e-5; linear/BH parity must be below 1e-6. Every case must exercise tangent assembly. Focused regressions cover Jacobi diagonals, nonfinite inputs, nonpositive curvature and the frozen WEDGE true-residual/budget contract. No skips are accepted. Each role used two CPU threads; no numerical gate was changed.

Each public export verifies current-machine identity and numerical gates before removing identifying runtime metadata. Records retain private-original hashes, a shared verification-batch identifier, distinct opaque machine tokens and canonical public-record digests. These are maintainer attestations, not remote hardware authentication. Originals and the batch key remain in private release storage.

Durations are diagnostic provenance, not performance acceptance. These records were executed anew with the exact-tag 5.3.3 wheel, rather than reused from an earlier release or a workflow-dispatch build.
