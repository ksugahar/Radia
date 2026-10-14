# Radia 5.3.1 exact-wheel acceptance

Source: `6634e31917f494628b2d4b90ee45af1fdaea07fd`, tag `v5.3.1`.
[Exact-tag native CI](https://github.com/ksugahar/Radia/actions/runs/38215376138) passed 501 tests.

- Wheel SHA-256: `d3c9bb4b76201c6ddc55556764c23c1464ca1d08927fe03ff0c520908ef88605`.
- `_radia_pybind.pyd` SHA-256: `7871a3df0172de256235cf6eeeba31bd3269d2b1305d0632c7fa2df4c1e31275`.
- All 263 packaged Python source files match the release source; fresh imports and installed RECORD hashes were checked in isolated wheel environments.

| Execution role | Strong nonlinear cases | Focused regressions | Record |
| --- | --- | --- | --- |
| Build | 6 passed | 8 passed | [acceptance](build/acceptance.json), [full6](build/full6.json), [focused](build/focused.xml) |
| Independent validation | 6 passed | 8 passed | [acceptance](validation/acceptance.json), [full6](validation/full6.json), [focused](validation/focused.xml) |

The six cases cover TET, HEX and WEDGE at orders 1 and 2, with 300000 A/m excitation. Nonlinear true residual remains at most 2e-5; linear/BH parity remains below 1e-6. Each case must exercise tangent assembly. The focused suite includes invalid/valid Jacobi diagonals, nonfinite diagonals, nonpositive curvature and the frozen WEDGE true-residual/budget regression. No skips are accepted.

The required structured two-element frozen WEDGE fixture was missing from Git. It is restored byte-for-byte under `candidate_d0d0bc4b5/validation-v5-frozen.npz`; SHA-256 `f4a679c782663e5bf0d3ead5ea0c2295517fea673915e4191aad03a5f7967758` matches the original tracked JSON record. The numerical input and tolerances were not changed. This fixture repair does not change the tagged solver or accepted wheel.

The export verifies the original runtime identity and every numerical gate before removing identifying runtime metadata. Public records retain hashes of the private originals, a shared verification-batch identifier, distinct opaque machine tokens, and canonical public-record digests. These are maintainer attestations, not remote hardware authentication. Raw records and the private batch key remain in private release storage.

Durations in the records are diagnostic provenance, not performance acceptance. Raw installation and execution logs are retained outside tracked source.

Publication remains held pending the exact-binary IH study and final release gates.
