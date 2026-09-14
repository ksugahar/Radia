# CI candidate: ESRF6 quadrature evidence

Status: the order-2 mixed-Omega bonus 12/16 plateau passed. The same-wheel
nonlinear three-engine run is still pending; this record does not approve a
release, BDM2, IMA, or general production acceptance.

The candidate is CI run 34800439509, source
`59b094d8ed2c19e35967fd7631f3a1473eb3d200`, Radia 4.95.91.
Wheel SHA-256:
`46e3aa1ddd89010e21419e1d28f6e44e403cf95014d7e5486109eca4232145b9`.
Native SHA-256:
`cbc9bbbc61b5cd0ee7035c642197522b2615a871e3d4c49c2887184f4f9040ec`.
The version alone does not distinguish this wheel from earlier candidates.

## Recovered results

Both files were recovered from hibino and their SHA-256 hashes were compared
with the remote originals before recording them here.

| File | SHA-256 |
| --- | --- |
| omega_algebraic_bonus12_16.json | 7cc8befa4eeb1f7ae056c816e2c36da5e07ce58bc36e99a0af7a65c8b343ca15 |
| omega_plateau_bonus12_16.json | fa1577ebf844f56e81ae68362756ae84e88301a00146fbe5b2b7506df6a99a18 |

| Metric | Observed | Predeclared limit |
| --- | ---: | ---: |
| Relative field RMS change | 2.8595794e-5 | 1e-3 |
| Maximum field change / high-rule RMS field | 5.1607991e-5 | 5e-3 |
| Relative change in Hodge harmonic norm | 9.5356089e-5 | 1e-4 |

The last metric is a relative **change**, not an absolute harmonic-norm bound.
Passing two quadrature rules is observed stability, not a rigorous integration
error bound. Projection and assembly use the same bonus within each solve.
The replay test recomputes the complete assessment from the raw JSON.

The existing hibino job directory remains owned by the HDiv acceptance task
until its running three-engine computation finishes. Cleanup follows recovery,
hash verification and commit of the final evidence, not this partial record.
