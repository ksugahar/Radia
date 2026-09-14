# CI candidate: ESRF6 quadrature evidence

Status: PASS for nominal ESRF6 BDM1 same-CI-wheel three-engine acceptance
and the order-2 mixed-Omega bonus 12/16 plateau. This is not release approval,
BDM2/IMA qualification, or an absolute-error certificate.

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

## Final nonlinear comparison

All three engines converged. HDiv used 8 Newton iterations, residual
1.6618566e-5 against 2e-5, one backtrack and no exhausted line search.
Reduced-A used 31 iterations (relative change 1.8592208e-5); mixed Omega
used 33 (relative B change 1.7435746e-5). FEM changes are stopping measures,
not true-error bounds. Mixed source projection and assembly use bonus 12.

Core 27-point vector RMS differences: HDiv/reduced-A 0.334891%, HDiv/mixed
0.669497%, reduced-A/mixed 0.781862%; the predeclared limit is 3%.
This acceptance uses full iron, BDM1, no IMA, nominal excitation and FEM order 2.
It does not qualify strong-saturation material accuracy from bore fields alone.

Final JSON SHA-256:
`b00ba5d7a6532721e23eefe6e8b1c0b8889b9f52e1f4240bbd2f2a99570d722f`.
Recovery archive (inputs, outputs, drivers, wheel and pip inventory):
`S:/Radia/validation_artifacts/esrf6_ci_acceptance_20260915/radia-hdiv-production-34800439509-recovery.tar.gz`.
Archive SHA-256:
`66371be72e10d1a4cd5efe1d8f78426513792b96838344f9d97342cb3a85828c`.
All 34 entries in `recovery_manifest.json` were verified after extraction.
The HDiv task owns cleanup of the completed hibino directory after this evidence
commit; no unrelated compute or shared runtime is included.
