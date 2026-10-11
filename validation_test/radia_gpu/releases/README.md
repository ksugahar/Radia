# Stable GPU release qualification

The release publisher uses the wheel from one successful `main` GPU CI run;
it does not rebuild it. CPU `radia`, `radia-mcp`, and `radia-optuna` releases
remain independent. Prereleases keep the installed 100 acceptance gate.

For a stable version, commit `<version>.json` here only after independent
review and final acceptance on SSH host `NGSolve`. This record may follow the
accepted build commit: the release gate still requires the wheel's exact source
commit to be an ancestor of `main`. Do not create placeholder acceptance records.

The qualification record uses schema `radia.gpu.stable-qualification.v1` and
contains `version`, `source_commit`, and `wheel_sha256` matching the actual
main CI wheel. Preserve the full completed review under `independent_review`:
the same `source_commit`, `outcome: accepted`, the independent `reviewer`,
the completed text in `report`, and an empty `unresolved_findings` list.
Resolve findings and review the final source before accepting it.

`a100_acceptance` contains the real installed acceptance record with schema
`radia.gpu.installed-acceptance.v1`, the same source/hash/version, `passed: true`,
`host: NGSolve`, the measured A100 `device`, at least two
`simulink_repeated_runs`, and `simulink_failure_cleanup: true`.
Set `cpu_parity_checked: true` only after testing each published CUDA solver
against its Python CPU reference. Retain these comparisons in `cpu_parity`,
keyed by capability, with measured `relative_l2` and `true_relative_residual`.
Use the solver IDs derived from the verified wheel's CUDA MEX filenames:
`csr_pcg` and, when present, `ams`. Missing or extra solver records fail;
Poisson alone cannot qualify a wheel that also publishes AMS.
The gate rejects empty results, nonfinite values, relative difference above
1e-8 or true residual above 1e-6. Retain the detailed numerical, runtime,
timing and memory evidence alongside the corresponding validation sources.

100 CI verifies LTspice; it is not required on the NGSolve host. NGSolve
starts from checked `.vol` files and does not install or run Cubit. Do not
deploy `radia-gpu` to 102. Commercial reference comparisons remain private
and never enter qualification records or release artifacts.

Stable GitHub releases omit the prerelease flag only after this qualification
passes. A stable package version does not certify unevaluated HDiv-MMM or other
capabilities; preserve the published per-feature support and limitation statements.
