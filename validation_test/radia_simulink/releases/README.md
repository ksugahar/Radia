# Stable GPU release qualification

The release publisher uses the wheel from one successful `main` GPU CI run;
it does not rebuild it. CPU `radia`, `radia-mcp`, and `radia-optuna` releases
remain independent. Every release, prerelease or stable, requires one installed
acceptance record (`radia.gpu.installed-acceptance.v1`) for the exact wheel in
the accepted artifact: `passed: true`, `ltspice_checked: true`, the same source,
hash, payload manifest and version, at least two `simulink_repeated_runs` and
`simulink_failure_cleanup: true`. A prerelease (`aN`, `bN`, `rcN`) needs only
that record; a stable version additionally needs the qualification record below.

For a stable version, commit `<version>.json` here only after independent
review and final installed acceptance on the Windows GPU validation host
(NVIDIA A100-SXM4-40GB). This record may follow the
accepted build commit: the release gate still requires the wheel's exact source
commit to be an ancestor of `main`. Do not create placeholder acceptance records.

The qualification record uses schema `radia.gpu.stable-qualification.v1` and
contains `version`, `source_commit`, and `wheel_sha256` matching the actual
main CI wheel. Preserve the full completed review under `independent_review`:
the same `source_commit`, `outcome: accepted`, the independent `reviewer`,
the completed text in `report`, and an empty `unresolved_findings` list.
Resolve findings and review the final source before accepting it.

`windows_gpu_acceptance` contains the real installed acceptance record with schema
`radia.gpu.installed-acceptance.v1`, the same source/hash/version, `passed: true`,
the payload `manifest_sha256` matching the verified wheel,
`ltspice_checked: true`, `execution_role: windows-gpu`, the measured
`device: NVIDIA A100-SXM4-40GB` (whitespace-normalized exact match), at least two
`simulink_repeated_runs`, and `simulink_failure_cleanup: true`.
`run_installed.py` records the actual OS hostname and MATLAB GPU device,
compute capability and total memory. Verify the validation-host identity in the raw
operational evidence retained outside Git. Public qualification records retain
the execution role and measured GPU fields, omitting identifying OS hostnames.
Set `cpu_parity_checked: true` only after testing each published CUDA solver
against its Python CPU reference. Retain these comparisons in `cpu_parity`,
keyed by capability, with measured `relative_l2` and `true_relative_residual`.
Use the solver IDs derived from the verified wheel's CUDA MEX filenames:
`csr_pcg` and, when present, `ams`. Missing or extra solver records fail;
Poisson alone cannot qualify a wheel that also publishes AMS.
The gate rejects empty results, nonfinite values, relative difference above
1e-8 or true residual above 1e-6. Retain the detailed numerical, runtime,
timing and memory evidence alongside the corresponding validation sources.

The A100 installed acceptance is the release gate; the CPU build host compiles no CUDA and runs no GPU acceptance. Records 0.1.0 to 0.1.4 in this directory are historical: they were accepted under the earlier build-host gate and keep their original fields unchanged. NGSolve
starts from checked `.vol` files and does not install or run Cubit. Do not
deploy `radia-gpu` to 102. Commercial reference comparisons remain private
and never enter qualification records or release artifacts.

Stable GitHub releases omit the prerelease flag only after this qualification
passes. A stable package version does not certify unevaluated HDiv-MMM or other
capabilities; preserve the published per-feature support and limitation statements.
