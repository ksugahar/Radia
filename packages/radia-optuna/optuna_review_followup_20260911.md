# Opus review follow-up: fail-closed evidence attribution

Baseline: `01297fc0b`, branch `codex/optuna-5-oracle`. This is the first
safety correction, not completion of the 15 findings or a release approval.
Existing GOT/Simulink and Optuna 5 numerical implementations are retained.
No optimizer production source is changed.

## Resolved or contained in this change

- Terminal-name matching is removed. Only explicit `optuna.*` attribute paths
  count as identified source references. `server.stop`, a `"stop"` string or
  another object's `params` cannot identify an Optuna API. Unit tests cover
  unrelated and correctly qualified receivers. This deliberately does not
  infer types for local variables or dynamic dispatch.
- Consequently, required coverage is now 79 identified mappings and 322
  asserted mappings out of 401; full compatibility is false. The 322 are not
  newly absent implementations. Release health must remain blocked until
  owner-aware evidence is supplied. Static reference attribution itself is
  not execution tracing or proof of output assertions.
- MCP oracle audit checks the normalized oracle digest, with a stale-digest
  regression. Health/release tests check actual incomplete coverage as well
  as synthetic complete evidence; mocked success is never production evidence.
- Standalone CI triggers include fixture, Optuna MATLAB test and precision
  helper changes in both push and pull-request filters.
- Test-manifest regeneration and canonical LF-normalized coverage serialization
  are checked in normal package tests. Digest comparisons use the same LF
  normalization as the generator.
- Private helpers are excluded and @Class/external-method discovery is handled
  without treating method filenames as top-level symbols. Synthetic layout
  regression covers both. The actual nsgaii namespace is recorded as a
  resolvable package rather than the former flattened label.
- MATLAB metadata checks cover all recorded classes/members, public method
  access, public property read access and enum constants. A missing class is
  reported with verifyNotEmpty and skipped locally, not an assertion that
  aborts the remaining audit.
- Precision is rerun in `results_optuna_precision_lab_20260911.json`. It records
  the exact UTF-8 test-source SHA256, so its lines are tied to that source.
  The September 10 JSON remains historical evidence for c90ab01cb.

## Open work, not waived

1. Supply owner-aware instance/dynamic contracts for the 322 asserted required
   entries, starting with Study/Trial operations. Do not restore suffix matches
   or hand-assert complete coverage to get a green release gate.
2. Complete package-qualified and continuation-aware inheritance handling.
   Current basename ambiguity fails closed instead of certifying an ambiguous
   public mapping.
3. Remove the generator's dependence on the executing Python version for
   inherited built-in member classification. Same-runtime determinism is not
   cross-version determinism.
4. Rebase and finish the separate constructor-default audit after its owner
   commits the seven-file WIP. `C:/temp/radia-optuna-ctor-audit` was inspected
   read-only and remains untouched.

Acceptance for this phase: 76 MATLAB upstream-suite tests plus precision
diagnostic self-test; 13 radia-optuna package tests; 8 focused MCP quality/policy
tests and the MCP compatibility/audit regression. These passing tests establish
honest incomplete reporting and unchanged tested numeric behavior, **not** a
passed release gate. No merge, push, tag, PyPI or release-quad operation occurred.
