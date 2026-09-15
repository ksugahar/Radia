# Second MCP review reconciliation

Review baseline: main `afb8c621c`, radia-mcp 1.4.54. This records the supplied
review summary, not an independent reproduction of every finding in the full
external report. No shared checkout, editable installation or live client is
changed by this maintenance patch.

| Finding | Disposition and evidence |
| --- | --- |
| Non-Windows `RADIA_MCP_TEMP` ignored | Fixed both native journal and headless batch paths through one helper. Mocked Windows/Linux/macOS cases check configured roots and platform defaults. The native journal suite passed twice consecutively (20 cases each); production Cubit was not started. Existing journal overwrite refusal remains in place. |
| README `49 servers` | Replaced current fixed counts with catalog discovery; historical changelog counts are not rewritten. |
| Missing 1.4.54 changes | Added SDK floor, argument binding/validation and restricted reload entries. New scratch-path changes are under Unreleased, not falsely described as published. |
| Blanket aliases and CBCR | Naming guidance explicitly excludes retired GUI tools; CBCR carries a historical-design notice and no longer claims current production status. Supported headless checkpoint operations remain available. |
| NumPy in minimal tests | Clarified runtime versus test dependencies. Mandatory runtime remains MCP; scoped CSV tests install NumPy, while minimum-SDK registration/dispatch does not add it. This does not establish all tools work without their optional dependencies. |
| 29 package tests with monorepo dependencies (review count) | Boundary refactored against the current main, rather than matching the old count. Repository-dependent tests now run in `tests/mcp_integration`, hybrid modules retain their unit tests, and four package-source paths no longer traverse the monorepo. The separate CI selection and package-only copy probe protect this boundary. See the scope and outstanding native evidence below. |
| Old S: source | An old checkout alone does not identify any live client. Use the original client's provenance; do not pull/repoint the shared tree based on its name. Integration and live acceptance remain separate. |

## Package-test boundary follow-up

- The current diff relocates 36 existing test functions: 35 repository
  contracts and one native-acceptance check. The SDK metadata/workflow test
  is split rather than removed. The duplicated workflow-discovery assertion
  is consolidated into `tests/test_ci_execution_policy.py`. Parameterized cases and new CI-boundary
  regressions explain why this is not the review's historical count of 29.
- Repository contracts cover saved field/motor/HDiv evidence, MATLAB/Optuna
  source manifests, repository MSH assets, Cubit golden scripts, the legacy
  md2html CLI, and CI policy. Missing required fixtures and Git inventory now
  fail instead of silently skipping.
- SDK behavior, MATLAB code generation, policy-string comparisons, Cubit
  session behavior, md2html conversion and pure force helpers remain package
  tests. A copy containing only the package executed 90 focused cases on LAB;
  repository integration executed 44 cases with no skips. These are not a
  full-package audit or numerical acceptance claim.
- The native motor-angle source-freshness test was preserved in
  `validation_test/radia_mcp/test_motor_angle_source_freshness.py`, unchanged
  in its requirement that recorded hashes match the candidate. It currently
  fails: the 2026-09-02 MATLAB evidence binds a different `radia_mex.cpp` hash.
  New native/MATLAB validation evidence is required before claiming current
  native acceptance. No evidence hash, threshold or native source was changed.
- CI runs the integration lane independently of the package's optional-import
  collection filter, then probes the focused package-only copy. A test migration
  is not permission to drop its failure signal.
- The Linux package-only probe exposed an unnecessary import-time NGSolve
  dependency in pure force helpers. NGSolve imports now occur only inside FEM
  operations; formulas and signatures are unchanged. A cold-process test blocks
  NGSolve explicitly and verifies pure helpers work while FEM calls fail loudly.
- The surface-load unit test now uses the independent analytic pressure
  `B^2/(2*mu)` as its expected value, instead of calling the separate Radia
  distribution's pressure adapter. That adapter still requires Radia; its
  misleading dependency-free docstring is corrected, not its behavior.

This patch closes the identified source-tree coupling, not package-wide
completion, a release, or the outstanding native-evidence refresh.

## Native follow-up

The 2026-09-15 clean MEX rebuild and dedicated LAB Engine run completed 87
tests (84 passed, 3 failed, none incomplete). The failure artifact is preserved
separately; historical passing evidence was not restamped. Two HEX directional
derivative comparisons require solver-side investigation; the application
failure-message check also needs its Python native runtime in the isolated
validation environment. See
`validation_test/radia_mcp/native_evidence_refresh_20260915.md` for exact
build identity, failures and execution boundaries. Native acceptance remains
open, rather than being hidden by the MCP test-boundary cleanup.
