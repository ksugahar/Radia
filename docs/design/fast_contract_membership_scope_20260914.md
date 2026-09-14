# Fast-contract profile membership scope

## Failure and cause

Main commit `1cf6c724bdebfa3585f40e0c08dbd54b51aae8e5`, Actions run
34772319443, passed 518 tests (one existing NGSolve-module skip), but exceeded
the unchanged 60-second tier budget: 67.56 seconds. Its preceding PR run
34772128198 passed the same count in 59.14 seconds, leaving little margin.
The first BH/PCHIP test took 11.66 seconds on main versus 6.64 seconds on the PR;
this includes first-use SciPy work, not a separately measured import benchmark.

The expansion to 31 files was not caused directly by editing radia-fast.yml.
Adding a BH test to fast-contracts.paths changed the profile structure.
`changed_impact_tests` therefore returned unknown and selected every registered
impact, including unrelated solver, build, release-packaging and subprocess tests.

## Change

- A profile paths-only change selects added and removed memberships, together
  with ordinary changed-source/test impacts and changed impact-rule owners.
  The unchanged baseline profile still runs. Reordering has no membership delta.
- Every other profile key is compared. Budget, inheritance, description, schema,
  unknown keys, added/removed profiles, invalid paths and unavailable comparison
  evidence retain conservative broad selection.
- A retired test absent both from the checkout and all current manifest
  references is reported and not handed to pytest. Still-referenced or newly
  misspelled missing paths fail configuration; renamed tests execute at the new
  path. This is not pytest skipping an existing test.
- BH is impact-selected rather than run for unrelated changes. It is explicitly
  retained in native-smoke, and selected by the actual helper, production
  `_nonlinear.py`, shipped BH table, native fixture, audit runner, tracked
  qualified result, the BH test itself, and fast-CI workflow changes.
- radia-fast.yml changes explicitly select the BH dependency check and the
  shared CI, profile-policy and selector contracts. SciPy remains declared in
  that existing isolated environment. No new lane, budget increase, solver
  change, numerical fixture edit or test duplication is introduced.

## Same-condition replay

Both runs used LAB, the same standalone minimal dependency venv
`C:/temp/bh-fast-ci-minimal-20260914`, the same checkout directory, comparison
base `6f447f4df388d38a1b796021440e1bd8fdaefc52`, and Git HEAD
`1cf6c724bdebfa3585f40e0c08dbd54b51aae8e5`. The after run applied the proposed
working-tree selector/manifest/regression changes; HEAD and changed-file query
were held fixed. Each invocation started a fresh pytest process. OS caches and
host load were not controlled, so this is not a per-test speedup benchmark.

| Replay | Selected files | Passed | Existing module skips | Tier elapsed |
|---|---:|---:|---:|---:|
| Before | 31 | 518 | 1 | 29.55 s |
| After | 13 | 167 | 1 | 5.34 s |

After includes 22 new selector cases at that measurement point. One additional
missing-new-member regression was subsequently added (78 focused tests PASS).
All eight BH contracts ran in both measurements. The existing NGSolve-only
module skip is unchanged. The before local time differs from CI; it does not
reproduce the CI host's cold-start timing, only its selection and test results.

Before JUnit totals include 5.56s for Simulink release packaging, 5.46s for
sparsesolv MATLAB contracts, and 2.28s for mdx preflight contracts. Those groups
are excluded from this unrelated change but retain their existing source/test
impact rules. Structural profile semantics changes still select broadly.

Command in both runs (only JUnit/output-log filenames differ):

```powershell
C:/temp/bh-fast-ci-minimal-20260914/Scripts/python.exe tools/run_test_tier.py --profile fast-contracts --since 6f447f4df388d38a1b796021440e1bd8fdaefc52 --junitxml C:/temp/fast-budget-before.xml
```

Raw logs and JUnit evidence are retained at
`S:/Radia/validation_artifacts/fast_contract_budget_20260914/`.
Formal acceptance still requires the new PR's exact-SHA CI result. Main merge
and release remain the management task's responsibility.
