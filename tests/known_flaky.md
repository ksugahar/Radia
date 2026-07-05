# Known-flaky tests

Tests that intermittently fail and then pass on rerun.  The CI workflows
run with `--reruns 2 --reruns-delay 1` (`build-test.yml` "Run basic tests"
+ panel step) and `--reruns 2` (`radia-mcp-matrix.yml`), and surface reruns
in the log via `-r aR`.  So a transient flake does **not** redden the run --
but the rerun is invisible unless you read the log.

This file is the **root-cause registry** so flakiness is *tracked*, not
silently absorbed by the reruns.  (16 of the last 80 failed CI runs were in
the `Run basic tests` step; some were genuine, some were flakes that a rerun
would have saved -- the point of this registry is to tell the two apart.)

| test | first seen | symptom | suspected root cause | status |
|---|---|---|---|---|
| `test_omega_reduced_omega` | a61cecdf | intermittent assertion failure | iterative solver tolerance near a band edge | rerun-masked; not yet root-caused |
| `test_B_accuracy_inside_iron` | 2026-05-20 | field value just outside the golden band | retired moment-path dipole-in-material sampling sensitivity | rerun-masked; not yet root-caused |

## Policy

1. **A test listed here MUST also be covered by `--reruns` in the workflow
   that runs it** (otherwise a flake reddens CI for everyone).
2. **Listing a test here is a STOPGAP, not a fix.**  The goal is to
   root-cause the flakiness (tighten the tolerance, widen the golden band
   with justification, pin the RNG seed, stabilize the mesh) and then
   **delete the row**.
3. **Never add a test here to silence a REAL regression.**  A test that
   fails *deterministically* is a bug to fix, not a flake to mask.  If a
   "flake" reproduces every run locally, it is not a flake.
4. When CI fails in `Run basic tests`, check this registry FIRST: if the
   failing test is listed and the run was a single attempt (rerun would have
   saved it), it is a known flake; otherwise treat it as a real failure.

See bug pattern `flaky-test-rerun-masked-no-rootcause`
(`bug_patterns_lookup(topic="ci")`).
