# Heavy CI execution and release artifacts

Keep the two existing heavy runners and their `100` labels. Heavy jobs use
run-local environments. Routing and receipt verification use hosted jobs.

On a main push, each heavy distribution lane can reuse a successful PR run
only when the merged PR belongs to this repository, its recorded merge-ref
commit has the PR head as its second parent, and its Git tree equals the main
commit tree. The original run, workflow path, PR association, required branch
checks (including pinned check applications), heavy jobs and unexpired identity
artifact are verified through GitHub APIs. Classic branch protection and
ruleset required checks are combined. On unprotected branches every observed
PR-head check must still be complete and non-failing, and the lane's heavy jobs
must have executed successfully when selected. A recorded non-impact decision
requires the corresponding jobs to be skipped; that same-tree exemption also
prevents a merge from creating new heavy work absent from PR policy.
Actions attaches checks to the PR head; the
archived checkout identity separately proves the tested merge-ref tree.
Missing evidence, API permissions,
failed/pending checks, different trees and direct pushes execute the lane.
An identity artifact alone is insufficient; the original run must have
completed successfully. Every decision records the tree, commit and run ID in
the workflow summary and `<lane>-tree-identity` artifact. Reuse adds the PR,
validated commit, required checks, source artifact ID and contributing run IDs.

## Impact closure

`tools/ci_distribution_manifest.json` owns workflow projections and distribution
boundaries. Root build/configuration changes remain conservative; unknown
content selects every owner. CPU native sources/headers and CMake toolchains
remain inputs because the AMS job builds Radia, SparseSolv, axifem and PEEC.
For PRs, `tools/ci_ams_inputs.py` follows Python imports from the full AMS,
sparsesolv, `solver-numerics` and solver-heavy profiles (including inherited
profiles), resolves subprocess CLI/file loaders and package adapters, and adds a
floor of every `src/radia/*.py` module and the declared fixture directories.
Any `src/radia/` change selects the lane. Changes to package metadata,
requirements, CI install steps, workflows, CMake/toolchain files, `Build.ps1`
or `tests/test_tier_manifest.json` run the full `solver-numerics` profile on the
PR. An unknown selector or invalid input fails closed. Direct main pushes retain
the full numerical gate.

The GPU lane retains the native quote-include closure already checked by
`tests/policy/test_workflow_owner_scope.py`, package sources and shared build
files. CUDA validation directories select executable helpers, binary fixtures and
JSON inputs, including GPU release qualification JSON. Optuna consumes its independent package,
oracle/test fixtures, CUDA validation drivers and root build configuration;
it does not consume unrelated application packages. Redundant package-subtree
patterns are unnecessary underneath the package wildcard.

## Numerical audits and publication

The SparseSolv AMS workflow is also the nightly and manually dispatched full
audit. Nightly runs execute the unchanged 91-file `solver-numerics` profile,
the AMS profile and forced solver-heavy validation. Failures fail the workflow;
they are never converted into successful receipts. PRs that select the lane
run AMS, the full `solver-numerics` profile and forced solver-heavy validation.
Tree receipts record the executed profile and are reused on main only when it
matches what main would run. No numerical test is deleted.

Tree reuse is queue optimization, not artifact promotion. Before release,
dispatch the relevant distribution workflow on the exact main release SHA.
Manual dispatch always runs the full lane without tree reuse. Optuna and GPU
publication accept successful main push or main manual runs, while retaining
exact commit, successful heavy-job, wheel hash and installed-wheel acceptance
checks. A receipt-only main run has no acceptable release wheel/heavy jobs.
Never substitute the PR wheel for an exact-main artifact.

`Radia Native Release` always executes the complete numerical profile and
forces solver-heavy validation on its exact release commit. Its existing
publication/promotion checks remain mandatory. Scheduled evidence from a
different SHA cannot qualify a release. No automatic deterministic retry or
additional runner is introduced.
