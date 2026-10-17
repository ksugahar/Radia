# Solver-heavy validated-source gate

Normal compact contracts remain ordinary `tests/` CI. The CPU AMS lane keeps
its native build, AMS tests and `solver-numerics` tests mandatory when its
existing impact selector requires them. Only the separate `solver-heavy`
validation step gains verified reuse. This does **not** make the remaining
native build or numerical suite cheap and does not waive release physics gates.

For solver-heavy, the gate either executes all required cases without skips or
verifies a prior PASS for exactly the same conservative input fingerprint:

- CPU sources, validation scripts, saved fixtures and inputs; production package
  sources; numerical test helpers; root build/dependency declarations;
- the relevant validation profile and inheritance, the gate and tier runner;
- actual Python/platform/processor, installed distribution versions, MSVC and
  Windows SDK identity, and fixed one-thread/non-partial environment contract.

Source and input file bytes are hashed. Relevant dirty/missing/linked inputs or
unknown compiler/SDK/environment facts disable reuse and execute validation.
Workflow routing and unrelated CI-only tests are not numerical inputs. The
closure is intentionally conservative: unrelated validation evidence changes
can cause execution; untracked runtime inputs are not claimed to be covered.
Native builds and their source/provenance checks continue to protect current
artifacts. Reusing source evidence never certifies a new wheel or native binary.

## Trust and failure behavior

A receipt is accepted only with authoritative GitHub metadata identifying a
successful, completed `main` push/manual run of the same CPU workflow. Its
commit must be an ancestor of fetched `origin/main`, its run/attempt and source
SHA must agree with the receipt, and its live artifact must be unexpired and at
most 30 days old. The API archive SHA-256, exact archive entries, evidence hash,
profile fingerprint, case count and no-failure/no-skip JUnit checks must agree.
PR/fork receipts, self-asserted JSON without this provenance, and receipts from
other workflows cannot authorize reuse. Archive redirects do not receive the
GitHub authorization token outside the GitHub API host.

Any cache miss, unknown condition, unavailable token/API, corrupted evidence or
changed fingerprint executes the required validation. Failed execution remains
failed. A verified reuse copies the original receipt and JUnit bytes unchanged
and saves a separate `verification.json` with the current source SHA, time and
original provenance. It does not relabel old results as current calculations.
Only actual execution uploads a new reusable receipt; reuse cannot manufacture
a newer authoritative proof or extend the original evidence's lifetime.

Manual workflow dispatch defaults `force_validation=true`, bypassing reuse.
The CLI also supports `--force`. Release/native workflow gates are unchanged.
No runner, service, label, resource-budget or currently running CI change is
performed by this implementation.

## Acceptance and remaining checks

Focused tests cover changed source/mesh/build/dependency conditions, routing-only
changes, missing/failed/expired/tampered/wrong-profile/wrong-run proofs, main
ancestry, immutable original identity, forced/missed/unknown execution and failed
validation. These are offline contract tests, not a claim of a live GitHub cache
hit. First deployment must execute on main to create an eligible receipt; a
later unchanged-source run must demonstrate a verified hit with the same
fingerprint. Confirm actual CMake compiler/SDK extraction and GitHub v4 artifact
digest metadata in that runner context. Unknown values safely execute instead.

Before integration with the capability-routing PR, merge its reviewed main,
retain all workflow projections, add this gate to the light lane's reviewed
dynamic-loader closure (its tests load this file), and repeat the light F2 audit
before updating that closure digest. Do not automatically regenerate signatures
as a substitute for reviewing execution changes.


Untracked production/validation inputs, including ignored inputs, prevent reuse;
ordinary generated compiler and Python cache outputs do not certify a new
artifact. The known `_radia_pybind.pyd.build.json` sidecar is accepted only after
checking its clean current-commit provenance, binary name, length and SHA-256.
Package test fixtures are included in the source fingerprint. This is separate
from the mandatory native build and solver-numerics checks; a matching source
receipt never certifies a newly built binary. Live main cache-hit and cache-miss
acceptance remains required before declaring the reuse feature operational.
