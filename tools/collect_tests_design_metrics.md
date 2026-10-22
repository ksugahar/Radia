# Private CI observation collector

`python tools/collect_tests_design_metrics.py --repo OWNER/REPO --merged-from UTC --merged-until UTC --out ABSOLUTE_PRIVATE_PATH`

Use an authenticated manager-owned `gh` context. The collector invokes only
`gh api --method GET`, never reads or prints credentials, and never writes to
GitHub. Do not put empirical snapshots in a public repository. The output must
be a fresh absolute directory outside Git worktrees, without junctions/symlinks.
Raw API bodies may contain private actors/runner details: review before sharing.

The merge window is UTC, lower-inclusive/upper-exclusive. Closed PRs and runs
are fully paginated, without search API queries (and their 1,000-result limit).
This first implementation scans repository history and can be expensive: it
has no resume mode. Use a fresh output directory on a later manager-approved
collection, preserving partial results. Independent pages can race new runs or
reruns; repeated IDs, missing attempts and count mismatches make collection
incomplete, with nonzero CLI exit. Request timeouts are bounded; no retry.

`summary.json` links every normalized observation to retained response endpoints
and hashes. `collection_complete` means retrieval checks passed, not that all
historical PR head associations were recoverable. Links use API PR associations
or exact final-head/merge SHA only; absent earlier links remain unresolved, never
inferred from branch names. PR/push events and attempt IDs remain separate.

Open-to-merge is elapsed calendar time, including review/human delay. Queue time
is reported only for attempt 1: the original run creation timestamp is not a
reliable rerun enqueue timestamp. `updated_at` is retained but never treated as
completion. Run wall time remains null. Per-job durations and their observed sum
are resource-time observations, not DAG critical path or user waiting time.
Unknown durations remain null. Cancelled jobs can have partial execution.
No causal throughput improvement, regression catch, lifecycle promotion or
release qualification is inferred. Catch/mutation/coverage joins and C5 receipt
reading are separate slices; this tool activates no schedule or calculation.

Contract tests use synthetic data, no token or API. Manager performs real private
collections. GitHub API contracts: https://docs.github.com/en/rest/actions/workflow-runs
and https://docs.github.com/en/rest/actions/workflow-jobs.


Job duration quality is independent of API retrieval completeness. Negative,
malformed or missing job intervals produce `duration_s: null` and per-job
`duration_unavailable_reason` (`negative_interval`, `invalid_timestamp` or
`missing_timestamp`). They increment unknown-duration counts and contribute no
value to the explicitly observed sum; they are never clipped to zero or repaired.
Equal valid timestamps are an observed zero interval, with no invented rounding
correction. These quality gaps alone do not set `collection_complete=false`;
API/page/identity failures and invalid PR/run intervals still do. No CI/physics
verdict is changed by the duration-quality metadata.
