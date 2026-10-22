# V-N execution receipts (C5a)

This explicitly loaded pytest plugin records execution completeness. It does not
schedule or build anything, certify physical accuracy, qualify a release, or
promote a validation record. It is not a regression baseline or lifecycle registry.

The execution producer supplies a reviewed, UTF-8 JSON plan with exactly:

| Field | Contract |
| --- | --- |
| `schema` | `radia.validation-execution-plan.v1` |
| `source_sha`, `source_tree` | Full Git object identities; checked against the actual clean checkout at entry and exit |
| `run_id`, `attempt` | Nonempty run identity and positive integer attempt |
| `profile` | `name` and SHA256 of the canonical execution projection; declared here, producer verifies the projection |
| `expected_nodeids` | Nonempty, unique, reviewed pytest nodeids under `validation_test/`; explicit parameter IDs |

Do not obtain expected IDs from the same live collection and automatically approve
them: an import skip would then disappear from both sides of the comparison.
Commission expected cases separately on an admitted environment and review the
plan against source/profile and declared parameter cases. No inventory is added
here for numerical models that have not been commissioned.

Load `validation_execution_receipt` from `tools` on PYTHONPATH and pass all four:
`--validation-plan`, `--validation-source-root`, `--validation-receipt`, and
`--validation-output-root`. Receipt and numerical output scopes must be separate,
outside the source checkout, and must not encompass it. Use a fresh receipt path
for each attempt. No options means the explicitly loaded plugin stays inactive.

The receipt schema identifier is `radia.validation-execution.v1`. `identity` is
the plan, bound by `plan_sha256`; `collected`, `executed`, `missing_collection`,
`missing_execution`, `unexpected`, and phase outcomes/durations describe actual
pytest observations. `attempted` includes call reports; `executed` excludes skipped/xfail/xpass calls.
`passed_pytest_cases` requires successful setup, call and teardown, and is not an
independent numerical-accuracy claim. Expected validation cases may run alongside ordinary source
tests in an existing tier; a failure elsewhere in that run still makes the run
fail. Unreviewed or duplicate validation cases block execution. Missing cases,
skips, xfail/xpass, collect-only, deselection, or duplicate phases are incomplete.
Assertions/collection errors, source changes, or changed source identity fail.
Neither elapsed time nor numerical thresholds are acceptance criteria here.

`state` remains initializing/collected/running and `status=incomplete` in atomic
checkpoints. Only a finished, complete session with successful phase/run results
can record `status=passed`. A timeout/kill leaves the last incomplete checkpoint;
the outer producer must mark the run interrupted/timed-out and report it, never
reuse that checkpoint as success. Final incomplete status forces nonzero exit
when pytest itself would return zero. Existing nonzero exits are preserved.

`numerical_outputs` lists relative JSON paths, byte counts and streaming SHA256
under the explicit output root. The files are unmodified observations, not
independent oracles; `none_observed` is explicit when no JSON was emitted. The
producer must declare which cases provide numerical JSON versus outcome-only
receipts. Resources, bounded child lifetime, fresh native import/source-build
provenance, Python/package versions, timeouts, safe publication/redaction, private
archive recovery and CI artifact transport remain producer responsibilities.
The plugin does not add local absolute paths or failure stack traces. The
producer must also review identifiers themselves for safe publication, disable
unreviewed plugin/retry autoload, and run one sequential phase per receipt.

C5b must wire this into the approved nightly and manual routes, preserve exact
release gates, and commission the full MCP inventory. This slice alone does not
claim that any numerical validation is scheduled or has passed.
