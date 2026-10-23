# Component/detector incidence (manual, private output)

`python tools/coverage_map.py --request <private-request.json> --out <new-absolute-private-map.json>`

This manual CLI reads existing A4 `radia.coverage-evidence.v1` exports; it runs
no solver, detector, GitHub API or coverage collection. No adoption/accuracy or
retirement decision is made. Actual inputs/maps belong outside Git. Output is
exclusive-create and rejected inside a Git worktree. Links/junctions, reference
traversal, hash mismatch, duplicate JSON keys and incompatible identities fail.

## Explicit input projections

Request fields (all required): `schema: radia.coverage-map.request.v1`,
`source_sha` (40 lowercase hex), `artifacts` (name to SHA256), `inventory`
(reference), and `detectors` (array). References have exactly `path` (normalized
relative to the request directory) and `sha256`. Inventory, test source, capture
and attribution bytes are all hash-checked. These are reader projections, NOT
a new canonical adoption-record or mutation schema; their reviewed declarations
still require producer/owner review. No real adoption IDs are invented.

Inventory fields: `schema: radia.coverage-map.inventory.v1`, `source_sha`,
nonempty `components`. Each component declares `id`, `kind`, and a nonempty
unique `expected_detector_ids` reporting denominator. Python components also
declare production `path` and `sha256` (the A4 captured source bytes):

| kind | Extra fields | Observation |
|---|---|---|
| python_module | none | One or more executed module lines |
| python_line | positive `line` | That line observed |
| python_branch | `arc: [from,to]` | That branch observed with branch tracing |
| python_function | `qualified_name` | UNKNOWN: A4 v1 does not attest function entry |
| native_entry | `entry`, `artifact`, `configuration` | Positive invocation of same artifact/configuration |
| native_internal | `entry`, `artifact`, `configuration` | UNKNOWN: entry is not internal coverage |

Arcs use A4 conventions, including negative nonzero entry/exit lines. A missing
file/capture is UNKNOWN, not zero. File hashes must match. No AST-span inference
is made from definition/import line hits. The inventory is explicit and may be
partial: it is not a discovery of all Radia APIs or executable code.

Each detector has `id` (canonical E2E identity supplied by the owner), `nodeid`,
`test_source` (reference), `metadata` (purpose/quantity/model/reference/declared
route where available), and `captures`. Each capture has `export` and
`attribution` references. Attribution fields, all required:

* `schema: radia.coverage-map.attribution.v1`, `source_sha`, `export_sha256`,
  `test_source_sha256` bind exact inputs;
* `detector_ids` and corresponding `nodeids` identify the capture's test scope;
* Boolean `reviewed` and `capture_complete` record owner-reviewed attribution.

Only a reviewed, single-detector attribution assigns hits. An aggregate export
cannot lend all its hits to each test. Unreviewed/missing bindings remain UNKNOWN.
Do not set `reviewed` just to obtain a non-UNKNOWN map: the supplying reviewer
must establish the identity, nodeid/test source, capture scope and duplicate/alias
resolution. Different runs, shards, heavy/light levels and repetitions of one
detector belong in its `captures` array, not invented detector IDs. Fully
duplicated aliases must be resolved by producer/reviewer before this reader.

All captures use the request's exact source and artifact identities; differing
routes/models may remain distinct detectors but incompatible source or native
versions cannot be added together. Native identity/configuration mismatch blocks.
A4's existing heavy/light `compare_exports` remains separately necessary when
qualifying a surrogate path; this map does not replace that comparison.

## Interpretation

The output includes both component→detector IDs and detector→component IDs.
`observed_test_count` is a lower bound, deduplicated by detector ID; metadata is
retained for independence review. Complete observations give categories 0/1/2+.
An authoritative zero requires every expected detector to have an applicable,
reviewed complete capture, passed original gate and adequate tracing. Exclusions,
unverified required subprocess capture, missing branch tracing or missing scope
prevent completeness. Known hits in incomplete captures remain lower bounds,
with `category: UNKNOWN` and the separate `observed_category` visible.

A-B-C and A-D-C give A/C two path detectors, B/D one, and inventoried E zero
only with complete no-hit captures from its expected detectors. Incidence does
not reconstruct chronological call order or prove a detector kills a mutation.
No statistical 90/99-percent reliability is inferred. Overlap is welcome and is
never retirement evidence. `qualified` and `automatic_retirement` are always
false. Mutation corroboration awaits the canonical schema; function entry,
native internals, history and management activity stay explicitly UNKNOWN in
this first slice. No thresholds, baselines or adoption criteria are changed.
