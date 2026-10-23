# Historical validation records

This directory holds one-time evidence, not a numerical release lane. Package-owned records belong in `packages/<package>/validation_records/`.

Each record has `record.json`, a hash-checked result JSON and a local entry script with SHA256. Use `tools/validation_record.schema.json`. Preserve original source, neutral runtime role, units, limits and timings. Keep failures and never rewrite raw results or claim stale evidence as current.

The stdlib checker supports only the fixed keywords used by this repository-owned schema and rejects unsupported keywords; it is not a general JSON Schema validator.

Run `python tools/check_validation_records.py --report <owned-output>/records.json` for read-only schema, path, hash, covered-content staleness checks. `repo-contracts` runs the same gate without executing recorded physics. Privacy patterns are a guard, not publication clearance; source/data still require confidentiality review. Creating this directory authorizes no record migration.

Each record is a development milestone using the latest `radia.validation-record.v3` schema. Required metadata includes `milestone`, a one-line `claim`, calendar `date` (YYYY-MM-DD), `radia_version`, `source.provenance`, `source.repo_sha`, and hash-checked `result` and entry artifacts. `conditions` declares the problem, mesh and neutral hardware description; `environment.threads` records the thread count. `comparison` names the reference and its independence (including any limitation, such as a same-solver regression lock). This metadata does not establish oracle independence by itself.

Champion categories are `accuracy-champion`, `speed` and `scale`; the other
fixed roles are `adoption` and `milestone` (no champion claim). `supersedes` lists prior same-distribution record IDs, a relation
(`improves`, `corrects`, or `replaces`) and its reason; use an empty list for a
first record. Targets must exist and relations must be acyclic. This preserves
the earlier evidence. It does not automatically rank performance, prove that
conditions match, or establish an unrestricted speed record.

## Adoption: criterion, recorded decision and protection

An optional `adoption` block joins the production `feature`, `criterion`, user
`approval`, recorded `result` and `protected_by` tests paths. No block means no
adoption claim. A valid FAIL or UNRESOLVED can be preserved; only PASS with
current status yields `adoption_claim=true`. The checker does not authenticate
the user: actual consent and its scope require independent review. Public
approval metadata contains only the approval date and identifier; keep private
consent locations outside the public record. Never invent or backdate approval.

The criterion includes reference type/citation, a quantity definition and unit,
convergence variables/explicit levels, and named series that hold other
variables fixed. Each series names its axis/direction, strict error decrease,
rate model, approved minimum rate, finest-level threshold/operator and expected
theory (norm, source, assumptions, limitations). Record an explicit failure
definition, the true-residual threshold and `saturation_policy="blocks-pass"`.
An empirical lower bound must not be called a theorem. No ordinary point RMS
is automatically an energy norm. An unexplained plateau is not an exception.
At least two levels are needed for a rate; any stronger level-count requirement
belongs in the user-approved criterion. A p-then-h sequence uses separate series
with shared endpoints, not a rate fitted across simultaneously changing axes.

Analytic references require a DOI and a hash-checked `symbolic_verification`
artifact (`path`, `sha256`, relative to the record). An independent route requires
an `independence` rationale and `shared_components`, even when the list is empty.
Measurements and benchmarks require an uncertainty statement. The checker checks
presence and bytes, not the scientific adequacy of these statements or DOI
resolution; a reviewer must verify the reference and symbolic calculation.

Approval contains `approved_by="user"`, calendar `approved_on`, public
`approval_ref` and the SHA256 returned by
`tools.validation_adoption.criterion_hash(criterion)`. This binds provenance to
the whole proposed criterion, not just its finest error threshold. Changing a
criterion needs renewed actual approval, even if someone can recompute its hash.

`adoption.result` has the exact execution source SHA, the enclosing result
artifact SHA256, a verdict, and every level's parameters, finite nonnegative
error, true relative residual, mesh SHA256 and descriptive runtime. The result
JSON contains an `adoption_result` object equal to this declaration except for
`artifact_sha256` (a hash cannot include itself). All execution source identities
must agree. Pre-rewrite/source-null evidence cannot create an exact-source
adoption certificate. Older source objects need not exist in the local Git
history. Source, artifact, mesh and parameter identity are provenance, not proof
that the described computation was actually executed; review its originals.

For adjacent errors E1,E2, `log-ratio` measures
`(log(E1)-log(E2))/abs(log(x1)-log(x2))`; `log-per-unit` uses
`(log(E1)-log(E2))/abs(x2-x1)`. The axis direction is explicit and other variables
must remain fixed. Missing or duplicate levels, changed source/parameters, an
inconsistent recorded verdict, or nonfinite data are rejected. Zero error does
not produce an infinite PASS rate: it makes the rate unresolved unless another
required gate already fails. Record outcome must agree with PASS/FAIL/UNRESOLVED.
Runtime is metadata, never a rerun or time assertion in CI.

`threshold_verification.mode="reviewer"` needs a public `review_ref` and explicitly
leaves code/criterion consistency to review. `literal-constants` additionally
names a hash-checked local validation source snapshot and nonempty bindings
from constant names to criterion JSON pointers. Only direct top-level literal
assignments are read using AST; no import or evaluation occurs. Each value/type
must match. This verifies only the declared bindings, not whether a real gate
uses that constant or whether other thresholds remain unbound. Gate wiring,
binding completeness, theory applicability and human consent always require
independent review; the report preserves that requirement.

The light regression lock is called a **surrogate test**, not a surrogate model.
It compares the adopted convergence sequence's coarsest level with its own
saved earlier output, not the analytic reference or finest solution. Its
`Assumes:` cites `surrogate test of adoption record <ID> (convergence level 1)`
with an actual approved ID. Missing protection is an explicit empty list, not
coverage. No documentation or MCP projection may infer verified status from a
missing, stale, failed or unapproved record.

`status` is `current` or `stale`; `covered_paths` declares the source paths whose changes invalidate current status. Preserve the stale reason and original evidence. `timing_s` values are descriptive records only and must never determine acceptance; the checker does not compare recorded times with limits.

Staleness uses `covered_sha256`, a nonempty map from covered repository paths to source-time SHA256. UTF-8 text normalizes CRLF and CR to LF; other bytes are hashed unchanged. Store a hash per matched file, with dependency paths selected by the record author when relevant to the claim; no root or package pyproject is forced. New, removed or changed covered files invalidate `current`; history ancestry is not required. An exact source SHA is provenance and need not exist as a local Git object; confirmed pre-rewrite identity loss follows the separate rule below. Optional `environment.ngsolve` records the actual experiment version, including pre-pin trials; it is not forced on unrelated packages. Only explicitly covered dependency metadata changes invalidate current status; the author is responsible for selecting a sufficient covered scope. Bundle the local script with the record so checking it needs no historical Git objects.

`source.source_only=true` explicitly identifies a source-only milestone. Otherwise at least one native or wheel SHA256 is required, in addition to the always-required result and entry artifact hashes. `superseded_by`, when set, must resolve to an existing different record in the same distribution. The checker never executes entry scripts or uses recorded timings to accept/reject a numerical result.

## Source provenance after history rewriting

`source.provenance="exact"` requires the known full execution `repo_sha`. A missing
local Git object or shallow checkout does not invalidate a known source identity;
the checker does not infer rewriting from Git object availability.

For evidence whose old execution SHA was lost through a confirmed history
rewrite, use `source.provenance="pre-rewrite"`, `source.repo_sha=null`, and
`status="stale"`. The reason must start with `pre-rewrite:` followed by an
explanation. This classification is an owner-reviewed provenance decision, not
an automatic fallback. It stays stale even when covered content matches today.
The old execution SHA and original-to-export mapping belong only in private
storage keyed by the stable record ID. Do not add a public `old_sha` field or
substitute today's HEAD for the historical execution. CI needs no private map.

Retain raw originals privately. Do not edit an old result or script in place to
make its provenance appear source-null. A public derived export must be separately
identified, hashed and reviewed; source-null does not authorize publishing its
originals. Embedded old execution identities and other private material require
owner review beyond this structural guard. Existing script/result/artifact,
clean-state and covered-source hash checks still apply to stale records. Do not
invent unavailable hashes, versions or dirty-state evidence to pass them.

There is no active v1/v2 parser or migration fallback. Update wrapper metadata to
v3 only with reviewed provenance and leave original result bytes untouched. A
fresh numerical run is a new exact-source record, not a relabeling of a stale one.

### Judgment rationale and reconsideration

An adoption includes `rationale`: a scoped inference, ID-bearing assumptions,
their evidence identifiers, invalidation conditions and recheck actions. Every
assumption must appear exactly once in `depends_on`; IDs are scoped by record ID.
Approval binds the canonical rationale hash as well as the criterion hash.
The checker returns assumption-to-record dependencies, enabling callers to join
(record ID, assumption ID) when a reference, coverage or mutation receipt fails.
It does not infer scientific truth or monitor external evidence automatically.
An invalidated premise requires independent reconsideration of every dependent
record; schema PASS alone cannot keep a scientific claim current. Evidence refs
are public identifiers, never private approval paths. No record is automatically
re-approved, re-run or relabelled by this structural checker.

## Three layers: adoption, protection, champion

```text
tests              coarsest convergence level -> saved-self surrogate test (every PR)
validation_test    approved convergence series -> adoption decision (trigger/release)
validation_records beyond adoption -> finest/highest-order one-off champion,
                   or a speed/scale record (schema-only CI; no scheduled replay)
```

A champion references its separate adoption decision by records-root-relative
ID and hash of that record.json, never by copying a second set of criteria.
Categories `accuracy-champion`, `speed`, and `scale` require a `champion` block;
an adoption decision is a separate non-champion record. It must contain a valid
PASS adoption, with matching feature; a stale origin remains historical evidence,
not current adoption certification. The reference hash prevents unnoticed edits.
Updating an origin to stale requires explicitly refreshing reviewed references.

The block names quantity/units, the selected level and selection reason, plus a
convergence table through that level: parameters, per-level source SHA, mesh
hash, result artifact path/hash, error and descriptive runtime_s. Each per-level
artifact contains `champion_level` equal to its declaration excluding `artifact`;
the enclosing result JSON contains `champion_result` equal to the block excluding
`adoption_record`. Thus no self-hashing artifact is required. Different levels
may legitimately have different source SHAs; review their comparability. Schema
checks identity and existence, not whether these are scientifically the best
results, higher resolution, or globally fastest/largest. Review the conditions,
accuracy and independent evidence before accepting a champion claim. No elapsed
time threshold, speed ranking or automatic adoption is inferred by the checker.

`supersedes` retains the predecessor and states improves/corrects/replaces with
reason. Covered-content hashes still require stale marking after code changes;
never delete an earlier champion because it is stale or superseded. Preserve
failed/partial runs outside a claimed best result. Do not periodically replay
champions: expensive one-off work beyond 30 minutes is not a PR or scheduled gate.
Trigger-driven validation uses validation_test, not re-running the champion.

A one-off extension may add a single new level: the linked adoption record retains
the preceding convergence sequence. Do not fabricate a second heavy run to fill
a table; provenance and comparability still require review.

## Reviewed decision bindings and reference resolution

Categories are the fixed roles milestone, adoption, accuracy-champion, speed,
and scale; free-form rankings such as fastest-ever are rejected. Category
adoption requires the adoption block, and an adoption block requires that role.
Approval additionally binds scope_sha256 = scope_hash(record), covering feature
and record ID; criterion/rationale hashes alone cannot transfer consent to
another entry or decision. approved_on must not be later than the record date.
These hashes and date checks still do not authenticate actual user consent.

Assumption evidence_refs resolve only to artifact:result or artifact:entry
(the record's hash-checked local artifacts), approval:<declared approval_ref>,
review:<declared threshold-verification review_ref>, or record:<existing ID>
in the same distribution's records root. Other references fail loudly. Record
references report observed bytes/hash and ID, not scientific validity or full
certification of that target. Review/approval identifiers remain declared
provenance requiring independent confirmation, never authentication by string.
Each break_conditions item declares an id and condition; invalidated_by lists
those declared IDs rather than unconstrained strings. The checker cannot decide
whether a failure condition has happened: external evidence and the inference
itself still require review. Namespace resolution is not proof of a premise.

For each protected_by Python module the checker parses its module docstring
without importing/executing it: Assumes must contain the exact phrase
surrogate test of adoption record <this ID> (convergence level 1). This checks
the citation, not production path coverage, mutation sensitivity or runtime.
Supersedes and superseded_by must agree in both directions; retain both records
and independently review the reason and comparability.

A rate passes its bound only with the direct comparison observed_rate >=
minimum_rate (no epsilon or rounding waiver). Floating-point boundary rounding
can conservatively reject mathematical equality; inspect the original evidence,
never silently loosen the human-owned bound to turn it into a pass.

### Champion and lifecycle fail-closed checks

A current champion requires a current origin; a stale origin requires the
champion to be marked stale too. Its original data remains retained. A champion
requires outcome=pass, the same quantity name/unit as its origin criterion, and
exactly the origin's variable keys at every level. Accuracy-champion errors must
be nonincreasing in declared table order and no worse than the finest origin
error (the best terminal error across the origin's declared convergence series).
This is a consistency guard, not a
new convergence-rate threshold, proof of globally best accuracy, or a timing
gate; scientific comparability and best-result selection still need review.

Adoption records require at least one protected_by lock. Each must be a .py
module with the exact record citation in Assumes. An adoption cannot use its
own record ID as record: evidence. The checker enforces these links without
importing tests or proving that they are collected, traverse the right code,
or kill a mutation. Those remain independent evidence/review obligations.
Contract checks must be able to fail on changed public behaviour; merely
asserting fixture existence or a constant flag is not evidence.
Future record dates relative to wall-clock
time are not checked (optional L5); approval-versus-record chronology and
calendar validity remain checked. Audit the actual dated provenance in review.
